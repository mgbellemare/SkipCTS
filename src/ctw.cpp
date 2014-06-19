/******************************
      Author: Joel Veness
        Date: 2011
******************************/

#include "ctw.hpp"

#include <vector>
#include <cassert>
#include <stack>
#include <iostream>
#include <cmath>

// boost includes
#include <boost/utility.hpp>


// enable both options below for better compression performance on text sources.
// disable both for vanilla CTW.

// do we use the zero redundancy estimator instead of the KT estimator?
static const bool UseZeroRedundancy = false;

// do we only perform weighting at byte boundaries in factored mode?
static const bool UseWeightingOnlyAtByteBoundaries = false;


// precompute some common logarithms
static const double log_point_five = std::log(0.5);
static const double log_quarter    = std::log(0.25);


/* create a new context node */
CTNode::CTNode() :
    m_log_prob_est(0.0),
    m_log_prob_weighted(0.0)
{
    m_count[0] = 0;    m_count[1] = 0;
    m_child[0] = NULL; m_child[1] = NULL;
}


/* update the weighted probabilities */
void CTNode::updateWeighted() {

    // computes P_w = log{0.5 * [P_kt + P_w0*P_w1]}
    double log_prob_on  = child(1) ? child(1)->logProbWeighted() : 0.0;
    double log_prob_off = child(0) ? child(0)->logProbWeighted() : 0.0;
    double log_one_plus_exp = log_prob_off + log_prob_on - logProbEstimated();

    // NOTE: no need to compute the log(1+e^x) if x is large, plus it avoids overflows
    if (log_one_plus_exp < 100.0) log_one_plus_exp = std::log(1.0 + std::exp(log_one_plus_exp));

    m_log_prob_weighted = log_point_five + logProbEstimated() + log_one_plus_exp;
}


/* process a new binary symbol */
void CTNode::update(bit_t b, bool skip) {

    // update the KT estimate and counts
    double log_kt_mul = logKTMul(b);
    m_log_prob_est += log_kt_mul;
    m_count[b]++;

    if (isLeaf()) {
        m_log_prob_weighted = logProbEstimated();
    } else {
        if (skip) {
            double log_prob_on  = child(1) ? child(1)->logProbWeighted() : 0.0;
            double log_prob_off = child(0) ? child(0)->logProbWeighted() : 0.0;
            m_log_prob_weighted = log_prob_on + log_prob_off;
        } else {
            updateWeighted();
        }
    }
}


/* is the current node a leaf node? */
bool CTNode::isLeaf() const {

    return child(0) == NULL && child(1) == NULL;
}


/* Krichevski-Trofimov estimated log probability accessor */
weight_t CTNode::logProbEstimated() const {

    if (UseZeroRedundancy) {
        if (m_count[0]+m_count[1] == 0) return 0.0;
        double rval = log_point_five + m_log_prob_est;
        if (m_count[0] == 0) rval = logAdd(log_quarter, rval);
        if (m_count[1] == 0) rval = logAdd(log_quarter, rval);
        return rval;
    }

    return m_log_prob_est;
}


/* logarithmic weighted probability estimate accessor */
weight_t CTNode::logProbWeighted() const {
    return m_log_prob_weighted;
}


/* child corresponding to a particular symbol */
const CTNode *CTNode::child(bit_t b) const {
    return m_child[b];
}


/* the number of times this context been visited */
int CTNode::visits() const {
    return m_count[0] + m_count[1];
}


/* compute the logarithm of the KT-estimator update multiplier */
double CTNode::logKTMul(bit_t b) const {

    static const double alpha = 0.5;
    static const double alpha2 = 2.0 * alpha;

    double kt_mul_numer = double(m_count[b]) + alpha;
    double kt_mul_denom = double(visits()) + alpha2;

    return std::log(kt_mul_numer / kt_mul_denom);
}


/* number of descendents of a node in the context tree */
size_t CTNode::size() const {

    size_t rval = 1;
    rval += child(0) ? child(0)->size() : 0;
    rval += child(1) ? child(1)->size() : 0;
    return rval;
}


/* create (if necessary) all of the nodes in the current context */
void ContextTree::createNodesInCurrentContext(const context_t &context) {

    CTNode **ctn = &m_root;

    for (size_t i = 0; i < context.size(); i++) {
        ctn = &((*ctn)->m_child[context[i]]);
        if (*ctn == NULL) {
            void *p = m_ctnode_pool.malloc();
            assert(p != NULL);  // TODO: make more robust
            *ctn = new (p) CTNode();
        }
    }
}


/* create a context tree of specified maximum depth and size */
ContextTree::ContextTree(history_t &history, size_t depth, int phase/*=-1*/) :
    m_ctnode_pool(sizeof(CTNode)),
    m_root(new (m_ctnode_pool.malloc()) CTNode()),
    m_phase(phase),
    m_depth(depth),
    m_history(history)
{
}


/* delete the context tree */
ContextTree::~ContextTree(void) {
    deleteCT(m_root);
}


/* recursively deletes the nodes in a context tree */
void ContextTree::deleteCT(CTNode *n) {

    if (n == NULL) return;

    if (n->m_child[0] != NULL) deleteCT(n->m_child[0]);
    if (n->m_child[1] != NULL) deleteCT(n->m_child[1]);

    m_ctnode_pool.free(n);
}


/* compute the current binary context */
void ContextTree::getContext(const history_t &h, context_t &context) const {

    context.clear();
    for (size_t i=0; i < m_depth; ++i) {
        context.push_back(h[h.size()-i-1]);
    }
}


/* updates the context tree with a single bit */
void ContextTree::update(bit_t b) {

    // compute the current context
    context_t context;
    context.reserve(m_depth);
    getContext(m_history, context);

    // 1. create new nodes in the context tree (if necessary)
    createNodesInCurrentContext(context);

    // 2. walk down the tree to the relevant leaf, saving the path as we go
    std::stack<CTNode *, std::vector<CTNode *> > path;
    path.push(m_root); // add the empty context
    CTNode *ctn = m_root;
    for (size_t i = 0; i < context.size(); i++) {
        ctn = ctn->m_child[context[i]];
        path.push(ctn);
    }

    // 3. update the probability estimates from the leaf node back up to the root
    int index = static_cast<int>(m_depth);
    for (; !path.empty(); path.pop()) {
        bool skip = UseWeightingOnlyAtByteBoundaries && m_phase > -1 &&
                    (index % 8) != m_phase && index != 0;
        path.top()->update(b, skip);
        index--;
    }

    // 4. update the history
    m_history.push_back(b != 0);
}


/* the probability of seeing a particular symbol next */
double ContextTree::prob(bit_t b) {

    typedef std::pair<CTNode *, CTNode> ctpair_t;

    double before = logBlockProbability();

    // compute the current context
    context_t context;
    getContext(m_history, context);

    // 1. record newly added or modified nodes
    std::vector<CTNode *> created;
    std::vector<ctpair_t> modified;

    CTNode **ctnp = &m_root;
    modified.push_back(ctpair_t(m_root, *m_root));
    for (size_t i = 0; i < context.size(); i++) {
        ctnp = &((*ctnp)->m_child[context[i]]);
        if (*ctnp == NULL) {
            void *p = m_ctnode_pool.malloc();
            assert(p != NULL);  // TODO: make more robust
            *ctnp = new (p) CTNode();
            created.push_back(*ctnp);
        } else {
            modified.push_back(ctpair_t(*ctnp, **ctnp));
        }
    }

    // 2. walk down the tree to the relevant leaf, saving the path as we go
    std::stack<CTNode *, std::vector<CTNode *> > path;
    path.push(m_root); // add the empty context
    CTNode *ctn = m_root;
    for (size_t i = 0; i < context.size(); i++) {
        ctn = ctn->m_child[context[i]];
        path.push(ctn);
    }

    // 3. update the probability estimates from the leaf node back up to the root
    int index = static_cast<int>(m_depth);
    for (; !path.empty(); path.pop()) {
        bool skip = UseWeightingOnlyAtByteBoundaries && m_phase > -1 &&
                    (index % 8) != m_phase && index != 0;
        path.top()->update(b, skip);
        index--;
    }

    double rval = std::exp(logBlockProbability() - before);

    // now revert the changes
    for (size_t i=0; i < created.size(); ++i) m_ctnode_pool.free(created[i]);
    for (size_t i=0; i < modified.size(); ++i) *modified[i].first = modified[i].second;

    return rval;
}


/* the depth of the context tree */
size_t ContextTree::depth() const {

    return m_depth;
}


/* number of nodes in the context tree */
size_t ContextTree::size(void) const {

    return m_root->size();
}


/* recover the memory used by a node */
void ContextTree::reclaimMemory(CTNode *n) {

    m_ctnode_pool.free(n);
}


/* the logarithm of the block probability of the whole sequence */
double ContextTree::logBlockProbability(void) const {

    return m_root->logProbWeighted();
}


