/******************************
      Author: Joel Veness
        Date: 2011
******************************/

#include "cts.hpp"
#include "fastmath.hpp"

#include <vector>
#include <cassert>
#include <stack>
#include <iostream>
#include <cmath>

// boost includes
#include <boost/utility.hpp>

// do we use discounting with the KT-estimator?
static const bool UseDiscounting = true;
static const count_t Gamma       = count_t(0.02);
static const count_t Discount    = count_t(1.0) - Gamma;
static const count_t KT_Alpha = 0.0625f;
static const count_t KT_Alpha2 = KT_Alpha + KT_Alpha;

// initialise weight for switching model
static const double SwitchPrior        = 0.925;

// precompute some commonly used logarithms
static const double log_switch_prior   = std::log(SwitchPrior);
static const double log_kt_prior       = std::log(1.0 - SwitchPrior);

// do we use unique path pruning?
static const bool UseUniquePathPruning       = true;
static const bool StrictModePathPruning      = true;

// do we use fast numerical approximations
static const bool UseFastLog         = true;
static const bool UseFastExp         = true;
static const bool UseFastJacobianLog = true;


/* compute a natural logarithm */
inline double ctsLog(double x) {
    return UseFastLog ? fast_log(x) : std::log(x);
}


/* compute an exponential */
inline double ctsExp(double x) {
    return UseFastExp ? fast_exp(x) : std::exp(x);
}


/* given log(x) and log(y), compute log(x+y). uses the following identity:
   log(x + y) = log(x) + log(1 + y/x) = log(x) + log(1+exp(log(y)-log(x))) */
inline double ctsLogAdd(double log_x, double log_y) {

    if (UseFastJacobianLog) {

        if (log_x < log_y) {
            return fast_jacoblog(log_y - log_x) + log_x;
        } else {
            return fast_jacoblog(log_x - log_y) + log_y;
        }

    } else {
        if (log_x < log_y) {
            return std::log(1.0 + std::exp(log_y - log_x)) + log_x;
        } else {
            return std::log(1.0 + std::exp(log_x - log_y)) + log_y;
        }
    }
}


/* display an SNode */
std::ostream &operator<<(std::ostream &o, const SNode &sn) {

    o << "(est: " << sn.m_log_prob_est
      << ", weighted: " << sn.m_log_prob_weighted
      << ", b: " << sn.m_log_b
      << ", s: " << sn.m_log_s
      << ", counts: " << sn.m_count[0] << "/" << sn.m_count[1]
      << ", pidx: " << sn.m_pruned_idx
      << ", children: " << sn.m_child[0] << "/" << sn.m_child[1]
      << ")";

    return o;
}


/* create a new switching node */
SNode::SNode(int depth) :
    m_log_prob_est(0.0),
    m_log_prob_weighted(0.0),
    m_log_b(log_kt_prior),
    m_log_s(log_switch_prior),
    m_pruned_idx(-1)
{
    m_count[0] = 0;    m_count[1] = 0;
    m_child[0] = NULL; m_child[1] = NULL;
}


/* create a new switching node from a pruned node */
SNode::SNode(const SNode &rhs, int pindx) :
    m_log_prob_est(rhs.m_log_prob_est),
    m_log_prob_weighted(rhs.m_log_prob_weighted),
    m_log_b(rhs.m_log_b),
    m_log_s(rhs.m_log_s),
    m_pruned_idx(pindx)
{
    m_count[0] = rhs.m_count[0];
    m_count[1] = rhs.m_count[1];
    m_child[0] = NULL;
    m_child[1] = NULL;
}


/* process a new binary symbol, with switching rate alpha, and blend 1-2*alpha */
void SNode::update(bit_t b, double log_alpha, double log_blend, double log_split_mul) {

    // update the KT estimate and counts
    double log_est_mul = logKTMul(b);
    m_log_prob_est += log_est_mul;
    if (UseDiscounting) {
        m_count[0] *= Discount;
        m_count[1] *= Discount;
    }
    m_count[b]++;

    if (isLeaf()) {

        if (m_pruned_idx >= 0) {
            // if we have pruned, we know weighted_pr == estimated_pr
            m_log_prob_weighted = logProbEstimated();
            m_log_b = ctsLogAdd(log_alpha + m_log_prob_weighted, log_blend + m_log_b + log_est_mul);
            m_log_s = ctsLogAdd(log_alpha + m_log_prob_weighted, log_blend + m_log_s + log_est_mul);
        } else {
            m_log_prob_weighted = logProbEstimated();
        }

        return;
    }

    m_log_prob_weighted = ctsLogAdd(m_log_b + log_est_mul, m_log_s + log_split_mul);
    m_log_b = ctsLogAdd(log_alpha + m_log_prob_weighted, log_blend + m_log_b + log_est_mul);
    m_log_s = ctsLogAdd(log_alpha + m_log_prob_weighted, log_blend + m_log_s + log_split_mul);
}


/* compute the result of an update call non-destructively */
double SNode::updateNonDestructive(bit_t b, double c_weighted, double c_weighted_old) const {

    // compute the KT estimate
    double log_est_mul = logKTMul(b);
    double log_prob_est = m_log_prob_est;
    log_prob_est += log_est_mul;

    if (isLeaf()) return log_prob_est;

    double log_split_mul = c_weighted - c_weighted_old;
    return ctsLogAdd(m_log_b + log_est_mul, m_log_s + log_split_mul);
}


/* is the current node a leaf node? */
bool SNode::isLeaf() const {

    return child(0) == NULL && child(1) == NULL;
}


/* Krichevski-Trofimov estimated log probability accessor */
weight_t SNode::logProbEstimated() const {

    return m_log_prob_est;
}


/* logarithmic weighted probability estimate accessor */
weight_t SNode::logProbWeighted() const {
    return m_log_prob_weighted;
}


/* child corresponding to a particular symbol */
const SNode *SNode::child(bit_t b) const {

    return m_child[b];
}


/* the number of times this context been visited */
count_t SNode::visits() const {

    return m_count[0] + m_count[1];
}


/* compute the logarithm of the KT-estimator update multiplier */
inline double SNode::logKTMul(bit_t b) const {

    count_t kt_mul_numer = m_count[b] + KT_Alpha;
    count_t kt_mul_denom = m_count[0] + m_count[1] + KT_Alpha2;

    return ctsLog(kt_mul_numer / kt_mul_denom);
}


/* number of descendents of a node in the context tree */
size_t SNode::size() const {

    size_t rval = 1;
    rval += child(0) ? child(0)->size() : 0;
    rval += child(1) ? child(1)->size() : 0;
    return rval;
}


/* determine whether two contexts are identical */
bool SwitchingTree::contextsEqual(const context_t &lhs, const context_t &rhs) {

    assert(lhs.size() == rhs.size());

    for (size_t i=0; i < lhs.size(); i++) {
        if (lhs[i] != rhs[i]) return false;
    }

    return true;
}


/* create (if necessary) all of the nodes in the current context */
void SwitchingTree::createNodesInCurrentContext(const context_t &context) {

    SNode **ctn = &m_root;

    if (!UseUniquePathPruning) {

        for (size_t i = 0; i < context.size(); i++) {
            ctn = &((*ctn)->m_child[context[i]]);
            if (*ctn == NULL) {
                void *p = m_ctnode_pool.malloc();
                assert(p != NULL);  // TODO: make more robust
                *ctn = new (p) SNode(static_cast<int>(i));
            }
        }
        return;
    }

    // unique path pruning - only create nodes that are needed!
    for (size_t i = 0; i < context.size(); i++) {

        SNode *n = *ctn;
        // if we encountered a node with pruning, restore the statistics
        if (n->m_pruned_idx >= 0) {

            // get the pruned context
            getContext(m_history, m_pcontext, n->m_pruned_idx);

            // strict unique path pruning check: if contexts are identical,
            // don't create _any_ more new nodes!
            if (StrictModePathPruning && contextsEqual(m_context, m_pcontext)) return;

            // now expand the old context out till it is unique again,
            // copying in the old relevant information
            SNode **pctn = ctn;
            int pidx = n->m_pruned_idx;
            for (size_t j=i; j < m_pcontext.size(); j++) {

                (*pctn)->m_pruned_idx = -1;

                pctn = &((*pctn)->m_child[m_pcontext[j]]);
                void *p = m_ctnode_pool.malloc();
                assert(p != NULL);  // TODO: make more robust
                *pctn = new (p) SNode(*n, j == m_pcontext.size()-1 ? -1 : pidx);

                if (m_pcontext[j] != context[j]) break;
            }
        }

        // create new node
        ctn = &((*ctn)->m_child[context[i]]);
        if (*ctn == NULL) {
            void *p = m_ctnode_pool.malloc();
            assert(p != NULL);  // TODO: make more robust
            *ctn = new (p) SNode(static_cast<int>(i));
            if (i+1 < m_context.size())
                (*ctn)->m_pruned_idx = static_cast<int>(m_history.size());
            break;
        }
    }
}


/* create a context tree of specified maximum depth and size */
SwitchingTree::SwitchingTree(history_t &history, size_t depth, int phase/*=-1*/) :
    m_ctnode_pool(sizeof(SNode)),
    m_root(new (m_ctnode_pool.malloc()) SNode(0)),
    m_phase(phase),
    m_depth(depth),
    m_history(history),
    m_prob_cache(-1)
{
}


/* delete the context tree */
SwitchingTree::~SwitchingTree(void) {
    deleteCT(m_root);
}


/* recursively deletes the nodes in a context tree */
void SwitchingTree::deleteCT(SNode *n) {

    if (n == NULL) return;

    if (n->m_child[0] != NULL) deleteCT(n->m_child[0]);
    if (n->m_child[1] != NULL) deleteCT(n->m_child[1]);

    m_ctnode_pool.free(n);
}


/* compute the current binary context */
void SwitchingTree::getContext(const history_t &h, context_t &context, int idx /* = -1 */) {

    size_t offset = idx < 0 ? h.size() : size_t(idx);
    context.clear();

    for (size_t i=0; i < m_depth; ++i) {
        context.push_back(h[offset-i-1]);
    }
}


/* computes the context, creates relevant nodes and determine the path to update */
void SwitchingTree::makeContextAndPath() {

    // compute the current context
    getContext(m_history, m_context);

    // 1. create new nodes in the context tree (if necessary)
    createNodesInCurrentContext(m_context);

    // 2. walk down the tree to the relevant leaf, saving the path as we go
    m_path.clear(); m_log_old_weights.clear();
    SNode *ctn = m_root;
    m_log_old_weights.push_back(m_root->logProbWeighted());
    m_path.push_back(m_root); // add the empty context

    for (size_t i = 0; i < m_context.size(); i++) {
        ctn = ctn->m_child[m_context[i]];
        if (ctn == NULL) break;
        m_log_old_weights.push_back(ctn->logProbWeighted());
        m_path.push_back(ctn);
    }
}


/* compute the switching rate for a given time t */
double SwitchingTree::switchRate(size_t t) const {

    return 1.0 / double(t - m_depth + 3);
}


/* updates the context tree with a single bit */
void SwitchingTree::update(bit_t b) {

    // avoid recomputing context and path if prob() just called
    if (m_prob_cache != static_cast<int>(m_history.size())) makeContextAndPath();

    // 3. update the probability estimates from the leaf node back up to the root
    double alpha = switchRate(m_history.size());

    double log_alpha = ctsLog(alpha);
    double log_blend = ctsLog(1.0 - 2.0*alpha);
    double log_split_mul = 0.0;

    SNode **sn = &m_path[m_path.size()-1];
    SNode *snc = NULL;
    size_t c = m_path.size()-1;

    while (!m_path.empty()) {
        (*sn)->update(b, log_alpha, log_blend, log_split_mul);
        log_split_mul = (*sn)->logProbWeighted() - m_log_old_weights[c];
        snc = *sn;
        sn--;
        c--;
        m_path.pop_back();
    }

    // 4. update the history
    m_history.push_back(b != 0);
}


/* the probability of seeing a particular symbol next */
double SwitchingTree::prob(bit_t b) {

    double before = logBlockProbability();

    makeContextAndPath();

    // remember for later so we can avoid calling makeContextAndPath
    // if update() is called immediately after this
    m_prob_cache = static_cast<int>(m_history.size());

     // 3. compute the probability estimates from the leaf node back up to the root
    double c_weighted = 0.0, c_old_weighted = 0.0;
    SNode *snc = NULL;
    SNode **sn = &m_path[m_path.size()-1];
    size_t c = 0;
    while (c < m_path.size()) {
        c_weighted = (*sn)->updateNonDestructive(b, c_weighted, c_old_weighted);
        c_old_weighted =  m_log_old_weights[m_path.size()-c-1];
        snc = *sn;
        sn--;
        c++;
    }

    return ctsExp(c_weighted - before);
}


/* the depth of the context tree */
size_t SwitchingTree::depth() const {

    return m_depth;
}


/* number of nodes in the context tree */
size_t SwitchingTree::size(void) const {

    return m_root->size();
}


/* recover the memory used by a node */
void SwitchingTree::reclaimMemory(SNode *n) {

    m_ctnode_pool.free(n);
}


/* the logarithm of the block probability of the whole sequence */
double SwitchingTree::logBlockProbability(void) const {

    return m_root->logProbWeighted();
}




