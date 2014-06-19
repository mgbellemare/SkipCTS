/******************************
      Author: Joel Veness
        Date: 2013
******************************/

#include "skipcts.hpp"
#include "fastmath.hpp"

#include <boost/random.hpp>


// adaptive KT estimator parameters
static const bool UseDiscounting = true;
static const count_t Gamma       = 0.02f;
static const count_t Discount    = 1.0f - Gamma;
static const count_t KT_Alpha = 0.0625f;
static const count_t KT_Alpha2 = KT_Alpha + KT_Alpha;


// initialise weights for switching model
static const double SplitPrior        = 0.925;
static const double LogSplitPrior   = std::log(SplitPrior);
static const double LogStopPrior    = std::log(1.0 - SplitPrior);

static const double LogOneHalf = std::log(0.5);

struct prior_t {
    double stop;
    double split;
    double skip;
};


// prior weights for either splitting, stopping or skipping
static const prior_t SplitSkipPrior[] = { 
    { LogStopPrior,    LogSplitPrior,   std::log(0.0)   },
    { std::log(0.075), std::log(0.85),  std::log(0.075) },
    { std::log(0.075), std::log(0.85),  std::log(0.075) },
    { std::log(0.075), std::log(0.85),  std::log(0.075) },
    { std::log(0.075), std::log(0.85),  std::log(0.075) },
    { std::log(0.075), std::log(0.85),  std::log(0.075) },
};


// zobrist key generation for context hashing
typedef boost::mt19937 randsrc_t;
typedef boost::variate_generator<randsrc_t &, boost::uniform_int<uint64_t> > zobrng_t;
static randsrc_t randsrc(666);
static boost::uniform_int<uint64_t> uint64_random_range(0, std::numeric_limits<uint64_t>::max());
static zobrng_t zobrng(randsrc, uint64_random_range);


// precomputed static tables
zobhash_t SkipCTS::s_zobtbl[MaxDepth][2];
weight_t SkipCTS::s_log_tbl[LogTblSize];


/* adds two numbers represented in the logarithmic domain */
inline static double fast_logadd(double log_x, double log_y) {

    if (log_x < log_y) {
        return fast_jacoblog(log_y - log_x) + log_x;
    } else {
        return fast_jacoblog(log_x - log_y) + log_y;
    }
}


/* skip nodes */
SkipNode::SkipNode() :
    m_log_prob_est(LogStopPrior),
    m_log_prob_split(LogSplitPrior),
    m_log_prob_weighted(0.0),
    m_log_skip_lik(NULL),
    m_depth(-1),
    m_buf(0.0)
{
    m_count[0] = 0.0f;
    m_count[1] = 0.0f;
}


SkipNode::~SkipNode() {
    delete [] m_log_skip_lik;
}


/* compute the logarithm of the KT-estimator update multiplier */
inline double SkipNode::logKTMul(bit_t b) const {

    return fast_log(ktMul(b));
}


/* compute the logarithm of the KT-estimator update multiplier */
inline double SkipNode::ktMul(bit_t b) const {

    count_t kt_mul_numer = m_count[b] + KT_Alpha;
    count_t kt_mul_denom = m_count[0] + m_count[1] + KT_Alpha2;

    return kt_mul_numer / kt_mul_denom;
}


/* update the KT estimates */
inline void SkipNode::updateKT(bit_t b, double log_est_mul) {

    if (UseDiscounting) {
        m_count[0] *= Discount;
        m_count[1] *= Discount;
    }
    m_count[b]++;
}

       
/* precompute the zobrist hash keys used for context hashing. */
void SkipCTS::initZobrist() {

    static bool init = false;

    if (!init) {
        for (int i=0; i < SkipCTS::MaxDepth; i++) {
            s_zobtbl[i][0] = zobrng();
            s_zobtbl[i][1] = zobrng();
        }
        init = true;
    }
}


/* initialise a table of precomputed logarithms */
void SkipCTS::initLogTbl() {

    static bool init = false;

    if (!init) {
        for (int i=0; i < SkipCTS::LogTblSize; i++) {
            s_log_tbl[i] = std::log(static_cast<double>(i));
        }
        init = true;
    }
}


/* initialise the hash deltas to save having to recompompute each 
   zobrist hash key from scratch. */
void SkipCTS::initHashDeltas() {
    
    // compute the hash deltas
    indices_t delta;

    for (int i=m_depth; i >= 0; i--) {
        
        indices_list_t &il = m_indices[i];
        
        for (int j=0; j < il.size(); j++) {

            indices_t &idxs = il[j];                
            indices_t diff(m_depth);

            indices_t::iterator it = std::set_symmetric_difference(
                delta.begin(), delta.end(), idxs.begin(),
                idxs.end(), diff.begin()
            );
            diff.resize(it - diff.begin());
            delta = idxs;
            idxs = diff;
        }
    }
}


/* initialise auxilary information such as number of remaining skips,
   and the depth of the context. */
void SkipCTS::initAuxInfo() {
    
    // initialise the auxilary information
    for (int i=0; i <= m_depth; i++) {
        
        m_auxinfo.push_back(aux_info_list_t());
        indices_list_t &il = m_indices[i];
        
        for (int j=0; j < il.size(); j++) {
            
            m_auxinfo.back().push_back(aux_info_t());
            
            // precompute last index
            m_auxinfo.back().back().last_idx   = -1;
            if (!il[j].empty()) 
                m_auxinfo.back().back().last_idx = static_cast<int>(il[j].back());
            
            // precompute number of skips remaining
            m_auxinfo.back().back().skips_left = m_skips;
            for (size_t k=0; k < il[j].size(); k++) {
                if ((k == 0 &&  il[j][k] > 0) || (k > 0 &&  il[j][k] >  il[j][k-1]+1)) 
                    m_auxinfo.back().back().skips_left--;
            }
        }
    }
}


/* compute the indices for the zobrist hashes, on a per depth basis */
void SkipCTS::initIndices() {
    
    for (size_t i=0; i <= m_depth; i++)
        m_indices.push_back(indices_list_t());

    skip_contexts_t::const_iterator it = m_skip_contexts.begin();
    for ( ; it != m_skip_contexts.end(); ++it) {
  
        const std::string &s = *it;
        int depth = static_cast<int>(std::count(s.begin(), s.end(), 'b'));
        m_indices[depth].push_back(indices_t());
        
        for (size_t j=0; j < s.length(); j++) {
            if (s[j] == 'b') m_indices[depth].back().push_back(j);
        }
    }
}


/* precomputes the context indices */
void SkipCTS::initContexts() {

    // generate contexts
    m_skip_contexts.clear();
    genBoundedSkipContexts("", m_depth, m_skips);

    initIndices();

    m_skip_contexts.clear();
}


/* create the skip context tree */
SkipCTS::SkipCTS(history_t &history, int depth, int max_skips, size_t log2_slots) :
    m_history(history),
    m_depth(depth),
    m_skips(max_skips),
    m_log_skip_preds(MaxDepth),
    m_nodes(new SkipNode[size_t(1) << log2_slots]),
    m_mask((size_t(1) << log2_slots)-1)
{
    initZobrist();
    initContexts();
    initAuxInfo();
    initHashDeltas();
    initLogTbl();
}


SkipCTS::~SkipCTS() {
    delete [] m_nodes;
}


/* the logarithm of the probability of all processed experience */
double SkipCTS::logBlockProbability() const {

    return getNode(0, 0, numSubmodels(-1, m_skips)).logProbWeighted();
}


/* compute the switching rate for a given time t */
double SkipCTS::switchRate(size_t t) const {

    return 1.0 / double(t - m_depth + 3);
}


int SkipCTS::numSubmodels(int position, int skips) const {

    int n = m_depth - position;
    return (skips == 0) ? std::min(n, 2) : n;
}


/* the probability of seeing a particular symbol next */
double SkipCTS::prob(bit_t b) {

    // We proceed as with update(), except that we keep track of the symbol probability at
    // each node instead of actually updating parameters
    getContext();

    zobhash_t hash = 0;
    int skips_left, last_idx;
    double symbolLogProb = LogOneHalf; 

    // propagate the symbol probability from leaves to root 
    for (int i=m_depth; i >= 0; i--) {
        
        // update the KT statistics, then the weighted 
        // probability for every node on this level
        const indices_list_t &il = m_indices[i];
        
        for (int j=0; j < il.size(); j++) {

            getContextInfo(hash, il[j]);
            skips_left = m_auxinfo[i][j].skips_left;
            last_idx   = m_auxinfo[i][j].last_idx;

            // update the node
            int n_submodels = numSubmodels(last_idx, skips_left);

            SkipNode &n = getNode(hash, i, n_submodels);

            // handle the stop case
            double log_est_mul = n.logKTMul(b);
            if (n_submodels == 1) {
                n.m_buf = symbolLogProb = log_est_mul;
                continue;
            }
                 
            // Here we rely on the property that log_prob_est and the like are unnormalized
            // posteriors. we accumulate the symbol log probability under 'symbolLogProb' 
            symbolLogProb = n.m_log_prob_est + log_est_mul;
 
            // handle the split case
            zobhash_t delta = s_zobtbl[last_idx+1][m_context[last_idx+1]];
            const SkipNode &nn = getNode(hash ^ delta, i+1, numSubmodels(last_idx+1, skips_left));
            // recall that m_buf contains the symbol probability at the child node
            double log_split_pred = nn.m_buf;
            symbolLogProb = fast_logadd(symbolLogProb, n.m_log_prob_split + log_split_pred);

            // handle the skipping case
            if (n_submodels > 2) {
                
                // if we did not yet allocate these models, this node must never have been
                // updated. We assume (perhaps incorrectly) that none of this node's children
                // exist and pretend they return a symbol probability of 0.5 
                if (n.m_log_skip_lik == NULL) {

                    const prior_t &p = SplitSkipPrior[skips_left];
                    symbolLogProb = fast_logadd(symbolLogProb, p.skip + LogOneHalf); 
                }
                
                // mix in the symbol probability from the skipping models
                else for (int k=last_idx+2; k < m_depth; k++) { 
                    
                    zobhash_t h = hash ^ s_zobtbl[k][m_context[k]];
                    SkipNode &sn = getNode(h, i+1, numSubmodels(k, skips_left - 1));
                    
                    double log_skip_pred = sn.m_buf;
                    int z = k - last_idx - 2;
                    symbolLogProb = fast_logadd(symbolLogProb, n.m_log_skip_lik[z] + log_skip_pred);
                }
            }

            // Finally we normalize by the mixture probability at this node
            symbolLogProb -= n.m_log_prob_weighted;
            n.m_buf = symbolLogProb;

            assert(n.m_buf < 0.0);
        }
    }
   
    // our scheme assumes that the last node processed is the root; the variable 'symbolLogProb'
    // contains its symbol probability 
    return fast_exp(symbolLogProb);
}


/* gets the node's index into the hash table */
inline SkipNode &SkipCTS::getNode(zobhash_t hash, int depth, int submodels) const {

    size_t key = static_cast<size_t>(hash & m_mask);

    // do a linear scan till we find either an empty slot, or
    // a populated slot with matching depth
    do {
        if (m_nodes[key].m_depth == -1) break;
        if (m_nodes[key].m_depth == depth && m_nodes[key].m_submodels == submodels) break;
        key = (key + 1) & m_mask;
    } while (true);

    // mark the slot as used
    m_nodes[key].m_depth = depth;
    m_nodes[key].m_submodels = submodels;

    return m_nodes[key];
}


/* performs the update operation to maintain a switching posterior. */
void SkipCTS::posteriorUpdate(
    SkipNode &n, double log_scale, double log_alpha, 
    double log_K, double log_mul, double &log_post
) const {

    log_post = log_scale + 
        fast_logadd(
            log_alpha + n.m_log_prob_weighted,
            log_K + log_post + log_mul
        );
}


/* lazy allocation of skipping posterior weights. */
void SkipCTS::lazyAllocate(SkipNode &n, int n_submodels, int skips_left) const {

    n.m_log_skip_lik = new weight_t[n_submodels-2];

    const prior_t &p = SplitSkipPrior[skips_left];
    n.m_log_prob_est    = p.stop;
    n.m_log_prob_split  = p.split;
    for (int k=0; k < n_submodels-2; k++) {
        n.m_log_skip_lik[k] = p.skip - s_log_tbl[n_submodels-2];  
    }
}


/* compute the information needed to update the current context stats. */
void SkipCTS::getContextInfo(zobhash_t &hash, const indices_t &idxs) const {

    // compute the hash
    for (int k=0; k < idxs.size(); k++) {
        size_t x = idxs[k];
        hash ^= s_zobtbl[x][m_context[x]];
    }
}


/* process a new piece of sensory experience */
void SkipCTS::update(bit_t b) {

    getContext();

    double alpha = switchRate(m_history.size());
    double log_alpha = fast_log(alpha);
    
    zobhash_t hash = 0;
    int skips_left, last_idx;

    // update nodes from deepest to shallowest
    for (int i=m_depth; i >= 0; i--) {
        
        // update the KT statistics, then the weighted 
        // probability for every node on this level
        const indices_list_t &il = m_indices[i];
        
        for (int j=0; j < il.size(); j++) {

            m_log_skip_preds.clear();

            getContextInfo(hash, il[j]);
            skips_left = m_auxinfo[i][j].skips_left;
            last_idx   = m_auxinfo[i][j].last_idx;

            // update the node
            int n_submodels = numSubmodels(last_idx, skips_left);

            SkipNode &n = getNode(hash, i, n_submodels);
            n.m_buf = n.m_log_prob_weighted;

            // lazy allocation of skipping prior weights
            if (n_submodels > 2 && n.m_log_skip_lik == NULL)
                lazyAllocate(n, n_submodels, skips_left);
    
            // handle the stop case
            double log_est_mul = n.logKTMul(b);
            if (n_submodels == 1) {
                n.updateKT(b, log_est_mul);
                n.m_log_prob_weighted += log_est_mul;
                n.m_buf = log_est_mul;
                continue;
            }
                 
            double log_acc = n.m_log_prob_est + log_est_mul;
            n.updateKT(b, log_est_mul);
                
            // handle the split case
            zobhash_t delta = s_zobtbl[last_idx+1][m_context[last_idx+1]];
            const SkipNode &nn = getNode(hash ^ delta, i+1, numSubmodels(last_idx+1, skips_left));
            double log_split_pred = nn.m_buf;
            log_acc = fast_logadd(log_acc, n.m_log_prob_split + log_split_pred);

            // handle the skipping case
            if (n_submodels > 2) {
                
                // update the skipping models
                for (int k=last_idx+2; k < m_depth; k++) { 
                    
                    zobhash_t h = hash ^ s_zobtbl[k][m_context[k]];
                    SkipNode &sn = getNode(h, i+1, numSubmodels(k, skips_left - 1));
                    
                    double log_skip_pred = sn.m_buf;
                    m_log_skip_preds.push_back(log_skip_pred);
                    int z = k - last_idx - 2;
                    log_acc = fast_logadd(log_acc, n.m_log_skip_lik[z] + log_skip_pred);
                }
            }

            // store the weighted probability
            n.m_log_prob_weighted = log_acc;

            assert(n.m_log_prob_weighted < n.m_buf);
            // Store the *difference* in log probability in m_buf
            n.m_buf = n.m_log_prob_weighted - n.m_buf;

            updatePosteriors(n, n_submodels, alpha, log_alpha, log_est_mul, log_split_pred);
        }
    }

    m_history.push_back(b != 0);
}


/* update the switching posterior weights */
void SkipCTS::updatePosteriors(SkipNode &n, int n_submodels, double alpha, 
    double log_alpha, double log_stop_mul, double log_split_mul) {

    // update switching log-posteriors
    double dn        =  static_cast<double>(n_submodels);
    double K         =  (1.0 - alpha) * dn - 1.0;
    double log_K     =  fast_log(K);
    double log_scale = -s_log_tbl[n_submodels-1];

    posteriorUpdate(n, log_scale, log_alpha, log_K, log_stop_mul, n.m_log_prob_est);
    posteriorUpdate(n, log_scale, log_alpha, log_K, log_split_mul, n.m_log_prob_split);
    for (int k=0; k < m_log_skip_preds.size(); k++) {
        posteriorUpdate(n, log_scale, log_alpha, log_K, m_log_skip_preds[k], n.m_log_skip_lik[k]);
    }
}


/* generate possible context strings with bounded # of skips */
void SkipCTS::genBoundedSkipContexts(std::string buf, int depth, int skips) {

    m_skip_contexts.push_back(buf);

    if (depth == 0) return;

    if (skips > 0) {
        std::string skipstr = "*";
        for (int i=0; i < depth-1; i++) {
            genBoundedSkipContexts(buf + skipstr + "b", depth - (i+1) - 1, skips-1);
            skipstr.append("*");
        }
    }

    genBoundedSkipContexts(buf + "b", depth - 1, skips);
}


/* prints a set of skips contexts */
void SkipCTS::printContexts(skip_contexts_t &l) const {

    skip_contexts_t::const_iterator it = l.begin();
    for (; it != l.end(); ++it) {
        std::cout << *it << std::endl;
    }
}


/* compute the current binary context */
inline void SkipCTS::getContext() {

    size_t offset = m_history.size();
    m_context.clear();

    for (size_t i=0; i < m_depth; ++i) {
        m_context.push_back(m_history[offset-i-1]);
    }
}

