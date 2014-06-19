#ifndef __SKIPCTS_HPP__
#define __SKIPCTS_HPP__

/******************************
      Author: Joel Veness
        Date: 2013
******************************/

#include "common.hpp"

#include <set>
#include <string>
#include <map>

// boost includes
#include <boost/utility.hpp>

typedef uint64_t zobhash_t;
typedef std::vector<std::string> skip_contexts_t;
typedef std::vector<size_t> indices_t;
typedef std::vector<indices_t> indices_list_t;

struct aux_info_t {
    int skips_left;
    int last_idx;
    /** The total number of context variables split on */
    int num_splits;
};
typedef std::vector<aux_info_t> aux_info_list_t;

// skipping context tree node
class SkipNode {

    friend class SkipCTS;
    friend std::ostream &operator<<(std::ostream &o, const SkipNode &);

    public:

        SkipNode();
        ~SkipNode();

        // update the KT estimates
        void updateKT(bit_t b, double log_est_mul);

        // update weighted estimates at leaf nodes
        void updateWeightedLeaf();

        // log weighted blocked probability
        weight_t logProbWeighted() const { return m_log_prob_weighted; }

        // log KT estimated probability
        weight_t logProbEstimated() const { return m_log_prob_est; }

        // log Split probability
        weight_t logProbSplit() const { return m_log_prob_split; }

    private:

        // compute the KT-estimator update multiplier
        double logKTMul(bit_t b) const;
        double ktMul(bit_t b) const;

        weight_t m_log_prob_est;
        weight_t m_log_prob_split;
        weight_t m_log_prob_weighted;
        weight_t m_buf;

        weight_t *m_log_skip_lik;

        // one slot for each binary value
        count_t m_count[2];

        short m_depth;
        short m_submodels;
};


// skip context tree switching
class SkipCTS : public Compressor, boost::noncopyable {

    static const int MaxLogDepth = 9;
    static const int MaxDepth = (1 << MaxLogDepth);
    static const int LogTblSize = MaxDepth * 2;

    public:

        SkipCTS(history_t &history, int depth, int max_skips=1, size_t log2_slots=24);
        ~SkipCTS();

        // file extension
        virtual const char *fileExtension() const { return "skipcts"; }

        // the logarithm of the probability of all processed experience
        virtual double logBlockProbability() const;

        // the probability of seeing a particular symbol next
        virtual double prob(bit_t b);

        // process a new piece of sensory experience
        virtual void update(bit_t b);

    private:

        // update the switching posterior weights 
        void updatePosteriors(SkipNode &n, int n_submodels, double alpha, 
            double log_alpha, double log_stop_mul, double log_split_mul);

        // compute the switching rate for a given time t
        double switchRate(size_t t) const;

        // gets the node's index into the hash table
        SkipNode &getNode(zobhash_t hash, int depth, int submodels) const;

        // Returns how many submodels a node at a given context position, with a number of
        // skips left, has.
        int numSubmodels(int position, int skips) const;

        // precomputes the zobrist hash keys
        static void initZobrist();

        // update the weighted estimates
        void updateWeighted(SkipNode &n, zobhash_t hash, int end_idx, int depth);

        // performs the update operation to maintain a switching posterior
        void posteriorUpdate(SkipNode &n, double log_scale, double log_alpha, 
            double log_K, double log_mul, double &log_post) const;

        // lazy allocation of skipping posterior weights
        void lazyAllocate(SkipNode &n, int n_submodels, int skips_left) const;

        // compute the information needed to update the current context stats
        void getContextInfo(zobhash_t &hash, const indices_t &idxs) const;

        // functions for precomputing the context indices
        void initContexts();
        void genBoundedSkipContexts(std::string buf, int depth, int skips);
        void printContexts(skip_contexts_t &l) const;

        // miscellaneous precomputations
        void initLogTbl();
        void initHashDeltas();
        void initAuxInfo();
        void initIndices();

        // compute the current binary context
        void getContext();

        // static members
        static zobhash_t s_zobtbl[MaxDepth][2];
        static weight_t s_log_tbl[LogTblSize];
        
        std::vector<aux_info_list_t> m_auxinfo;
        history_t &m_history;
        int m_depth;
        int m_skips;
        skip_contexts_t m_skip_contexts;
        std::vector<indices_list_t> m_indices;
        std::vector<size_t> m_context;
        std::vector<weight_t> m_log_skip_preds;
        SkipNode *m_nodes;
        zobhash_t m_mask;
};


#endif // __SKIPCTS_HPP__






