#ifndef __CTS_HPP__
#define __CTS_HPP__

/****************************************************************
      Author: Joel Veness
        Date: 2011
        Info: A Context Tree Switching (CTS) implementation.
*****************************************************************/

#include "common.hpp"

// boost includes
#include <boost/utility.hpp>
#include <boost/pool/pool.hpp>


// context tree node
class SNode {

    friend class SwitchingTree;
    friend std::ostream &operator<<(std::ostream &o, const SNode &);

    public:

        SNode(int depth);
        explicit SNode(const SNode &rhs, int pindx);

        /// process a new binary symbol, with switching rate alpha, and blend 1-2*alpha
        void update(bit_t b, double log_alpha, double log_blend, double log_split_mul);

        /// log weighted blocked probability
        weight_t logProbWeighted() const;

        /// log KT estimated probability
        weight_t logProbEstimated() const;

        /// child corresponding to a particular symbol
        const SNode *child(bit_t b) const;

        /// the number of times this context been visited
        count_t visits() const;

        /// number of descendents
        size_t size() const;

    private:

        // compute the result of an update call non-destructively
        double updateNonDestructive(bit_t b, double c_weighted, double c_weighted_old) const;

        // is the current node a leaf node?
        bool isLeaf() const;

        // compute the logarithm of the KT-estimator update multiplier
        double logKTMul(bit_t b) const;

        weight_t m_log_prob_est;
        weight_t m_log_prob_weighted;

        // switching weights
        weight_t m_log_b;
        weight_t m_log_s;

        // one slot for each binary value
        count_t m_count[2];
        SNode *m_child[2];

        // m_backidx >= 0 when a previous symbol's context was pruned
        int m_pruned_idx;
};

extern std::ostream &operator<<(std::ostream &o, const SNode &sn);


// a context tree used for CTW mixing
class SwitchingTree : public Compressor, boost::noncopyable {

    typedef std::pair<SNode *, SNode> ctpair_t;

    public:

        /// create a context tree of specified maximum depth and size
        SwitchingTree(history_t &history, size_t depth, int phase=-1);

        /// delete the context tree
        ~SwitchingTree();

        /// file extension
        const char *fileExtension() const { return "cts"; }

        /// the logarithm of the probability of all processed experience
        double logBlockProbability() const;

        /// the probability of seeing a particular symbol next
        double prob(bit_t b);

        /// process a new piece of sensory experience
        void update(bit_t b);

        /// the depth of the context tree
        size_t depth() const;

        /// number of nodes in the context tree
        size_t size() const;

    private:

        // compute the switching rate for a given time t
        double switchRate(size_t t) const;

        // recover the memory used by a node
        void reclaimMemory(SNode *n);

        // computes the context, creates relevant nodes and determine the path to update
        void makeContextAndPath();

        // determine whether two contexts are identical
        static bool contextsEqual(const context_t &lhs, const context_t &rhs);

        // compute the current context
        void getContext(const history_t &h, context_t &context, int idx = -1);

        // create (if necessary) all of the nodes in the current context
        void createNodesInCurrentContext(const context_t &context);

        // recursively deletes the nodes in a context tree
        void deleteCT(SNode *root);

        boost::pool<> m_ctnode_pool;

        SNode *m_root;
        int m_phase;
        size_t m_depth;
        context_t m_context;
        context_t m_pcontext;
        std::vector<SNode *> m_created;
        std::vector<ctpair_t> m_modified;
        std::vector<SNode *> m_path;
        std::vector<weight_t> m_log_old_weights;
        history_t &m_history;
        int m_prob_cache;
};

#endif // __CTS_HPP__

