#ifndef __FACTOR_HPP__
#define __FACTOR_HPP__

/******************************
      Author: Joel Veness
        Date: 2011
******************************/

#include "common.hpp"

// a generic structure for combining compressors for byte oriented data
template <typename T, size_t N>
class Factor : public Compressor {

    public:

        /// create a context tree of specified maximum depth and size
        Factor(history_t &history, size_t depth);

        /// delete the factored context tree
        ~Factor();

        /// file extension
        const char *fileExtension() const;

        /// the logarithm of the probability of all processed experience
        double logBlockProbability() const;

        // the probability of seeing a particular symbol next
        double prob(bit_t b);

        /// process a new piece of sensory experience
        void update(bit_t b);

        /// the depth of the context tree
        size_t depth() const;

        /// number of nodes in the context tree
        size_t size() const;

    private:

        /// copy contructor / assignment operator disabled
        Factor(const Factor &rhs);
        const Factor &operator=(const Factor &rhs);

        T *m_models[N];
        size_t m_depth;
        history_t &m_history;
};


/* create the factored context tree */
template <typename T, size_t N>
Factor<T,N>::Factor(history_t &history, size_t depth) :
    m_depth(depth),
    m_history(history)
{
    for (size_t i=0; i < N; i++) {
        //m_models[i] = new T(history, depth+i, static_cast<int>(i));
        m_models[i] = new T(history, depth+i);
    }
}


/* delete the factored context tree */
template <typename T, size_t N>
Factor<T,N>::~Factor() {

    for (size_t i=0; i < N; i++) {
        delete m_models[i];
    }
}


/* the logarithm of the probability of all processed experience */
template <typename T, size_t N>
double Factor<T,N>::logBlockProbability() const {

    double sum = 0.0;
    for (size_t i=0; i < N; i++) {
        sum += m_models[i]->logBlockProbability();
    }
    return sum;
}


/* the probability of seeing a particular symbol next */
template <typename T, size_t N>
double Factor<T,N>::prob(bit_t b) {

    size_t idx = (m_history.size() - m_depth) % N;
    return m_models[idx]->prob(b);
}


/* process a new piece of data */
template <typename T, size_t N>
void Factor<T,N>::update(bit_t b) {

    size_t idx = (m_history.size() - m_depth) % N;
    m_models[idx]->update(b);
}


/* the depth of the context tree */
template <typename T, size_t N>
size_t Factor<T,N>::depth() const {

    return m_depth;
}


/* number of nodes in the factored context tree */
template <typename T, size_t N>
size_t Factor<T,N>::size() const {

    size_t sum = 0;
    for (size_t i=0; i < N; i++) {
        sum += m_models[i]->size();
    }
    return sum;
}


/* file extension */
template <typename T, size_t N>
const char *Factor<T, N>::fileExtension() const {

    static const std::string ext = std::string("fac") + m_models[0]->fileExtension();
    return ext.c_str();
}


#endif // __FACTOR_HPP__

