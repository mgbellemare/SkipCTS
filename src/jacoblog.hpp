#ifndef __JACOBLOG_HPP__
#define __JACOBLOG_HPP__

/******************************
      Author: Joel Veness
        Date: 2011
******************************/

#include <cassert>


/* builds a cache for a fast Jacobian Logarithm approximation */
class JacobianLogTable {

    public:

        // precompute cache
        JacobianLogTable(size_t entries);

        // release cache memory
        ~JacobianLogTable();

        // compute the Jacobian Logarithm log(1+exp(x))
        double jacobianLog(double x) const;

    private:

        double *m_tbl, *m_tbl2;
        double m_index_mul;
        double m_step;
};


/* create a table with a given number of entries */
JacobianLogTable::JacobianLogTable(size_t entries) :
    m_index_mul(0.01 * double(entries)),
    m_step(100.0 / double(entries))
{
    assert(entries >= 32);

    m_tbl  = new double[entries];
    m_tbl2 = new double[entries];
    assert(m_tbl != NULL && m_tbl2 != NULL);

    const double StepR = 1.0 / m_step;

    for (size_t i=0; i < entries; i++) {
        double x = m_step * i;
        m_tbl[i] = std::log(1.0 + std::exp(x)) - x;
    }

    // precompute (y1-y0) * StepR
    for (size_t i=0; i < entries-1; i++) {
        m_tbl2[i] = (m_tbl[i+1] - m_tbl[i]) * StepR;
    }
    m_tbl2[entries-1] = 0.0;
}


/* release the cache memory */
inline JacobianLogTable::~JacobianLogTable() {

    delete[] m_tbl;
    delete[] m_tbl2;
}


/* compute the Jacobian Logarithm with log(1+exp(x)) = x + corr */
inline double JacobianLogTable::jacobianLog(double x) const {

    if (x >= 100.0) return x;

    // use linear interpolation to construct an estimate
    int i = int(x * m_index_mul);
    double x0 = m_step * i;
    return x + m_tbl[i] + (x - x0) * m_tbl2[i];
}


#endif // __JACOBLOG_HPP__

