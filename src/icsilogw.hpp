#ifndef __ICSILOGW_HPP__
#define __ICSILOGW_HPP__

/******************************
      Author: Joel Veness
        Date: 2011
******************************/

#include "icsilog.h"

#include <cassert>

// a simple C++ wrapper for ICSILog
class ICSILog {

    public:

        ICSILog(unsigned int n);
        ~ICSILog() { free(m_table); }

        // fast approximation to std::log
        float log(float x) { return icsi_log(x, m_table, m_n); }

    private:

        // disable copy constructor and assignment operator
        ICSILog(const ICSILog &rhs);
        ICSILog &operator=(const ICSILog &rhs);

        float *m_table;
        int m_n;
};


/* precompute the ICSILog tables */
inline ICSILog::ICSILog(unsigned int n) {

    assert(n <= 64);
    m_n = int(n);
    m_table = (float*) malloc(int(1 << n) * sizeof(float));
    fill_icsi_log_table(n, m_table);
}


#endif // __ICSILOGW_HPP__

