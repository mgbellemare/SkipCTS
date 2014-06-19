/******************************
      Author: Joel Veness
        Date: 2011
*******************************/

#include "ac.hpp"

#include <limits>
#include <iostream>
#include <fstream>
#include <cassert>
#include <stdexcept>


// arithmetic encoder/decoder constants
static const size_t CodeValueBits = 32;
static const uint64_t TopValue    = (uint64_t(1) << CodeValueBits) - 1;
static const uint64_t FirstQtr    = TopValue / 4 + 1;
static const uint64_t Half        = 2 * FirstQtr;
static const uint64_t ThirdQtr    = 3 * FirstQtr;

// precomputed values
static const double TopValueDbl = double(TopValue);


/* force probabilities recieved from the statistical modeller to be far enough
   from 0 and 1 so that our arithmetic encoder can use them */
static inline void clipProbabilities(double &p) {

    p = std::max(p, 0.0001);
    p = std::min(p, 0.9999);
}


/* create the arithmetic encoder */
BinaryArithmeticEncoder::BinaryArithmeticEncoder(const std::string &filename) {

    m_ofs.open(filename.c_str(), std::ios::out | std::ios::binary);
    if (!m_ofs) {
        std::string msg = "BinaryArithmeticEncoder cannot open ";
        msg += filename;
        throw std::runtime_error(msg);
    }

    // reserve space for header (just the file size in bits)
    char buf[sizeof(m_size)];
    m_ofs.write(buf, sizeof(m_size));

    m_low   = 0;
    m_high  = TopValue;
    m_fbits = 0;

    m_upto    = 0;
    m_buffer  = 0;
    m_size    = 0;
}


/* write the next bit, being careful to record previous middle subdivisions */
void BinaryArithmeticEncoder::bit_plus_follow(int bit) {

    write(bit);
    while (m_fbits > 0)  {
        write(!bit);
        m_fbits -= 1;
    }
}


/* send a bit to the encoder */
void BinaryArithmeticEncoder::encode(int bit, double pr_one) {

    assert(m_ofs.good());

    clipProbabilities(pr_one);

    double pr_zero = 1.0 - pr_one;
    uint64_t range = (m_high - m_low) + 1;
    uint64_t split = m_low + range * uint64_t(pr_zero * TopValueDbl) / TopValue;

    // narrow range
    if (bit == 1) {
        m_low  = split;
    } else {
        m_high = split - 1;
    }

    for (;;)  {

        if (m_high < Half)  {
            bit_plus_follow(0);
        }  else if (m_low >= Half)  {
            bit_plus_follow(1);
            m_low -= Half;
            m_high -= Half;
        }  else if (m_low >= FirstQtr && m_high < ThirdQtr)  {
            m_fbits += 1;
            m_low  -= FirstQtr;
            m_high -= FirstQtr;
        }  else break;

        m_low  = 2 * m_low;
        m_high = 2 * m_high + 1;
    }

    m_size++;
}


/* sends a bit to the stream */
void BinaryArithmeticEncoder::write(uint32_t bit) {

    m_buffer |= char(bit) << m_upto;
    m_upto++;

    if (m_upto == 8) flush();
}


/* flushes the current write buffer */
void BinaryArithmeticEncoder::flush() {

    m_ofs.put(m_buffer);
    m_buffer = 0;
    m_upto = 0;
}


/* finalises the arithmetic encoded file */
void BinaryArithmeticEncoder::close() {

    m_fbits += 1;
    if (m_low < FirstQtr) bit_plus_follow(0);
    else bit_plus_follow(1);

    // flush the buffer
    if (m_upto > 0) flush();

    // write the final file size to the header
    m_ofs.seekp(0, std::ios_base::beg);
    m_ofs.write((const char *) &m_size, sizeof(m_size));
    m_ofs.close();
}


/* opens a new decoder that can read bits from a file */
BinaryArithmeticDecoder::BinaryArithmeticDecoder(const std::string &filename) {

    m_ifs.open(filename.c_str(), std::ios::in | std::ios::binary);
    if (!m_ifs) {
        std::string msg = "BinaryArithmeticDecoder cannot open ";
        msg += filename;
        throw std::runtime_error(msg);
    }

    // determine the size of the file
    m_ifs.read((char *) &m_size, sizeof(m_size));

    m_bits_left = 0;
    m_value = 0;
    for (size_t i=1; i <= CodeValueBits; ++i) {
        m_value = 2 * m_value + readBit();
    }
    m_low = 0;
    m_high = TopValue;
}


/* read a new bit */
int BinaryArithmeticDecoder::decode(double pr_one) {

    clipProbabilities(pr_one);

    double pr_zero = 1.0 - pr_one;
    uint64_t range = (m_high - m_low) + 1;
    uint64_t split = m_low + range * uint64_t(pr_zero * TopValueDbl) / TopValue;

    // determine bit
    int bit = m_value < split ? 0 : 1;

    // narrow range
    if (bit == 1) {
        m_low  = split;
    } else {
        m_high = split - 1;
    }

    // rescale interval
    for (;;)  {

        if (m_high < Half)  {
            // do nothing
        }  else if (m_low >= Half)  {
            m_value -= Half;
            m_low   -= Half;
            m_high  -= Half;
        }  else if (m_low >= FirstQtr && m_high < ThirdQtr)  {
            m_value -= FirstQtr;
            m_low   -= FirstQtr;
            m_high  -= FirstQtr;
        } else break;

        m_low   = 2 * m_low;
        m_high  = 2 * m_high + 1;
        m_value = 2 * m_value + readBit();
    }

    return bit;
}


/* reads another bit from the file buffer */
uint32_t BinaryArithmeticDecoder::readBit() {

    // load some more data if needed
    if (m_bits_left == 0) {

        // if no data is left, pad with zeros
        if (!m_ifs.good()) return 0;

        m_buffer = m_ifs.get();
        m_bits_left = 8;
    }

    uint32_t rval = m_buffer & 1;
    m_buffer >>= 1;
    m_bits_left--;

    return rval;
}


/* finalise the decoding process */
void BinaryArithmeticDecoder::close() {

    m_ifs.close();
}



