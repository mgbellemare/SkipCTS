#ifndef __AC_HPP__
#define __AC_HPP__

/***************************************************************
      Author: Joel Veness
        Date: 2011
        Info: A binary arithmetic decoder/encoder.

        Thanks to the Arithmetic Coding Library by
        Fred Wheeler (2000) for inspiration, and also to
        the accompanying source code provided with the article
        "Arithmetic Coding for Data Compression", by
        Witten et al.
*****************************************************************/

#include <fstream>
#include <iostream>

#include <boost/cstdint.hpp>


// fixed size types
typedef boost::uint32_t uint32_t;
typedef boost::uint64_t uint64_t;
typedef boost::int8_t int8_t;


// binary arithemtic encoder
class BinaryArithmeticEncoder {

    public:

        // opens a new encoder that outputs to filename
        BinaryArithmeticEncoder(const std::string &filename);

        // send a bit to the encoder and the probability of it being a 1
        void encode(int bit, double pr_one);

        // finalise the encoding process
        void close();

    private:

        // write the next bit, being careful to record previous middle subdivisions
        void bit_plus_follow(int bit);

        // writes a bit to the internal buffer
        void write(uint32_t bit);

        // flushes the internal buffer to file
        void flush();

        std::ofstream m_ofs;

        uint64_t m_low;
        uint64_t m_high;
        uint64_t m_fbits;

        char m_buffer;
        char m_upto;
        uint32_t m_size;
};


// binary arithmetic decoder
class BinaryArithmeticDecoder {

    public:

        // opens a new decoder that can read bits from a file
        BinaryArithmeticDecoder(const std::string &filename);

        // read a new bit
        int decode(double pr_one);

        // size of decompressed file
        uint32_t size() const { return m_size; }

        // finalise the decoding process
        void close();

    private:

        uint32_t readBit();

        std::ifstream m_ifs;

        uint64_t m_low;
        uint64_t m_high;
        uint64_t m_value;

        uint32_t m_size;

        int8_t m_bits_left;
        int8_t m_buffer;
};


#endif // __AC_HPP__


