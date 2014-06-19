/*  Lookup table generator for the ICSI logarithm function
    Version 0.6 beta
    Build date: November 13th, 2007

    Copyright (C) 2007 International Computer Science Institute
    1947 Center Street. Suite 600
    Berkeley, CA 94704

    Contact information:
    Oriol Vinyals   vinyals@icsi.berkeley.edu
    Gerald Friedland    fractor@icsi.berkeley.edu

    Acknowledgements:
    Thanks to Harrison Ainsworth (hxa7241@gmail.com) for his idea that
    doubled the accuracy.

    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License along
    with this program; if not, write to the Free Software Foundation, Inc.,
    51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
*/

#include "icsilog.h"

#include <math.h>


double log2(double x) {
    return log(x) / log(2.0);
}


/*
This method fills a given array of floats with the information necessary to compute the icsi_log. This method has to be called before any call to icsi_log.
Parameters:
    n is the number of bits used from the mantissa (0<=n<=23). Higher n means higher accuracy but slower execution. We found that a good value for n is 14.
    lookup_table requires a float* pointing to a continuous (preallocated) memory array of 2^n*sizeof(float) bytes.
Return values: void
*/
void fill_icsi_log_table(const int n, float *lookup_table)
{
    float numlog;
    int incr,i,p;
    int *const exp_ptr = ((int*)&numlog);
    int x = *exp_ptr; /*x is the float treated as an integer*/
    x = 0x3F800000; /*set the exponent to 0 so numlog=1.0*/
        *exp_ptr = x;
    incr = 1 << (23-n); /*amount to increase the mantissa*/
    p = 1 << n;
    for(i=0;i<p;i++)
    {
        lookup_table[i] = (float) log2(numlog); /*save the log of the value*/
        x += incr;
        *exp_ptr = x; /*update the float value*/
    }
}


/* ICSIlog V 2.0 */
void fill_icsi_log_table2(const unsigned precision, float* const   pTable)
{
    /* step along table elements and x-axis positions
      (start with extra half increment, so the steps intersect at their midpoints.) */
    float oneToTwo = 1.0f + (1.0f / (float)( 1 <<(precision + 1) ));
    int i;
    for(i = 0;  i < (1 << precision);  ++i )
    {
        // make y-axis value for table element
        pTable[i] = logf(oneToTwo) / 0.69314718055995f;

        oneToTwo += 1.0f / (float)( 1 << precision );
    }
}
