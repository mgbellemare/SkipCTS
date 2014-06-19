SkipCTS
=======

Skip Context Tree Switching - Reference Implementation

Authors: Marc G. Bellemare, Joel Veness, Erik Talvitie

This is an implementation of the Skip Context Tree Switching (SkipCTS) algorithm based on
Joel Veness' original source for Context Tree Switching (CTS). SkipCTS is described in
the "Skip Context Tree Switching" by the authors, presented at ICML 2014. This paper can
be found at

http://www.github.com/mgbellemare/SkipCTS

along with the latest version of the present source code.

In this code, the implementations for CTS and SkipCTS differ in a number of ways:
    * The SkipCTS tree is stored within a hash table, whose size is specified with --slots
    * SkipCTS is additionally parametrized by K, the number of allowed skips (--skips) 
    * SkipCTS is also parametrized by a new set of skipping prior parameters (see skipcts.cpp)


This implementation of SkipCTS requires
    * Boost 1.54.0
    * CMake


See the included readme.txt for more details. 
