******************************* 
  Skip Context Tree Switching
*******************************


Authors: Marc G. Bellemare, Joel Veness, Erik Talvitie
Date: 17/06/2014

Introduction
------------

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


Installation Requirements
-------------------------

This implementation of SkipCTS requires
    * Boost 1.54.0
    * CMake

A typical installation on a Unix system proceeds as follows: 

> mkdir build
> cd build
> cmake ../src
> make -j 4
> cd ..
> ln -s build/skipcts .

It is possible to specify the Boost library paths for CMake, by invoking cmake as

> export BOOST_ROOT=<path-to-boost>
> cmake -DBOOST_NO_SYSTEM_PATHS=1 ../src


SkipCTS may then be run as, e.g.,

> ./skipcts --method skipcts --file <file> --depth 48 --slots 26 --skips 1 


Bugs, Caveats, Features
-----------------------

Speed. The per-step cost of SkipCTS is O(D^{2K+1}). As such, bit-level SkipCTS can run slowly. 
Restricting where the skips can occur is one way to speed things up. Intuitively, we should get 
most of the benefits of skipping by considering skips at byte boundaries only. This extension 
should be fairly straightforward -- let us know if you try it!

SkipCTS::prob(). The prob() method currently assumes that a NULL child node implies a NULL 
subtree. This is not strictly true when hashing is involved, and could in theory result in 
incorrect compression/decompression. One possible fix would be to create the queried nodes, then 
delete them after the call to prob(). 

FacSkipCTS. The bit-factored SkipCTS is not currently fully implemented.

Hashing. The hashing scheme should never hash together nodes which split over a different
number of symbols (the 'depth' field), or have a different number of submodels. This could
result in greater-than-one probabilities. Hence the convolutions of the getNode() method. 




**************************  
  Context Tree Switching
**************************

Author:    Joel Veness
Date:      14/11/2011

Introduction:
-------------

This is an implementation of the Context Tree Switching algorithm, 
described in the accompanying technical report at 
    http://arxiv.org/abs/1111.3182v1. 

It is an extension of the Context Tree Weighting algorithm, and works by replacing
the recursive weighting step by a computationally efficienct switching 
technique. The probabilities estimated by this technique then drive a 
standard binary arithmetic encoder, which produces the compressed file.

This is a proof-of-concept implementation, and is currently quite memory hungry. 
However, provided you have a modern machine with 2 gigs or more of RAM, the default
settings should work fine for files less than 10meg.

Program Usage:
---------------

Run cts.exe --help for program options. It should display something like:

Usage:
  --help
  --depth arg (=48)      Maximum depth to look back in bits. Higher values use more RAM.
  --method arg (=faccts) Compression method. (ctw, cts, facctw, faccts)
  --file arg             File to compress/decompress. A compressed file with the same name plus an
                         extension (e.g. foo.cts) will be produced. Decompression is chosen
                         automatically if the file extension is .ctw, .cts, .facctw, or .faccts.

A sample usage might go like this:

> cts.exe --depth 128 --file foo

Later, to recover foo, run

> cts.exe --depth 128 --file foo.faccts

About the different methods:
----------------------------

ctw:	A bare-bones implementation of vanilla CTW, with no enhancements.
cts:    An implementation of cts. Treats everything as a raw binary stream.
facctw: CTW plus a binary decomposition to exploit byte oriented data.
faccts: CTS plus a binary decomposition to exploit byte oriented data

There are some optional flags that can be set in the source code.
Much better CTW performance can be obtained by enabling the zero-redundancy 
estimator and disabling weighting at byte boundaries. For the best CTW
implementation I am aware of, see http://www.ele.tue.nl/ctw/.

The most important parameter in CTS seems to be the choice of prior weight to assign
to a switch (controlled by the log_switch_prior and log_kt_prior parameters in the code). 
The analysis and results reported in the technical report for vanilla CTS use a value of 0.5. 
Much better across the board performance appears when log_switch_prior is set to 0.925.
The code here uses this better performing constant as the default.

Additional Credits:
--------------------

In addition to my co-authors on the attached technical report, I'd like to mention:

   - Frans Willems, Yuri Shtarkov, Tjalling Tjalkens, and Paul Volf for their fantastic papers!
     Their work got me interested in data compression, and in particular techniques that mix theory and practice.
     
   - The CTW 0.1 website (http://www.ele.tue.nl/ctw/) contained lots of helpful ideas.
     In particular, the "Strict Unique Path Pruning" idea can luckily be adapted to CTS,
     which kept the memory requirements to a much more reasonable level.
     
Musings/Todo:
--------------
   
   - The settings used to compress and decompress need to match. If they don't, all hell will break loose.
     In the future this should be done properly by storing appropriate information in the file header.
     
   - The memory requirements can be reduced considerably by packing the information
     required at each node. This would be good to do in the future.
     
   - This implementation could be improved by converting the logic to use only integer arithmetic.
   
   - There are probably many ways to speed this implementation up considerably. 
     I'd welcome some suggestions, since this is my first data compression program.
     
   - It's quite likely that compression performance can be improved by tuning some parameters inside cts.cpp. 
     Changing the global log_kt_prior constant to a depth dependent constant should be investigated.
   
   - For maximum compression, disable the fast math routines in cts.cpp. By default these settings
     are enabled since the gain in speed is large. This more than offsets the loss in compression,
     which I found to be within a margin of ~= 0.01 bits per character on the Calgary Corpus.
     
     
