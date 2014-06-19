#ifndef __FACMODELS_HPP__
#define __FACMODELS_HPP__

/******************************
      Author: Joel Veness
        Date: 2011
******************************/

#include "common.hpp"
#include "ctw.hpp"
#include "cts.hpp"
#include "skipcts.hpp"
#include "factor.hpp"

typedef Factor<ContextTree,8>   FactoredContextTree;
typedef Factor<SwitchingTree,8> FactoredSwitchingTree;
typedef Factor<SkipCTS, 8>      FactoredSkipCTS;

#endif // __FACMODELS_HPP__


