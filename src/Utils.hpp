/*
 * Authors: Gustavo V. Barroso
 * Created: 08/12/2022
 * Last modified: 08/12/2022
 *
 */


#ifndef _UTILS_H_
#define _UTILS_H_

#include <iostream>
#include <cmath>
#include <cstring>
#include <vector>
#include <algorithm>
#include <cstdlib>
#include <utility>

#include <Bpp/Numeric/Constraints.h>
#include <Bpp/Numeric/ParameterList.h>

#include "Model.hpp"
#include "OptionsContainer.hpp"

class Utils
{

private:

public:
  Utils():
  { }

public:
  std::vector<bpp::ParameterList> fetchParamsTable(const std::string& fileName);

};

#endif
