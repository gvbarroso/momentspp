/*
 * Authors: Gustavo V. Barroso
 * Created: 29/07/2022
 * Last modified: 15/09/2022
 *
 */


#ifndef _OPTIMIZATIONWRAPPER_H_
#define _OPTIMIZATIONWRAPPER_H_

#include <string>
#include <vector>
#include <limits>

#include <Bpp/Numeric/Function/ReparametrizationFunctionWrapper.h>
#include <Bpp/Numeric/Function/PowellMultiDimensions.h>
#include <Bpp/Numeric/Function/BfgsMultiDimensions.h>
#include <Bpp/Numeric/Function/ThreePointsNumericalDerivative.h>
#include <Bpp/Numeric/Matrix/MatrixTools.h>
//#include <Bpp/Numeric/Random/RandomTools.h>
#include <Bpp/Numeric/ParameterList.h>
#include <Bpp/App/BppApplication.h>
#include <Bpp/App/ApplicationTools.h>
#include <Bpp/Phyl/PseudoNewtonOptimizer.h>

#include "OptionsContainer.hpp"
#include "Model.hpp"


class OptimizationWrapper
{
    
private:
  // population index -> index of 2 parental populations in immediately previous epoch (can be the same, eg (3, (1,1)))
  std::vector<std::map<size_t, std::pair<size_t, size_t>>> popMaps_;
  OptionsContainer options_;
  
public:
  OptimizationWrapper(const OptionsContainer& options):
  popMaps_(),
  options_(options)
  { }

  OptimizationWrapper(const std::map<std::string, std::string>& options):
  popMaps_(),
  options_(options)
  { }
  
public:
  void optimize(); // optimze from scratch

  void parsePopsFile(const std::string& name);

private:
  void fitModel_(Model* model);

  void writeEstimatesToFile_(Model* model);
  
};

#endif
