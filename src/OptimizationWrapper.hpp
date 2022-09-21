/*
 * Authors: Gustavo V. Barroso
 * Created: 29/07/2022
 * Last modified: 21/09/2022
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

#include "PolymorphismData.hpp"
#include "OptionsContainer.hpp"
#include "Model.hpp"


class OptimizationWrapper
{
    
private:
  std::vector<std::vector<size_t>> popList_; // population indices per epoch
  OptionsContainer options_;
  
public:
  OptimizationWrapper(const OptionsContainer& options):
  popList_(),
  options_(options)
  { }

  OptimizationWrapper(const std::map<std::string, std::string>& options):
  popList_(),
  options_(options)
  { }
  
public:
  void optimize(const PolymorphismData& data);

  void parsePopsFile(const std::string& name);

private:
  void fitModel_(Model* model);

  void writeEstimatesToFile_(Model* model);
  
};

#endif
