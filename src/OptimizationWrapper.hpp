/*
 * Authors: Gustavo V. Barroso
 * Created: 29/07/2022
 * Last modified: 29/08/2022
 *
 */


#ifndef _OPTIMIZATIONWRAPPER_H_
#define _OPTIMIZATIONWRAPPER_H_

#include <string>
#include <vector>
#include <limits>

#include <Bpp/Numeric/Function/ReparametrizationFunctionWrapper.h>
#include <Bpp/Numeric/Function/PowellMultiDimensions.h>
#include <Bpp/Numeric/Function/Bfgs.h>
#include <Bpp/Numeric/Function/ThreePointsNumericalDerivative.h>
#include <Bpp/Numeric/Matrix/MatrixTools.h>
#include <Bpp/Numeric/Random/RandomTools.h>
#include <Bpp/Numeric/ParameterList.h>
#include <Bpp/App/BppApplication.h>
#include <Bpp/App/ApplicationTools.h>
#include <Bpp/Phyl/PseudoNewtonOptimizer.h>

#include "Model.hpp"
#include "OptionsContainer.hpp"


class OptimizationWrapper
{
    
private:
  OptionsContainer options_;
  
public:
  OptimizationWrapper(const OptionsContainer& options):
  options_(options),
  { }

  OptimizationWrapper(const std::map<std::string, std::string>& options):
  options_(options),
  { }
  
public:
  void optimize();
  
  void optimize(const bpp::ParameterList& backupParams); // to resume optimisation after a crash:
  
private:
  void fitModel_(Model* model);

  void writeEstimatesToFile_(Model* model);
  
};

#endif
