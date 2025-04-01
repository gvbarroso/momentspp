/*
 * Authors: Gustavo V. Barroso
 * Created: 29/07/2022
 * Last modified: 01/04/2025
 *
 */


#ifndef _OPTIMIZATIONWRAPPER_H_
#define _OPTIMIZATIONWRAPPER_H_

#include <map>
#include <string>
#include <vector>
#include <limits>

#include <Bpp/Numeric/Function/ReparametrizationFunctionWrapper.h>
#include <Bpp/Numeric/Function/PowellMultiDimensions.h>
#include <Bpp/Numeric/Function/BfgsMultiDimensions.h>
#include <Bpp/Numeric/Function/ThreePointsNumericalDerivative.h>
#include <Bpp/Numeric/Matrix/MatrixTools.h>
#include <Bpp/Numeric/ParameterList.h>
#include <Bpp/App/BppApplication.h>
#include <Bpp/App/ApplicationTools.h>

#include "Population.hpp"
#include "Data.hpp"
#include "Demes.hpp"
#include "OptionsContainer.hpp"
#include "Model.hpp"


class OptimizationWrapper
{
    
private:
  OptionsContainer options_;
  
public:
  OptimizationWrapper(const OptionsContainer& opt):
  options_(opt)
  { }
  
public:
  void fitModel(std::shared_ptr<Model> model);

private:
  void writeEstimatesToFile_(std::shared_ptr<Model> model);
  
};

#endif
