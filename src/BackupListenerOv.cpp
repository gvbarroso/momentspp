/*
 * Authors: Gustavo V. Barroso and Julien Y. Dutheil
 * Created: 24/08/2022
 * Last modified: 24/08/2022
 *
 */


#include "BackupListenerOv.hpp"

//prints list of current best parameter values (in original space) to backup file
void BackupListenerOv::optimizationStepPerformed(const bpp::OptimizationEvent& event)
{
  const bpp::Function* reparamFun = event.getOptimizer()->getFunction();
  bpp::ParameterList params = dynamic_cast<const bpp::ReparametrizationFunctionWrapper*>(reparamFun)->getFunction().getParameters();

  std::ofstream file(backupFile_.c_str(), std::ios::out);
  double aic = 2. * params.size() + 2. * event.getOptimizer()->getFunction()->getValue();
  file << "AIC = " << aic << std::endl << std::endl;

  for(size_t i = 0; i < params.size(); ++i)
    file << params[i].getName() << "\t" << params[i].getValue() << std::endl;

  file.close();
}
