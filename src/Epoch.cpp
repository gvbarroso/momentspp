/*
 * Authors: Gustavo V. Barroso
 * Created: 31/08/2022
 * Last modified: 08/03/2023
 *
 */

#include <ios>

#include "Log.hpp"
#include "Epoch.hpp"

void Epoch::fireParameterChanged(const bpp::ParameterList& params)
{
  if(matchParametersValues(params))
    updateOperators_(params);

  Eigen::SparseMatrix<double> mat = operators_[0]->getTransitionMatrix(); // init mat

  for(size_t i = 1; i < operators_.size(); ++i)
    mat = mat * operators_[i]->getTransitionMatrix();

  transitionMatrix_ = mat; // converts to dense format
}

void Epoch::computeExpectedSumStats(Eigen::VectorXd& y)
{
  y = transitionMatrix_.pow(duration()) * y; // heavier Eigen3 linear algebra takes place
}

void Epoch::transferStatistics(Eigen::VectorXd& y) // y comes from previous Epoch
{
  Eigen::VectorXd tmp(ssl_.getMoments().size()); // y and tmp have potentially different sizes
  tmp.setZero();

  // for each Moment in *this Epoch, we assign its value from its parental Moment from the previous Epoch
  for(int i = 0; i < tmp.size(); ++i)
  {
    int idx = ssl_.getMoments()[i]->getParent()->getPosition();
    tmp(i) = y(idx);
  }

  y = tmp; // swap
}

void Epoch::updateMoments(const Eigen::VectorXd& y)
{
  assert(y.size() == static_cast<int>(ssl_.getMoments().size()));

  for(int i = 0; i < y.size(); ++i)
    ssl_.getMoments()[i]->setValue(y(i));
}

void Epoch::printRecursions(std::ostream& stream)
{
  for(size_t i = 0; i < ssl_.getMoments().size(); ++i)
  {
    if(ssl_.getMoments()[i]->getName() != "I")
    {
      int pos = static_cast<int>(ssl_.getMoments()[i]->getPosition()); // row in delta matrix
      stream << "delta_" << ssl_.getMoments()[i]->getName() << " = ";

      for(size_t j = 0; j < operators_.size(); ++j)
      {
        for(size_t k = 0; k < operators_[j]->getParameters().size(); ++k)
        {
          bpp::Parameter param = operators_[j]->getParameters()[k];
          std::string name = param.getName();
          auto mat = operators_[j]->getMatrix(k); // hard copy delta matrix
          mat = mat / param.getValue(); // convert back to coefficients

          for(int l = 0; l < mat.cols(); ++l)
          {
            if(mat.coeffRef(pos, l) != 0)
            {
              if(mat.coeffRef(pos, l) > 0)
                stream << "+";

              stream << bpp::TextTools::toString(mat.coeffRef(pos, l)) + "*" + name + "*" + ssl_.getMoments()[l]->getName() + " ";
            }
          }
        }
      }

      stream << "\n";
    }
  }
}

void Epoch::computeSteadyState_()
{
  #ifdef VERBOSE
  Log logger;
  logger.openFile(getName() + "_matrices.txt");
  Eigen::SparseMatrix<double> test(ssl_.getNumStats(), ssl_.getNumStats());
  test.setIdentity();

  for(size_t i = 0; i < operators_.size(); ++i)
  {
    Eigen::SparseMatrix<double> tmp(ssl_.getNumStats(), ssl_.getNumStats());
    tmp.setZero();
    for(size_t j = 0; j < operators_[i]->getMatrices().size(); ++j)
    {
      logger.getLogFile() << "\n\nsum of entries (delta matrix " << j << ") = " << std::setprecision(1e-12) << std::scientific << operators_[i]->getMatrices()[j].sum() << "\n";
      tmp += operators_[i]->getMatrices()[j];
      bpp::ParameterList pl;
      pl.addParameter(operators_[i]->getParameters()[j]);
      pl.printParameters(logger.getLogFile());
      logger.getLogFile() << std::setprecision(1e-12) << std::scientific << operators_[i]->getMatrices()[j] << "\n";
    }

    operators_[i]->getParameters().printParameters(logger.getLogFile());
    logger.getLogFile() << "\n\nsum of entries (operator combined delta matrix) = " << std::setprecision(1e-12) << std::scientific << tmp.sum() << "\n";
    logger.getLogFile() << std::setprecision(1e-12) << tmp << "\n";
    logger.getLogFile() << "operator transition matrix:\n";
    logger.getLogFile() << operators_[i]->getTransitionMatrix() << "\n\n";

    logger.getLogFile() << "accumulated transition matrix:\n";
    test = operators_[i]->getTransitionMatrix() * test;
    logger.getLogFile() << std::setprecision(1e-12) << test << "\n";
  }
  #endif

  Eigen::SparseMatrix<double> mat = operators_[0]->getTransitionMatrix(); // init mat

  for(size_t i = 1; i < operators_.size(); ++i)
    mat = operators_[i]->getTransitionMatrix() * mat;

  transitionMatrix_ = mat; // converts to dense format

  // we find the eigenvector associated with thr leading eigenvalue (== 1) in transitionMatrix_
  Eigen::EigenSolver<Eigen::MatrixXd> es(transitionMatrix_);

  int idx = 0;
  for(int i = 0; i < es.eigenvalues().size(); ++i)
  {
    // finding the maximum value (should be == 1., but searching for equality is problematic due to precision)
    if(es.eigenvalues().real()(i) > es.eigenvalues().real()(idx))
      idx = i;
  }

  steadYstate_ = es.eigenvectors().col(idx).real();
  steadYstate_ /= steadYstate_(ssl_.getDummyIndex()); // divide by I moment, which embodies scaling constant used for Eigen decomposition

  updateMoments(steadYstate_);
}
