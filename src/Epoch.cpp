/*
 * Authors: Gustavo V. Barroso
 * Created: 31/08/2022
 * Last modified: 05/12/2022
 *
 */

#include <ios>

#include "Log.hpp"
#include "Epoch.hpp"

void Epoch::fireParameterChanged(const bpp::ParameterList& params)
{
  if(matchParametersValues(params))
    updateOperators_(params);

  Eigen::SparseMatrix<double> mat = operators_[0]->fetchCombinedMatrix(); // init mat

  for(size_t i = 1; i < operators_.size(); ++i)
    mat = mat * operators_[i]->fetchCombinedMatrix();

  transitionMatrix_ = mat; // converts to dense format
}

// this method is where the heavier Eigen3 linear algebra takes place
void Epoch::computeExpectedSumStats(Eigen::VectorXd& y)
{
  Log timer;
  timer.openFile("timing.txt");

  Eigen::EigenSolver<Eigen::MatrixXd> es_(transitionMatrix_);

  timer.start_timer();
  auto tmp1 = transitionMatrix_.pow(duration()) * y;
  timer.stop_timer();

  timer.start_timer();
  auto tmp2 = es_.eigenvectors().real() * es_.eigenvalues().real().array().pow(duration()).matrix().asDiagonal() * es_.eigenvectors().real().inverse() * y;
  timer.stop_timer();

  timer.start_timer();
  y = es_.pseudoEigenvectors() * es_.pseudoEigenvalueMatrix().pow(duration()) * es_.pseudoEigenvectors().inverse() * y;
  timer.stop_timer();
}

void Epoch::transferStatistics(Eigen::VectorXd& y)
{
  Eigen::VectorXd tmp(ssl_.getMoments().size()); // y and tmp have potentially different sizes
  tmp.setZero();

  // for each Moment in *this epoch, we assign its value from its parental Moment from the previous epoch (from which y comes)
  for(int i = 0; i < tmp.size(); ++i)
  {
    int idx = ssl_.getMoments()[i].getParent()->getPosition();
    tmp(i) = y(idx);
  }

  y = tmp; // swap
}

void Epoch::updateMoments(const Eigen::VectorXd& y)
{
  if(y.size() != ssl_.getMoments().size())
    throw bpp::Exception("Epoch::attempted to update moments from vector of different size!");

  for(int i = 0; i < y.size(); ++i)
    ssl_.getMoments()[i].setValue(y(i));
}

void Epoch::computeSteadyState_()
{
  #ifdef VERBOSE
  Log logger;
  logger.openFile("matrices.txt");
  Eigen::SparseMatrix<double> test(ssl_.getNumStats(), ssl_.getNumStats());
  test.setIdentity();

  for(size_t i = 0; i < operators_.size(); ++i)
  {
    Eigen::SparseMatrix<double> tmp(ssl_.getNumStats(), ssl_.getNumStats());
    tmp.setZero();
    for(size_t j = 0; j < operators_[i]->getMatrices().size(); ++j)
    {
      tmp += operators_[i]->getMatrices()[j];
      bpp::ParameterList pl;
      pl.addParameter(operators_[i]->getParameters()[j]);
      pl.printParameters(logger.getLogFile());
      logger.getLogFile() << std::setprecision(1e-12) << std::scientific << operators_[i]->getMatrices()[j] << "\n";
    }

    operators_[i]->getParameters().printParameters(logger.getLogFile());
    logger.getLogFile() << "\n\nsum of entries = " << std::setprecision(1e-12) << std::scientific << tmp.sum() << "\n";
    logger.getLogFile() << std::scientific << tmp << "\n";

    logger.getLogFile() << "accumulated transition matrix:\n";
    test = operators_[i]->fetchCombinedMatrix() * test;
    logger.getLogFile() << std::setprecision(1e-12) << std::scientific << test << "\n";
  }
  #endif

  Eigen::SparseMatrix<double> mat = operators_[0]->fetchCombinedMatrix(); // init mat

  for(size_t i = 1; i < operators_.size(); ++i)
    mat = operators_[i]->fetchCombinedMatrix() * mat;

  transitionMatrix_ = mat; // converts to dense format

  // we find the eigenvector associated with (leading) eigenvalue == 1 in transitionMatrix_
  Eigen::EigenSolver<Eigen::MatrixXd> es(transitionMatrix_);

  int idx = 0;
  for(int i = 0; i < es.eigenvalues().size(); ++i)
  {
    // finding the maximum value (should be == 1., but searching for equality is problematic due to precision issues)
    if(es.eigenvalues().real()(i) > es.eigenvalues().real()(idx))
      idx = i;
  }

  steadYstate_ = es.eigenvectors().col(idx).real();
  steadYstate_ /= steadYstate_(ssl_.getDummyIndex()); // divide by I moment, which embodies constant used for Eigen decomposition

  updateMoments(steadYstate_);
}
