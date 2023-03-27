/*
 * Authors: Gustavo V. Barroso
 * Created: 10/08/2022
 * Last modified: 27/03/2022
 *
 */


#ifndef _MIGRATION_H_
#define _MIGRATION_H_

#include "AbstractOperator.hpp"

class Migration: public AbstractOperator
{

private:
  Eigen::MatrixXd littleMigMat_; // P x P

public:
  Migration(const bpp::ParameterList migParams, const SumStatsLibrary& sslib):
  AbstractOperator(),
  littleMigMat_()
  {
    includeParameters_(migParams);
    prevParams_.addParameters(getParameters()); // inits list of "previous" parameters
    setLittleMat_();
    //testStationary_();
    setUpMatrices_(sslib);
  }

  Migration(const std::vector<double>& initValues, std::shared_ptr<bpp::IntervalConstraint> ic, const SumStatsLibrary& sslib):
  AbstractOperator(),
  littleMigMat_()
  {
    size_t idx = 0;
    for(auto itI = std::begin(sslib.getPopIndices()); itI != std::end(sslib.getPopIndices()); ++itI) // for each population modeled in epoch i
    {
      for(auto itJ = std::begin(sslib.getPopIndices()); itJ != std::end(sslib.getPopIndices()); ++itJ)
      {
        if((*itI) != (*itJ)) // if population indices are different
        {
          //std::shared_ptr<bpp::Parameter> param = std::make_shared<bpp::Parameter>("m_" + bpp::TextTools::toString((*itI)) + bpp::TextTools::toString((*itJ)), initValue, ic);
          //addParameter_(param.get());
          addParameter_(new bpp::Parameter("m_" + bpp::TextTools::toString((*itI)) + "_" + bpp::TextTools::toString((*itJ)), initValues[idx], ic));
          ++idx;
        }
      }
    }

    prevParams_.addParameters(getParameters()); // inits list of "previous" parameters
    setLittleMat_();
    //testStationary_();
    setUpMatrices_(sslib);
  }

  virtual Migration* clone() const override
  {
    return new Migration(*this);
  }

  const Eigen::MatrixXd& getLittleMigMat()
  {
    return littleMigMat_;
  }

private:
  void setUpMatrices_(const SumStatsLibrary& sslib) override;

  void updateMatrices_() override;

  void setLittleMat_();

  //void testStationary_();

  size_t fetchNumPops_();

};

#endif
