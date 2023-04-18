/*
 * Authors: Gustavo V. Barroso
 * Created: 10/08/2022
 * Last modified: 18/04/2022
 *
 */


#ifndef _MIGRATION_H_
#define _MIGRATION_H_

#include "AbstractOperator.hpp"
#include "Graph.hpp"

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
    prevParams_.addParameters(getParameters());
    setLittleMat_();
    testFlow_();
    setUpMatrices_(sslib);
  }

  Migration(const Eigen::MatrixXd& migMat, std::shared_ptr<bpp::IntervalConstraint> ic, const SumStatsLibrary& sslib):
  AbstractOperator(),
  littleMigMat_(migMat)
  {
    for(size_t i = 0; i < sslib.getPopIndices().size(); ++i) // for each population modeled in epoch i
    {
      for(size_t j = 0; j < sslib.getPopIndices().size(); ++j)
      {
        if((i != j) && (littleMigMat_(i, j) != 0.))
        {
          //std::shared_ptr<bpp::Parameter> param = std::make_shared<bpp::Parameter>("m_" + bpp::TextTools::toString((*itI)) + bpp::TextTools::toString((*itJ)), initValues[idx], ic);
          //addParameter_(param.get());
          addParameter_(new bpp::Parameter("m_" + bpp::TextTools::toString(i) + "_" + bpp::TextTools::toString(j), littleMigMat_(i, j), ic));
        }
      }
    }

    prevParams_.addParameters(getParameters());
    testFlow_();
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

  void testFlow_();

  size_t fetchNumPops_();

};

#endif
