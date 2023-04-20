/*
 * Authors: Gustavo V. Barroso
 * Created: 10/08/2022
 * Last modified: 20/04/2022
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
    setUpMatrices_(sslib);
  }

  Migration(const Eigen::MatrixXd& migMat, std::shared_ptr<bpp::IntervalConstraint> ic, const SumStatsLibrary& sslib):
  AbstractOperator(),
  littleMigMat_(migMat)
  {
    // for each pair of populations modeled in the epoch to which *this operator belongs
    for(size_t i = 0; i < sslib.getPopIndices().size(); ++i)
    {
      for(size_t j = 0; j < sslib.getPopIndices().size(); ++j)
      {
        if(i != j)
          addParameter_(new bpp::Parameter("m_" + bpp::TextTools::toString(i) + "_" + bpp::TextTools::toString(j), littleMigMat_(i, j), ic));
      }
    }

    prevParams_.addParameters(getParameters());
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

  void testFlow();

private:
  void setUpMatrices_(const SumStatsLibrary& sslib) override;

  void updateMatrices_() override;

  void setLittleMat_();

  size_t fetchNumPops_();

};

#endif
