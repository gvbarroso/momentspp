/*
 * Authors: Gustavo V. Barroso
 * Created: 10/08/2022
 * Last modified: 30/05/2023
 *
 */


#ifndef _MIGRATION_H_
#define _MIGRATION_H_

#include "AbstractOperator.hpp"
#include "Graph.hpp"

class Migration: public AbstractOperator
{

private:
  Eigen::Matrix<long double, Eigen::Dynamic, Eigen::Dynamic> littleMigMat_; // P x P

public:
  Migration(const bpp::ParameterList migParams, const SumStatsLibrary& sslib):
  AbstractOperator(sslib.getPopIndices()),
  littleMigMat_()
  {
    includeParameters_(migParams);
    prevParams_.addParameters(getParameters());
    setLittleMat_();
    setUpMatrices_(sslib);
  }

  Migration(const Eigen::Matrix<long double, Eigen::Dynamic, Eigen::Dynamic>& migMat, std::shared_ptr<bpp::IntervalConstraint> ic, const SumStatsLibrary& sslib):
  AbstractOperator(sslib.getPopIndices()),
  littleMigMat_(migMat)
  {
    // for each pair of populations modeled in the epoch to which *this operator belongs
    for(size_t i = 0; i < popIndices_.size(); ++i)
    {
      size_t id = popIndices_[i];

      for(size_t j = 0; j < popIndices_.size(); ++j)
      {
        size_t jd = popIndices_[j];

        if(id != jd)
          addParameter_(new bpp::Parameter("m_" + bpp::TextTools::toString(id) + "_" + bpp::TextTools::toString(jd), littleMigMat_(i, j), ic));
      }
    }

    prevParams_.addParameters(getParameters());
    setUpMatrices_(sslib);
  }

  virtual Migration* clone() const override
  {
    return new Migration(*this);
  }

  const Eigen::Matrix<long double, Eigen::Dynamic, Eigen::Dynamic>& getLittleMigMat()
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
