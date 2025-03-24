/*
 * Authors: Gustavo V. Barroso
 * Created: 21/03/2025
 * Last modified: 24/04/2025
 *
 */


#ifndef _NEUTRAL_MIGRATION_H_
#define _NEUTRAL_MIGRATION_H_

#include "AbstractOperator.hpp"
#include "Graph.hpp"

class NeutralMigration: public AbstractOperator
{

private:
  Eigen::MatrixXd littleMigMat_; // P x P

public:
  NeutralMigration(const bpp::ParameterList migParams, const SumStatsLibrary& sslib):
  AbstractOperator(sslib.getPopIndices()),
  littleMigMat_()
  {
    includeParameters_(migParams);
    prevParams_.addParameters(getParameters());
    setLittleMat_();
    setUpMatrices_(sslib);
  }

  NeutralMigration(const Eigen::MatrixXd& migMat, std::shared_ptr<bpp::IntervalConstraint> ic, const SumStatsLibrary& sslib):
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
