/*
 * Authors: Gustavo V. Barroso
 * Created: 31/10/2022
 * Last modified: 18/04/2023
 *
 */


#ifndef _DEMES_H_
#define _DEMES_H_

#include <fstream>
#include <vector>
#include <memory>
#include <utility>
#include <algorithm>
#include <map>
#include <ios>
#include <limits>

#include <yaml-cpp/node/parse.h>
#include <yaml-cpp/node/ptr.h>
#include <yaml-cpp/node/emit.h>
#include <yaml-cpp/node/convert.h>
#include <yaml-cpp/node/impl.h>
#include <yaml-cpp/node/node.h>
#include <yaml-cpp/node/iterator.h>
#include <yaml-cpp/node/type.h>
#include <yaml-cpp/node/detail/node_ref.h>
#include <yaml-cpp/node/detail/node.h>
#include <yaml-cpp/node/detail/memory.h>
#include <yaml-cpp/node/detail/iterator.h>
#include <yaml-cpp/node/detail/impl.h>
#include <yaml-cpp/node/detail/iterator_fwd.h>
#include <yaml-cpp/node/detail/node_iterator.h>
#include <yaml-cpp/node/detail/node_data.h>
#include <yaml-cpp/yaml.h>
#include <yaml-cpp/parser.h>
#include <yaml-cpp/emitter.h>
#include <yaml-cpp/depthguard.h>
#include <yaml-cpp/noexcept.h>
#include <yaml-cpp/anchor.h>
#include <yaml-cpp/stlemitter.h>
#include <yaml-cpp/binary.h>
#include <yaml-cpp/emitterstyle.h>
#include <yaml-cpp/mark.h>
#include <yaml-cpp/ostream_wrapper.h>
#include <yaml-cpp/emitfromevents.h>
#include <yaml-cpp/traits.h>
#include <yaml-cpp/eventhandler.h>
#include <yaml-cpp/dll.h>
#include <yaml-cpp/emitterdef.h>
#include <yaml-cpp/exceptions.h>
#include <yaml-cpp/null.h>

#include <Bpp/Exceptions.h>

#include "Population.hpp"

class Demes
{

private:
  // demes file describes the demographic structure allowed
  // param values are optimized by moments++ but no "new param" is included (eg, admix event)
  YAML::Node model_;

  // pop-id -> Population*, one per epoch
  std::vector<std::vector<std::shared_ptr<Population>>> pops_;

public:
  Demes(const std::string& file):
  model_(),
  pops_(0)
  {
    parse_(file);
  }

  Demes():
  model_(),
  pops_(0)
  { }

public:
  ~Demes()
  { }

  const YAML::Node& getModel_()
  {
    return model_;
  }

  const std::vector<std::vector<std::shared_ptr<Population>>>& getPopMaps()
  {
    return pops_;
  }

  const std::vector<std::vector<std::shared_ptr<Population>>>& getPopMaps() const
  {
    return pops_;
  }

  size_t getNumEpochs()
  {
    return pops_.size();
  }

  size_t getNumEpochs() const
  {
    return pops_.size();
  }

  size_t getNumPops(size_t epoch)
  {
    return pops_[epoch].size();
  }

  void write(const std::string& fileName);

private:
  void parse_(const std::string& fileName);

};

#endif
