/*
 * Authors: Gustavo V. Barroso
 * Created: 31/10/2022
 * Last modified: 11/04/2023
 *
 */

#include "Demes.hpp"

void Demes::parse_(const std::string& fileName)
{
  std::cout << "\nparsing " << fileName << "\n";
  model_ = YAML::LoadFile(fileName);

  if(model_["description"])
    std::cout << "model: " << model_["description"] << "\n\n";

  if(model_["time_units"].as<std::string>() != "generations")
    throw bpp::Exception("moments++ demes model must specify time in units of generations! [time_units]");

  for(YAML::const_iterator it = model_.begin(); it != model_.end(); ++it)
  {
    if(it->first.as<std::string>() == "demes")
    {
      assert(it->second.size() > 0);

      for(size_t i = 0; i < it->second.size(); ++i)
      {
        std::cout << it->second[i];

        if(it->second[i].IsMap())
          std::cout << "map in\n";

        if(i < it->second.size() - 1)
          std::cout << "\t";
      }

      std::cout << "\n";
    }

    else if(it->first.as<std::string>() == "migrations")
    {
      std::cout << "TODO\n";
    }

    else if(it->first.as<std::string>() == "pulsess")
    {
      std::cout << "TODO\n";
    }
  }

  std::cout << "\n";
}
