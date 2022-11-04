/*
 * Authors: Gustavo V. Barroso
 * Created: 31/10/2022
 * Last modified: 02/10/2022
 *
 */

#include "Demes.hpp"

void Demes::parse(const std::string& fileName)
{
  std::cout << "parsing " << fileName << "\n";
  YAML::Node test = YAML::LoadFile(fileName);

  std::ofstream fout("test.yaml");
  fout << test;
  fout.close();

  for(YAML::const_iterator it = test.begin(); it != test.end(); ++it)
  {
    YAML::Node tmp = *it;

    std::cout << it->first.as<std::string>() << "\t";

    if(it->second.size() == 0)
      std::cout << it->second;

    else
    {
      for(size_t i = 0; i < it->second.size(); ++i)
      {
        std::cout << it->second[i];

        if(i < it->second.size() - 1)
          std::cout << "\t";
      }
    }

    std::cout << "\n";
  }

}
