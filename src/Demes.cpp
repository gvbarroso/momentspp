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

  for(YAML::const_iterator it=test.begin();it!=test.end();++it)
    std::cout << "node: " << it->first.as<std::string>() << " is\n ";// << it->second.as<std::string>() << "\n";

}
