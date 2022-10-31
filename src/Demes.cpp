/*
 * Authors: Gustavo V. Barroso
 * Created: 31/10/2022
 * Last modified: 31/10/2022
 *
 */

#include "Demes.hpp"

void parse(const std::string& fileName)
{
  //parser.HandleNextDocument(test);
  YAML::Node test = YAML::LoadFile(fileName);


  std::ofstream fout("test_out.yaml");
  fout << test;
}
