/*
 * Authors: Gustavo V. Barroso
 * Created: 08/12/2022
 * Last modified: 08/12/2022
 *
 */


#include "Utils.hpp"

std::vector<bpp::ParameterList> Utils::fetchParamsTable(const std::string& fileName)
{
  std::vector<bpp::ParameterList> listOfModels(0); // each combination of parameters defines a "model"
  std::vector<std::string> paramNames(0);

  std::ifstream paramsFile;
  paramsFile.open(fileName);

  if(paramsFile.is_open())
  {
    std::string line = "";
    std::vector<std::string> splitLine(0);

    getline(paramsFile, line); // header

    boost::split(splitLine, line, [=](char c) { return c == ','; }); // comma-separated values
    paramNames = splitLine;

    size_t numPops = std::sqrt(splitLine.size() - 3); // integer, # of N_i's + m_ij's

    while(getline(paramsFile, line))
    {
      bpp::ParameterList pl;
      boost::split(splitLine, line, [=](char c) { return c == ','; });

      for(size_t i = 0; i < splitLine.size(); ++i)
        pl.addParameter(new bpp::Parameter(paramNames[i], splitLine[i]));

      listOfModels.push_back(pl);
    }

    paramsFile.close();
  }

  else
    throw bpp::Exception("UTILS::Could not open " + fileName);
}

// TODO when models are specified by YAML Demes files, we can have a file with a list of names of Demes files:
// model_1.yaml
// model_2.yaml
// etc
void Utils::parseListOfModelsFile(const std::string& fileName);
