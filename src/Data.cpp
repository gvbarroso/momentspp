/*
 * Authors: Gustavo V. Barroso
 * Created: 22/09/2022
 * Last modified: 01/11/2023
 *
 */


#include "Data.hpp"


void Data::parse_(const std::string& file)
{
  ssl_.readStatsFromFile(file);
}

