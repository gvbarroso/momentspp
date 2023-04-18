/*
 * Authors: Gustavo V. Barroso
 * Created: 31/10/2022
 * Last modified: 13/04/2023
 *
 */

#include <Eigen/Core>
#include <Eigen/Dense>

#include "Demes.hpp"

void Demes::parse_(const std::string& fileName)
{
  std::cout << "\nparsing " << fileName << "\n";
  model_ = YAML::LoadFile(fileName);

  if(model_["description"])
    std::cout << "model: " << model_["description"] << "\n\n";

  if(model_["defaults"])
    throw bpp::Exception("moments++ does not allow Demes defaults!");

  if(model_["time_units"].as<std::string>() != "generations")
    throw bpp::Exception("moments++ requires Demes [time_units] to be \"generations\"!");

  std::vector<std::vector<std::shared_ptr<Population>>> popsInv(0);

  for(YAML::const_iterator it = model_.begin(); it != model_.end(); ++it)
  {
    if(it->first.as<std::string>() == "demes")
    {
      YAML::Node pops = it->second;

      for(size_t i = 0; i < pops.size(); ++i) // deme by deme
      {
        std::string name = pops[i]["name"].as<std::string>(); // stats file should match pop names in Demes file
        std::string des = "none";

        if(pops[i]["description"])
          des = pops[i]["description"].as<std::string>();

        std::vector<std::shared_ptr<Population>> singlePopOverTime(0); // pop i is represented by a series of constant-size populations
        YAML::Node popEpochs = pops[i]["epochs"];

        for(size_t j = 0; j < popEpochs.size(); ++j) // epochs of focal pop i as they appear in the Demes (YAML) file
        {
          size_t startTime = std::numeric_limits<int>::max();
          if(popEpochs[j]["start_time"] && popEpochs[j]["start_time"].as<std::string>() != ".inf")
            startTime = popEpochs[j]["start_time"].as<int>();

          size_t endTime = 0;
          if(popEpochs[j]["end_time"])
            endTime = popEpochs[j]["end_time"].as<int>();

          size_t size = 0;
          if(popEpochs[j]["start_size"])
          {
            size = popEpochs[j]["start_size"].as<int>();

            if(popEpochs[j]["end_size"] && popEpochs[j]["end_size"] != popEpochs[j]["start_size"])
              throw bpp::Exception("moments++ requires a deme's [start_size] and [end_size] to be equal within each epoch!");
          }

          // each instance of pop i (one per epoch) is treated as a different Population object in moments++
          std::shared_ptr<Population> p = std::make_shared<Population>(name, des, i, startTime, endTime, size, false);
          singlePopOverTime.push_back(p);
        }

        if(pops[i]["ancestors"])
        {
          if(i == 0)
            throw bpp::Exception("ancestor(s) specified for first population of Demes file! Please list them in chronological order!");

          std::shared_ptr<Population> child = singlePopOverTime.front(); // most ancient instance of pop i

          if(pops[i]["ancestors"].size() == 1)
          {
            std::string ancName = pops[i]["ancestors"][0].as<std::string>();

            for(size_t k = 0; k < popsInv.size(); ++k)
            {
              if(popsInv[k].front()->getName() == ancName)
              {
                std::shared_ptr<Population> parent = popsInv[k].front();

                child->setLeftParent(parent);
                child->setRightParent(parent);

                if(child->getStartTime() == std::numeric_limits<int>::max())
                  child->setStartTime(parent->getEndTime());

                if(child->getSize() == 0)
                  child->setSize(parent->getSize());

                break;
              }
            }
          }

          else if(pops[i]["ancestors"].size() == 2)
          {
            std::string ancNameFirst = pops[i]["ancestors"][0].as<std::string>();
            std::string ancNameSecond = pops[i]["ancestors"][1].as<std::string>();

            for(size_t k = 0; k < popsInv.size(); ++k)
            {
              if(popsInv[k].front()->getName() == ancNameFirst)
              {
                std::shared_ptr<Population> parent = popsInv[k].front();
                child->setLeftParent(parent);
              }

              else if(popsInv[k].front()->getName() == ancNameSecond)
              {
                std::shared_ptr<Population> parent = popsInv[k].front();
                child->setRightParent(parent);
              }
            }

            double f = pops[i]["proportions"][0].as<double>();
            double g = pops[i]["proportions"][1].as<double>();

            // TODO forward f and g to Admixture operator

            if((f + g) != 1.)
              throw bpp::Exception("admixture proportions in Demes file don't sum to 1.0!");

            if(child->getSize() == 0)
              throw bpp::Exception("demes with two ancestors must have specified start_sizes!");

            if(child->getStartTime() == std::numeric_limits<int>::max())
              std::cout << "Fix this\n";
          }

          else if(pops[i]["ancestors"].size() > 2)
            throw bpp::Exception("more than two ancestors for a single population!");
        }

        for(size_t k = 1; k < singlePopOverTime.size(); ++k)
        {
          singlePopOverTime[k]->setLeftParent(singlePopOverTime[k - 1]);
          singlePopOverTime[k]->setRightParent(singlePopOverTime[k - 1]);
          singlePopOverTime[k]->setStartTime(singlePopOverTime[k - 1]->getEndTime());
        }

        popsInv.push_back(singlePopOverTime);
      }

      // slice time into epochs based on populations time boundaries
      std::vector<size_t> timeBounds(0);
      for(size_t i = 0; i < popsInv.size(); ++i)
      {
        for(size_t j = 0; j < popsInv[i].size(); ++j)
        {
          timeBounds.push_back(popsInv[i][j]->getStartTime());
          timeBounds.push_back(popsInv[i][j]->getEndTime());
        }
      }

      std::sort(std::begin(timeBounds), std::end(timeBounds), std::greater<int>()); // descending order
      timeBounds.erase(std::unique(std::begin(timeBounds), std::end(timeBounds)), std::end(timeBounds));

      // first pass: split populations that span more than one epoch
      for(size_t i = 1; i < timeBounds.size(); ++i)
      {
        size_t epochStart = timeBounds[i - 1];
        size_t epochEnd = timeBounds[i];

        for(size_t j = 0; j < popsInv.size(); ++j)
        {
          for(auto itPop = std::begin(popsInv[j]); itPop < std::end(popsInv[j]); ++itPop)
          {
            size_t popStart = (*itPop)->getStartTime();
            size_t popEnd = (*itPop)->getEndTime();

            if(popStart == epochStart && popEnd < epochEnd)
            {
              std::shared_ptr<Population> splitLeft = std::make_shared<Population>(*(*itPop).get());
              std::shared_ptr<Population> splitRight = std::make_shared<Population>(*(*itPop).get());

              splitLeft->setEndTime(epochEnd);

              splitRight->setStartTime(splitLeft->getEndTime());
              splitRight->setEndTime(popEnd);
              splitRight->setLeftParent(splitLeft);
              splitRight->setRightParent(splitLeft);

              // pseudo code
              itPop = popsInv[j].erase(itPop);
              itPop = popsInv[j].insert(itPop, splitRight);
              itPop = popsInv[j].insert(itPop, splitLeft);
              itPop = std::next(itPop, 1);
            }
          }
        }
      }

      // second pass: organize pops within epochs
      pops_.resize(timeBounds.size() - 1);
      for(size_t i = 1; i < timeBounds.size(); ++i)
      {
        size_t epochStart = timeBounds[i - 1];
        size_t epochEnd = timeBounds[i];

        for(size_t j = 0; j < popsInv.size(); ++j)
        {
          for(auto itPop = std::begin(popsInv[j]); itPop < std::end(popsInv[j]); ++itPop)
          {
            size_t popStart = (*itPop)->getStartTime();
            size_t popEnd = (*itPop)->getEndTime();

            if(popStart == epochStart && popEnd == epochEnd)
            {
              pops_[i - 1].push_back(*itPop);
              break;
            }
          }
        }
      }
    }

    else if(it->first.as<std::string>() == "migrations")
    {
      std::cout << "build littleMigMat\n"; // NOTE: it is epoch-specific
    }

    else if(it->first.as<std::string>() == "pulses")
    {
      std::cout << "TODO\n"; // NOTE: it is epoch-specific
    }

    else if(it->first.as<std::string>() == "metadata")
    {
      std::cout << "extract u, r...\n";
      //double u = 0.;
      //double r = 0.;
    }
  }

  for(auto it = std::begin(pops_); it != std::end(pops_); ++it)
  {
    for(size_t x = 0; x < it->size(); ++x)
      (*it)[x]->printAttributes(std::cout);
  }
}
