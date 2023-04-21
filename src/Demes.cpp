/*
 * Authors: Gustavo V. Barroso
 * Created: 31/10/2022
 * Last modified: 21/04/2023
 *
 */

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

  // two vectors that must live throughout this method
  std::vector<size_t> timeBounds(0);
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
              throw bpp::Exception("deme's [start_size] and [end_size] must be equal within each epoch!");
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

            double f = pops[i]["proportions"][0].as<double>();
            double g = pops[i]["proportions"][1].as<double>();

            // forward f and g to Admixture operator

            if((f + g) != 1.)
              throw bpp::Exception("ancestral admixture proportions in Demes file don't sum to 1.0!");

            if(child->getSize() == 0)
              throw bpp::Exception("demes with two ancestors must have specified start_size's!");

            if(child->getStartTime() == std::numeric_limits<int>::max())
              throw bpp::Exception("demes with two ancestors must have specified start_time's!");

            for(size_t k = 0; k < popsInv.size(); ++k)
            {
              if(popsInv[k].front()->getName() == ancNameFirst)
              {
                for(size_t l = 0; l < popsInv[k].size(); ++l)
                {
                  if(popsInv[k][l]->getEndTime() == child->getStartTime())
                  {
                    std::shared_ptr<Population> parent = popsInv[k][l];
                    child->setLeftParent(parent);
                  }
                }
              }

              else if(popsInv[k].front()->getName() == ancNameSecond)
              {
                for(size_t l = 0; l < popsInv[k].size(); ++l)
                {
                  if(popsInv[k][l]->getEndTime() == child->getStartTime())
                  {
                    std::shared_ptr<Population> parent = popsInv[k][l];
                    child->setRightParent(parent);
                  }
                }
              }
            }

            if(child->getLeftParent() == nullptr || child->getRightParent() == nullptr)
              throw bpp::Exception("could not pinpoint parents of pop " + child->getName() + " in time.");
          }

          else if(pops[i]["ancestors"].size() > 2)
            throw bpp::Exception("more than two ancestors for a single population!");
        }

        for(size_t k = 1; k < singlePopOverTime.size(); ++k)
        {
          singlePopOverTime[k]->setLeftParent(singlePopOverTime[k - 1]);
          singlePopOverTime[k]->setRightParent(singlePopOverTime[k - 1]);
          singlePopOverTime[k]->setStartTime(singlePopOverTime[k - 1]->getEndTime());

          if(singlePopOverTime[k]->getSize() == 0)
            singlePopOverTime[k]->setSize(singlePopOverTime[k - 1]->getSize());
        }

        popsInv.push_back(singlePopOverTime);
      }

      // slice time into epochs based on populations time boundaries
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

      size_t numEpochs = timeBounds.size() - 1;

      pops_.resize(numEpochs);
      migRates_.resize(numEpochs);
      mutRates_.reserve(numEpochs);
      recRates_.reserve(numEpochs);

      // second pass: organize pops within epochs
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

      for(size_t i = 0; i < numEpochs; ++i)
      {
        size_t p = pops_[i].size();
        Eigen::MatrixXd mat(p, p);
        mat.setZero();

        migRates_[i] = mat;
      }

      // sets defaults to mutation and recombination rates
      double rate = 1e-8;

      for(size_t i = 0; i < numEpochs; ++i)
      {
        mutRates_.emplace_back(rate);
        recRates_.emplace_back(rate);
      }
    }

    else if(it->first.as<std::string>() == "migrations")
    {
      YAML::Node migs = it->second;

      std::string source = "";
      std::string dest = "";
      double rate = 0.;
      size_t startTime = 0;
      size_t endTime = 0;

      for(size_t i = 0; i < migs.size(); ++i) // mig period by mig period
      {
        if(migs[i]["source"])
          source = migs[i]["source"].as<std::string>();

        else
          throw bpp::Exception("'migrations' field in Demes file must explicitly speficy 'source'!");

        if(migs[i]["dest"])
          dest = migs[i]["dest"].as<std::string>();

        else
          throw bpp::Exception("'migrations' field in Demes file must explicitly speficy 'dest'!");

        if(migs[i]["rate"])
          rate = migs[i]["rate"].as<double>();

        else
          throw bpp::Exception("'migrations' field in Demes file must explicitly speficy 'rate'!");

        if(migs[i]["start_time"])
          startTime = migs[i]["start_time"].as<size_t>();

        else
          throw bpp::Exception("'migrations' field in Demes file must explicitly speficy 'start_time'!");

        if(migs[i]["end_time"])
          endTime = migs[i]["end_time"].as<size_t>();

        else
          throw bpp::Exception("'migrations' field in Demes file must explicitly speficy 'end_time'!");

        bool match = 0;
        for(size_t j = 1; j < timeBounds.size(); ++j)
        {
          if(startTime == timeBounds[j - 1] && endTime == timeBounds[j])
          {
            int row = -1;
            int col = -1;

            for(size_t k = 0; k < pops_[j - 1].size(); ++ k)
            {
              if(pops_[j - 1][k]->getName() == source)
                row = k;

              else if(pops_[j - 1][k]->getName() == dest)
                col = k;
            }

            migRates_[j - 1](row, col) = rate;

            match = 1;
            break;
          }
        }

        if(!match)
          throw bpp::Exception("start_time and end_time of 'migrations' in Demes file do not match the span of any epoch!");

      }
    }

    else if(it->first.as<std::string>() == "pulses")
    {
      YAML::Node pulses = it->second;

      for(size_t i = 0; i < pulses.size(); ++i) // pulse by pulse
      {
        YAML::Node sources = pulses[i]["sources"];
        YAML::Node fs = pulses[i]["proportions"];

        std::string source = "";
        double f = -1.;

        if(sources.size() == 1)
          source = sources[0].as<std::string>();

        else
          throw bpp::Exception("only a single 'source' per admixture 'pulse' is allowed in Demes model");

        if(fs.size() == 1)
          f = fs[0].as<double>(); // proportion of source ancestry

        else
          throw bpp::Exception("only a single 'proportion' per admixture 'pulse' is allowed in Demes model");

        std::string dest = pulses[i]["dest"].as<std::string>();
        size_t time = pulses[i]["time"].as<size_t>();

        std::cout << source << "~~~~~>" << dest << " [" << f << "] at " << time << " gens. ago\n";
      }
    }

    // move to mutation and recombination rates (optional, defaults give above)
    if(it->first.as<std::string>() == "mutation")
    {
      YAML::Node muts = it->second;

      double rate = 0.;
      size_t startTime = timeBounds.front();
      size_t endTime = 0;

      for(size_t i = 0; i < muts.size(); ++i) // mut period by mut period
      {
        if(muts[i]["rate"])
          rate = muts[i]["rate"].as<double>();

        if(muts[i]["start_time"])
          startTime = muts[i]["start_time"].as<size_t>();

        if(muts[i]["end_time"])
          endTime = muts[i]["end_time"].as<size_t>();

        bool match = 1;
        for(size_t j = 1; j < timeBounds.size(); ++j)
        {
          if((startTime == timeBounds[j - 1] && endTime == timeBounds[j]) || (startTime == timeBounds.front() && endTime == 0))
            mutRates_[j - 1] = rate;

          else
            match = 0;
        }

        if(!match)
          throw bpp::Exception("start_time and end_time of 'mutation' in Demes file do not match the span of any epoch!");
      }
    }

    else if(it->first.as<std::string>() == "recombination")
    {
      YAML::Node recs = it->second;

      double rate = 1e-8;
      size_t startTime = timeBounds.front();
      size_t endTime = 0;

      for(size_t i = 0; i < recs.size(); ++i) // rec period by rec period
      {
        if(recs[i]["rate"])
          rate = recs[i]["rate"].as<double>();

        if(recs[i]["start_time"])
          startTime = recs[i]["start_time"].as<size_t>();

        if(recs[i]["end_time"])
          endTime = recs[i]["end_time"].as<size_t>();

        bool match = 1;
        for(size_t j = 1; j < timeBounds.size(); ++j)
        {
          if((startTime == timeBounds[j - 1] && endTime == timeBounds[j]) || (startTime == timeBounds.front() && endTime == 0))
            recRates_[j - 1] = rate;

          else
            match = 0;
        }

        if(!match)
          throw bpp::Exception("start_time and end_time of 'recombination' in Demes file do not match the span of any epoch!");
      }
    }
  }

  /*for(auto it = std::begin(pops_); it != std::end(pops_); ++it)
  {
    for(size_t x = 0; x < it->size(); ++x)
      (*it)[x]->printAttributes(std::cout);
  }*/
}
