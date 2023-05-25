/*
 * Authors: Gustavo V. Barroso
 * Created: 31/10/2022
 * Last modified: 24/05/2023
 *
 */

#include "Demes.hpp"

void Demes::parse_(const std::string& fileName)
{
  std::cout << "\nParsing " << fileName << "\n";
  model_ = YAML::LoadFile(fileName);

  if(model_["description"])
    std::cout << "model: " << model_["description"] << "\n\n";

  if(model_["defaults"])
    throw bpp::Exception("Demes::defaults are not allowed!");

  if(model_["time_units"].as<std::string>() != "generations")
    throw bpp::Exception("Demes:: [time_units] must be \"generations\"!");

  // two vectors that must live throughout this method
  std::vector<size_t> timeBounds(0);
  std::vector<std::vector<std::shared_ptr<Population>>> popsOverTime(0);

  if(model_["pulses"]) // first pass on "pulses" to store times
  {
    YAML::Node pulses = model_["pulses"];

    for(size_t i = 0; i < pulses.size(); ++i)
      timeBounds.push_back(pulses[i]["time"].as<size_t>() - 1); // -1 to allow 1-gen epochs
  }

  for(YAML::const_iterator it = model_.begin(); it != model_.end(); ++it)
  {
    if(it->first.as<std::string>() == "demes")
    {
      YAML::Node pops = it->second;

      for(size_t i = 0; i < pops.size(); ++i) // deme by deme, pop index (i) is fixed by order of listing demes in Demes file
      {
        // NOTE: stats file should match pop names and/or indices in Demes file
        std::string name = pops[i]["name"].as<std::string>();
        std::string des = "none";

        if(pops[i]["description"])
          des = pops[i]["description"].as<std::string>();

        size_t startTime = std::numeric_limits<int>::max();
        if(pops[i]["start_time"] && pops[i]["start_time"].as<std::string>() != ".inf")
          startTime = pops[i]["start_time"].as<int>();

        std::vector<std::shared_ptr<Population>> singlePopOverTime(0); // deme i is represented by a series of populations of piece-wise constant Ne
        YAML::Node popEpochs = pops[i]["epochs"];

        for(size_t j = 0; j < popEpochs.size(); ++j) // epochs of focal pop i as they appear in the Demes file
        {
          size_t endTime = 0;
          if(popEpochs[j]["end_time"])
            endTime = popEpochs[j]["end_time"].as<int>();

          size_t size = 0;
          if(popEpochs[j]["start_size"])
          {
            size = popEpochs[j]["start_size"].as<int>();

            if(popEpochs[j]["end_size"] && popEpochs[j]["end_size"] != popEpochs[j]["start_size"])
              throw bpp::Exception("Demes::[start_size] and [end_size] must be equal within each epoch!");
          }

          // each instance of pop i (one per epoch) is treated as a different Population object in moments++
          singlePopOverTime.push_back(std::make_shared<Population>(name, des, i, startTime, endTime, size, false));
        }

        if(pops[i]["ancestors"])
        {
          std::shared_ptr<Population> child = singlePopOverTime.front(); // most ancient instance of pop i

          if(pops[i]["ancestors"].size() == 1)
          {
            std::string ancName = pops[i]["ancestors"][0].as<std::string>();

            for(size_t k = 0; k < popsOverTime.size(); ++k)
            {
              if(popsOverTime[k].back()->getName() == ancName)
              {
                if(startTime == std::numeric_limits<int>::max())
                {
                  std::shared_ptr<Population> parent = popsOverTime[k].back();

                  child->setLeftParent(parent);
                  child->setRightParent(parent);
                  child->setStartTime(parent->getEndTime());

                  if(child->getSize() == 0)
                    child->setSize(parent->getSize());
                }

                else
                {
                  bool match = 0;
                  for(size_t l = 0; l < popsOverTime[k].size(); ++l)
                  {
                    if(startTime == popsOverTime[k][l]->getEndTime())
                    {
                      std::shared_ptr<Population> parent = popsOverTime[k][l];

                      child->setLeftParent(parent);
                      child->setRightParent(parent);

                      if(child->getSize() == 0)
                        child->setSize(parent->getSize());

                      match = 1;
                    }
                  }

                  if(match == 0)
                    throw bpp::Exception("Demes::could not find ancestor of pop " + name + "at specified start_time!");
                }

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

            if((f + g) != 1.)
              throw bpp::Exception("Demes::ancestral admixture proportions don't sum to 1.0!");

            if(child->getSize() == 0)
              throw bpp::Exception("Demes::demes with two ancestors must have specified start_size's!");

            if(child->getStartTime() == std::numeric_limits<int>::max())
              throw bpp::Exception("Demes::demes with two ancestors must have specified start_time's!");

            child->setProportions(std::make_pair(f, g));

            for(size_t k = 0; k < popsOverTime.size(); ++k)
            {
              if(popsOverTime[k].front()->getName() == ancNameFirst)
              {
                for(size_t l = 0; l < popsOverTime[k].size(); ++l)
                {
                  if(popsOverTime[k][l]->getEndTime() == child->getStartTime())
                    child->setLeftParent(popsOverTime[k][l]); // treated as source population in pulse of admixture
                }
              }

              else if(popsOverTime[k].front()->getName() == ancNameSecond)
              {
                for(size_t l = 0; l < popsOverTime[k].size(); ++l)
                {
                  if(popsOverTime[k][l]->getEndTime() == child->getStartTime())
                    child->setRightParent(popsOverTime[k][l]); // child will "copy" stats from right parent
                }
              }
            }

            if(child->getLeftParent() == nullptr || child->getRightParent() == nullptr)
              throw bpp::Exception("Demes::could not find ancestors of pop " + child->getName() + " (both must have an epoch's end_time matching " + child->getName() + "'s start_time)");
          }

          else if(pops[i]["ancestors"].size() > 2)
            throw bpp::Exception("Demes::more than two ancestors specified for a single population!");
        } // exists 'ancestors' field

        for(size_t k = 1; k < singlePopOverTime.size(); ++k)
        {
          singlePopOverTime[k]->setLeftParent(singlePopOverTime[k - 1]);
          singlePopOverTime[k]->setRightParent(singlePopOverTime[k - 1]);
          singlePopOverTime[k]->setStartTime(singlePopOverTime[k - 1]->getEndTime());

          if(singlePopOverTime[k]->getSize() == 0)
            singlePopOverTime[k]->setSize(singlePopOverTime[k - 1]->getSize());
        }

        popsOverTime.push_back(singlePopOverTime);
      } // ends loop over demes

      // slice time into epochs based on populations time boundaries
      for(size_t i = 0; i < popsOverTime.size(); ++i)
      {
        for(size_t j = 0; j < popsOverTime[i].size(); ++j)
        {
          timeBounds.push_back(popsOverTime[i][j]->getStartTime());
          timeBounds.push_back(popsOverTime[i][j]->getEndTime());

          // NOTE introducing 1-generation epoch to handle admixture that forms new population
          if(popsOverTime[i][j]->hasDistinctParents())
            timeBounds.push_back(popsOverTime[i][j]->getStartTime() - 1);
        }
      }

      std::sort(std::begin(timeBounds), std::end(timeBounds), std::greater<int>()); // descending order
      timeBounds.erase(std::unique(std::begin(timeBounds), std::end(timeBounds)), std::end(timeBounds));

      // first pass: split populations that span more than one epoch
      for(size_t i = 1; i < timeBounds.size(); ++i)
      {
        size_t epochStart = timeBounds[i - 1];
        size_t epochEnd = timeBounds[i];

        //std::cout << "\nepoch start: " << epochStart << ", epoch end: " << epochEnd << "\n";

        for(size_t j = 0; j < popsOverTime.size(); ++j) // for each "deme"
        {
          //std::cout << popsOverTime[j].front()->getName() << "\n";

          for(auto itPop = std::begin(popsOverTime[j]); itPop < std::end(popsOverTime[j]); ++itPop)
          {
            size_t popStart = (*itPop)->getStartTime();
            size_t popEnd = (*itPop)->getEndTime();

            //std::cout << "pop start: " << popStart << ", pop end: " << popEnd << "\n";

            if(popStart == epochStart && popEnd < epochEnd) // must split
            {
              std::shared_ptr<Population> splitLeft = std::make_shared<Population>(*(*itPop).get()); // more ancient instance of pop
              std::shared_ptr<Population> splitRight = std::make_shared<Population>(*(*itPop).get()); // more recent instance of pop

              splitLeft->setEndTime(epochEnd);

              splitRight->setStartTime(splitLeft->getEndTime());
              splitRight->setEndTime(popEnd);
              splitRight->setLeftParent(splitLeft);
              splitRight->setRightParent(splitLeft);

              itPop = popsOverTime[j].erase(itPop);
              itPop = popsOverTime[j].insert(itPop, splitRight);
              itPop = popsOverTime[j].insert(itPop, splitLeft);
              itPop = std::next(itPop, 1);
            }
          }
        }
      }

      /*for(size_t i = 0; i < popsOverTime.size(); ++i)
        for(size_t j = 0; j < popsOverTime[i].size(); ++j)
          popsOverTime[i][j]->printAttributes(std::cout);*/

      // epochs are formalized (Epoch objects are instantiated) in the main function
      // their skeleton are prepared here:
      size_t numEpochs = timeBounds.size() - 1;

      pops_.resize(numEpochs);
      migRates_.resize(numEpochs);
      pulses_.resize(numEpochs);
      mutRates_.reserve(numEpochs);
      recRates_.reserve(numEpochs);

      // second pass: organize pops within epochs
      for(size_t i = 0; i < numEpochs; ++i)
      {
        size_t epochStart = timeBounds[i];
        size_t epochEnd = timeBounds[i + 1];

        for(size_t j = 0; j < popsOverTime.size(); ++j)
        {
          for(auto itPop = std::begin(popsOverTime[j]); itPop < std::end(popsOverTime[j]); ++itPop)
          {
            size_t popStart = (*itPop)->getStartTime();
            size_t popEnd = (*itPop)->getEndTime();

            if(popStart == epochStart && popEnd == epochEnd)
            {
              pops_[i].push_back(*itPop);
              break;
            }
          }
        }
      }

      // inits parameters for the different operators
      for(size_t i = 0; i < numEpochs; ++i)
      {
        double rate = 1e-8; // default mu and r

        mutRates_.emplace_back(rate);
        recRates_.emplace_back(rate);

        size_t p = pops_[i].size();
        Eigen::MatrixXd mat(p, p);
        mat.setZero();

        migRates_[i] = mat; // littleMigMat_ inside Migration class
        pulses_[i] = mat; // littleAdmixMat_ inside Admixture class
      }

      // search for populations with two ancestors: pick one of them to copy moments from,
      // then apply Admixture as if it were a pulse (c.f. Model::linkMoments_())
      for(size_t j = 1; j < numEpochs; ++j)
      {
        for(size_t k = 0; k < pops_[j].size(); ++ k)
        {
          int row = -1;
          int col = k;

          if(pops_[j][k]->hasDistinctParents()) // admixture forms a new population
          {
            // search for (left) ancestral population in the current (1-generation) epoch
            for(size_t l = 0; l < pops_[j].size(); ++ l)
              if(pops_[j][l]->getName() == pops_[j][k]->getLeftParent()->getName())
                row = l;

            pulses_[j](row, col) = pops_[j][k]->getProportions().first;
          }
        }
      }
    } // exits 'demes' field of Demes file

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
          throw bpp::Exception("Demes::'migrations' field must explicitly speficy 'source'!");

        if(migs[i]["dest"])
          dest = migs[i]["dest"].as<std::string>();

        else
          throw bpp::Exception("Demes::'migrations' field must explicitly speficy 'dest'!");

        if(migs[i]["rate"])
          rate = migs[i]["rate"].as<double>();

        else
          throw bpp::Exception("Demes::'migrations' field must explicitly speficy 'rate'!");

        if(migs[i]["start_time"])
        {
          if(migs[i]["start_time"].as<std::string>() == ".inf")
            startTime = std::numeric_limits<int>::max();

          else
            startTime = migs[i]["start_time"].as<size_t>();
        }

        else
          throw bpp::Exception("Demes::'migrations' field must explicitly specify 'start_time'!");

        if(migs[i]["end_time"])
          endTime = migs[i]["end_time"].as<size_t>();

        else
          throw bpp::Exception("Demes::'migrations' field must explicitly speficy 'end_time'!");

        bool match = 0;
        for(size_t j = 1; j < timeBounds.size(); ++j)
        {
          // minding 1-gen epochs introduced by admixture
          if((startTime == timeBounds[j - 1] || (startTime - 1) == timeBounds[j - 1]) && endTime == timeBounds[j])
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
          throw bpp::Exception("Demes::start_time and end_time of 'migrations' in Demes file do not match the span of any epoch!");
      }
    } // exits 'migrations' field of Demes file

    else if(it->first.as<std::string>() == "pulses")
    {
      YAML::Node pulses = it->second;

      for(size_t i = 0; i < pulses.size(); ++i) // pulse by pulse
      {
        YAML::Node sources = pulses[i]["sources"];
        YAML::Node proportions = pulses[i]["proportions"];

        std::string source = "";
        double f = -1.;

        if(sources.size() == 1)
          source = sources[0].as<std::string>();

        else
          throw bpp::Exception("Demes::only a single 'source' pop. (specified within brackets) per admixture 'pulse' is allowed!");

        if(proportions.size() == 1)
          f = proportions[0].as<double>();

        else
          throw bpp::Exception("Demes::only a single 'proportion' (specified within brackets) per admixture 'pulse' is allowed!");

        std::string dest = pulses[i]["dest"].as<std::string>();
        size_t time = pulses[i]["time"].as<size_t>() - 1; // -1 to help define 1-gen epochs for admixture

        bool valid = 0;
        for(size_t j = 1; j < timeBounds.size(); ++j)
        {
          if(time == timeBounds[j])
          {
            valid = 1;

            int row = -1;
            int col = -1;

            for(size_t k = 0; k < pops_[j - 1].size(); ++ k)
            {
              // ancestral populations are found in the previous epoch
              if(pops_[j - 1][k]->getName() == dest)
                row = k;

              else if(pops_[j - 1][k]->getName() == source)
                col = k;
            }

            if(row == -1 || col == -1)
              throw bpp::Exception("Demes::could not find admixing populations " + source + " & " + dest + " in epoch " + bpp::TextTools::toString(j - 1) + "!");

            pulses_[j - 1](row, col) = f;  // "from" (f), "to" (1-f, ommited)
          }
        }

        if(!valid)
          throw bpp::Exception("Demes::time of admixture pulse must match the start of an epoch!");
      }
    } // exits 'pulses' field of Demes file

    // move to custom fields for moments++ (mutation, recombination, selection)
    // NOTE: this field must be placed last in the demes file
    else if(it->first.as<std::string>() == "metadata")
    {
      YAML::Node meta = it->second;

      for(size_t i = 0; i < meta.size(); ++i)
      {
        if(meta["mutation"])
        {
          YAML::Node muts = meta["mutation"];

          double rate = 1e-8;
          size_t startTime = timeBounds.front();
          size_t endTime = 0;

          for(size_t j = 0; j < muts.size(); ++j)
          {
            if(muts["rate"])
              rate = muts["rate"].as<double>();

            if(muts["start_time"])
              startTime = muts["start_time"].as<size_t>();

            if(muts["end_time"])
              endTime = muts["end_time"].as<size_t>();

            bool match = 1;
            for(size_t k = 1; k < timeBounds.size(); ++k)
            {
              if((startTime == timeBounds[k - 1] && endTime == timeBounds[k]) || (startTime == timeBounds.front() && endTime == 0))
                mutRates_[k - 1] = rate;

              else
                match = 0;
            }

            if(!match)
              throw bpp::Exception("Demes::start_time and end_time of 'mutation' do not match the span of any epoch!");
          }
        }

        else if(meta["recombination"])
        {
          YAML::Node recs = meta["recombination"];

          double rate = 1e-8;
          size_t startTime = timeBounds.front();
          size_t endTime = 0;

          for(size_t j = 0; j < recs.size(); ++j)
          {
            if(recs["rate"])
              rate = recs["rate"].as<double>();

            if(recs["start_time"])
              startTime = recs["start_time"].as<size_t>();

            if(recs["end_time"])
              endTime = recs["end_time"].as<size_t>();

            bool match = 1;
            for(size_t k = 1; k < timeBounds.size(); ++k)
            {
              if((startTime == timeBounds[k - 1] && endTime == timeBounds[k]) || (startTime == timeBounds.front() && endTime == 0))
                recRates_[k - 1] = rate;

              else
                match = 0;
            }

            if(!match)
              throw bpp::Exception("Demes::start_time and end_time of 'recombination' do not match the span of any epoch!");
          }
        }

        else if(meta["selection"])
        {
          YAML::Node sel = meta["selection"];

          std::cout << "TODO: list demes which experience selection on the left locus\n"; // then pop->setSelectiveConstraint(bool)
          double s = 1e-3;
          size_t startTime = timeBounds.front();
          size_t endTime = 0;

          for(size_t j = 0; j < sel.size(); ++j) // sel period by sel period
          {
            if(sel["s"])
              s = sel["s"].as<double>();

            if(sel[j]["start_time"])
              startTime = sel["start_time"].as<size_t>();

            if(sel[j]["end_time"])
              endTime = sel["end_time"].as<size_t>();

            bool match = 1;
            for(size_t k = 1; k < timeBounds.size(); ++k)
            {
              if((startTime == timeBounds[k - 1] && endTime == timeBounds[k]) || (startTime == timeBounds.front() && endTime == 0))
                selCoeffs_[k - 1] = s;

              else
                match = 0;
            }

            if(!match)
              throw bpp::Exception("Demes::start_time and end_time of 'selection' do not match the span of any epoch!");
          }
        }
      }
    } // exits 'metadata' field of Demes file
  }
}
