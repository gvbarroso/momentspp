/*
 * Authors: Gustavo V. Barroso 
 * Created: 21/10/2022
 * Last modified: 21/10/2022
 *
 */

#ifndef _LOG_
#define _LOG_

#include <iomanip>
#include <chrono>
#include <string>
#include <iostream>
#include <fstream>


class Log {

private:
  std::chrono::time_point<std::chrono::high_resolution_clock> startTimePoint_;
  std::chrono::time_point<std::chrono::high_resolution_clock> endTimePoint_;

  std::ofstream logFile_;

public:
  Log():
  startTimePoint_(),
  endTimePoint_(),
  logFile_()
  { }

  ~Log()
  {
    if(logFile_.is_open())
      logFile_.close();
  }

  void stop_timer(double mult, const std::string& task, const std::string& unit);

  void start_timer()
  {
    startTimePoint_ = std::chrono::high_resolution_clock::now();
  }

  void openFile(const std::string& name)
  {
    if(!logFile_.is_open())
      logFile_.open(name);
  }

  void closeFile()
  {
    if(logFile_.is_open())
      logFile_.close();
  }

  std::ofstream& getLogFile()
  {
    return logFile_;
  }
  
};

#endif
