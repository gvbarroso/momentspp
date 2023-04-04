/*
 * Authors: Gustavo V. Barroso 
 * Created: 21/10/2022
 * Last modified: 04/04/2023
 *
 */

#ifndef _LOG_
#define _LOG_

#include <iomanip>
#include <chrono>
#include <string>
#include <iostream>

class Log
{

private:
  std::chrono::time_point<std::chrono::high_resolution_clock> startTimePoint_;
  std::chrono::time_point<std::chrono::high_resolution_clock> endTimePoint_;

public:
  Log():
  startTimePoint_(),
  endTimePoint_()
  { }

  void stop_timer(std::ostream& stream);

  void start_timer()
  {
    startTimePoint_ = std::chrono::high_resolution_clock::now();
  }
  
};

#endif
