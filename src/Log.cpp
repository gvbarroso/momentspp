/*
 * Authors: Gustavo V. Barroso 
 * Created: 21/10/2022
 * Last modified: 21/10/2022
 *
 */


#include "Log.hpp"

void Log::stop_timer(double mult, const std::string& task, const std::string& unit)
{
  endTimePoint_ = std::chrono::high_resolution_clock::now();

  auto start = std::chrono::time_point_cast<std::chrono::microseconds>(startTimePoint_).time_since_epoch().count();
  auto end = std::chrono::time_point_cast<std::chrono::microseconds>(endTimePoint_).time_since_epoch().count();

  auto duration = end - start;
  double conv = duration / mult;
    
  logFile_ << "Duration of " << task << " = " << std::setprecision(3) << conv << " (" << unit << ")\n";
}

