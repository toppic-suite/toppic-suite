#include "base/file_util.hpp"
#include "console/summary.hpp"

namespace prot {

void Summary::outputSummary(Argument arguments, 
                            boost::posix_time::ptime start_time,
                            boost::posix_time::ptime stop_time) {
  
  std::string spec_file_name = arguments.getArguments()["spectrumFileName"];
  std::string base_name = FileUtil::basename(spec_file_name);
  std::string output_file_name = base_name + "." + "SUMMARY";

  std::ofstream output; 
  output.open(output_file_name, std::ios::out | std::ios::app);

  boost::posix_time::time_duration duration;
  duration = stop_time - start_time;
  output << "1. Time" << std::endl;
  output << "Start time: " << boost::posix_time::to_simple_string(start_time) << std::endl;
  output << "Stop time: " << boost::posix_time::to_simple_string(stop_time) << std::endl;
  output << "Running time: " << boost::posix_time::to_simple_string(duration) << std::endl;


  Argument::outputArguments(output, arguments.getArguments());
  output.close();
}

} /* namespace prot */
