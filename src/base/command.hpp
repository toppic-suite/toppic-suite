/*
 * command.hpp
 *
 *  Created on: Apr 24, 2014
 *      Author: xunlikun
 */

#ifndef COMMAND_HPP_
#define COMMAND_HPP_

#include <map>
#include <fstream>
#include <algorithm>
#include <memory>

#include "base/logger.hpp"
#include "base/xml_dom_document.hpp"
#include "base/string_util.hpp"

namespace prot {

class command {
 public:
  command();
  virtual ~command();

  int setCommandValue(std::string command,std::string value);
  bool validateCommandValue(std::string command,std::string value);
  std::map<std::string,std::string> getArgument(){return arguments_;}
  void setArgumentsWithConfigFile(std::string filename);

 private:
  void initOption();
  void initArguments();

  std::map<std::string,std::string> options_;
  std::map<std::string,std::string> arguments_;
};

typedef std::shared_ptr<command> cmdPtr;

void showCommondList();
bool existFile(std::string filename);
std::string getExePath(std::string & command_path);

int getOS();
int runCommand(std::string cmd);
int runCommand(std::string cmd,std::string mod);


class SystemInfo{
 public:
  static void initSystemInfo(std::string argument_zero);
  static std::string getExeFilePath(){return exe_path_;};
  static std::string getOSString(){return os_;};
  static void setXmlPath(std::string xml_path){xml_path_=xml_path;};
  static void setHtmlPath(std::string html_path){html_path_=html_path;};
  static std::string getXmlPath(){return xml_path_;};
  static std::string getHtmlPath(){return html_path_;};
 private:
  static std::string exe_path_;
  static std::string xml_path_;
  static std::string html_path_;
  static std::string os_;
};
} /* namespace prot */

#endif /* COMMAND_HPP_ */
