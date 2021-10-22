//
// Created by abbash on 10/11/2021.
//

#include <cmath>
#include <Python.h>
#include <vector>
#include <memory>
#include "common/util/logger.hpp"

#ifndef TOPPIC_PYTHON_CALLER_H
#define TOPPIC_PYTHON_CALLER_H

namespace toppic {

class PythonCaller {
 public:
  PythonCaller();
  void close_connection(){ Py_Finalize(); };
  std::vector<double> process(std::vector<double> xic, int ref_scan);

 private:
  PyObject *sysmodule_;
  PyObject *syspath_;
  PyObject *mymodule_;
  PyObject *myfunc_;
};

typedef std::shared_ptr<PythonCaller> PythonCallerPtr;

}


#endif //TOPPIC_PYTHON_CALLER_H
