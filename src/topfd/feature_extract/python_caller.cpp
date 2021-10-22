//
// Created by abbash on 10/11/2021.
//
//#include <cmath>
//#include <Python.h>
#include "python_caller.h"

namespace toppic {

PythonCaller::PythonCaller() {
  Py_Initialize();
  sysmodule_ = PyImport_ImportModule("sys");
  syspath_ = PyObject_GetAttrString(sysmodule_, "path");
  PyList_Append(syspath_, PyUnicode_FromString("./toppic_resources"));

  mymodule_ = PyImport_ImportModule("script");
  myfunc_ = PyObject_GetAttrString(mymodule_, "fit_vogit");
//  PyRun_SimpleString("from lmfit.models import SkewedVoigtModel");
//  PyRun_SimpleString("import numpy");

}

std::vector<double> PythonCaller::process(std::vector<double> xic, int ref_scan) {
  PyObject *pXVec = PyTuple_New(xic.size());
  for (size_t i = 0; i < xic.size(); ++i)
    PyTuple_SetItem(pXVec, i, PyFloat_FromDouble(xic[i]));

//  LOG_ERROR("Tuple packed " << ref_scan << " - " << xic[ref_scan]);
  PyObject *pArgs = PyTuple_Pack(2, PyLong_FromLong(ref_scan), pXVec);
  PyObject *pValue = PyObject_CallObject(myfunc_, pArgs);
//  LOG_ERROR("Python Returned!!!");
  std::vector<double> fitted_arr;
  for (int i = 0; i < PyList_GET_SIZE(pValue); i++)
    fitted_arr.push_back(PyFloat_AsDouble(PyList_GetItem(pValue, i)));
  return fitted_arr;
}

}