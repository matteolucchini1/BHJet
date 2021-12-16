// Provide the Python interface to the BHJet model

#include <iostream>

// Allow the model to be called from Python the same way
// XSPEC models are.
//
#include "sherpa/astro/xspec_extension.hh"

extern "C" {

  void xspec_jetinterp(const double *ear, int ne, const double *param, int spectrumNumber,
		       double *photeng, double *photspec,
		       const char* initStr);
}

// Wrapper

static PyMethodDef Wrappers[] = {
  XSPECMODELFCT_C_NORM(xspec_jetinterp, 28),
  { NULL, NULL, 0, NULL }
};

// Create the Python module

static struct PyModuleDef wrapper_module = {
  PyModuleDef_HEAD_INIT,
  "_models",
  NULL,
  -1,
  Wrappers,
};

// Technically we should ensure that the Sherpa XSPEC module is
// initialized. However, we should not need it (famous last
// words!).
//
PyMODINIT_FUNC PyInit__models(void) {
  import_array();
  return PyModule_Create(&wrapper_module);
}
