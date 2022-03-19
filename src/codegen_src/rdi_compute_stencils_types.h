// Copyright 2022 The NumGeom Group, Stony Brook University
//
// rdi_compute_stencils_types.h
//
// Code generation for function 'rdi_compute_stencils'
//

#ifndef RDI_COMPUTE_STENCILS_TYPES_H
#define RDI_COMPUTE_STENCILS_TYPES_H

// Include files
#include "rtwtypes.h"

// Type Definitions
namespace rdi_kernel {
struct RdiParams {
  int dim;
  double ring;
  int maxStclSize;
  int verbose;
  int degree;
  int surfType;
  double epsBeta;
  double cGlobal;
  double cLocal;
  double kappa1;
  double kappa0;
  boolean_T markNearDis;
  boolean_T wlsInterp0;
  boolean_T wlsUseDag;
  int nThreads;
  boolean_T parTask;
};

} // namespace rdi_kernel

#endif
// End of code generation (rdi_compute_stencils_types.h)
