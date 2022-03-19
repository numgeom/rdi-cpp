// Copyright 2022 The NumGeom Group, Stony Brook University
//
// librdi_types.h
//
// Code generation for function 'rdi_default_params'
//

#ifndef LIBRDI_TYPES_H
#define LIBRDI_TYPES_H

// Include files
#include "rtwtypes.h"
#include "coder_array.h"

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

struct Omp4mPart {
  int nparts;
  ::coder::array<int, 1U> part_ptr;
  ::coder::array<int, 1U> part_list;
  ::coder::array<int, 1U> shared_ents;
};

struct RdiMesh {
  ::coder::array<double, 2U> xs;
  ::coder::array<int, 2U> conn;
  int topoDim;
  int surfType;
  ::coder::array<double, 2U> dirs;
  ::coder::array<signed char, 1U> features;
  ::coder::array<double, 1U> cellSizes;
  double hGlobal;
  ::coder::array<double, 1U> cellWeights;
  ::coder::array<int, 1U> n2cPtr;
  ::coder::array<int, 1U> n2cList;
  ::coder::array<int, 1U> n2nPtr;
  ::coder::array<int, 1U> n2nList;
  ::coder::array<Omp4mPart, 1U> parts;
};

} // namespace rdi_kernel

#endif
// End of code generation (librdi_types.h)
