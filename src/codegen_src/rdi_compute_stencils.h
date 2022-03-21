// Copyright 2022 The NumGeom Group, Stony Brook University
//
// rdi_compute_stencils.h
//
// Code generation for function 'rdi_compute_stencils'
//

#ifndef RDI_COMPUTE_STENCILS_H
#define RDI_COMPUTE_STENCILS_H

// Include files
#include "rtwtypes.h"
#include "m2c_lib.h"
#include "rdi_compute_stencils_types.h"
#include "coder_array.h"
#include "rdi_params.hpp"
#include <cstddef>
#include <cstdlib>

// Function Declarations
namespace rdi_stencils {
static inline void rdi_compute_stencils(int n, const ::coder::array<int, 2U> &conn,
                                  const RdiParams *params,
                                  const ::coder::array<int, 1U> &nrange,
                                  ::coder::array<int, 2U> &stcls);

static inline void rdi_compute_stencils(int n, const ::coder::array<int, 2U> &conn,
                                  const RdiParams *params,
                                  ::coder::array<int, 2U> &stcls);

} // namespace rdi_stencils

#include "rdi_compute_stencils.cpp"
#endif
// End of code generation (rdi_compute_stencils.h)
