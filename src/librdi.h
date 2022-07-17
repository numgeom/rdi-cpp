// Copyright 2022 The NumGeom Group, Stony Brook University
//
// librdi.h
//
// Code generation for function 'librdi'
//

#ifndef LIBRDI_H
#define LIBRDI_H

// Include files
#include "rtwtypes.h"
#include "m2c_lib.h"
#include "librdi_types.h"
#include "coder_array.h"
#include <cstddef>
#include <cstdlib>

// Type Declarations
namespace rdi_kernel {
struct RdiObject;

struct WlsMesh;

} // namespace rdi_kernel


// Function Declarations
namespace rdi_kernel {

static inline void rdi_apply(const RdiObject *rdi, const WlsMesh *mesh,
                       const ::coder::array<real_T, 2U> &fs,
                       const ::coder::array<real_T, 2U> &df,
                       ::coder::array<real_T, 2U> &alpha,
                       ::coder::array<real_T, 2U> &beta,
                       ::coder::array<int8_T, 2U> &distags,
                       ::coder::array<real_T, 2U> &alphanode);

static inline void rdi_apply(const RdiObject *rdi, const WlsMesh *mesh,
                       const ::coder::array<real_T, 2U> &fs,
                       const ::coder::array<real_T, 2U> &df, ::coder::SizeType varargin_2,
                       ::coder::array<real_T, 2U> &alpha,
                       ::coder::array<real_T, 2U> &beta,
                       ::coder::array<int8_T, 2U> &distags,
                       ::coder::array<real_T, 2U> &alphanode);

static inline void rdi_apply3(const RdiObject *rdi, const WlsMesh *mesh,
                       const ::coder::array<real_T, 2U> &fs,
                       const ::coder::array<real_T, 2U> &df,
                       ::coder::array<real_T, 2U> &alpha,
                       ::coder::array<real_T, 2U> &beta,
                       ::coder::array<int8_T, 2U> &distags,
                       ::coder::array<real_T, 2U> &alphanode);

static inline void rdi_compute_osusop(RdiObject *rdi, WlsMesh *mesh);

static inline void rdi_compute_osusop2(RdiObject *rdi, WlsMesh *mesh,
                                boolean_T interp0);

static inline void rdi_compute_osusop(RdiObject *rdi, WlsMesh *mesh,
                                boolean_T interp0, ::coder::SizeType varargin_2);

static inline void rdi_compute_osusop4(RdiObject *rdi, WlsMesh *mesh,
                                boolean_T interp0);

static inline void rdi_create(WlsMesh *mesh, RdiObject *rdi);

static inline void rdi_create(WlsMesh *mesh, ::coder::SizeType degree, RdiObject *rdi);

static inline void rdi_create(WlsMesh *mesh, ::coder::SizeType degree, ::coder::SizeType nparts,
                        RdiObject *rdi);

static inline void rdi_create(WlsMesh *mesh, ::coder::SizeType degree, ::coder::SizeType nparts,
                        ::coder::SizeType nrmid, RdiObject *rdi);

static inline void rdi_create5(WlsMesh *mesh, ::coder::SizeType degree, ::coder::SizeType nparts,
                        ::coder::SizeType nrmid, real_T ring, RdiObject *rdi);

static inline void rdi_create(WlsMesh *mesh, ::coder::SizeType degree, ::coder::SizeType nparts,
                        ::coder::SizeType nrmid, real_T ring, ::coder::SizeType varargin_2,
                        RdiObject *rdi);

static inline void rdi_create7(WlsMesh *mesh, ::coder::SizeType degree, ::coder::SizeType nparts,
                        ::coder::SizeType nrmid, real_T ring, RdiObject *rdi);

static inline void rdi_postproc_markers(::coder::array<int8_T, 2U> &distags,
                                  const WlsMesh *mesh);

static inline void rdi_postproc_markers(::coder::array<int8_T, 2U> &distags,
                                  const WlsMesh *mesh, ::coder::SizeType nlayers);

} // namespace rdi_kernel

#include "librdi.cpp"
#endif
// End of code generation (librdi.h)
