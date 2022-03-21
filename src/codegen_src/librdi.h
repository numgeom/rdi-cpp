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
#include "rdi_params.hpp"
#include <cstddef>
#include <cstdlib>

// Type Declarations
namespace rdi_kernel {
struct RdiMesh;

}

// Function Declarations
namespace rdi_kernel {

static inline void rdi_assemble_osusop(
    const RdiMesh *mesh, const ::coder::array<int, 2U> &stcls,
    const RdiParams *params, ::coder::array<int, 1U> &rowPtr,
    ::coder::array<int, 1U> &colInd, ::coder::array<double, 1U> &vals,
    ::coder::array<int, 1U> &nnzPr, ::coder::array<int, 1U> &rdNodes);

static inline void rdi_assemble_osusop2(
    const RdiMesh *mesh, const ::coder::array<int, 2U> &stcls,
    const RdiParams *params, ::coder::array<int, 1U> &rowPtr,
    ::coder::array<int, 1U> &colInd, ::coder::array<double, 1U> &vals,
    ::coder::array<int, 1U> &nnzPr, ::coder::array<int, 1U> &rdNodes);

static inline void rdi_build_node2cell(int n, const ::coder::array<int, 2U> &conn,
                                 ::coder::array<int, 1U> &n2cPtr,
                                 ::coder::array<int, 1U> &n2cList);

static inline void rdi_build_node2node(int n, const ::coder::array<int, 2U> &conn,
                                 int dim, const ::coder::array<int, 1U> &n2cPtr,
                                 const ::coder::array<int, 1U> &n2cList,
                                 ::coder::array<int, 1U> &n2nPtr,
                                 ::coder::array<int, 1U> &n2nList);

static inline void rdi_compute_cellsizes(const ::coder::array<double, 2U> &xs,
                                   const ::coder::array<int, 2U> &conn,
                                   const ::coder::array<int, 1U> &n2nPtr,
                                   const ::coder::array<int, 1U> &n2nList,
                                   const RdiParams *params,
                                   const ::coder::array<double, 2U> &nrms,
                                   ::coder::array<double, 1U> &h);

static inline void rdi_compute_cellsizes2(const ::coder::array<double, 2U> &xs,
                                   const ::coder::array<int, 2U> &conn,
                                   const ::coder::array<int, 1U> &n2nPtr,
                                   const ::coder::array<int, 1U> &n2nList,
                                   const RdiParams *params,
                                   const ::coder::array<double, 2U> &nrms,
                                   ::coder::array<double, 1U> &h);

static inline void rdi_compute_cellsizes(const ::coder::array<double, 2U> &xs,
                                   const ::coder::array<int, 2U> &conn,
                                   const ::coder::array<int, 1U> &n2nPtr,
                                   const ::coder::array<int, 1U> &n2nList,
                                   const RdiParams *params,
                                   ::coder::array<double, 1U> &h);

static inline void rdi_compute_cellweights(const ::coder::array<double, 2U> &xs,
                                     const ::coder::array<int, 2U> &conn,
                                     const RdiParams *params,
                                     ::coder::array<double, 1U> &w);

static inline void rdi_compute_cellweights2(const ::coder::array<double, 2U> &xs,
                                     const ::coder::array<int, 2U> &conn,
                                     const RdiParams *params,
                                     ::coder::array<double, 1U> &w);

static inline void rdi_compute_inds(const ::coder::array<int, 1U> &rowPtr,
                              const ::coder::array<int, 1U> &colInd,
                              const ::coder::array<double, 1U> &vals,
                              const ::coder::array<double, 2U> &fs,
                              const ::coder::array<double, 2U> &dfGlobal,
                              const RdiMesh *mesh, const RdiParams *params,
                              ::coder::array<double, 2U> &alphaCell,
                              ::coder::array<double, 2U> &alphaNode,
                              ::coder::array<double, 2U> &beta);

static inline void rdi_compute_inds2(const ::coder::array<int, 1U> &rowPtr,
                              const ::coder::array<int, 1U> &colInd,
                              const ::coder::array<double, 1U> &vals,
                              const ::coder::array<double, 2U> &fs,
                              const ::coder::array<double, 2U> &dfGlobal,
                              const RdiMesh *mesh, const RdiParams *params,
                              ::coder::array<double, 2U> &alphaCell,
                              ::coder::array<double, 2U> &alphaNode,
                              ::coder::array<double, 2U> &beta);

static inline void rdi_compute_oscind(const ::coder::array<double, 2U> &dfGlobal,
                                const ::coder::array<double, 2U> &alphaCell,
                                const RdiMesh *mesh, const RdiParams *params,
                                ::coder::array<double, 2U> &beta);

static inline void rdi_compute_oscind2(const ::coder::array<double, 2U> &dfGlobal,
                                const ::coder::array<double, 2U> &alphaCell,
                                const RdiMesh *mesh, const RdiParams *params,
                                ::coder::array<double, 2U> &beta);

static inline void rdi_compute_osusind(const ::coder::array<int, 1U> &rowPtr,
                                 const ::coder::array<int, 1U> &colInd,
                                 const ::coder::array<double, 1U> &vals,
                                 const ::coder::array<double, 2U> &fs,
                                 const RdiMesh *mesh, const RdiParams *params,
                                 ::coder::array<double, 2U> &alphaCell,
                                 ::coder::array<double, 2U> &alphaNode);

static inline void rdi_compute_osusind2(const ::coder::array<int, 1U> &rowPtr,
                                 const ::coder::array<int, 1U> &colInd,
                                 const ::coder::array<double, 1U> &vals,
                                 const ::coder::array<double, 2U> &fs,
                                 const RdiMesh *mesh, const RdiParams *params,
                                 ::coder::array<double, 2U> &alphaCell,
                                 ::coder::array<double, 2U> &alphaNode);

static inline void rdi_default_params(int dim, int nThreads, RdiParams *params);

static inline void rdi_default_params2(int dim, int nThreads, RdiParams *params);

static inline void rdi_default_params(int dim, RdiParams *params);

static inline void
rdi_mark_discontinuities(const ::coder::array<double, 2U> &fs,
                          const ::coder::array<double, 2U> &alphaCell,
                          const ::coder::array<double, 2U> &beta,
                          const ::coder::array<double, 2U> &dfGlobal,
                          const RdiMesh *mesh, const RdiParams *params,
                          ::coder::array<signed char, 2U> &disTags);

static inline void rdi_partition(const RdiParams *params, RdiMesh *mesh);

} // namespace rdi_kernel

#include "librdi.cpp"
#endif
// End of code generation (librdi.h)
