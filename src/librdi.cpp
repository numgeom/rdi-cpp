// Copyright 2022 The NumGeom Group, Stony Brook University
//
// librdi.cpp
//
// Code generation for function 'librdi'
//

// Include files
#include "librdi.h"
#include "m2c_lib.h"
#include "librdi_types.h"
#include "coder_array.h"
#include "m2c_lib.h"
#ifdef OMP4M_HAS_METIS
#include "metis.h"
#endif
#include "sfe_quadrules_rowmajor.h"
#include "sfe_shapefuncs_rowmajor.h"
#include "wls_lapack.hpp"
#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstring>
#ifdef _OPENMP
#include <omp.h>
#endif
#include <stdexcept>
#include <stdio.h>

// Type Definitions
namespace rdi_kernel {
struct SfeObject {
  int32_T etypes[2];
  int32_T nnodes[2];
  int32_T geom_dim;
  int32_T topo_dim;
  int8_T facetid;
  int32_T nqp;
  ::coder::array<real_T, 1U> ws;
  ::coder::array<real_T, 2U> cs;
  ::coder::array<real_T, 2U> shapes_sol;
  ::coder::array<real_T, 2U> shapes_geom;
  ::coder::array<real_T, 3U> derivs_sol;
  ::coder::array<real_T, 3U> derivs_geom;
  ::coder::array<real_T, 2U> cs_phy;
  ::coder::array<real_T, 2U> grads_sol;
  ::coder::array<real_T, 2U> grads_geom;
  ::coder::array<real_T, 2U> jacTs;
  ::coder::array<real_T, 1U> wdetJ;
  ::coder::array<real_T, 2U> dwork1;
  ::coder::array<real_T, 2U> dwork2;
  ::coder::array<real_T, 2U> xswork;
  ::coder::array<int32_T, 2U> iwork;
};
struct b_WlsObject {
  int32_T nstpnts;
  int32_T degree;
  int32_T order;
  boolean_T unimono;
  int32_T interp0;
  int32_T stride;
  ::coder::array<real_T, 2U> us;
  ::coder::bounded_array<real_T, 3U, 2U> origin;
  ::coder::array<real_T, 1U> rweights;
  ::coder::bounded_array<real_T, 3U, 2U> hs_inv;
  ::coder::array<real_T, 2U> V;
  ::coder::array<real_T, 2U> QR;
  ::coder::array<real_T, 2U> rhs;
  int32_T nevpnts;
  int32_T nrows;
  int32_T ncols;
  int32_T rank;
  boolean_T fullrank;
  ::coder::array<int32_T, 1U> jpvt;
  ::coder::array<real_T, 1U> work;
  boolean_T rowmajor;
  ::coder::array<real_T, 2U> QRt;
  ::coder::bounded_array<real_T, 4U, 1U> runtimes;
};
struct c_WlsObject {
  int32_T nstpnts;
  int32_T degree;
  int32_T order;
  boolean_T unimono;
  int32_T interp0;
  int32_T stride;
  ::coder::array<real_T, 2U> us;
  ::coder::bounded_array<real_T, 3U, 2U> origin;
  ::coder::array<real_T, 1U> rweights;
  ::coder::bounded_array<real_T, 3U, 2U> hs_inv;
  ::coder::array<real_T, 2U> V;
  ::coder::array<real_T, 2U> QR;
  ::coder::array<real_T, 2U> rhs;
  int32_T nevpnts;
  int32_T nrows;
  int32_T ncols;
  int32_T rank;
  boolean_T fullrank;
  ::coder::array<int32_T, 1U> jpvt;
  ::coder::array<real_T, 1U> work;
  boolean_T rowmajor;
  ::coder::array<real_T, 2U> QRt;
  ::coder::bounded_array<real_T, 4U, 1U> runtimes;
};

struct d_WlsObject {
  int32_T nstpnts;
  int32_T degree;
  int32_T order;
  boolean_T unimono;
  int32_T interp0;
  int32_T stride;
  ::coder::array<real_T, 2U> us;
  ::coder::bounded_array<real_T, 3U, 2U> origin;
  ::coder::array<real_T, 1U> rweights;
  ::coder::array<real_T, 2U> hs_inv;
  ::coder::array<real_T, 2U> V;
  ::coder::array<real_T, 2U> QR;
  ::coder::array<real_T, 2U> rhs;
  int32_T nevpnts;
  int32_T nrows;
  int32_T ncols;
  int32_T rank;
  boolean_T fullrank;
  ::coder::array<int32_T, 1U> jpvt;
  ::coder::array<real_T, 1U> work;
  boolean_T rowmajor;
  ::coder::array<real_T, 2U> QRt;
  ::coder::bounded_array<real_T, 4U, 1U> runtimes;
};

struct e_WlsObject {
  int32_T nstpnts;
  int32_T degree;
  int32_T order;
  boolean_T unimono;
  int32_T interp0;
  int32_T stride;
  ::coder::array<real_T, 2U> us;
  ::coder::bounded_array<real_T, 3U, 2U> origin;
  ::coder::array<real_T, 1U> rweights;
  ::coder::array<real_T, 2U> hs_inv;
  ::coder::array<real_T, 2U> V;
  ::coder::array<real_T, 2U> QR;
  ::coder::array<real_T, 2U> rhs;
  int32_T nevpnts;
  int32_T nrows;
  int32_T ncols;
  int32_T rank;
  boolean_T fullrank;
  ::coder::array<int32_T, 1U> jpvt;
  ::coder::array<real_T, 1U> work;
  boolean_T rowmajor;
  ::coder::array<real_T, 2U> QRt;
  ::coder::bounded_array<real_T, 4U, 1U> runtimes;
};

} // namespace rdi_kernel

// Variable Definitions
namespace rdi_kernel {
static const int16_T iv[250]{
    1,  0,   0,   0,  0,  0,   0,   0,   0,   0,   0,   0,   0,  0,  0,  0,  0,
    0,  0,   0,   0,  0,  0,   0,   0,   0,   0,   0,   0,   0,  0,  0,  0,  0,
    0,  2,   0,   0,  0,  3,   0,   0,   0,   4,   4,   0,   0,  5,  5,  0,  0,
    6,  6,   0,   0,  7,  7,   0,   0,   0,   0,   0,   0,   0,  0,  0,  0,  3,
    0,  0,   0,   6,  0,  0,   0,   10,  10,  0,   0,   15,  15, 15, 0,  21, 21,
    21, 0,   28,  28, 28, 0,   0,   0,   0,   0,   0,   0,   0,  0,  4,  0,  0,
    0,  9,   0,   0,  0,  16,  16,  0,   0,   25,  25,  0,   0,  36, 36, 0,  0,
    49, 49,  0,   0,  0,  0,   0,   0,   0,   0,   0,   0,   4,  0,  0,  0,  10,
    0,  0,   0,   20, 20, 0,   0,   35,  35,  35,  0,   56,  56, 56, 0,  84, 84,
    84, 0,   0,   0,  0,  0,   0,   0,   0,   0,   5,   0,   0,  0,  14, 0,  0,
    0,  30,  30,  0,  0,  55,  55,  55,  0,   91,  0,   0,   0,  0,  0,  0,  0,
    0,  0,   0,   0,  0,  0,   0,   0,   6,   0,   0,   0,   18, 0,  0,  0,  40,
    40, 0,   0,   75, 75, 75,  0,   126, 126, 0,   0,   196, 0,  0,  0,  0,  0,
    0,  0,   0,   0,  0,  0,   8,   0,   0,   0,   27,  0,   0,  0,  64, 64, 0,
    0,  125, 125, 0,  0,  216, 216, 0,   0,   343, 343, -1};

} // namespace rdi_kernel

// Function Declarations
namespace rdi_kernel {
static inline
::coder::SizeType append_wlsmesh_kring(WlsMesh *mesh);

static inline
void
assemble_body(RdiObject *rdi, const ::coder::array<real_T, 2U> &mesh_coords,
              const ::coder::array<ConnData, 1U> &mesh_elemtables,
              const ::coder::array<uint64_T, 1U> &mesh_teids,
              const ::coder::array<int64_T, 1U> &mesh_node2elems_row_ptr,
              const ::coder::array<int32_T, 1U> &mesh_node2elems_col_ind,
              const ::coder::array<Stencils, 1U> &mesh_stencils, ::coder::SizeType stclid,
              boolean_T interp0);

static inline
void
assemble_body_range(const ::coder::array<int64_T, 1U> &row_ptr,
                    const ::coder::array<int32_T, 1U> &col_ind,
                    const ::coder::array<real_T, 2U> &mesh_coords,
                    const ::coder::array<ConnData, 1U> &mesh_elemtables,
                    const ::coder::array<uint64_T, 1U> &mesh_teids,
                    const ::coder::array<int64_T, 1U> &stcl_row_ptr,
                    const ::coder::array<int32_T, 1U> &stcl_col_ind,
                    const ::coder::array<int64_T, 1U> &n2e_row_ptr,
                    const ::coder::array<int32_T, 1U> &n2e_col_ind,
                    ::coder::SizeType degree, boolean_T interp0,
                    const ::coder::array<int32_T, 1U> &nrange, ::coder::SizeType istart,
                    ::coder::SizeType iend, ::coder::array<real_T, 1U> &val,
                    ::coder::array<boolean_T, 1U> &rdtags, boolean_T *fullrank);

static inline
void assemble_body_task(
    RdiObject *rdi, const ::coder::array<real_T, 2U> &mesh_coords,
    const ::coder::array<ConnData, 1U> &mesh_elemtables,
    const ::coder::array<uint64_T, 1U> &mesh_teids,
    const ::coder::array<int64_T, 1U> &mesh_node2elems_row_ptr,
    const ::coder::array<int32_T, 1U> &mesh_node2elems_col_ind,
    const ::coder::array<Omp4mPart, 1U> &mesh_nodeparts,
    const ::coder::array<int64_T, 1U> &stcl_row_ptr,
    const ::coder::array<int32_T, 1U> &stcl_col_ind, boolean_T interp0);

static inline
void
assemble_surf(RdiObject *rdi, const ::coder::array<real_T, 2U> &mesh_coords,
              const ::coder::array<ConnData, 1U> &mesh_elemtables,
              const ::coder::array<uint64_T, 1U> &mesh_teids,
              const ::coder::array<NormalsData, 1U> &mesh_nrmstables,
              const ::coder::array<int64_T, 1U> &mesh_node2elems_row_ptr,
              const ::coder::array<int32_T, 1U> &mesh_node2elems_col_ind,
              const ::coder::array<Stencils, 1U> &mesh_stencils, ::coder::SizeType stclid,
              boolean_T interp0);

static inline
void
assemble_surf_range(const ::coder::array<int64_T, 1U> &row_ptr,
                    const ::coder::array<int32_T, 1U> &col_ind,
                    const ::coder::array<real_T, 2U> &mesh_coords,
                    const ::coder::array<ConnData, 1U> &mesh_elemtables,
                    const ::coder::array<uint64_T, 1U> &mesh_teids,
                    const ::coder::array<int64_T, 1U> &stcl_row_ptr,
                    const ::coder::array<int32_T, 1U> &stcl_col_ind,
                    const ::coder::array<int64_T, 1U> &n2e_row_ptr,
                    const ::coder::array<int32_T, 1U> &n2e_col_ind,
                    const ::coder::array<real_T, 2U> &nrms, ::coder::SizeType degree,
                    boolean_T interp0,
                    const ::coder::array<int32_T, 1U> &nrange, ::coder::SizeType istart,
                    ::coder::SizeType iend, ::coder::array<real_T, 1U> &val,
                    ::coder::array<boolean_T, 1U> &rdtags, boolean_T *fullrank);

static inline
void
assemble_surf_task(RdiObject *rdi,
                   const ::coder::array<real_T, 2U> &mesh_coords,
                   const ::coder::array<ConnData, 1U> &mesh_elemtables,
                   const ::coder::array<uint64_T, 1U> &mesh_teids,
                   const ::coder::array<int64_T, 1U> &mesh_node2elems_row_ptr,
                   const ::coder::array<int32_T, 1U> &mesh_node2elems_col_ind,
                   const ::coder::array<Omp4mPart, 1U> &mesh_nodeparts,
                   const ::coder::array<int64_T, 1U> &stcl_row_ptr,
                   const ::coder::array<int32_T, 1U> &stcl_col_ind,
                   const ::coder::array<real_T, 2U> &nrms, boolean_T interp0);

static inline
::coder::SizeType b_append_wlsmesh_kring(WlsMesh *mesh);

static inline
void
b_assemble_body(RdiObject *rdi, const ::coder::array<real_T, 2U> &mesh_coords,
                const ::coder::array<ConnData, 1U> &mesh_elemtables,
                const ::coder::array<uint64_T, 1U> &mesh_teids,
                const ::coder::array<int64_T, 1U> &mesh_node2elems_row_ptr,
                const ::coder::array<int32_T, 1U> &mesh_node2elems_col_ind,
                const ::coder::array<Stencils, 1U> &mesh_stencils,
                const ::coder::array<Omp4mPart, 1U> &mesh_nodeparts,
                ::coder::SizeType stclid, boolean_T interp0, ::coder::SizeType varargin_2);

static inline
boolean_T b_assemble_body_range(
    const ::coder::array<int64_T, 1U> &row_ptr,
    const ::coder::array<int32_T, 1U> &col_ind,
    const ::coder::array<real_T, 2U> &mesh_coords,
    const ::coder::array<ConnData, 1U> &mesh_elemtables,
    const ::coder::array<uint64_T, 1U> &mesh_teids,
    const ::coder::array<int64_T, 1U> &stcl_row_ptr,
    const ::coder::array<int32_T, 1U> &stcl_col_ind,
    const ::coder::array<int64_T, 1U> &n2e_row_ptr,
    const ::coder::array<int32_T, 1U> &n2e_col_ind, ::coder::SizeType degree,
    boolean_T interp0, const ::coder::array<int32_T, 1U> &nrange, ::coder::SizeType iend,
    ::coder::array<real_T, 1U> &val, ::coder::array<boolean_T, 1U> &rdtags);

static inline
void
b_assemble_surf(RdiObject *rdi, const ::coder::array<real_T, 2U> &mesh_coords,
                const ::coder::array<ConnData, 1U> &mesh_elemtables,
                const ::coder::array<uint64_T, 1U> &mesh_teids,
                const ::coder::array<NormalsData, 1U> &mesh_nrmstables,
                const ::coder::array<int64_T, 1U> &mesh_node2elems_row_ptr,
                const ::coder::array<int32_T, 1U> &mesh_node2elems_col_ind,
                const ::coder::array<Stencils, 1U> &mesh_stencils,
                const ::coder::array<Omp4mPart, 1U> &mesh_nodeparts,
                ::coder::SizeType stclid, boolean_T interp0, ::coder::SizeType varargin_2);

static inline
boolean_T b_assemble_surf_range(
    const ::coder::array<int64_T, 1U> &row_ptr,
    const ::coder::array<int32_T, 1U> &col_ind,
    const ::coder::array<real_T, 2U> &mesh_coords,
    const ::coder::array<ConnData, 1U> &mesh_elemtables,
    const ::coder::array<uint64_T, 1U> &mesh_teids,
    const ::coder::array<int64_T, 1U> &stcl_row_ptr,
    const ::coder::array<int32_T, 1U> &stcl_col_ind,
    const ::coder::array<int64_T, 1U> &n2e_row_ptr,
    const ::coder::array<int32_T, 1U> &n2e_col_ind,
    const ::coder::array<real_T, 2U> &nrms, ::coder::SizeType degree, boolean_T interp0,
    const ::coder::array<int32_T, 1U> &nrange, ::coder::SizeType iend,
    ::coder::array<real_T, 1U> &val, ::coder::array<boolean_T, 1U> &rdtags);

static inline
void b_compute_stencils_1d(WlsMesh *mesh, ::coder::SizeType stclidx,
                                  const real_T krings[3]);

static inline
void b_compute_stencils_2d(WlsMesh *mesh, ::coder::SizeType stclidx,
                                  const real_T krings[3]);

static inline
void b_compute_stencils_kernel_1d(WlsMesh *mesh, ::coder::SizeType stclidx,
                                         real_T ring, const int32_T bounds[2]);

static inline
void b_compute_stencils_kernel_2d(WlsMesh *mesh, ::coder::SizeType stclidx,
                                         real_T ring, const int32_T bounds[2]);

static inline
void b_wlsmesh_compute_1ring(WlsMesh *mesh);

static inline
void b_wlsmesh_compute_meshprop(WlsMesh *mesh, ::coder::SizeType nrmidx);

static inline
void b_wlsmesh_compute_stencils(WlsMesh *mesh, ::coder::SizeType stclidx,
                                       real_T krings);

static inline
void bar_quadrules(::coder::SizeType degree, ::coder::array<real_T, 2U> &cs,
                          ::coder::array<real_T, 1U> &ws);

static inline
void build_part(::coder::SizeType nParts,
                       const ::coder::array<int32_T, 1U> &nparts,
                       const ::coder::array<int32_T, 1U> &cparts,
                       const ::coder::array<int32_T, 1U> &eptr,
                       const ::coder::array<int32_T, 1U> &eind,
                       ::coder::array<boolean_T, 1U> &ctags,
                       ::coder::array<int32_T, 1U> &iwork,
                       ::coder::array<int32_T, 1U> &partptr,
                       ::coder::array<int32_T, 1U> &partlist,
                       ::coder::array<int32_T, 1U> &sharedents);

static inline
::coder::SizeType c_append_wlsmesh_kring(WlsMesh *mesh);

static inline
void
c_assemble_body(RdiObject *rdi, const ::coder::array<real_T, 2U> &mesh_coords,
                const ::coder::array<ConnData, 1U> &mesh_elemtables,
                const ::coder::array<uint64_T, 1U> &mesh_teids,
                const ::coder::array<int64_T, 1U> &mesh_node2elems_row_ptr,
                const ::coder::array<int32_T, 1U> &mesh_node2elems_col_ind,
                const ::coder::array<Stencils, 1U> &mesh_stencils,
                const ::coder::array<Omp4mPart, 1U> &mesh_nodeparts,
                ::coder::SizeType stclid, boolean_T interp0);

static inline
boolean_T
c_assemble_body_range(const ::coder::array<int64_T, 1U> &row_ptr,
                      const ::coder::array<int32_T, 1U> &col_ind,
                      const ::coder::array<real_T, 2U> &mesh_coords,
                      const ::coder::array<ConnData, 1U> &mesh_elemtables,
                      const ::coder::array<uint64_T, 1U> &mesh_teids,
                      const ::coder::array<int64_T, 1U> &stcl_row_ptr,
                      const ::coder::array<int32_T, 1U> &stcl_col_ind,
                      const ::coder::array<int64_T, 1U> &n2e_row_ptr,
                      const ::coder::array<int32_T, 1U> &n2e_col_ind,
                      ::coder::SizeType degree, boolean_T interp0, ::coder::SizeType iend,
                      ::coder::array<real_T, 1U> &val,
                      ::coder::array<boolean_T, 1U> &rdtags);

static inline
void
c_assemble_surf(RdiObject *rdi, const ::coder::array<real_T, 2U> &mesh_coords,
                const ::coder::array<ConnData, 1U> &mesh_elemtables,
                const ::coder::array<uint64_T, 1U> &mesh_teids,
                const ::coder::array<NormalsData, 1U> &mesh_nrmstables,
                const ::coder::array<int64_T, 1U> &mesh_node2elems_row_ptr,
                const ::coder::array<int32_T, 1U> &mesh_node2elems_col_ind,
                const ::coder::array<Stencils, 1U> &mesh_stencils,
                const ::coder::array<Omp4mPart, 1U> &mesh_nodeparts,
                ::coder::SizeType stclid, boolean_T interp0);

static inline
boolean_T
c_assemble_surf_range(const ::coder::array<int64_T, 1U> &row_ptr,
                      const ::coder::array<int32_T, 1U> &col_ind,
                      const ::coder::array<real_T, 2U> &mesh_coords,
                      const ::coder::array<ConnData, 1U> &mesh_elemtables,
                      const ::coder::array<uint64_T, 1U> &mesh_teids,
                      const ::coder::array<int64_T, 1U> &stcl_row_ptr,
                      const ::coder::array<int32_T, 1U> &stcl_col_ind,
                      const ::coder::array<int64_T, 1U> &n2e_row_ptr,
                      const ::coder::array<int32_T, 1U> &n2e_col_ind,
                      const ::coder::array<real_T, 2U> &nrms, ::coder::SizeType degree,
                      boolean_T interp0, ::coder::SizeType iend,
                      ::coder::array<real_T, 1U> &val,
                      ::coder::array<boolean_T, 1U> &rdtags);

static inline
void call_metis_mesh(int32_T n, ::coder::array<int32_T, 1U> &eptr,
                            ::coder::array<int32_T, 1U> &eind, int32_T nParts,
                            ::coder::array<int32_T, 1U> &nparts,
                            ::coder::array<int32_T, 1U> &cparts);

namespace coder {
static inline
real_T sum(const ::coder::array<real_T, 1U> &x);

}
static inline
void
compute_beta_kernel(const ::coder::array<real_T, 2U> &mesh_coords,
                    const ::coder::array<int64_T, 1U> &mesh_node2elems_row_ptr,
                    const ::coder::array<int32_T, 1U> &mesh_node2elems_col_ind,
                    const ::coder::array<real_T, 1U> &mesh_elemmeas,
                    real_T mesh_globalh, const ::coder::array<real_T, 2U> &df,
                    const ::coder::array<real_T, 2U> &alpha, real_T epsbeta,
                    ::coder::array<real_T, 2U> &beta);

static inline
::coder::SizeType compute_connected_components(
    ::coder::array<boolean_T, 1U> &visited, ::coder::array<int32_T, 1U> &iwork,
    const ::coder::array<int64_T, 1U> &G_row_ptr,
    const ::coder::array<int32_T, 1U> &G_col_ind, ::coder::SizeType G_ncols);

static inline
void
compute_fconn_graph(::coder::array<boolean_T, 1U> &visited,
                    ::coder::array<int32_T, 1U> &iwork,
                    const ::coder::array<ConnData, 1U> &mesh_elemtables,
                    const ::coder::array<uint64_T, 1U> &mesh_teids,
                    const ::coder::array<int64_T, 1U> &mesh_node2elems_row_ptr,
                    const ::coder::array<int32_T, 1U> &mesh_node2elems_col_ind,
                    ::coder::SizeType nid, const ::coder::array<int8_T, 2U> &distags,
                    ::coder::SizeType col, ::coder::array<int64_T, 1U> &G_row_ptr,
                    ::coder::array<int32_T, 1U> &G_col_ind, int32_T *G_ncols);

static inline
void
compute_measure_kernel(const ::coder::array<real_T, 2U> &mesh_coords,
                       const ::coder::array<ConnData, 1U> &mesh_elemtables,
                       const ::coder::array<uint64_T, 1U> &mesh_teids,
                       ::coder::array<real_T, 1U> &m);

static inline
void compute_meshsizes_kernel(
    const ::coder::array<real_T, 2U> &mesh_coords,
    const ::coder::array<ConnData, 1U> &mesh_elemtables,
    const ::coder::array<uint64_T, 1U> &mesh_teids,
    const ::coder::array<int64_T, 1U> &mesh_node2nodes_row_ptr,
    const ::coder::array<int32_T, 1U> &mesh_node2nodes_col_ind,
    const ::coder::array<real_T, 2U> &nrms, ::coder::array<real_T, 1U> &elemh,
    ::coder::array<real_T, 1U> &nodeh, real_T *globalh);

static inline
void compute_meshsizes_kernel(
    const ::coder::array<real_T, 2U> &mesh_coords,
    const ::coder::array<ConnData, 1U> &mesh_elemtables,
    const ::coder::array<uint64_T, 1U> &mesh_teids,
    const ::coder::array<int64_T, 1U> &mesh_node2nodes_row_ptr,
    const ::coder::array<int32_T, 1U> &mesh_node2nodes_col_ind,
    ::coder::array<real_T, 1U> &elemh, ::coder::array<real_T, 1U> &nodeh,
    real_T *globalh);

static inline
void
compute_nodal_alpha(const ::coder::array<real_T, 2U> &alphacell,
                    const ::coder::array<real_T, 2U> &mesh_coords,
                    const ::coder::array<int64_T, 1U> &mesh_node2elems_row_ptr,
                    const ::coder::array<int32_T, 1U> &mesh_node2elems_col_ind,
                    ::coder::array<real_T, 2U> &alphanode);

static inline
void compute_stencils_1d(WlsMesh *mesh, ::coder::SizeType stclidx,
                                const real_T krings[3]);

static inline
void compute_stencils_2d(WlsMesh *mesh, ::coder::SizeType stclidx,
                                const real_T krings[3]);

static inline
void compute_stencils_kernel_1d(WlsMesh *mesh, ::coder::SizeType stclidx,
                                       real_T ring, const int32_T bounds[2]);

static inline
void compute_stencils_kernel_2d(WlsMesh *mesh, ::coder::SizeType stclidx,
                                       real_T ring, const int32_T bounds[2]);

static inline
void crsAx_kernel(const ::coder::array<int64_T, 1U> &row_ptr,
                         const ::coder::array<int32_T, 1U> &col_ind,
                         const ::coder::array<real_T, 1U> &val,
                         const ::coder::array<real_T, 2U> &x,
                         ::coder::array<real_T, 2U> &b);

static inline
void crsCompress(CrsMatrix *A, const ::coder::array<int32_T, 1U> &nnzs);

static inline
void crsProdMatVec(const ::coder::array<int64_T, 1U> &A_row_ptr,
                          const ::coder::array<int32_T, 1U> &A_col_ind,
                          const ::coder::array<real_T, 1U> &A_val,
                          const ::coder::array<real_T, 2U> &x,
                          ::coder::array<real_T, 2U> &b);

static inline
void crs_compress(::coder::array<int64_T, 1U> &row_ptr,
                         ::coder::array<int32_T, 1U> &col_ind,
                         const ::coder::array<int32_T, 1U> &nnzs);

static inline
void crs_prod_mat_vec(const ::coder::array<int64_T, 1U> &A_rowptr,
                             const ::coder::array<int32_T, 1U> &A_colind,
                             const ::coder::array<real_T, 1U> &A_val,
                             const ::coder::array<real_T, 2U> &x,
                             ::coder::array<real_T, 2U> &b);

static inline
void determine_rdnodes(boolean_T fullrank,
                              ::coder::array<boolean_T, 1U> &rdtags,
                              ::coder::array<int32_T, 1U> &rdnodes);

static inline
void extract_sub(::coder::SizeType n, const ::coder::array<int32_T, 1U> &crange,
                        const ::coder::array<int32_T, 1U> &eptr,
                        const ::coder::array<int32_T, 1U> &eind,
                        ::coder::array<int32_T, 1U> &iwork,
                        ::coder::array<boolean_T, 1U> &ntags,
                        ::coder::array<int32_T, 1U> &eptrloc,
                        ::coder::array<int32_T, 1U> &eindloc, int32_T *nnodes);

static inline
void f_WlsObject(::coder::SizeType degree, b_WlsObject *wls);

static inline
real_T find_kth_shortest_dist(::coder::array<real_T, 1U> &arr, ::coder::SizeType k,
                                     ::coder::SizeType l, ::coder::SizeType r);

static inline
void gen_vander(const ::coder::array<real_T, 2U> &us, ::coder::SizeType npoints,
                       ::coder::SizeType degree, ::coder::array<real_T, 2U> &V);

static inline
void gen_vander(const ::coder::array<real_T, 2U> &us, ::coder::SizeType npoints,
                       ::coder::SizeType degree,
                       const ::coder::array<real_T, 1U> &weights,
                       ::coder::array<real_T, 2U> &V);

static inline
void gen_vander_2d(const ::coder::array<real_T, 2U> &us, ::coder::SizeType npoints,
                          ::coder::SizeType degree, ::coder::array<real_T, 2U> &V);

static inline
void gen_vander_2d(const ::coder::array<real_T, 2U> &us, ::coder::SizeType npoints,
                          ::coder::SizeType degree,
                          const ::coder::array<real_T, 1U> &weights,
                          ::coder::array<real_T, 2U> &V);

static inline
void gen_vander_3d(const ::coder::array<real_T, 2U> &us, ::coder::SizeType npoints,
                          ::coder::SizeType degree, ::coder::array<real_T, 2U> &V);

static inline
void gen_vander_3d(const ::coder::array<real_T, 2U> &us, ::coder::SizeType npoints,
                          ::coder::SizeType degree,
                          const ::coder::array<real_T, 1U> &weights,
                          ::coder::array<real_T, 2U> &V);

static inline
void hexa_125(real_T xi, real_T eta, real_T zeta, real_T sfvals[125],
                     real_T sdvals[375]);

static inline
void hexa_216(real_T xi, real_T eta, real_T zeta, real_T sfvals[216],
                     real_T sdvals[648]);

static inline
void hexa_343(real_T xi, real_T eta, real_T zeta, real_T sfvals[343],
                     real_T sdvals[1029]);

static inline
void hexa_64(real_T xi, real_T eta, real_T zeta, real_T sfvals[64],
                    real_T sdvals[192]);

static inline
void hexa_gl_125(real_T xi, real_T eta, real_T zeta, real_T sfvals[125],
                        real_T sdvals[375]);

static inline
void hexa_gl_216(real_T xi, real_T eta, real_T zeta, real_T sfvals[216],
                        real_T sdvals[648]);

static inline
void hexa_gl_343(real_T xi, real_T eta, real_T zeta, real_T sfvals[343],
                        real_T sdvals[1029]);

static inline
void hexa_gl_64(real_T xi, real_T eta, real_T zeta, real_T sfvals[64],
                       real_T sdvals[192]);

static inline
void init_osusop(const ::coder::array<real_T, 2U> &mesh_coords,
                        const ::coder::array<ConnData, 1U> &mesh_elemtables,
                        const ::coder::array<uint64_T, 1U> &mesh_teids,
                        const ::coder::array<Stencils, 1U> &mesh_stencils,
                        ::coder::SizeType stclid, ::coder::array<int64_T, 1U> &A_row_ptr,
                        ::coder::array<int32_T, 1U> &A_col_ind,
                        ::coder::array<real_T, 1U> &A_val, int32_T *A_ncols,
                        ::coder::array<int32_T, 1U> &nnzs);

static inline
void insert_mem_crs(::coder::SizeType i, ::coder::array<int64_T, 1U> &row_ptr,
                           ::coder::array<int32_T, 1U> &col_ind,
                           ::coder::array<real_T, 1U> &val);

static inline
int64_T m2cFind(const ::coder::array<int32_T, 1U> &keys, ::coder::SizeType key,
                       int64_T b_first, int64_T last);

static inline
void m2cSort(::coder::array<int32_T, 1U> &keys, ::coder::SizeType b_first,
                    ::coder::SizeType last);

static inline
void majorityTransform(const b_WlsObject *r, c_WlsObject *r1);

static inline
void majorityTransform(const d_WlsObject *r, e_WlsObject *r1);

static inline
void mark_discontinuities_kernel(
    const ::coder::array<real_T, 2U> &mesh_coords,
    const ::coder::array<ConnData, 1U> &mesh_elemtables,
    const ::coder::array<uint64_T, 1U> &mesh_teids,
    const ::coder::array<int64_T, 1U> &mesh_node2elems_row_ptr,
    const ::coder::array<int32_T, 1U> &mesh_node2elems_col_ind,
    const ::coder::array<real_T, 1U> &mesh_elemh, real_T mesh_globalh,
    const ::coder::array<real_T, 2U> &fs, const ::coder::array<real_T, 2U> &df,
    const ::coder::array<real_T, 2U> &alpha,
    const ::coder::array<real_T, 2U> &beta, real_T cglobal, real_T clocal,
    real_T kappa1, real_T kappa0, ::coder::array<int8_T, 2U> &distags);

static inline
void obtain_nring_1d(const ::coder::array<ConnData, 1U> &mesh_elemtables,
                            const ::coder::array<uint64_T, 1U> &mesh_teids,
                            const ::coder::array<uint64_T, 2U> &mesh_sibhfs,
                            const ::coder::array<uint64_T, 1U> &mesh_v2hfid,
                            ::coder::SizeType vid, real_T ring, ::coder::SizeType maxnpnts,
                            ::coder::array<boolean_T, 1U> &vtags,
                            ::coder::array<boolean_T, 1U> &ftags,
                            ::coder::array<int32_T, 1U> &ngbvs, int32_T *nverts,
                            ::coder::array<int32_T, 1U> &ngbfs, int32_T *nfaces,
                            ::coder::array<uint64_T, 1U> &hebuf,
                            boolean_T *reflected, boolean_T *overflow);

static inline
void obtain_nring_2d(const ::coder::array<ConnData, 1U> &mesh_elemtables,
                            const ::coder::array<uint64_T, 1U> &mesh_teids,
                            const ::coder::array<uint64_T, 2U> &mesh_sibhfs,
                            const ::coder::array<uint64_T, 1U> &mesh_v2hfid,
                            ::coder::SizeType vid, real_T ring, ::coder::SizeType maxnpnts,
                            ::coder::array<boolean_T, 1U> &vtags,
                            ::coder::array<boolean_T, 1U> &ftags,
                            const ::coder::array<int32_T, 1U> &bridges,
                            ::coder::array<int32_T, 1U> &ngbvs, int32_T *nverts,
                            ::coder::array<int32_T, 1U> &ngbfs, int32_T *nfaces,
                            ::coder::array<uint64_T, 1U> &hebuf,
                            boolean_T *reflected, boolean_T *overflow);

static inline
void omp4mRecurPartMesh(::coder::SizeType n,
                               const ::coder::array<int32_T, 2U> &cells,
                               ::coder::SizeType dim, ::coder::SizeType nLevels, ::coder::SizeType nParts,
                               ::coder::array<Omp4mPart, 1U> &parts);

static inline
void prism_126(real_T xi, real_T eta, real_T zeta, real_T sfvals[126],
                      real_T sdvals[378]);

static inline
void prism_196(real_T xi, real_T eta, real_T zeta, real_T sfvals[196],
                      real_T sdvals[588]);

static inline
void prism_40(real_T xi, real_T eta, real_T zeta, real_T sfvals[40],
                     real_T sdvals[120]);

static inline
void prism_75(real_T xi, real_T eta, real_T zeta, real_T sfvals[75],
                     real_T sdvals[225]);

static inline
void prism_gl_126(real_T xi, real_T eta, real_T zeta, real_T sfvals[126],
                         real_T sdvals[378]);

static inline
void prism_gl_40(real_T xi, real_T eta, real_T zeta, real_T sfvals[40],
                        real_T sdvals[120]);

static inline
void prism_gl_75(real_T xi, real_T eta, real_T zeta, real_T sfvals[75],
                        real_T sdvals[225]);

static inline
void pyra_30(real_T xi, real_T eta, real_T zeta, real_T sfvals[30],
                    real_T sdvals[90]);

static inline
void pyra_55(real_T xi, real_T eta, real_T zeta, real_T sfvals[55],
                    real_T sdvals[165]);

static inline
void pyra_gl_30(real_T xi, real_T eta, real_T zeta, real_T sfvals[30],
                       real_T sdvals[90]);

static inline
void pyra_gl_55(real_T xi, real_T eta, real_T zeta, real_T sfvals[55],
                       real_T sdvals[165]);

static inline
void pyra_quadrules(::coder::SizeType degree, ::coder::array<real_T, 2U> &cs,
                           ::coder::array<real_T, 1U> &ws);

static inline
void quad_25(real_T xi, real_T eta, real_T sfvals[25],
                    real_T sdvals[50]);

static inline
void quad_36(real_T xi, real_T eta, real_T sfvals[36],
                    real_T sdvals[72]);

static inline
void quad_49(real_T xi, real_T eta, real_T sfvals[49],
                    real_T sdvals[98]);

static inline
void quad_gl_25(real_T xi, real_T eta, real_T sfvals[25],
                       real_T sdvals[50]);

static inline
void quad_gl_36(real_T xi, real_T eta, real_T sfvals[36],
                       real_T sdvals[72]);

static inline
void quad_gl_49(real_T xi, real_T eta, real_T sfvals[49],
                       real_T sdvals[98]);

static inline
void rdi_compute_oscind(
    real_T rdi_epsbeta, const ::coder::array<real_T, 2U> &mesh_coords,
    const ::coder::array<int64_T, 1U> &mesh_node2elems_row_ptr,
    const ::coder::array<int32_T, 1U> &mesh_node2elems_col_ind,
    const ::coder::array<real_T, 1U> &mesh_elemmeas, real_T mesh_globalh,
    const ::coder::array<real_T, 2U> &df,
    const ::coder::array<real_T, 2U> &alpha, ::coder::array<real_T, 2U> &beta);

static inline
void rdi_compute_osusind(
    const ::coder::array<int64_T, 1U> &rdi_A_row_ptr,
    const ::coder::array<int32_T, 1U> &rdi_A_col_ind,
    const ::coder::array<real_T, 1U> &rdi_A_val, boolean_T rdi_fullrank,
    const ::coder::array<real_T, 2U> &mesh_coords,
    const ::coder::array<int64_T, 1U> &mesh_node2elems_row_ptr,
    const ::coder::array<int32_T, 1U> &mesh_node2elems_col_ind,
    const ::coder::array<real_T, 2U> &fs, ::coder::array<real_T, 2U> &alphacell,
    ::coder::array<real_T, 2U> &alphanode);

static inline
void
rdi_contract_markers(::coder::array<int8_T, 2U> &distags,
                     const ::coder::array<real_T, 2U> &mesh_coords,
                     const ::coder::array<ConnData, 1U> &mesh_elemtables,
                     const ::coder::array<uint64_T, 1U> &mesh_teids,
                     const ::coder::array<int64_T, 1U> &mesh_node2nodes_row_ptr,
                     const ::coder::array<int32_T, 1U> &mesh_node2nodes_col_ind,
                     const ::coder::array<int64_T, 1U> &mesh_node2elems_row_ptr,
                     const ::coder::array<int32_T, 1U> &mesh_node2elems_col_ind,
                     ::coder::SizeType nlayers);

static inline
void
rdi_expand_markers(::coder::array<int8_T, 2U> &distags,
                   const ::coder::array<int64_T, 1U> &mesh_node2nodes_row_ptr,
                   const ::coder::array<int32_T, 1U> &mesh_node2nodes_col_ind,
                   ::coder::SizeType nlayers);

static inline
void rdi_mark_discontinuities(
    real_T rdi_cglobal, real_T rdi_clocal, real_T rdi_kappa0, real_T rdi_kappa1,
    const ::coder::array<real_T, 2U> &mesh_coords,
    const ::coder::array<ConnData, 1U> &mesh_elemtables,
    const ::coder::array<uint64_T, 1U> &mesh_teids,
    const ::coder::array<int64_T, 1U> &mesh_node2elems_row_ptr,
    const ::coder::array<int32_T, 1U> &mesh_node2elems_col_ind,
    const ::coder::array<real_T, 1U> &mesh_elemh, real_T mesh_globalh,
    const ::coder::array<real_T, 2U> &fs, const ::coder::array<real_T, 2U> &df,
    const ::coder::array<real_T, 2U> &alpha,
    const ::coder::array<real_T, 2U> &beta,
    ::coder::array<int8_T, 2U> &distags);

static inline
void rdi_update_osusop(RdiObject *rdi, WlsMesh *mesh, boolean_T interp0);

static inline
void rrqr_factor(const ::coder::array<real_T, 2U> &A, real_T thres,
                        ::coder::SizeType rowoffset, ::coder::SizeType coloffset, ::coder::SizeType m,
                        ::coder::SizeType n, ::coder::array<real_T, 2U> &QR,
                        ::coder::array<int32_T, 1U> &p, int32_T *rank,
                        ::coder::array<real_T, 1U> &work);

static inline
void rrqr_qmulti(const ::coder::array<real_T, 2U> &QR, ::coder::SizeType m,
                        ::coder::SizeType n, ::coder::SizeType rank, ::coder::array<real_T, 2U> &bs,
                        ::coder::SizeType nrhs, ::coder::array<real_T, 1U> &work);

static inline
void rrqr_rtsolve(const ::coder::array<real_T, 2U> &QR, ::coder::SizeType n,
                         ::coder::SizeType rank, ::coder::array<real_T, 2U> &bs,
                         ::coder::SizeType nrhs);

static inline
void sfe1_tabulate_shapefuncs(::coder::SizeType etype,
                                     const ::coder::array<real_T, 2U> &cs,
                                     ::coder::array<real_T, 2U> &sfvals,
                                     ::coder::array<real_T, 3U> &sdvals);

static inline
void sfe2_tabulate_equi_quad(::coder::SizeType etype,
                                    const ::coder::array<real_T, 2U> &cs,
                                    ::coder::array<real_T, 2U> &sfvals,
                                    ::coder::array<real_T, 3U> &sdvals);

static inline
void sfe2_tabulate_equi_tri(::coder::SizeType etype,
                                   const ::coder::array<real_T, 2U> &cs,
                                   ::coder::array<real_T, 2U> &sfvals,
                                   ::coder::array<real_T, 3U> &sdvals);

static inline
void sfe2_tabulate_fek_tri(::coder::SizeType etype,
                                  const ::coder::array<real_T, 2U> &cs,
                                  ::coder::array<real_T, 2U> &sfvals,
                                  ::coder::array<real_T, 3U> &sdvals);

static inline
void sfe2_tabulate_gl_quad(::coder::SizeType etype,
                                  const ::coder::array<real_T, 2U> &cs,
                                  ::coder::array<real_T, 2U> &sfvals,
                                  ::coder::array<real_T, 3U> &sdvals);

static inline
void sfe2_tabulate_gl_tri(::coder::SizeType etype,
                                 const ::coder::array<real_T, 2U> &cs,
                                 ::coder::array<real_T, 2U> &sfvals,
                                 ::coder::array<real_T, 3U> &sdvals);

static inline
void sfe2_tabulate_shapefuncs(::coder::SizeType etype,
                                     const ::coder::array<real_T, 2U> &cs,
                                     ::coder::array<real_T, 2U> &sfvals,
                                     ::coder::array<real_T, 3U> &sdvals);

static inline
void sfe3_tabulate_equi_hexa(::coder::SizeType etype,
                                    const ::coder::array<real_T, 2U> &cs,
                                    ::coder::array<real_T, 2U> &sfvals,
                                    ::coder::array<real_T, 3U> &sdvals);

static inline
void sfe3_tabulate_equi_prism(::coder::SizeType etype,
                                     const ::coder::array<real_T, 2U> &cs,
                                     ::coder::array<real_T, 2U> &sfvals,
                                     ::coder::array<real_T, 3U> &sdvals);

static inline
void sfe3_tabulate_equi_pyra(::coder::SizeType etype,
                                    const ::coder::array<real_T, 2U> &cs,
                                    ::coder::array<real_T, 2U> &sfvals,
                                    ::coder::array<real_T, 3U> &sdvals);

static inline
void sfe3_tabulate_equi_tet(::coder::SizeType etype,
                                   const ::coder::array<real_T, 2U> &cs,
                                   ::coder::array<real_T, 2U> &sfvals,
                                   ::coder::array<real_T, 3U> &sdvals);

static inline
void sfe3_tabulate_gl_hexa(::coder::SizeType etype,
                                  const ::coder::array<real_T, 2U> &cs,
                                  ::coder::array<real_T, 2U> &sfvals,
                                  ::coder::array<real_T, 3U> &sdvals);

static inline
void sfe3_tabulate_gl_prism(::coder::SizeType etype,
                                   const ::coder::array<real_T, 2U> &cs,
                                   ::coder::array<real_T, 2U> &sfvals,
                                   ::coder::array<real_T, 3U> &sdvals);

static inline
void sfe3_tabulate_gl_pyra(::coder::SizeType etype,
                                  const ::coder::array<real_T, 2U> &cs,
                                  ::coder::array<real_T, 2U> &sfvals,
                                  ::coder::array<real_T, 3U> &sdvals);

static inline
void sfe3_tabulate_gl_tet(::coder::SizeType etype,
                                 const ::coder::array<real_T, 2U> &cs,
                                 ::coder::array<real_T, 2U> &sfvals,
                                 ::coder::array<real_T, 3U> &sdvals);

static inline
void sfe3_tabulate_shapefuncs(::coder::SizeType etype,
                                     const ::coder::array<real_T, 2U> &cs,
                                     ::coder::array<real_T, 2U> &sfvals,
                                     ::coder::array<real_T, 3U> &sdvals);

static inline
void sfe_init(SfeObject *sfe, ::coder::SizeType etypes,
                     const ::coder::array<real_T, 2U> &xs);

static inline
void sfe_init(SfeObject *sfe, const ::coder::array<real_T, 2U> &xs);

static inline
void tabulate_quadratures(::coder::SizeType etype, ::coder::SizeType qd,
                                 ::coder::array<real_T, 2U> &cs,
                                 ::coder::array<real_T, 1U> &ws);

static inline
void tabulate_shapefuncs(::coder::SizeType etype,
                                const ::coder::array<real_T, 2U> &cs,
                                ::coder::array<real_T, 2U> &sfvals,
                                ::coder::array<real_T, 3U> &sdvals);

static inline
void tet_20(real_T xi, real_T eta, real_T zeta, real_T sfvals[20],
                   real_T sdvals[60]);

static inline
void tet_35(real_T xi, real_T eta, real_T zeta, real_T sfvals[35],
                   real_T sdvals[105]);

static inline
void tet_56(real_T xi, real_T eta, real_T zeta, real_T sfvals[56],
                   real_T sdvals[168]);

static inline
void tet_84(real_T xi, real_T eta, real_T zeta, real_T sfvals[84],
                   real_T sdvals[252]);

static inline
void tet_gl_20(real_T xi, real_T eta, real_T zeta, real_T sfvals[20],
                      real_T sdvals[60]);

static inline
void tet_gl_35(real_T xi, real_T eta, real_T zeta, real_T sfvals[35],
                      real_T sdvals[105]);

static inline
void tri_21(real_T xi, real_T eta, real_T sfvals[21], real_T sdvals[42]);

static inline
void tri_28(real_T xi, real_T eta, real_T sfvals[28], real_T sdvals[56]);

static inline
void tri_fek_15(real_T xi, real_T eta, real_T sfvals[15],
                       real_T sdvals[30]);

static inline
void tri_fek_21(real_T xi, real_T eta, real_T sfvals[21],
                       real_T sdvals[42]);

static inline
void tri_fek_28(real_T xi, real_T eta, real_T sfvals[28],
                       real_T sdvals[56]);

static inline
void tri_gl_21(real_T xi, real_T eta, real_T sfvals[21],
                      real_T sdvals[42]);

static inline
void tri_quadrules(::coder::SizeType degree, ::coder::array<real_T, 2U> &cs,
                          ::coder::array<real_T, 1U> &ws);

static inline
void
update_osusop(CrsMatrix *A, ::coder::array<int32_T, 1U> &nnzs,
              const ::coder::array<real_T, 2U> &mesh_coords,
              const ::coder::array<int64_T, 1U> &mesh_node2elems_row_ptr,
              const ::coder::array<int32_T, 1U> &mesh_node2elems_col_ind,
              const ::coder::array<Stencils, 1U> &mesh_stencils,
              ::coder::SizeType stclid);

static inline
void wls_buhmann_weights(const ::coder::array<real_T, 2U> &us,
                                ::coder::SizeType npoints, ::coder::SizeType degree,
                                const ::coder::array<real_T, 1U> &params_sh,
                                const ::coder::array<real_T, 2U> &params_pw,
                                ::coder::array<real_T, 1U> &ws);

static inline
void wls_eno_weights(const ::coder::array<real_T, 2U> &us,
                            ::coder::SizeType npoints, ::coder::SizeType degree,
                            const ::coder::array<real_T, 2U> &us_unscaled,
                            const ::coder::array<real_T, 1U> &params_sh,
                            const ::coder::array<real_T, 2U> &params_pw,
                            ::coder::array<real_T, 1U> &ws);

static inline
void wls_func(e_WlsObject *wls,
                     const ::coder::array<real_T, 2U> &eval_pnts,
                     ::coder::array<real_T, 2U> &varargout_1);

static inline
void wls_init(e_WlsObject *wls, const ::coder::array<real_T, 2U> &us,
                     const char_T weight_name_data[],
                     const ::coder::array<real_T, 1U> &weight_params_shared,
                     const ::coder::array<real_T, 2U> &weight_params_pointwise,
                     const ::coder::array<boolean_T, 1U> &weight_omit_rows,
                     ::coder::SizeType degree, ::coder::SizeType interp0, ::coder::SizeType nstpnts);

static inline
void wls_invdist_weights(const ::coder::array<real_T, 2U> &us,
                                ::coder::SizeType npoints, ::coder::SizeType degree,
                                const ::coder::array<real_T, 1U> &params_sh,
                                const ::coder::array<real_T, 2U> &params_pw,
                                ::coder::array<real_T, 1U> &ws);

static inline
void wls_invdist_weights(const ::coder::array<real_T, 2U> &us,
                                ::coder::SizeType npoints, real_T degree,
                                ::coder::array<real_T, 1U> &ws);

static inline
void wls_kernel(e_WlsObject *wls,
                       const ::coder::array<real_T, 2U> &eval_pnts,
                       ::coder::array<real_T, 2U> &vdops);

static inline
void wls_resize(e_WlsObject *wls, ::coder::SizeType dim, ::coder::SizeType nstpnts,
                       ::coder::SizeType degree);

static inline
void wls_solve_sys(e_WlsObject *wls, ::coder::array<real_T, 2U> &vdops);

static inline
void wls_update_rhs(e_WlsObject *wls);

static inline
void wlsmesh_compute_1ring(WlsMesh *mesh);

static inline
void wlsmesh_compute_meshprop(WlsMesh *mesh, ::coder::SizeType nrmidx);

static inline
void wlsmesh_compute_nodeparts(WlsMesh *mesh, ::coder::SizeType nparts);

static inline
void wlsmesh_compute_stencils(WlsMesh *mesh, ::coder::SizeType stclidx);

static inline
void wlsmesh_compute_stencils(WlsMesh *mesh, ::coder::SizeType stclidx,
                                     real_T krings);

} // namespace rdi_kernel

// Function Definitions
namespace rdi_kernel {
static ::coder::SizeType append_wlsmesh_kring(WlsMesh *mesh)
{
  static const char_T name[7]{'O', 'n', 'e', 'R', 'i', 'n', 'g'};
  ::coder::SizeType stclidx;
  stclidx = mesh->stencils.size(0) + 1;
  // Stencils - Object containing the stencils of vertices in a WlsMesh based on
  mesh->stencils.set_size(mesh->stencils.size(0) + 1);
  mesh->stencils[stclidx - 1].reflected.set_size(0);
  mesh->stencils[stclidx - 1].vidmap.set_size(0);
  mesh->stencils[stclidx - 1].ngbverts.row_ptr.set_size(0);
  mesh->stencils[stclidx - 1].ngbverts.ncols = 0;
  //  Construct A.row_ptr
  mesh->stencils[stclidx - 1].ngbverts.col_ind.set_size(0);
  mesh->stencils[stclidx - 1].ngbelems.row_ptr.set_size(0);
  mesh->stencils[stclidx - 1].ngbelems.ncols = 0;
  //  Construct A.row_ptr
  mesh->stencils[stclidx - 1].ngbelems.col_ind.set_size(0);
  mesh->stencils[stclidx - 1].name.set_size(1, 7);
  for (::coder::SizeType i{0}; i < 7; i++) {
    mesh->stencils[stclidx - 1].name[i] = name[i];
  }
  return stclidx;
}

// assemble_body - Interior body assembly for OSUS operator
static inline
void
assemble_body(RdiObject *rdi, const ::coder::array<real_T, 2U> &mesh_coords,
              const ::coder::array<ConnData, 1U> &mesh_elemtables,
              const ::coder::array<uint64_T, 1U> &mesh_teids,
              const ::coder::array<int64_T, 1U> &mesh_node2elems_row_ptr,
              const ::coder::array<int32_T, 1U> &mesh_node2elems_col_ind,
              const ::coder::array<Stencils, 1U> &mesh_stencils, ::coder::SizeType stclid,
              boolean_T interp0)
{
  if ((stclid != rdi->stclid) && (stclid != rdi->extstclid)) {
    m2cErrMsgIdAndTxt("assemble_body:badStencilID", "bad stencil ID (index) %d",
                      (int)stclid);
  }
  //  Set fullrank to be true upon input
  if (mesh_stencils[stclid - 1].vidmap.size(0) != 0) {
    //  For extended stencils (resolving rank deficiencies), must be serial
    rdi->fullrank = b_assemble_body_range(
        rdi->A.row_ptr, rdi->A.col_ind, mesh_coords, mesh_elemtables,
        mesh_teids, mesh_stencils[stclid - 1].ngbverts.row_ptr,
        mesh_stencils[stclid - 1].ngbverts.col_ind, mesh_node2elems_row_ptr,
        mesh_node2elems_col_ind, rdi->degree, interp0,
        mesh_stencils[stclid - 1].vidmap,
        mesh_stencils[stclid - 1].vidmap.size(0), rdi->A.val, rdi->rdtags);
  } else {
    //  Must be serial
    rdi->fullrank = c_assemble_body_range(
        rdi->A.row_ptr, rdi->A.col_ind, mesh_coords, mesh_elemtables,
        mesh_teids, mesh_stencils[stclid - 1].ngbverts.row_ptr,
        mesh_stencils[stclid - 1].ngbverts.col_ind, mesh_node2elems_row_ptr,
        mesh_node2elems_col_ind, rdi->degree, interp0, mesh_coords.size(0),
        rdi->A.val, rdi->rdtags);
  }
}

static void
assemble_body_range(const ::coder::array<int64_T, 1U> &row_ptr,
                    const ::coder::array<int32_T, 1U> &col_ind,
                    const ::coder::array<real_T, 2U> &mesh_coords,
                    const ::coder::array<ConnData, 1U> &mesh_elemtables,
                    const ::coder::array<uint64_T, 1U> &mesh_teids,
                    const ::coder::array<int64_T, 1U> &stcl_row_ptr,
                    const ::coder::array<int32_T, 1U> &stcl_col_ind,
                    const ::coder::array<int64_T, 1U> &n2e_row_ptr,
                    const ::coder::array<int32_T, 1U> &n2e_col_ind,
                    ::coder::SizeType degree, boolean_T interp0,
                    const ::coder::array<int32_T, 1U> &nrange, ::coder::SizeType istart,
                    ::coder::SizeType iend, ::coder::array<real_T, 1U> &val,
                    ::coder::array<boolean_T, 1U> &rdtags, boolean_T *fullrank)
{
  static const char_T name[7]{'B', 'u', 'h', 'm', 'a', 'n', 'n'};
  ::coder::array<real_T, 2U> coeffs_;
  ::coder::array<real_T, 2U> w_params_pointwise;
  ::coder::array<real_T, 2U> xs_;
  ::coder::array<real_T, 1U> Ns_;
  ::coder::array<real_T, 1U> w_params_shared;
  ::coder::array<int32_T, 1U> eids_;
  ::coder::array<int32_T, 1U> nodes_;
  ::coder::array<boolean_T, 1U> w_omit_rows;
  b_WlsObject expl_temp;
  e_WlsObject wls_;
  ::coder::SizeType c_loop_ub;
  ::coder::SizeType dim;
  ::coder::SizeType loop_ub;
  char_T wgts__name_data[7];
  f_WlsObject(degree, &expl_temp);
  wls_.runtimes.size[0] = expl_temp.runtimes.size[0];
  loop_ub = expl_temp.runtimes.size[0];
  if (loop_ub - 1 >= 0) {
    std::copy(&expl_temp.runtimes.data[0], &expl_temp.runtimes.data[loop_ub],
              &wls_.runtimes.data[0]);
  }
  wls_.QRt.set_size(expl_temp.QRt.size(0), expl_temp.QRt.size(1));
  loop_ub = expl_temp.QRt.size(1) * expl_temp.QRt.size(0);
  for (::coder::SizeType i{0}; i < loop_ub; i++) {
    wls_.QRt[i] = expl_temp.QRt[i];
  }
  wls_.rowmajor = expl_temp.rowmajor;
  wls_.work.set_size(expl_temp.work.size(0));
  loop_ub = expl_temp.work.size(0);
  for (::coder::SizeType i{0}; i < loop_ub; i++) {
    wls_.work[i] = expl_temp.work[i];
  }
  wls_.jpvt.set_size(expl_temp.jpvt.size(0));
  loop_ub = expl_temp.jpvt.size(0);
  for (::coder::SizeType i{0}; i < loop_ub; i++) {
    wls_.jpvt[i] = expl_temp.jpvt[i];
  }
  wls_.fullrank = expl_temp.fullrank;
  wls_.rank = expl_temp.rank;
  wls_.ncols = expl_temp.ncols;
  wls_.nrows = expl_temp.nrows;
  wls_.nevpnts = expl_temp.nevpnts;
  wls_.rhs.set_size(expl_temp.rhs.size(0), expl_temp.rhs.size(1));
  loop_ub = expl_temp.rhs.size(1) * expl_temp.rhs.size(0);
  for (::coder::SizeType i{0}; i < loop_ub; i++) {
    wls_.rhs[i] = expl_temp.rhs[i];
  }
  wls_.QR.set_size(expl_temp.QR.size(0), expl_temp.QR.size(1));
  loop_ub = expl_temp.QR.size(1) * expl_temp.QR.size(0);
  for (::coder::SizeType i{0}; i < loop_ub; i++) {
    wls_.QR[i] = expl_temp.QR[i];
  }
  wls_.V.set_size(expl_temp.V.size(0), expl_temp.V.size(1));
  loop_ub = expl_temp.V.size(1) * expl_temp.V.size(0);
  for (::coder::SizeType i{0}; i < loop_ub; i++) {
    wls_.V[i] = expl_temp.V[i];
  }
  wls_.hs_inv.set_size(1, expl_temp.hs_inv.size[1]);
  loop_ub = expl_temp.hs_inv.size[1];
  for (::coder::SizeType i{0}; i < loop_ub; i++) {
    wls_.hs_inv[i] = expl_temp.hs_inv.data[i];
  }
  wls_.rweights.set_size(expl_temp.rweights.size(0));
  loop_ub = expl_temp.rweights.size(0);
  for (::coder::SizeType i{0}; i < loop_ub; i++) {
    wls_.rweights[i] = expl_temp.rweights[i];
  }
  wls_.origin.size[1] = expl_temp.origin.size[1];
  wls_.origin.size[0] = 1;
  loop_ub = expl_temp.origin.size[1];
  if (loop_ub - 1 >= 0) {
    std::copy(&expl_temp.origin.data[0], &expl_temp.origin.data[loop_ub],
              &wls_.origin.data[0]);
  }
  wls_.us.set_size(expl_temp.us.size(0), expl_temp.us.size(1));
  loop_ub = expl_temp.us.size(0);
  for (::coder::SizeType i{0}; i < loop_ub; i++) {
    ::coder::SizeType b_loop_ub;
    b_loop_ub = expl_temp.us.size(1);
    for (::coder::SizeType i1{0}; i1 < b_loop_ub; i1++) {
      wls_.us[i1 + wls_.us.size(1) * i] =
          expl_temp.us[i1 + expl_temp.us.size(1) * i];
    }
  }
  wls_.stride = expl_temp.stride;
  wls_.interp0 = expl_temp.interp0;
  wls_.unimono = expl_temp.unimono;
  wls_.order = expl_temp.order;
  wls_.degree = expl_temp.degree;
  wls_.nstpnts = expl_temp.nstpnts;
  w_params_shared.set_size(0);
  w_params_pointwise.set_size(::coder::SizeType(0), ::coder::SizeType(0));
  w_omit_rows.set_size(0);
  for (::coder::SizeType i{0}; i < 7; i++) {
    wgts__name_data[i] = name[i];
  }
  // Local buffers with patterns "(\w+)?params_pointwise" are legitimate
  if (istart <= iend) {
    dim = mesh_coords.size(1);
    c_loop_ub = mesh_coords.size(1);
  }
  for (::coder::SizeType b_i{istart}; b_i <= iend; b_i++) {
    int64_T c_i;
    int64_T n;
    ::coder::SizeType k;
    ::coder::SizeType nid;
    if (nrange.size(0) == 0) {
      nid = b_i - 1;
    } else {
      nid = nrange[b_i - 1] - 1;
    }
    //  Fetch local data
    n = stcl_row_ptr[nid + 1] - stcl_row_ptr[nid];
    nodes_.set_size((n + 1L));
    xs_.set_size((n + 1L), mesh_coords.size(1));
    nodes_[0] = nid + 1;
    for (::coder::SizeType i{0}; i < c_loop_ub; i++) {
      xs_[i] = 0.0;
    }
    for (c_i = 2L; c_i - 1L <= n; c_i++) {
      k = stcl_col_ind[static_cast<::coder::SizeType>((stcl_row_ptr[nid] + c_i) - 2L) -
                       1];
      nodes_[c_i - 1] = k;
      for (::coder::SizeType j{0}; j < dim; j++) {
        xs_[j + xs_.size(1) * (c_i - 1)] =
            mesh_coords[j + mesh_coords.size(1) * (k - 1)] -
            mesh_coords[j + mesh_coords.size(1) * nid];
      }
    }
    //  Compute wls
    wls_init(&wls_, xs_, wgts__name_data, w_params_shared, w_params_pointwise,
             w_omit_rows, degree, interp0, xs_.size(0));
    if (wls_.rank < 0) {
      //  LAPACK error
      m2cErrMsgIdAndTxt("assemble_body_range:badLapack",
                        "LAPACK error code %d for node %d", wls_.rank, (int)nid + 1);
    }
    if (!wls_.fullrank) {
      //  Not full rank
      rdtags[nid] = true;
      *fullrank = false;
    } else {
      int64_T m;
      ::coder::SizeType b_dim;
      ::coder::SizeType b_m;
      ::coder::SizeType b_n;
      //  Get centers and shape function Ns at centers
      b_dim = mesh_coords.size(1) - 1;
      //  # of nearby elements
      m = n2e_row_ptr[nid + 1] - n2e_row_ptr[nid];
      eids_.set_size(m);
      xs_.set_size(m, mesh_coords.size(1));
      Ns_.set_size(m);
      for (int64_T b_j{1L}; b_j <= m; b_j++) {
        uint64_T c;
        ::coder::SizeType eid;
        ::coder::SizeType leid;
        ::coder::SizeType npe;
        eid = n2e_col_ind[static_cast<::coder::SizeType>((n2e_row_ptr[nid] + b_j) - 1L) -
                          1] -
              1;
        eids_[b_j - 1] = eid + 1;
        c = mesh_teids[eid] & 255UL;
        leid = (mesh_teids[eid] >> 8) - 1;
        npe =
            mesh_elemtables
                [static_cast<::coder::SizeType>(
                     mesh_teids[n2e_col_ind[static_cast<::coder::SizeType>(
                                                (n2e_row_ptr[nid] + b_j) - 1L) -
                                            1] -
                                1] &
                     255UL) -
                 1]
                    .conn.size(1);
        Ns_[b_j - 1] =
            1.0 /
            static_cast<real_T>(
                mesh_elemtables[c - 1].conn.size(1));
        //  Compute localized center
        for (k = 0; k <= b_dim; k++) {
          xs_[k + xs_.size(1) * (b_j - 1)] = mesh_coords
              [k + mesh_coords.size(1) *
                       (mesh_elemtables[c - 1]
                            .conn[mesh_elemtables[c - 1]
                                      .conn.size(1) *
                                  leid] -
                        1)];
        }
        for (::coder::SizeType d_i{2}; d_i <= npe; d_i++) {
          for (k = 0; k <= b_dim; k++) {
            xs_[k + xs_.size(1) * (b_j - 1)] =
                xs_[k + xs_.size(1) * (b_j - 1)] +
                mesh_coords
                    [k + mesh_coords.size(1) *
                             (mesh_elemtables[c - 1].conn
                                  [(d_i +
                                    mesh_elemtables[c - 1]
                                            .conn.size(1) *
                                        leid) -
                                   1] -
                              1)];
          }
        }
        for (k = 0; k <= b_dim; k++) {
          xs_[k + xs_.size(1) * (b_j - 1)] =
              xs_[k + xs_.size(1) * (b_j - 1)] *
                  Ns_[b_j - 1] -
              mesh_coords[k + mesh_coords.size(1) * nid];
        }
      }
      //  Compute WLS fitting at centers
      wls_func(&wls_, xs_, coeffs_);
      //  Compute local OSUS coefficients
      b_m = Ns_.size(0) - 1;
      b_n = coeffs_.size(0);
      //  WALF fitting
      for (k = 0; k < b_n; k++) {
        for (::coder::SizeType j{0}; j <= b_m; j++) {
          coeffs_[j + coeffs_.size(1) * k] =
              coeffs_[j + coeffs_.size(1) * k] * Ns_[j];
        }
      }
      //  Substract linear interp
      for (::coder::SizeType j{0}; j <= b_m; j++) {
        coeffs_[j] = coeffs_[j] - Ns_[j];
      }
      //  Add to global CRS
      b_n = eids_.size(0);
      b_m = nodes_.size(0);
      for (::coder::SizeType j{0}; j < b_n; j++) {
        ::coder::SizeType r;
        r = eids_[j];
        for (::coder::SizeType d_i{0}; d_i < b_m; d_i++) {
          int64_T b_k;
          b_k = 0L;
          if (r < row_ptr.size(0)) {
            int64_T bounds_idx_1;
            //  Perform linear search
            bounds_idx_1 = row_ptr[r] - 1L;
            //  Perform linear search
            c_i = row_ptr[r - 1];
            while (c_i <= bounds_idx_1) {
              if (col_ind[c_i - 1] == nodes_[d_i]) {
                b_k = c_i;
                c_i = bounds_idx_1 + 1L;
              } else {
                c_i++;
              }
            }
          }
          val[b_k - 1] =
              val[b_k - 1] +
              coeffs_[j + coeffs_.size(1) * d_i];
        }
      }
    }
  }
}

static void assemble_body_task(
    RdiObject *rdi, const ::coder::array<real_T, 2U> &mesh_coords,
    const ::coder::array<ConnData, 1U> &mesh_elemtables,
    const ::coder::array<uint64_T, 1U> &mesh_teids,
    const ::coder::array<int64_T, 1U> &mesh_node2elems_row_ptr,
    const ::coder::array<int32_T, 1U> &mesh_node2elems_col_ind,
    const ::coder::array<Omp4mPart, 1U> &mesh_nodeparts,
    const ::coder::array<int64_T, 1U> &stcl_row_ptr,
    const ::coder::array<int32_T, 1U> &stcl_col_ind, boolean_T interp0)
{
  ::coder::SizeType lvl;
  boolean_T exitg1;
  //  kernel for task
  lvl = 0;
  exitg1 = false;
  while ((!exitg1) && (lvl <= mesh_nodeparts.size(0) - 1)) {
    ::coder::SizeType i;
#pragma omp single
    { // single
      i = mesh_nodeparts[lvl].nparts;
      for (::coder::SizeType part{0}; part < i; part++) {
#pragma omp task default(shared)
        { // task
          assemble_body_range(
              rdi->A.row_ptr, rdi->A.col_ind, mesh_coords, mesh_elemtables,
              mesh_teids, stcl_row_ptr, stcl_col_ind, mesh_node2elems_row_ptr,
              mesh_node2elems_col_ind, rdi->degree, interp0,
              mesh_nodeparts[lvl].part_list, mesh_nodeparts[lvl].part_ptr[part],
              mesh_nodeparts[lvl].part_ptr[part + 1] - 1, rdi->A.val,
              rdi->rdtags, &rdi->fullrank);
        } // task
      }
    } // single
    //  implicit barrier
    if (mesh_nodeparts[lvl].shared_ents.size(0) != 0) {
#pragma omp single nowait
      { // single
        assemble_body_range(rdi->A.row_ptr, rdi->A.col_ind, mesh_coords,
                            mesh_elemtables, mesh_teids, stcl_row_ptr,
                            stcl_col_ind, mesh_node2elems_row_ptr,
                            mesh_node2elems_col_ind, rdi->degree, interp0,
                            mesh_nodeparts[lvl].shared_ents, 1,
                            mesh_nodeparts[lvl].shared_ents.size(0), rdi->A.val,
                            rdi->rdtags, &rdi->fullrank);
      } // single
      exitg1 = true;
    } else {
      lvl++;
    }
  }
}

// assemble_surf - Surface assembler for OSUS operator
static void
assemble_surf(RdiObject *rdi, const ::coder::array<real_T, 2U> &mesh_coords,
              const ::coder::array<ConnData, 1U> &mesh_elemtables,
              const ::coder::array<uint64_T, 1U> &mesh_teids,
              const ::coder::array<NormalsData, 1U> &mesh_nrmstables,
              const ::coder::array<int64_T, 1U> &mesh_node2elems_row_ptr,
              const ::coder::array<int32_T, 1U> &mesh_node2elems_col_ind,
              const ::coder::array<Stencils, 1U> &mesh_stencils, ::coder::SizeType stclid,
              boolean_T interp0)
{
  ::coder::SizeType nrmid;
  if ((stclid != rdi->stclid) && (stclid != rdi->extstclid)) {
    m2cErrMsgIdAndTxt("assemble_surf:badStencilID", "bad stencil ID (index) %d",
                      (int)stclid);
  }
  //  Set fullrank to be true upon input
  nrmid = rdi->nrmid;
  if (mesh_stencils[stclid - 1].vidmap.size(0) != 0) {
    //  For extended stencils (resolving rank deficiencies), must be serial
    rdi->fullrank = b_assemble_surf_range(
        rdi->A.row_ptr, rdi->A.col_ind, mesh_coords, mesh_elemtables,
        mesh_teids, mesh_stencils[stclid - 1].ngbverts.row_ptr,
        mesh_stencils[stclid - 1].ngbverts.col_ind, mesh_node2elems_row_ptr,
        mesh_node2elems_col_ind, mesh_nrmstables[nrmid - 1].normals,
        rdi->degree, interp0, mesh_stencils[stclid - 1].vidmap,
        mesh_stencils[stclid - 1].vidmap.size(0), rdi->A.val, rdi->rdtags);
  } else {
    //  Must be serial
    rdi->fullrank = c_assemble_surf_range(
        rdi->A.row_ptr, rdi->A.col_ind, mesh_coords, mesh_elemtables,
        mesh_teids, mesh_stencils[stclid - 1].ngbverts.row_ptr,
        mesh_stencils[stclid - 1].ngbverts.col_ind, mesh_node2elems_row_ptr,
        mesh_node2elems_col_ind, mesh_nrmstables[nrmid - 1].normals,
        rdi->degree, interp0, mesh_coords.size(0), rdi->A.val, rdi->rdtags);
  }
}

static void
assemble_surf_range(const ::coder::array<int64_T, 1U> &row_ptr,
                    const ::coder::array<int32_T, 1U> &col_ind,
                    const ::coder::array<real_T, 2U> &mesh_coords,
                    const ::coder::array<ConnData, 1U> &mesh_elemtables,
                    const ::coder::array<uint64_T, 1U> &mesh_teids,
                    const ::coder::array<int64_T, 1U> &stcl_row_ptr,
                    const ::coder::array<int32_T, 1U> &stcl_col_ind,
                    const ::coder::array<int64_T, 1U> &n2e_row_ptr,
                    const ::coder::array<int32_T, 1U> &n2e_col_ind,
                    const ::coder::array<real_T, 2U> &nrms, ::coder::SizeType degree,
                    boolean_T interp0,
                    const ::coder::array<int32_T, 1U> &nrange, ::coder::SizeType istart,
                    ::coder::SizeType iend, ::coder::array<real_T, 1U> &val,
                    ::coder::array<boolean_T, 1U> &rdtags, boolean_T *fullrank)
{
  static const char_T name[7]{'B', 'u', 'h', 'm', 'a', 'n', 'n'};
  ::coder::array<real_T, 2U> coeffs_;
  ::coder::array<real_T, 2U> us_;
  ::coder::array<real_T, 2U> wgts__params_pointwise;
  ::coder::array<real_T, 2U> xs_;
  ::coder::array<real_T, 1U> Ns_;
  ::coder::array<real_T, 1U> w_params_shared;
  ::coder::array<int32_T, 1U> eids_;
  ::coder::array<int32_T, 1U> nodes_;
  ::coder::array<boolean_T, 1U> w_omit_rows;
  b_WlsObject r;
  c_WlsObject expl_temp;
  d_WlsObject b_expl_temp;
  e_WlsObject wls_;
  real_T t_data[6];
  ::coder::SizeType b_loop_ub;
  ::coder::SizeType dim;
  ::coder::SizeType i;
  ::coder::SizeType i1;
  ::coder::SizeType loop_ub;
  ::coder::SizeType t_size_idx_1;
  char_T wgts__name_data[7];
  f_WlsObject(degree, &r);
  majorityTransform(&r, &expl_temp);
  b_expl_temp.us.set_size(expl_temp.us.size(0), expl_temp.us.size(1));
  loop_ub = expl_temp.us.size(1) * expl_temp.us.size(0);
  for (i = 0; i < loop_ub; i++) {
    b_expl_temp.us[i] = expl_temp.us[i];
  }
  b_expl_temp.hs_inv.set_size(expl_temp.hs_inv.size[0], 1);
  loop_ub = expl_temp.hs_inv.size[0];
  for (i = 0; i < loop_ub; i++) {
    b_expl_temp.hs_inv[i] = expl_temp.hs_inv.data[i];
  }
  b_expl_temp.runtimes.size[0] = expl_temp.runtimes.size[0];
  loop_ub = expl_temp.runtimes.size[0];
  if (loop_ub - 1 >= 0) {
    std::copy(&expl_temp.runtimes.data[0], &expl_temp.runtimes.data[loop_ub],
              &b_expl_temp.runtimes.data[0]);
  }
  b_expl_temp.QRt.set_size(expl_temp.QRt.size(0), expl_temp.QRt.size(1));
  loop_ub = expl_temp.QRt.size(1) * expl_temp.QRt.size(0);
  for (i = 0; i < loop_ub; i++) {
    b_expl_temp.QRt[i] = expl_temp.QRt[i];
  }
  b_expl_temp.rowmajor = expl_temp.rowmajor;
  b_expl_temp.work.set_size(expl_temp.work.size(0));
  loop_ub = expl_temp.work.size(0);
  for (i = 0; i < loop_ub; i++) {
    b_expl_temp.work[i] = expl_temp.work[i];
  }
  b_expl_temp.jpvt.set_size(expl_temp.jpvt.size(0));
  loop_ub = expl_temp.jpvt.size(0);
  for (i = 0; i < loop_ub; i++) {
    b_expl_temp.jpvt[i] = expl_temp.jpvt[i];
  }
  b_expl_temp.fullrank = expl_temp.fullrank;
  b_expl_temp.rank = expl_temp.rank;
  b_expl_temp.ncols = expl_temp.ncols;
  b_expl_temp.nrows = expl_temp.nrows;
  b_expl_temp.nevpnts = expl_temp.nevpnts;
  b_expl_temp.rhs.set_size(expl_temp.rhs.size(0), expl_temp.rhs.size(1));
  loop_ub = expl_temp.rhs.size(1) * expl_temp.rhs.size(0);
  for (i = 0; i < loop_ub; i++) {
    b_expl_temp.rhs[i] = expl_temp.rhs[i];
  }
  b_expl_temp.QR.set_size(expl_temp.QR.size(0), expl_temp.QR.size(1));
  loop_ub = expl_temp.QR.size(1) * expl_temp.QR.size(0);
  for (i = 0; i < loop_ub; i++) {
    b_expl_temp.QR[i] = expl_temp.QR[i];
  }
  b_expl_temp.V.set_size(expl_temp.V.size(0), expl_temp.V.size(1));
  loop_ub = expl_temp.V.size(1) * expl_temp.V.size(0);
  for (i = 0; i < loop_ub; i++) {
    b_expl_temp.V[i] = expl_temp.V[i];
  }
  b_expl_temp.rweights.set_size(expl_temp.rweights.size(0));
  loop_ub = expl_temp.rweights.size(0);
  for (i = 0; i < loop_ub; i++) {
    b_expl_temp.rweights[i] = expl_temp.rweights[i];
  }
  b_expl_temp.origin.size[1] = 1;
  b_expl_temp.origin.size[0] = expl_temp.origin.size[0];
  loop_ub = expl_temp.origin.size[0];
  if (loop_ub - 1 >= 0) {
    std::copy(&expl_temp.origin.data[0], &expl_temp.origin.data[loop_ub],
              &b_expl_temp.origin.data[0]);
  }
  b_expl_temp.stride = expl_temp.stride;
  b_expl_temp.interp0 = expl_temp.interp0;
  b_expl_temp.unimono = expl_temp.unimono;
  b_expl_temp.order = expl_temp.order;
  b_expl_temp.degree = expl_temp.degree;
  b_expl_temp.nstpnts = expl_temp.nstpnts;
  majorityTransform(&b_expl_temp, &wls_);
  w_params_shared.set_size(0);
  w_omit_rows.set_size(0);
  for (i = 0; i < 7; i++) {
    wgts__name_data[i] = name[i];
  }
  // Local buffers with patterns "wgts__params_pointwise" are legitimate
  if (istart <= iend) {
    dim = mesh_coords.size(1);
    b_loop_ub = mesh_coords.size(1);
    t_size_idx_1 = nrms.size(1) - 1;
    i1 = nrms.size(1);
  }
  for (::coder::SizeType b_i{istart}; b_i <= iend; b_i++) {
    int64_T n;
    ::coder::SizeType b_dim;
    ::coder::SizeType b_npoints;
    ::coder::SizeType k;
    ::coder::SizeType nid;
    ::coder::SizeType npoints;
    boolean_T b;
    boolean_T b1;
    if (nrange.size(0) == 0) {
      nid = b_i - 1;
    } else {
      nid = nrange[b_i - 1] - 1;
    }
    //  Fetch local data
    n = stcl_row_ptr[nid + 1] - stcl_row_ptr[nid];
    nodes_.set_size((n + 1L));
    xs_.set_size((n + 1L), mesh_coords.size(1));
    nodes_[0] = nid + 1;
    for (i = 0; i < b_loop_ub; i++) {
      xs_[i] = 0.0;
    }
    for (int64_T c_i{2L}; c_i - 1L <= n; c_i++) {
      k = stcl_col_ind[static_cast<::coder::SizeType>((stcl_row_ptr[nid] + c_i) - 2L) -
                       1];
      nodes_[c_i - 1] = k;
      for (::coder::SizeType j{0}; j < dim; j++) {
        xs_[j + xs_.size(1) * (c_i - 1)] =
            mesh_coords[j + mesh_coords.size(1) * (k - 1)] -
            mesh_coords[j + mesh_coords.size(1) * nid];
      }
    }
    npoints = xs_.size(0);
    //  Get normal direction
    if (i1 == 2) {
      t_data[0] = -nrms[nrms.size(1) * nid + 1];
      t_data[t_size_idx_1] = nrms[nrms.size(1) * nid];
    } else {
      real_T a;
      real_T a_tmp;
      a_tmp = std::abs(nrms[nrms.size(1) * nid]);
      if ((a_tmp > std::abs(nrms[nrms.size(1) * nid + 1])) &&
          (a_tmp > std::abs(nrms[nrms.size(1) * nid + 2]))) {
        t_data[0] = -nrms[nrms.size(1) * nid] * nrms[nrms.size(1) * nid + 1];
        t_data[t_size_idx_1] =
            1.0 - nrms[nrms.size(1) * nid + 1] * nrms[nrms.size(1) * nid + 1];
        t_data[t_size_idx_1 * 2] =
            -nrms[nrms.size(1) * nid + 1] * nrms[nrms.size(1) * nid + 2];
      } else {
        t_data[0] = 1.0 - nrms[nrms.size(1) * nid] * nrms[nrms.size(1) * nid];
        t_data[t_size_idx_1] =
            -nrms[nrms.size(1) * nid] * nrms[nrms.size(1) * nid + 1];
        t_data[t_size_idx_1 * 2] =
            -nrms[nrms.size(1) * nid] * nrms[nrms.size(1) * nid + 2];
      }
      a_tmp = t_data[t_size_idx_1 * 2];
      a = std::sqrt((t_data[0] * t_data[0] +
                     t_data[t_size_idx_1] * t_data[t_size_idx_1]) +
                    a_tmp * a_tmp);
      t_data[0] /= a;
      t_data[t_size_idx_1] /= a;
      t_data[t_size_idx_1 * 2] /= a;
      //  cross
      t_data[1] = t_data[t_size_idx_1 * 2] * nrms[nrms.size(1) * nid + 1] -
                  t_data[t_size_idx_1] * nrms[nrms.size(1) * nid + 2];
      t_data[t_size_idx_1 + 1] =
          t_data[0] * nrms[nrms.size(1) * nid + 2] -
          nrms[nrms.size(1) * nid] * t_data[t_size_idx_1 * 2];
      t_data[t_size_idx_1 * 2 + 1] =
          nrms[nrms.size(1) * nid] * t_data[t_size_idx_1] -
          t_data[0] * nrms[nrms.size(1) * nid + 1];
    }
    //  Project onto tangent plane
    b_npoints = xs_.size(0) - 1;
    b_dim = xs_.size(1) - 2;
    us_.set_size(xs_.size(0), xs_.size(1) - 1);
    for (k = 0; k <= b_npoints; k++) {
      for (::coder::SizeType j{0}; j <= b_dim; j++) {
        us_[j + us_.size(1) * k] = 0.0;
      }
    }
    //  matrix-matrix, using mem-efficient loop
    for (::coder::SizeType ii{0}; ii <= b_npoints; ii++) {
      for (k = 0; k <= b_dim + 1; k++) {
        for (::coder::SizeType j{0}; j <= b_dim; j++) {
          us_[j + us_.size(1) * ii] =
              us_[j + us_.size(1) * ii] +
              xs_[k + xs_.size(1) * ii] * t_data[j + t_size_idx_1 * k];
        }
      }
    }
    //  Compute normal matches as additional weights to WLS
    wgts__params_pointwise.set_size(xs_.size(0), 1);
    wgts__params_pointwise[0] = 1.0;
    b = true;
    i = wgts__params_pointwise.size(0);
    b1 = wgts__params_pointwise.size(0) <= 0;
    loop_ub = 0;
    for (::coder::SizeType j{2}; j <= npoints; j++) {
      real_T sigma;
      ::coder::SizeType sigma_tmp;
      if (b1 || (j - 1 >= i)) {
        loop_ub = 0;
        b = true;
      } else if (b) {
        b = false;
        loop_ub = (j - 1) % wgts__params_pointwise.size(0) +
                  (j - 1) / wgts__params_pointwise.size(0);
      } else if (loop_ub > 2147483646) {
        loop_ub = (j - 1) % wgts__params_pointwise.size(0) +
                  (j - 1) / wgts__params_pointwise.size(0);
      } else {
        loop_ub++;
        if (loop_ub > wgts__params_pointwise.size(0) - 1) {
          loop_ub = (loop_ub - wgts__params_pointwise.size(0)) + 1;
        }
      }
      sigma_tmp = nodes_[j - 1] - 1;
      sigma = nrms[nrms.size(1) * nid] * nrms[nrms.size(1) * sigma_tmp] +
              nrms[nrms.size(1) * nid + 1] * nrms[nrms.size(1) * sigma_tmp + 1];
      if (nrms.size(1) > 2) {
        sigma +=
            nrms[nrms.size(1) * nid + 2] * nrms[nrms.size(1) * sigma_tmp + 2];
      }
      wgts__params_pointwise[loop_ub] = sigma;
    }
    //  Compute wls
    wls_init(&wls_, us_, wgts__name_data, w_params_shared,
             wgts__params_pointwise, w_omit_rows, degree,
             interp0, xs_.size(0));
    if (wls_.rank < 0) {
      //  LAPACK error
      m2cErrMsgIdAndTxt("assemble_surf_range:badLapack",
                        "LAPACK error code %d for node %d", wls_.rank, (int)nid + 1);
    }
    if (!wls_.fullrank) {
      //  Not full rank
      rdtags[nid] = true;
      *fullrank = false;
    } else {
      int64_T m;
      ::coder::SizeType b_m;
      ::coder::SizeType b_n;
      //  Get centers and shape function Ns at centers
      b_dim = mesh_coords.size(1) - 1;
      //  # of nearby elements
      m = n2e_row_ptr[nid + 1] - n2e_row_ptr[nid];
      eids_.set_size(m);
      xs_.set_size(m, mesh_coords.size(1));
      Ns_.set_size(m);
      for (int64_T b_j{1L}; b_j <= m; b_j++) {
        uint64_T c;
        ::coder::SizeType eid;
        ::coder::SizeType leid;
        ::coder::SizeType npe;
        eid = n2e_col_ind[static_cast<::coder::SizeType>((n2e_row_ptr[nid] + b_j) - 1L) -
                          1] -
              1;
        eids_[b_j - 1] = eid + 1;
        c = mesh_teids[eid] & 255UL;
        leid = (mesh_teids[eid] >> 8) - 1;
        npe =
            mesh_elemtables
                [static_cast<::coder::SizeType>(
                     mesh_teids[n2e_col_ind[static_cast<::coder::SizeType>(
                                                (n2e_row_ptr[nid] + b_j) - 1L) -
                                            1] -
                                1] &
                     255UL) -
                 1]
                    .conn.size(1);
        Ns_[b_j - 1] =
            1.0 /
            static_cast<real_T>(
                mesh_elemtables[c - 1].conn.size(1));
        //  Compute localized center
        for (k = 0; k <= b_dim; k++) {
          xs_[k + xs_.size(1) * (b_j - 1)] = mesh_coords
              [k + mesh_coords.size(1) *
                       (mesh_elemtables[c - 1]
                            .conn[mesh_elemtables[c - 1]
                                      .conn.size(1) *
                                  leid] -
                        1)];
        }
        for (::coder::SizeType d_i{2}; d_i <= npe; d_i++) {
          for (k = 0; k <= b_dim; k++) {
            xs_[k + xs_.size(1) * (b_j - 1)] =
                xs_[k + xs_.size(1) * (b_j - 1)] +
                mesh_coords
                    [k + mesh_coords.size(1) *
                             (mesh_elemtables[c - 1].conn
                                  [(d_i +
                                    mesh_elemtables[c - 1]
                                            .conn.size(1) *
                                        leid) -
                                   1] -
                              1)];
          }
        }
        for (k = 0; k <= b_dim; k++) {
          xs_[k + xs_.size(1) * (b_j - 1)] =
              xs_[k + xs_.size(1) * (b_j - 1)] *
                  Ns_[b_j - 1] -
              mesh_coords[k + mesh_coords.size(1) * nid];
        }
      }
      //  Project centers onto tangent plane
      npoints = xs_.size(0) - 1;
      b_dim = xs_.size(1) - 2;
      us_.set_size(xs_.size(0), xs_.size(1) - 1);
      for (k = 0; k <= npoints; k++) {
        for (::coder::SizeType j{0}; j <= b_dim; j++) {
          us_[j + us_.size(1) * k] = 0.0;
        }
      }
      //  matrix-matrix, using mem-efficient loop
      for (::coder::SizeType ii{0}; ii <= npoints; ii++) {
        for (k = 0; k <= b_dim + 1; k++) {
          for (::coder::SizeType j{0}; j <= b_dim; j++) {
            us_[j + us_.size(1) * ii] =
                us_[j + us_.size(1) * ii] +
                xs_[k + xs_.size(1) * ii] * t_data[j + t_size_idx_1 * k];
          }
        }
      }
      //  Compute WLS fitting at centers
      wls_func(&wls_, us_, coeffs_);
      //  Compute local OSUS coefficients
      b_m = Ns_.size(0) - 1;
      b_n = coeffs_.size(0);
      //  WALF fitting
      for (k = 0; k < b_n; k++) {
        for (::coder::SizeType j{0}; j <= b_m; j++) {
          coeffs_[j + coeffs_.size(1) * k] =
              coeffs_[j + coeffs_.size(1) * k] * Ns_[j];
        }
      }
      //  Substract linear interp
      for (::coder::SizeType j{0}; j <= b_m; j++) {
        coeffs_[j] = coeffs_[j] - Ns_[j];
      }
      //  Add to global CRS
      b_n = eids_.size(0);
      b_m = nodes_.size(0);
      for (::coder::SizeType j{0}; j < b_n; j++) {
        ::coder::SizeType b_r;
        b_r = eids_[j];
        for (::coder::SizeType d_i{0}; d_i < b_m; d_i++) {
          int64_T b_k;
          b_k = 0L;
          if (b_r < row_ptr.size(0)) {
            //  Perform linear search
            b_k = m2cFind(col_ind, nodes_[d_i], row_ptr[b_r - 1],
                          row_ptr[b_r] - 1L);
          }
          val[b_k - 1] =
              val[b_k - 1] +
              coeffs_[j + coeffs_.size(1) * d_i];
        }
      }
    }
  }
}

static void
assemble_surf_task(RdiObject *rdi,
                   const ::coder::array<real_T, 2U> &mesh_coords,
                   const ::coder::array<ConnData, 1U> &mesh_elemtables,
                   const ::coder::array<uint64_T, 1U> &mesh_teids,
                   const ::coder::array<int64_T, 1U> &mesh_node2elems_row_ptr,
                   const ::coder::array<int32_T, 1U> &mesh_node2elems_col_ind,
                   const ::coder::array<Omp4mPart, 1U> &mesh_nodeparts,
                   const ::coder::array<int64_T, 1U> &stcl_row_ptr,
                   const ::coder::array<int32_T, 1U> &stcl_col_ind,
                   const ::coder::array<real_T, 2U> &nrms, boolean_T interp0)
{
  ::coder::SizeType lvl;
  boolean_T exitg1;
  //  kernel for task
  lvl = 0;
  exitg1 = false;
  while ((!exitg1) && (lvl <= mesh_nodeparts.size(0) - 1)) {
    ::coder::SizeType i;
#pragma omp single
    { // single
      i = mesh_nodeparts[lvl].nparts;
      for (::coder::SizeType part{0}; part < i; part++) {
#pragma omp task default(shared)
        { // task
          assemble_surf_range(
              rdi->A.row_ptr, rdi->A.col_ind, mesh_coords, mesh_elemtables,
              mesh_teids, stcl_row_ptr, stcl_col_ind, mesh_node2elems_row_ptr,
              mesh_node2elems_col_ind, nrms, rdi->degree, interp0,
              mesh_nodeparts[lvl].part_list, mesh_nodeparts[lvl].part_ptr[part],
              mesh_nodeparts[lvl].part_ptr[part + 1] - 1, rdi->A.val,
              rdi->rdtags, &rdi->fullrank);
        } // task
      }
    } // single
    //  implicit barrier
    if (mesh_nodeparts[lvl].shared_ents.size(0) != 0) {
#pragma omp single nowait
      { // single
        assemble_surf_range(rdi->A.row_ptr, rdi->A.col_ind, mesh_coords,
                            mesh_elemtables, mesh_teids, stcl_row_ptr,
                            stcl_col_ind, mesh_node2elems_row_ptr,
                            mesh_node2elems_col_ind, nrms, rdi->degree, interp0,
                            mesh_nodeparts[lvl].shared_ents, 1,
                            mesh_nodeparts[lvl].shared_ents.size(0), rdi->A.val,
                            rdi->rdtags, &rdi->fullrank);
      } // single
      exitg1 = true;
    } else {
      lvl++;
    }
  }
}

static ::coder::SizeType b_append_wlsmesh_kring(WlsMesh *mesh)
{
  static const char_T name[18]{'R', 'd', 'i', 'P', 'r', 'i', 'm', 'a', 'r',
                               'y', 'S', 't', 'e', 'n', 'c', 'i', 'l', 's'};
  ::coder::SizeType stclidx;
  stclidx = mesh->stencils.size(0) + 1;
  // Stencils - Object containing the stencils of vertices in a WlsMesh based on
  mesh->stencils.set_size(mesh->stencils.size(0) + 1);
  mesh->stencils[stclidx - 1].reflected.set_size(0);
  mesh->stencils[stclidx - 1].vidmap.set_size(0);
  mesh->stencils[stclidx - 1].ngbverts.row_ptr.set_size(0);
  mesh->stencils[stclidx - 1].ngbverts.ncols = 0;
  //  Construct A.row_ptr
  mesh->stencils[stclidx - 1].ngbverts.col_ind.set_size(0);
  mesh->stencils[stclidx - 1].ngbelems.row_ptr.set_size(0);
  mesh->stencils[stclidx - 1].ngbelems.ncols = 0;
  //  Construct A.row_ptr
  mesh->stencils[stclidx - 1].ngbelems.col_ind.set_size(0);
  mesh->stencils[stclidx - 1].name.set_size(1, 18);
  for (::coder::SizeType i{0}; i < 18; i++) {
    mesh->stencils[stclidx - 1].name[i] = name[i];
  }
  return stclidx;
}

static inline
void
b_assemble_body(RdiObject *rdi, const ::coder::array<real_T, 2U> &mesh_coords,
                const ::coder::array<ConnData, 1U> &mesh_elemtables,
                const ::coder::array<uint64_T, 1U> &mesh_teids,
                const ::coder::array<int64_T, 1U> &mesh_node2elems_row_ptr,
                const ::coder::array<int32_T, 1U> &mesh_node2elems_col_ind,
                const ::coder::array<Stencils, 1U> &mesh_stencils,
                const ::coder::array<Omp4mPart, 1U> &mesh_nodeparts,
                ::coder::SizeType stclid, boolean_T interp0, ::coder::SizeType varargin_2)
{
  if ((stclid != rdi->stclid) && (stclid != rdi->extstclid)) {
    m2cErrMsgIdAndTxt("assemble_body:badStencilID", "bad stencil ID (index) %d",
                      (int)stclid);
  }
  //  Set fullrank to be true upon input
  rdi->fullrank = true;
  if (mesh_stencils[stclid - 1].vidmap.size(0) != 0) {
    //  For extended stencils (resolving rank deficiencies), must be serial
    rdi->fullrank = b_assemble_body_range(
        rdi->A.row_ptr, rdi->A.col_ind, mesh_coords, mesh_elemtables,
        mesh_teids, mesh_stencils[stclid - 1].ngbverts.row_ptr,
        mesh_stencils[stclid - 1].ngbverts.col_ind, mesh_node2elems_row_ptr,
        mesh_node2elems_col_ind, rdi->degree, interp0,
        mesh_stencils[stclid - 1].vidmap,
        mesh_stencils[stclid - 1].vidmap.size(0), rdi->A.val, rdi->rdtags);
  } else if (mesh_nodeparts.size(0) == 0) {
    //  Must be serial
    rdi->fullrank = c_assemble_body_range(
        rdi->A.row_ptr, rdi->A.col_ind, mesh_coords, mesh_elemtables,
        mesh_teids, mesh_stencils[stclid - 1].ngbverts.row_ptr,
        mesh_stencils[stclid - 1].ngbverts.col_ind, mesh_node2elems_row_ptr,
        mesh_node2elems_col_ind, rdi->degree, interp0, mesh_coords.size(0),
        rdi->A.val, rdi->rdtags);
  } else {
    boolean_T m2cTryBlkErrFlag;
    ::coder::SizeType nthreads;
    if (varargin_2 <= 0) {
      nthreads = 1;
#ifdef _OPENMP
      nthreads = omp_get_max_threads();
#endif // _OPENMP
    } else {
      nthreads = varargin_2;
    }
    m2cTryBlkErrFlag = 0;
#pragma omp parallel default(shared) num_threads(nthreads)
    try { // try
      assemble_body_task(rdi, mesh_coords, mesh_elemtables, mesh_teids,
                         mesh_node2elems_row_ptr, mesh_node2elems_col_ind,
                         mesh_nodeparts,
                         mesh_stencils[stclid - 1].ngbverts.row_ptr,
                         mesh_stencils[stclid - 1].ngbverts.col_ind, interp0);
    } catch (const std::runtime_error &m2cExc) {
      m2cTryBlkErrFlag = 1;
      m2cPrintError("runtime_error %s\n", m2cExc.what());
    } catch (const std::logic_error &m2cExc) {
      m2cTryBlkErrFlag = 1;
      m2cPrintError("logic_error %s\n", m2cExc.what());
    } catch (const std::exception &m2cExc) {
      m2cTryBlkErrFlag = 1;
      m2cPrintError("exception %s\n", m2cExc.what());
    } catch (...) {
      m2cTryBlkErrFlag = 1;
      m2cPrintError("Unknown error detected from C++ exceptions\n");
      fflush(stderr);
    } // end try
    if ((int32_T)m2cTryBlkErrFlag != 0) {
      throw std::runtime_error("omp4m:runtimeErrorInThread");
    }
  }
}

static boolean_T b_assemble_body_range(
    const ::coder::array<int64_T, 1U> &row_ptr,
    const ::coder::array<int32_T, 1U> &col_ind,
    const ::coder::array<real_T, 2U> &mesh_coords,
    const ::coder::array<ConnData, 1U> &mesh_elemtables,
    const ::coder::array<uint64_T, 1U> &mesh_teids,
    const ::coder::array<int64_T, 1U> &stcl_row_ptr,
    const ::coder::array<int32_T, 1U> &stcl_col_ind,
    const ::coder::array<int64_T, 1U> &n2e_row_ptr,
    const ::coder::array<int32_T, 1U> &n2e_col_ind, ::coder::SizeType degree,
    boolean_T interp0, const ::coder::array<int32_T, 1U> &nrange, ::coder::SizeType iend,
    ::coder::array<real_T, 1U> &val, ::coder::array<boolean_T, 1U> &rdtags)
{
  static const char_T name[7]{'B', 'u', 'h', 'm', 'a', 'n', 'n'};
  ::coder::array<real_T, 2U> coeffs_;
  ::coder::array<real_T, 2U> w_params_pointwise;
  ::coder::array<real_T, 2U> xs_;
  ::coder::array<real_T, 1U> Ns_;
  ::coder::array<real_T, 1U> w_params_shared;
  ::coder::array<int32_T, 1U> eids_;
  ::coder::array<int32_T, 1U> nodes_;
  ::coder::array<boolean_T, 1U> w_omit_rows;
  b_WlsObject expl_temp;
  e_WlsObject wls_;
  ::coder::SizeType c_loop_ub;
  ::coder::SizeType dim;
  ::coder::SizeType loop_ub;
  char_T wgts__name_data[7];
  boolean_T fullrank;
  fullrank = true;
  f_WlsObject(degree, &expl_temp);
  wls_.runtimes.size[0] = expl_temp.runtimes.size[0];
  loop_ub = expl_temp.runtimes.size[0];
  if (loop_ub - 1 >= 0) {
    std::copy(&expl_temp.runtimes.data[0], &expl_temp.runtimes.data[loop_ub],
              &wls_.runtimes.data[0]);
  }
  wls_.QRt.set_size(expl_temp.QRt.size(0), expl_temp.QRt.size(1));
  loop_ub = expl_temp.QRt.size(1) * expl_temp.QRt.size(0);
  for (::coder::SizeType i{0}; i < loop_ub; i++) {
    wls_.QRt[i] = expl_temp.QRt[i];
  }
  wls_.rowmajor = expl_temp.rowmajor;
  wls_.work.set_size(expl_temp.work.size(0));
  loop_ub = expl_temp.work.size(0);
  for (::coder::SizeType i{0}; i < loop_ub; i++) {
    wls_.work[i] = expl_temp.work[i];
  }
  wls_.jpvt.set_size(expl_temp.jpvt.size(0));
  loop_ub = expl_temp.jpvt.size(0);
  for (::coder::SizeType i{0}; i < loop_ub; i++) {
    wls_.jpvt[i] = expl_temp.jpvt[i];
  }
  wls_.fullrank = expl_temp.fullrank;
  wls_.rank = expl_temp.rank;
  wls_.ncols = expl_temp.ncols;
  wls_.nrows = expl_temp.nrows;
  wls_.nevpnts = expl_temp.nevpnts;
  wls_.rhs.set_size(expl_temp.rhs.size(0), expl_temp.rhs.size(1));
  loop_ub = expl_temp.rhs.size(1) * expl_temp.rhs.size(0);
  for (::coder::SizeType i{0}; i < loop_ub; i++) {
    wls_.rhs[i] = expl_temp.rhs[i];
  }
  wls_.QR.set_size(expl_temp.QR.size(0), expl_temp.QR.size(1));
  loop_ub = expl_temp.QR.size(1) * expl_temp.QR.size(0);
  for (::coder::SizeType i{0}; i < loop_ub; i++) {
    wls_.QR[i] = expl_temp.QR[i];
  }
  wls_.V.set_size(expl_temp.V.size(0), expl_temp.V.size(1));
  loop_ub = expl_temp.V.size(1) * expl_temp.V.size(0);
  for (::coder::SizeType i{0}; i < loop_ub; i++) {
    wls_.V[i] = expl_temp.V[i];
  }
  wls_.hs_inv.set_size(1, expl_temp.hs_inv.size[1]);
  loop_ub = expl_temp.hs_inv.size[1];
  for (::coder::SizeType i{0}; i < loop_ub; i++) {
    wls_.hs_inv[i] = expl_temp.hs_inv.data[i];
  }
  wls_.rweights.set_size(expl_temp.rweights.size(0));
  loop_ub = expl_temp.rweights.size(0);
  for (::coder::SizeType i{0}; i < loop_ub; i++) {
    wls_.rweights[i] = expl_temp.rweights[i];
  }
  wls_.origin.size[1] = expl_temp.origin.size[1];
  wls_.origin.size[0] = 1;
  loop_ub = expl_temp.origin.size[1];
  if (loop_ub - 1 >= 0) {
    std::copy(&expl_temp.origin.data[0], &expl_temp.origin.data[loop_ub],
              &wls_.origin.data[0]);
  }
  wls_.us.set_size(expl_temp.us.size(0), expl_temp.us.size(1));
  loop_ub = expl_temp.us.size(0);
  for (::coder::SizeType i{0}; i < loop_ub; i++) {
    ::coder::SizeType b_loop_ub;
    b_loop_ub = expl_temp.us.size(1);
    for (::coder::SizeType i1{0}; i1 < b_loop_ub; i1++) {
      wls_.us[i1 + wls_.us.size(1) * i] =
          expl_temp.us[i1 + expl_temp.us.size(1) * i];
    }
  }
  wls_.stride = expl_temp.stride;
  wls_.interp0 = expl_temp.interp0;
  wls_.unimono = expl_temp.unimono;
  wls_.order = expl_temp.order;
  wls_.degree = expl_temp.degree;
  wls_.nstpnts = expl_temp.nstpnts;
  w_params_shared.set_size(0);
  w_params_pointwise.set_size(::coder::SizeType(0), ::coder::SizeType(0));
  w_omit_rows.set_size(0);
  for (::coder::SizeType i{0}; i < 7; i++) {
    wgts__name_data[i] = name[i];
  }
  // Local buffers with patterns "(\w+)?params_pointwise" are legitimate
  if (iend >= 1) {
    dim = mesh_coords.size(1);
    c_loop_ub = mesh_coords.size(1);
  }
  for (::coder::SizeType b_i{1}; b_i <= iend; b_i++) {
    int64_T c_i;
    int64_T n;
    int64_T n_tmp;
    ::coder::SizeType k;
    ::coder::SizeType nid;
    if (nrange.size(0) == 0) {
      nid = b_i;
    } else {
      nid = nrange[b_i - 1];
    }
    //  Fetch local data
    n_tmp = stcl_row_ptr[b_i - 1];
    n = stcl_row_ptr[b_i] - n_tmp;
    nodes_.set_size((n + 1L));
    xs_.set_size((n + 1L), mesh_coords.size(1));
    nodes_[0] = nid;
    for (::coder::SizeType i{0}; i < c_loop_ub; i++) {
      xs_[i] = 0.0;
    }
    for (c_i = 2L; c_i - 1L <= n; c_i++) {
      k = stcl_col_ind[static_cast<::coder::SizeType>((n_tmp + c_i) - 2L) - 1];
      nodes_[c_i - 1] = k;
      for (::coder::SizeType j{0}; j < dim; j++) {
        xs_[j + xs_.size(1) * (c_i - 1)] =
            mesh_coords[j + mesh_coords.size(1) * (k - 1)] -
            mesh_coords[j + mesh_coords.size(1) * (nid - 1)];
      }
    }
    //  Compute wls
    wls_init(&wls_, xs_, wgts__name_data, w_params_shared, w_params_pointwise,
             w_omit_rows, degree, interp0, xs_.size(0));
    if (wls_.rank < 0) {
      //  LAPACK error
      m2cErrMsgIdAndTxt("assemble_body_range:badLapack",
                        "LAPACK error code %d for node %d", wls_.rank, (int)nid);
    }
    if (!wls_.fullrank) {
      //  Not full rank
      rdtags[nid - 1] = true;
      fullrank = false;
    } else {
      int64_T m;
      ::coder::SizeType b_dim;
      ::coder::SizeType b_m;
      ::coder::SizeType b_n;
      //  Get centers and shape function Ns at centers
      b_dim = mesh_coords.size(1) - 1;
      //  # of nearby elements
      n_tmp = n2e_row_ptr[nid - 1];
      m = n2e_row_ptr[nid] - n_tmp;
      eids_.set_size(m);
      xs_.set_size(m, mesh_coords.size(1));
      Ns_.set_size(m);
      for (int64_T b_j{1L}; b_j <= m; b_j++) {
        uint64_T c;
        ::coder::SizeType eid;
        ::coder::SizeType leid;
        ::coder::SizeType npe;
        eid = n2e_col_ind[static_cast<::coder::SizeType>((n_tmp + b_j) - 1L) - 1] - 1;
        eids_[b_j - 1] = eid + 1;
        c = mesh_teids[eid] & 255UL;
        leid = (mesh_teids[eid] >> 8) - 1;
        npe = mesh_elemtables
                  [static_cast<::coder::SizeType>(
                       mesh_teids[n2e_col_ind[static_cast<::coder::SizeType>(
                                                  (n2e_row_ptr[nid - 1] + b_j) -
                                                  1L) -
                                              1] -
                                  1] &
                       255UL) -
                   1]
                      .conn.size(1);
        Ns_[b_j - 1] =
            1.0 /
            static_cast<real_T>(
                mesh_elemtables[c - 1].conn.size(1));
        //  Compute localized center
        for (k = 0; k <= b_dim; k++) {
          xs_[k + xs_.size(1) * (b_j - 1)] = mesh_coords
              [k + mesh_coords.size(1) *
                       (mesh_elemtables[c - 1]
                            .conn[mesh_elemtables[c - 1]
                                      .conn.size(1) *
                                  leid] -
                        1)];
        }
        for (::coder::SizeType d_i{2}; d_i <= npe; d_i++) {
          for (k = 0; k <= b_dim; k++) {
            xs_[k + xs_.size(1) * (b_j - 1)] =
                xs_[k + xs_.size(1) * (b_j - 1)] +
                mesh_coords
                    [k + mesh_coords.size(1) *
                             (mesh_elemtables[c - 1].conn
                                  [(d_i +
                                    mesh_elemtables[c - 1]
                                            .conn.size(1) *
                                        leid) -
                                   1] -
                              1)];
          }
        }
        for (k = 0; k <= b_dim; k++) {
          xs_[k + xs_.size(1) * (b_j - 1)] =
              xs_[k + xs_.size(1) * (b_j - 1)] *
                  Ns_[b_j - 1] -
              mesh_coords[k + mesh_coords.size(1) * (nid - 1)];
        }
      }
      //  Compute WLS fitting at centers
      wls_func(&wls_, xs_, coeffs_);
      //  Compute local OSUS coefficients
      b_m = Ns_.size(0) - 1;
      b_n = coeffs_.size(0);
      //  WALF fitting
      for (k = 0; k < b_n; k++) {
        for (::coder::SizeType j{0}; j <= b_m; j++) {
          coeffs_[j + coeffs_.size(1) * k] =
              coeffs_[j + coeffs_.size(1) * k] * Ns_[j];
        }
      }
      //  Substract linear interp
      for (::coder::SizeType j{0}; j <= b_m; j++) {
        coeffs_[j] = coeffs_[j] - Ns_[j];
      }
      //  Add to global CRS
      b_n = eids_.size(0);
      b_m = nodes_.size(0);
      for (::coder::SizeType j{0}; j < b_n; j++) {
        ::coder::SizeType r;
        r = eids_[j];
        for (::coder::SizeType d_i{0}; d_i < b_m; d_i++) {
          int64_T b_k;
          b_k = 0L;
          if (r < row_ptr.size(0)) {
            int64_T bounds_idx_1;
            //  Perform linear search
            bounds_idx_1 = row_ptr[r] - 1L;
            //  Perform linear search
            c_i = row_ptr[r - 1];
            while (c_i <= bounds_idx_1) {
              if (col_ind[c_i - 1] == nodes_[d_i]) {
                b_k = c_i;
                c_i = bounds_idx_1 + 1L;
              } else {
                c_i++;
              }
            }
          }
          val[b_k - 1] =
              val[b_k - 1] +
              coeffs_[j + coeffs_.size(1) * d_i];
        }
      }
    }
  }
  return fullrank;
}

static void
b_assemble_surf(RdiObject *rdi, const ::coder::array<real_T, 2U> &mesh_coords,
                const ::coder::array<ConnData, 1U> &mesh_elemtables,
                const ::coder::array<uint64_T, 1U> &mesh_teids,
                const ::coder::array<NormalsData, 1U> &mesh_nrmstables,
                const ::coder::array<int64_T, 1U> &mesh_node2elems_row_ptr,
                const ::coder::array<int32_T, 1U> &mesh_node2elems_col_ind,
                const ::coder::array<Stencils, 1U> &mesh_stencils,
                const ::coder::array<Omp4mPart, 1U> &mesh_nodeparts,
                ::coder::SizeType stclid, boolean_T interp0, ::coder::SizeType varargin_2)
{
  ::coder::SizeType nrmid;
  if ((stclid != rdi->stclid) && (stclid != rdi->extstclid)) {
    m2cErrMsgIdAndTxt("assemble_surf:badStencilID", "bad stencil ID (index) %d",
                      (int)stclid);
  }
  //  Set fullrank to be true upon input
  rdi->fullrank = true;
  //  Overcome a "bug" in Coder regarding creating copies of normal data
  nrmid = rdi->nrmid;
  if (mesh_stencils[stclid - 1].vidmap.size(0) != 0) {
    //  For extended stencils (resolving rank deficiencies), must be serial
    rdi->fullrank = b_assemble_surf_range(
        rdi->A.row_ptr, rdi->A.col_ind, mesh_coords, mesh_elemtables,
        mesh_teids, mesh_stencils[stclid - 1].ngbverts.row_ptr,
        mesh_stencils[stclid - 1].ngbverts.col_ind, mesh_node2elems_row_ptr,
        mesh_node2elems_col_ind, mesh_nrmstables[nrmid - 1].normals,
        rdi->degree, interp0, mesh_stencils[stclid - 1].vidmap,
        mesh_stencils[stclid - 1].vidmap.size(0), rdi->A.val, rdi->rdtags);
  } else if (mesh_nodeparts.size(0) == 0) {
    //  Must be serial
    rdi->fullrank = c_assemble_surf_range(
        rdi->A.row_ptr, rdi->A.col_ind, mesh_coords, mesh_elemtables,
        mesh_teids, mesh_stencils[stclid - 1].ngbverts.row_ptr,
        mesh_stencils[stclid - 1].ngbverts.col_ind, mesh_node2elems_row_ptr,
        mesh_node2elems_col_ind, mesh_nrmstables[nrmid - 1].normals,
        rdi->degree, interp0, mesh_coords.size(0), rdi->A.val, rdi->rdtags);
  } else {
    boolean_T m2cTryBlkErrFlag;
    ::coder::SizeType nthreads;
    if (varargin_2 <= 0) {
      nthreads = 1;
#ifdef _OPENMP
      nthreads = omp_get_max_threads();
#endif // _OPENMP
    } else {
      nthreads = varargin_2;
    }
    m2cTryBlkErrFlag = 0;
#pragma omp parallel default(shared) num_threads(nthreads)
    try { // try
      assemble_surf_task(rdi, mesh_coords, mesh_elemtables, mesh_teids,
                         mesh_node2elems_row_ptr, mesh_node2elems_col_ind,
                         mesh_nodeparts,
                         mesh_stencils[stclid - 1].ngbverts.row_ptr,
                         mesh_stencils[stclid - 1].ngbverts.col_ind,
                         mesh_nrmstables[nrmid - 1].normals, interp0);
    } catch (const std::runtime_error &m2cExc) {
      m2cTryBlkErrFlag = 1;
      m2cPrintError("runtime_error %s\n", m2cExc.what());
    } catch (const std::logic_error &m2cExc) {
      m2cTryBlkErrFlag = 1;
      m2cPrintError("logic_error %s\n", m2cExc.what());
    } catch (const std::exception &m2cExc) {
      m2cTryBlkErrFlag = 1;
      m2cPrintError("exception %s\n", m2cExc.what());
    } catch (...) {
      m2cTryBlkErrFlag = 1;
      m2cPrintError("Unknown error detected from C++ exceptions\n");
      fflush(stderr);
    } // end try
    if ((int32_T)m2cTryBlkErrFlag != 0) {
      throw std::runtime_error("omp4m:runtimeErrorInThread");
    }
  }
}

static boolean_T b_assemble_surf_range(
    const ::coder::array<int64_T, 1U> &row_ptr,
    const ::coder::array<int32_T, 1U> &col_ind,
    const ::coder::array<real_T, 2U> &mesh_coords,
    const ::coder::array<ConnData, 1U> &mesh_elemtables,
    const ::coder::array<uint64_T, 1U> &mesh_teids,
    const ::coder::array<int64_T, 1U> &stcl_row_ptr,
    const ::coder::array<int32_T, 1U> &stcl_col_ind,
    const ::coder::array<int64_T, 1U> &n2e_row_ptr,
    const ::coder::array<int32_T, 1U> &n2e_col_ind,
    const ::coder::array<real_T, 2U> &nrms, ::coder::SizeType degree, boolean_T interp0,
    const ::coder::array<int32_T, 1U> &nrange, ::coder::SizeType iend,
    ::coder::array<real_T, 1U> &val, ::coder::array<boolean_T, 1U> &rdtags)
{
  static const char_T name[7]{'B', 'u', 'h', 'm', 'a', 'n', 'n'};
  ::coder::array<real_T, 2U> coeffs_;
  ::coder::array<real_T, 2U> us_;
  ::coder::array<real_T, 2U> wgts__params_pointwise;
  ::coder::array<real_T, 2U> xs_;
  ::coder::array<real_T, 1U> Ns_;
  ::coder::array<real_T, 1U> w_params_shared;
  ::coder::array<int32_T, 1U> eids_;
  ::coder::array<int32_T, 1U> nodes_;
  ::coder::array<boolean_T, 1U> w_omit_rows;
  b_WlsObject r;
  c_WlsObject expl_temp;
  d_WlsObject b_expl_temp;
  e_WlsObject wls_;
  real_T t_data[6];
  ::coder::SizeType b_loop_ub;
  ::coder::SizeType dim;
  ::coder::SizeType i;
  ::coder::SizeType i1;
  ::coder::SizeType loop_ub;
  ::coder::SizeType t_size_idx_1;
  char_T wgts__name_data[7];
  boolean_T fullrank;
  fullrank = true;
  f_WlsObject(degree, &r);
  majorityTransform(&r, &expl_temp);
  b_expl_temp.us.set_size(expl_temp.us.size(0), expl_temp.us.size(1));
  loop_ub = expl_temp.us.size(1) * expl_temp.us.size(0);
  for (i = 0; i < loop_ub; i++) {
    b_expl_temp.us[i] = expl_temp.us[i];
  }
  b_expl_temp.hs_inv.set_size(expl_temp.hs_inv.size[0], 1);
  loop_ub = expl_temp.hs_inv.size[0];
  for (i = 0; i < loop_ub; i++) {
    b_expl_temp.hs_inv[i] = expl_temp.hs_inv.data[i];
  }
  b_expl_temp.runtimes.size[0] = expl_temp.runtimes.size[0];
  loop_ub = expl_temp.runtimes.size[0];
  if (loop_ub - 1 >= 0) {
    std::copy(&expl_temp.runtimes.data[0], &expl_temp.runtimes.data[loop_ub],
              &b_expl_temp.runtimes.data[0]);
  }
  b_expl_temp.QRt.set_size(expl_temp.QRt.size(0), expl_temp.QRt.size(1));
  loop_ub = expl_temp.QRt.size(1) * expl_temp.QRt.size(0);
  for (i = 0; i < loop_ub; i++) {
    b_expl_temp.QRt[i] = expl_temp.QRt[i];
  }
  b_expl_temp.rowmajor = expl_temp.rowmajor;
  b_expl_temp.work.set_size(expl_temp.work.size(0));
  loop_ub = expl_temp.work.size(0);
  for (i = 0; i < loop_ub; i++) {
    b_expl_temp.work[i] = expl_temp.work[i];
  }
  b_expl_temp.jpvt.set_size(expl_temp.jpvt.size(0));
  loop_ub = expl_temp.jpvt.size(0);
  for (i = 0; i < loop_ub; i++) {
    b_expl_temp.jpvt[i] = expl_temp.jpvt[i];
  }
  b_expl_temp.fullrank = expl_temp.fullrank;
  b_expl_temp.rank = expl_temp.rank;
  b_expl_temp.ncols = expl_temp.ncols;
  b_expl_temp.nrows = expl_temp.nrows;
  b_expl_temp.nevpnts = expl_temp.nevpnts;
  b_expl_temp.rhs.set_size(expl_temp.rhs.size(0), expl_temp.rhs.size(1));
  loop_ub = expl_temp.rhs.size(1) * expl_temp.rhs.size(0);
  for (i = 0; i < loop_ub; i++) {
    b_expl_temp.rhs[i] = expl_temp.rhs[i];
  }
  b_expl_temp.QR.set_size(expl_temp.QR.size(0), expl_temp.QR.size(1));
  loop_ub = expl_temp.QR.size(1) * expl_temp.QR.size(0);
  for (i = 0; i < loop_ub; i++) {
    b_expl_temp.QR[i] = expl_temp.QR[i];
  }
  b_expl_temp.V.set_size(expl_temp.V.size(0), expl_temp.V.size(1));
  loop_ub = expl_temp.V.size(1) * expl_temp.V.size(0);
  for (i = 0; i < loop_ub; i++) {
    b_expl_temp.V[i] = expl_temp.V[i];
  }
  b_expl_temp.rweights.set_size(expl_temp.rweights.size(0));
  loop_ub = expl_temp.rweights.size(0);
  for (i = 0; i < loop_ub; i++) {
    b_expl_temp.rweights[i] = expl_temp.rweights[i];
  }
  b_expl_temp.origin.size[1] = 1;
  b_expl_temp.origin.size[0] = expl_temp.origin.size[0];
  loop_ub = expl_temp.origin.size[0];
  if (loop_ub - 1 >= 0) {
    std::copy(&expl_temp.origin.data[0], &expl_temp.origin.data[loop_ub],
              &b_expl_temp.origin.data[0]);
  }
  b_expl_temp.stride = expl_temp.stride;
  b_expl_temp.interp0 = expl_temp.interp0;
  b_expl_temp.unimono = expl_temp.unimono;
  b_expl_temp.order = expl_temp.order;
  b_expl_temp.degree = expl_temp.degree;
  b_expl_temp.nstpnts = expl_temp.nstpnts;
  majorityTransform(&b_expl_temp, &wls_);
  w_params_shared.set_size(0);
  w_omit_rows.set_size(0);
  for (i = 0; i < 7; i++) {
    wgts__name_data[i] = name[i];
  }
  // Local buffers with patterns "wgts__params_pointwise" are legitimate
  if (iend >= 1) {
    dim = mesh_coords.size(1);
    b_loop_ub = mesh_coords.size(1);
    t_size_idx_1 = nrms.size(1) - 1;
    i1 = nrms.size(1);
  }
  for (::coder::SizeType b_i{1}; b_i <= iend; b_i++) {
    int64_T n;
    int64_T n_tmp;
    ::coder::SizeType b_dim;
    ::coder::SizeType b_npoints;
    ::coder::SizeType k;
    ::coder::SizeType nid;
    ::coder::SizeType npoints;
    boolean_T b;
    boolean_T b1;
    if (nrange.size(0) == 0) {
      nid = b_i - 1;
    } else {
      nid = nrange[b_i - 1] - 1;
    }
    //  Fetch local data
    n_tmp = stcl_row_ptr[b_i - 1];
    n = stcl_row_ptr[b_i] - n_tmp;
    nodes_.set_size((n + 1L));
    xs_.set_size((n + 1L), mesh_coords.size(1));
    nodes_[0] = nid + 1;
    for (i = 0; i < b_loop_ub; i++) {
      xs_[i] = 0.0;
    }
    for (int64_T c_i{2L}; c_i - 1L <= n; c_i++) {
      k = stcl_col_ind[static_cast<::coder::SizeType>((n_tmp + c_i) - 2L) - 1];
      nodes_[c_i - 1] = k;
      for (::coder::SizeType j{0}; j < dim; j++) {
        xs_[j + xs_.size(1) * (c_i - 1)] =
            mesh_coords[j + mesh_coords.size(1) * (k - 1)] -
            mesh_coords[j + mesh_coords.size(1) * nid];
      }
    }
    npoints = xs_.size(0);
    //  Get normal direction
    if (i1 == 2) {
      t_data[0] = -nrms[nrms.size(1) * nid + 1];
      t_data[t_size_idx_1] = nrms[nrms.size(1) * nid];
    } else {
      real_T a;
      real_T a_tmp;
      a_tmp = std::abs(nrms[nrms.size(1) * nid]);
      if ((a_tmp > std::abs(nrms[nrms.size(1) * nid + 1])) &&
          (a_tmp > std::abs(nrms[nrms.size(1) * nid + 2]))) {
        t_data[0] = -nrms[nrms.size(1) * nid] * nrms[nrms.size(1) * nid + 1];
        t_data[t_size_idx_1] =
            1.0 - nrms[nrms.size(1) * nid + 1] * nrms[nrms.size(1) * nid + 1];
        t_data[t_size_idx_1 * 2] =
            -nrms[nrms.size(1) * nid + 1] * nrms[nrms.size(1) * nid + 2];
      } else {
        t_data[0] = 1.0 - nrms[nrms.size(1) * nid] * nrms[nrms.size(1) * nid];
        t_data[t_size_idx_1] =
            -nrms[nrms.size(1) * nid] * nrms[nrms.size(1) * nid + 1];
        t_data[t_size_idx_1 * 2] =
            -nrms[nrms.size(1) * nid] * nrms[nrms.size(1) * nid + 2];
      }
      a_tmp = t_data[t_size_idx_1 * 2];
      a = std::sqrt((t_data[0] * t_data[0] +
                     t_data[t_size_idx_1] * t_data[t_size_idx_1]) +
                    a_tmp * a_tmp);
      t_data[0] /= a;
      t_data[t_size_idx_1] /= a;
      t_data[t_size_idx_1 * 2] /= a;
      //  cross
      t_data[1] = t_data[t_size_idx_1 * 2] * nrms[nrms.size(1) * nid + 1] -
                  t_data[t_size_idx_1] * nrms[nrms.size(1) * nid + 2];
      t_data[t_size_idx_1 + 1] =
          t_data[0] * nrms[nrms.size(1) * nid + 2] -
          nrms[nrms.size(1) * nid] * t_data[t_size_idx_1 * 2];
      t_data[t_size_idx_1 * 2 + 1] =
          nrms[nrms.size(1) * nid] * t_data[t_size_idx_1] -
          t_data[0] * nrms[nrms.size(1) * nid + 1];
    }
    //  Project onto tangent plane
    b_npoints = xs_.size(0) - 1;
    b_dim = xs_.size(1) - 2;
    us_.set_size(xs_.size(0), xs_.size(1) - 1);
    for (k = 0; k <= b_npoints; k++) {
      for (::coder::SizeType j{0}; j <= b_dim; j++) {
        us_[j + us_.size(1) * k] = 0.0;
      }
    }
    //  matrix-matrix, using mem-efficient loop
    for (::coder::SizeType ii{0}; ii <= b_npoints; ii++) {
      for (k = 0; k <= b_dim + 1; k++) {
        for (::coder::SizeType j{0}; j <= b_dim; j++) {
          us_[j + us_.size(1) * ii] =
              us_[j + us_.size(1) * ii] +
              xs_[k + xs_.size(1) * ii] * t_data[j + t_size_idx_1 * k];
        }
      }
    }
    //  Compute normal matches as additional weights to WLS
    wgts__params_pointwise.set_size(xs_.size(0), 1);
    wgts__params_pointwise[0] = 1.0;
    b = true;
    i = wgts__params_pointwise.size(0);
    b1 = wgts__params_pointwise.size(0) <= 0;
    loop_ub = 0;
    for (::coder::SizeType j{2}; j <= npoints; j++) {
      real_T sigma;
      ::coder::SizeType sigma_tmp;
      if (b1 || (j - 1 >= i)) {
        loop_ub = 0;
        b = true;
      } else if (b) {
        b = false;
        loop_ub = (j - 1) % wgts__params_pointwise.size(0) +
                  (j - 1) / wgts__params_pointwise.size(0);
      } else if (loop_ub > 2147483646) {
        loop_ub = (j - 1) % wgts__params_pointwise.size(0) +
                  (j - 1) / wgts__params_pointwise.size(0);
      } else {
        loop_ub++;
        if (loop_ub > wgts__params_pointwise.size(0) - 1) {
          loop_ub = (loop_ub - wgts__params_pointwise.size(0)) + 1;
        }
      }
      sigma_tmp = nodes_[j - 1] - 1;
      sigma = nrms[nrms.size(1) * nid] * nrms[nrms.size(1) * sigma_tmp] +
              nrms[nrms.size(1) * nid + 1] * nrms[nrms.size(1) * sigma_tmp + 1];
      if (nrms.size(1) > 2) {
        sigma +=
            nrms[nrms.size(1) * nid + 2] * nrms[nrms.size(1) * sigma_tmp + 2];
      }
      wgts__params_pointwise[loop_ub] = sigma;
    }
    //  Compute wls
    wls_init(&wls_, us_, wgts__name_data, w_params_shared,
             wgts__params_pointwise, w_omit_rows, degree,
             interp0, xs_.size(0));
    if (wls_.rank < 0) {
      //  LAPACK error
      m2cErrMsgIdAndTxt("assemble_surf_range:badLapack",
                        "LAPACK error code %d for node %d", wls_.rank, (int)nid + 1);
    }
    if (!wls_.fullrank) {
      //  Not full rank
      rdtags[nid] = true;
      fullrank = false;
    } else {
      int64_T m;
      ::coder::SizeType b_m;
      ::coder::SizeType b_n;
      //  Get centers and shape function Ns at centers
      b_dim = mesh_coords.size(1) - 1;
      //  # of nearby elements
      m = n2e_row_ptr[nid + 1] - n2e_row_ptr[nid];
      eids_.set_size(m);
      xs_.set_size(m, mesh_coords.size(1));
      Ns_.set_size(m);
      for (int64_T b_j{1L}; b_j <= m; b_j++) {
        uint64_T c;
        ::coder::SizeType eid;
        ::coder::SizeType leid;
        ::coder::SizeType npe;
        eid = n2e_col_ind[static_cast<::coder::SizeType>((n2e_row_ptr[nid] + b_j) - 1L) -
                          1] -
              1;
        eids_[b_j - 1] = eid + 1;
        c = mesh_teids[eid] & 255UL;
        leid = (mesh_teids[eid] >> 8) - 1;
        npe =
            mesh_elemtables
                [static_cast<::coder::SizeType>(
                     mesh_teids[n2e_col_ind[static_cast<::coder::SizeType>(
                                                (n2e_row_ptr[nid] + b_j) - 1L) -
                                            1] -
                                1] &
                     255UL) -
                 1]
                    .conn.size(1);
        Ns_[b_j - 1] =
            1.0 /
            static_cast<real_T>(
                mesh_elemtables[c - 1].conn.size(1));
        //  Compute localized center
        for (k = 0; k <= b_dim; k++) {
          xs_[k + xs_.size(1) * (b_j - 1)] = mesh_coords
              [k + mesh_coords.size(1) *
                       (mesh_elemtables[c - 1]
                            .conn[mesh_elemtables[c - 1]
                                      .conn.size(1) *
                                  leid] -
                        1)];
        }
        for (::coder::SizeType d_i{2}; d_i <= npe; d_i++) {
          for (k = 0; k <= b_dim; k++) {
            xs_[k + xs_.size(1) * (b_j - 1)] =
                xs_[k + xs_.size(1) * (b_j - 1)] +
                mesh_coords
                    [k + mesh_coords.size(1) *
                             (mesh_elemtables[c - 1].conn
                                  [(d_i +
                                    mesh_elemtables[c - 1]
                                            .conn.size(1) *
                                        leid) -
                                   1] -
                              1)];
          }
        }
        for (k = 0; k <= b_dim; k++) {
          xs_[k + xs_.size(1) * (b_j - 1)] =
              xs_[k + xs_.size(1) * (b_j - 1)] *
                  Ns_[b_j - 1] -
              mesh_coords[k + mesh_coords.size(1) * nid];
        }
      }
      //  Project centers onto tangent plane
      npoints = xs_.size(0) - 1;
      b_dim = xs_.size(1) - 2;
      us_.set_size(xs_.size(0), xs_.size(1) - 1);
      for (k = 0; k <= npoints; k++) {
        for (::coder::SizeType j{0}; j <= b_dim; j++) {
          us_[j + us_.size(1) * k] = 0.0;
        }
      }
      //  matrix-matrix, using mem-efficient loop
      for (::coder::SizeType ii{0}; ii <= npoints; ii++) {
        for (k = 0; k <= b_dim + 1; k++) {
          for (::coder::SizeType j{0}; j <= b_dim; j++) {
            us_[j + us_.size(1) * ii] =
                us_[j + us_.size(1) * ii] +
                xs_[k + xs_.size(1) * ii] * t_data[j + t_size_idx_1 * k];
          }
        }
      }
      //  Compute WLS fitting at centers
      wls_func(&wls_, us_, coeffs_);
      //  Compute local OSUS coefficients
      b_m = Ns_.size(0) - 1;
      b_n = coeffs_.size(0);
      //  WALF fitting
      for (k = 0; k < b_n; k++) {
        for (::coder::SizeType j{0}; j <= b_m; j++) {
          coeffs_[j + coeffs_.size(1) * k] =
              coeffs_[j + coeffs_.size(1) * k] * Ns_[j];
        }
      }
      //  Substract linear interp
      for (::coder::SizeType j{0}; j <= b_m; j++) {
        coeffs_[j] = coeffs_[j] - Ns_[j];
      }
      //  Add to global CRS
      b_n = eids_.size(0);
      b_m = nodes_.size(0);
      for (::coder::SizeType j{0}; j < b_n; j++) {
        ::coder::SizeType b_r;
        b_r = eids_[j];
        for (::coder::SizeType d_i{0}; d_i < b_m; d_i++) {
          int64_T b_k;
          b_k = 0L;
          if (b_r < row_ptr.size(0)) {
            //  Perform linear search
            b_k = m2cFind(col_ind, nodes_[d_i], row_ptr[b_r - 1],
                          row_ptr[b_r] - 1L);
          }
          val[b_k - 1] =
              val[b_k - 1] +
              coeffs_[j + coeffs_.size(1) * d_i];
        }
      }
    }
  }
  return fullrank;
}

static void b_compute_stencils_1d(WlsMesh *mesh, ::coder::SizeType stclidx,
                                  const real_T krings[3])
{
  int32_T bounds[2];
  bounds[0] = 0;
  //  set default upper bound
  bounds[1] = static_cast<::coder::SizeType>(std::round(2.0 * krings[0]));
  //  Ensure that max is no smaller than corresponding min
  if (bounds[1] <= 0) {
    bounds[1] = 0;
  }
  //  Parallel
  b_compute_stencils_kernel_1d(mesh, stclidx, krings[0], bounds);
}

static void b_compute_stencils_2d(WlsMesh *mesh, ::coder::SizeType stclidx,
                                  const real_T krings[3])
{
  int32_T bounds[2];
  bounds[0] = 0;
  //  set default upper bound
  bounds[1] = static_cast<::coder::SizeType>(
      std::round(12.0 * (krings[0] + 0.5) * (krings[0] + 1.0)));
  //  Ensure that max is no smaller than corresponding min
  if (bounds[1] <= 0) {
    bounds[1] = 0;
  }
  //  Parallel
  b_compute_stencils_kernel_2d(mesh, stclidx, krings[0], bounds);
}

static void b_compute_stencils_kernel_1d(WlsMesh *mesh, ::coder::SizeType stclidx,
                                         real_T ring, const int32_T bounds[2])
{
  ::coder::array<uint64_T, 1U> hebuf_;
  ::coder::array<int32_T, 1U> ngbfs_;
  ::coder::array<int32_T, 1U> ngbvs_;
  ::coder::array<boolean_T, 1U> ftags_;
  ::coder::array<boolean_T, 1U> vtags_;
  ::coder::SizeType iend;
  ::coder::SizeType istart;
  ::coder::SizeType n;
  int32_T nfaces;
  ::coder::SizeType nthreads;
  int32_T nverts;
  ::coder::SizeType u0;
  ::coder::SizeType u1;
  boolean_T hadoverflow;
  boolean_T overflow;
  boolean_T reflected;
  //  Determine total number of vertices
  if (mesh->stencils[stclidx - 1].vidmap.size(0) == 0) {
    n = mesh->coords.size(0);
  } else {
    n = mesh->stencils[stclidx - 1].vidmap.size(0);
  }
  //  Determine partitioning
  nthreads = 1;
#ifdef _OPENMP
  nthreads = omp_get_num_threads();
#endif // _OPENMP
  if (nthreads == 1) {
    istart = 1;
    iend = n;
  } else {
    ::coder::SizeType b_remainder;
    ::coder::SizeType chunk;
    ::coder::SizeType threadID;
    threadID = 0;
#ifdef _OPENMP
    threadID = omp_get_thread_num();
#endif // _OPENMP
    chunk = n / nthreads;
    b_remainder = n - nthreads * chunk;
    u1 = threadID;
    if (b_remainder <= threadID) {
      u1 = b_remainder;
    }
    istart = (threadID * chunk + u1) + 1;
    iend = ((istart + chunk) + (threadID < b_remainder)) - 1;
  }
#pragma omp single
  { // single
    mesh->stencils[stclidx - 1].reflected.set_size(n);
    mesh->stencils[stclidx - 1].ngbverts.row_ptr.set_size(n + 1);
    mesh->stencils[stclidx - 1].ngbverts.row_ptr[0] = 1L;
    mesh->stencils[stclidx - 1].ngbverts.ncols = mesh->coords.size(0);
    mesh->stencils[stclidx - 1].ngbelems.row_ptr.set_size(n + 1);
    mesh->stencils[stclidx - 1].ngbelems.row_ptr[0] = 1L;
    mesh->stencils[stclidx - 1].ngbelems.ncols = mesh->teids.size(0);
  } // single
  //  Assemble rowptrs
  for (::coder::SizeType i{istart}; i <= iend; i++) {
    //  Get vertex ID
    mesh->stencils[stclidx - 1].ngbverts.row_ptr[i] = bounds[1];
    mesh->stencils[stclidx - 1].ngbelems.row_ptr[i] = bounds[1] << 1;
  }
#pragma omp barrier
#pragma omp single nowait
  { // single
    for (::coder::SizeType i{0}; i < n; i++) {
      mesh->stencils[stclidx - 1].ngbverts.row_ptr[i + 1] =
          mesh->stencils[stclidx - 1].ngbverts.row_ptr[i + 1] +
          mesh->stencils[stclidx - 1].ngbverts.row_ptr[i];
    }
  } // single
#pragma omp single
  { // single
    for (::coder::SizeType i{0}; i < n; i++) {
      mesh->stencils[stclidx - 1].ngbelems.row_ptr[i + 1] =
          mesh->stencils[stclidx - 1].ngbelems.row_ptr[i + 1] +
          mesh->stencils[stclidx - 1].ngbelems.row_ptr[i];
    }
  } // single
#pragma omp single
  { // single
    mesh->stencils[stclidx - 1].ngbverts.col_ind.set_size((
        mesh->stencils[stclidx - 1].ngbverts.row_ptr[n] - 1L));
    mesh->stencils[stclidx - 1].ngbelems.col_ind.set_size((
        mesh->stencils[stclidx - 1].ngbelems.row_ptr[n] - 1L));
    mesh->iwork1.set_size(n);
    mesh->iwork2.set_size(n);
  } // single
  //  Local buffers
  vtags_.set_size(mesh->coords.size(0));
  u0 = mesh->coords.size(0);
  for (u1 = 0; u1 < u0; u1++) {
    vtags_[u1] = false;
  }
  ftags_.set_size(mesh->teids.size(0));
  u0 = mesh->teids.size(0);
  for (u1 = 0; u1 < u0; u1++) {
    ftags_[u1] = false;
  }
  hadoverflow = false;
  //  Loop begins
  for (::coder::SizeType i{istart}; i <= iend; i++) {
    int64_T estart;
    int64_T vstart;
    ::coder::SizeType vid;
    vstart = mesh->stencils[stclidx - 1].ngbverts.row_ptr[i - 1] - 1L;
    estart = mesh->stencils[stclidx - 1].ngbelems.row_ptr[i - 1] - 1L;
    //  Prevent Coder from creating copies of rowptrs
    //  Get vertex ID
    if (mesh->stencils[stclidx - 1].vidmap.size(0) == 0) {
      vid = i;
    } else {
      vid = mesh->stencils[stclidx - 1].vidmap[i - 1];
    }
    obtain_nring_1d(mesh->elemtables, mesh->teids, mesh->sibhfs, mesh->v2hfid,
                    vid, ring, bounds[1], vtags_, ftags_, ngbvs_, &nverts,
                    ngbfs_, &nfaces, hebuf_, &reflected, &overflow);
    if ((!hadoverflow) && overflow) {
      m2cWarnMsgIdAndTxt("obtain_nring_2d:overflow",
                         "Buffers are too small to contain neighborhood");
      hadoverflow = true;
    }
    mesh->stencils[stclidx - 1].reflected[i - 1] = reflected;
    mesh->iwork1[i - 1] = nverts;
    for (int64_T j{1L}; j <= nverts; j++) {
      mesh->stencils[stclidx - 1]
          .ngbverts.col_ind[(vstart + j) - 1] =
          ngbvs_[j - 1];
    }
    mesh->iwork2[i - 1] = nfaces;
    for (int64_T j{1L}; j <= nfaces; j++) {
      mesh->stencils[stclidx - 1]
          .ngbelems.col_ind[(estart + j) - 1] =
          ngbfs_[j - 1];
    }
  }
#pragma omp barrier
#pragma omp single nowait
  { // single
    crs_compress(mesh->stencils[stclidx - 1].ngbverts.row_ptr,
                 mesh->stencils[stclidx - 1].ngbverts.col_ind, mesh->iwork1);
  } // single
#pragma omp single
  { // single
    crs_compress(mesh->stencils[stclidx - 1].ngbelems.row_ptr,
                 mesh->stencils[stclidx - 1].ngbelems.col_ind, mesh->iwork2);
  } // single
}

static void b_compute_stencils_kernel_2d(WlsMesh *mesh, ::coder::SizeType stclidx,
                                         real_T ring, const int32_T bounds[2])
{
  ::coder::array<uint64_T, 1U> hebuf_;
  ::coder::array<int32_T, 1U> ngbfs_;
  ::coder::array<int32_T, 1U> ngbvs_;
  ::coder::array<boolean_T, 1U> ftags_;
  ::coder::array<boolean_T, 1U> vtags_;
  ::coder::SizeType iend;
  ::coder::SizeType istart;
  ::coder::SizeType n;
  int32_T nfaces;
  ::coder::SizeType nthreads;
  int32_T nverts;
  ::coder::SizeType u0;
  ::coder::SizeType u1;
  boolean_T hadoverflow;
  boolean_T overflow;
  boolean_T reflected;
  //  Determine total number of vertices
  if (mesh->stencils[stclidx - 1].vidmap.size(0) == 0) {
    n = mesh->coords.size(0);
  } else {
    n = mesh->stencils[stclidx - 1].vidmap.size(0);
  }
  //  Determine partitioning
  nthreads = 1;
#ifdef _OPENMP
  nthreads = omp_get_num_threads();
#endif // _OPENMP
  if (nthreads == 1) {
    istart = 1;
    iend = n;
  } else {
    ::coder::SizeType b_remainder;
    ::coder::SizeType chunk;
    ::coder::SizeType threadID;
    threadID = 0;
#ifdef _OPENMP
    threadID = omp_get_thread_num();
#endif // _OPENMP
    chunk = n / nthreads;
    b_remainder = n - nthreads * chunk;
    u1 = threadID;
    if (b_remainder <= threadID) {
      u1 = b_remainder;
    }
    istart = (threadID * chunk + u1) + 1;
    iend = ((istart + chunk) + (threadID < b_remainder)) - 1;
  }
#pragma omp single
  { // single
    mesh->stencils[stclidx - 1].reflected.set_size(n);
    mesh->stencils[stclidx - 1].ngbverts.row_ptr.set_size(n + 1);
    mesh->stencils[stclidx - 1].ngbverts.row_ptr[0] = 1L;
    mesh->stencils[stclidx - 1].ngbverts.ncols = mesh->coords.size(0);
    mesh->stencils[stclidx - 1].ngbelems.row_ptr.set_size(n + 1);
    mesh->stencils[stclidx - 1].ngbelems.row_ptr[0] = 1L;
    mesh->stencils[stclidx - 1].ngbelems.ncols = mesh->teids.size(0);
  } // single
  //  Assemble rowptrs
  for (::coder::SizeType i{istart}; i <= iend; i++) {
    //  Get vertex ID
    mesh->stencils[stclidx - 1].ngbverts.row_ptr[i] = bounds[1];
    mesh->stencils[stclidx - 1].ngbelems.row_ptr[i] = bounds[1] << 1;
  }
#pragma omp barrier
#pragma omp single nowait
  { // single
    for (::coder::SizeType i{0}; i < n; i++) {
      mesh->stencils[stclidx - 1].ngbverts.row_ptr[i + 1] =
          mesh->stencils[stclidx - 1].ngbverts.row_ptr[i + 1] +
          mesh->stencils[stclidx - 1].ngbverts.row_ptr[i];
    }
  } // single
#pragma omp single
  { // single
    for (::coder::SizeType i{0}; i < n; i++) {
      mesh->stencils[stclidx - 1].ngbelems.row_ptr[i + 1] =
          mesh->stencils[stclidx - 1].ngbelems.row_ptr[i + 1] +
          mesh->stencils[stclidx - 1].ngbelems.row_ptr[i];
    }
  } // single
#pragma omp single
  { // single
    mesh->stencils[stclidx - 1].ngbverts.col_ind.set_size((
        mesh->stencils[stclidx - 1].ngbverts.row_ptr[n] - 1L));
    mesh->stencils[stclidx - 1].ngbelems.col_ind.set_size((
        mesh->stencils[stclidx - 1].ngbelems.row_ptr[n] - 1L));
    mesh->iwork1.set_size(n);
    mesh->iwork2.set_size(n);
  } // single
  //  Local buffers
  vtags_.set_size(mesh->coords.size(0));
  u0 = mesh->coords.size(0);
  for (u1 = 0; u1 < u0; u1++) {
    vtags_[u1] = false;
  }
  ftags_.set_size(mesh->teids.size(0));
  u0 = mesh->teids.size(0);
  for (u1 = 0; u1 < u0; u1++) {
    ftags_[u1] = false;
  }
  hadoverflow = false;
  //  Loop begins
  for (::coder::SizeType i{istart}; i <= iend; i++) {
    int64_T estart;
    int64_T vstart;
    ::coder::SizeType vid;
    vstart = mesh->stencils[stclidx - 1].ngbverts.row_ptr[i - 1] - 1L;
    estart = mesh->stencils[stclidx - 1].ngbelems.row_ptr[i - 1] - 1L;
    //  Prevent Coder from creating copies of rowptrs
    //  Get vertex ID
    if (mesh->stencils[stclidx - 1].vidmap.size(0) == 0) {
      vid = i;
    } else {
      vid = mesh->stencils[stclidx - 1].vidmap[i - 1];
    }
    obtain_nring_2d(mesh->elemtables, mesh->teids, mesh->sibhfs, mesh->v2hfid,
                    vid, ring, bounds[1], vtags_, ftags_, mesh->bridges, ngbvs_,
                    &nverts, ngbfs_, &nfaces, hebuf_, &reflected, &overflow);
    if ((!hadoverflow) && overflow) {
      m2cWarnMsgIdAndTxt("obtain_nring_2d:overflow",
                         "Buffers are too small to contain neighborhood");
      hadoverflow = true;
    }
    mesh->stencils[stclidx - 1].reflected[i - 1] = reflected;
    mesh->iwork1[i - 1] = nverts;
    for (int64_T j{1L}; j <= nverts; j++) {
      mesh->stencils[stclidx - 1]
          .ngbverts.col_ind[(vstart + j) - 1] =
          ngbvs_[j - 1];
    }
    mesh->iwork2[i - 1] = nfaces;
    for (int64_T j{1L}; j <= nfaces; j++) {
      mesh->stencils[stclidx - 1]
          .ngbelems.col_ind[(estart + j) - 1] =
          ngbfs_[j - 1];
    }
  }
#pragma omp barrier
#pragma omp single nowait
  { // single
    crs_compress(mesh->stencils[stclidx - 1].ngbverts.row_ptr,
                 mesh->stencils[stclidx - 1].ngbverts.col_ind, mesh->iwork1);
  } // single
#pragma omp single
  { // single
    crs_compress(mesh->stencils[stclidx - 1].ngbelems.row_ptr,
                 mesh->stencils[stclidx - 1].ngbelems.col_ind, mesh->iwork2);
  } // single
}

static inline
void b_wlsmesh_compute_1ring(WlsMesh *mesh)
{
#pragma omp single
  { // single
    append_wlsmesh_kring(mesh);
  } // single
  //  Next, compute 1-ring
  wlsmesh_compute_stencils(mesh, mesh->stencils.size(0));
#pragma omp barrier
#pragma omp single
  { // single
    mesh->node2nodes = mesh->stencils[mesh->stencils.size(0) - 1].ngbverts;
    mesh->node2elems = mesh->stencils[mesh->stencils.size(0) - 1].ngbelems;
  } // single
}

static inline
void b_wlsmesh_compute_meshprop(WlsMesh *mesh, ::coder::SizeType nrmidx)
{
  if (mesh->teids.size(0) == 0) {
    m2cErrMsgIdAndTxt("wlsmesh_compute_meshprop:missingSetup",
                      "must call setup first");
  }
  if ((mesh->coords.size(1) == mesh->topo_ndims + 1) &&
      ((nrmidx < 1) || (nrmidx > mesh->nrmstables.size(0)))) {
    m2cErrMsgIdAndTxt("wlsmesh_compute_meshprop:missingNrms",
                      "missing normals for surface");
  }
  //  Measure
  compute_measure_kernel(mesh->coords, mesh->elemtables, mesh->teids,
                         mesh->elemmeas);
  //  Sizes
  if (mesh->coords.size(1) == mesh->topo_ndims + 1) {
    //  Surface
    compute_meshsizes_kernel(mesh->coords, mesh->elemtables, mesh->teids,
                             mesh->node2nodes.row_ptr, mesh->node2nodes.col_ind,
                             mesh->nrmstables[nrmidx - 1].normals, mesh->elemh,
                             mesh->nodeh, &mesh->globalh);
  } else {
    compute_meshsizes_kernel(mesh->coords, mesh->elemtables, mesh->teids,
                             mesh->node2nodes.row_ptr, mesh->node2nodes.col_ind,
                             mesh->elemh, mesh->nodeh, &mesh->globalh);
  }
}

static void b_wlsmesh_compute_stencils(WlsMesh *mesh, ::coder::SizeType stclidx,
                                       real_T krings)
{
  real_T mykrings[3];
  if (stclidx > mesh->stencils.size(0)) {
    m2cErrMsgIdAndTxt("wlsmesh_compute_stencils:badIdx",
                      "Invalid kring stencil index %d", (int)stclidx);
  }
  //  Determine ring
  mykrings[1] = 0.0;
  mykrings[2] = 0.0;
  mykrings[0] = krings;
  if (krings < 1.0) {
    mykrings[0] = 1.0;
  }
  //  default 1-ring
  if (mesh->topo_ndims == 1) {
    b_compute_stencils_1d(mesh, stclidx, mykrings);
  } else if (mesh->topo_ndims == 2) {
    b_compute_stencils_2d(mesh, stclidx, mykrings);
  } else {
    m2cAssert(false, "Not impl");
  }
}

//  bar_quadrules - Obtain quadrature points and weights of a bar element.
static void bar_quadrules(::coder::SizeType degree, ::coder::array<real_T, 2U> &cs,
                          ::coder::array<real_T, 1U> &ws)
{
  if (degree <= 1) {
    ::coder::SizeType nqp;
    nqp = ::sfe_qrules::bar_deg1_qrule();
    cs.set_size(nqp, 1);
    ws.set_size(nqp);
    ::sfe_qrules::bar_deg1_qrule(&cs[0], &(ws.data())[0]);
  } else if (degree <= 3) {
    ::coder::SizeType nqp;
    nqp = ::sfe_qrules::bar_deg3_qrule();
    cs.set_size(nqp, 1);
    ws.set_size(nqp);
    ::sfe_qrules::bar_deg3_qrule(&cs[0], &(ws.data())[0]);
  } else if (degree <= 5) {
    ::coder::SizeType nqp;
    nqp = ::sfe_qrules::bar_deg5_qrule();
    cs.set_size(nqp, 1);
    ws.set_size(nqp);
    ::sfe_qrules::bar_deg5_qrule(&cs[0], &(ws.data())[0]);
  } else if (degree <= 7) {
    ::coder::SizeType nqp;
    nqp = ::sfe_qrules::bar_deg7_qrule();
    cs.set_size(nqp, 1);
    ws.set_size(nqp);
    ::sfe_qrules::bar_deg7_qrule(&cs[0], &(ws.data())[0]);
  } else if (degree <= 9) {
    ::coder::SizeType nqp;
    nqp = ::sfe_qrules::bar_deg9_qrule();
    cs.set_size(nqp, 1);
    ws.set_size(nqp);
    ::sfe_qrules::bar_deg9_qrule(&cs[0], &(ws.data())[0]);
  } else if (degree <= 11) {
    ::coder::SizeType nqp;
    nqp = ::sfe_qrules::bar_deg11_qrule();
    cs.set_size(nqp, 1);
    ws.set_size(nqp);
    ::sfe_qrules::bar_deg11_qrule(&cs[0], &(ws.data())[0]);
  } else if (degree <= 13) {
    ::coder::SizeType nqp;
    nqp = ::sfe_qrules::bar_deg13_qrule();
    cs.set_size(nqp, 1);
    ws.set_size(nqp);
    ::sfe_qrules::bar_deg13_qrule(&cs[0], &(ws.data())[0]);
  } else {
    ::coder::SizeType nqp;
    if (degree > 15) {
      m2cWarnMsgIdAndTxt("bar_quadrules:UnsupportedDegree",
                         "Only support up to degree 15");
    }
    nqp = ::sfe_qrules::bar_deg15_qrule();
    cs.set_size(nqp, 1);
    ws.set_size(nqp);
    ::sfe_qrules::bar_deg15_qrule(&cs[0], &(ws.data())[0]);
  }
}

static void build_part(::coder::SizeType nParts,
                       const ::coder::array<int32_T, 1U> &nparts,
                       const ::coder::array<int32_T, 1U> &cparts,
                       const ::coder::array<int32_T, 1U> &eptr,
                       const ::coder::array<int32_T, 1U> &eind,
                       ::coder::array<boolean_T, 1U> &ctags,
                       ::coder::array<int32_T, 1U> &iwork,
                       ::coder::array<int32_T, 1U> &partptr,
                       ::coder::array<int32_T, 1U> &partlist,
                       ::coder::array<int32_T, 1U> &sharedents)
{
  ::coder::SizeType i;
  ::coder::SizeType j;
  ::coder::SizeType k;
  ::coder::SizeType m;
  ::coder::SizeType n0;
  //  local function to compute the partitioning
  m = cparts.size(0) - 1;
  partptr.set_size(nParts + 1);
  for (i = 0; i <= nParts; i++) {
    partptr[i] = 0;
  }
  partptr[0] = 1;
  for (::coder::SizeType b_i{0}; b_i <= m; b_i++) {
    partptr[cparts[b_i]] = partptr[cparts[b_i]] + 1;
  }
  //  accumulate
  for (::coder::SizeType b_i{0}; b_i < nParts; b_i++) {
    partptr[b_i + 1] = partptr[b_i + 1] + partptr[b_i];
  }
  if (iwork.size(0) < partptr[nParts] - 1) {
    iwork.set_size(partptr[nParts] - 1);
  }
  for (::coder::SizeType b_i{0}; b_i <= m; b_i++) {
    i = partptr[cparts[b_i] - 1];
    iwork[i - 1] = b_i + 1;
    partptr[cparts[b_i] - 1] = i + 1;
  }
  //  reset
  for (::coder::SizeType b_i{nParts}; b_i >= 1; b_i--) {
    partptr[b_i] = partptr[b_i - 1];
  }
  partptr[0] = 1;
  //  determine shared region
  for (::coder::SizeType part{0}; part < nParts; part++) {
    ::coder::SizeType i1;
    i = partptr[part];
    i1 = partptr[part + 1] - 1;
    for (::coder::SizeType b_i{i}; b_i <= i1; b_i++) {
      ::coder::SizeType j_tmp;
      boolean_T exitg1;
      boolean_T shared;
      shared = false;
      //  loop through all nodes in cell
      j_tmp = iwork[b_i - 1];
      j = eptr[j_tmp - 1];
      exitg1 = false;
      while ((!exitg1) && (j <= eptr[j_tmp] - 1)) {
        if (nparts[eind[j - 1] - 1] != part + 1) {
          shared = true;
          exitg1 = true;
        } else {
          j++;
        }
      }
      if (shared) {
        ctags[j_tmp - 1] = true;
      }
    }
  }
  //  remove shared cells
  k = 1;
  n0 = 0;
  for (::coder::SizeType b_i{0}; b_i < nParts; b_i++) {
    i = partptr[b_i + 1] - 1;
    for (j = k; j <= i; j++) {
      if (ctags[iwork[j - 1] - 1]) {
        n0++;
      }
    }
    //  record next position
    k = partptr[b_i + 1];
    partptr[b_i + 1] = partptr[b_i + 1] - n0;
  }
  partlist.set_size(partptr[nParts] - 1);
  k = -1;
  for (::coder::SizeType b_i{0}; b_i <= m; b_i++) {
    if (!ctags[iwork[b_i] - 1]) {
      k++;
      partlist[k] = iwork[b_i];
    }
  }
  //  get shared cells
  sharedents.set_size(n0);
  k = -1;
  for (::coder::SizeType b_i{0}; b_i <= m; b_i++) {
    if (ctags[b_i]) {
      k++;
      sharedents[k] = b_i + 1;
      ctags[b_i] = false;
      //  reset
    }
  }
}

static ::coder::SizeType c_append_wlsmesh_kring(WlsMesh *mesh)
{
  static const char_T name[19]{'R', 'd', 'i', 'E', 'x', 't', 'e', 'n', 'd', 'e',
                               'd', 'S', 't', 'e', 'n', 'c', 'i', 'l', 's'};
  ::coder::SizeType stclidx;
  stclidx = mesh->stencils.size(0) + 1;
  // Stencils - Object containing the stencils of vertices in a WlsMesh based on
  mesh->stencils.set_size(mesh->stencils.size(0) + 1);
  mesh->stencils[stclidx - 1].reflected.set_size(0);
  mesh->stencils[stclidx - 1].vidmap.set_size(0);
  mesh->stencils[stclidx - 1].ngbverts.row_ptr.set_size(0);
  mesh->stencils[stclidx - 1].ngbverts.ncols = 0;
  //  Construct A.row_ptr
  mesh->stencils[stclidx - 1].ngbverts.col_ind.set_size(0);
  mesh->stencils[stclidx - 1].ngbelems.row_ptr.set_size(0);
  mesh->stencils[stclidx - 1].ngbelems.ncols = 0;
  //  Construct A.row_ptr
  mesh->stencils[stclidx - 1].ngbelems.col_ind.set_size(0);
  mesh->stencils[stclidx - 1].name.set_size(1, 19);
  for (::coder::SizeType i{0}; i < 19; i++) {
    mesh->stencils[stclidx - 1].name[i] = name[i];
  }
  return stclidx;
}

static inline
void
c_assemble_body(RdiObject *rdi, const ::coder::array<real_T, 2U> &mesh_coords,
                const ::coder::array<ConnData, 1U> &mesh_elemtables,
                const ::coder::array<uint64_T, 1U> &mesh_teids,
                const ::coder::array<int64_T, 1U> &mesh_node2elems_row_ptr,
                const ::coder::array<int32_T, 1U> &mesh_node2elems_col_ind,
                const ::coder::array<Stencils, 1U> &mesh_stencils,
                const ::coder::array<Omp4mPart, 1U> &mesh_nodeparts,
                ::coder::SizeType stclid, boolean_T interp0)
{
  if ((stclid != rdi->stclid) && (stclid != rdi->extstclid)) {
    m2cErrMsgIdAndTxt("assemble_body:badStencilID", "bad stencil ID (index) %d",
                      (int)stclid);
  }
#pragma omp single
  { // single
    rdi->fullrank = true;
  } // single
  if (mesh_stencils[stclid - 1].vidmap.size(0) != 0) {
#pragma omp single
    { // single
      rdi->fullrank = b_assemble_body_range(
          rdi->A.row_ptr, rdi->A.col_ind, mesh_coords, mesh_elemtables,
          mesh_teids, mesh_stencils[stclid - 1].ngbverts.row_ptr,
          mesh_stencils[stclid - 1].ngbverts.col_ind, mesh_node2elems_row_ptr,
          mesh_node2elems_col_ind, rdi->degree, interp0,
          mesh_stencils[stclid - 1].vidmap,
          mesh_stencils[stclid - 1].vidmap.size(0), rdi->A.val, rdi->rdtags);
    } // single
  } else if (mesh_nodeparts.size(0) == 0) {
#pragma omp single
    { // single
      rdi->fullrank = c_assemble_body_range(
          rdi->A.row_ptr, rdi->A.col_ind, mesh_coords, mesh_elemtables,
          mesh_teids, mesh_stencils[stclid - 1].ngbverts.row_ptr,
          mesh_stencils[stclid - 1].ngbverts.col_ind, mesh_node2elems_row_ptr,
          mesh_node2elems_col_ind, rdi->degree, interp0, mesh_coords.size(0),
          rdi->A.val, rdi->rdtags);
    } // single
  } else {
    assemble_body_task(rdi, mesh_coords, mesh_elemtables, mesh_teids,
                       mesh_node2elems_row_ptr, mesh_node2elems_col_ind,
                       mesh_nodeparts,
                       mesh_stencils[stclid - 1].ngbverts.row_ptr,
                       mesh_stencils[stclid - 1].ngbverts.col_ind, interp0);
  }
}

static boolean_T
c_assemble_body_range(const ::coder::array<int64_T, 1U> &row_ptr,
                      const ::coder::array<int32_T, 1U> &col_ind,
                      const ::coder::array<real_T, 2U> &mesh_coords,
                      const ::coder::array<ConnData, 1U> &mesh_elemtables,
                      const ::coder::array<uint64_T, 1U> &mesh_teids,
                      const ::coder::array<int64_T, 1U> &stcl_row_ptr,
                      const ::coder::array<int32_T, 1U> &stcl_col_ind,
                      const ::coder::array<int64_T, 1U> &n2e_row_ptr,
                      const ::coder::array<int32_T, 1U> &n2e_col_ind,
                      ::coder::SizeType degree, boolean_T interp0, ::coder::SizeType iend,
                      ::coder::array<real_T, 1U> &val,
                      ::coder::array<boolean_T, 1U> &rdtags)
{
  static const char_T name[7]{'B', 'u', 'h', 'm', 'a', 'n', 'n'};
  ::coder::array<real_T, 2U> coeffs_;
  ::coder::array<real_T, 2U> w_params_pointwise;
  ::coder::array<real_T, 2U> xs_;
  ::coder::array<real_T, 1U> Ns_;
  ::coder::array<real_T, 1U> w_params_shared;
  ::coder::array<int32_T, 1U> eids_;
  ::coder::array<int32_T, 1U> nodes_;
  ::coder::array<boolean_T, 1U> w_omit_rows;
  b_WlsObject expl_temp;
  e_WlsObject wls_;
  ::coder::SizeType c_loop_ub;
  ::coder::SizeType dim;
  ::coder::SizeType loop_ub;
  char_T wgts__name_data[7];
  boolean_T fullrank;
  fullrank = true;
  f_WlsObject(degree, &expl_temp);
  wls_.runtimes.size[0] = expl_temp.runtimes.size[0];
  loop_ub = expl_temp.runtimes.size[0];
  if (loop_ub - 1 >= 0) {
    std::copy(&expl_temp.runtimes.data[0], &expl_temp.runtimes.data[loop_ub],
              &wls_.runtimes.data[0]);
  }
  wls_.QRt.set_size(expl_temp.QRt.size(0), expl_temp.QRt.size(1));
  loop_ub = expl_temp.QRt.size(1) * expl_temp.QRt.size(0);
  for (::coder::SizeType i{0}; i < loop_ub; i++) {
    wls_.QRt[i] = expl_temp.QRt[i];
  }
  wls_.rowmajor = expl_temp.rowmajor;
  wls_.work.set_size(expl_temp.work.size(0));
  loop_ub = expl_temp.work.size(0);
  for (::coder::SizeType i{0}; i < loop_ub; i++) {
    wls_.work[i] = expl_temp.work[i];
  }
  wls_.jpvt.set_size(expl_temp.jpvt.size(0));
  loop_ub = expl_temp.jpvt.size(0);
  for (::coder::SizeType i{0}; i < loop_ub; i++) {
    wls_.jpvt[i] = expl_temp.jpvt[i];
  }
  wls_.fullrank = expl_temp.fullrank;
  wls_.rank = expl_temp.rank;
  wls_.ncols = expl_temp.ncols;
  wls_.nrows = expl_temp.nrows;
  wls_.nevpnts = expl_temp.nevpnts;
  wls_.rhs.set_size(expl_temp.rhs.size(0), expl_temp.rhs.size(1));
  loop_ub = expl_temp.rhs.size(1) * expl_temp.rhs.size(0);
  for (::coder::SizeType i{0}; i < loop_ub; i++) {
    wls_.rhs[i] = expl_temp.rhs[i];
  }
  wls_.QR.set_size(expl_temp.QR.size(0), expl_temp.QR.size(1));
  loop_ub = expl_temp.QR.size(1) * expl_temp.QR.size(0);
  for (::coder::SizeType i{0}; i < loop_ub; i++) {
    wls_.QR[i] = expl_temp.QR[i];
  }
  wls_.V.set_size(expl_temp.V.size(0), expl_temp.V.size(1));
  loop_ub = expl_temp.V.size(1) * expl_temp.V.size(0);
  for (::coder::SizeType i{0}; i < loop_ub; i++) {
    wls_.V[i] = expl_temp.V[i];
  }
  wls_.hs_inv.set_size(1, expl_temp.hs_inv.size[1]);
  loop_ub = expl_temp.hs_inv.size[1];
  for (::coder::SizeType i{0}; i < loop_ub; i++) {
    wls_.hs_inv[i] = expl_temp.hs_inv.data[i];
  }
  wls_.rweights.set_size(expl_temp.rweights.size(0));
  loop_ub = expl_temp.rweights.size(0);
  for (::coder::SizeType i{0}; i < loop_ub; i++) {
    wls_.rweights[i] = expl_temp.rweights[i];
  }
  wls_.origin.size[1] = expl_temp.origin.size[1];
  wls_.origin.size[0] = 1;
  loop_ub = expl_temp.origin.size[1];
  if (loop_ub - 1 >= 0) {
    std::copy(&expl_temp.origin.data[0], &expl_temp.origin.data[loop_ub],
              &wls_.origin.data[0]);
  }
  wls_.us.set_size(expl_temp.us.size(0), expl_temp.us.size(1));
  loop_ub = expl_temp.us.size(0);
  for (::coder::SizeType i{0}; i < loop_ub; i++) {
    ::coder::SizeType b_loop_ub;
    b_loop_ub = expl_temp.us.size(1);
    for (::coder::SizeType i1{0}; i1 < b_loop_ub; i1++) {
      wls_.us[i1 + wls_.us.size(1) * i] =
          expl_temp.us[i1 + expl_temp.us.size(1) * i];
    }
  }
  wls_.stride = expl_temp.stride;
  wls_.interp0 = expl_temp.interp0;
  wls_.unimono = expl_temp.unimono;
  wls_.order = expl_temp.order;
  wls_.degree = expl_temp.degree;
  wls_.nstpnts = expl_temp.nstpnts;
  w_params_shared.set_size(0);
  w_params_pointwise.set_size(::coder::SizeType(0), ::coder::SizeType(0));
  w_omit_rows.set_size(0);
  for (::coder::SizeType i{0}; i < 7; i++) {
    wgts__name_data[i] = name[i];
  }
  // Local buffers with patterns "(\w+)?params_pointwise" are legitimate
  if (iend >= 1) {
    dim = mesh_coords.size(1);
    c_loop_ub = mesh_coords.size(1);
  }
  for (::coder::SizeType b_i{1}; b_i <= iend; b_i++) {
    int64_T c_i;
    int64_T n;
    int64_T n_tmp;
    ::coder::SizeType k;
    //  Fetch local data
    n_tmp = stcl_row_ptr[b_i - 1];
    n = stcl_row_ptr[b_i] - n_tmp;
    nodes_.set_size((n + 1L));
    xs_.set_size((n + 1L), mesh_coords.size(1));
    nodes_[0] = b_i;
    for (::coder::SizeType i{0}; i < c_loop_ub; i++) {
      xs_[i] = 0.0;
    }
    for (c_i = 2L; c_i - 1L <= n; c_i++) {
      k = stcl_col_ind[static_cast<::coder::SizeType>((n_tmp + c_i) - 2L) - 1];
      nodes_[c_i - 1] = k;
      for (::coder::SizeType j{0}; j < dim; j++) {
        xs_[j + xs_.size(1) * (c_i - 1)] =
            mesh_coords[j + mesh_coords.size(1) * (k - 1)] -
            mesh_coords[j + mesh_coords.size(1) * (b_i - 1)];
      }
    }
    //  Compute wls
    wls_init(&wls_, xs_, wgts__name_data, w_params_shared, w_params_pointwise,
             w_omit_rows, degree, interp0, xs_.size(0));
    if (wls_.rank < 0) {
      //  LAPACK error
      m2cErrMsgIdAndTxt("assemble_body_range:badLapack",
                        "LAPACK error code %d for node %d", wls_.rank, (int)b_i);
    }
    if (!wls_.fullrank) {
      //  Not full rank
      rdtags[b_i - 1] = true;
      fullrank = false;
    } else {
      int64_T m;
      ::coder::SizeType b_dim;
      ::coder::SizeType b_m;
      ::coder::SizeType b_n;
      //  Get centers and shape function Ns at centers
      b_dim = mesh_coords.size(1) - 1;
      //  # of nearby elements
      n_tmp = n2e_row_ptr[b_i - 1];
      m = n2e_row_ptr[b_i] - n_tmp;
      eids_.set_size(m);
      xs_.set_size(m, mesh_coords.size(1));
      Ns_.set_size(m);
      for (int64_T b_j{1L}; b_j <= m; b_j++) {
        uint64_T c;
        ::coder::SizeType eid;
        ::coder::SizeType leid;
        ::coder::SizeType npe;
        eid = n2e_col_ind[static_cast<::coder::SizeType>((n_tmp + b_j) - 1L) - 1] - 1;
        eids_[b_j - 1] = eid + 1;
        c = mesh_teids[eid] & 255UL;
        leid = (mesh_teids[eid] >> 8) - 1;
        npe = mesh_elemtables
                  [static_cast<::coder::SizeType>(
                       mesh_teids[n2e_col_ind[static_cast<::coder::SizeType>(
                                                  (n2e_row_ptr[b_i - 1] + b_j) -
                                                  1L) -
                                              1] -
                                  1] &
                       255UL) -
                   1]
                      .conn.size(1);
        Ns_[b_j - 1] =
            1.0 /
            static_cast<real_T>(
                mesh_elemtables[c - 1].conn.size(1));
        //  Compute localized center
        for (k = 0; k <= b_dim; k++) {
          xs_[k + xs_.size(1) * (b_j - 1)] = mesh_coords
              [k + mesh_coords.size(1) *
                       (mesh_elemtables[c - 1]
                            .conn[mesh_elemtables[c - 1]
                                      .conn.size(1) *
                                  leid] -
                        1)];
        }
        for (::coder::SizeType d_i{2}; d_i <= npe; d_i++) {
          for (k = 0; k <= b_dim; k++) {
            xs_[k + xs_.size(1) * (b_j - 1)] =
                xs_[k + xs_.size(1) * (b_j - 1)] +
                mesh_coords
                    [k + mesh_coords.size(1) *
                             (mesh_elemtables[c - 1].conn
                                  [(d_i +
                                    mesh_elemtables[c - 1]
                                            .conn.size(1) *
                                        leid) -
                                   1] -
                              1)];
          }
        }
        for (k = 0; k <= b_dim; k++) {
          xs_[k + xs_.size(1) * (b_j - 1)] =
              xs_[k + xs_.size(1) * (b_j - 1)] *
                  Ns_[b_j - 1] -
              mesh_coords[k + mesh_coords.size(1) * (b_i - 1)];
        }
      }
      //  Compute WLS fitting at centers
      wls_func(&wls_, xs_, coeffs_);
      //  Compute local OSUS coefficients
      b_m = Ns_.size(0) - 1;
      b_n = coeffs_.size(0);
      //  WALF fitting
      for (k = 0; k < b_n; k++) {
        for (::coder::SizeType j{0}; j <= b_m; j++) {
          coeffs_[j + coeffs_.size(1) * k] =
              coeffs_[j + coeffs_.size(1) * k] * Ns_[j];
        }
      }
      //  Substract linear interp
      for (::coder::SizeType j{0}; j <= b_m; j++) {
        coeffs_[j] = coeffs_[j] - Ns_[j];
      }
      //  Add to global CRS
      b_n = eids_.size(0);
      b_m = nodes_.size(0);
      for (::coder::SizeType j{0}; j < b_n; j++) {
        ::coder::SizeType r;
        r = eids_[j];
        for (::coder::SizeType d_i{0}; d_i < b_m; d_i++) {
          int64_T b_k;
          b_k = 0L;
          if (r < row_ptr.size(0)) {
            int64_T bounds_idx_1;
            //  Perform linear search
            bounds_idx_1 = row_ptr[r] - 1L;
            //  Perform linear search
            c_i = row_ptr[r - 1];
            while (c_i <= bounds_idx_1) {
              if (col_ind[c_i - 1] == nodes_[d_i]) {
                b_k = c_i;
                c_i = bounds_idx_1 + 1L;
              } else {
                c_i++;
              }
            }
          }
          val[b_k - 1] =
              val[b_k - 1] +
              coeffs_[j + coeffs_.size(1) * d_i];
        }
      }
    }
  }
  return fullrank;
}

static void
c_assemble_surf(RdiObject *rdi, const ::coder::array<real_T, 2U> &mesh_coords,
                const ::coder::array<ConnData, 1U> &mesh_elemtables,
                const ::coder::array<uint64_T, 1U> &mesh_teids,
                const ::coder::array<NormalsData, 1U> &mesh_nrmstables,
                const ::coder::array<int64_T, 1U> &mesh_node2elems_row_ptr,
                const ::coder::array<int32_T, 1U> &mesh_node2elems_col_ind,
                const ::coder::array<Stencils, 1U> &mesh_stencils,
                const ::coder::array<Omp4mPart, 1U> &mesh_nodeparts,
                ::coder::SizeType stclid, boolean_T interp0)
{
  ::coder::SizeType nrmid;
  if ((stclid != rdi->stclid) && (stclid != rdi->extstclid)) {
    m2cErrMsgIdAndTxt("assemble_surf:badStencilID", "bad stencil ID (index) %d",
                      (int)stclid);
  }
#pragma omp single
  { // single
    rdi->fullrank = true;
  } // single
  //  Overcome a "bug" in Coder regarding creating copies of normal data
  nrmid = rdi->nrmid;
  if (mesh_stencils[stclid - 1].vidmap.size(0) != 0) {
#pragma omp single
    { // single
      rdi->fullrank = b_assemble_surf_range(
          rdi->A.row_ptr, rdi->A.col_ind, mesh_coords, mesh_elemtables,
          mesh_teids, mesh_stencils[stclid - 1].ngbverts.row_ptr,
          mesh_stencils[stclid - 1].ngbverts.col_ind, mesh_node2elems_row_ptr,
          mesh_node2elems_col_ind, mesh_nrmstables[nrmid - 1].normals,
          rdi->degree, interp0, mesh_stencils[stclid - 1].vidmap,
          mesh_stencils[stclid - 1].vidmap.size(0), rdi->A.val, rdi->rdtags);
    } // single
  } else if (mesh_nodeparts.size(0) == 0) {
#pragma omp single
    { // single
      rdi->fullrank = c_assemble_surf_range(
          rdi->A.row_ptr, rdi->A.col_ind, mesh_coords, mesh_elemtables,
          mesh_teids, mesh_stencils[stclid - 1].ngbverts.row_ptr,
          mesh_stencils[stclid - 1].ngbverts.col_ind, mesh_node2elems_row_ptr,
          mesh_node2elems_col_ind, mesh_nrmstables[nrmid - 1].normals,
          rdi->degree, interp0, mesh_coords.size(0), rdi->A.val, rdi->rdtags);
    } // single
  } else {
    assemble_surf_task(rdi, mesh_coords, mesh_elemtables, mesh_teids,
                       mesh_node2elems_row_ptr, mesh_node2elems_col_ind,
                       mesh_nodeparts,
                       mesh_stencils[stclid - 1].ngbverts.row_ptr,
                       mesh_stencils[stclid - 1].ngbverts.col_ind,
                       mesh_nrmstables[nrmid - 1].normals, interp0);
  }
}

static boolean_T
c_assemble_surf_range(const ::coder::array<int64_T, 1U> &row_ptr,
                      const ::coder::array<int32_T, 1U> &col_ind,
                      const ::coder::array<real_T, 2U> &mesh_coords,
                      const ::coder::array<ConnData, 1U> &mesh_elemtables,
                      const ::coder::array<uint64_T, 1U> &mesh_teids,
                      const ::coder::array<int64_T, 1U> &stcl_row_ptr,
                      const ::coder::array<int32_T, 1U> &stcl_col_ind,
                      const ::coder::array<int64_T, 1U> &n2e_row_ptr,
                      const ::coder::array<int32_T, 1U> &n2e_col_ind,
                      const ::coder::array<real_T, 2U> &nrms, ::coder::SizeType degree,
                      boolean_T interp0, ::coder::SizeType iend,
                      ::coder::array<real_T, 1U> &val,
                      ::coder::array<boolean_T, 1U> &rdtags)
{
  static const char_T name[7]{'B', 'u', 'h', 'm', 'a', 'n', 'n'};
  ::coder::array<real_T, 2U> coeffs_;
  ::coder::array<real_T, 2U> us_;
  ::coder::array<real_T, 2U> wgts__params_pointwise;
  ::coder::array<real_T, 2U> xs_;
  ::coder::array<real_T, 1U> Ns_;
  ::coder::array<real_T, 1U> w_params_shared;
  ::coder::array<int32_T, 1U> eids_;
  ::coder::array<int32_T, 1U> nodes_;
  ::coder::array<boolean_T, 1U> w_omit_rows;
  b_WlsObject r;
  c_WlsObject expl_temp;
  d_WlsObject b_expl_temp;
  e_WlsObject wls_;
  real_T t_data[6];
  ::coder::SizeType b_loop_ub;
  ::coder::SizeType dim;
  ::coder::SizeType i;
  ::coder::SizeType i1;
  ::coder::SizeType loop_ub;
  ::coder::SizeType t_size_idx_1;
  char_T wgts__name_data[7];
  boolean_T fullrank;
  fullrank = true;
  f_WlsObject(degree, &r);
  majorityTransform(&r, &expl_temp);
  b_expl_temp.us.set_size(expl_temp.us.size(0), expl_temp.us.size(1));
  loop_ub = expl_temp.us.size(1) * expl_temp.us.size(0);
  for (i = 0; i < loop_ub; i++) {
    b_expl_temp.us[i] = expl_temp.us[i];
  }
  b_expl_temp.hs_inv.set_size(expl_temp.hs_inv.size[0], 1);
  loop_ub = expl_temp.hs_inv.size[0];
  for (i = 0; i < loop_ub; i++) {
    b_expl_temp.hs_inv[i] = expl_temp.hs_inv.data[i];
  }
  b_expl_temp.runtimes.size[0] = expl_temp.runtimes.size[0];
  loop_ub = expl_temp.runtimes.size[0];
  if (loop_ub - 1 >= 0) {
    std::copy(&expl_temp.runtimes.data[0], &expl_temp.runtimes.data[loop_ub],
              &b_expl_temp.runtimes.data[0]);
  }
  b_expl_temp.QRt.set_size(expl_temp.QRt.size(0), expl_temp.QRt.size(1));
  loop_ub = expl_temp.QRt.size(1) * expl_temp.QRt.size(0);
  for (i = 0; i < loop_ub; i++) {
    b_expl_temp.QRt[i] = expl_temp.QRt[i];
  }
  b_expl_temp.rowmajor = expl_temp.rowmajor;
  b_expl_temp.work.set_size(expl_temp.work.size(0));
  loop_ub = expl_temp.work.size(0);
  for (i = 0; i < loop_ub; i++) {
    b_expl_temp.work[i] = expl_temp.work[i];
  }
  b_expl_temp.jpvt.set_size(expl_temp.jpvt.size(0));
  loop_ub = expl_temp.jpvt.size(0);
  for (i = 0; i < loop_ub; i++) {
    b_expl_temp.jpvt[i] = expl_temp.jpvt[i];
  }
  b_expl_temp.fullrank = expl_temp.fullrank;
  b_expl_temp.rank = expl_temp.rank;
  b_expl_temp.ncols = expl_temp.ncols;
  b_expl_temp.nrows = expl_temp.nrows;
  b_expl_temp.nevpnts = expl_temp.nevpnts;
  b_expl_temp.rhs.set_size(expl_temp.rhs.size(0), expl_temp.rhs.size(1));
  loop_ub = expl_temp.rhs.size(1) * expl_temp.rhs.size(0);
  for (i = 0; i < loop_ub; i++) {
    b_expl_temp.rhs[i] = expl_temp.rhs[i];
  }
  b_expl_temp.QR.set_size(expl_temp.QR.size(0), expl_temp.QR.size(1));
  loop_ub = expl_temp.QR.size(1) * expl_temp.QR.size(0);
  for (i = 0; i < loop_ub; i++) {
    b_expl_temp.QR[i] = expl_temp.QR[i];
  }
  b_expl_temp.V.set_size(expl_temp.V.size(0), expl_temp.V.size(1));
  loop_ub = expl_temp.V.size(1) * expl_temp.V.size(0);
  for (i = 0; i < loop_ub; i++) {
    b_expl_temp.V[i] = expl_temp.V[i];
  }
  b_expl_temp.rweights.set_size(expl_temp.rweights.size(0));
  loop_ub = expl_temp.rweights.size(0);
  for (i = 0; i < loop_ub; i++) {
    b_expl_temp.rweights[i] = expl_temp.rweights[i];
  }
  b_expl_temp.origin.size[1] = 1;
  b_expl_temp.origin.size[0] = expl_temp.origin.size[0];
  loop_ub = expl_temp.origin.size[0];
  if (loop_ub - 1 >= 0) {
    std::copy(&expl_temp.origin.data[0], &expl_temp.origin.data[loop_ub],
              &b_expl_temp.origin.data[0]);
  }
  b_expl_temp.stride = expl_temp.stride;
  b_expl_temp.interp0 = expl_temp.interp0;
  b_expl_temp.unimono = expl_temp.unimono;
  b_expl_temp.order = expl_temp.order;
  b_expl_temp.degree = expl_temp.degree;
  b_expl_temp.nstpnts = expl_temp.nstpnts;
  majorityTransform(&b_expl_temp, &wls_);
  w_params_shared.set_size(0);
  w_omit_rows.set_size(0);
  for (i = 0; i < 7; i++) {
    wgts__name_data[i] = name[i];
  }
  // Local buffers with patterns "wgts__params_pointwise" are legitimate
  if (iend >= 1) {
    dim = mesh_coords.size(1);
    b_loop_ub = mesh_coords.size(1);
    t_size_idx_1 = nrms.size(1) - 1;
    i1 = nrms.size(1);
  }
  for (::coder::SizeType b_i{1}; b_i <= iend; b_i++) {
    int64_T n;
    int64_T n_tmp;
    ::coder::SizeType b_dim;
    ::coder::SizeType b_npoints;
    ::coder::SizeType k;
    ::coder::SizeType npoints;
    boolean_T b;
    boolean_T b1;
    //  Fetch local data
    n_tmp = stcl_row_ptr[b_i - 1];
    n = stcl_row_ptr[b_i] - n_tmp;
    nodes_.set_size((n + 1L));
    xs_.set_size((n + 1L), mesh_coords.size(1));
    nodes_[0] = b_i;
    for (i = 0; i < b_loop_ub; i++) {
      xs_[i] = 0.0;
    }
    for (int64_T c_i{2L}; c_i - 1L <= n; c_i++) {
      k = stcl_col_ind[static_cast<::coder::SizeType>((n_tmp + c_i) - 2L) - 1];
      nodes_[c_i - 1] = k;
      for (::coder::SizeType j{0}; j < dim; j++) {
        xs_[j + xs_.size(1) * (c_i - 1)] =
            mesh_coords[j + mesh_coords.size(1) * (k - 1)] -
            mesh_coords[j + mesh_coords.size(1) * (b_i - 1)];
      }
    }
    npoints = xs_.size(0);
    //  Get normal direction
    if (i1 == 2) {
      t_data[0] = -nrms[nrms.size(1) * (b_i - 1) + 1];
      t_data[t_size_idx_1] = nrms[nrms.size(1) * (b_i - 1)];
    } else {
      real_T a;
      real_T a_tmp;
      real_T d;
      real_T d1;
      boolean_T guard1{false};
      d = nrms[nrms.size(1) * (b_i - 1)];
      d1 = nrms[nrms.size(1) * (b_i - 1) + 1];
      a_tmp = std::abs(d);
      guard1 = false;
      if (a_tmp > std::abs(d1)) {
        real_T d2;
        d2 = nrms[nrms.size(1) * (b_i - 1) + 2];
        if (a_tmp > std::abs(d2)) {
          t_data[0] = -d * d1;
          t_data[t_size_idx_1] = 1.0 - d1 * d1;
          t_data[t_size_idx_1 * 2] = -d1 * d2;
        } else {
          guard1 = true;
        }
      } else {
        guard1 = true;
      }
      if (guard1) {
        t_data[0] = 1.0 - d * d;
        t_data[t_size_idx_1] = -d * d1;
        t_data[t_size_idx_1 * 2] = -d * nrms[nrms.size(1) * (b_i - 1) + 2];
      }
      a_tmp = t_data[t_size_idx_1 * 2];
      a = std::sqrt((t_data[0] * t_data[0] +
                     t_data[t_size_idx_1] * t_data[t_size_idx_1]) +
                    a_tmp * a_tmp);
      t_data[0] /= a;
      t_data[t_size_idx_1] /= a;
      t_data[t_size_idx_1 * 2] /= a;
      //  cross
      a_tmp = nrms[nrms.size(1) * (b_i - 1) + 2];
      t_data[1] = t_data[t_size_idx_1 * 2] * d1 - t_data[t_size_idx_1] * a_tmp;
      t_data[t_size_idx_1 + 1] =
          t_data[0] * a_tmp - d * t_data[t_size_idx_1 * 2];
      t_data[t_size_idx_1 * 2 + 1] = d * t_data[t_size_idx_1] - t_data[0] * d1;
    }
    //  Project onto tangent plane
    b_npoints = xs_.size(0) - 1;
    b_dim = xs_.size(1) - 2;
    us_.set_size(xs_.size(0), xs_.size(1) - 1);
    for (k = 0; k <= b_npoints; k++) {
      for (::coder::SizeType j{0}; j <= b_dim; j++) {
        us_[j + us_.size(1) * k] = 0.0;
      }
    }
    //  matrix-matrix, using mem-efficient loop
    for (::coder::SizeType ii{0}; ii <= b_npoints; ii++) {
      for (k = 0; k <= b_dim + 1; k++) {
        for (::coder::SizeType j{0}; j <= b_dim; j++) {
          us_[j + us_.size(1) * ii] =
              us_[j + us_.size(1) * ii] +
              xs_[k + xs_.size(1) * ii] * t_data[j + t_size_idx_1 * k];
        }
      }
    }
    //  Compute normal matches as additional weights to WLS
    wgts__params_pointwise.set_size(xs_.size(0), 1);
    wgts__params_pointwise[0] = 1.0;
    b = true;
    i = wgts__params_pointwise.size(0);
    b1 = wgts__params_pointwise.size(0) <= 0;
    loop_ub = 0;
    for (::coder::SizeType j{2}; j <= npoints; j++) {
      real_T sigma;
      ::coder::SizeType sigma_tmp;
      if (b1 || (j - 1 >= i)) {
        loop_ub = 0;
        b = true;
      } else if (b) {
        b = false;
        loop_ub = (j - 1) % wgts__params_pointwise.size(0) +
                  (j - 1) / wgts__params_pointwise.size(0);
      } else if (loop_ub > 2147483646) {
        loop_ub = (j - 1) % wgts__params_pointwise.size(0) +
                  (j - 1) / wgts__params_pointwise.size(0);
      } else {
        loop_ub++;
        if (loop_ub > wgts__params_pointwise.size(0) - 1) {
          loop_ub = (loop_ub - wgts__params_pointwise.size(0)) + 1;
        }
      }
      sigma_tmp = nodes_[j - 1] - 1;
      sigma = nrms[nrms.size(1) * (b_i - 1)] * nrms[nrms.size(1) * sigma_tmp] +
              nrms[nrms.size(1) * (b_i - 1) + 1] *
                  nrms[nrms.size(1) * sigma_tmp + 1];
      if (nrms.size(1) > 2) {
        sigma += nrms[nrms.size(1) * (b_i - 1) + 2] *
                 nrms[nrms.size(1) * sigma_tmp + 2];
      }
      wgts__params_pointwise[loop_ub] = sigma;
    }
    //  Compute wls
    wls_init(&wls_, us_, wgts__name_data, w_params_shared,
             wgts__params_pointwise, w_omit_rows, degree,
             interp0, xs_.size(0));
    if (wls_.rank < 0) {
      //  LAPACK error
      m2cErrMsgIdAndTxt("assemble_surf_range:badLapack",
                        "LAPACK error code %d for node %d", wls_.rank, (int)b_i);
    }
    if (!wls_.fullrank) {
      //  Not full rank
      rdtags[b_i - 1] = true;
      fullrank = false;
    } else {
      int64_T m;
      ::coder::SizeType b_m;
      ::coder::SizeType b_n;
      //  Get centers and shape function Ns at centers
      b_dim = mesh_coords.size(1) - 1;
      //  # of nearby elements
      n_tmp = n2e_row_ptr[b_i - 1];
      m = n2e_row_ptr[b_i] - n_tmp;
      eids_.set_size(m);
      xs_.set_size(m, mesh_coords.size(1));
      Ns_.set_size(m);
      for (int64_T b_j{1L}; b_j <= m; b_j++) {
        uint64_T c;
        ::coder::SizeType eid;
        ::coder::SizeType leid;
        ::coder::SizeType npe;
        eid = n2e_col_ind[static_cast<::coder::SizeType>((n_tmp + b_j) - 1L) - 1] - 1;
        eids_[b_j - 1] = eid + 1;
        c = mesh_teids[eid] & 255UL;
        leid = (mesh_teids[eid] >> 8) - 1;
        npe = mesh_elemtables
                  [static_cast<::coder::SizeType>(
                       mesh_teids[n2e_col_ind[static_cast<::coder::SizeType>(
                                                  (n2e_row_ptr[b_i - 1] + b_j) -
                                                  1L) -
                                              1] -
                                  1] &
                       255UL) -
                   1]
                      .conn.size(1);
        Ns_[b_j - 1] =
            1.0 /
            static_cast<real_T>(
                mesh_elemtables[c - 1].conn.size(1));
        //  Compute localized center
        for (k = 0; k <= b_dim; k++) {
          xs_[k + xs_.size(1) * (b_j - 1)] = mesh_coords
              [k + mesh_coords.size(1) *
                       (mesh_elemtables[c - 1]
                            .conn[mesh_elemtables[c - 1]
                                      .conn.size(1) *
                                  leid] -
                        1)];
        }
        for (::coder::SizeType d_i{2}; d_i <= npe; d_i++) {
          for (k = 0; k <= b_dim; k++) {
            xs_[k + xs_.size(1) * (b_j - 1)] =
                xs_[k + xs_.size(1) * (b_j - 1)] +
                mesh_coords
                    [k + mesh_coords.size(1) *
                             (mesh_elemtables[c - 1].conn
                                  [(d_i +
                                    mesh_elemtables[c - 1]
                                            .conn.size(1) *
                                        leid) -
                                   1] -
                              1)];
          }
        }
        for (k = 0; k <= b_dim; k++) {
          xs_[k + xs_.size(1) * (b_j - 1)] =
              xs_[k + xs_.size(1) * (b_j - 1)] *
                  Ns_[b_j - 1] -
              mesh_coords[k + mesh_coords.size(1) * (b_i - 1)];
        }
      }
      //  Project centers onto tangent plane
      npoints = xs_.size(0) - 1;
      b_dim = xs_.size(1) - 2;
      us_.set_size(xs_.size(0), xs_.size(1) - 1);
      for (k = 0; k <= npoints; k++) {
        for (::coder::SizeType j{0}; j <= b_dim; j++) {
          us_[j + us_.size(1) * k] = 0.0;
        }
      }
      //  matrix-matrix, using mem-efficient loop
      for (::coder::SizeType ii{0}; ii <= npoints; ii++) {
        for (k = 0; k <= b_dim + 1; k++) {
          for (::coder::SizeType j{0}; j <= b_dim; j++) {
            us_[j + us_.size(1) * ii] =
                us_[j + us_.size(1) * ii] +
                xs_[k + xs_.size(1) * ii] * t_data[j + t_size_idx_1 * k];
          }
        }
      }
      //  Compute WLS fitting at centers
      wls_func(&wls_, us_, coeffs_);
      //  Compute local OSUS coefficients
      b_m = Ns_.size(0) - 1;
      b_n = coeffs_.size(0);
      //  WALF fitting
      for (k = 0; k < b_n; k++) {
        for (::coder::SizeType j{0}; j <= b_m; j++) {
          coeffs_[j + coeffs_.size(1) * k] =
              coeffs_[j + coeffs_.size(1) * k] * Ns_[j];
        }
      }
      //  Substract linear interp
      for (::coder::SizeType j{0}; j <= b_m; j++) {
        coeffs_[j] = coeffs_[j] - Ns_[j];
      }
      //  Add to global CRS
      b_n = eids_.size(0);
      b_m = nodes_.size(0);
      for (::coder::SizeType j{0}; j < b_n; j++) {
        ::coder::SizeType b_r;
        b_r = eids_[j];
        for (::coder::SizeType d_i{0}; d_i < b_m; d_i++) {
          int64_T b_k;
          b_k = 0L;
          if (b_r < row_ptr.size(0)) {
            //  Perform linear search
            b_k = m2cFind(col_ind, nodes_[d_i], row_ptr[b_r - 1],
                          row_ptr[b_r] - 1L);
          }
          val[b_k - 1] =
              val[b_k - 1] +
              coeffs_[j + coeffs_.size(1) * d_i];
        }
      }
    }
  }
  return fullrank;
}

static void call_metis_mesh(int32_T n, ::coder::array<int32_T, 1U> &eptr,
                            ::coder::array<int32_T, 1U> &eind, int32_T nParts,
                            ::coder::array<int32_T, 1U> &nparts,
                            ::coder::array<int32_T, 1U> &cparts)
{
  int32_T opts_data[40];
  ::coder::SizeType loop_ub;
  int32_T m;
  ::coder::SizeType status;
  //  this function calls METIS
  m = eptr.size(0) - 1;
  //  num. of cells
  nparts.set_size(n);
  for (::coder::SizeType i{0}; i < n; i++) {
    nparts[i] = 1;
  }
  cparts.set_size(eptr.size(0) - 1);
  loop_ub = eptr.size(0);
  for (::coder::SizeType i{0}; i <= loop_ub - 2; i++) {
    cparts[i] = 1;
  }
  int32_t ncuts;
#ifdef OMP4M_HAS_METIS
  //  if we have METIS
  m2cAssert((int32_T)sizeof(idx_t) == 4,
            "Assuming unsigned 32-bit integer-type for mtmetis_vtx_type.");
  //  setup options
  m2cAssert((int32_T)METIS_NOPTIONS <= 40, "");
  METIS_SetDefaultOptions(&opts_data[0]);
  opts_data[(int32_T)METIS_OPTION_NUMBERING] = 1;
  //  1-based
  status = METIS_PartMeshNodal(&m, &n, &(eptr.data())[0], &(eind.data())[0],
                               NULL, NULL, &nParts, NULL, &opts_data[0], &ncuts,
                               &(cparts.data())[0], &(nparts.data())[0]);
  if (status != (int32_T)METIS_OK) {
    m2cErrMsgIdAndTxt("metis:MetisError", "METIS returned error %d.", (int)status);
  }
#else
  //  if we do NOT have METIS
  m2cWarnMsgIdAndTxt(
      "metis:missingMetis",
      "METIS is not available. All nodes are assigned to partition 1.");
#endif
}

namespace coder {
static real_T sum(const ::coder::array<real_T, 1U> &x)
{
  real_T y;
  ::coder::SizeType vlen;
  vlen = x.size(0);
  if (x.size(0) == 0) {
    y = 0.0;
  } else {
    ::coder::SizeType firstBlockLength;
    ::coder::SizeType lastBlockLength;
    ::coder::SizeType nblocks;
    if (x.size(0) <= 1024) {
      firstBlockLength = x.size(0);
      lastBlockLength = 0;
      nblocks = 1;
    } else {
      firstBlockLength = 1024;
      nblocks = x.size(0) / 1024;
      lastBlockLength = x.size(0) - (nblocks << 10);
      if (lastBlockLength > 0) {
        nblocks++;
      } else {
        lastBlockLength = 1024;
      }
    }
    y = x[0];
    for (::coder::SizeType k{2}; k <= firstBlockLength; k++) {
      if (vlen >= 2) {
        y += x[k - 1];
      }
    }
    for (::coder::SizeType ib{2}; ib <= nblocks; ib++) {
      real_T bsum;
      ::coder::SizeType hi;
      firstBlockLength = (ib - 1) << 10;
      bsum = x[firstBlockLength];
      if (ib == nblocks) {
        hi = lastBlockLength;
      } else {
        hi = 1024;
      }
      for (::coder::SizeType k{2}; k <= hi; k++) {
        if (vlen >= 2) {
          bsum += x[(firstBlockLength + k) - 1];
        }
      }
      y += bsum;
    }
  }
  return y;
}

} // namespace coder
static void
compute_beta_kernel(const ::coder::array<real_T, 2U> &mesh_coords,
                    const ::coder::array<int64_T, 1U> &mesh_node2elems_row_ptr,
                    const ::coder::array<int32_T, 1U> &mesh_node2elems_col_ind,
                    const ::coder::array<real_T, 1U> &mesh_elemmeas,
                    real_T mesh_globalh, const ::coder::array<real_T, 2U> &df,
                    const ::coder::array<real_T, 2U> &alpha, real_T epsbeta,
                    ::coder::array<real_T, 2U> &beta)
{
  real_T epsh2;
  ::coder::SizeType iend;
  ::coder::SizeType istart;
  ::coder::SizeType nrhs;
  ::coder::SizeType nthreads;
  ::coder::SizeType u1;
  nrhs = df.size(1);
#pragma omp single
  { // single
    beta.set_size(mesh_coords.size(0), df.size(1));
  } // single
  nthreads = 1;
#ifdef _OPENMP
  nthreads = omp_get_num_threads();
#endif // _OPENMP
  if (nthreads == 1) {
    istart = 0;
    iend = mesh_coords.size(0);
  } else {
    ::coder::SizeType b_remainder;
    ::coder::SizeType chunk;
    ::coder::SizeType threadID;
    threadID = 0;
#ifdef _OPENMP
    threadID = omp_get_thread_num();
#endif // _OPENMP
    chunk = mesh_coords.size(0) / nthreads;
    b_remainder = mesh_coords.size(0) - nthreads * chunk;
    u1 = threadID;
    if (b_remainder <= threadID) {
      u1 = b_remainder;
    }
    istart = threadID * chunk + u1;
    iend = (istart + chunk) + (threadID < b_remainder);
  }
  epsh2 = epsbeta * mesh_globalh * mesh_globalh;
  //  C++
  for (::coder::SizeType i{istart + 1}; i <= iend; i++) {
    for (::coder::SizeType k{0}; k < nrhs; k++) {
      real_T abar;
      real_T asum;
      real_T atop;
      real_T wbar;
      asum = 0.0;
      wbar = 0.0;
      for (int64_T j = mesh_node2elems_row_ptr[i - 1];
           j < mesh_node2elems_row_ptr[i]; j++) {
        real_T asum_tmp;
        u1 = mesh_node2elems_col_ind[j - 1] - 1;
        asum_tmp = mesh_elemmeas[u1];
        asum += asum_tmp * alpha[k + alpha.size(1) * u1];
        wbar += asum_tmp;
      }
      //  sum(w_e*a_e)/sum(w_e)
      abar = asum / wbar;
      atop = 0.0;
      for (int64_T j = mesh_node2elems_row_ptr[i - 1];
           j < mesh_node2elems_row_ptr[i]; j++) {
        u1 = mesh_node2elems_col_ind[j - 1] - 1;
        atop +=
            mesh_elemmeas[u1] * std::abs(alpha[k + alpha.size(1) * u1] - abar);
      }
      beta[k + beta.size(1) * (i - 1)] =
          atop /
          ((std::abs(asum) + wbar * (df[k] * epsh2)) + 2.2250738585072014E-308);
    }
  }
}

static ::coder::SizeType compute_connected_components(
    ::coder::array<boolean_T, 1U> &visited, ::coder::array<int32_T, 1U> &iwork,
    const ::coder::array<int64_T, 1U> &G_row_ptr,
    const ::coder::array<int32_T, 1U> &G_col_ind, ::coder::SizeType G_ncols)
{
  ::coder::SizeType nc;
  //  Kernel for computing #. of connected components
  nc = 0;
  for (::coder::SizeType i{0}; i < G_ncols; i++) {
    if (!visited[i]) {
      ::coder::SizeType stkptr;
      stkptr = 1;
      iwork[0] = i + 1;
      while (stkptr != 0) {
        ::coder::SizeType v;
        v = iwork[stkptr - 1] - 1;
        iwork[stkptr - 1] = 0;
        stkptr--;
        if (!visited[v]) {
          visited[v] = true;
          for (int64_T b_i = G_row_ptr[v]; b_i < G_row_ptr[v + 1]; b_i++) {
            ::coder::SizeType c_i;
            c_i = G_col_ind[b_i - 1];
            if (!visited[c_i - 1]) {
              iwork[stkptr] = c_i;
              stkptr++;
            }
          }
        }
      }
      nc++;
    }
  }
  //  Reset
  for (::coder::SizeType i{0}; i < G_ncols; i++) {
    visited[i] = false;
  }
  return nc;
}

static void
compute_fconn_graph(::coder::array<boolean_T, 1U> &visited,
                    ::coder::array<int32_T, 1U> &iwork,
                    const ::coder::array<ConnData, 1U> &mesh_elemtables,
                    const ::coder::array<uint64_T, 1U> &mesh_teids,
                    const ::coder::array<int64_T, 1U> &mesh_node2elems_row_ptr,
                    const ::coder::array<int32_T, 1U> &mesh_node2elems_col_ind,
                    ::coder::SizeType nid, const ::coder::array<int8_T, 2U> &distags,
                    ::coder::SizeType col, ::coder::array<int64_T, 1U> &G_row_ptr,
                    ::coder::array<int32_T, 1U> &G_col_ind, int32_T *G_ncols)
{
  uint64_T leid;
  ::coder::SizeType b_iv;
  ::coder::SizeType c_i;
  ::coder::SizeType etable;
  ::coder::SizeType iend;
  ::coder::SizeType loop_ub;
  ::coder::SizeType ndn;
  ::coder::SizeType npe;
  ::coder::SizeType offset;
  ::coder::SizeType start;
  ::coder::SizeType v;
  //  Compute conn graph
  ndn = -1;
  for (int64_T i = mesh_node2elems_row_ptr[nid - 1];
       i < mesh_node2elems_row_ptr[nid]; i++) {
    etable =
        (
            mesh_teids[mesh_node2elems_col_ind[i - 1] -
                       1] &
            255UL) -
        1;
    leid =
        mesh_teids[mesh_node2elems_col_ind[i - 1] - 1] >>
        8;
    npe = mesh_elemtables[etable].conn.size(1);
    for (::coder::SizeType j{0}; j < npe; j++) {
      v = mesh_elemtables[etable]
              .conn[j + mesh_elemtables[etable].conn.size(1) *
                            (leid - 1)] -
          1;
      if ((distags[(col + distags.size(1) * v) - 1] != 0) && (!visited[v])) {
        ndn++;
        iwork[ndn] = v + 1;
        visited[v] = true;
      }
    }
  }
  //  Reset visited tags
  for (::coder::SizeType b_i{0}; b_i <= ndn; b_i++) {
    visited[iwork[b_i] - 1] = false;
  }
  G_row_ptr.set_size(ndn + 2);
  for (c_i = 0; c_i <= ndn + 1; c_i++) {
    G_row_ptr[c_i] = 0L;
  }
  G_row_ptr[0] = 1L;
  //  Use G.col_ind as buffer to compute inverse index mapping
  if (ndn + 1 < 1) {
    loop_ub = -1;
  } else {
    loop_ub = ndn;
  }
  G_col_ind.set_size(loop_ub + 1);
  for (c_i = 0; c_i <= loop_ub; c_i++) {
    G_col_ind[c_i] = iwork[c_i];
  }
  for (::coder::SizeType b_i{0}; b_i <= ndn; b_i++) {
    iwork[G_col_ind[b_i] - 1] = b_i + 1;
  }
  //  Determine total edges
  for (int64_T i = mesh_node2elems_row_ptr[nid - 1];
       i < mesh_node2elems_row_ptr[nid]; i++) {
    uint64_T c_tmp;
    ::coder::SizeType b_leid;
    c_tmp =
        mesh_teids[mesh_node2elems_col_ind[i - 1] - 1] &
        255UL;
    leid =
        mesh_teids[mesh_node2elems_col_ind[i - 1] - 1] >>
        8;
    b_leid = (
        mesh_teids[mesh_node2elems_col_ind[i - 1] - 1] >>
        8);
    npe = mesh_elemtables[c_tmp - 1].conn.size(1) - 1;
    for (::coder::SizeType j{0}; j <= npe; j++) {
      v = mesh_elemtables[c_tmp - 1]
              .conn[j + mesh_elemtables[c_tmp - 1]
                                .conn.size(1) *
                            (leid - 1)] -
          1;
      if (distags[(col + distags.size(1) * v) - 1] != 0) {
        b_iv = iwork[v];
        for (::coder::SizeType k{0}; k <= npe; k++) {
          if ((k != j) &&
              (distags
                   [(col +
                     distags.size(1) *
                         (mesh_elemtables[c_tmp - 1].conn
                              [k +
                               mesh_elemtables[c_tmp - 1]
                                       .conn.size(1) *
                                   (b_leid - 1)] -
                          1)) -
                    1] != 0)) {
            G_row_ptr[b_iv] = G_row_ptr[b_iv] + 1L;
          }
        }
      }
    }
  }
  for (::coder::SizeType b_i{0}; b_i <= ndn; b_i++) {
    G_row_ptr[b_i + 1] = G_row_ptr[b_i + 1] + G_row_ptr[b_i];
  }
  G_col_ind.set_size(
      (G_row_ptr[G_row_ptr.size(0) - 1] - 1L));
  for (int64_T i = mesh_node2elems_row_ptr[nid - 1];
       i < mesh_node2elems_row_ptr[nid]; i++) {
    etable =
        (
            mesh_teids[mesh_node2elems_col_ind[i - 1] -
                       1] &
            255UL) -
        1;
    leid =
        mesh_teids[mesh_node2elems_col_ind[i - 1] - 1] >>
        8;
    npe = mesh_elemtables[etable].conn.size(1) - 1;
    for (::coder::SizeType j{0}; j <= npe; j++) {
      v = mesh_elemtables[etable]
              .conn[j + mesh_elemtables[etable].conn.size(1) *
                            (leid - 1)] -
          1;
      if (distags[(col + distags.size(1) * v) - 1] != 0) {
        b_iv = iwork[v] - 1;
        for (::coder::SizeType k{0}; k <= npe; k++) {
          if (k != j) {
            ::coder::SizeType v0;
            v0 = mesh_elemtables[etable]
                     .conn[k + mesh_elemtables[etable].conn.size(1) *
                                   (leid - 1)] -
                 1;
            if (distags[(col + distags.size(1) * v0) - 1] != 0) {
              G_col_ind[(G_row_ptr[b_iv]) - 1] = iwork[v0];
              G_row_ptr[b_iv] = G_row_ptr[b_iv] + 1L;
            }
          }
        }
      }
    }
  }
  for (::coder::SizeType b_i{ndn + 1}; b_i >= 2; b_i--) {
    G_row_ptr[b_i - 1] = G_row_ptr[b_i - 2];
  }
  G_row_ptr[0] = 1L;
  iend = G_row_ptr.size(0);
  for (::coder::SizeType b_i{0}; b_i <= iend - 2; b_i++) {
    int64_T b_j;
    boolean_T ascend;
    boolean_T exitg1;
    ascend = true;
    b_j = G_row_ptr[b_i] + 1L;
    exitg1 = false;
    while ((!exitg1) && (b_j <= G_row_ptr[b_i + 1] - 1L)) {
      if (G_col_ind[b_j - 1] <
          G_col_ind[(b_j - 1L) - 1]) {
        ascend = false;
        exitg1 = true;
      } else {
        b_j++;
      }
    }
    if (!ascend) {
      //  sort in place
      m2cSort(G_col_ind, (G_row_ptr[b_i]),
              (G_row_ptr[b_i + 1] - 1L));
    }
  }
  offset = 0;
  start = 0;
  c_i = G_row_ptr.size(0);
  for (::coder::SizeType b_i{0}; b_i <= c_i - 2; b_i++) {
    int64_T i1;
    ::coder::SizeType i2;
    if (offset != 0) {
      G_col_ind[start - offset] = G_col_ind[start];
    }
    loop_ub = start + 2;
    i1 = G_row_ptr[b_i + 1];
    i2 = i1 - 1;
    for (::coder::SizeType j{loop_ub}; j <= i2; j++) {
      ::coder::SizeType i3;
      ::coder::SizeType i4;
      i3 = G_col_ind[j - 1];
      i4 = j - offset;
      if (i3 == G_col_ind[i4 - 2]) {
        offset++;
      } else if (offset != 0) {
        G_col_ind[i4 - 1] = i3;
      }
    }
    start = i1 - 1;
    G_row_ptr[b_i + 1] = i1 - offset;
  }
  if (offset != 0) {
    ::coder::SizeType newlen;
    newlen = G_col_ind.size(0) - offset;
    if (newlen < 1) {
      newlen = 0;
    }
    G_col_ind.set_size(newlen);
  }
  *G_ncols = ndn + 1;
}

static void
compute_measure_kernel(const ::coder::array<real_T, 2U> &mesh_coords,
                       const ::coder::array<ConnData, 1U> &mesh_elemtables,
                       const ::coder::array<uint64_T, 1U> &mesh_teids,
                       ::coder::array<real_T, 1U> &m)
{
  ::coder::array<SfeObject, 1U> sfes_;
  ::coder::array<real_T, 2U> xs_;
  SfeObject sfe;
  ::coder::SizeType iend;
  ::coder::SizeType istart;
  ::coder::SizeType loop_ub;
  ::coder::SizeType nthreads;
  ::coder::SizeType u1;
  sfe.etypes[0] = 0;
  sfe.nnodes[0] = 0;
  sfe.etypes[1] = 0;
  sfe.nnodes[1] = 0;
  sfe.geom_dim = 0;
  sfe.topo_dim = 0;
  sfe.facetid = 0;
  sfe.nqp = 0;
  sfe.ws.set_size(0);
  sfe.cs.set_size(::coder::SizeType(0), ::coder::SizeType(0));
  sfe.shapes_sol.set_size(::coder::SizeType(0), ::coder::SizeType(0));
  sfe.shapes_geom.set_size(::coder::SizeType(0), ::coder::SizeType(0));
  sfe.derivs_sol.set_size(::coder::SizeType(0), ::coder::SizeType(0), 0);
  sfe.derivs_geom.set_size(::coder::SizeType(0), ::coder::SizeType(0), 0);
  sfe.cs_phy.set_size(::coder::SizeType(0), ::coder::SizeType(0));
  sfe.grads_sol.set_size(::coder::SizeType(0), ::coder::SizeType(0));
  sfe.grads_geom.set_size(::coder::SizeType(0), ::coder::SizeType(0));
  sfe.jacTs.set_size(0, 3);
  sfe.wdetJ.set_size(0);
  sfe.dwork1.set_size(::coder::SizeType(0), ::coder::SizeType(0));
  sfe.dwork2.set_size(::coder::SizeType(0), ::coder::SizeType(0));
  sfe.xswork.set_size(::coder::SizeType(0), ::coder::SizeType(0));
  sfe.iwork.set_size(::coder::SizeType(0), ::coder::SizeType(0));
  sfes_.set_size(mesh_elemtables.size(0));
  loop_ub = mesh_elemtables.size(0);
  for (::coder::SizeType i{0}; i < loop_ub; i++) {
    sfes_[i] = sfe;
  }
#pragma omp single
  { // single
    m.set_size(mesh_teids.size(0));
  } // single
  nthreads = 1;
#ifdef _OPENMP
  nthreads = omp_get_num_threads();
#endif // _OPENMP
  if (nthreads == 1) {
    istart = 0;
    iend = mesh_teids.size(0);
  } else {
    ::coder::SizeType b_remainder;
    ::coder::SizeType chunk;
    ::coder::SizeType threadID;
    threadID = 0;
#ifdef _OPENMP
    threadID = omp_get_thread_num();
#endif // _OPENMP
    chunk = mesh_teids.size(0) / nthreads;
    b_remainder = mesh_teids.size(0) - nthreads * chunk;
    u1 = threadID;
    if (b_remainder <= threadID) {
      u1 = b_remainder;
    }
    istart = threadID * chunk + u1;
    iend = (istart + chunk) + (threadID < b_remainder);
  }
  //  Loop begins here
  for (::coder::SizeType b_i{istart + 1}; b_i <= iend; b_i++) {
    uint64_T c;
    ::coder::SizeType eid;
    c = mesh_teids[b_i - 1] & 255UL;
    eid = (mesh_teids[b_i - 1] >> 8);
    loop_ub = mesh_coords.size(1);
    xs_.set_size(mesh_elemtables[c - 1].conn.size(1),
                 mesh_coords.size(1));
    u1 = mesh_elemtables[c - 1].conn.size(1);
    for (::coder::SizeType i{0}; i < u1; i++) {
      for (::coder::SizeType i1{0}; i1 < loop_ub; i1++) {
        xs_[i1 + xs_.size(1) * i] = mesh_coords
            [i1 +
             mesh_coords.size(1) *
                 (mesh_elemtables[c - 1]
                      .conn[i + mesh_elemtables[c - 1]
                                        .conn.size(1) *
                                    (eid - 1)] -
                  1)];
      }
    }
    if (sfes_[c - 1].etypes[0] <= 0) {
      sfe_init(&sfes_[c - 1],
               mesh_elemtables[c - 1].etype, xs_);
    } else {
      sfe_init(&sfes_[c - 1], xs_);
    }
    m[b_i - 1] = coder::sum(sfes_[c - 1].wdetJ);
  }
}

static void compute_meshsizes_kernel(
    const ::coder::array<real_T, 2U> &mesh_coords,
    const ::coder::array<ConnData, 1U> &mesh_elemtables,
    const ::coder::array<uint64_T, 1U> &mesh_teids,
    const ::coder::array<int64_T, 1U> &mesh_node2nodes_row_ptr,
    const ::coder::array<int32_T, 1U> &mesh_node2nodes_col_ind,
    ::coder::array<real_T, 1U> &elemh, ::coder::array<real_T, 1U> &nodeh,
    real_T *globalh)
{
  real_T d;
  real_T hmax;
  ::coder::SizeType b_remainder;
  ::coder::SizeType c_i;
  ::coder::SizeType chunk;
  ::coder::SizeType dim;
  ::coder::SizeType iend;
  ::coder::SizeType istart;
  ::coder::SizeType loop_ub;
  ::coder::SizeType nthreads;
  ::coder::SizeType threadID;
  ::coder::SizeType u1;
  dim = mesh_coords.size(1);
  nthreads = 1;
#ifdef _OPENMP
  nthreads = omp_get_num_threads();
#endif // _OPENMP
  if (nthreads == 1) {
    istart = 0;
    iend = mesh_coords.size(0);
  } else {
    threadID = 0;
#ifdef _OPENMP
    threadID = omp_get_thread_num();
#endif // _OPENMP
    chunk = mesh_coords.size(0) / nthreads;
    b_remainder = mesh_coords.size(0) - nthreads * chunk;
    u1 = threadID;
    if (b_remainder <= threadID) {
      u1 = b_remainder;
    }
    istart = threadID * chunk + u1;
    iend = (istart + chunk) + (threadID < b_remainder);
  }
#pragma omp single
  { // single
    nodeh.set_size(mesh_coords.size(0));
  } // single
  //  Nodal h
  for (::coder::SizeType i{istart + 1}; i <= iend; i++) {
    int64_T b_i;
    nodeh[i - 1] = 0.0;
    b_i = mesh_node2nodes_row_ptr[i - 1];
    if (b_i <= mesh_node2nodes_row_ptr[i] - 1L) {
      c_i = i;
      loop_ub = mesh_coords.size(1);
    }
    for (int64_T k = mesh_node2nodes_row_ptr[i - 1];
         k < mesh_node2nodes_row_ptr[i]; k++) {
      real_T v_data[3];
      real_T len;
      for (u1 = 0; u1 < loop_ub; u1++) {
        v_data[u1] = mesh_coords[u1 + mesh_coords.size(1) *
                                          (mesh_node2nodes_col_ind
                                               [k - 1] -
                                           1)] -
                     mesh_coords[u1 + mesh_coords.size(1) * (c_i - 1)];
      }
      len = 0.0;
      for (::coder::SizeType ii{0}; ii < dim; ii++) {
        d = v_data[ii];
        len += std::sqrt(d * d);
      }
      nodeh[i - 1] = nodeh[i - 1] + len;
    }
    //  average
    nodeh[i - 1] =
        nodeh[i - 1] / static_cast<real_T>(mesh_node2nodes_row_ptr[i] - b_i);
  }
#pragma omp barrier
#pragma omp single
  { // single
    elemh.set_size(mesh_teids.size(0));
  } // single
  nthreads = 1;
#ifdef _OPENMP
  nthreads = omp_get_num_threads();
#endif // _OPENMP
  if (nthreads == 1) {
    istart = 0;
    iend = mesh_teids.size(0);
  } else {
    threadID = 0;
#ifdef _OPENMP
    threadID = omp_get_thread_num();
#endif // _OPENMP
    chunk = mesh_teids.size(0) / nthreads;
    b_remainder = mesh_teids.size(0) - nthreads * chunk;
    u1 = threadID;
    if (b_remainder <= threadID) {
      u1 = b_remainder;
    }
    istart = threadID * chunk + u1;
    iend = (istart + chunk) + (threadID < b_remainder);
  }
  //  Elemental h
  hmax = 0.0;
  for (::coder::SizeType i{istart + 1}; i <= iend; i++) {
    real_T h0;
    ::coder::SizeType etable;
    etable = (mesh_teids[i - 1] & 255UL) - 1;
    h0 = 0.0;
    u1 = mesh_elemtables[etable].conn.size(1) - 1;
    for (::coder::SizeType j{0}; j <= u1; j++) {
      h0 += nodeh[mesh_elemtables[etable].conn
                      [j +
                       mesh_elemtables[etable].conn.size(1) *
                           ((mesh_teids[i - 1] >> 8) - 1)] -
                  1];
    }
    //  Average
    d = h0 / static_cast<real_T>(mesh_elemtables[etable].conn.size(1));
    elemh[i - 1] = d;
    hmax = std::fmax(hmax, d);
  }
#pragma omp barrier
#pragma omp critical
  { // critical
  } // critical
  *globalh = hmax;
}

static void compute_meshsizes_kernel(
    const ::coder::array<real_T, 2U> &mesh_coords,
    const ::coder::array<ConnData, 1U> &mesh_elemtables,
    const ::coder::array<uint64_T, 1U> &mesh_teids,
    const ::coder::array<int64_T, 1U> &mesh_node2nodes_row_ptr,
    const ::coder::array<int32_T, 1U> &mesh_node2nodes_col_ind,
    const ::coder::array<real_T, 2U> &nrms, ::coder::array<real_T, 1U> &elemh,
    ::coder::array<real_T, 1U> &nodeh, real_T *globalh)
{
  real_T t_data[6];
  real_T d;
  real_T hmax;
  ::coder::SizeType b_remainder;
  ::coder::SizeType c_i;
  ::coder::SizeType chunk;
  ::coder::SizeType dim;
  ::coder::SizeType iend;
  ::coder::SizeType istart;
  ::coder::SizeType loop_ub;
  ::coder::SizeType nthreads;
  ::coder::SizeType threadID;
  ::coder::SizeType u1;
  dim = mesh_coords.size(1) - 1;
  nthreads = 1;
#ifdef _OPENMP
  nthreads = omp_get_num_threads();
#endif // _OPENMP
  if (nthreads == 1) {
    istart = 0;
    iend = mesh_coords.size(0);
  } else {
    threadID = 0;
#ifdef _OPENMP
    threadID = omp_get_thread_num();
#endif // _OPENMP
    chunk = mesh_coords.size(0) / nthreads;
    b_remainder = mesh_coords.size(0) - nthreads * chunk;
    u1 = threadID;
    if (b_remainder <= threadID) {
      u1 = b_remainder;
    }
    istart = threadID * chunk + u1;
    iend = (istart + chunk) + (threadID < b_remainder);
  }
#pragma omp single
  { // single
    nodeh.set_size(mesh_coords.size(0));
  } // single
  //  Nodal h
  for (::coder::SizeType i{istart + 1}; i <= iend; i++) {
    int64_T b_i;
    nodeh[i - 1] = 0.0;
    b_i = mesh_node2nodes_row_ptr[i - 1];
    if (b_i <= mesh_node2nodes_row_ptr[i] - 1L) {
      c_i = i;
      loop_ub = mesh_coords.size(1);
    }
    for (int64_T k = mesh_node2nodes_row_ptr[i - 1];
         k < mesh_node2nodes_row_ptr[i]; k++) {
      real_T v_data[3];
      real_T len;
      for (u1 = 0; u1 < loop_ub; u1++) {
        v_data[u1] = mesh_coords[u1 + mesh_coords.size(1) *
                                          (mesh_node2nodes_col_ind
                                               [k - 1] -
                                           1)] -
                     mesh_coords[u1 + mesh_coords.size(1) * (c_i - 1)];
      }
      if ((nrms.size(0) != 0) && (nrms.size(1) != 0)) {
        real_T u_data[2];
        //  Surface
        u1 = nrms.size(1) - 1;
        if (nrms.size(1) == 2) {
          //  2D
          t_data[0] = -nrms[nrms.size(1) * (i - 1) + 1];
          t_data[u1] = nrms[nrms.size(1) * (i - 1)];
        } else {
          real_T a;
          real_T a_tmp;
          real_T d1;
          boolean_T guard1{false};
          d = nrms[nrms.size(1) * (i - 1)];
          a_tmp = std::abs(d);
          d1 = nrms[nrms.size(1) * (i - 1) + 1];
          guard1 = false;
          if (a_tmp > std::abs(d1)) {
            real_T d2;
            d2 = nrms[nrms.size(1) * (i - 1) + 2];
            if (a_tmp > std::abs(d2)) {
              t_data[0] = -d * d1;
              t_data[u1] = 1.0 - d1 * d1;
              t_data[u1 * 2] = -nrms[nrms.size(1) * (i - 1) + 1] * d2;
            } else {
              guard1 = true;
            }
          } else {
            guard1 = true;
          }
          if (guard1) {
            t_data[0] = 1.0 - d * d;
            t_data[u1] = -nrms[nrms.size(1) * (i - 1)] *
                         nrms[nrms.size(1) * (i - 1) + 1];
            t_data[u1 * 2] = -nrms[nrms.size(1) * (i - 1)] *
                             nrms[nrms.size(1) * (i - 1) + 2];
          }
          a_tmp = t_data[u1 * 2];
          a = std::sqrt((t_data[0] * t_data[0] + t_data[u1] * t_data[u1]) +
                        a_tmp * a_tmp);
          t_data[0] /= a;
          t_data[u1] /= a;
          t_data[u1 * 2] /= a;
          //  cross
          a_tmp = nrms[nrms.size(1) * (i - 1) + 2];
          t_data[1] = t_data[u1 * 2] * d1 - t_data[u1] * a_tmp;
          t_data[u1 + 1] = t_data[0] * a_tmp - d * t_data[u1 * 2];
          t_data[u1 * 2 + 1] = d * t_data[u1] - t_data[0] * d1;
        }
        //  Project
        for (::coder::SizeType ii{0}; ii < dim; ii++) {
          u_data[ii] = 0.0;
          for (::coder::SizeType jj{0}; jj <= dim; jj++) {
            u_data[ii] += v_data[jj] * t_data[ii + u1 * jj];
          }
        }
        if (dim + 1 == 2) {
          len = std::abs(u_data[0]);
        } else {
          len = std::sqrt(u_data[0] * u_data[0] + u_data[1] * u_data[1]);
        }
      } else {
        len = 0.0;
        for (::coder::SizeType ii{0}; ii <= dim; ii++) {
          d = v_data[ii];
          len += std::sqrt(d * d);
        }
      }
      nodeh[i - 1] = nodeh[i - 1] + len;
    }
    //  average
    nodeh[i - 1] =
        nodeh[i - 1] / static_cast<real_T>(mesh_node2nodes_row_ptr[i] - b_i);
  }
#pragma omp barrier
#pragma omp single
  { // single
    elemh.set_size(mesh_teids.size(0));
  } // single
  nthreads = 1;
#ifdef _OPENMP
  nthreads = omp_get_num_threads();
#endif // _OPENMP
  if (nthreads == 1) {
    istart = 0;
    iend = mesh_teids.size(0);
  } else {
    threadID = 0;
#ifdef _OPENMP
    threadID = omp_get_thread_num();
#endif // _OPENMP
    chunk = mesh_teids.size(0) / nthreads;
    b_remainder = mesh_teids.size(0) - nthreads * chunk;
    u1 = threadID;
    if (b_remainder <= threadID) {
      u1 = b_remainder;
    }
    istart = threadID * chunk + u1;
    iend = (istart + chunk) + (threadID < b_remainder);
  }
  //  Elemental h
  hmax = 0.0;
  for (::coder::SizeType i{istart + 1}; i <= iend; i++) {
    real_T h0;
    ::coder::SizeType etable;
    etable = (mesh_teids[i - 1] & 255UL) - 1;
    h0 = 0.0;
    u1 = mesh_elemtables[etable].conn.size(1) - 1;
    for (::coder::SizeType j{0}; j <= u1; j++) {
      h0 += nodeh[mesh_elemtables[etable].conn
                      [j +
                       mesh_elemtables[etable].conn.size(1) *
                           ((mesh_teids[i - 1] >> 8) - 1)] -
                  1];
    }
    //  Average
    d = h0 / static_cast<real_T>(mesh_elemtables[etable].conn.size(1));
    elemh[i - 1] = d;
    hmax = std::fmax(hmax, d);
  }
#pragma omp barrier
#pragma omp critical
  { // critical
  } // critical
  *globalh = hmax;
}

static void
compute_nodal_alpha(const ::coder::array<real_T, 2U> &alphacell,
                    const ::coder::array<real_T, 2U> &mesh_coords,
                    const ::coder::array<int64_T, 1U> &mesh_node2elems_row_ptr,
                    const ::coder::array<int32_T, 1U> &mesh_node2elems_col_ind,
                    ::coder::array<real_T, 2U> &alphanode)
{
  ::coder::SizeType iend;
  ::coder::SizeType istart;
  ::coder::SizeType nrhs;
  ::coder::SizeType nthreads;
  nrhs = alphacell.size(1);
#pragma omp single
  { // single
    alphanode.set_size(mesh_coords.size(0), alphacell.size(1));
  } // single
  //  implicit barrier
  nthreads = 1;
#ifdef _OPENMP
  nthreads = omp_get_num_threads();
#endif // _OPENMP
  if (nthreads == 1) {
    istart = 0;
    iend = mesh_coords.size(0);
  } else {
    ::coder::SizeType b_remainder;
    ::coder::SizeType chunk;
    ::coder::SizeType threadID;
    ::coder::SizeType u1;
    threadID = 0;
#ifdef _OPENMP
    threadID = omp_get_thread_num();
#endif // _OPENMP
    chunk = mesh_coords.size(0) / nthreads;
    b_remainder = mesh_coords.size(0) - nthreads * chunk;
    u1 = threadID;
    if (b_remainder <= threadID) {
      u1 = b_remainder;
    }
    istart = threadID * chunk + u1;
    iend = (istart + chunk) + (threadID < b_remainder);
  }
  //  Compute nodal alpha values
  for (::coder::SizeType i{istart + 1}; i <= iend; i++) {
    for (::coder::SizeType k{0}; k < nrhs; k++) {
      real_T a;
      a = 0.0;
      for (int64_T j = mesh_node2elems_row_ptr[i - 1];
           j < mesh_node2elems_row_ptr[i]; j++) {
        real_T d;
        d = std::abs(
            alphacell[k + alphacell.size(1) *
                              (mesh_node2elems_col_ind[j -
                                                       1] -
                               1)]);
        if (a < d) {
          a = d;
        }
      }
      alphanode[k + alphanode.size(1) * (i - 1)] = a;
    }
  }
}

// compute_stencils_1d - Assembly of stencils for 1D
static void compute_stencils_1d(WlsMesh *mesh, ::coder::SizeType stclidx,
                                const real_T krings[3])
{
  int32_T bounds[2];
  bounds[0] = 0;
  //  set default upper bound
  bounds[1] = static_cast<::coder::SizeType>(std::round(2.0 * krings[0]));
  //  Ensure that max is no smaller than corresponding min
  if (bounds[1] <= 0) {
    bounds[1] = 0;
  }
  //  Serial
  compute_stencils_kernel_1d(mesh, stclidx, krings[0], bounds);
}

// compute_stencils_2d - Assembly of stencils for 2D
static void compute_stencils_2d(WlsMesh *mesh, ::coder::SizeType stclidx,
                                const real_T krings[3])
{
  int32_T bounds[2];
  bounds[0] = 0;
  //  set default upper bound
  bounds[1] = static_cast<::coder::SizeType>(
      std::round(12.0 * (krings[0] + 0.5) * (krings[0] + 1.0)));
  //  Ensure that max is no smaller than corresponding min
  if (bounds[1] <= 0) {
    bounds[1] = 0;
  }
  //  Serial
  compute_stencils_kernel_2d(mesh, stclidx, krings[0], bounds);
}

static void compute_stencils_kernel_1d(WlsMesh *mesh, ::coder::SizeType stclidx,
                                       real_T ring, const int32_T bounds[2])
{
  ::coder::array<uint64_T, 1U> hebuf_;
  ::coder::array<int32_T, 1U> ngbfs_;
  ::coder::array<int32_T, 1U> ngbvs_;
  ::coder::array<boolean_T, 1U> ftags_;
  ::coder::array<boolean_T, 1U> vtags_;
  ::coder::SizeType loop_ub;
  ::coder::SizeType n;
  int32_T nfaces;
  int32_T nverts;
  boolean_T hadoverflow;
  boolean_T overflow;
  boolean_T reflected;
  //  Determine total number of vertices
  if (mesh->stencils[stclidx - 1].vidmap.size(0) == 0) {
    n = mesh->coords.size(0);
  } else {
    n = mesh->stencils[stclidx - 1].vidmap.size(0);
  }
  //  Determine partitioning
  mesh->stencils[stclidx - 1].reflected.set_size(n);
  mesh->stencils[stclidx - 1].ngbverts.row_ptr.set_size(n + 1);
  mesh->stencils[stclidx - 1].ngbverts.row_ptr[0] = 1L;
  mesh->stencils[stclidx - 1].ngbverts.ncols = mesh->coords.size(0);
  mesh->stencils[stclidx - 1].ngbelems.row_ptr.set_size(n + 1);
  mesh->stencils[stclidx - 1].ngbelems.row_ptr[0] = 1L;
  mesh->stencils[stclidx - 1].ngbelems.ncols = mesh->teids.size(0);
  //  Assemble rowptrs
  for (::coder::SizeType i{0}; i < n; i++) {
    //  Get vertex ID
    mesh->stencils[stclidx - 1].ngbverts.row_ptr[i + 1] = bounds[1];
    mesh->stencils[stclidx - 1].ngbelems.row_ptr[i + 1] = bounds[1] << 1;
  }
  for (::coder::SizeType i{0}; i < n; i++) {
    mesh->stencils[stclidx - 1].ngbverts.row_ptr[i + 1] =
        mesh->stencils[stclidx - 1].ngbverts.row_ptr[i + 1] +
        mesh->stencils[stclidx - 1].ngbverts.row_ptr[i];
    mesh->stencils[stclidx - 1].ngbelems.row_ptr[i + 1] =
        mesh->stencils[stclidx - 1].ngbelems.row_ptr[i + 1] +
        mesh->stencils[stclidx - 1].ngbelems.row_ptr[i];
  }
  //  Allocate space for col_inds and actual nnzs
  mesh->stencils[stclidx - 1].ngbverts.col_ind.set_size((
      mesh->stencils[stclidx - 1].ngbverts.row_ptr[n] - 1L));
  mesh->stencils[stclidx - 1].ngbelems.col_ind.set_size((
      mesh->stencils[stclidx - 1].ngbelems.row_ptr[n] - 1L));
  mesh->iwork1.set_size(n);
  mesh->iwork2.set_size(n);
  vtags_.set_size(mesh->coords.size(0));
  loop_ub = mesh->coords.size(0);
  for (::coder::SizeType b_i{0}; b_i < loop_ub; b_i++) {
    vtags_[b_i] = false;
  }
  ftags_.set_size(mesh->teids.size(0));
  loop_ub = mesh->teids.size(0);
  for (::coder::SizeType b_i{0}; b_i < loop_ub; b_i++) {
    ftags_[b_i] = false;
  }
  hadoverflow = false;
  //  Loop begins
  for (::coder::SizeType i{0}; i < n; i++) {
    int64_T estart;
    int64_T vstart;
    ::coder::SizeType vid;
    vstart = mesh->stencils[stclidx - 1].ngbverts.row_ptr[i] - 1L;
    estart = mesh->stencils[stclidx - 1].ngbelems.row_ptr[i] - 1L;
    //  Prevent Coder from creating copies of rowptrs
    //  Get vertex ID
    if (mesh->stencils[stclidx - 1].vidmap.size(0) == 0) {
      vid = i + 1;
    } else {
      vid = mesh->stencils[stclidx - 1].vidmap[i];
    }
    obtain_nring_1d(mesh->elemtables, mesh->teids, mesh->sibhfs, mesh->v2hfid,
                    vid, ring, bounds[1], vtags_, ftags_, ngbvs_, &nverts,
                    ngbfs_, &nfaces, hebuf_, &reflected, &overflow);
    if ((!hadoverflow) && overflow) {
      m2cWarnMsgIdAndTxt("obtain_nring_2d:overflow",
                         "Buffers are too small to contain neighborhood");
      hadoverflow = true;
    }
    mesh->stencils[stclidx - 1].reflected[i] = reflected;
    mesh->iwork1[i] = nverts;
    for (int64_T j{1L}; j <= nverts; j++) {
      mesh->stencils[stclidx - 1]
          .ngbverts.col_ind[(vstart + j) - 1] =
          ngbvs_[j - 1];
    }
    mesh->iwork2[i] = nfaces;
    for (int64_T j{1L}; j <= nfaces; j++) {
      mesh->stencils[stclidx - 1]
          .ngbelems.col_ind[(estart + j) - 1] =
          ngbfs_[j - 1];
    }
  }
  //  Compress
  crs_compress(mesh->stencils[stclidx - 1].ngbverts.row_ptr,
               mesh->stencils[stclidx - 1].ngbverts.col_ind, mesh->iwork1);
  crs_compress(mesh->stencils[stclidx - 1].ngbelems.row_ptr,
               mesh->stencils[stclidx - 1].ngbelems.col_ind, mesh->iwork2);
}

static void compute_stencils_kernel_2d(WlsMesh *mesh, ::coder::SizeType stclidx,
                                       real_T ring, const int32_T bounds[2])
{
  ::coder::array<uint64_T, 1U> hebuf_;
  ::coder::array<int32_T, 1U> ngbfs_;
  ::coder::array<int32_T, 1U> ngbvs_;
  ::coder::array<boolean_T, 1U> ftags_;
  ::coder::array<boolean_T, 1U> vtags_;
  ::coder::SizeType loop_ub;
  ::coder::SizeType n;
  int32_T nfaces;
  int32_T nverts;
  boolean_T hadoverflow;
  boolean_T overflow;
  boolean_T reflected;
  //  Determine total number of vertices
  if (mesh->stencils[stclidx - 1].vidmap.size(0) == 0) {
    n = mesh->coords.size(0);
  } else {
    n = mesh->stencils[stclidx - 1].vidmap.size(0);
  }
  //  Determine partitioning
  mesh->stencils[stclidx - 1].reflected.set_size(n);
  mesh->stencils[stclidx - 1].ngbverts.row_ptr.set_size(n + 1);
  mesh->stencils[stclidx - 1].ngbverts.row_ptr[0] = 1L;
  mesh->stencils[stclidx - 1].ngbverts.ncols = mesh->coords.size(0);
  mesh->stencils[stclidx - 1].ngbelems.row_ptr.set_size(n + 1);
  mesh->stencils[stclidx - 1].ngbelems.row_ptr[0] = 1L;
  mesh->stencils[stclidx - 1].ngbelems.ncols = mesh->teids.size(0);
  //  Assemble rowptrs
  for (::coder::SizeType i{0}; i < n; i++) {
    //  Get vertex ID
    mesh->stencils[stclidx - 1].ngbverts.row_ptr[i + 1] = bounds[1];
    mesh->stencils[stclidx - 1].ngbelems.row_ptr[i + 1] = bounds[1] << 1;
  }
  for (::coder::SizeType i{0}; i < n; i++) {
    mesh->stencils[stclidx - 1].ngbverts.row_ptr[i + 1] =
        mesh->stencils[stclidx - 1].ngbverts.row_ptr[i + 1] +
        mesh->stencils[stclidx - 1].ngbverts.row_ptr[i];
    mesh->stencils[stclidx - 1].ngbelems.row_ptr[i + 1] =
        mesh->stencils[stclidx - 1].ngbelems.row_ptr[i + 1] +
        mesh->stencils[stclidx - 1].ngbelems.row_ptr[i];
  }
  //  Allocate space for col_inds and actual nnzs
  mesh->stencils[stclidx - 1].ngbverts.col_ind.set_size((
      mesh->stencils[stclidx - 1].ngbverts.row_ptr[n] - 1L));
  mesh->stencils[stclidx - 1].ngbelems.col_ind.set_size((
      mesh->stencils[stclidx - 1].ngbelems.row_ptr[n] - 1L));
  mesh->iwork1.set_size(n);
  mesh->iwork2.set_size(n);
  vtags_.set_size(mesh->coords.size(0));
  loop_ub = mesh->coords.size(0);
  for (::coder::SizeType b_i{0}; b_i < loop_ub; b_i++) {
    vtags_[b_i] = false;
  }
  ftags_.set_size(mesh->teids.size(0));
  loop_ub = mesh->teids.size(0);
  for (::coder::SizeType b_i{0}; b_i < loop_ub; b_i++) {
    ftags_[b_i] = false;
  }
  hadoverflow = false;
  //  Loop begins
  for (::coder::SizeType i{0}; i < n; i++) {
    int64_T estart;
    int64_T vstart;
    ::coder::SizeType vid;
    vstart = mesh->stencils[stclidx - 1].ngbverts.row_ptr[i] - 1L;
    estart = mesh->stencils[stclidx - 1].ngbelems.row_ptr[i] - 1L;
    //  Prevent Coder from creating copies of rowptrs
    //  Get vertex ID
    if (mesh->stencils[stclidx - 1].vidmap.size(0) == 0) {
      vid = i + 1;
    } else {
      vid = mesh->stencils[stclidx - 1].vidmap[i];
    }
    obtain_nring_2d(mesh->elemtables, mesh->teids, mesh->sibhfs, mesh->v2hfid,
                    vid, ring, bounds[1], vtags_, ftags_, mesh->bridges, ngbvs_,
                    &nverts, ngbfs_, &nfaces, hebuf_, &reflected, &overflow);
    if ((!hadoverflow) && overflow) {
      m2cWarnMsgIdAndTxt("obtain_nring_2d:overflow",
                         "Buffers are too small to contain neighborhood");
      hadoverflow = true;
    }
    mesh->stencils[stclidx - 1].reflected[i] = reflected;
    mesh->iwork1[i] = nverts;
    for (int64_T j{1L}; j <= nverts; j++) {
      mesh->stencils[stclidx - 1]
          .ngbverts.col_ind[(vstart + j) - 1] =
          ngbvs_[j - 1];
    }
    mesh->iwork2[i] = nfaces;
    for (int64_T j{1L}; j <= nfaces; j++) {
      mesh->stencils[stclidx - 1]
          .ngbelems.col_ind[(estart + j) - 1] =
          ngbfs_[j - 1];
    }
  }
  //  Compress
  crs_compress(mesh->stencils[stclidx - 1].ngbverts.row_ptr,
               mesh->stencils[stclidx - 1].ngbverts.col_ind, mesh->iwork1);
  crs_compress(mesh->stencils[stclidx - 1].ngbelems.row_ptr,
               mesh->stencils[stclidx - 1].ngbelems.col_ind, mesh->iwork2);
}

//  crsAx_kernel - Kernel function for evaluating A*x in each thread
static void crsAx_kernel(const ::coder::array<int64_T, 1U> &row_ptr,
                         const ::coder::array<int32_T, 1U> &col_ind,
                         const ::coder::array<real_T, 1U> &val,
                         const ::coder::array<real_T, 2U> &x,
                         ::coder::array<real_T, 2U> &b)
{
  ::coder::SizeType iend;
  ::coder::SizeType istart;
  ::coder::SizeType nrhs;
  ::coder::SizeType nthreads;
  nrhs = x.size(1) - 1;
  nthreads = 1;
#ifdef _OPENMP
  nthreads = omp_get_num_threads();
#endif // _OPENMP
  if (nthreads == 1) {
    istart = 0;
    iend = row_ptr.size(0) - 1;
  } else {
    ::coder::SizeType b_remainder;
    ::coder::SizeType chunk;
    ::coder::SizeType threadID;
    ::coder::SizeType u1;
    threadID = 0;
#ifdef _OPENMP
    threadID = omp_get_thread_num();
#endif // _OPENMP
    chunk = (row_ptr.size(0) - 1) / nthreads;
    b_remainder = (row_ptr.size(0) - nthreads * chunk) - 1;
    u1 = threadID;
    if (b_remainder <= threadID) {
      u1 = b_remainder;
    }
    istart = threadID * chunk + u1;
    iend = (istart + chunk) + (threadID < b_remainder);
  }
#pragma omp single
  { // single
    b.set_size(row_ptr.size(0) - 1, x.size(1));
  } // single
  //  Optimized for row-major
  for (::coder::SizeType i{istart + 1}; i <= iend; i++) {
    for (::coder::SizeType k{0}; k <= nrhs; k++) {
      b[k + b.size(1) * (i - 1)] = 0.0;
    }
    for (int64_T j = row_ptr[i - 1]; j < row_ptr[i]; j++) {
      for (::coder::SizeType k{0}; k <= nrhs; k++) {
        b[k + b.size(1) * (i - 1)] =
            b[k + b.size(1) * (i - 1)] +
            val[j - 1] *
                x[k + x.size(1) * (col_ind[j - 1] - 1)];
      }
    }
  }
}

// crsCompress - Compresses a CRS graph or matrix with (optionally) nnz in
static void crsCompress(CrsMatrix *A, const ::coder::array<int32_T, 1U> &nnzs)
{
  ::coder::array<real_T, 1U> val;
  ::coder::array<int32_T, 1U> col_ind;
  ::coder::array<int32_T, 1U> r;
  ::coder::SizeType istart;
  ::coder::SizeType j;
  ::coder::SizeType loop_ub;
  ::coder::SizeType nrows;
  ::coder::SizeType u1;
  r.set_size(A->col_ind.size(0));
  loop_ub = A->col_ind.size(0);
  for (u1 = 0; u1 < loop_ub; u1++) {
    r[u1] = A->col_ind[u1];
  }
  // Local buffers with patterns "(\w+)?col_ind(\w+)?|r" are legitimate
  // Local buffers with patterns "(\w+)?val(\w+)?" are legitimate
  nrows = A->row_ptr.size(0);
  m2cAssert(A->row_ptr.size(0) - 1 == nnzs.size(0),
            "Length of nnzs must be equal to the number of rows.");
  //  Get the size upon input
  j = 0;
  istart = (A->row_ptr[0]);
  for (::coder::SizeType i{0}; i <= nrows - 2; i++) {
    ::coder::SizeType b_j0;
    b_j0 = j;
    loop_ub = nnzs[i];
    u1 = (A->row_ptr[i + 1]) - istart;
    if (loop_ub <= u1) {
      u1 = loop_ub;
    }
    u1 = (istart + u1) - 1;
    for (::coder::SizeType k{istart}; k <= u1; k++) {
      real_T d;
      d = A->val[k - 1];
      if (d != 0.0) {
        //  Skip zero entries
        r[j] = r[k - 1];
        A->val[j] = d;
        j++;
      }
    }
    istart = (A->row_ptr[i + 1]);
    A->row_ptr[i + 1] = A->row_ptr[i] + (j - b_j0);
  }
  //  resize inplace
  if (j < 1) {
    u1 = 0;
  } else {
    u1 = j;
  }
  r.set_size(u1);
  if (j < 1) {
    loop_ub = 0;
  } else {
    loop_ub = j;
  }
  A->val.set_size(loop_ub);
  if (j > A->col_ind.size(0) / 2) {
    //  make copies
    col_ind.set_size(u1);
    for (::coder::SizeType i{0}; i < u1; i++) {
      col_ind[i] = r[i];
    }
    val.set_size(A->val.size(0));
    for (::coder::SizeType i{0}; i < loop_ub; i++) {
      val[i] = A->val[i];
    }
    r.set_size(col_ind.size(0));
    loop_ub = col_ind.size(0);
    for (u1 = 0; u1 < loop_ub; u1++) {
      r[u1] = col_ind[u1];
    }
    A->val.set_size(val.size(0));
    loop_ub = val.size(0);
    for (u1 = 0; u1 < loop_ub; u1++) {
      A->val[u1] = val[u1];
    }
  }
  A->col_ind.set_size(r.size(0));
  loop_ub = r.size(0);
  for (u1 = 0; u1 < loop_ub; u1++) {
    A->col_ind[u1] = r[u1];
  }
}

//  crsProdMatVec - Compute b=A*x for CRS matrix A and vector(s) x.
static inline
void crsProdMatVec(const ::coder::array<int64_T, 1U> &A_row_ptr,
                          const ::coder::array<int32_T, 1U> &A_col_ind,
                          const ::coder::array<real_T, 1U> &A_val,
                          const ::coder::array<real_T, 2U> &x,
                          ::coder::array<real_T, 2U> &b)
{
  crs_prod_mat_vec(A_row_ptr, A_col_ind, A_val, x, b);
}

// crs_compress - Compresses a padded CRS matrix given the nnz in each row.
static void crs_compress(::coder::array<int64_T, 1U> &row_ptr,
                         ::coder::array<int32_T, 1U> &col_ind,
                         const ::coder::array<int32_T, 1U> &nnzs)
{
  ::coder::array<int32_T, 1U> col_ind0;
  ::coder::SizeType a;
  ::coder::SizeType b_i;
  ::coder::SizeType istart;
  ::coder::SizeType j;
  ::coder::SizeType nrows;
  ::coder::SizeType u0;
  ::coder::SizeType u1;
  // Local buffers with patterns "(\w+)?col_ind(\w+)?|r" are legitimate
  nrows = row_ptr.size(0);
  m2cAssert(row_ptr.size(0) - 1 == nnzs.size(0),
            "Length of nnzs must be equal to the number of rows.");
  //  Get the size upon input
  a = col_ind.size(0);
  j = 0;
  istart = (row_ptr[0]);
  for (::coder::SizeType i{0}; i <= nrows - 2; i++) {
    ::coder::SizeType b_j0;
    ::coder::SizeType u1_tmp;
    u1_tmp = (row_ptr[i + 1]);
    b_j0 = j;
    u0 = nnzs[i];
    u1 = u1_tmp - istart;
    if (u0 <= u1) {
      u1 = u0;
    }
    b_i = (istart + u1) - 1;
    for (::coder::SizeType k{istart}; k <= b_i; k++) {
      col_ind[(j + k) - istart] = col_ind[k - 1];
    }
    j = ((j + b_i) - istart) + 1;
    istart = u1_tmp;
    row_ptr[i + 1] = row_ptr[i] + (j - b_j0);
  }
  //  resize inplace
  if (j < 1) {
    b_i = 0;
  } else {
    b_i = j;
  }
  col_ind.set_size(b_i);
  if (j > a / 2) {
    //  make copies
    col_ind0.set_size(col_ind.size(0));
    u0 = col_ind.size(0);
    for (u1 = 0; u1 < u0; u1++) {
      col_ind0[u1] = col_ind[u1];
    }
    for (::coder::SizeType i{0}; i < b_i; i++) {
      col_ind[i] = col_ind0[i];
    }
  }
}

//  crs_prod_mat_vec - Compute b=A*x for CRS matrix A and vector(s) x.
static inline
void crs_prod_mat_vec(const ::coder::array<int64_T, 1U> &A_rowptr,
                             const ::coder::array<int32_T, 1U> &A_colind,
                             const ::coder::array<real_T, 1U> &A_val,
                             const ::coder::array<real_T, 2U> &x,
                             ::coder::array<real_T, 2U> &b)
{
  crsAx_kernel(A_rowptr, A_colind, A_val, x, b);
}

// determine_rdnodes - Determine rank-deficient (RD) nodes
static inline
void determine_rdnodes(boolean_T fullrank,
                              ::coder::array<boolean_T, 1U> &rdtags,
                              ::coder::array<int32_T, 1U> &rdnodes)
{
  if (fullrank) {
    rdnodes.set_size(0);
  } else {
    ::coder::SizeType j;
    ::coder::SizeType n;
    ::coder::SizeType rds;
    rds = 0;
    n = rdtags.size(0) - 1;
    for (::coder::SizeType i{0}; i <= n; i++) {
      if (rdtags[i]) {
        rds++;
      }
    }
    rdnodes.set_size(rds);
    j = 0;
    for (::coder::SizeType i{0}; i <= n; i++) {
      if (rdtags[i]) {
        rdnodes[j] = i + 1;
        j++;
        rdtags[i] = false;
        //  reset
      }
    }
  }
}

static void extract_sub(::coder::SizeType n, const ::coder::array<int32_T, 1U> &crange,
                        const ::coder::array<int32_T, 1U> &eptr,
                        const ::coder::array<int32_T, 1U> &eind,
                        ::coder::array<int32_T, 1U> &iwork,
                        ::coder::array<boolean_T, 1U> &ntags,
                        ::coder::array<int32_T, 1U> &eptrloc,
                        ::coder::array<int32_T, 1U> &eindloc, int32_T *nnodes)
{
  ::coder::array<int32_T, 1U> l2g_;
  ::coder::SizeType j;
  ::coder::SizeType k;
  ::coder::SizeType m;
  ::coder::SizeType u0;
  ::coder::SizeType y;
  //  local function to extract and localize a subgraph
  m = crange.size(0) - 1;
  u0 = crange.size(0);
  k = n;
  if (u0 >= n) {
    k = u0;
  }
  if (iwork.size(0) < k) {
    u0 = crange.size(0);
    k = n;
    if (u0 >= n) {
      k = u0;
    }
    iwork.set_size(k);
  }
  eptrloc.set_size(crange.size(0) + 1);
  eptrloc[0] = 1;
  for (::coder::SizeType i{0}; i <= m; i++) {
    eptrloc[i + 1] = (eptrloc[i] + eptr[crange[i]]) - eptr[crange[i] - 1];
  }
  eindloc.set_size(eptrloc[crange.size(0)] - 1);
  for (::coder::SizeType i{0}; i <= m; i++) {
    ::coder::SizeType n0;
    n0 = eptrloc[i + 1] - eptrloc[i];
    for (j = 0; j < n0; j++) {
      eindloc[(eptrloc[i] + j) - 1] = eind[(eptr[crange[i] - 1] + j) - 1];
    }
    k = eptr[crange[i] - 1];
    u0 = eptr[crange[i]] - 1;
    for (j = k; j <= u0; j++) {
      ntags[eind[j - 1] - 1] = true;
    }
  }
  u0 = ntags.size(0);
  y = ntags[0];
  for (k = 2; k <= u0; k++) {
    if (u0 >= 2) {
      y += ntags[k - 1];
    }
  }
  l2g_.set_size(y);
  j = -1;
  for (::coder::SizeType i{0}; i < n; i++) {
    if (ntags[i]) {
      j++;
      l2g_[j] = i + 1;
      ntags[i] = false;
      //  reset
    }
  }
  //  global to local
  for (::coder::SizeType i{0}; i < y; i++) {
    iwork[l2g_[i] - 1] = i + 1;
  }
  k = eindloc.size(0);
  for (::coder::SizeType i{0}; i < k; i++) {
    eindloc[i] = iwork[eindloc[i] - 1];
  }
  *nnodes = y;
}

static void f_WlsObject(::coder::SizeType degree, b_WlsObject *wls)
{
  wls->nstpnts = 0;
  wls->degree = degree;
  wls->order = 0;
  wls->unimono = false;
  wls->interp0 = 0;
  wls->stride = 0;
  wls->us.set_size(0, 3);
  wls->origin.size[1] = 3;
  wls->origin.size[0] = 1;
  wls->origin.data[0] = 0.0;
  wls->origin.data[1] = 0.0;
  wls->origin.data[2] = 0.0;
  wls->rweights.set_size(0);
  wls->hs_inv.size[1] = 0;
  wls->hs_inv.size[0] = 1;
  wls->V.set_size(::coder::SizeType(0), ::coder::SizeType(0));
  wls->QR.set_size(::coder::SizeType(0), ::coder::SizeType(0));
  wls->rhs.set_size(::coder::SizeType(0), ::coder::SizeType(0));
  wls->nevpnts = 0;
  wls->nrows = 0;
  wls->ncols = 0;
  wls->rank = 0;
  wls->fullrank = false;
  wls->jpvt.set_size(0);
  wls->work.set_size(0);
  wls->rowmajor = true;
  wls->QRt.set_size(::coder::SizeType(0), ::coder::SizeType(0));
  wls->runtimes.size[0] = 0;
}

static real_T find_kth_shortest_dist(::coder::array<real_T, 1U> &arr, ::coder::SizeType k,
                                     ::coder::SizeType l, ::coder::SizeType r)
{
  real_T dist;
  real_T val;
  ::coder::SizeType i;
  ::coder::SizeType j;
  //  Find the kth smallest number in arr(l:r).
  if (k < l) {
    k = l;
  }
  if (k > r) {
    k = r;
  }
  val = arr[l - 1];
  i = l;
  j = r;
  while (i <= j) {
    real_T d;
    real_T d1;
    ::coder::SizeType exitg1;
    do {
      exitg1 = 0;
      d = arr[i - 1];
      if (d < val) {
        i++;
      } else {
        exitg1 = 1;
      }
    } while (exitg1 == 0);
    do {
      exitg1 = 0;
      d1 = arr[j - 1];
      if (d1 > val) {
        j--;
      } else {
        exitg1 = 1;
      }
    } while (exitg1 == 0);
    if (i <= j) {
      arr[i - 1] = d1;
      arr[j - 1] = d;
      i++;
      j--;
    }
  }
  if (k <= j) {
    dist = find_kth_shortest_dist(arr, k, l, j);
  } else if (k >= i) {
    dist = find_kth_shortest_dist(arr, k, i, r);
  } else {
    dist = val;
  }
  return dist;
}

//  gen_vander  Wrapper function for computing confluent Vandermonde matrix in
static void gen_vander(const ::coder::array<real_T, 2U> &us, ::coder::SizeType npoints,
                       ::coder::SizeType degree, ::coder::array<real_T, 2U> &V)
{
  switch (us.size(1)) {
  case 1: {
    ::coder::SizeType b_npoints;
    ::coder::SizeType i;
    ::coder::SizeType i1;
    ::coder::SizeType i2;
    boolean_T b;
    boolean_T b1;
    b_npoints = npoints - 1;
    m2cAssert(us.size(1) == 1, "");
    //  Handle input arguments
    if (npoints == 0) {
      b_npoints = us.size(0) - 1;
    } else {
      m2cAssert(npoints <= us.size(0), "Input us is too small.");
    }
    m2cAssert(degree >= 0, "Degree must be nonnegative");
    //  Number of row blocks
    V.set_size(degree + 1, us.size(0));
    //  Compute rows corresponding to function values
    if (degree != 0) {
      b = true;
      b1 = ((us.size(1) <= 0) || (us.size(0) <= 0));
      i = us.size(1) * us.size(0);
      i1 = 0;
      for (::coder::SizeType iPnt{0}; iPnt <= b_npoints; iPnt++) {
        if (b1 || (iPnt >= i)) {
          i1 = 0;
          b = true;
        } else if (b) {
          b = false;
          i1 = iPnt % us.size(0) * us.size(1) + iPnt / us.size(0);
        } else {
          i2 = us.size(1) * us.size(0) - 1;
          if (i1 > MAX_int32_T - us.size(1)) {
            i1 = iPnt % us.size(0) * us.size(1) + iPnt / us.size(0);
          } else {
            i1 += us.size(1);
            if (i1 > i2) {
              i1 -= i2;
            }
          }
        }
        V[iPnt] = 1.0;
        V[iPnt + V.size(1)] = us[i1];
      }
    } else {
      for (::coder::SizeType iPnt{0}; iPnt <= b_npoints; iPnt++) {
        V[iPnt] = 1.0;
      }
    }
    i = degree + 1;
    for (::coder::SizeType ii{2}; ii <= i; ii++) {
      b = true;
      b1 = ((us.size(1) <= 0) || (us.size(0) <= 0));
      i1 = us.size(1) * us.size(0);
      i2 = 0;
      for (::coder::SizeType iPnt{0}; iPnt <= b_npoints; iPnt++) {
        if (b1 || (iPnt >= i1)) {
          i2 = 0;
          b = true;
        } else if (b) {
          b = false;
          i2 = iPnt % us.size(0) * us.size(1) + iPnt / us.size(0);
        } else {
          ::coder::SizeType i3;
          i3 = us.size(1) * us.size(0) - 1;
          if (i2 > MAX_int32_T - us.size(1)) {
            i2 = iPnt % us.size(0) * us.size(1) + iPnt / us.size(0);
          } else {
            i2 += us.size(1);
            if (i2 > i3) {
              i2 -= i3;
            }
          }
        }
        V[iPnt + V.size(1) * (ii - 1)] =
            V[iPnt + V.size(1) * (ii - 2)] * us[i2];
      }
    }
    //  Add row blocks corresponding to kth derivatives
  } break;
  case 2:
    gen_vander_2d(us, npoints, degree, V);
    break;
  default:
    gen_vander_3d(us, npoints, degree, V);
    break;
  }
}

//  gen_vander  Wrapper function for computing confluent Vandermonde matrix in
static void gen_vander(const ::coder::array<real_T, 2U> &us, ::coder::SizeType npoints,
                       ::coder::SizeType degree,
                       const ::coder::array<real_T, 1U> &weights,
                       ::coder::array<real_T, 2U> &V)
{
  switch (us.size(1)) {
  case 1: {
    ::coder::SizeType i;
    ::coder::SizeType i1;
    ::coder::SizeType i2;
    boolean_T b;
    boolean_T b1;
    m2cAssert(us.size(1) == 1, "");
    //  Handle input arguments
    m2cAssert(npoints <= us.size(0), "Input us is too small.");
    m2cAssert(degree >= 0, "Degree must be nonnegative");
    //  Number of row blocks
    V.set_size(degree + 1, us.size(0));
    //  Compute rows corresponding to function values
    if (weights.size(0) == 0) {
      if (degree != 0) {
        b = true;
        b1 = ((us.size(1) <= 0) || (us.size(0) <= 0));
        i = us.size(1) * us.size(0);
        i1 = 0;
        for (::coder::SizeType iPnt{0}; iPnt < npoints; iPnt++) {
          if (b1 || (iPnt >= i)) {
            i1 = 0;
            b = true;
          } else if (b) {
            b = false;
            i1 = iPnt % us.size(0) * us.size(1) + iPnt / us.size(0);
          } else {
            i2 = us.size(1) * us.size(0) - 1;
            if (i1 > MAX_int32_T - us.size(1)) {
              i1 = iPnt % us.size(0) * us.size(1) + iPnt / us.size(0);
            } else {
              i1 += us.size(1);
              if (i1 > i2) {
                i1 -= i2;
              }
            }
          }
          V[iPnt] = 1.0;
          V[iPnt + V.size(1)] = us[i1];
        }
      } else {
        for (::coder::SizeType iPnt{0}; iPnt < npoints; iPnt++) {
          V[iPnt] = 1.0;
        }
      }
    } else if (degree != 0) {
      b = true;
      b1 = ((us.size(1) <= 0) || (us.size(0) <= 0));
      i = us.size(1) * us.size(0);
      i1 = 0;
      for (::coder::SizeType iPnt{0}; iPnt < npoints; iPnt++) {
        if (b1 || (iPnt >= i)) {
          i1 = 0;
          b = true;
        } else if (b) {
          b = false;
          i1 = iPnt % us.size(0) * us.size(1) + iPnt / us.size(0);
        } else {
          i2 = us.size(1) * us.size(0) - 1;
          if (i1 > MAX_int32_T - us.size(1)) {
            i1 = iPnt % us.size(0) * us.size(1) + iPnt / us.size(0);
          } else {
            i1 += us.size(1);
            if (i1 > i2) {
              i1 -= i2;
            }
          }
        }
        V[iPnt] = weights[iPnt];
        V[iPnt + V.size(1)] = us[i1] * weights[iPnt];
      }
    } else {
      for (::coder::SizeType iPnt{0}; iPnt < npoints; iPnt++) {
        V[iPnt] = weights[iPnt];
      }
    }
    i = degree + 1;
    for (::coder::SizeType ii{2}; ii <= i; ii++) {
      b = true;
      b1 = ((us.size(1) <= 0) || (us.size(0) <= 0));
      i1 = us.size(1) * us.size(0);
      i2 = 0;
      for (::coder::SizeType iPnt{0}; iPnt < npoints; iPnt++) {
        if (b1 || (iPnt >= i1)) {
          i2 = 0;
          b = true;
        } else if (b) {
          b = false;
          i2 = iPnt % us.size(0) * us.size(1) + iPnt / us.size(0);
        } else {
          ::coder::SizeType i3;
          i3 = us.size(1) * us.size(0) - 1;
          if (i2 > MAX_int32_T - us.size(1)) {
            i2 = iPnt % us.size(0) * us.size(1) + iPnt / us.size(0);
          } else {
            i2 += us.size(1);
            if (i2 > i3) {
              i2 -= i3;
            }
          }
        }
        V[iPnt + V.size(1) * (ii - 1)] =
            V[iPnt + V.size(1) * (ii - 2)] * us[i2];
      }
    }
    //  Add row blocks corresponding to kth derivatives
  } break;
  case 2:
    gen_vander_2d(us, npoints, degree, weights, V);
    break;
  default:
    gen_vander_3d(us, npoints, degree, weights, V);
    break;
  }
}

//  gen_vander_2d  Generate generalized/confluent Vandermonde matrix in 2D.
static void gen_vander_2d(const ::coder::array<real_T, 2U> &us, ::coder::SizeType npoints,
                          ::coder::SizeType degree,
                          const ::coder::array<real_T, 1U> &weights,
                          ::coder::array<real_T, 2U> &V)
{
  ::coder::SizeType b_degree;
  ::coder::SizeType c;
  if (npoints > us.size(0)) {
    m2cErrMsgIdAndTxt("wlslib:BufferTooSmall", "Input us is too small.");
  }
  //  Number of row blocks
  if (degree >= 0) {
    b_degree = (degree + 1) * (degree + 2) / 2;
  } else {
    b_degree = (1 - degree) * (1 - degree);
  }
  V.set_size(b_degree, us.size(0));
  //  compute 0th order generalized Vandermonde matrix
  if (weights.size(0) == 0) {
    if (degree != 0) {
      for (::coder::SizeType iPnt{0}; iPnt < npoints; iPnt++) {
        V[iPnt] = 1.0;
        V[iPnt + V.size(1)] = us[us.size(1) * iPnt];
        V[iPnt + V.size(1) * 2] = us[us.size(1) * iPnt + 1];
      }
    } else {
      for (::coder::SizeType iPnt{0}; iPnt < npoints; iPnt++) {
        V[iPnt] = 1.0;
      }
    }
  } else if (degree != 0) {
    for (::coder::SizeType iPnt{0}; iPnt < npoints; iPnt++) {
      V[iPnt] = weights[iPnt];
      V[iPnt + V.size(1)] = us[us.size(1) * iPnt] * weights[iPnt];
      V[iPnt + V.size(1) * 2] = us[us.size(1) * iPnt + 1] * weights[iPnt];
    }
  } else {
    for (::coder::SizeType iPnt{0}; iPnt < npoints; iPnt++) {
      V[iPnt] = weights[iPnt];
    }
  }
  c = 3;
  if (degree < 0) {
    b_degree = -degree;
  } else {
    b_degree = degree;
  }
  for (::coder::SizeType deg{2}; deg <= b_degree; deg++) {
    for (::coder::SizeType j{0}; j < deg; j++) {
      for (::coder::SizeType iPnt{0}; iPnt < npoints; iPnt++) {
        V[iPnt + V.size(1) * c] =
            V[iPnt + V.size(1) * (c - deg)] * us[us.size(1) * iPnt];
      }
      c++;
    }
    for (::coder::SizeType iPnt{0}; iPnt < npoints; iPnt++) {
      V[iPnt + V.size(1) * c] =
          V[iPnt + V.size(1) * ((c - deg) - 1)] * us[us.size(1) * iPnt + 1];
    }
    c++;
  }
  //  Compute the bi-degree terms if degree<0
  b_degree = -degree;
  for (::coder::SizeType deg{b_degree}; deg >= 1; deg--) {
    for (::coder::SizeType k{0}; k < deg; k++) {
      for (::coder::SizeType iPnt{0}; iPnt < npoints; iPnt++) {
        V[iPnt + V.size(1) * c] =
            V[iPnt + V.size(1) * (c - deg)] * us[us.size(1) * iPnt];
      }
      c++;
    }
  }
  //  compute higher order confluent Vandermonde matrix blocks incrementally
}

//  gen_vander_2d  Generate generalized/confluent Vandermonde matrix in 2D.
static void gen_vander_2d(const ::coder::array<real_T, 2U> &us, ::coder::SizeType npoints,
                          ::coder::SizeType degree, ::coder::array<real_T, 2U> &V)
{
  ::coder::SizeType b_degree;
  ::coder::SizeType c;
  if (npoints == 0) {
    npoints = us.size(0);
  } else if (npoints > us.size(0)) {
    m2cErrMsgIdAndTxt("wlslib:BufferTooSmall", "Input us is too small.");
  }
  //  Number of row blocks
  if (degree >= 0) {
    b_degree = (degree + 1) * (degree + 2) / 2;
  } else {
    b_degree = (1 - degree) * (1 - degree);
  }
  V.set_size(b_degree, us.size(0));
  //  compute 0th order generalized Vandermonde matrix
  if (degree != 0) {
    for (::coder::SizeType iPnt{0}; iPnt < npoints; iPnt++) {
      V[iPnt] = 1.0;
      V[iPnt + V.size(1)] = us[us.size(1) * iPnt];
      V[iPnt + V.size(1) * 2] = us[us.size(1) * iPnt + 1];
    }
  } else {
    for (::coder::SizeType iPnt{0}; iPnt < npoints; iPnt++) {
      V[iPnt] = 1.0;
    }
  }
  c = 3;
  if (degree < 0) {
    b_degree = -degree;
  } else {
    b_degree = degree;
  }
  for (::coder::SizeType deg{2}; deg <= b_degree; deg++) {
    for (::coder::SizeType j{0}; j < deg; j++) {
      for (::coder::SizeType iPnt{0}; iPnt < npoints; iPnt++) {
        V[iPnt + V.size(1) * c] =
            V[iPnt + V.size(1) * (c - deg)] * us[us.size(1) * iPnt];
      }
      c++;
    }
    for (::coder::SizeType iPnt{0}; iPnt < npoints; iPnt++) {
      V[iPnt + V.size(1) * c] =
          V[iPnt + V.size(1) * ((c - deg) - 1)] * us[us.size(1) * iPnt + 1];
    }
    c++;
  }
  //  Compute the bi-degree terms if degree<0
  b_degree = -degree;
  for (::coder::SizeType deg{b_degree}; deg >= 1; deg--) {
    for (::coder::SizeType k{0}; k < deg; k++) {
      for (::coder::SizeType iPnt{0}; iPnt < npoints; iPnt++) {
        V[iPnt + V.size(1) * c] =
            V[iPnt + V.size(1) * (c - deg)] * us[us.size(1) * iPnt];
      }
      c++;
    }
  }
  //  compute higher order confluent Vandermonde matrix blocks incrementally
}

//  gen_vander_3d  Generate generalized/confluent Vandermonde matrix in 3D.
static void gen_vander_3d(const ::coder::array<real_T, 2U> &us, ::coder::SizeType npoints,
                          ::coder::SizeType degree,
                          const ::coder::array<real_T, 1U> &weights,
                          ::coder::array<real_T, 2U> &V)
{
  ::coder::SizeType b_degree;
  ::coder::SizeType c;
  ::coder::SizeType d;
  ::coder::SizeType deg;
  ::coder::SizeType i;
  if (npoints > us.size(0)) {
    m2cErrMsgIdAndTxt("wlslib:BufferTooSmall", "Input us is too small.");
  }
  //  Allocate storage for V
  if (degree >= 0) {
    b_degree = (degree + 1) * (degree + 2) * (degree + 3) / 6;
  } else {
    b_degree = (1 - degree) * (1 - degree) * (1 - degree);
  }
  V.set_size(b_degree, us.size(0));
  //  compute 0th order generalized Vandermonde matrix
  if (weights.size(0) == 0) {
    if (degree != 0) {
      for (::coder::SizeType iPnt{0}; iPnt < npoints; iPnt++) {
        V[iPnt] = 1.0;
        V[iPnt + V.size(1)] = us[us.size(1) * iPnt];
        V[iPnt + V.size(1) * 2] = us[us.size(1) * iPnt + 1];
        V[iPnt + V.size(1) * 3] = us[us.size(1) * iPnt + 2];
      }
    } else {
      for (::coder::SizeType iPnt{0}; iPnt < npoints; iPnt++) {
        V[iPnt] = 1.0;
      }
    }
  } else if (degree != 0) {
    for (::coder::SizeType iPnt{0}; iPnt < npoints; iPnt++) {
      V[iPnt] = weights[iPnt];
      V[iPnt + V.size(1)] = us[us.size(1) * iPnt] * weights[iPnt];
      V[iPnt + V.size(1) * 2] = us[us.size(1) * iPnt + 1] * weights[iPnt];
      V[iPnt + V.size(1) * 3] = us[us.size(1) * iPnt + 2] * weights[iPnt];
    }
  } else {
    for (::coder::SizeType iPnt{0}; iPnt < npoints; iPnt++) {
      V[iPnt] = weights[iPnt];
    }
  }
  c = 4;
  d = 3;
  if (degree < 0) {
    i = -degree;
  } else {
    i = degree;
  }
  for (deg = 2; deg <= i; deg++) {
    //  Within each level, use convention of Pascal triangle with x^deg at peak
    for (::coder::SizeType j{0}; j < deg; j++) {
      for (::coder::SizeType iPnt{0}; iPnt < npoints; iPnt++) {
        V[iPnt + V.size(1) * c] =
            V[iPnt + V.size(1) * (c - d)] * us[us.size(1) * iPnt];
      }
      c++;
    }
    for (::coder::SizeType iPnt{0}; iPnt < npoints; iPnt++) {
      V[iPnt + V.size(1) * c] =
          V[iPnt + V.size(1) * ((c - d) - 1)] * us[us.size(1) * iPnt + 1];
    }
    c++;
    for (::coder::SizeType j{0}; j < d; j++) {
      for (::coder::SizeType iPnt{0}; iPnt < npoints; iPnt++) {
        V[iPnt + V.size(1) * c] = V[iPnt + V.size(1) * (((c - d) - deg) - 1)] *
                                  us[us.size(1) * iPnt + 2];
      }
      c++;
    }
    d = (d + deg) + 1;
  }
  //  Compute the tri-degree terms if degree<0
  if (degree < 0) {
    ::coder::SizeType cornerTriangle;
    ::coder::SizeType excess;
    ::coder::SizeType maxLayers;
    ::coder::SizeType nTermsInLayer;
    deg = -degree;
    maxLayers = -degree * 3;
    // max number of layers needed in the Pascal tetrahedron
    cornerTriangle = 0;
    // number of elements subtracted in each corner Pascal triangle
    nTermsInLayer = d;
    // initializing number of elements in layer
    excess = 0;
    // excess based on overlapping of growing Pascal triangles
    i = 1 - degree;
    for (::coder::SizeType p{i}; p <= maxLayers; p++) {
      ::coder::SizeType counterBottomRow;
      ::coder::SizeType gap;
      ::coder::SizeType nTermsInPrevLayer;
      //  Within each level, x^deg is at the peak of Pascal triangle
      cornerTriangle = (cornerTriangle + p) + degree;
      counterBottomRow = 1;
      // counter for the bottom row to be subtracted later
      for (::coder::SizeType k{0}; k < deg; k++) {
        for (::coder::SizeType iPnt{0}; iPnt < npoints; iPnt++) {
          V[iPnt + V.size(1) * c] = V[iPnt + V.size(1) * (c - nTermsInLayer)] *
                                    us[us.size(1) * iPnt + 1];
        }
        c++;
        counterBottomRow++;
      }
      deg--;
      b_degree = ((degree + degree) + p) - 1;
      if (b_degree < 0) {
        b_degree = 0;
      }
      excess += b_degree;
      d = (d + p) + 1;
      // number of terms in Pascal tetrahedron
      nTermsInPrevLayer = nTermsInLayer;
      nTermsInLayer = d + 3 * (excess - cornerTriangle);
      gap = (nTermsInPrevLayer + counterBottomRow) - 1;
      b_degree = nTermsInLayer - counterBottomRow;
      for (::coder::SizeType j{0}; j <= b_degree; j++) {
        for (::coder::SizeType iPnt{0}; iPnt < npoints; iPnt++) {
          V[iPnt + V.size(1) * c] =
              V[iPnt + V.size(1) * (c - gap)] * us[us.size(1) * iPnt + 2];
        }
        c++;
      }
    }
  }
  m2cAssert(true, "");
}

//  gen_vander_3d  Generate generalized/confluent Vandermonde matrix in 3D.
static void gen_vander_3d(const ::coder::array<real_T, 2U> &us, ::coder::SizeType npoints,
                          ::coder::SizeType degree, ::coder::array<real_T, 2U> &V)
{
  ::coder::SizeType b_degree;
  ::coder::SizeType c;
  ::coder::SizeType d;
  ::coder::SizeType deg;
  ::coder::SizeType i;
  if (npoints == 0) {
    npoints = us.size(0);
  } else if (npoints > us.size(0)) {
    m2cErrMsgIdAndTxt("wlslib:BufferTooSmall", "Input us is too small.");
  }
  //  Allocate storage for V
  if (degree >= 0) {
    b_degree = (degree + 1) * (degree + 2) * (degree + 3) / 6;
  } else {
    b_degree = (1 - degree) * (1 - degree) * (1 - degree);
  }
  V.set_size(b_degree, us.size(0));
  //  compute 0th order generalized Vandermonde matrix
  if (degree != 0) {
    for (::coder::SizeType iPnt{0}; iPnt < npoints; iPnt++) {
      V[iPnt] = 1.0;
      V[iPnt + V.size(1)] = us[us.size(1) * iPnt];
      V[iPnt + V.size(1) * 2] = us[us.size(1) * iPnt + 1];
      V[iPnt + V.size(1) * 3] = us[us.size(1) * iPnt + 2];
    }
  } else {
    for (::coder::SizeType iPnt{0}; iPnt < npoints; iPnt++) {
      V[iPnt] = 1.0;
    }
  }
  c = 4;
  d = 3;
  if (degree < 0) {
    i = -degree;
  } else {
    i = degree;
  }
  for (deg = 2; deg <= i; deg++) {
    //  Within each level, use convention of Pascal triangle with x^deg at peak
    for (::coder::SizeType j{0}; j < deg; j++) {
      for (::coder::SizeType iPnt{0}; iPnt < npoints; iPnt++) {
        V[iPnt + V.size(1) * c] =
            V[iPnt + V.size(1) * (c - d)] * us[us.size(1) * iPnt];
      }
      c++;
    }
    for (::coder::SizeType iPnt{0}; iPnt < npoints; iPnt++) {
      V[iPnt + V.size(1) * c] =
          V[iPnt + V.size(1) * ((c - d) - 1)] * us[us.size(1) * iPnt + 1];
    }
    c++;
    for (::coder::SizeType j{0}; j < d; j++) {
      for (::coder::SizeType iPnt{0}; iPnt < npoints; iPnt++) {
        V[iPnt + V.size(1) * c] = V[iPnt + V.size(1) * (((c - d) - deg) - 1)] *
                                  us[us.size(1) * iPnt + 2];
      }
      c++;
    }
    d = (d + deg) + 1;
  }
  //  Compute the tri-degree terms if degree<0
  if (degree < 0) {
    ::coder::SizeType cornerTriangle;
    ::coder::SizeType excess;
    ::coder::SizeType maxLayers;
    ::coder::SizeType nTermsInLayer;
    deg = -degree;
    maxLayers = -degree * 3;
    // max number of layers needed in the Pascal tetrahedron
    cornerTriangle = 0;
    // number of elements subtracted in each corner Pascal triangle
    nTermsInLayer = d;
    // initializing number of elements in layer
    excess = 0;
    // excess based on overlapping of growing Pascal triangles
    i = 1 - degree;
    for (::coder::SizeType p{i}; p <= maxLayers; p++) {
      ::coder::SizeType counterBottomRow;
      ::coder::SizeType gap;
      ::coder::SizeType nTermsInPrevLayer;
      //  Within each level, x^deg is at the peak of Pascal triangle
      cornerTriangle = (cornerTriangle + p) + degree;
      counterBottomRow = 1;
      // counter for the bottom row to be subtracted later
      for (::coder::SizeType k{0}; k < deg; k++) {
        for (::coder::SizeType iPnt{0}; iPnt < npoints; iPnt++) {
          V[iPnt + V.size(1) * c] = V[iPnt + V.size(1) * (c - nTermsInLayer)] *
                                    us[us.size(1) * iPnt + 1];
        }
        c++;
        counterBottomRow++;
      }
      deg--;
      b_degree = ((degree + degree) + p) - 1;
      if (b_degree < 0) {
        b_degree = 0;
      }
      excess += b_degree;
      d = (d + p) + 1;
      // number of terms in Pascal tetrahedron
      nTermsInPrevLayer = nTermsInLayer;
      nTermsInLayer = d + 3 * (excess - cornerTriangle);
      gap = (nTermsInPrevLayer + counterBottomRow) - 1;
      b_degree = nTermsInLayer - counterBottomRow;
      for (::coder::SizeType j{0}; j <= b_degree; j++) {
        for (::coder::SizeType iPnt{0}; iPnt < npoints; iPnt++) {
          V[iPnt + V.size(1) * c] =
              V[iPnt + V.size(1) * (c - gap)] * us[us.size(1) * iPnt + 2];
        }
        c++;
      }
    }
  }
  m2cAssert(true, "");
}

// hexa_125 - Triquartic hexahedral element with equidistant points
static inline
void hexa_125(real_T xi, real_T eta, real_T zeta, real_T sfvals[125],
                     real_T sdvals[375])
{
  ::sfe_sfuncs::hexa_125_sfunc(xi, eta, zeta, &sfvals[0], &sdvals[0]);
}

// hexa_216 - Triquintic hexahedral element with equidistant points
static inline
void hexa_216(real_T xi, real_T eta, real_T zeta, real_T sfvals[216],
                     real_T sdvals[648])
{
  ::sfe_sfuncs::hexa_216_sfunc(xi, eta, zeta, &sfvals[0], &sdvals[0]);
}

// hexa_343 - Trisextic hexahedral element with equidistant points
static inline
void hexa_343(real_T xi, real_T eta, real_T zeta, real_T sfvals[343],
                     real_T sdvals[1029])
{
  ::sfe_sfuncs::hexa_343_sfunc(xi, eta, zeta, &sfvals[0], &sdvals[0]);
}

// hexa_64 - Tricubic hexahedral element with equidistant nodes
static inline
void hexa_64(real_T xi, real_T eta, real_T zeta, real_T sfvals[64],
                    real_T sdvals[192])
{
  ::sfe_sfuncs::hexa_64_sfunc(xi, eta, zeta, &sfvals[0], &sdvals[0]);
}

// hexa_gl_125 - Triquartic hexahedral element with Gauss-Lobatto points
static inline
void hexa_gl_125(real_T xi, real_T eta, real_T zeta, real_T sfvals[125],
                        real_T sdvals[375])
{
  ::sfe_sfuncs::hexa_gl_125_sfunc(xi, eta, zeta, &sfvals[0], &sdvals[0]);
}

// hexa_gl_216 - Triquintic hexahedral element with equidistant points
static inline
void hexa_gl_216(real_T xi, real_T eta, real_T zeta, real_T sfvals[216],
                        real_T sdvals[648])
{
  ::sfe_sfuncs::hexa_gl_216_sfunc(xi, eta, zeta, &sfvals[0], &sdvals[0]);
}

// hexa_gl_343 - Trisextic hexahedral element with equidistant points
static inline
void hexa_gl_343(real_T xi, real_T eta, real_T zeta, real_T sfvals[343],
                        real_T sdvals[1029])
{
  ::sfe_sfuncs::hexa_gl_343_sfunc(xi, eta, zeta, &sfvals[0], &sdvals[0]);
}

// hexa_gl_64 - Tricubic hexahedral element with Gauss-Lobatto nodes
static inline
void hexa_gl_64(real_T xi, real_T eta, real_T zeta, real_T sfvals[64],
                       real_T sdvals[192])
{
  ::sfe_sfuncs::hexa_gl_64_sfunc(xi, eta, zeta, &sfvals[0], &sdvals[0]);
}

// init_osusop - Initialize CRS for OSUS operator
static void init_osusop(const ::coder::array<real_T, 2U> &mesh_coords,
                        const ::coder::array<ConnData, 1U> &mesh_elemtables,
                        const ::coder::array<uint64_T, 1U> &mesh_teids,
                        const ::coder::array<Stencils, 1U> &mesh_stencils,
                        ::coder::SizeType stclid, ::coder::array<int64_T, 1U> &A_row_ptr,
                        ::coder::array<int32_T, 1U> &A_col_ind,
                        ::coder::array<real_T, 1U> &A_val, int32_T *A_ncols,
                        ::coder::array<int32_T, 1U> &nnzs)
{
  ::coder::array<int32_T, 1U> nodes_;
  ::coder::array<boolean_T, 1U> visited_;
  ::coder::SizeType i;
  ::coder::SizeType i1;
  ::coder::SizeType j;
  ::coder::SizeType loop_ub;
  ::coder::SizeType nid;
  ::coder::SizeType npe;
  A_row_ptr.set_size(mesh_teids.size(0) + 1);
  A_row_ptr[0] = 1L;
  *A_ncols = mesh_coords.size(0);
  //  Allocate for nnz per row
  nnzs.set_size(mesh_teids.size(0));
  visited_.set_size(mesh_coords.size(0));
  loop_ub = mesh_coords.size(0);
  for (i = 0; i < loop_ub; i++) {
    visited_[i] = false;
  }
  nodes_.set_size(mesh_coords.size(0));
  i = mesh_elemtables.size(0) - 1;
  for (::coder::SizeType etable{0}; etable <= i; etable++) {
    npe = mesh_elemtables[etable].conn.size(1);
    i1 = mesh_elemtables[etable].conn.size(0);
    for (::coder::SizeType e{0}; e < i1; e++) {
      ::coder::SizeType eid;
      //  Element ID
      eid = mesh_elemtables[etable].istart + e;
      j = 0;
      for (::coder::SizeType b_i{0}; b_i < npe; b_i++) {
        nid = mesh_elemtables[etable]
                  .conn[b_i + mesh_elemtables[etable].conn.size(1) * e];
        //  Stencil associated with nid
        for (int64_T k = mesh_stencils[stclid - 1].ngbverts.row_ptr[nid - 1];
             k < mesh_stencils[stclid - 1].ngbverts.row_ptr[nid]; k++) {
          //  Check if already visited
          loop_ub = mesh_stencils[stclid - 1]
                        .ngbverts.col_ind[k - 1];
          if (!visited_[loop_ub - 1]) {
            j++;
            nodes_[j - 1] = loop_ub;
            visited_[loop_ub - 1] = true;
          }
        }
      }
      //  Reset visited nodes
      for (::coder::SizeType b_i{0}; b_i < j; b_i++) {
        visited_[nodes_[b_i] - 1] = false;
      }
      //  Actual nnz
      nnzs[eid - 1] = j;
      //  Nnz with extra reserved space
      A_row_ptr[eid] =
          A_row_ptr[eid - 1] +
          static_cast<int64_T>(std::round(static_cast<real_T>(j) * 1.2));
    }
  }
  //  Allocate col_ind and val
  A_col_ind.set_size(
      (A_row_ptr[A_row_ptr.size(0) - 1] - 1L));
  A_val.set_size((A_row_ptr[A_row_ptr.size(0) - 1] - 1L));
  loop_ub = (A_row_ptr[A_row_ptr.size(0) - 1] - 1L);
  for (i1 = 0; i1 < loop_ub; i1++) {
    A_val[i1] = 0.0;
  }
  //  Loop again to put actual IDs to A
  for (::coder::SizeType etable{0}; etable <= i; etable++) {
    npe = mesh_elemtables[etable].conn.size(1);
    i1 = mesh_elemtables[etable].conn.size(0);
    for (::coder::SizeType e{0}; e < i1; e++) {
      int64_T istart;
      //  Element ID
      j = -1;
      for (::coder::SizeType b_i{0}; b_i < npe; b_i++) {
        nid = mesh_elemtables[etable]
                  .conn[b_i + mesh_elemtables[etable].conn.size(1) * e];
        //  Stencil associated with nid
        for (int64_T k = mesh_stencils[stclid - 1].ngbverts.row_ptr[nid - 1];
             k < mesh_stencils[stclid - 1].ngbverts.row_ptr[nid]; k++) {
          //  Check if already visited
          loop_ub = mesh_stencils[stclid - 1]
                        .ngbverts.col_ind[k - 1];
          if (!visited_[loop_ub - 1]) {
            j++;
            nodes_[j] = loop_ub;
            visited_[loop_ub - 1] = true;
          }
        }
      }
      //  Put to A and reset visited nodes
      istart = A_row_ptr[(mesh_elemtables[etable].istart + e) - 1];
      for (::coder::SizeType b_i{0}; b_i <= j; b_i++) {
        visited_[nodes_[b_i] - 1] = false;
        A_col_ind[static_cast<::coder::SizeType>((istart + (b_i + 1)) - 1L) - 1] =
            nodes_[b_i];
      }
    }
  }
}

static void insert_mem_crs(::coder::SizeType i, ::coder::array<int64_T, 1U> &row_ptr,
                           ::coder::array<int32_T, 1U> &col_ind,
                           ::coder::array<real_T, 1U> &val)
{
  ::coder::array<real_T, 1U> b_val;
  ::coder::array<int32_T, 1U> b_col_ind;
  int64_T dif;
  int64_T nnzcur;
  int64_T x;
  ::coder::SizeType b_i;
  ::coder::SizeType b_loop_ub;
  ::coder::SizeType i1;
  ::coder::SizeType loop_ub;
  ::coder::SizeType n;
  // Local buffers with patterns "\w+" are legitimate
  nnzcur = row_ptr[i] - row_ptr[i - 1];
  //  20% more
  x = static_cast<int64_T>(std::round(static_cast<real_T>(nnzcur) * 1.2));
  if (x < nnzcur + 1L) {
    x = nnzcur + 1L;
  }
  dif = x - nnzcur;
  if (row_ptr[i] - 1L < 1L) {
    loop_ub = 0;
  } else {
    loop_ub = (row_ptr[i] - 1L);
  }
  if (row_ptr[i] > col_ind.size(0)) {
    b_i = 0;
    i1 = 0;
  } else {
    b_i = (row_ptr[i]) - 1;
    i1 = col_ind.size(0);
  }
  b_col_ind.set_size(((loop_ub + dif) + i1) - b_i);
  for (::coder::SizeType i2{0}; i2 < loop_ub; i2++) {
    b_col_ind[i2] = col_ind[i2];
  }
  b_loop_ub = dif;
  for (::coder::SizeType i2{0}; i2 < b_loop_ub; i2++) {
    b_col_ind[i2 + loop_ub] = 0;
  }
  b_loop_ub = i1 - b_i;
  for (i1 = 0; i1 < b_loop_ub; i1++) {
    b_col_ind[(i1 + loop_ub) + dif] = col_ind[b_i + i1];
  }
  col_ind.set_size(b_col_ind.size(0));
  loop_ub = b_col_ind.size(0);
  for (b_i = 0; b_i < loop_ub; b_i++) {
    col_ind[b_i] = b_col_ind[b_i];
  }
  if (row_ptr[i] - 1L < 1L) {
    loop_ub = 0;
  } else {
    loop_ub = (row_ptr[i] - 1L);
  }
  if (row_ptr[i] > val.size(0)) {
    b_i = 0;
    i1 = 0;
  } else {
    b_i = (row_ptr[i]) - 1;
    i1 = val.size(0);
  }
  b_val.set_size(((loop_ub + dif) + i1) - b_i);
  for (::coder::SizeType i2{0}; i2 < loop_ub; i2++) {
    b_val[i2] = val[i2];
  }
  b_loop_ub = dif;
  for (::coder::SizeType i2{0}; i2 < b_loop_ub; i2++) {
    b_val[i2 + loop_ub] = 0.0;
  }
  b_loop_ub = i1 - b_i;
  for (i1 = 0; i1 < b_loop_ub; i1++) {
    b_val[(i1 + loop_ub) + dif] = val[b_i + i1];
  }
  val.set_size(b_val.size(0));
  loop_ub = b_val.size(0);
  for (b_i = 0; b_i < loop_ub; b_i++) {
    val[b_i] = b_val[b_i];
  }
  n = row_ptr.size(0);
  b_i = i + 1;
  for (::coder::SizeType j{b_i}; j <= n; j++) {
    row_ptr[j - 1] = row_ptr[j - 1] + dif;
  }
}

// m2cFind - Search the position (index) of a `key` in a given `keys`
static int64_T m2cFind(const ::coder::array<int32_T, 1U> &keys, ::coder::SizeType key,
                       int64_T b_first, int64_T last)
{
  int64_T i;
  int64_T k;
  k = 0L;
  //  Perform linear search
  i = b_first;
  while (i <= last) {
    if (keys[i - 1] == key) {
      k = i;
      i = last + 1L;
    } else {
      i++;
    }
  }
  return k;
}

//  m2cSort - Sort keys in an array using std::sort in C++ STL.
static inline
void m2cSort(::coder::array<int32_T, 1U> &keys, ::coder::SizeType b_first,
                    ::coder::SizeType last)
{
  if (last > b_first) {
    auto last_ptr = 1 + (&keys[last - 1]);
    std::sort(&keys[b_first - 1], last_ptr);
  }
}

static void majorityTransform(const b_WlsObject *r, c_WlsObject *r1)
{
  ::coder::SizeType b_loop_ub;
  ::coder::SizeType loop_ub;
  r1->runtimes.size[0] = r->runtimes.size[0];
  loop_ub = r->runtimes.size[0];
  if (loop_ub - 1 >= 0) {
    std::copy(&r->runtimes.data[0], &r->runtimes.data[loop_ub],
              &r1->runtimes.data[0]);
  }
  r1->QRt.set_size(r->QRt.size(1), r->QRt.size(0));
  loop_ub = r->QRt.size(1);
  for (::coder::SizeType i{0}; i < loop_ub; i++) {
    b_loop_ub = r->QRt.size(0);
    for (::coder::SizeType i1{0}; i1 < b_loop_ub; i1++) {
      r1->QRt[i1 + r1->QRt.size(1) * i] = r->QRt[i + r->QRt.size(1) * i1];
    }
  }
  r1->rowmajor = r->rowmajor;
  r1->work.set_size(r->work.size(0));
  loop_ub = r->work.size(0);
  for (::coder::SizeType i{0}; i < loop_ub; i++) {
    r1->work[i] = r->work[i];
  }
  r1->jpvt.set_size(r->jpvt.size(0));
  loop_ub = r->jpvt.size(0);
  for (::coder::SizeType i{0}; i < loop_ub; i++) {
    r1->jpvt[i] = r->jpvt[i];
  }
  r1->fullrank = r->fullrank;
  r1->rank = r->rank;
  r1->ncols = r->ncols;
  r1->nrows = r->nrows;
  r1->nevpnts = r->nevpnts;
  r1->rhs.set_size(r->rhs.size(1), r->rhs.size(0));
  loop_ub = r->rhs.size(1);
  for (::coder::SizeType i{0}; i < loop_ub; i++) {
    b_loop_ub = r->rhs.size(0);
    for (::coder::SizeType i1{0}; i1 < b_loop_ub; i1++) {
      r1->rhs[i1 + r1->rhs.size(1) * i] = r->rhs[i + r->rhs.size(1) * i1];
    }
  }
  r1->QR.set_size(r->QR.size(1), r->QR.size(0));
  loop_ub = r->QR.size(1);
  for (::coder::SizeType i{0}; i < loop_ub; i++) {
    b_loop_ub = r->QR.size(0);
    for (::coder::SizeType i1{0}; i1 < b_loop_ub; i1++) {
      r1->QR[i1 + r1->QR.size(1) * i] = r->QR[i + r->QR.size(1) * i1];
    }
  }
  r1->V.set_size(r->V.size(1), r->V.size(0));
  loop_ub = r->V.size(1);
  for (::coder::SizeType i{0}; i < loop_ub; i++) {
    b_loop_ub = r->V.size(0);
    for (::coder::SizeType i1{0}; i1 < b_loop_ub; i1++) {
      r1->V[i1 + r1->V.size(1) * i] = r->V[i + r->V.size(1) * i1];
    }
  }
  r1->hs_inv.size[1] = 1;
  r1->hs_inv.size[0] = r->hs_inv.size[1];
  loop_ub = r->hs_inv.size[1];
  if (loop_ub - 1 >= 0) {
    std::copy(&r->hs_inv.data[0], &r->hs_inv.data[loop_ub],
              &r1->hs_inv.data[0]);
  }
  r1->rweights.set_size(r->rweights.size(0));
  loop_ub = r->rweights.size(0);
  for (::coder::SizeType i{0}; i < loop_ub; i++) {
    r1->rweights[i] = r->rweights[i];
  }
  r1->origin.size[1] = 1;
  r1->origin.size[0] = r->origin.size[1];
  loop_ub = r->origin.size[1];
  if (loop_ub - 1 >= 0) {
    std::copy(&r->origin.data[0], &r->origin.data[loop_ub],
              &r1->origin.data[0]);
  }
  r1->us.set_size(r->us.size(1), r->us.size(0));
  loop_ub = r->us.size(1);
  for (::coder::SizeType i{0}; i < loop_ub; i++) {
    b_loop_ub = r->us.size(0);
    for (::coder::SizeType i1{0}; i1 < b_loop_ub; i1++) {
      r1->us[i1 + r1->us.size(1) * i] = r->us[i + r->us.size(1) * i1];
    }
  }
  r1->stride = r->stride;
  r1->interp0 = r->interp0;
  r1->unimono = r->unimono;
  r1->order = r->order;
  r1->degree = r->degree;
  r1->nstpnts = r->nstpnts;
}

static void majorityTransform(const d_WlsObject *r, e_WlsObject *r1)
{
  ::coder::SizeType b_loop_ub;
  ::coder::SizeType loop_ub;
  r1->runtimes.size[0] = r->runtimes.size[0];
  loop_ub = r->runtimes.size[0];
  if (loop_ub - 1 >= 0) {
    std::copy(&r->runtimes.data[0], &r->runtimes.data[loop_ub],
              &r1->runtimes.data[0]);
  }
  r1->QRt.set_size(r->QRt.size(1), r->QRt.size(0));
  loop_ub = r->QRt.size(1);
  for (::coder::SizeType i{0}; i < loop_ub; i++) {
    b_loop_ub = r->QRt.size(0);
    for (::coder::SizeType i1{0}; i1 < b_loop_ub; i1++) {
      r1->QRt[i1 + r1->QRt.size(1) * i] = r->QRt[i + r->QRt.size(1) * i1];
    }
  }
  r1->rowmajor = r->rowmajor;
  r1->work.set_size(r->work.size(0));
  loop_ub = r->work.size(0);
  for (::coder::SizeType i{0}; i < loop_ub; i++) {
    r1->work[i] = r->work[i];
  }
  r1->jpvt.set_size(r->jpvt.size(0));
  loop_ub = r->jpvt.size(0);
  for (::coder::SizeType i{0}; i < loop_ub; i++) {
    r1->jpvt[i] = r->jpvt[i];
  }
  r1->fullrank = r->fullrank;
  r1->rank = r->rank;
  r1->ncols = r->ncols;
  r1->nrows = r->nrows;
  r1->nevpnts = r->nevpnts;
  r1->rhs.set_size(r->rhs.size(1), r->rhs.size(0));
  loop_ub = r->rhs.size(1);
  for (::coder::SizeType i{0}; i < loop_ub; i++) {
    b_loop_ub = r->rhs.size(0);
    for (::coder::SizeType i1{0}; i1 < b_loop_ub; i1++) {
      r1->rhs[i1 + r1->rhs.size(1) * i] = r->rhs[i + r->rhs.size(1) * i1];
    }
  }
  r1->QR.set_size(r->QR.size(1), r->QR.size(0));
  loop_ub = r->QR.size(1);
  for (::coder::SizeType i{0}; i < loop_ub; i++) {
    b_loop_ub = r->QR.size(0);
    for (::coder::SizeType i1{0}; i1 < b_loop_ub; i1++) {
      r1->QR[i1 + r1->QR.size(1) * i] = r->QR[i + r->QR.size(1) * i1];
    }
  }
  r1->V.set_size(r->V.size(1), r->V.size(0));
  loop_ub = r->V.size(1);
  for (::coder::SizeType i{0}; i < loop_ub; i++) {
    b_loop_ub = r->V.size(0);
    for (::coder::SizeType i1{0}; i1 < b_loop_ub; i1++) {
      r1->V[i1 + r1->V.size(1) * i] = r->V[i + r->V.size(1) * i1];
    }
  }
  r1->hs_inv.set_size(1, r->hs_inv.size(0));
  loop_ub = r->hs_inv.size(0);
  for (::coder::SizeType i{0}; i < loop_ub; i++) {
    r1->hs_inv[i] = r->hs_inv[i];
  }
  r1->rweights.set_size(r->rweights.size(0));
  loop_ub = r->rweights.size(0);
  for (::coder::SizeType i{0}; i < loop_ub; i++) {
    r1->rweights[i] = r->rweights[i];
  }
  r1->origin.size[1] = r->origin.size[0];
  r1->origin.size[0] = 1;
  loop_ub = r->origin.size[0];
  if (loop_ub - 1 >= 0) {
    std::copy(&r->origin.data[0], &r->origin.data[loop_ub],
              &r1->origin.data[0]);
  }
  r1->us.set_size(r->us.size(1), r->us.size(0));
  loop_ub = r->us.size(1);
  for (::coder::SizeType i{0}; i < loop_ub; i++) {
    b_loop_ub = r->us.size(0);
    for (::coder::SizeType i1{0}; i1 < b_loop_ub; i1++) {
      r1->us[i1 + r1->us.size(1) * i] = r->us[i + r->us.size(1) * i1];
    }
  }
  r1->stride = r->stride;
  r1->interp0 = r->interp0;
  r1->unimono = r->unimono;
  r1->order = r->order;
  r1->degree = r->degree;
  r1->nstpnts = r->nstpnts;
}

static void mark_discontinuities_kernel(
    const ::coder::array<real_T, 2U> &mesh_coords,
    const ::coder::array<ConnData, 1U> &mesh_elemtables,
    const ::coder::array<uint64_T, 1U> &mesh_teids,
    const ::coder::array<int64_T, 1U> &mesh_node2elems_row_ptr,
    const ::coder::array<int32_T, 1U> &mesh_node2elems_col_ind,
    const ::coder::array<real_T, 1U> &mesh_elemh, real_T mesh_globalh,
    const ::coder::array<real_T, 2U> &fs, const ::coder::array<real_T, 2U> &df,
    const ::coder::array<real_T, 2U> &alpha,
    const ::coder::array<real_T, 2U> &beta, real_T cglobal, real_T clocal,
    real_T kappa1, real_T kappa0, ::coder::array<int8_T, 2U> &distags)
{
  real_T thres;
  ::coder::SizeType iend;
  ::coder::SizeType istart;
  ::coder::SizeType nrhs;
  ::coder::SizeType nthreads;
  ::coder::SizeType u1;
  nrhs = alpha.size(1);
#pragma omp single
  { // single
    distags.set_size(mesh_coords.size(0), alpha.size(1));
  } // single
  nthreads = 1;
#ifdef _OPENMP
  nthreads = omp_get_num_threads();
#endif // _OPENMP
  if (nthreads == 1) {
    istart = 0;
    iend = mesh_coords.size(0);
  } else {
    ::coder::SizeType b_remainder;
    ::coder::SizeType chunk;
    ::coder::SizeType threadID;
    threadID = 0;
#ifdef _OPENMP
    threadID = omp_get_thread_num();
#endif // _OPENMP
    chunk = mesh_coords.size(0) / nthreads;
    b_remainder = mesh_coords.size(0) - nthreads * chunk;
    u1 = threadID;
    if (b_remainder <= threadID) {
      u1 = b_remainder;
    }
    istart = threadID * chunk + u1;
    iend = (istart + chunk) + (threadID < b_remainder);
  }
  thres = cglobal * std::pow(mesh_globalh, 1.5);
  for (::coder::SizeType i{istart + 1}; i <= iend; i++) {
    for (::coder::SizeType k{0}; k < nrhs; k++) {
      int64_T j;
      real_T tauglobal;
      boolean_T discell;
      boolean_T exitg1;
      distags[k + distags.size(1) * (i - 1)] = 0;
      tauglobal = thres * df[k];
      discell = false;
      j = mesh_node2elems_row_ptr[i - 1];
      exitg1 = false;
      while ((!exitg1) && (j <= mesh_node2elems_row_ptr[i] - 1L)) {
        real_T b_fmax;
        real_T b_fmin;
        uint64_T c;
        ::coder::SizeType leid;
        ::coder::SizeType u0;
        u0 = mesh_node2elems_col_ind[j - 1] - 1;
        c = mesh_teids[u0] & 255UL;
        leid = (
            mesh_teids[mesh_node2elems_col_ind[j - 1] -
                       1] >>
            8);
        b_fmax = -1.7976931348623157E+308;
        b_fmin = 1.7976931348623157E+308;
        u1 = mesh_elemtables[(
                                 mesh_teids[mesh_node2elems_col_ind
                                                [j - 1] -
                                            1] &
                                 255UL) -
                             1]
                 .conn.size(1);
        for (::coder::SizeType ii{0}; ii < u1; ii++) {
          real_T fvalue;
          fvalue =
              fs[k + fs.size(1) *
                         (mesh_elemtables[c - 1].conn
                              [ii + mesh_elemtables[c - 1]
                                            .conn.size(1) *
                                        (leid - 1)] -
                          1)];
          if (b_fmax < fvalue) {
            b_fmax = fvalue;
          }
          if (b_fmin > fvalue) {
            b_fmin = fvalue;
          }
        }
        //  Compute local df and thres
        if (std::abs(alpha[k + alpha.size(1) * u0]) >
            std::fmax(tauglobal,
                      clocal * (b_fmax - b_fmin) * std::sqrt(mesh_elemh[u0]))) {
          //  Dis cell
          discell = true;
          exitg1 = true;
        } else {
          j++;
        }
      }
      if (discell) {
        real_T d;
        d = beta[k + beta.size(1) * (i - 1)];
        if (d > kappa1) {
          distags[k + distags.size(1) * (i - 1)] = 2;
          //  C1
          if (d > kappa0) {
            distags[k + distags.size(1) * (i - 1)] = 1;
          }
          //  C0
        }
      }
    }
  }
}

// obtain_nring_1d - Collect n-ring vertices of a 1D mesh
static void obtain_nring_1d(const ::coder::array<ConnData, 1U> &mesh_elemtables,
                            const ::coder::array<uint64_T, 1U> &mesh_teids,
                            const ::coder::array<uint64_T, 2U> &mesh_sibhfs,
                            const ::coder::array<uint64_T, 1U> &mesh_v2hfid,
                            ::coder::SizeType vid, real_T ring, ::coder::SizeType maxnpnts,
                            ::coder::array<boolean_T, 1U> &vtags,
                            ::coder::array<boolean_T, 1U> &ftags,
                            ::coder::array<int32_T, 1U> &ngbvs, int32_T *nverts,
                            ::coder::array<int32_T, 1U> &ngbfs, int32_T *nfaces,
                            ::coder::array<uint64_T, 1U> &hebuf,
                            boolean_T *reflected, boolean_T *overflow)
{
  uint64_T c_tmp;
  uint64_T fid;
  uint64_T oppfid;
  boolean_T oneringonly;
  fid = mesh_v2hfid[vid - 1] >> 8;
  c_tmp = mesh_v2hfid[vid - 1] & 255UL;
  *nfaces = 0;
  *overflow = false;
  //  Create output and buffers
  ngbvs.set_size(maxnpnts);
  ngbfs.set_size(maxnpnts);
  hebuf.set_size(maxnpnts);
  oneringonly = ring == 1.0;
  //  One ring comes from elements on both sides of vertex
  if (mesh_sibhfs[c_tmp +
                  mesh_sibhfs.size(1) * (fid - 1)] !=
      0UL) {
    uint64_T c;
    *nverts = 2;
    *nfaces = 2;
    ngbvs[0] =
        mesh_elemtables[(
                            mesh_teids[fid - 1] & 255UL) -
                        1]
            .conn[(mesh_elemtables[(
                                       mesh_teids[fid -
                                                  1] &
                                       255UL) -
                                   1]
                           .conn.size(1) *
                       ((
                            mesh_teids[fid - 1] >> 8) -
                        1) -
                   c_tmp) +
                  1];
    ngbfs[0] = fid;
    oppfid =
        mesh_sibhfs[c_tmp +
                    mesh_sibhfs.size(1) * (fid - 1)] >>
        8;
    c = mesh_sibhfs[c_tmp +
                    mesh_sibhfs.size(1) * (fid - 1)] &
        255UL;
    ngbvs[1] =
        mesh_elemtables[(
                            mesh_teids[oppfid - 1] &
                            255UL) -
                        1]
            .conn[(mesh_elemtables[(
                                       mesh_teids[oppfid -
                                                  1] &
                                       255UL) -
                                   1]
                           .conn.size(1) *
                       ((
                            mesh_teids[oppfid - 1] >> 8) -
                        1) -
                   c) +
                  1];
    ngbfs[1] = oppfid;
    if (!oneringonly) {
      hebuf[0] = ((fid << 8) + (2 - c_tmp)) - 1UL;
      hebuf[1] = ((oppfid << 8) + (2 - c)) - 1UL;
    }
  } else {
    //  If vertex is border edge, insert its incident border vertex.
    *nverts = 1;
    ngbvs[0] =
        mesh_elemtables[(
                            mesh_teids[fid - 1] & 255UL) -
                        1]
            .conn[(mesh_elemtables[(
                                       mesh_teids[fid -
                                                  1] &
                                       255UL) -
                                   1]
                           .conn.size(1) *
                       ((
                            mesh_teids[fid - 1] >> 8) -
                        1) -
                   c_tmp) +
                  1];
    if (!oneringonly) {
      hebuf[0] = ((fid << 8) + (2 - c_tmp)) - 1UL;
      ring *= 2.0;
    }
  }
  *reflected = false;
  if (ring != 1.0) {
    real_T cur_ring;
    real_T orig_ring;
    ::coder::SizeType nverts_pre;
    vtags[vid - 1] = true;
    for (::coder::SizeType i{0}; i < *nverts; i++) {
      vtags[ngbvs[i] - 1] = true;
    }
    for (::coder::SizeType i{0}; i < *nfaces; i++) {
      ftags[ngbfs[i] - 1] = true;
    }
    //  Define buffers and prepare tags for further processing
    nverts_pre = 0;
    orig_ring = ring;
    //  store original ring size
    cur_ring = 1.0;
    ::coder::SizeType b_i;
    ::coder::SizeType exitg1;
    ::coder::SizeType nfaces_pre;
    ::coder::SizeType nverts_last;
    do {
      exitg1 = 0;
      //  Collect next level of ring
      nverts_last = *nverts;
      nfaces_pre = *nfaces;
      b_i = nverts_pre + 1;
      for (::coder::SizeType ii{b_i}; ii <= nverts_last; ii++) {
        c_tmp = hebuf[ii - 1];
        if (c_tmp != 0UL) {
          ::coder::SizeType lid_tmp;
          fid = hebuf[ii - 1] >> 8;
          lid_tmp = (c_tmp & 255UL);
          oppfid = mesh_sibhfs[lid_tmp + mesh_sibhfs.size(1) *
                                             (fid - 1)] >>
                   8;
          if ((*overflow) ||
              ((!vtags[ngbvs[ii - 1] - 1]) && (*nverts >= maxnpnts))) {
            *overflow = true;
          } else {
            *overflow = false;
          }
          if ((oppfid != 0UL) && (!*overflow)) {
            ::coder::SizeType vopp;
            lid_tmp = static_cast<::coder::SizeType>(
                mesh_sibhfs[lid_tmp + mesh_sibhfs.size(1) *
                                          (fid - 1)] &
                255UL);
            vopp =
                mesh_elemtables[(
                                    mesh_teids[oppfid -
                                               1] &
                                    255UL) -
                                1]
                    .conn[(mesh_elemtables[static_cast<::coder::SizeType>(
                                               mesh_teids[(
                                                              oppfid) -
                                                          1] &
                                               255UL) -
                                           1]
                                   .conn.size(1) *
                               ((
                                    mesh_teids[oppfid -
                                               1] >>
                                    8) -
                                1) -
                           lid_tmp) +
                          1];
            (*nverts)++;
            (*nfaces)++;
            ngbvs[*nverts - 1] = vopp;
            ngbfs[*nfaces - 1] = oppfid;
            hebuf[*nverts - 1] = ((oppfid << 8) + (2 - lid_tmp)) - 1UL;
            vtags[vopp - 1] = true;
          } else if ((oppfid == 0UL) && (!*overflow)) {
            *reflected = true;
            ring = orig_ring + (orig_ring - cur_ring);
          }
        }
      }
      cur_ring++;
      if ((cur_ring >= ring) || (*nfaces == nfaces_pre) || (*overflow)) {
        exitg1 = 1;
      } else {
        nverts_pre = nverts_last;
      }
    } while (exitg1 == 0);
    //  Reset flags
    vtags[vid - 1] = false;
    for (::coder::SizeType i{0}; i < *nverts; i++) {
      vtags[ngbvs[i] - 1] = false;
    }
    for (::coder::SizeType i{0}; i < *nfaces; i++) {
      ftags[ngbfs[i] - 1] = false;
    }
  }
}

// obtain_nring_2d - Collect n-ring vertices of a 2D mesh
static void obtain_nring_2d(const ::coder::array<ConnData, 1U> &mesh_elemtables,
                            const ::coder::array<uint64_T, 1U> &mesh_teids,
                            const ::coder::array<uint64_T, 2U> &mesh_sibhfs,
                            const ::coder::array<uint64_T, 1U> &mesh_v2hfid,
                            ::coder::SizeType vid, real_T ring, ::coder::SizeType maxnpnts,
                            ::coder::array<boolean_T, 1U> &vtags,
                            ::coder::array<boolean_T, 1U> &ftags,
                            const ::coder::array<int32_T, 1U> &bridges,
                            ::coder::array<int32_T, 1U> &ngbvs, int32_T *nverts,
                            ::coder::array<int32_T, 1U> &ngbfs, int32_T *nfaces,
                            ::coder::array<uint64_T, 1U> &hebuf,
                            boolean_T *reflected, boolean_T *overflow)
{
  static const int8_T nxt[8]{2, 3, 1, 0, 2, 3, 4, 1};
  static const int8_T prv[8]{3, 1, 2, 0, 4, 1, 2, 3};
  uint64_T c_tmp;
  uint64_T fid;
  uint64_T fid_in;
  uint64_T opp;
  ::coder::SizeType exitg1;
  ::coder::SizeType i;
  ::coder::SizeType lid;
  ::coder::SizeType lid_prv;
  ::coder::SizeType maxnf;
  ::coder::SizeType ngbvs_tmp;
  boolean_T oneringonly;
  maxnf = maxnpnts << 1;
  //  element space
  fid = mesh_v2hfid[vid - 1] >> 8;
  c_tmp = mesh_v2hfid[vid - 1] & 255UL;
  lid = c_tmp;
  *nverts = 0;
  *nfaces = 0;
  *overflow = false;
  //  Create output and buffers
  ngbvs.set_size(maxnpnts);
  ngbfs.set_size(maxnf);
  hebuf.set_size(maxnpnts);
  oneringonly = ring == 1.0;
  if (mesh_sibhfs[c_tmp +
                  mesh_sibhfs.size(1) * (fid - 1)] !=
      0UL) {
    fid_in = fid;
  } else {
    fid_in = 0UL;
    //  Get shape of the element
    *nverts = 1;
    ngbvs_tmp = (mesh_teids[fid - 1] &
                                     255UL) -
                1;
    ngbvs[0] =
        mesh_elemtables[ngbvs_tmp]
            .conn[(nxt[c_tmp +
                       ((1 - ((mesh_elemtables[ngbvs_tmp].etype >> 5 & 7) == 2))
                        << 2)] +
                   mesh_elemtables[ngbvs_tmp].conn.size(1) *
                       ((
                            mesh_teids[fid - 1] >> 8) -
                        1)) -
                  1];
    if (!oneringonly) {
      hebuf[0] = 0UL;
    }
  }
  //  Rotate counterclockwise order around vertex and insert vertices
  do {
    exitg1 = 0;
    //  Get shape of the element
    ngbvs_tmp = (mesh_teids[fid - 1] &
                                     255UL) -
                1;
    lid_prv =
        prv[lid +
            ((1 - ((mesh_elemtables[ngbvs_tmp].etype >> 5 & 7) == 2)) << 2)];
    if ((*nverts < maxnpnts) && (*nfaces < maxnf)) {
      (*nverts)++;
      ngbvs[*nverts - 1] =
          mesh_elemtables[ngbvs_tmp]
              .conn[(lid_prv +
                     mesh_elemtables[ngbvs_tmp].conn.size(1) *
                         ((
                              mesh_teids[fid - 1] >> 8) -
                          1)) -
                    1];
      (*nfaces)++;
      ngbfs[*nfaces - 1] = fid;
      if (!oneringonly) {
        //  Save a starting edge for newly inserted vertex to allow early
        hebuf[*nverts - 1] = ((fid << 8) + lid_prv) - 1UL;
      }
    } else {
      *overflow = true;
    }
    //  Opposite face ID
    opp = mesh_sibhfs[(lid_prv +
                       mesh_sibhfs.size(1) * (fid - 1)) -
                      1];
    fid = opp >> 8;
    if (fid == fid_in) {
      exitg1 = 1;
    } else {
      lid = (opp & 255UL);
    }
  } while (exitg1 == 0);
  //  Handling bridges
  if ((bridges.size(0) != 0) && (*nverts < maxnpnts)) {
    i = bridges[vid - 1];
    if (i > 0) {
      //  use tags
      vtags[vid - 1] = true;
      for (::coder::SizeType b_i{0}; b_i < *nverts; b_i++) {
        vtags[ngbvs[b_i] - 1] = true;
      }
      if (!vtags[i - 1]) {
        //  reset tags
        vtags[vid - 1] = false;
        for (::coder::SizeType b_i{0}; b_i < *nverts; b_i++) {
          vtags[ngbvs[b_i] - 1] = false;
        }
        (*nverts)++;
        ngbvs[*nverts - 1] = i;
      }
    }
  }
  *reflected = false;
  if (ring != 1.0) {
    real_T cur_ring;
    real_T orig_ring;
    real_T ring_full;
    ::coder::SizeType nfaces_pre;
    ::coder::SizeType nverts_pre;
    orig_ring = ring;
    //  store original ring size
    vtags[vid - 1] = true;
    for (::coder::SizeType b_i{0}; b_i < *nverts; b_i++) {
      vtags[ngbvs[b_i] - 1] = true;
    }
    for (::coder::SizeType b_i{0}; b_i < *nfaces; b_i++) {
      ftags[ngbfs[b_i] - 1] = true;
    }
    //  Define buffers and prepare tags for further processing
    nverts_pre = 0;
    nfaces_pre = 0;
    //  Second, build full-size ring
    ring_full = std::trunc(ring);
    cur_ring = 1.0;
    boolean_T guard1{false};
    do {
      ::coder::SizeType lid_tmp;
      ::coder::SizeType nverts_last;
      ::coder::SizeType v;
      exitg1 = 0;
      guard1 = false;
      if ((cur_ring > ring_full) ||
          ((cur_ring == ring_full) && (ring_full != ring))) {
        ::coder::SizeType nfaces_last;
        //  Collect halfring
        nfaces_last = *nfaces;
        nverts_last = *nverts;
        i = nfaces_pre + 1;
        for (::coder::SizeType ii{i}; ii <= nfaces_last; ii++) {
          ngbvs_tmp = ngbfs[ii - 1] - 1;
          lid_tmp = (mesh_teids[ngbvs_tmp] & 255UL) - 1;
          if ((mesh_elemtables[lid_tmp].etype >> 5 & 7) == 2) {
            uint64_T oppe;
            boolean_T b_guard1{false};
            //  take opposite vertex in opposite face of triangle
            oppe = mesh_sibhfs[mesh_sibhfs.size(1) * ngbvs_tmp];
            c_tmp = oppe >> 8;
            b_guard1 = false;
            if (oppe != 0UL) {
              boolean_T b;
              b = ftags[c_tmp - 1];
              if (!b) {
                ngbvs_tmp =
                    (
                        mesh_teids[c_tmp - 1] & 255UL) -
                    1;
                v = mesh_elemtables[ngbvs_tmp].conn
                        [(prv[(oppe & 255UL) +
                              ((1 - ((mesh_elemtables[ngbvs_tmp].etype >> 5 &
                                      7) == 2))
                               << 2)] +
                          mesh_elemtables[ngbvs_tmp].conn.size(1) *
                              ((
                                   mesh_teids[c_tmp -
                                              1] >>
                                   8) -
                               1)) -
                         1] -
                    1;
                if ((*overflow) || ((!vtags[v]) && (*nverts >= maxnpnts)) ||
                    (*nfaces >= maxnf)) {
                  *overflow = true;
                } else {
                  *overflow = false;
                }
                if (!*overflow) {
                  (*nfaces)++;
                  ngbfs[*nfaces - 1] = c_tmp;
                  ftags[c_tmp - 1] = true;
                }
                if ((!vtags[v]) && (!*overflow)) {
                  (*nverts)++;
                  ngbvs[*nverts - 1] = v + 1;
                  vtags[v] = true;
                }
              } else {
                b_guard1 = true;
              }
            } else {
              b_guard1 = true;
            }
            if (b_guard1) {
              if ((oppe == 0UL) && (orig_ring == ring)) {
                //  Reflect rings to enlarge the number of remaining rings
                *reflected = true;
                ring = orig_ring + (orig_ring - cur_ring);
              }
              oppe = mesh_sibhfs[mesh_sibhfs.size(1) * ngbvs_tmp + 1];
              c_tmp = oppe >> 8;
              if ((oppe != 0UL) && (!ftags[c_tmp - 1])) {
                ngbvs_tmp =
                    (
                        mesh_teids[c_tmp - 1] & 255UL) -
                    1;
                v = mesh_elemtables[ngbvs_tmp].conn
                        [(prv[(oppe & 255UL) +
                              ((1 - ((mesh_elemtables[ngbvs_tmp].etype >> 5 &
                                      7) == 2))
                               << 2)] +
                          mesh_elemtables[ngbvs_tmp].conn.size(1) *
                              ((
                                   mesh_teids[c_tmp -
                                              1] >>
                                   8) -
                               1)) -
                         1] -
                    1;
                if ((*overflow) || ((!vtags[v]) && (*nverts >= maxnpnts)) ||
                    ((!ftags[c_tmp - 1]) &&
                     (*nfaces >= maxnf))) {
                  *overflow = true;
                } else {
                  *overflow = false;
                }
                if ((!ftags[c_tmp - 1]) && (!*overflow)) {
                  (*nfaces)++;
                  ngbfs[*nfaces - 1] = c_tmp;
                  ftags[c_tmp - 1] = true;
                }
                if ((!vtags[v]) && (!*overflow)) {
                  (*nverts)++;
                  ngbvs[*nverts - 1] = v + 1;
                  vtags[v] = true;
                }
              } else {
                if ((oppe == 0UL) && (orig_ring == ring)) {
                  //  Reflect rings to enlarge the number of remaining rings
                  *reflected = true;
                  ring = orig_ring + (orig_ring - cur_ring);
                }
                oppe = mesh_sibhfs[mesh_sibhfs.size(1) * ngbvs_tmp + 2];
                c_tmp = oppe >> 8;
                if ((oppe != 0UL) &&
                    (!ftags[c_tmp - 1])) {
                  ngbvs_tmp =
                      (
                          mesh_teids[c_tmp - 1] & 255UL) -
                      1;
                  v = mesh_elemtables[ngbvs_tmp].conn
                          [(prv[(oppe & 255UL) +
                                ((1 - ((mesh_elemtables[ngbvs_tmp].etype >> 5 &
                                        7) == 2))
                                 << 2)] +
                            mesh_elemtables[ngbvs_tmp].conn.size(1) *
                                ((
                                     mesh_teids[c_tmp -
                                                1] >>
                                     8) -
                                 1)) -
                           1] -
                      1;
                  if ((*overflow) || ((!vtags[v]) && (*nverts >= maxnpnts)) ||
                      ((!ftags[c_tmp - 1]) &&
                       (*nfaces >= maxnf))) {
                    *overflow = true;
                  } else {
                    *overflow = false;
                  }
                  if ((!ftags[c_tmp - 1]) &&
                      (!*overflow)) {
                    (*nfaces)++;
                    ngbfs[*nfaces - 1] = c_tmp;
                    ftags[c_tmp - 1] = true;
                  }
                  if ((!vtags[v]) && (!*overflow)) {
                    (*nverts)++;
                    ngbvs[*nverts - 1] = v + 1;
                    vtags[v] = true;
                  }
                } else if ((oppe == 0UL) && (orig_ring == ring)) {
                  //  Reflect rings to enlarge the number of remaining rings
                  *reflected = true;
                  ring = orig_ring + (orig_ring - cur_ring);
                }
              }
            }
          } else {
            //  Insert missing vertices in quads
            v = mesh_elemtables[lid_tmp]
                    .conn[mesh_elemtables[lid_tmp].conn.size(1) *
                          ((mesh_teids[ngbfs[ii - 1] - 1] >>
                                                8) -
                           1)];
            if (!vtags[v - 1]) {
              if (*nverts >= maxnpnts) {
                *overflow = true;
              } else {
                (*nverts)++;
                ngbvs[*nverts - 1] = v;
                vtags[v - 1] = true;
              }
            } else {
              v = mesh_elemtables[lid_tmp]
                      .conn[mesh_elemtables[lid_tmp].conn.size(1) *
                                ((
                                     mesh_teids[ngbfs[ii - 1] - 1] >> 8) -
                                 1) +
                            1];
              if (!vtags[v - 1]) {
                if (*nverts >= maxnpnts) {
                  *overflow = true;
                } else {
                  (*nverts)++;
                  ngbvs[*nverts - 1] = v;
                  vtags[v - 1] = true;
                }
              } else {
                v = mesh_elemtables[lid_tmp]
                        .conn[mesh_elemtables[lid_tmp].conn.size(1) *
                                  ((
                                       mesh_teids[ngbfs[ii - 1] - 1] >> 8) -
                                   1) +
                              2];
                if (!vtags[v - 1]) {
                  if (*nverts >= maxnpnts) {
                    *overflow = true;
                  } else {
                    (*nverts)++;
                    ngbvs[*nverts - 1] = v;
                    vtags[v - 1] = true;
                  }
                } else {
                  v = mesh_elemtables[lid_tmp]
                          .conn[mesh_elemtables[lid_tmp].conn.size(1) *
                                    ((
                                         mesh_teids[ngbfs[ii - 1] - 1] >> 8) -
                                     1) +
                                3];
                  if (!vtags[v - 1]) {
                    if (*nverts >= maxnpnts) {
                      *overflow = true;
                    } else {
                      (*nverts)++;
                      ngbvs[*nverts - 1] = v;
                      vtags[v - 1] = true;
                    }
                  }
                }
              }
            }
          }
        }
        if ((cur_ring + 0.5 >= ring) || (*nverts >= maxnpnts) ||
            (*nfaces >= maxnf)) {
          exitg1 = 1;
        } else {
          //  If needs to expand, then undo the last half ring
          i = nverts_last + 1;
          for (::coder::SizeType b_i{i}; b_i <= *nverts; b_i++) {
            vtags[ngbvs[b_i - 1] - 1] = false;
          }
          *nverts = nverts_last;
          i = nfaces_last + 1;
          for (::coder::SizeType b_i{i}; b_i <= *nfaces; b_i++) {
            ftags[ngbfs[b_i - 1] - 1] = false;
          }
          *nfaces = nfaces_last;
          guard1 = true;
        }
      } else {
        guard1 = true;
      }
      if (guard1) {
        //  Collect next full level of ring
        nverts_last = *nverts;
        nfaces_pre = *nfaces;
        i = nverts_pre + 1;
        for (::coder::SizeType ii{i}; ii <= nverts_last; ii++) {
          ::coder::SizeType pos;
          boolean_T allow_early_term;
          boolean_T isfirst;
          fid = mesh_v2hfid[ngbvs[ii - 1] - 1] >> 8;
          lid_tmp = ngbvs[ii - 1] - 1;
          ngbvs_tmp = (mesh_v2hfid[lid_tmp] & 255UL);
          lid = ngbvs_tmp;
          pos =
              1 -
              ((mesh_elemtables[(
                                    mesh_teids[fid - 1] &
                                    255UL) -
                                1]
                        .etype >>
                    5 &
                7) == 2);
          //  Allow early termination of the loop if an incident halfedge
          c_tmp = hebuf[ii - 1];
          if ((c_tmp != 0UL) &&
              (mesh_sibhfs[ngbvs_tmp + mesh_sibhfs.size(1) *
                                           (fid - 1)] !=
               0UL)) {
            allow_early_term = true;
            fid = hebuf[ii - 1] >> 8;
            lid = (c_tmp & 255UL);
            pos =
                1 - ((mesh_elemtables[(
                                          mesh_teids[fid -
                                                     1] &
                                          255UL) -
                                      1]
                              .etype >>
                          5 &
                      7) == 2);
          } else {
            allow_early_term = false;
          }
          //  Starting point of counterclockwise rotation
          if (mesh_sibhfs[lid + mesh_sibhfs.size(1) *
                                    (fid - 1)] != 0UL) {
            fid_in = fid;
          } else {
            fid_in = 0UL;
            if (((bridges.size(0) == 0) || (bridges[lid_tmp] == 0)) &&
                (orig_ring == ring)) {
              //  Reflect rings to enlarge the number of remaining rings
              *reflected = true;
              ring = orig_ring + (orig_ring - cur_ring);
            }
          }
          ngbvs_tmp = lid + (pos << 2);
          v = mesh_elemtables[(
                                  mesh_teids[fid - 1] &
                                  255UL) -
                              1]
                  .conn
                      [(nxt[ngbvs_tmp] +
                        mesh_elemtables
                                [(
                                     mesh_teids[fid - 1] &
                                     255UL) -
                                 1]
                                    .conn.size(1) *
                            ((
                                 mesh_teids[fid - 1] >>
                                 8) -
                             1)) -
                       1] -
              1;
          if ((*overflow) || ((!vtags[v]) && (*nverts >= maxnpnts))) {
            *overflow = true;
          } else {
            *overflow = false;
          }
          if ((!*overflow) && (!vtags[v])) {
            (*nverts)++;
            ngbvs[*nverts - 1] = v + 1;
            vtags[v] = true;
            //  Save starting position for next vertex
            hebuf[*nverts - 1] = ((fid << 8) + nxt[ngbvs_tmp]) - 1UL;
          }
          //  Rotate counterclockwise around the vertex.
          isfirst = true;
          ::coder::SizeType exitg2;
          boolean_T guard2{false};
          do {
            exitg2 = 0;
            //  Insert vertx into list
            ngbvs_tmp = (
                            mesh_teids[fid - 1] & 255UL) -
                        1;
            lid_prv =
                prv[lid +
                    ((1 - ((mesh_elemtables[ngbvs_tmp].etype >> 5 & 7) == 2))
                     << 2)];
            //  Insert face into list
            guard2 = false;
            if (ftags[fid - 1]) {
              if (allow_early_term && (!isfirst)) {
                exitg2 = 1;
              } else {
                guard2 = true;
              }
            } else {
              //  If the face has already been inserted, then the vertex
              v = mesh_elemtables[ngbvs_tmp].conn
                      [(lid_prv +
                        mesh_elemtables[ngbvs_tmp].conn.size(1) *
                            ((
                                 mesh_teids[fid - 1] >>
                                 8) -
                             1)) -
                       1] -
                  1;
              if ((*overflow) || ((!vtags[v]) && (*nverts >= maxnpnts)) ||
                  ((!ftags[fid - 1]) &&
                   (*nfaces >= maxnf))) {
                *overflow = true;
              } else {
                *overflow = false;
              }
              if ((!vtags[v]) && (!*overflow)) {
                (*nverts)++;
                ngbvs[*nverts - 1] = v + 1;
                vtags[v] = true;
                //  Save starting position for next ring
                hebuf[*nverts - 1] = ((fid << 8) + lid_prv) - 1UL;
              }
              if ((!ftags[fid - 1]) && (!*overflow)) {
                (*nfaces)++;
                ngbfs[*nfaces - 1] = fid;
                ftags[fid - 1] = true;
              }
              isfirst = false;
              guard2 = true;
            }
            if (guard2) {
              opp =
                  mesh_sibhfs[(lid_prv + mesh_sibhfs.size(1) *
                                             (fid - 1)) -
                              1];
              fid = opp >> 8;
              if (fid == fid_in) {
                exitg2 = 1;
              } else {
                lid = (opp & 255UL);
              }
            }
          } while (exitg2 == 0);
          //  Handle bridges
          if (bridges.size(0) != 0) {
            ngbvs_tmp = bridges[ngbvs[ii - 1] - 1];
            if ((ngbvs_tmp > 0) && (!vtags[ngbvs_tmp - 1]) &&
                (*nverts < maxnpnts)) {
              (*nverts)++;
              ngbvs[*nverts - 1] = ngbvs_tmp;
              vtags[ngbvs_tmp - 1] = true;
            }
          }
        }
        cur_ring++;
        if ((cur_ring >= ring) || (*nfaces == nfaces_pre) || (*overflow)) {
          exitg1 = 1;
        } else {
          nverts_pre = nverts_last;
        }
      }
    } while (exitg1 == 0);
    //  Reset flags
    vtags[vid - 1] = false;
    for (::coder::SizeType b_i{0}; b_i < *nverts; b_i++) {
      vtags[ngbvs[b_i] - 1] = false;
    }
    for (::coder::SizeType b_i{0}; b_i < *nfaces; b_i++) {
      ftags[ngbfs[b_i] - 1] = false;
    }
  }
}

// omp4mRecurPartMesh Recursively partition an unstructured mesh
static void omp4mRecurPartMesh(::coder::SizeType n,
                               const ::coder::array<int32_T, 2U> &cells,
                               ::coder::SizeType dim, ::coder::SizeType nLevels, ::coder::SizeType nParts,
                               ::coder::array<Omp4mPart, 1U> &parts)
{
  ::coder::array<int32_T, 1U> cparts;
  ::coder::array<int32_T, 1U> eind;
  ::coder::array<int32_T, 1U> eindloc;
  ::coder::array<int32_T, 1U> eptr;
  ::coder::array<int32_T, 1U> eptrloc;
  ::coder::array<int32_T, 1U> iwork;
  ::coder::array<int32_T, 1U> l2gmap_;
  ::coder::array<int32_T, 1U> r;
  ::coder::array<int32_T, 1U> r1;
  ::coder::array<int32_T, 1U> vparts;
  ::coder::array<boolean_T, 1U> visited;
  Omp4mPart b_parts;
  int32_T nnodes;
  if ((n < 1) || ((cells.size(0) == 0) || (cells.size(1) == 0)) || (dim < 1) ||
      (dim > 3)) {
    parts.set_size(0);
  } else {
    //  check inputs
    if (nLevels < 1) {
      nLevels = dim;
    }
    //  quick return if possible
    if (nParts == 1) {
      b_parts.nparts = 1;
      b_parts.part_ptr.set_size(2);
      b_parts.part_ptr[0] = 1;
      b_parts.part_ptr[1] = n + 1;
      b_parts.part_list.set_size(0);
      b_parts.shared_ents.set_size(0);
      parts.set_size(1);
      parts[0] = b_parts;
      parts[0].part_list.set_size(n);
      for (::coder::SizeType i{0}; i < n; i++) {
        parts[0].part_list[i] = i + 1;
      }
    } else {
      ::coder::SizeType b_i;
      ::coder::SizeType idx;
      ::coder::SizeType m;
      ::coder::SizeType npts;
      ::coder::SizeType u0;
      ::coder::SizeType u1;
      parts.set_size(nLevels);
      //  buffer
      u0 = cells.size(0);
      u1 = n;
      if (u0 >= n) {
        u1 = u0;
      }
      iwork.set_size(0);
      eptr.set_size(n + 1);
      for (b_i = 0; b_i <= n; b_i++) {
        eptr[b_i] = 0;
      }
      eptr[0] = 1;
      //  determine number of incident cells
      m = cells.size(0) - 1;
      //  number of cells
      for (::coder::SizeType cid{0}; cid <= m; cid++) {
        //  local function to get number of points per cell
        for (npts = cells.size(1) - 1; cells[npts + cells.size(1) * cid] <= 0;
             npts--) {
        }
        for (::coder::SizeType i{0}; i <= npts; i++) {
          idx = cells[i + cells.size(1) * cid];
          eptr[idx] = eptr[idx] + 1;
        }
      }
      for (::coder::SizeType i{0}; i < n; i++) {
        eptr[i + 1] = eptr[i + 1] + eptr[i];
      }
      //  allocate n2cList
      eind.set_size(eptr[n] - 1);
      for (::coder::SizeType cid{0}; cid <= m; cid++) {
        //  local function to get number of points per cell
        for (npts = cells.size(1) - 1; cells[npts + cells.size(1) * cid] <= 0;
             npts--) {
        }
        for (::coder::SizeType i{0}; i <= npts; i++) {
          idx = cells[i + cells.size(1) * cid] - 1;
          eind[eptr[idx] - 1] = cid + 1;
          eptr[idx] = eptr[idx] + 1;
        }
      }
      for (::coder::SizeType i{n}; i >= 1; i--) {
        eptr[i] = eptr[i - 1];
      }
      eptr[0] = 1;
      //  call metis
      call_metis_mesh(cells.size(0), eptr, eind, nParts, vparts, cparts);
      //  build partition, note if nodal-based, then nparts is for cell
      visited.set_size(u1);
      for (b_i = 0; b_i < u1; b_i++) {
        visited[b_i] = false;
      }
      build_part(nParts, vparts, cparts, eptr, eind, visited, iwork,
                 parts[0].part_ptr, parts[0].part_list, parts[0].shared_ents);
      parts[0].nparts = nParts;
      if ((nLevels == 1) || (parts[0].shared_ents.size(0) == 0)) {
        if (nLevels > 1) {
          b_parts = parts[0];
          parts.set_size(1);
          parts[0] = b_parts;
        }
      } else {
        ::coder::SizeType lvl;
        boolean_T exitg1;
        //  for more than 2 levels, recursive loop begins
        lvl = 1;
        exitg1 = false;
        while ((!exitg1) && (lvl + 1 <= nLevels)) {
          l2gmap_.set_size(parts[lvl - 1].shared_ents.size(0));
          u0 = parts[lvl - 1].shared_ents.size(0);
          for (b_i = 0; b_i < u0; b_i++) {
            l2gmap_[b_i] = parts[lvl - 1].shared_ents[b_i];
          }
          parts[lvl - 1].shared_ents.set_size(0);
          //  first, extract and localize subgraph
          extract_sub(cells.size(0), l2gmap_, eptr, eind, iwork, visited,
                      eptrloc, eindloc, &nnodes);
          //  call metis
          call_metis_mesh(nnodes, eptrloc, eindloc, nParts, vparts, cparts);
          //  build partition, note if nodal-based, then nparts is for cell
          build_part(nParts, vparts, cparts, eptrloc, eindloc, visited, iwork,
                     parts[lvl].part_ptr, r, r1);
          parts[lvl].part_list.set_size(r.size(0));
          u0 = r.size(0);
          for (b_i = 0; b_i < u0; b_i++) {
            parts[lvl].part_list[b_i] = r[b_i];
          }
          parts[lvl].shared_ents.set_size(r1.size(0));
          u0 = r1.size(0);
          for (b_i = 0; b_i < u0; b_i++) {
            parts[lvl].shared_ents[b_i] = r1[b_i];
          }
          parts[lvl].nparts = nParts;
          //  from localized IDs to global IDs
          b_i = parts[lvl].part_list.size(0);
          for (::coder::SizeType i{0}; i < b_i; i++) {
            parts[lvl].part_list[i] = l2gmap_[parts[lvl].part_list[i] - 1];
          }
          b_i = parts[lvl].shared_ents.size(0);
          for (::coder::SizeType i{0}; i < b_i; i++) {
            parts[lvl].shared_ents[i] = l2gmap_[parts[lvl].shared_ents[i] - 1];
          }
          //  check
          if (parts[lvl].shared_ents.size(0) == 0) {
            if (lvl + 1 < nLevels) {
              parts.set_size(lvl + 1);
            }
            exitg1 = true;
          } else {
            lvl++;
          }
        }
        // Local buffers with patterns "\w+" are legitimate
      }
    }
  }
}

// prism_126 - Quintic prismatic element with equidistant nodes
static inline
void prism_126(real_T xi, real_T eta, real_T zeta, real_T sfvals[126],
                      real_T sdvals[378])
{
  ::sfe_sfuncs::prism_126_sfunc(xi, eta, zeta, &sfvals[0], &sdvals[0]);
}

// prism_196 - Sextic prismatic element with equidistant nodes
static inline
void prism_196(real_T xi, real_T eta, real_T zeta, real_T sfvals[196],
                      real_T sdvals[588])
{
  ::sfe_sfuncs::prism_196_sfunc(xi, eta, zeta, &sfvals[0], &sdvals[0]);
}

// prism_40 - Cubic prismatic element
static inline
void prism_40(real_T xi, real_T eta, real_T zeta, real_T sfvals[40],
                     real_T sdvals[120])
{
  ::sfe_sfuncs::prism_40_sfunc(xi, eta, zeta, &sfvals[0], &sdvals[0]);
}

// prism_75 - Quartic prismatic element with equidistant nodes
static inline
void prism_75(real_T xi, real_T eta, real_T zeta, real_T sfvals[75],
                     real_T sdvals[225])
{
  ::sfe_sfuncs::prism_75_sfunc(xi, eta, zeta, &sfvals[0], &sdvals[0]);
}

// prism_gl_126 - Quintic prismatic element with equidistant nodes
static inline
void prism_gl_126(real_T xi, real_T eta, real_T zeta, real_T sfvals[126],
                         real_T sdvals[378])
{
  ::sfe_sfuncs::prism_gl_126_sfunc(xi, eta, zeta, &sfvals[0], &sdvals[0]);
}

// prism_gl_40 - Quadratic prismatic element with Gauss-Lobatto nodes
static inline
void prism_gl_40(real_T xi, real_T eta, real_T zeta, real_T sfvals[40],
                        real_T sdvals[120])
{
  ::sfe_sfuncs::prism_gl_40_sfunc(xi, eta, zeta, &sfvals[0], &sdvals[0]);
}

// prism_gl_75 - Quartic prismatic element with Gauss-Lobatto nodes
static inline
void prism_gl_75(real_T xi, real_T eta, real_T zeta, real_T sfvals[75],
                        real_T sdvals[225])
{
  ::sfe_sfuncs::prism_gl_75_sfunc(xi, eta, zeta, &sfvals[0], &sdvals[0]);
}

// pyra_30 - Compute shape functions and their derivatives of pyra_30
static inline
void pyra_30(real_T xi, real_T eta, real_T zeta, real_T sfvals[30],
                    real_T sdvals[90])
{
  ::sfe_sfuncs::pyra_30_sfunc(xi, eta, zeta, &sfvals[0], &sdvals[0]);
}

// pyra_55 - Compute shape functions and their derivatives of pyra_55
static inline
void pyra_55(real_T xi, real_T eta, real_T zeta, real_T sfvals[55],
                    real_T sdvals[165])
{
  ::sfe_sfuncs::pyra_55_sfunc(xi, eta, zeta, &sfvals[0], &sdvals[0]);
}

// pyra_gl_30 - Compute shape functions and their derivatives of pyra_gl_30
static inline
void pyra_gl_30(real_T xi, real_T eta, real_T zeta, real_T sfvals[30],
                       real_T sdvals[90])
{
  ::sfe_sfuncs::pyra_gl_30_sfunc(xi, eta, zeta, &sfvals[0], &sdvals[0]);
}

// pyra_gl_55 - Compute shape functions and their derivatives of pyra_gl_55
static inline
void pyra_gl_55(real_T xi, real_T eta, real_T zeta, real_T sfvals[55],
                       real_T sdvals[165])
{
  ::sfe_sfuncs::pyra_gl_55_sfunc(xi, eta, zeta, &sfvals[0], &sdvals[0]);
}

//  pyra_quadrules - Obtain quadrature points and weights of a pyramidal
static void pyra_quadrules(::coder::SizeType degree, ::coder::array<real_T, 2U> &cs,
                           ::coder::array<real_T, 1U> &ws)
{
  if (degree <= 1) {
    ::coder::SizeType nqp;
    nqp = ::sfe_qrules::pyra_deg1_qrule();
    cs.set_size(nqp, 3);
    ws.set_size(nqp);
    ::sfe_qrules::pyra_deg1_qrule(&cs[0], &(ws.data())[0]);
  } else if (degree <= 2) {
    ::coder::SizeType nqp;
    nqp = ::sfe_qrules::pyra_deg2_qrule();
    cs.set_size(nqp, 3);
    ws.set_size(nqp);
    ::sfe_qrules::pyra_deg2_qrule(&cs[0], &(ws.data())[0]);
  } else if (degree <= 3) {
    ::coder::SizeType nqp;
    nqp = ::sfe_qrules::pyra_deg3_qrule();
    cs.set_size(nqp, 3);
    ws.set_size(nqp);
    ::sfe_qrules::pyra_deg3_qrule(&cs[0], &(ws.data())[0]);
  } else if (degree <= 5) {
    ::coder::SizeType nqp;
    nqp = ::sfe_qrules::pyra_deg5_qrule();
    cs.set_size(nqp, 3);
    ws.set_size(nqp);
    ::sfe_qrules::pyra_deg5_qrule(&cs[0], &(ws.data())[0]);
  } else if (degree <= 6) {
    ::coder::SizeType nqp;
    nqp = ::sfe_qrules::pyra_deg6_qrule();
    cs.set_size(nqp, 3);
    ws.set_size(nqp);
    ::sfe_qrules::pyra_deg6_qrule(&cs[0], &(ws.data())[0]);
  } else if (degree <= 7) {
    ::coder::SizeType nqp;
    nqp = ::sfe_qrules::pyra_deg7_qrule();
    cs.set_size(nqp, 3);
    ws.set_size(nqp);
    ::sfe_qrules::pyra_deg7_qrule(&cs[0], &(ws.data())[0]);
  } else if (degree <= 9) {
    ::coder::SizeType nqp;
    nqp = ::sfe_qrules::pyra_deg9_qrule();
    cs.set_size(nqp, 3);
    ws.set_size(nqp);
    ::sfe_qrules::pyra_deg9_qrule(&cs[0], &(ws.data())[0]);
  } else if (degree <= 11) {
    ::coder::SizeType nqp;
    nqp = ::sfe_qrules::pyra_deg11_qrule();
    cs.set_size(nqp, 3);
    ws.set_size(nqp);
    ::sfe_qrules::pyra_deg11_qrule(&cs[0], &(ws.data())[0]);
  } else {
    ::coder::SizeType nqp;
    if (degree > 13) {
      m2cWarnMsgIdAndTxt("pyra_quadrules:UnsupportedDegree",
                         "Only support up to degree 13");
    }
    nqp = ::sfe_qrules::pyra_deg13_qrule();
    cs.set_size(nqp, 3);
    ws.set_size(nqp);
    ::sfe_qrules::pyra_deg13_qrule(&cs[0], &(ws.data())[0]);
  }
}

// quad_25 - Biquartic quadrilateral element with equidistant points
static inline
void quad_25(real_T xi, real_T eta, real_T sfvals[25], real_T sdvals[50])
{
  ::sfe_sfuncs::quad_25_sfunc(xi, eta, &sfvals[0], &sdvals[0]);
}

// quad_36   Biquintic quadrilateral element with equidistant points
static inline
void quad_36(real_T xi, real_T eta, real_T sfvals[36], real_T sdvals[72])
{
  ::sfe_sfuncs::quad_36_sfunc(xi, eta, &sfvals[0], &sdvals[0]);
}

// quad_49 - Bisextic quadrilateral element with equidistant points
static inline
void quad_49(real_T xi, real_T eta, real_T sfvals[49], real_T sdvals[98])
{
  ::sfe_sfuncs::quad_49_sfunc(xi, eta, &sfvals[0], &sdvals[0]);
}

// quad_gl_25 - Biquartic quadrilateral element with Gauss-Lobatto points
static inline
void quad_gl_25(real_T xi, real_T eta, real_T sfvals[25],
                       real_T sdvals[50])
{
  ::sfe_sfuncs::quad_gl_25_sfunc(xi, eta, &sfvals[0], &sdvals[0]);
}

// quad_gl_36 - Biquintic quadrilateral element with equidistant points
static inline
void quad_gl_36(real_T xi, real_T eta, real_T sfvals[36],
                       real_T sdvals[72])
{
  ::sfe_sfuncs::quad_gl_36_sfunc(xi, eta, &sfvals[0], &sdvals[0]);
}

// quad_gl_49 - Bisextic quadrilateral element with equidistant points
static inline
void quad_gl_49(real_T xi, real_T eta, real_T sfvals[49],
                       real_T sdvals[98])
{
  ::sfe_sfuncs::quad_gl_49_sfunc(xi, eta, &sfvals[0], &sdvals[0]);
}

// rdi_compute_oscind - Compute oscillation indicators (beta values)
static void rdi_compute_oscind(
    real_T rdi_epsbeta, const ::coder::array<real_T, 2U> &mesh_coords,
    const ::coder::array<int64_T, 1U> &mesh_node2elems_row_ptr,
    const ::coder::array<int32_T, 1U> &mesh_node2elems_col_ind,
    const ::coder::array<real_T, 1U> &mesh_elemmeas, real_T mesh_globalh,
    const ::coder::array<real_T, 2U> &df,
    const ::coder::array<real_T, 2U> &alpha, ::coder::array<real_T, 2U> &beta)
{
  real_T epsbeta;
  if (alpha.size(1) != df.size(1)) {
    m2cErrMsgIdAndTxt("rdi_compute_oscind:badShape",
                      "Unmatched # of functions");
  }
  epsbeta = rdi_epsbeta;
  if (rdi_epsbeta <= 0.0) {
    epsbeta = 0.0002;
  }
  compute_beta_kernel(mesh_coords, mesh_node2elems_row_ptr,
                      mesh_node2elems_col_ind, mesh_elemmeas, mesh_globalh, df,
                      alpha, epsbeta, beta);
}

// rdi_compute_osusind - Compute over-/under-shoot indicators
static inline
void rdi_compute_osusind(
    const ::coder::array<int64_T, 1U> &rdi_A_row_ptr,
    const ::coder::array<int32_T, 1U> &rdi_A_col_ind,
    const ::coder::array<real_T, 1U> &rdi_A_val, boolean_T rdi_fullrank,
    const ::coder::array<real_T, 2U> &mesh_coords,
    const ::coder::array<int64_T, 1U> &mesh_node2elems_row_ptr,
    const ::coder::array<int32_T, 1U> &mesh_node2elems_col_ind,
    const ::coder::array<real_T, 2U> &fs, ::coder::array<real_T, 2U> &alphacell,
    ::coder::array<real_T, 2U> &alphanode)
{
  if (!rdi_fullrank) {
    m2cErrMsgIdAndTxt("rdi_compute_osusind:badA",
                      "OSUS operator must be full-rank");
  }
  if (mesh_coords.size(0) != fs.size(0)) {
    m2cErrMsgIdAndTxt("rdi_compute_osusind:badShape",
                      "Unmatched input fs size");
  }
  //  Using advanced OpenMP mode to compute cell-based alpha value
  crsProdMatVec(rdi_A_row_ptr, rdi_A_col_ind, rdi_A_val, fs, alphacell);
#pragma omp barrier
  compute_nodal_alpha(alphacell, mesh_coords, mesh_node2elems_row_ptr,
                      mesh_node2elems_col_ind, alphanode);
}

static void
rdi_contract_markers(::coder::array<int8_T, 2U> &distags,
                     const ::coder::array<real_T, 2U> &mesh_coords,
                     const ::coder::array<ConnData, 1U> &mesh_elemtables,
                     const ::coder::array<uint64_T, 1U> &mesh_teids,
                     const ::coder::array<int64_T, 1U> &mesh_node2nodes_row_ptr,
                     const ::coder::array<int32_T, 1U> &mesh_node2nodes_col_ind,
                     const ::coder::array<int64_T, 1U> &mesh_node2elems_row_ptr,
                     const ::coder::array<int32_T, 1U> &mesh_node2elems_col_ind,
                     ::coder::SizeType nlayers)
{
  ::coder::array<int64_T, 1U> lgraph__row_ptr;
  ::coder::array<int32_T, 1U> bdnodes_;
  ::coder::array<int32_T, 1U> iwork_;
  ::coder::array<int32_T, 1U> lgraph__col_ind;
  ::coder::array<boolean_T, 1U> visited_;
  int32_T lgraph__ncols;
  if (nlayers > 0) {
    ::coder::SizeType loop_ub;
    ::coder::SizeType nrhs;
    ::coder::SizeType u1;
    if (distags.size(0) != mesh_coords.size(0)) {
      m2cErrMsgIdAndTxt("rdi_contract_markers:badSize", "unmatched nnodes");
    }
    nrhs = distags.size(1);
    //  Create workspace
    visited_.set_size(distags.size(0));
    loop_ub = distags.size(0);
    for (u1 = 0; u1 < loop_ub; u1++) {
      visited_[u1] = false;
    }
    loop_ub = mesh_coords.size(0);
    u1 = mesh_teids.size(0);
    if (loop_ub >= u1) {
      u1 = loop_ub;
    }
    iwork_.set_size(u1);
    for (::coder::SizeType col{0}; col < nrhs; col++) {
      for (::coder::SizeType layer{0}; layer < nlayers; layer++) {
        ::coder::SizeType nbp;
        ::coder::SizeType npts;
        //  Kernel for determine boundary nodes
        npts = distags.size(0);
        bdnodes_.set_size(distags.size(0));
        nbp = 0;
        for (::coder::SizeType i{0}; i < npts; i++) {
          if (distags[col + distags.size(1) * i] != 0) {
            int64_T j;
            boolean_T exitg1;
            j = mesh_node2nodes_row_ptr[i];
            exitg1 = false;
            while ((!exitg1) && (j <= mesh_node2nodes_row_ptr[i + 1] - 1L)) {
              if (distags[col +
                          distags.size(1) *
                              (mesh_node2nodes_col_ind[j -
                                                       1] -
                               1)] == 0) {
                nbp++;
                bdnodes_[nbp - 1] = i + 1;
                exitg1 = true;
              } else {
                j++;
              }
            }
          }
        }
        if (nbp < 1) {
          u1 = 0;
        } else {
          u1 = nbp;
        }
        bdnodes_.set_size(u1);
        for (::coder::SizeType i{0}; i < u1; i++) {
          int8_T b_i;
          b_i = distags[col + distags.size(1) * (bdnodes_[i] - 1)];
          if (b_i < 0) {
            ::coder::SizeType n1;
            ::coder::SizeType n2;
            compute_fconn_graph(visited_, iwork_, mesh_elemtables, mesh_teids,
                                mesh_node2elems_row_ptr,
                                mesh_node2elems_col_ind, bdnodes_[i], distags,
                                col + 1, lgraph__row_ptr, lgraph__col_ind,
                                &lgraph__ncols);
            n1 = compute_connected_components(visited_, iwork_, lgraph__row_ptr,
                                              lgraph__col_ind, lgraph__ncols);
            distags[col + distags.size(1) * (bdnodes_[i] - 1)] = 0;
            compute_fconn_graph(visited_, iwork_, mesh_elemtables, mesh_teids,
                                mesh_node2elems_row_ptr,
                                mesh_node2elems_col_ind, bdnodes_[i], distags,
                                col + 1, lgraph__row_ptr, lgraph__col_ind,
                                &lgraph__ncols);
            n2 = compute_connected_components(visited_, iwork_, lgraph__row_ptr,
                                              lgraph__col_ind, lgraph__ncols);
            // Local buffers with patterns "(\w+_)?(row_ptr|col_ind)" are legitimate
            if (n1 != n2) {
              distags[col + distags.size(1) * (bdnodes_[i] - 1)] = b_i;
            }
          }
        }
      }
    }
  }
}

static void
rdi_expand_markers(::coder::array<int8_T, 2U> &distags,
                   const ::coder::array<int64_T, 1U> &mesh_node2nodes_row_ptr,
                   const ::coder::array<int32_T, 1U> &mesh_node2nodes_col_ind,
                   ::coder::SizeType nlayers)
{
  ::coder::array<int32_T, 1U> expmarkers_;
  ::coder::array<int32_T, 1U> lvlptr_;
  if (nlayers > 0) {
    ::coder::SizeType nrhs;
    nrhs = distags.size(1);
    expmarkers_.set_size(distags.size(0));
    lvlptr_.set_size(nlayers + 1);
    lvlptr_[0] = 1;
    for (::coder::SizeType col{0}; col < nrhs; col++) {
      int64_T j;
      ::coder::SizeType i;
      ::coder::SizeType i1;
      ::coder::SizeType tag_tmp_tmp;
      int8_T tag;
      //  Resets the counter if nrhs > 1
      lvlptr_[1] = lvlptr_[0];
      i = distags.size(0);
      for (::coder::SizeType b_i{0}; b_i < i; b_i++) {
        tag = distags[col + distags.size(1) * b_i];
        if (tag > 0) {
          //  Get connectivity of this node to get 1-ring, and mark them
          for (j = mesh_node2nodes_row_ptr[b_i];
               j < mesh_node2nodes_row_ptr[b_i + 1]; j++) {
            //  Preserve originally marked nodes and prefer C_0 over C_1
            i1 = mesh_node2nodes_col_ind[j - 1];
            tag_tmp_tmp = distags[col + distags.size(1) * (i1 - 1)];
            if ((tag_tmp_tmp == -2) || (tag_tmp_tmp == 0)) {
              distags[col + distags.size(1) * (i1 - 1)] =
                  static_cast<int8_T>(-tag);
              expmarkers_[lvlptr_[1] - 1] = i1;
              lvlptr_[1] = lvlptr_[1] + 1;
            }
          }
        }
      }
      for (::coder::SizeType k{2}; k <= nlayers; k++) {
        lvlptr_[k] = lvlptr_[k - 1];
        i = lvlptr_[k - 2];
        i1 = lvlptr_[k - 1] - 1;
        for (::coder::SizeType b_i{i}; b_i <= i1; b_i++) {
          int64_T i2;
          tag_tmp_tmp = expmarkers_[b_i - 1];
          tag = distags[col + distags.size(1) * (tag_tmp_tmp - 1)];
          //  Get connectivity of this node to get 1-ring, and mark them
          j = mesh_node2nodes_row_ptr[tag_tmp_tmp - 1];
          i2 = mesh_node2nodes_row_ptr[tag_tmp_tmp] - 1L;
          while (j <= i2) {
            ::coder::SizeType i3;
            //  Prefer C_0 over C_1
            tag_tmp_tmp = mesh_node2nodes_col_ind[j - 1];
            i3 = distags[col + distags.size(1) * (tag_tmp_tmp - 1)];
            if ((i3 == -2) || (i3 == 0)) {
              if (tag < 0) {
                distags[col + distags.size(1) * (tag_tmp_tmp - 1)] =
                    static_cast<int8_T>(-static_cast<int8_T>(-tag));
              } else {
                distags[col + distags.size(1) * (tag_tmp_tmp - 1)] =
                    static_cast<int8_T>(-tag);
              }
              expmarkers_[lvlptr_[k] - 1] = tag_tmp_tmp;
              lvlptr_[k] = lvlptr_[k] + 1;
            }
            j++;
          }
        }
      }
    }
  }
}

// rdi_mark_discontinuities - Determine discontinuous nodes
static void rdi_mark_discontinuities(
    real_T rdi_cglobal, real_T rdi_clocal, real_T rdi_kappa0, real_T rdi_kappa1,
    const ::coder::array<real_T, 2U> &mesh_coords,
    const ::coder::array<ConnData, 1U> &mesh_elemtables,
    const ::coder::array<uint64_T, 1U> &mesh_teids,
    const ::coder::array<int64_T, 1U> &mesh_node2elems_row_ptr,
    const ::coder::array<int32_T, 1U> &mesh_node2elems_col_ind,
    const ::coder::array<real_T, 1U> &mesh_elemh, real_T mesh_globalh,
    const ::coder::array<real_T, 2U> &fs, const ::coder::array<real_T, 2U> &df,
    const ::coder::array<real_T, 2U> &alpha,
    const ::coder::array<real_T, 2U> &beta, ::coder::array<int8_T, 2U> &distags)
{
  real_T cglobal;
  real_T clocal;
  real_T kappa0;
  real_T kappa1;
  if (mesh_coords.size(0) != fs.size(0)) {
    m2cErrMsgIdAndTxt("rdi_mark_discontinuities:badShape",
                      "Unmatched input fs size");
  }
  if ((df.size(1) != fs.size(1)) || (alpha.size(1) != df.size(1)) ||
      (beta.size(1) != df.size(1))) {
    m2cErrMsgIdAndTxt("rdi_mark_discontinuities:badShape",
                      "Unmatched # of functions");
  }
  //  Get parameters
  cglobal = rdi_cglobal;
  if (rdi_cglobal < 0.0) {
    cglobal = 0.05;
  }
  clocal = rdi_clocal;
  if (rdi_clocal < 0.0) {
    clocal = 0.5;
  }
  kappa1 = rdi_kappa1;
  if (rdi_kappa1 <= 0.0) {
    kappa1 = 0.3;
  }
  //  C1 dis
  kappa0 = rdi_kappa0;
  if (rdi_kappa0 <= 0.0) {
    kappa0 = 1.0;
  }
  //  C0 dis
  if (kappa0 <= kappa1) {
    kappa0 = kappa1;
  }
  mark_discontinuities_kernel(mesh_coords, mesh_elemtables, mesh_teids,
                              mesh_node2elems_row_ptr, mesh_node2elems_col_ind,
                              mesh_elemh, mesh_globalh, fs, df, alpha, beta,
                              cglobal, clocal, kappa1, kappa0, distags);
}

// rdi_update_osusop - Update OSUS operator for rank-deficient nodes
static void rdi_update_osusop(RdiObject *rdi, WlsMesh *mesh, boolean_T interp0)
{
  ::coder::SizeType extstclid;
  m2cAssert(!rdi->fullrank, "Cannot be full-rank");
  //  Compute stencils with larger ring size for extended RD nodes
  rdi->ring += 1.0 / static_cast<real_T>(mesh->topo_ndims);
  wlsmesh_compute_stencils(mesh, rdi->extstclid, rdi->ring);
  //  Update the operator sparsity pattern
  update_osusop(&rdi->A, rdi->nnzs, mesh->coords, mesh->node2elems.row_ptr,
                mesh->node2elems.col_ind, mesh->stencils, rdi->extstclid);
  if (mesh->topo_ndims == mesh->coords.size(1)) {
    //  {1,2,3}-D body assembly
    assemble_body(rdi, mesh->coords, mesh->elemtables, mesh->teids,
                  mesh->node2elems.row_ptr, mesh->node2elems.col_ind,
                  mesh->stencils, rdi->extstclid, interp0);
  } else {
    if (mesh->topo_ndims + 1 != mesh->coords.size(1)) {
      m2cErrMsgIdAndTxt("rdi_assemble_osusop:badDim",
                        "Must be either body or surface");
    }
    assemble_surf(rdi, mesh->coords, mesh->elemtables, mesh->teids,
                  mesh->nrmstables, mesh->node2elems.row_ptr,
                  mesh->node2elems.col_ind, mesh->stencils, rdi->extstclid,
                  interp0);
  }
  //  Determine potential RD nodes and reset rdtags. RD nodes are stored in
  extstclid = rdi->extstclid;
  //  Overcome a "bug" in coder regarding creating buffers
  determine_rdnodes(rdi->fullrank, rdi->rdtags,
                    mesh->stencils[extstclid - 1].vidmap);
  if (rdi->fullrank) {
    crsCompress(&rdi->A, rdi->nnzs);
  }
}

//  rrqr_factor  Compute rank-revealing QR with column pivoting
static void rrqr_factor(const ::coder::array<real_T, 2U> &A, real_T thres,
                        ::coder::SizeType rowoffset, ::coder::SizeType coloffset, ::coder::SizeType m,
                        ::coder::SizeType n, ::coder::array<real_T, 2U> &QR,
                        ::coder::array<int32_T, 1U> &p, int32_T *rank,
                        ::coder::array<real_T, 1U> &work)
{
  ::coder::SizeType wsize;
  if (m == 0) {
    m = A.size(1) - rowoffset;
  } else {
    m2cAssert(m + rowoffset <= A.size(1),
              "Number of rows cannot exceed nrows(A).");
  }
  if (n == 0) {
    n = A.size(0) - coloffset;
  } else {
    m2cAssert(n + coloffset <= A.size(0),
              "Number of ncolumns cannot exceed ncols(A).");
  }
  //  Preallocate output arguments
  m2cAssert(QR.size(1) == A.size(1),
            "The number of rows in QR must be equal to that of A.");
  m2cAssert(QR.size(0) >= n + 1,
            "The number of columns in QR must be greater than that of A.");
  m2cAssert(p.size(0) >= n, "Length of permutation vector must be no smaller "
                            "than the number of columns.");
  //  Allocate work space if needed
  wsize = wls::query_work_size(m, n);
  work.set_size(wsize);
  p[0] = 0;
  *rank = wls::rrqr_factor(&A[rowoffset + A.size(1) * coloffset], thres, m, n,
                           &QR[0], &(p.data())[0], &(work.data())[0], wsize,
                           A.size(1));
}

//  rrqr_qmulti  Perform Q*bs, where Q is stored implicitly in QR
static void rrqr_qmulti(const ::coder::array<real_T, 2U> &QR, ::coder::SizeType m,
                        ::coder::SizeType n, ::coder::SizeType rank, ::coder::array<real_T, 2U> &bs,
                        ::coder::SizeType nrhs, ::coder::array<real_T, 1U> &work)
{
  ::coder::SizeType stride_bs;
  ::coder::SizeType u1;
  ::coder::SizeType wsize;
  stride_bs = bs.size(1);
  //  Obtain input arguments
  if (m == 0) {
    m = QR.size(1);
  }
  if (n == 0) {
    n = QR.size(0) - 1;
  }
  if (rank == 0) {
    rank = n;
  }
  u1 = n;
  if (m <= n) {
    u1 = m;
  }
  if ((rank > u1) || (rank < 1)) {
    m2cErrMsgIdAndTxt(
        "wlslib:WrongRank",
        "Rank %d must be a positive value no greater than min(%d, %d).", (int)rank,
        (int)m, (int)n);
  }
  if (nrhs == 0) {
    nrhs = bs.size(0);
  }
  //  Resize work space if needed
  wsize = wls::query_work_size(m, n);
  work.set_size(wsize);
  u1 = n + 1;
  for (::coder::SizeType i{u1}; i <= m; i++) {
    for (::coder::SizeType j{0}; j < nrhs; j++) {
      bs[(i + bs.size(1) * j) - 1] = 0.0;
    }
  }
  //  Invoke C++ function
  wls::rrqr_qmulti(&QR[0], m, n, rank, QR.size(1), nrhs, &bs[0], stride_bs,
                   &(work.data())[0], wsize);
}

//  rrqr_rtsolve  Perform forward substitution to compute bs=R'\bs, where R is
static void rrqr_rtsolve(const ::coder::array<real_T, 2U> &QR, ::coder::SizeType n,
                         ::coder::SizeType rank, ::coder::array<real_T, 2U> &bs,
                         ::coder::SizeType nrhs)
{
  ::coder::SizeType i;
  if (n == 0) {
    n = QR.size(0) - 1;
  }
  if (rank == 0) {
    rank = n;
  }
  if (QR.size(1) > n) {
    i = n;
  } else {
    i = QR.size(1);
  }
  if ((rank > i) || (rank < 1)) {
    m2cErrMsgIdAndTxt(
        "wlslib:WrongRank",
        "Rank %d must be a positive value no greater than min(%d, %d).", (int)rank,
        (int)QR.size(1), (int)n);
  }
  if (nrhs == 0) {
    nrhs = bs.size(0);
  }
  //  Obtain stride
  wls::rrqr_rtsolve(&QR[0], n, rank, QR.size(1), nrhs, &bs[0], bs.size(1));
}

static void sfe1_tabulate_shapefuncs(::coder::SizeType etype,
                                     const ::coder::array<real_T, 2U> &cs,
                                     ::coder::array<real_T, 2U> &sfvals,
                                     ::coder::array<real_T, 3U> &sdvals)
{
  ::coder::SizeType i;
  ::coder::SizeType nqp;
  //  Tabulate shape functions and derivative at given points.
  nqp = cs.size(0) - 1;
  i = iv[etype - 1];
  sfvals.set_size(cs.size(0), i);
  sdvals.set_size(cs.size(0), i, cs.size(1));
  if ((etype & 3) == 0) {
    switch (etype) {
    case 36: {
      ::coder::SizeType i1;
      boolean_T b;
      boolean_T b1;
      b = true;
      b1 = ((cs.size(1) <= 0) || (cs.size(0) <= 0));
      i = cs.size(1) * cs.size(0);
      i1 = 0;
      for (::coder::SizeType q{0}; q <= nqp; q++) {
        real_T N[2];
        real_T deriv[2];
        if (b1 || (q >= i)) {
          i1 = 0;
          b = true;
        } else if (b) {
          b = false;
          i1 = q % cs.size(0) * cs.size(1) + q / cs.size(0);
        } else {
          ::coder::SizeType i2;
          i2 = cs.size(1) * cs.size(0) - 1;
          if (i1 > MAX_int32_T - cs.size(1)) {
            i1 = q % cs.size(0) * cs.size(1) + q / cs.size(0);
          } else {
            i1 += cs.size(1);
            if (i1 > i2) {
              i1 -= i2;
            }
          }
        }
        ::sfe_sfuncs::bar_2_sfunc(cs[i1], &N[0], &deriv[0]);
        sfvals[sfvals.size(1) * q] = N[0];
        sdvals[sdvals.size(2) * sdvals.size(1) * q] = deriv[0];
        sfvals[sfvals.size(1) * q + 1] = N[1];
        sdvals[sdvals.size(2) + sdvals.size(2) * sdvals.size(1) * q] = deriv[1];
      }
    } break;
    case 40: {
      ::coder::SizeType i1;
      boolean_T b;
      boolean_T b1;
      b = true;
      b1 = ((cs.size(1) <= 0) || (cs.size(0) <= 0));
      i = cs.size(1) * cs.size(0);
      i1 = 0;
      for (::coder::SizeType q{0}; q <= nqp; q++) {
        real_T b_N[3];
        real_T b_deriv[3];
        if (b1 || (q >= i)) {
          i1 = 0;
          b = true;
        } else if (b) {
          b = false;
          i1 = q % cs.size(0) * cs.size(1) + q / cs.size(0);
        } else {
          ::coder::SizeType i2;
          i2 = cs.size(1) * cs.size(0) - 1;
          if (i1 > MAX_int32_T - cs.size(1)) {
            i1 = q % cs.size(0) * cs.size(1) + q / cs.size(0);
          } else {
            i1 += cs.size(1);
            if (i1 > i2) {
              i1 -= i2;
            }
          }
        }
        ::sfe_sfuncs::bar_3_sfunc(cs[i1], &b_N[0], &b_deriv[0]);
        sfvals[sfvals.size(1) * q] = b_N[0];
        sdvals[sdvals.size(2) * sdvals.size(1) * q] = b_deriv[0];
        sfvals[sfvals.size(1) * q + 1] = b_N[1];
        sdvals[sdvals.size(2) + sdvals.size(2) * sdvals.size(1) * q] =
            b_deriv[1];
        sfvals[sfvals.size(1) * q + 2] = b_N[2];
        sdvals[sdvals.size(2) * 2 + sdvals.size(2) * sdvals.size(1) * q] =
            b_deriv[2];
      }
    } break;
    case 44: {
      ::coder::SizeType i1;
      boolean_T b;
      boolean_T b1;
      b = true;
      b1 = ((cs.size(1) <= 0) || (cs.size(0) <= 0));
      i = cs.size(1) * cs.size(0);
      i1 = 0;
      for (::coder::SizeType q{0}; q <= nqp; q++) {
        real_T c_N[4];
        real_T c_deriv[4];
        if (b1 || (q >= i)) {
          i1 = 0;
          b = true;
        } else if (b) {
          b = false;
          i1 = q % cs.size(0) * cs.size(1) + q / cs.size(0);
        } else {
          ::coder::SizeType i2;
          i2 = cs.size(1) * cs.size(0) - 1;
          if (i1 > MAX_int32_T - cs.size(1)) {
            i1 = q % cs.size(0) * cs.size(1) + q / cs.size(0);
          } else {
            i1 += cs.size(1);
            if (i1 > i2) {
              i1 -= i2;
            }
          }
        }
        ::sfe_sfuncs::bar_4_sfunc(cs[i1], &c_N[0], &c_deriv[0]);
        sfvals[sfvals.size(1) * q] = c_N[0];
        sdvals[sdvals.size(2) * sdvals.size(1) * q] = c_deriv[0];
        sfvals[sfvals.size(1) * q + 1] = c_N[1];
        sdvals[sdvals.size(2) + sdvals.size(2) * sdvals.size(1) * q] =
            c_deriv[1];
        sfvals[sfvals.size(1) * q + 2] = c_N[2];
        sdvals[sdvals.size(2) * 2 + sdvals.size(2) * sdvals.size(1) * q] =
            c_deriv[2];
        sfvals[sfvals.size(1) * q + 3] = c_N[3];
        sdvals[sdvals.size(2) * 3 + sdvals.size(2) * sdvals.size(1) * q] =
            c_deriv[3];
      }
    } break;
    case 48: {
      ::coder::SizeType i1;
      boolean_T b;
      boolean_T b1;
      b = true;
      b1 = ((cs.size(1) <= 0) || (cs.size(0) <= 0));
      i = cs.size(1) * cs.size(0);
      i1 = 0;
      for (::coder::SizeType q{0}; q <= nqp; q++) {
        real_T d_N[5];
        real_T d_deriv[5];
        if (b1 || (q >= i)) {
          i1 = 0;
          b = true;
        } else if (b) {
          b = false;
          i1 = q % cs.size(0) * cs.size(1) + q / cs.size(0);
        } else {
          ::coder::SizeType i2;
          i2 = cs.size(1) * cs.size(0) - 1;
          if (i1 > MAX_int32_T - cs.size(1)) {
            i1 = q % cs.size(0) * cs.size(1) + q / cs.size(0);
          } else {
            i1 += cs.size(1);
            if (i1 > i2) {
              i1 -= i2;
            }
          }
        }
        ::sfe_sfuncs::bar_5_sfunc(cs[i1], &d_N[0], &d_deriv[0]);
        sfvals[sfvals.size(1) * q] = d_N[0];
        sdvals[sdvals.size(2) * sdvals.size(1) * q] = d_deriv[0];
        sfvals[sfvals.size(1) * q + 1] = d_N[1];
        sdvals[sdvals.size(2) + sdvals.size(2) * sdvals.size(1) * q] =
            d_deriv[1];
        sfvals[sfvals.size(1) * q + 2] = d_N[2];
        sdvals[sdvals.size(2) * 2 + sdvals.size(2) * sdvals.size(1) * q] =
            d_deriv[2];
        sfvals[sfvals.size(1) * q + 3] = d_N[3];
        sdvals[sdvals.size(2) * 3 + sdvals.size(2) * sdvals.size(1) * q] =
            d_deriv[3];
        sfvals[sfvals.size(1) * q + 4] = d_N[4];
        sdvals[sdvals.size(2) * 4 + sdvals.size(2) * sdvals.size(1) * q] =
            d_deriv[4];
      }
    } break;
    case 52: {
      ::coder::SizeType i1;
      boolean_T b;
      boolean_T b1;
      b = true;
      b1 = ((cs.size(1) <= 0) || (cs.size(0) <= 0));
      i = cs.size(1) * cs.size(0);
      i1 = 0;
      for (::coder::SizeType q{0}; q <= nqp; q++) {
        real_T e_N[6];
        real_T e_deriv[6];
        if (b1 || (q >= i)) {
          i1 = 0;
          b = true;
        } else if (b) {
          b = false;
          i1 = q % cs.size(0) * cs.size(1) + q / cs.size(0);
        } else {
          ::coder::SizeType i2;
          i2 = cs.size(1) * cs.size(0) - 1;
          if (i1 > MAX_int32_T - cs.size(1)) {
            i1 = q % cs.size(0) * cs.size(1) + q / cs.size(0);
          } else {
            i1 += cs.size(1);
            if (i1 > i2) {
              i1 -= i2;
            }
          }
        }
        ::sfe_sfuncs::bar_6_sfunc(cs[i1], &e_N[0], &e_deriv[0]);
        sfvals[sfvals.size(1) * q] = e_N[0];
        sdvals[sdvals.size(2) * sdvals.size(1) * q] = e_deriv[0];
        sfvals[sfvals.size(1) * q + 1] = e_N[1];
        sdvals[sdvals.size(2) + sdvals.size(2) * sdvals.size(1) * q] =
            e_deriv[1];
        sfvals[sfvals.size(1) * q + 2] = e_N[2];
        sdvals[sdvals.size(2) * 2 + sdvals.size(2) * sdvals.size(1) * q] =
            e_deriv[2];
        sfvals[sfvals.size(1) * q + 3] = e_N[3];
        sdvals[sdvals.size(2) * 3 + sdvals.size(2) * sdvals.size(1) * q] =
            e_deriv[3];
        sfvals[sfvals.size(1) * q + 4] = e_N[4];
        sdvals[sdvals.size(2) * 4 + sdvals.size(2) * sdvals.size(1) * q] =
            e_deriv[4];
        sfvals[sfvals.size(1) * q + 5] = e_N[5];
        sdvals[sdvals.size(2) * 5 + sdvals.size(2) * sdvals.size(1) * q] =
            e_deriv[5];
      }
    } break;
    default: {
      ::coder::SizeType i1;
      boolean_T b;
      boolean_T b1;
      m2cAssert(etype == 56, "Only support up to sextic.");
      b = true;
      b1 = ((cs.size(1) <= 0) || (cs.size(0) <= 0));
      i = cs.size(1) * cs.size(0);
      i1 = 0;
      for (::coder::SizeType q{0}; q <= nqp; q++) {
        real_T f_N[7];
        real_T f_deriv[7];
        if (b1 || (q >= i)) {
          i1 = 0;
          b = true;
        } else if (b) {
          b = false;
          i1 = q % cs.size(0) * cs.size(1) + q / cs.size(0);
        } else {
          ::coder::SizeType i2;
          i2 = cs.size(1) * cs.size(0) - 1;
          if (i1 > MAX_int32_T - cs.size(1)) {
            i1 = q % cs.size(0) * cs.size(1) + q / cs.size(0);
          } else {
            i1 += cs.size(1);
            if (i1 > i2) {
              i1 -= i2;
            }
          }
        }
        ::sfe_sfuncs::bar_7_sfunc(cs[i1], &f_N[0], &f_deriv[0]);
        sfvals[sfvals.size(1) * q] = f_N[0];
        sdvals[sdvals.size(2) * sdvals.size(1) * q] = f_deriv[0];
        sfvals[sfvals.size(1) * q + 1] = f_N[1];
        sdvals[sdvals.size(2) + sdvals.size(2) * sdvals.size(1) * q] =
            f_deriv[1];
        sfvals[sfvals.size(1) * q + 2] = f_N[2];
        sdvals[sdvals.size(2) * 2 + sdvals.size(2) * sdvals.size(1) * q] =
            f_deriv[2];
        sfvals[sfvals.size(1) * q + 3] = f_N[3];
        sdvals[sdvals.size(2) * 3 + sdvals.size(2) * sdvals.size(1) * q] =
            f_deriv[3];
        sfvals[sfvals.size(1) * q + 4] = f_N[4];
        sdvals[sdvals.size(2) * 4 + sdvals.size(2) * sdvals.size(1) * q] =
            f_deriv[4];
        sfvals[sfvals.size(1) * q + 5] = f_N[5];
        sdvals[sdvals.size(2) * 5 + sdvals.size(2) * sdvals.size(1) * q] =
            f_deriv[5];
        sfvals[sfvals.size(1) * q + 6] = f_N[6];
        sdvals[sdvals.size(2) * 6 + sdvals.size(2) * sdvals.size(1) * q] =
            f_deriv[6];
      }
    } break;
    }
  } else {
    //  GL
    switch (etype) {
    case 45: {
      ::coder::SizeType i1;
      boolean_T b;
      boolean_T b1;
      b = true;
      b1 = ((cs.size(1) <= 0) || (cs.size(0) <= 0));
      i = cs.size(1) * cs.size(0);
      i1 = 0;
      for (::coder::SizeType q{0}; q <= nqp; q++) {
        real_T c_N[4];
        real_T c_deriv[4];
        if (b1 || (q >= i)) {
          i1 = 0;
          b = true;
        } else if (b) {
          b = false;
          i1 = q % cs.size(0) * cs.size(1) + q / cs.size(0);
        } else {
          ::coder::SizeType i2;
          i2 = cs.size(1) * cs.size(0) - 1;
          if (i1 > MAX_int32_T - cs.size(1)) {
            i1 = q % cs.size(0) * cs.size(1) + q / cs.size(0);
          } else {
            i1 += cs.size(1);
            if (i1 > i2) {
              i1 -= i2;
            }
          }
        }
        ::sfe_sfuncs::bar_gl_4_sfunc(cs[i1], &c_N[0], &c_deriv[0]);
        sfvals[sfvals.size(1) * q] = c_N[0];
        sdvals[sdvals.size(2) * sdvals.size(1) * q] = c_deriv[0];
        sfvals[sfvals.size(1) * q + 1] = c_N[1];
        sdvals[sdvals.size(2) + sdvals.size(2) * sdvals.size(1) * q] =
            c_deriv[1];
        sfvals[sfvals.size(1) * q + 2] = c_N[2];
        sdvals[sdvals.size(2) * 2 + sdvals.size(2) * sdvals.size(1) * q] =
            c_deriv[2];
        sfvals[sfvals.size(1) * q + 3] = c_N[3];
        sdvals[sdvals.size(2) * 3 + sdvals.size(2) * sdvals.size(1) * q] =
            c_deriv[3];
      }
    } break;
    case 49: {
      ::coder::SizeType i1;
      boolean_T b;
      boolean_T b1;
      b = true;
      b1 = ((cs.size(1) <= 0) || (cs.size(0) <= 0));
      i = cs.size(1) * cs.size(0);
      i1 = 0;
      for (::coder::SizeType q{0}; q <= nqp; q++) {
        real_T d_N[5];
        real_T d_deriv[5];
        if (b1 || (q >= i)) {
          i1 = 0;
          b = true;
        } else if (b) {
          b = false;
          i1 = q % cs.size(0) * cs.size(1) + q / cs.size(0);
        } else {
          ::coder::SizeType i2;
          i2 = cs.size(1) * cs.size(0) - 1;
          if (i1 > MAX_int32_T - cs.size(1)) {
            i1 = q % cs.size(0) * cs.size(1) + q / cs.size(0);
          } else {
            i1 += cs.size(1);
            if (i1 > i2) {
              i1 -= i2;
            }
          }
        }
        ::sfe_sfuncs::bar_gl_5_sfunc(cs[i1], &d_N[0], &d_deriv[0]);
        sfvals[sfvals.size(1) * q] = d_N[0];
        sdvals[sdvals.size(2) * sdvals.size(1) * q] = d_deriv[0];
        sfvals[sfvals.size(1) * q + 1] = d_N[1];
        sdvals[sdvals.size(2) + sdvals.size(2) * sdvals.size(1) * q] =
            d_deriv[1];
        sfvals[sfvals.size(1) * q + 2] = d_N[2];
        sdvals[sdvals.size(2) * 2 + sdvals.size(2) * sdvals.size(1) * q] =
            d_deriv[2];
        sfvals[sfvals.size(1) * q + 3] = d_N[3];
        sdvals[sdvals.size(2) * 3 + sdvals.size(2) * sdvals.size(1) * q] =
            d_deriv[3];
        sfvals[sfvals.size(1) * q + 4] = d_N[4];
        sdvals[sdvals.size(2) * 4 + sdvals.size(2) * sdvals.size(1) * q] =
            d_deriv[4];
      }
    } break;
    case 53: {
      ::coder::SizeType i1;
      boolean_T b;
      boolean_T b1;
      b = true;
      b1 = ((cs.size(1) <= 0) || (cs.size(0) <= 0));
      i = cs.size(1) * cs.size(0);
      i1 = 0;
      for (::coder::SizeType q{0}; q <= nqp; q++) {
        real_T e_N[6];
        real_T e_deriv[6];
        if (b1 || (q >= i)) {
          i1 = 0;
          b = true;
        } else if (b) {
          b = false;
          i1 = q % cs.size(0) * cs.size(1) + q / cs.size(0);
        } else {
          ::coder::SizeType i2;
          i2 = cs.size(1) * cs.size(0) - 1;
          if (i1 > MAX_int32_T - cs.size(1)) {
            i1 = q % cs.size(0) * cs.size(1) + q / cs.size(0);
          } else {
            i1 += cs.size(1);
            if (i1 > i2) {
              i1 -= i2;
            }
          }
        }
        ::sfe_sfuncs::bar_gl_6_sfunc(cs[i1], &e_N[0], &e_deriv[0]);
        sfvals[sfvals.size(1) * q] = e_N[0];
        sdvals[sdvals.size(2) * sdvals.size(1) * q] = e_deriv[0];
        sfvals[sfvals.size(1) * q + 1] = e_N[1];
        sdvals[sdvals.size(2) + sdvals.size(2) * sdvals.size(1) * q] =
            e_deriv[1];
        sfvals[sfvals.size(1) * q + 2] = e_N[2];
        sdvals[sdvals.size(2) * 2 + sdvals.size(2) * sdvals.size(1) * q] =
            e_deriv[2];
        sfvals[sfvals.size(1) * q + 3] = e_N[3];
        sdvals[sdvals.size(2) * 3 + sdvals.size(2) * sdvals.size(1) * q] =
            e_deriv[3];
        sfvals[sfvals.size(1) * q + 4] = e_N[4];
        sdvals[sdvals.size(2) * 4 + sdvals.size(2) * sdvals.size(1) * q] =
            e_deriv[4];
        sfvals[sfvals.size(1) * q + 5] = e_N[5];
        sdvals[sdvals.size(2) * 5 + sdvals.size(2) * sdvals.size(1) * q] =
            e_deriv[5];
      }
    } break;
    default: {
      ::coder::SizeType i1;
      boolean_T b;
      boolean_T b1;
      m2cAssert(etype == 57, "Only support up to sextic.");
      b = true;
      b1 = ((cs.size(1) <= 0) || (cs.size(0) <= 0));
      i = cs.size(1) * cs.size(0);
      i1 = 0;
      for (::coder::SizeType q{0}; q <= nqp; q++) {
        real_T f_N[7];
        real_T f_deriv[7];
        if (b1 || (q >= i)) {
          i1 = 0;
          b = true;
        } else if (b) {
          b = false;
          i1 = q % cs.size(0) * cs.size(1) + q / cs.size(0);
        } else {
          ::coder::SizeType i2;
          i2 = cs.size(1) * cs.size(0) - 1;
          if (i1 > MAX_int32_T - cs.size(1)) {
            i1 = q % cs.size(0) * cs.size(1) + q / cs.size(0);
          } else {
            i1 += cs.size(1);
            if (i1 > i2) {
              i1 -= i2;
            }
          }
        }
        ::sfe_sfuncs::bar_gl_7_sfunc(cs[i1], &f_N[0], &f_deriv[0]);
        sfvals[sfvals.size(1) * q] = f_N[0];
        sdvals[sdvals.size(2) * sdvals.size(1) * q] = f_deriv[0];
        sfvals[sfvals.size(1) * q + 1] = f_N[1];
        sdvals[sdvals.size(2) + sdvals.size(2) * sdvals.size(1) * q] =
            f_deriv[1];
        sfvals[sfvals.size(1) * q + 2] = f_N[2];
        sdvals[sdvals.size(2) * 2 + sdvals.size(2) * sdvals.size(1) * q] =
            f_deriv[2];
        sfvals[sfvals.size(1) * q + 3] = f_N[3];
        sdvals[sdvals.size(2) * 3 + sdvals.size(2) * sdvals.size(1) * q] =
            f_deriv[3];
        sfvals[sfvals.size(1) * q + 4] = f_N[4];
        sdvals[sdvals.size(2) * 4 + sdvals.size(2) * sdvals.size(1) * q] =
            f_deriv[4];
        sfvals[sfvals.size(1) * q + 5] = f_N[5];
        sdvals[sdvals.size(2) * 5 + sdvals.size(2) * sdvals.size(1) * q] =
            f_deriv[5];
        sfvals[sfvals.size(1) * q + 6] = f_N[6];
        sdvals[sdvals.size(2) * 6 + sdvals.size(2) * sdvals.size(1) * q] =
            f_deriv[6];
      }
    } break;
    }
  }
}

static void sfe2_tabulate_equi_quad(::coder::SizeType etype,
                                    const ::coder::array<real_T, 2U> &cs,
                                    ::coder::array<real_T, 2U> &sfvals,
                                    ::coder::array<real_T, 3U> &sdvals)
{
  real_T b_tmp_data[1029];
  real_T tmp_data[1029];
  real_T dv10[98];
  real_T dv7[72];
  real_T dv6[50];
  real_T dv11[49];
  real_T dv9[36];
  real_T dv8[25];
  ::coder::SizeType b_tmp_size_idx_1;
  ::coder::SizeType b_tmp_size_idx_2;
  ::coder::SizeType i;
  ::coder::SizeType i4;
  ::coder::SizeType i5;
  ::coder::SizeType i6;
  ::coder::SizeType i7;
  ::coder::SizeType loop_ub;
  ::coder::SizeType nqp;
  int16_T b_unnamed_idx_1;
  int16_T b_unnamed_idx_2;
  //  triangular
  nqp = cs.size(0) - 1;
  i = iv[etype - 1];
  sfvals.set_size(cs.size(0), i);
  sdvals.set_size(cs.size(0), i, cs.size(1));
  switch (etype) {
  case 100: {
    for (::coder::SizeType q{0}; q <= nqp; q++) {
      real_T dv3[8];
      real_T dv[4];
      ::coder::SizeType i1;
      ::coder::SizeType i2;
      ::coder::SizeType tmp_size_idx_1;
      ::coder::SizeType tmp_size_idx_2;
      ::coder::SizeType ub_loop;
      int16_T unnamed_idx_1;
      int16_T unnamed_idx_2;
      ::sfe_sfuncs::quad_4_sfunc(cs[cs.size(1) * q], cs[cs.size(1) * q + 1],
                                 &dv[0], &dv3[0]);
      ub_loop = sfvals.size(1);
      for (i = 0; i < ub_loop; i++) {
        sfvals[i + sfvals.size(1) * q] = dv[i];
      }
      unnamed_idx_1 = static_cast<int16_T>(sdvals.size(1));
      unnamed_idx_2 = static_cast<int16_T>(sdvals.size(2));
      i = 0;
      i1 = 0;
      ub_loop = 0;
      i2 = 0;
      tmp_size_idx_2 = sdvals.size(2);
      tmp_size_idx_1 = sdvals.size(1);
      for (::coder::SizeType i3{0}; i3 < unnamed_idx_1 * unnamed_idx_2; i3++) {
        tmp_data[i1 + tmp_size_idx_2 * i] = dv3[i2 + (ub_loop << 1)];
        ub_loop++;
        i++;
        if (i > tmp_size_idx_1 - 1) {
          i = 0;
          i1++;
        }
        if (ub_loop > 3) {
          ub_loop = 0;
          i2++;
        }
      }
      for (i = 0; i < tmp_size_idx_1; i++) {
        for (i1 = 0; i1 < tmp_size_idx_2; i1++) {
          sdvals[(i1 + sdvals.size(2) * i) +
                 sdvals.size(2) * sdvals.size(1) * q] =
              tmp_data[i1 + tmp_size_idx_2 * i];
        }
      }
    }
  } break;
  case 104: {
    for (::coder::SizeType q{0}; q <= nqp; q++) {
      real_T dv4[18];
      real_T dv1[9];
      ::coder::SizeType i1;
      ::coder::SizeType i2;
      ::coder::SizeType tmp_size_idx_1;
      ::coder::SizeType tmp_size_idx_2;
      ::coder::SizeType ub_loop;
      int16_T unnamed_idx_1;
      int16_T unnamed_idx_2;
      ::sfe_sfuncs::quad_9_sfunc(cs[cs.size(1) * q], cs[cs.size(1) * q + 1],
                                 &dv1[0], &dv4[0]);
      ub_loop = sfvals.size(1);
      for (i = 0; i < ub_loop; i++) {
        sfvals[i + sfvals.size(1) * q] = dv1[i];
      }
      unnamed_idx_1 = static_cast<int16_T>(sdvals.size(1));
      unnamed_idx_2 = static_cast<int16_T>(sdvals.size(2));
      i = 0;
      i1 = 0;
      ub_loop = 0;
      i2 = 0;
      tmp_size_idx_2 = sdvals.size(2);
      tmp_size_idx_1 = sdvals.size(1);
      for (::coder::SizeType i3{0}; i3 < unnamed_idx_1 * unnamed_idx_2; i3++) {
        tmp_data[i1 + tmp_size_idx_2 * i] = dv4[i2 + (ub_loop << 1)];
        ub_loop++;
        i++;
        if (i > tmp_size_idx_1 - 1) {
          i = 0;
          i1++;
        }
        if (ub_loop > 8) {
          ub_loop = 0;
          i2++;
        }
      }
      for (i = 0; i < tmp_size_idx_1; i++) {
        for (i1 = 0; i1 < tmp_size_idx_2; i1++) {
          sdvals[(i1 + sdvals.size(2) * i) +
                 sdvals.size(2) * sdvals.size(1) * q] =
              tmp_data[i1 + tmp_size_idx_2 * i];
        }
      }
    }
  } break;
  case 108: {
    for (::coder::SizeType q{0}; q <= nqp; q++) {
      real_T dv5[32];
      real_T dv2[16];
      ::coder::SizeType i1;
      ::coder::SizeType i2;
      ::coder::SizeType tmp_size_idx_1;
      ::coder::SizeType tmp_size_idx_2;
      ::coder::SizeType ub_loop;
      int16_T unnamed_idx_1;
      int16_T unnamed_idx_2;
      ::sfe_sfuncs::quad_16_sfunc(cs[cs.size(1) * q], cs[cs.size(1) * q + 1],
                                  &dv2[0], &dv5[0]);
      ub_loop = sfvals.size(1);
      for (i = 0; i < ub_loop; i++) {
        sfvals[i + sfvals.size(1) * q] = dv2[i];
      }
      unnamed_idx_1 = static_cast<int16_T>(sdvals.size(1));
      unnamed_idx_2 = static_cast<int16_T>(sdvals.size(2));
      i = 0;
      i1 = 0;
      ub_loop = 0;
      i2 = 0;
      tmp_size_idx_2 = sdvals.size(2);
      tmp_size_idx_1 = sdvals.size(1);
      for (::coder::SizeType i3{0}; i3 < unnamed_idx_1 * unnamed_idx_2; i3++) {
        tmp_data[i1 + tmp_size_idx_2 * i] = dv5[i2 + (ub_loop << 1)];
        ub_loop++;
        i++;
        if (i > tmp_size_idx_1 - 1) {
          i = 0;
          i1++;
        }
        if (ub_loop > 15) {
          ub_loop = 0;
          i2++;
        }
      }
      for (i = 0; i < tmp_size_idx_1; i++) {
        for (i1 = 0; i1 < tmp_size_idx_2; i1++) {
          sdvals[(i1 + sdvals.size(2) * i) +
                 sdvals.size(2) * sdvals.size(1) * q] =
              tmp_data[i1 + tmp_size_idx_2 * i];
        }
      }
    }
  } break;
  case 112: {
    ::coder::SizeType ub_loop;
    //  Serial mode
    ub_loop = cs.size(0) - 1;
    for (::coder::SizeType b_q = 0; b_q <= ub_loop; b_q++) {
      quad_25(cs[cs.size(1) * b_q], cs[cs.size(1) * b_q + 1], dv8, dv6);
      loop_ub = sfvals.size(1);
      for (i4 = 0; i4 < loop_ub; i4++) {
        sfvals[i4 + sfvals.size(1) * b_q] = dv8[i4];
      }
      b_unnamed_idx_1 = static_cast<int16_T>(sdvals.size(1));
      b_unnamed_idx_2 = static_cast<int16_T>(sdvals.size(2));
      i4 = 0;
      i5 = 0;
      loop_ub = 0;
      i6 = 0;
      b_tmp_size_idx_2 = sdvals.size(2);
      b_tmp_size_idx_1 = sdvals.size(1);
      for (i7 = 0; i7 < b_unnamed_idx_1 * b_unnamed_idx_2; i7++) {
        b_tmp_data[i5 + b_tmp_size_idx_2 * i4] = dv6[i6 + (loop_ub << 1)];
        loop_ub++;
        i4++;
        if (i4 > b_tmp_size_idx_1 - 1) {
          i4 = 0;
          i5++;
        }
        if (loop_ub > 24) {
          loop_ub = 0;
          i6++;
        }
      }
      for (i4 = 0; i4 < b_tmp_size_idx_1; i4++) {
        for (i5 = 0; i5 < b_tmp_size_idx_2; i5++) {
          sdvals[(i5 + sdvals.size(2) * i4) +
                 sdvals.size(2) * sdvals.size(1) * b_q] =
              b_tmp_data[i5 + b_tmp_size_idx_2 * i4];
        }
      }
    }
  } break;
  case 116: {
    ::coder::SizeType ub_loop;
    //  Serial mode
    ub_loop = cs.size(0) - 1;
    for (::coder::SizeType b_q = 0; b_q <= ub_loop; b_q++) {
      quad_36(cs[cs.size(1) * b_q], cs[cs.size(1) * b_q + 1], dv9, dv7);
      loop_ub = sfvals.size(1);
      for (i4 = 0; i4 < loop_ub; i4++) {
        sfvals[i4 + sfvals.size(1) * b_q] = dv9[i4];
      }
      b_unnamed_idx_1 = static_cast<int16_T>(sdvals.size(1));
      b_unnamed_idx_2 = static_cast<int16_T>(sdvals.size(2));
      i4 = 0;
      i5 = 0;
      loop_ub = 0;
      i6 = 0;
      b_tmp_size_idx_2 = sdvals.size(2);
      b_tmp_size_idx_1 = sdvals.size(1);
      for (i7 = 0; i7 < b_unnamed_idx_1 * b_unnamed_idx_2; i7++) {
        b_tmp_data[i5 + b_tmp_size_idx_2 * i4] = dv7[i6 + (loop_ub << 1)];
        loop_ub++;
        i4++;
        if (i4 > b_tmp_size_idx_1 - 1) {
          i4 = 0;
          i5++;
        }
        if (loop_ub > 35) {
          loop_ub = 0;
          i6++;
        }
      }
      for (i4 = 0; i4 < b_tmp_size_idx_1; i4++) {
        for (i5 = 0; i5 < b_tmp_size_idx_2; i5++) {
          sdvals[(i5 + sdvals.size(2) * i4) +
                 sdvals.size(2) * sdvals.size(1) * b_q] =
              b_tmp_data[i5 + b_tmp_size_idx_2 * i4];
        }
      }
    }
  } break;
  default: {
    ::coder::SizeType ub_loop;
    m2cAssert(etype == 120, "Only supports up to sextic.");
    //  Serial mode
    ub_loop = cs.size(0) - 1;
    for (::coder::SizeType b_q = 0; b_q <= ub_loop; b_q++) {
      quad_49(cs[cs.size(1) * b_q], cs[cs.size(1) * b_q + 1], dv11, dv10);
      loop_ub = sfvals.size(1);
      for (i4 = 0; i4 < loop_ub; i4++) {
        sfvals[i4 + sfvals.size(1) * b_q] = dv11[i4];
      }
      b_unnamed_idx_1 = static_cast<int16_T>(sdvals.size(1));
      b_unnamed_idx_2 = static_cast<int16_T>(sdvals.size(2));
      i4 = 0;
      i5 = 0;
      loop_ub = 0;
      i6 = 0;
      b_tmp_size_idx_2 = sdvals.size(2);
      b_tmp_size_idx_1 = sdvals.size(1);
      for (i7 = 0; i7 < b_unnamed_idx_1 * b_unnamed_idx_2; i7++) {
        b_tmp_data[i5 + b_tmp_size_idx_2 * i4] = dv10[i6 + (loop_ub << 1)];
        loop_ub++;
        i4++;
        if (i4 > b_tmp_size_idx_1 - 1) {
          i4 = 0;
          i5++;
        }
        if (loop_ub > 48) {
          loop_ub = 0;
          i6++;
        }
      }
      for (i4 = 0; i4 < b_tmp_size_idx_1; i4++) {
        for (i5 = 0; i5 < b_tmp_size_idx_2; i5++) {
          sdvals[(i5 + sdvals.size(2) * i4) +
                 sdvals.size(2) * sdvals.size(1) * b_q] =
              b_tmp_data[i5 + b_tmp_size_idx_2 * i4];
        }
      }
    }
  } break;
  }
}

static void sfe2_tabulate_equi_tri(::coder::SizeType etype,
                                   const ::coder::array<real_T, 2U> &cs,
                                   ::coder::array<real_T, 2U> &sfvals,
                                   ::coder::array<real_T, 3U> &sdvals)
{
  real_T b_tmp_data[1029];
  real_T tmp_data[1029];
  real_T dv9[56];
  real_T dv7[42];
  real_T dv10[28];
  real_T dv8[21];
  ::coder::SizeType b_tmp_size_idx_1;
  ::coder::SizeType b_tmp_size_idx_2;
  ::coder::SizeType i;
  ::coder::SizeType i3;
  ::coder::SizeType i5;
  ::coder::SizeType i6;
  ::coder::SizeType i7;
  ::coder::SizeType loop_ub;
  ::coder::SizeType nqp;
  int16_T b_unnamed_idx_1;
  int16_T b_unnamed_idx_2;
  //  triangular
  nqp = cs.size(0) - 1;
  i = iv[etype - 1];
  sfvals.set_size(cs.size(0), i);
  sdvals.set_size(cs.size(0), i, cs.size(1));
  switch (etype) {
  case 68: {
    for (::coder::SizeType q{0}; q <= nqp; q++) {
      real_T dv1[6];
      real_T dv[3];
      ::coder::SizeType i1;
      ::coder::SizeType i2;
      ::coder::SizeType tmp_size_idx_1;
      ::coder::SizeType tmp_size_idx_2;
      ::coder::SizeType ub_loop;
      int16_T unnamed_idx_1;
      int16_T unnamed_idx_2;
      ::sfe_sfuncs::tri_3_sfunc(cs[cs.size(1) * q], cs[cs.size(1) * q + 1],
                                &dv[0], &dv1[0]);
      ub_loop = sfvals.size(1);
      for (i = 0; i < ub_loop; i++) {
        sfvals[i + sfvals.size(1) * q] = dv[i];
      }
      unnamed_idx_1 = static_cast<int16_T>(sdvals.size(1));
      unnamed_idx_2 = static_cast<int16_T>(sdvals.size(2));
      i = 0;
      i1 = 0;
      ub_loop = 0;
      i2 = 0;
      tmp_size_idx_2 = sdvals.size(2);
      tmp_size_idx_1 = sdvals.size(1);
      for (::coder::SizeType i4{0}; i4 < unnamed_idx_1 * unnamed_idx_2; i4++) {
        tmp_data[i1 + tmp_size_idx_2 * i] = dv1[i2 + (ub_loop << 1)];
        ub_loop++;
        i++;
        if (i > tmp_size_idx_1 - 1) {
          i = 0;
          i1++;
        }
        if (ub_loop > 2) {
          ub_loop = 0;
          i2++;
        }
      }
      for (i = 0; i < tmp_size_idx_1; i++) {
        for (i1 = 0; i1 < tmp_size_idx_2; i1++) {
          sdvals[(i1 + sdvals.size(2) * i) +
                 sdvals.size(2) * sdvals.size(1) * q] =
              tmp_data[i1 + tmp_size_idx_2 * i];
        }
      }
    }
  } break;
  case 72: {
    for (::coder::SizeType q{0}; q <= nqp; q++) {
      real_T dv4[12];
      real_T dv1[6];
      ::coder::SizeType i1;
      ::coder::SizeType i2;
      ::coder::SizeType tmp_size_idx_1;
      ::coder::SizeType tmp_size_idx_2;
      ::coder::SizeType ub_loop;
      int16_T unnamed_idx_1;
      int16_T unnamed_idx_2;
      ::sfe_sfuncs::tri_6_sfunc(cs[cs.size(1) * q], cs[cs.size(1) * q + 1],
                                &dv1[0], &dv4[0]);
      ub_loop = sfvals.size(1);
      for (i = 0; i < ub_loop; i++) {
        sfvals[i + sfvals.size(1) * q] = dv1[i];
      }
      unnamed_idx_1 = static_cast<int16_T>(sdvals.size(1));
      unnamed_idx_2 = static_cast<int16_T>(sdvals.size(2));
      i = 0;
      i1 = 0;
      ub_loop = 0;
      i2 = 0;
      tmp_size_idx_2 = sdvals.size(2);
      tmp_size_idx_1 = sdvals.size(1);
      for (::coder::SizeType i4{0}; i4 < unnamed_idx_1 * unnamed_idx_2; i4++) {
        tmp_data[i1 + tmp_size_idx_2 * i] = dv4[i2 + (ub_loop << 1)];
        ub_loop++;
        i++;
        if (i > tmp_size_idx_1 - 1) {
          i = 0;
          i1++;
        }
        if (ub_loop > 5) {
          ub_loop = 0;
          i2++;
        }
      }
      for (i = 0; i < tmp_size_idx_1; i++) {
        for (i1 = 0; i1 < tmp_size_idx_2; i1++) {
          sdvals[(i1 + sdvals.size(2) * i) +
                 sdvals.size(2) * sdvals.size(1) * q] =
              tmp_data[i1 + tmp_size_idx_2 * i];
        }
      }
    }
  } break;
  case 76: {
    for (::coder::SizeType q{0}; q <= nqp; q++) {
      real_T dv5[20];
      real_T dv2[10];
      ::coder::SizeType i1;
      ::coder::SizeType i2;
      ::coder::SizeType tmp_size_idx_1;
      ::coder::SizeType tmp_size_idx_2;
      ::coder::SizeType ub_loop;
      int16_T unnamed_idx_1;
      int16_T unnamed_idx_2;
      ::sfe_sfuncs::tri_10_sfunc(cs[cs.size(1) * q], cs[cs.size(1) * q + 1],
                                 &dv2[0], &dv5[0]);
      ub_loop = sfvals.size(1);
      for (i = 0; i < ub_loop; i++) {
        sfvals[i + sfvals.size(1) * q] = dv2[i];
      }
      unnamed_idx_1 = static_cast<int16_T>(sdvals.size(1));
      unnamed_idx_2 = static_cast<int16_T>(sdvals.size(2));
      i = 0;
      i1 = 0;
      ub_loop = 0;
      i2 = 0;
      tmp_size_idx_2 = sdvals.size(2);
      tmp_size_idx_1 = sdvals.size(1);
      for (::coder::SizeType i4{0}; i4 < unnamed_idx_1 * unnamed_idx_2; i4++) {
        tmp_data[i1 + tmp_size_idx_2 * i] = dv5[i2 + (ub_loop << 1)];
        ub_loop++;
        i++;
        if (i > tmp_size_idx_1 - 1) {
          i = 0;
          i1++;
        }
        if (ub_loop > 9) {
          ub_loop = 0;
          i2++;
        }
      }
      for (i = 0; i < tmp_size_idx_1; i++) {
        for (i1 = 0; i1 < tmp_size_idx_2; i1++) {
          sdvals[(i1 + sdvals.size(2) * i) +
                 sdvals.size(2) * sdvals.size(1) * q] =
              tmp_data[i1 + tmp_size_idx_2 * i];
        }
      }
    }
  } break;
  case 80: {
    for (::coder::SizeType q{0}; q <= nqp; q++) {
      real_T dv6[30];
      real_T dv3[15];
      ::coder::SizeType i1;
      ::coder::SizeType i2;
      ::coder::SizeType tmp_size_idx_1;
      ::coder::SizeType tmp_size_idx_2;
      ::coder::SizeType ub_loop;
      int16_T unnamed_idx_1;
      int16_T unnamed_idx_2;
      ::sfe_sfuncs::tri_15_sfunc(cs[cs.size(1) * q], cs[cs.size(1) * q + 1],
                                 &dv3[0], &dv6[0]);
      ub_loop = sfvals.size(1);
      for (i = 0; i < ub_loop; i++) {
        sfvals[i + sfvals.size(1) * q] = dv3[i];
      }
      unnamed_idx_1 = static_cast<int16_T>(sdvals.size(1));
      unnamed_idx_2 = static_cast<int16_T>(sdvals.size(2));
      i = 0;
      i1 = 0;
      ub_loop = 0;
      i2 = 0;
      tmp_size_idx_2 = sdvals.size(2);
      tmp_size_idx_1 = sdvals.size(1);
      for (::coder::SizeType i4{0}; i4 < unnamed_idx_1 * unnamed_idx_2; i4++) {
        tmp_data[i1 + tmp_size_idx_2 * i] = dv6[i2 + (ub_loop << 1)];
        ub_loop++;
        i++;
        if (i > tmp_size_idx_1 - 1) {
          i = 0;
          i1++;
        }
        if (ub_loop > 14) {
          ub_loop = 0;
          i2++;
        }
      }
      for (i = 0; i < tmp_size_idx_1; i++) {
        for (i1 = 0; i1 < tmp_size_idx_2; i1++) {
          sdvals[(i1 + sdvals.size(2) * i) +
                 sdvals.size(2) * sdvals.size(1) * q] =
              tmp_data[i1 + tmp_size_idx_2 * i];
        }
      }
    }
  } break;
  case 84: {
    ::coder::SizeType ub_loop;
    //  Serial mode
    ub_loop = cs.size(0) - 1;
    for (::coder::SizeType b_q = 0; b_q <= ub_loop; b_q++) {
      tri_21(cs[cs.size(1) * b_q], cs[cs.size(1) * b_q + 1], dv8, dv7);
      loop_ub = sfvals.size(1);
      for (i3 = 0; i3 < loop_ub; i3++) {
        sfvals[i3 + sfvals.size(1) * b_q] = dv8[i3];
      }
      b_unnamed_idx_1 = static_cast<int16_T>(sdvals.size(1));
      b_unnamed_idx_2 = static_cast<int16_T>(sdvals.size(2));
      i3 = 0;
      i5 = 0;
      loop_ub = 0;
      i6 = 0;
      b_tmp_size_idx_2 = sdvals.size(2);
      b_tmp_size_idx_1 = sdvals.size(1);
      for (i7 = 0; i7 < b_unnamed_idx_1 * b_unnamed_idx_2; i7++) {
        b_tmp_data[i5 + b_tmp_size_idx_2 * i3] = dv7[i6 + (loop_ub << 1)];
        loop_ub++;
        i3++;
        if (i3 > b_tmp_size_idx_1 - 1) {
          i3 = 0;
          i5++;
        }
        if (loop_ub > 20) {
          loop_ub = 0;
          i6++;
        }
      }
      for (i3 = 0; i3 < b_tmp_size_idx_1; i3++) {
        for (i5 = 0; i5 < b_tmp_size_idx_2; i5++) {
          sdvals[(i5 + sdvals.size(2) * i3) +
                 sdvals.size(2) * sdvals.size(1) * b_q] =
              b_tmp_data[i5 + b_tmp_size_idx_2 * i3];
        }
      }
    }
  } break;
  default: {
    ::coder::SizeType ub_loop;
    m2cAssert(etype == 88, "Only support up to sextic.");
    //  Serial mode
    ub_loop = cs.size(0) - 1;
    for (::coder::SizeType b_q = 0; b_q <= ub_loop; b_q++) {
      tri_28(cs[cs.size(1) * b_q], cs[cs.size(1) * b_q + 1], dv10, dv9);
      loop_ub = sfvals.size(1);
      for (i3 = 0; i3 < loop_ub; i3++) {
        sfvals[i3 + sfvals.size(1) * b_q] = dv10[i3];
      }
      b_unnamed_idx_1 = static_cast<int16_T>(sdvals.size(1));
      b_unnamed_idx_2 = static_cast<int16_T>(sdvals.size(2));
      i3 = 0;
      i5 = 0;
      loop_ub = 0;
      i6 = 0;
      b_tmp_size_idx_2 = sdvals.size(2);
      b_tmp_size_idx_1 = sdvals.size(1);
      for (i7 = 0; i7 < b_unnamed_idx_1 * b_unnamed_idx_2; i7++) {
        b_tmp_data[i5 + b_tmp_size_idx_2 * i3] = dv9[i6 + (loop_ub << 1)];
        loop_ub++;
        i3++;
        if (i3 > b_tmp_size_idx_1 - 1) {
          i3 = 0;
          i5++;
        }
        if (loop_ub > 27) {
          loop_ub = 0;
          i6++;
        }
      }
      for (i3 = 0; i3 < b_tmp_size_idx_1; i3++) {
        for (i5 = 0; i5 < b_tmp_size_idx_2; i5++) {
          sdvals[(i5 + sdvals.size(2) * i3) +
                 sdvals.size(2) * sdvals.size(1) * b_q] =
              b_tmp_data[i5 + b_tmp_size_idx_2 * i3];
        }
      }
    }
  } break;
  }
}

static void sfe2_tabulate_fek_tri(::coder::SizeType etype,
                                  const ::coder::array<real_T, 2U> &cs,
                                  ::coder::array<real_T, 2U> &sfvals,
                                  ::coder::array<real_T, 3U> &sdvals)
{
  real_T tmp_data[1029];
  real_T dv4[56];
  real_T dv1[42];
  real_T dv[30];
  real_T dv5[28];
  real_T dv3[21];
  real_T dv2[15];
  ::coder::SizeType i;
  ::coder::SizeType i1;
  ::coder::SizeType i2;
  ::coder::SizeType i3;
  ::coder::SizeType loop_ub;
  ::coder::SizeType tmp_size_idx_1;
  ::coder::SizeType tmp_size_idx_2;
  ::coder::SizeType ub_loop;
  int16_T unnamed_idx_1;
  int16_T unnamed_idx_2;
  //  triangular
  ub_loop = iv[etype - 1];
  sfvals.set_size(cs.size(0), ub_loop);
  sdvals.set_size(cs.size(0), ub_loop, cs.size(1));
  switch (etype) {
  case 82:
    //  Serial mode
    ub_loop = cs.size(0) - 1;
    for (::coder::SizeType q = 0; q <= ub_loop; q++) {
      tri_fek_15(cs[cs.size(1) * q], cs[cs.size(1) * q + 1], dv2, dv);
      loop_ub = sfvals.size(1);
      for (i = 0; i < loop_ub; i++) {
        sfvals[i + sfvals.size(1) * q] = dv2[i];
      }
      unnamed_idx_1 = static_cast<int16_T>(sdvals.size(1));
      unnamed_idx_2 = static_cast<int16_T>(sdvals.size(2));
      i = 0;
      i1 = 0;
      loop_ub = 0;
      i2 = 0;
      tmp_size_idx_2 = sdvals.size(2);
      tmp_size_idx_1 = sdvals.size(1);
      for (i3 = 0; i3 < unnamed_idx_1 * unnamed_idx_2; i3++) {
        tmp_data[i1 + tmp_size_idx_2 * i] = dv[i2 + (loop_ub << 1)];
        loop_ub++;
        i++;
        if (i > tmp_size_idx_1 - 1) {
          i = 0;
          i1++;
        }
        if (loop_ub > 14) {
          loop_ub = 0;
          i2++;
        }
      }
      for (i = 0; i < tmp_size_idx_1; i++) {
        for (i1 = 0; i1 < tmp_size_idx_2; i1++) {
          sdvals[(i1 + sdvals.size(2) * i) +
                 sdvals.size(2) * sdvals.size(1) * q] =
              tmp_data[i1 + tmp_size_idx_2 * i];
        }
      }
    }
    break;
  case 86:
    //  Serial mode
    ub_loop = cs.size(0) - 1;
    for (::coder::SizeType q = 0; q <= ub_loop; q++) {
      tri_fek_21(cs[cs.size(1) * q], cs[cs.size(1) * q + 1], dv3, dv1);
      loop_ub = sfvals.size(1);
      for (i = 0; i < loop_ub; i++) {
        sfvals[i + sfvals.size(1) * q] = dv3[i];
      }
      unnamed_idx_1 = static_cast<int16_T>(sdvals.size(1));
      unnamed_idx_2 = static_cast<int16_T>(sdvals.size(2));
      i = 0;
      i1 = 0;
      loop_ub = 0;
      i2 = 0;
      tmp_size_idx_2 = sdvals.size(2);
      tmp_size_idx_1 = sdvals.size(1);
      for (i3 = 0; i3 < unnamed_idx_1 * unnamed_idx_2; i3++) {
        tmp_data[i1 + tmp_size_idx_2 * i] = dv1[i2 + (loop_ub << 1)];
        loop_ub++;
        i++;
        if (i > tmp_size_idx_1 - 1) {
          i = 0;
          i1++;
        }
        if (loop_ub > 20) {
          loop_ub = 0;
          i2++;
        }
      }
      for (i = 0; i < tmp_size_idx_1; i++) {
        for (i1 = 0; i1 < tmp_size_idx_2; i1++) {
          sdvals[(i1 + sdvals.size(2) * i) +
                 sdvals.size(2) * sdvals.size(1) * q] =
              tmp_data[i1 + tmp_size_idx_2 * i];
        }
      }
    }
    break;
  default:
    m2cAssert(etype == 90, "Only supports up to sextic.");
    //  Serial mode
    ub_loop = cs.size(0) - 1;
    for (::coder::SizeType q = 0; q <= ub_loop; q++) {
      tri_fek_28(cs[cs.size(1) * q], cs[cs.size(1) * q + 1], dv5, dv4);
      loop_ub = sfvals.size(1);
      for (i = 0; i < loop_ub; i++) {
        sfvals[i + sfvals.size(1) * q] = dv5[i];
      }
      unnamed_idx_1 = static_cast<int16_T>(sdvals.size(1));
      unnamed_idx_2 = static_cast<int16_T>(sdvals.size(2));
      i = 0;
      i1 = 0;
      loop_ub = 0;
      i2 = 0;
      tmp_size_idx_2 = sdvals.size(2);
      tmp_size_idx_1 = sdvals.size(1);
      for (i3 = 0; i3 < unnamed_idx_1 * unnamed_idx_2; i3++) {
        tmp_data[i1 + tmp_size_idx_2 * i] = dv4[i2 + (loop_ub << 1)];
        loop_ub++;
        i++;
        if (i > tmp_size_idx_1 - 1) {
          i = 0;
          i1++;
        }
        if (loop_ub > 27) {
          loop_ub = 0;
          i2++;
        }
      }
      for (i = 0; i < tmp_size_idx_1; i++) {
        for (i1 = 0; i1 < tmp_size_idx_2; i1++) {
          sdvals[(i1 + sdvals.size(2) * i) +
                 sdvals.size(2) * sdvals.size(1) * q] =
              tmp_data[i1 + tmp_size_idx_2 * i];
        }
      }
    }
    break;
  }
}

static void sfe2_tabulate_gl_quad(::coder::SizeType etype,
                                  const ::coder::array<real_T, 2U> &cs,
                                  ::coder::array<real_T, 2U> &sfvals,
                                  ::coder::array<real_T, 3U> &sdvals)
{
  real_T b_tmp_data[1029];
  real_T tmp_data[1029];
  real_T dv6[98];
  real_T dv3[72];
  real_T dv2[50];
  real_T dv7[49];
  real_T dv5[36];
  real_T dv4[25];
  ::coder::SizeType b_tmp_size_idx_1;
  ::coder::SizeType b_tmp_size_idx_2;
  ::coder::SizeType i;
  ::coder::SizeType i4;
  ::coder::SizeType i5;
  ::coder::SizeType i6;
  ::coder::SizeType i7;
  ::coder::SizeType loop_ub;
  ::coder::SizeType nqp;
  int16_T b_unnamed_idx_1;
  int16_T b_unnamed_idx_2;
  //  quad
  nqp = cs.size(0);
  i = iv[etype - 1];
  sfvals.set_size(cs.size(0), i);
  sdvals.set_size(cs.size(0), i, cs.size(1));
  switch (etype) {
  case 109: {
    for (::coder::SizeType q{0}; q < nqp; q++) {
      real_T dv1[32];
      real_T dv[16];
      ::coder::SizeType i1;
      ::coder::SizeType i2;
      ::coder::SizeType tmp_size_idx_1;
      ::coder::SizeType tmp_size_idx_2;
      ::coder::SizeType ub_loop;
      int16_T unnamed_idx_1;
      int16_T unnamed_idx_2;
      ::sfe_sfuncs::quad_gl_16_sfunc(cs[cs.size(1) * q], cs[cs.size(1) * q + 1],
                                     &dv[0], &dv1[0]);
      ub_loop = sfvals.size(1);
      for (i = 0; i < ub_loop; i++) {
        sfvals[i + sfvals.size(1) * q] = dv[i];
      }
      unnamed_idx_1 = static_cast<int16_T>(sdvals.size(1));
      unnamed_idx_2 = static_cast<int16_T>(sdvals.size(2));
      i = 0;
      i1 = 0;
      ub_loop = 0;
      i2 = 0;
      tmp_size_idx_2 = sdvals.size(2);
      tmp_size_idx_1 = sdvals.size(1);
      for (::coder::SizeType i3{0}; i3 < unnamed_idx_1 * unnamed_idx_2; i3++) {
        tmp_data[i1 + tmp_size_idx_2 * i] = dv1[i2 + (ub_loop << 1)];
        ub_loop++;
        i++;
        if (i > tmp_size_idx_1 - 1) {
          i = 0;
          i1++;
        }
        if (ub_loop > 15) {
          ub_loop = 0;
          i2++;
        }
      }
      for (i = 0; i < tmp_size_idx_1; i++) {
        for (i1 = 0; i1 < tmp_size_idx_2; i1++) {
          sdvals[(i1 + sdvals.size(2) * i) +
                 sdvals.size(2) * sdvals.size(1) * q] =
              tmp_data[i1 + tmp_size_idx_2 * i];
        }
      }
    }
  } break;
  case 113: {
    ::coder::SizeType ub_loop;
    //  Serial mode
    ub_loop = cs.size(0) - 1;
    for (::coder::SizeType b_q = 0; b_q <= ub_loop; b_q++) {
      quad_gl_25(cs[cs.size(1) * b_q], cs[cs.size(1) * b_q + 1], dv4, dv2);
      loop_ub = sfvals.size(1);
      for (i4 = 0; i4 < loop_ub; i4++) {
        sfvals[i4 + sfvals.size(1) * b_q] = dv4[i4];
      }
      b_unnamed_idx_1 = static_cast<int16_T>(sdvals.size(1));
      b_unnamed_idx_2 = static_cast<int16_T>(sdvals.size(2));
      i4 = 0;
      i5 = 0;
      loop_ub = 0;
      i6 = 0;
      b_tmp_size_idx_2 = sdvals.size(2);
      b_tmp_size_idx_1 = sdvals.size(1);
      for (i7 = 0; i7 < b_unnamed_idx_1 * b_unnamed_idx_2; i7++) {
        b_tmp_data[i5 + b_tmp_size_idx_2 * i4] = dv2[i6 + (loop_ub << 1)];
        loop_ub++;
        i4++;
        if (i4 > b_tmp_size_idx_1 - 1) {
          i4 = 0;
          i5++;
        }
        if (loop_ub > 24) {
          loop_ub = 0;
          i6++;
        }
      }
      for (i4 = 0; i4 < b_tmp_size_idx_1; i4++) {
        for (i5 = 0; i5 < b_tmp_size_idx_2; i5++) {
          sdvals[(i5 + sdvals.size(2) * i4) +
                 sdvals.size(2) * sdvals.size(1) * b_q] =
              b_tmp_data[i5 + b_tmp_size_idx_2 * i4];
        }
      }
    }
  } break;
  case 117: {
    ::coder::SizeType ub_loop;
    //  Serial mode
    ub_loop = cs.size(0) - 1;
    for (::coder::SizeType b_q = 0; b_q <= ub_loop; b_q++) {
      quad_gl_36(cs[cs.size(1) * b_q], cs[cs.size(1) * b_q + 1], dv5, dv3);
      loop_ub = sfvals.size(1);
      for (i4 = 0; i4 < loop_ub; i4++) {
        sfvals[i4 + sfvals.size(1) * b_q] = dv5[i4];
      }
      b_unnamed_idx_1 = static_cast<int16_T>(sdvals.size(1));
      b_unnamed_idx_2 = static_cast<int16_T>(sdvals.size(2));
      i4 = 0;
      i5 = 0;
      loop_ub = 0;
      i6 = 0;
      b_tmp_size_idx_2 = sdvals.size(2);
      b_tmp_size_idx_1 = sdvals.size(1);
      for (i7 = 0; i7 < b_unnamed_idx_1 * b_unnamed_idx_2; i7++) {
        b_tmp_data[i5 + b_tmp_size_idx_2 * i4] = dv3[i6 + (loop_ub << 1)];
        loop_ub++;
        i4++;
        if (i4 > b_tmp_size_idx_1 - 1) {
          i4 = 0;
          i5++;
        }
        if (loop_ub > 35) {
          loop_ub = 0;
          i6++;
        }
      }
      for (i4 = 0; i4 < b_tmp_size_idx_1; i4++) {
        for (i5 = 0; i5 < b_tmp_size_idx_2; i5++) {
          sdvals[(i5 + sdvals.size(2) * i4) +
                 sdvals.size(2) * sdvals.size(1) * b_q] =
              b_tmp_data[i5 + b_tmp_size_idx_2 * i4];
        }
      }
    }
  } break;
  default: {
    ::coder::SizeType ub_loop;
    m2cAssert(etype == 121, "Only supports up to sextic.");
    //  Serial mode
    ub_loop = cs.size(0) - 1;
    for (::coder::SizeType b_q = 0; b_q <= ub_loop; b_q++) {
      quad_gl_49(cs[cs.size(1) * b_q], cs[cs.size(1) * b_q + 1], dv7, dv6);
      loop_ub = sfvals.size(1);
      for (i4 = 0; i4 < loop_ub; i4++) {
        sfvals[i4 + sfvals.size(1) * b_q] = dv7[i4];
      }
      b_unnamed_idx_1 = static_cast<int16_T>(sdvals.size(1));
      b_unnamed_idx_2 = static_cast<int16_T>(sdvals.size(2));
      i4 = 0;
      i5 = 0;
      loop_ub = 0;
      i6 = 0;
      b_tmp_size_idx_2 = sdvals.size(2);
      b_tmp_size_idx_1 = sdvals.size(1);
      for (i7 = 0; i7 < b_unnamed_idx_1 * b_unnamed_idx_2; i7++) {
        b_tmp_data[i5 + b_tmp_size_idx_2 * i4] = dv6[i6 + (loop_ub << 1)];
        loop_ub++;
        i4++;
        if (i4 > b_tmp_size_idx_1 - 1) {
          i4 = 0;
          i5++;
        }
        if (loop_ub > 48) {
          loop_ub = 0;
          i6++;
        }
      }
      for (i4 = 0; i4 < b_tmp_size_idx_1; i4++) {
        for (i5 = 0; i5 < b_tmp_size_idx_2; i5++) {
          sdvals[(i5 + sdvals.size(2) * i4) +
                 sdvals.size(2) * sdvals.size(1) * b_q] =
              b_tmp_data[i5 + b_tmp_size_idx_2 * i4];
        }
      }
    }
  } break;
  }
}

static void sfe2_tabulate_gl_tri(::coder::SizeType etype,
                                 const ::coder::array<real_T, 2U> &cs,
                                 ::coder::array<real_T, 2U> &sfvals,
                                 ::coder::array<real_T, 3U> &sdvals)
{
  real_T b_tmp_data[1029];
  real_T tmp_data[1029];
  real_T dv6[56];
  real_T dv4[42];
  real_T dv7[28];
  real_T dv5[21];
  ::coder::SizeType b_tmp_size_idx_1;
  ::coder::SizeType b_tmp_size_idx_2;
  ::coder::SizeType i;
  ::coder::SizeType i3;
  ::coder::SizeType i5;
  ::coder::SizeType i6;
  ::coder::SizeType i7;
  ::coder::SizeType loop_ub;
  ::coder::SizeType nqp;
  int16_T b_unnamed_idx_1;
  int16_T b_unnamed_idx_2;
  //  triangular
  nqp = cs.size(0) - 1;
  i = iv[etype - 1];
  sfvals.set_size(cs.size(0), i);
  sdvals.set_size(cs.size(0), i, cs.size(1));
  switch (etype) {
  case 77: {
    for (::coder::SizeType q{0}; q <= nqp; q++) {
      real_T dv2[20];
      real_T dv[10];
      ::coder::SizeType i1;
      ::coder::SizeType i2;
      ::coder::SizeType tmp_size_idx_1;
      ::coder::SizeType tmp_size_idx_2;
      ::coder::SizeType ub_loop;
      int16_T unnamed_idx_1;
      int16_T unnamed_idx_2;
      ::sfe_sfuncs::tri_gl_10_sfunc(cs[cs.size(1) * q], cs[cs.size(1) * q + 1],
                                    &dv[0], &dv2[0]);
      ub_loop = sfvals.size(1);
      for (i = 0; i < ub_loop; i++) {
        sfvals[i + sfvals.size(1) * q] = dv[i];
      }
      unnamed_idx_1 = static_cast<int16_T>(sdvals.size(1));
      unnamed_idx_2 = static_cast<int16_T>(sdvals.size(2));
      i = 0;
      i1 = 0;
      ub_loop = 0;
      i2 = 0;
      tmp_size_idx_2 = sdvals.size(2);
      tmp_size_idx_1 = sdvals.size(1);
      for (::coder::SizeType i4{0}; i4 < unnamed_idx_1 * unnamed_idx_2; i4++) {
        tmp_data[i1 + tmp_size_idx_2 * i] = dv2[i2 + (ub_loop << 1)];
        ub_loop++;
        i++;
        if (i > tmp_size_idx_1 - 1) {
          i = 0;
          i1++;
        }
        if (ub_loop > 9) {
          ub_loop = 0;
          i2++;
        }
      }
      for (i = 0; i < tmp_size_idx_1; i++) {
        for (i1 = 0; i1 < tmp_size_idx_2; i1++) {
          sdvals[(i1 + sdvals.size(2) * i) +
                 sdvals.size(2) * sdvals.size(1) * q] =
              tmp_data[i1 + tmp_size_idx_2 * i];
        }
      }
    }
  } break;
  case 81: {
    for (::coder::SizeType q{0}; q <= nqp; q++) {
      real_T dv3[30];
      real_T dv1[15];
      ::coder::SizeType i1;
      ::coder::SizeType i2;
      ::coder::SizeType tmp_size_idx_1;
      ::coder::SizeType tmp_size_idx_2;
      ::coder::SizeType ub_loop;
      int16_T unnamed_idx_1;
      int16_T unnamed_idx_2;
      ::sfe_sfuncs::tri_gl_15_sfunc(cs[cs.size(1) * q], cs[cs.size(1) * q + 1],
                                    &dv1[0], &dv3[0]);
      ub_loop = sfvals.size(1);
      for (i = 0; i < ub_loop; i++) {
        sfvals[i + sfvals.size(1) * q] = dv1[i];
      }
      unnamed_idx_1 = static_cast<int16_T>(sdvals.size(1));
      unnamed_idx_2 = static_cast<int16_T>(sdvals.size(2));
      i = 0;
      i1 = 0;
      ub_loop = 0;
      i2 = 0;
      tmp_size_idx_2 = sdvals.size(2);
      tmp_size_idx_1 = sdvals.size(1);
      for (::coder::SizeType i4{0}; i4 < unnamed_idx_1 * unnamed_idx_2; i4++) {
        tmp_data[i1 + tmp_size_idx_2 * i] = dv3[i2 + (ub_loop << 1)];
        ub_loop++;
        i++;
        if (i > tmp_size_idx_1 - 1) {
          i = 0;
          i1++;
        }
        if (ub_loop > 14) {
          ub_loop = 0;
          i2++;
        }
      }
      for (i = 0; i < tmp_size_idx_1; i++) {
        for (i1 = 0; i1 < tmp_size_idx_2; i1++) {
          sdvals[(i1 + sdvals.size(2) * i) +
                 sdvals.size(2) * sdvals.size(1) * q] =
              tmp_data[i1 + tmp_size_idx_2 * i];
        }
      }
    }
  } break;
  case 85: {
    ::coder::SizeType ub_loop;
    //  Serial mode
    ub_loop = cs.size(0) - 1;
    for (::coder::SizeType b_q = 0; b_q <= ub_loop; b_q++) {
      tri_gl_21(cs[cs.size(1) * b_q], cs[cs.size(1) * b_q + 1], dv5, dv4);
      loop_ub = sfvals.size(1);
      for (i3 = 0; i3 < loop_ub; i3++) {
        sfvals[i3 + sfvals.size(1) * b_q] = dv5[i3];
      }
      b_unnamed_idx_1 = static_cast<int16_T>(sdvals.size(1));
      b_unnamed_idx_2 = static_cast<int16_T>(sdvals.size(2));
      i3 = 0;
      i5 = 0;
      loop_ub = 0;
      i6 = 0;
      b_tmp_size_idx_2 = sdvals.size(2);
      b_tmp_size_idx_1 = sdvals.size(1);
      for (i7 = 0; i7 < b_unnamed_idx_1 * b_unnamed_idx_2; i7++) {
        b_tmp_data[i5 + b_tmp_size_idx_2 * i3] = dv4[i6 + (loop_ub << 1)];
        loop_ub++;
        i3++;
        if (i3 > b_tmp_size_idx_1 - 1) {
          i3 = 0;
          i5++;
        }
        if (loop_ub > 20) {
          loop_ub = 0;
          i6++;
        }
      }
      for (i3 = 0; i3 < b_tmp_size_idx_1; i3++) {
        for (i5 = 0; i5 < b_tmp_size_idx_2; i5++) {
          sdvals[(i5 + sdvals.size(2) * i3) +
                 sdvals.size(2) * sdvals.size(1) * b_q] =
              b_tmp_data[i5 + b_tmp_size_idx_2 * i3];
        }
      }
    }
  } break;
  default: {
    ::coder::SizeType ub_loop;
    m2cAssert(etype == 85, "Only support up to sextic");
    //  Serial mode
    ub_loop = cs.size(0) - 1;
    for (::coder::SizeType b_q = 0; b_q <= ub_loop; b_q++) {
      tri_fek_28(cs[cs.size(1) * b_q], cs[cs.size(1) * b_q + 1], dv7, dv6);
      loop_ub = sfvals.size(1);
      for (i3 = 0; i3 < loop_ub; i3++) {
        sfvals[i3 + sfvals.size(1) * b_q] = dv7[i3];
      }
      b_unnamed_idx_1 = static_cast<int16_T>(sdvals.size(1));
      b_unnamed_idx_2 = static_cast<int16_T>(sdvals.size(2));
      i3 = 0;
      i5 = 0;
      loop_ub = 0;
      i6 = 0;
      b_tmp_size_idx_2 = sdvals.size(2);
      b_tmp_size_idx_1 = sdvals.size(1);
      for (i7 = 0; i7 < b_unnamed_idx_1 * b_unnamed_idx_2; i7++) {
        b_tmp_data[i5 + b_tmp_size_idx_2 * i3] = dv6[i6 + (loop_ub << 1)];
        loop_ub++;
        i3++;
        if (i3 > b_tmp_size_idx_1 - 1) {
          i3 = 0;
          i5++;
        }
        if (loop_ub > 27) {
          loop_ub = 0;
          i6++;
        }
      }
      for (i3 = 0; i3 < b_tmp_size_idx_1; i3++) {
        for (i5 = 0; i5 < b_tmp_size_idx_2; i5++) {
          sdvals[(i5 + sdvals.size(2) * i3) +
                 sdvals.size(2) * sdvals.size(1) * b_q] =
              b_tmp_data[i5 + b_tmp_size_idx_2 * i3];
        }
      }
    }
  } break;
  }
}

// sfe2_tabulate_shapefuncs - Tabulate shape functions and sdvals at given
static inline
void sfe2_tabulate_shapefuncs(::coder::SizeType etype,
                                     const ::coder::array<real_T, 2U> &cs,
                                     ::coder::array<real_T, 2U> &sfvals,
                                     ::coder::array<real_T, 3U> &sdvals)
{
  switch (etype & 3) {
  case 0:
    //  equi kernel
    if ((etype >> 5 & 7) == 2) {
      sfe2_tabulate_equi_tri(etype, cs, sfvals, sdvals);
    } else {
      sfe2_tabulate_equi_quad(etype, cs, sfvals, sdvals);
    }
    break;
  case 1:
    //  GL kernel
    if ((etype >> 5 & 7) == 2) {
      sfe2_tabulate_gl_tri(etype, cs, sfvals, sdvals);
    } else {
      sfe2_tabulate_gl_quad(etype, cs, sfvals, sdvals);
    }
    break;
  default:
    //  FEK kernel
    sfe2_tabulate_fek_tri(etype, cs, sfvals, sdvals);
    break;
  }
}

static void sfe3_tabulate_equi_hexa(::coder::SizeType etype,
                                    const ::coder::array<real_T, 2U> &cs,
                                    ::coder::array<real_T, 2U> &sfvals,
                                    ::coder::array<real_T, 3U> &sdvals)
{
  real_T b_tmp_data[1029];
  real_T dv10[1029];
  real_T tmp_data[1029];
  real_T dv6[648];
  real_T dv5[375];
  real_T dv11[343];
  real_T dv9[216];
  real_T dv4[192];
  real_T dv8[125];
  real_T dv7[64];
  ::coder::SizeType b_tmp_size_idx_1;
  ::coder::SizeType b_tmp_size_idx_2;
  ::coder::SizeType i;
  ::coder::SizeType i4;
  ::coder::SizeType i5;
  ::coder::SizeType i6;
  ::coder::SizeType i7;
  ::coder::SizeType loop_ub;
  ::coder::SizeType nqp;
  int16_T b_unnamed_idx_1;
  int16_T b_unnamed_idx_2;
  //  hex
  nqp = cs.size(0) - 1;
  i = iv[etype - 1];
  sfvals.set_size(cs.size(0), i);
  sdvals.set_size(cs.size(0), i, cs.size(1));
  switch (etype) {
  case 228: {
    for (::coder::SizeType q{0}; q <= nqp; q++) {
      real_T dv2[24];
      real_T dv[8];
      ::coder::SizeType i1;
      ::coder::SizeType i2;
      ::coder::SizeType tmp_size_idx_1;
      ::coder::SizeType tmp_size_idx_2;
      ::coder::SizeType ub_loop;
      int16_T unnamed_idx_1;
      int16_T unnamed_idx_2;
      ::sfe_sfuncs::hexa_8_sfunc(cs[cs.size(1) * q], cs[cs.size(1) * q + 1],
                                 cs[cs.size(1) * q + 2], &dv[0], &dv2[0]);
      ub_loop = sfvals.size(1);
      for (i = 0; i < ub_loop; i++) {
        sfvals[i + sfvals.size(1) * q] = dv[i];
      }
      unnamed_idx_1 = static_cast<int16_T>(sdvals.size(1));
      unnamed_idx_2 = static_cast<int16_T>(sdvals.size(2));
      i = 0;
      i1 = 0;
      ub_loop = 0;
      i2 = 0;
      tmp_size_idx_2 = sdvals.size(2);
      tmp_size_idx_1 = sdvals.size(1);
      for (::coder::SizeType i3{0}; i3 < unnamed_idx_1 * unnamed_idx_2; i3++) {
        tmp_data[i1 + tmp_size_idx_2 * i] = dv2[i2 + 3 * ub_loop];
        ub_loop++;
        i++;
        if (i > tmp_size_idx_1 - 1) {
          i = 0;
          i1++;
        }
        if (ub_loop > 7) {
          ub_loop = 0;
          i2++;
        }
      }
      for (i = 0; i < tmp_size_idx_1; i++) {
        for (i1 = 0; i1 < tmp_size_idx_2; i1++) {
          sdvals[(i1 + sdvals.size(2) * i) +
                 sdvals.size(2) * sdvals.size(1) * q] =
              tmp_data[i1 + tmp_size_idx_2 * i];
        }
      }
    }
  } break;
  case 232: {
    for (::coder::SizeType q{0}; q <= nqp; q++) {
      real_T dv3[81];
      real_T dv1[27];
      ::coder::SizeType i1;
      ::coder::SizeType i2;
      ::coder::SizeType tmp_size_idx_1;
      ::coder::SizeType tmp_size_idx_2;
      ::coder::SizeType ub_loop;
      int16_T unnamed_idx_1;
      int16_T unnamed_idx_2;
      ::sfe_sfuncs::hexa_27_sfunc(cs[cs.size(1) * q], cs[cs.size(1) * q + 1],
                                  cs[cs.size(1) * q + 2], &dv1[0], &dv3[0]);
      ub_loop = sfvals.size(1);
      for (i = 0; i < ub_loop; i++) {
        sfvals[i + sfvals.size(1) * q] = dv1[i];
      }
      unnamed_idx_1 = static_cast<int16_T>(sdvals.size(1));
      unnamed_idx_2 = static_cast<int16_T>(sdvals.size(2));
      i = 0;
      i1 = 0;
      ub_loop = 0;
      i2 = 0;
      tmp_size_idx_2 = sdvals.size(2);
      tmp_size_idx_1 = sdvals.size(1);
      for (::coder::SizeType i3{0}; i3 < unnamed_idx_1 * unnamed_idx_2; i3++) {
        tmp_data[i1 + tmp_size_idx_2 * i] = dv3[i2 + 3 * ub_loop];
        ub_loop++;
        i++;
        if (i > tmp_size_idx_1 - 1) {
          i = 0;
          i1++;
        }
        if (ub_loop > 26) {
          ub_loop = 0;
          i2++;
        }
      }
      for (i = 0; i < tmp_size_idx_1; i++) {
        for (i1 = 0; i1 < tmp_size_idx_2; i1++) {
          sdvals[(i1 + sdvals.size(2) * i) +
                 sdvals.size(2) * sdvals.size(1) * q] =
              tmp_data[i1 + tmp_size_idx_2 * i];
        }
      }
    }
  } break;
  case 236: {
    ::coder::SizeType ub_loop;
    //  Serial mode
    ub_loop = cs.size(0) - 1;
    for (::coder::SizeType b_q = 0; b_q <= ub_loop; b_q++) {
      hexa_64(cs[cs.size(1) * b_q], cs[cs.size(1) * b_q + 1],
              cs[cs.size(1) * b_q + 2], dv7, dv4);
      loop_ub = sfvals.size(1);
      for (i4 = 0; i4 < loop_ub; i4++) {
        sfvals[i4 + sfvals.size(1) * b_q] = dv7[i4];
      }
      b_unnamed_idx_1 = static_cast<int16_T>(sdvals.size(1));
      b_unnamed_idx_2 = static_cast<int16_T>(sdvals.size(2));
      i4 = 0;
      i5 = 0;
      loop_ub = 0;
      i6 = 0;
      b_tmp_size_idx_2 = sdvals.size(2);
      b_tmp_size_idx_1 = sdvals.size(1);
      for (i7 = 0; i7 < b_unnamed_idx_1 * b_unnamed_idx_2; i7++) {
        b_tmp_data[i5 + b_tmp_size_idx_2 * i4] = dv4[i6 + 3 * loop_ub];
        loop_ub++;
        i4++;
        if (i4 > b_tmp_size_idx_1 - 1) {
          i4 = 0;
          i5++;
        }
        if (loop_ub > 63) {
          loop_ub = 0;
          i6++;
        }
      }
      for (i4 = 0; i4 < b_tmp_size_idx_1; i4++) {
        for (i5 = 0; i5 < b_tmp_size_idx_2; i5++) {
          sdvals[(i5 + sdvals.size(2) * i4) +
                 sdvals.size(2) * sdvals.size(1) * b_q] =
              b_tmp_data[i5 + b_tmp_size_idx_2 * i4];
        }
      }
    }
  } break;
  case 240: {
    ::coder::SizeType ub_loop;
    //  Serial mode
    ub_loop = cs.size(0) - 1;
    for (::coder::SizeType b_q = 0; b_q <= ub_loop; b_q++) {
      hexa_125(cs[cs.size(1) * b_q], cs[cs.size(1) * b_q + 1],
               cs[cs.size(1) * b_q + 2], dv8, dv5);
      loop_ub = sfvals.size(1);
      for (i4 = 0; i4 < loop_ub; i4++) {
        sfvals[i4 + sfvals.size(1) * b_q] = dv8[i4];
      }
      b_unnamed_idx_1 = static_cast<int16_T>(sdvals.size(1));
      b_unnamed_idx_2 = static_cast<int16_T>(sdvals.size(2));
      i4 = 0;
      i5 = 0;
      loop_ub = 0;
      i6 = 0;
      b_tmp_size_idx_2 = sdvals.size(2);
      b_tmp_size_idx_1 = sdvals.size(1);
      for (i7 = 0; i7 < b_unnamed_idx_1 * b_unnamed_idx_2; i7++) {
        b_tmp_data[i5 + b_tmp_size_idx_2 * i4] = dv5[i6 + 3 * loop_ub];
        loop_ub++;
        i4++;
        if (i4 > b_tmp_size_idx_1 - 1) {
          i4 = 0;
          i5++;
        }
        if (loop_ub > 124) {
          loop_ub = 0;
          i6++;
        }
      }
      for (i4 = 0; i4 < b_tmp_size_idx_1; i4++) {
        for (i5 = 0; i5 < b_tmp_size_idx_2; i5++) {
          sdvals[(i5 + sdvals.size(2) * i4) +
                 sdvals.size(2) * sdvals.size(1) * b_q] =
              b_tmp_data[i5 + b_tmp_size_idx_2 * i4];
        }
      }
    }
  } break;
  case 244: {
    ::coder::SizeType ub_loop;
    //  Serial mode
    ub_loop = cs.size(0) - 1;
    for (::coder::SizeType b_q = 0; b_q <= ub_loop; b_q++) {
      hexa_216(cs[cs.size(1) * b_q], cs[cs.size(1) * b_q + 1],
               cs[cs.size(1) * b_q + 2], dv9, dv6);
      loop_ub = sfvals.size(1);
      for (i4 = 0; i4 < loop_ub; i4++) {
        sfvals[i4 + sfvals.size(1) * b_q] = dv9[i4];
      }
      b_unnamed_idx_1 = static_cast<int16_T>(sdvals.size(1));
      b_unnamed_idx_2 = static_cast<int16_T>(sdvals.size(2));
      i4 = 0;
      i5 = 0;
      loop_ub = 0;
      i6 = 0;
      b_tmp_size_idx_2 = sdvals.size(2);
      b_tmp_size_idx_1 = sdvals.size(1);
      for (i7 = 0; i7 < b_unnamed_idx_1 * b_unnamed_idx_2; i7++) {
        b_tmp_data[i5 + b_tmp_size_idx_2 * i4] = dv6[i6 + 3 * loop_ub];
        loop_ub++;
        i4++;
        if (i4 > b_tmp_size_idx_1 - 1) {
          i4 = 0;
          i5++;
        }
        if (loop_ub > 215) {
          loop_ub = 0;
          i6++;
        }
      }
      for (i4 = 0; i4 < b_tmp_size_idx_1; i4++) {
        for (i5 = 0; i5 < b_tmp_size_idx_2; i5++) {
          sdvals[(i5 + sdvals.size(2) * i4) +
                 sdvals.size(2) * sdvals.size(1) * b_q] =
              b_tmp_data[i5 + b_tmp_size_idx_2 * i4];
        }
      }
    }
  } break;
  default: {
    ::coder::SizeType ub_loop;
    m2cAssert(etype == 248, "Hex elements supports up to sextic.");
    //  Serial mode
    ub_loop = cs.size(0) - 1;
    for (::coder::SizeType b_q = 0; b_q <= ub_loop; b_q++) {
      hexa_343(cs[cs.size(1) * b_q], cs[cs.size(1) * b_q + 1],
               cs[cs.size(1) * b_q + 2], dv11, dv10);
      loop_ub = sfvals.size(1);
      for (i4 = 0; i4 < loop_ub; i4++) {
        sfvals[i4 + sfvals.size(1) * b_q] = dv11[i4];
      }
      b_unnamed_idx_1 = static_cast<int16_T>(sdvals.size(1));
      b_unnamed_idx_2 = static_cast<int16_T>(sdvals.size(2));
      i4 = 0;
      i5 = 0;
      loop_ub = 0;
      i6 = 0;
      b_tmp_size_idx_2 = sdvals.size(2);
      b_tmp_size_idx_1 = sdvals.size(1);
      for (i7 = 0; i7 < b_unnamed_idx_1 * b_unnamed_idx_2; i7++) {
        b_tmp_data[i5 + b_tmp_size_idx_2 * i4] = dv10[i6 + 3 * loop_ub];
        loop_ub++;
        i4++;
        if (i4 > b_tmp_size_idx_1 - 1) {
          i4 = 0;
          i5++;
        }
        if (loop_ub > 342) {
          loop_ub = 0;
          i6++;
        }
      }
      for (i4 = 0; i4 < b_tmp_size_idx_1; i4++) {
        for (i5 = 0; i5 < b_tmp_size_idx_2; i5++) {
          sdvals[(i5 + sdvals.size(2) * i4) +
                 sdvals.size(2) * sdvals.size(1) * b_q] =
              b_tmp_data[i5 + b_tmp_size_idx_2 * i4];
        }
      }
    }
  } break;
  }
}

static void sfe3_tabulate_equi_prism(::coder::SizeType etype,
                                     const ::coder::array<real_T, 2U> &cs,
                                     ::coder::array<real_T, 2U> &sfvals,
                                     ::coder::array<real_T, 3U> &sdvals)
{
  real_T b_tmp_data[1029];
  real_T tmp_data[1029];
  real_T dv9[588];
  real_T dv5[378];
  real_T dv4[225];
  real_T dv10[196];
  real_T dv8[126];
  real_T dv3[120];
  real_T dv7[75];
  real_T dv6[40];
  ::coder::SizeType b_tmp_size_idx_1;
  ::coder::SizeType b_tmp_size_idx_2;
  ::coder::SizeType i;
  ::coder::SizeType i4;
  ::coder::SizeType i5;
  ::coder::SizeType i6;
  ::coder::SizeType i7;
  ::coder::SizeType loop_ub;
  ::coder::SizeType nqp;
  int16_T b_unnamed_idx_1;
  int16_T b_unnamed_idx_2;
  //  prisms
  nqp = cs.size(0) - 1;
  i = iv[etype - 1];
  sfvals.set_size(cs.size(0), i);
  sdvals.set_size(cs.size(0), i, cs.size(1));
  switch (etype) {
  case 196: {
    for (::coder::SizeType q{0}; q <= nqp; q++) {
      real_T dv1[18];
      real_T dv[6];
      ::coder::SizeType i1;
      ::coder::SizeType i2;
      ::coder::SizeType tmp_size_idx_1;
      ::coder::SizeType tmp_size_idx_2;
      ::coder::SizeType ub_loop;
      int16_T unnamed_idx_1;
      int16_T unnamed_idx_2;
      ::sfe_sfuncs::prism_6_sfunc(cs[cs.size(1) * q], cs[cs.size(1) * q + 1],
                                  cs[cs.size(1) * q + 2], &dv[0], &dv1[0]);
      ub_loop = sfvals.size(1);
      for (i = 0; i < ub_loop; i++) {
        sfvals[i + sfvals.size(1) * q] = dv[i];
      }
      unnamed_idx_1 = static_cast<int16_T>(sdvals.size(1));
      unnamed_idx_2 = static_cast<int16_T>(sdvals.size(2));
      i = 0;
      i1 = 0;
      ub_loop = 0;
      i2 = 0;
      tmp_size_idx_2 = sdvals.size(2);
      tmp_size_idx_1 = sdvals.size(1);
      for (::coder::SizeType i3{0}; i3 < unnamed_idx_1 * unnamed_idx_2; i3++) {
        tmp_data[i1 + tmp_size_idx_2 * i] = dv1[i2 + 3 * ub_loop];
        ub_loop++;
        i++;
        if (i > tmp_size_idx_1 - 1) {
          i = 0;
          i1++;
        }
        if (ub_loop > 5) {
          ub_loop = 0;
          i2++;
        }
      }
      for (i = 0; i < tmp_size_idx_1; i++) {
        for (i1 = 0; i1 < tmp_size_idx_2; i1++) {
          sdvals[(i1 + sdvals.size(2) * i) +
                 sdvals.size(2) * sdvals.size(1) * q] =
              tmp_data[i1 + tmp_size_idx_2 * i];
        }
      }
    }
  } break;
  case 200: {
    for (::coder::SizeType q{0}; q <= nqp; q++) {
      real_T dv2[54];
      real_T dv1[18];
      ::coder::SizeType i1;
      ::coder::SizeType i2;
      ::coder::SizeType tmp_size_idx_1;
      ::coder::SizeType tmp_size_idx_2;
      ::coder::SizeType ub_loop;
      int16_T unnamed_idx_1;
      int16_T unnamed_idx_2;
      ::sfe_sfuncs::prism_18_sfunc(cs[cs.size(1) * q], cs[cs.size(1) * q + 1],
                                   cs[cs.size(1) * q + 2], &dv1[0], &dv2[0]);
      ub_loop = sfvals.size(1);
      for (i = 0; i < ub_loop; i++) {
        sfvals[i + sfvals.size(1) * q] = dv1[i];
      }
      unnamed_idx_1 = static_cast<int16_T>(sdvals.size(1));
      unnamed_idx_2 = static_cast<int16_T>(sdvals.size(2));
      i = 0;
      i1 = 0;
      ub_loop = 0;
      i2 = 0;
      tmp_size_idx_2 = sdvals.size(2);
      tmp_size_idx_1 = sdvals.size(1);
      for (::coder::SizeType i3{0}; i3 < unnamed_idx_1 * unnamed_idx_2; i3++) {
        tmp_data[i1 + tmp_size_idx_2 * i] = dv2[i2 + 3 * ub_loop];
        ub_loop++;
        i++;
        if (i > tmp_size_idx_1 - 1) {
          i = 0;
          i1++;
        }
        if (ub_loop > 17) {
          ub_loop = 0;
          i2++;
        }
      }
      for (i = 0; i < tmp_size_idx_1; i++) {
        for (i1 = 0; i1 < tmp_size_idx_2; i1++) {
          sdvals[(i1 + sdvals.size(2) * i) +
                 sdvals.size(2) * sdvals.size(1) * q] =
              tmp_data[i1 + tmp_size_idx_2 * i];
        }
      }
    }
  } break;
  case 204: {
    ::coder::SizeType ub_loop;
    //  Serial mode
    ub_loop = cs.size(0) - 1;
    for (::coder::SizeType b_q = 0; b_q <= ub_loop; b_q++) {
      prism_40(cs[cs.size(1) * b_q], cs[cs.size(1) * b_q + 1],
               cs[cs.size(1) * b_q + 2], dv6, dv3);
      loop_ub = sfvals.size(1);
      for (i4 = 0; i4 < loop_ub; i4++) {
        sfvals[i4 + sfvals.size(1) * b_q] = dv6[i4];
      }
      b_unnamed_idx_1 = static_cast<int16_T>(sdvals.size(1));
      b_unnamed_idx_2 = static_cast<int16_T>(sdvals.size(2));
      i4 = 0;
      i5 = 0;
      loop_ub = 0;
      i6 = 0;
      b_tmp_size_idx_2 = sdvals.size(2);
      b_tmp_size_idx_1 = sdvals.size(1);
      for (i7 = 0; i7 < b_unnamed_idx_1 * b_unnamed_idx_2; i7++) {
        b_tmp_data[i5 + b_tmp_size_idx_2 * i4] = dv3[i6 + 3 * loop_ub];
        loop_ub++;
        i4++;
        if (i4 > b_tmp_size_idx_1 - 1) {
          i4 = 0;
          i5++;
        }
        if (loop_ub > 39) {
          loop_ub = 0;
          i6++;
        }
      }
      for (i4 = 0; i4 < b_tmp_size_idx_1; i4++) {
        for (i5 = 0; i5 < b_tmp_size_idx_2; i5++) {
          sdvals[(i5 + sdvals.size(2) * i4) +
                 sdvals.size(2) * sdvals.size(1) * b_q] =
              b_tmp_data[i5 + b_tmp_size_idx_2 * i4];
        }
      }
    }
  } break;
  case 208: {
    ::coder::SizeType ub_loop;
    //  Serial mode
    ub_loop = cs.size(0) - 1;
    for (::coder::SizeType b_q = 0; b_q <= ub_loop; b_q++) {
      prism_75(cs[cs.size(1) * b_q], cs[cs.size(1) * b_q + 1],
               cs[cs.size(1) * b_q + 2], dv7, dv4);
      loop_ub = sfvals.size(1);
      for (i4 = 0; i4 < loop_ub; i4++) {
        sfvals[i4 + sfvals.size(1) * b_q] = dv7[i4];
      }
      b_unnamed_idx_1 = static_cast<int16_T>(sdvals.size(1));
      b_unnamed_idx_2 = static_cast<int16_T>(sdvals.size(2));
      i4 = 0;
      i5 = 0;
      loop_ub = 0;
      i6 = 0;
      b_tmp_size_idx_2 = sdvals.size(2);
      b_tmp_size_idx_1 = sdvals.size(1);
      for (i7 = 0; i7 < b_unnamed_idx_1 * b_unnamed_idx_2; i7++) {
        b_tmp_data[i5 + b_tmp_size_idx_2 * i4] = dv4[i6 + 3 * loop_ub];
        loop_ub++;
        i4++;
        if (i4 > b_tmp_size_idx_1 - 1) {
          i4 = 0;
          i5++;
        }
        if (loop_ub > 74) {
          loop_ub = 0;
          i6++;
        }
      }
      for (i4 = 0; i4 < b_tmp_size_idx_1; i4++) {
        for (i5 = 0; i5 < b_tmp_size_idx_2; i5++) {
          sdvals[(i5 + sdvals.size(2) * i4) +
                 sdvals.size(2) * sdvals.size(1) * b_q] =
              b_tmp_data[i5 + b_tmp_size_idx_2 * i4];
        }
      }
    }
  } break;
  case 212: {
    ::coder::SizeType ub_loop;
    //  Serial mode
    ub_loop = cs.size(0) - 1;
    for (::coder::SizeType b_q = 0; b_q <= ub_loop; b_q++) {
      prism_126(cs[cs.size(1) * b_q], cs[cs.size(1) * b_q + 1],
                cs[cs.size(1) * b_q + 2], dv8, dv5);
      loop_ub = sfvals.size(1);
      for (i4 = 0; i4 < loop_ub; i4++) {
        sfvals[i4 + sfvals.size(1) * b_q] = dv8[i4];
      }
      b_unnamed_idx_1 = static_cast<int16_T>(sdvals.size(1));
      b_unnamed_idx_2 = static_cast<int16_T>(sdvals.size(2));
      i4 = 0;
      i5 = 0;
      loop_ub = 0;
      i6 = 0;
      b_tmp_size_idx_2 = sdvals.size(2);
      b_tmp_size_idx_1 = sdvals.size(1);
      for (i7 = 0; i7 < b_unnamed_idx_1 * b_unnamed_idx_2; i7++) {
        b_tmp_data[i5 + b_tmp_size_idx_2 * i4] = dv5[i6 + 3 * loop_ub];
        loop_ub++;
        i4++;
        if (i4 > b_tmp_size_idx_1 - 1) {
          i4 = 0;
          i5++;
        }
        if (loop_ub > 125) {
          loop_ub = 0;
          i6++;
        }
      }
      for (i4 = 0; i4 < b_tmp_size_idx_1; i4++) {
        for (i5 = 0; i5 < b_tmp_size_idx_2; i5++) {
          sdvals[(i5 + sdvals.size(2) * i4) +
                 sdvals.size(2) * sdvals.size(1) * b_q] =
              b_tmp_data[i5 + b_tmp_size_idx_2 * i4];
        }
      }
    }
  } break;
  default: {
    ::coder::SizeType ub_loop;
    m2cAssert(etype == 216, "prismatic elements supports up to sextic.");
    //  Serial mode
    ub_loop = cs.size(0) - 1;
    for (::coder::SizeType b_q = 0; b_q <= ub_loop; b_q++) {
      prism_196(cs[cs.size(1) * b_q], cs[cs.size(1) * b_q + 1],
                cs[cs.size(1) * b_q + 2], dv10, dv9);
      loop_ub = sfvals.size(1);
      for (i4 = 0; i4 < loop_ub; i4++) {
        sfvals[i4 + sfvals.size(1) * b_q] = dv10[i4];
      }
      b_unnamed_idx_1 = static_cast<int16_T>(sdvals.size(1));
      b_unnamed_idx_2 = static_cast<int16_T>(sdvals.size(2));
      i4 = 0;
      i5 = 0;
      loop_ub = 0;
      i6 = 0;
      b_tmp_size_idx_2 = sdvals.size(2);
      b_tmp_size_idx_1 = sdvals.size(1);
      for (i7 = 0; i7 < b_unnamed_idx_1 * b_unnamed_idx_2; i7++) {
        b_tmp_data[i5 + b_tmp_size_idx_2 * i4] = dv9[i6 + 3 * loop_ub];
        loop_ub++;
        i4++;
        if (i4 > b_tmp_size_idx_1 - 1) {
          i4 = 0;
          i5++;
        }
        if (loop_ub > 195) {
          loop_ub = 0;
          i6++;
        }
      }
      for (i4 = 0; i4 < b_tmp_size_idx_1; i4++) {
        for (i5 = 0; i5 < b_tmp_size_idx_2; i5++) {
          sdvals[(i5 + sdvals.size(2) * i4) +
                 sdvals.size(2) * sdvals.size(1) * b_q] =
              b_tmp_data[i5 + b_tmp_size_idx_2 * i4];
        }
      }
    }
  } break;
  }
}

static void sfe3_tabulate_equi_pyra(::coder::SizeType etype,
                                    const ::coder::array<real_T, 2U> &cs,
                                    ::coder::array<real_T, 2U> &sfvals,
                                    ::coder::array<real_T, 3U> &sdvals)
{
  real_T b_tmp_data[1029];
  real_T tmp_data[1029];
  real_T dv6[165];
  real_T dv4[90];
  real_T dv7[55];
  real_T dv5[30];
  ::coder::SizeType b_tmp_size_idx_1;
  ::coder::SizeType b_tmp_size_idx_2;
  ::coder::SizeType i;
  ::coder::SizeType i2;
  ::coder::SizeType i5;
  ::coder::SizeType i6;
  ::coder::SizeType i7;
  ::coder::SizeType loop_ub;
  ::coder::SizeType nqp;
  int16_T b_unnamed_idx_1;
  int16_T b_unnamed_idx_2;
  //  pyra
  nqp = cs.size(0) - 1;
  i = iv[etype - 1];
  sfvals.set_size(cs.size(0), i);
  sdvals.set_size(cs.size(0), i, cs.size(1));
  switch (etype) {
  case 164: {
    for (::coder::SizeType q{0}; q <= nqp; q++) {
      real_T dv2[15];
      real_T dv[5];
      ::coder::SizeType i1;
      ::coder::SizeType i3;
      ::coder::SizeType tmp_size_idx_1;
      ::coder::SizeType tmp_size_idx_2;
      ::coder::SizeType ub_loop;
      int16_T unnamed_idx_1;
      int16_T unnamed_idx_2;
      ::sfe_sfuncs::pyra_5_sfunc(cs[cs.size(1) * q], cs[cs.size(1) * q + 1],
                                 cs[cs.size(1) * q + 2], &dv[0], &dv2[0]);
      ub_loop = sfvals.size(1);
      for (i = 0; i < ub_loop; i++) {
        sfvals[i + sfvals.size(1) * q] = dv[i];
      }
      unnamed_idx_1 = static_cast<int16_T>(sdvals.size(1));
      unnamed_idx_2 = static_cast<int16_T>(sdvals.size(2));
      i = 0;
      i1 = 0;
      ub_loop = 0;
      i3 = 0;
      tmp_size_idx_2 = sdvals.size(2);
      tmp_size_idx_1 = sdvals.size(1);
      for (::coder::SizeType i4{0}; i4 < unnamed_idx_1 * unnamed_idx_2; i4++) {
        tmp_data[i1 + tmp_size_idx_2 * i] = dv2[i3 + 3 * ub_loop];
        ub_loop++;
        i++;
        if (i > tmp_size_idx_1 - 1) {
          i = 0;
          i1++;
        }
        if (ub_loop > 4) {
          ub_loop = 0;
          i3++;
        }
      }
      for (i = 0; i < tmp_size_idx_1; i++) {
        for (i1 = 0; i1 < tmp_size_idx_2; i1++) {
          sdvals[(i1 + sdvals.size(2) * i) +
                 sdvals.size(2) * sdvals.size(1) * q] =
              tmp_data[i1 + tmp_size_idx_2 * i];
        }
      }
    }
  } break;
  case 168: {
    for (::coder::SizeType q{0}; q <= nqp; q++) {
      real_T dv3[42];
      real_T dv1[14];
      ::coder::SizeType i1;
      ::coder::SizeType i3;
      ::coder::SizeType tmp_size_idx_1;
      ::coder::SizeType tmp_size_idx_2;
      ::coder::SizeType ub_loop;
      int16_T unnamed_idx_1;
      int16_T unnamed_idx_2;
      ::sfe_sfuncs::pyra_14_sfunc(cs[cs.size(1) * q], cs[cs.size(1) * q + 1],
                                  cs[cs.size(1) * q + 2], &dv1[0], &dv3[0]);
      ub_loop = sfvals.size(1);
      for (i = 0; i < ub_loop; i++) {
        sfvals[i + sfvals.size(1) * q] = dv1[i];
      }
      unnamed_idx_1 = static_cast<int16_T>(sdvals.size(1));
      unnamed_idx_2 = static_cast<int16_T>(sdvals.size(2));
      i = 0;
      i1 = 0;
      ub_loop = 0;
      i3 = 0;
      tmp_size_idx_2 = sdvals.size(2);
      tmp_size_idx_1 = sdvals.size(1);
      for (::coder::SizeType i4{0}; i4 < unnamed_idx_1 * unnamed_idx_2; i4++) {
        tmp_data[i1 + tmp_size_idx_2 * i] = dv3[i3 + 3 * ub_loop];
        ub_loop++;
        i++;
        if (i > tmp_size_idx_1 - 1) {
          i = 0;
          i1++;
        }
        if (ub_loop > 13) {
          ub_loop = 0;
          i3++;
        }
      }
      for (i = 0; i < tmp_size_idx_1; i++) {
        for (i1 = 0; i1 < tmp_size_idx_2; i1++) {
          sdvals[(i1 + sdvals.size(2) * i) +
                 sdvals.size(2) * sdvals.size(1) * q] =
              tmp_data[i1 + tmp_size_idx_2 * i];
        }
      }
    }
  } break;
  case 172: {
    ::coder::SizeType ub_loop;
    //  Serial mode
    ub_loop = cs.size(0) - 1;
    for (::coder::SizeType b_q = 0; b_q <= ub_loop; b_q++) {
      pyra_30(cs[cs.size(1) * b_q], cs[cs.size(1) * b_q + 1],
              cs[cs.size(1) * b_q + 2], dv5, dv4);
      loop_ub = sfvals.size(1);
      for (i2 = 0; i2 < loop_ub; i2++) {
        sfvals[i2 + sfvals.size(1) * b_q] = dv5[i2];
      }
      b_unnamed_idx_1 = static_cast<int16_T>(sdvals.size(1));
      b_unnamed_idx_2 = static_cast<int16_T>(sdvals.size(2));
      i2 = 0;
      i5 = 0;
      loop_ub = 0;
      i6 = 0;
      b_tmp_size_idx_2 = sdvals.size(2);
      b_tmp_size_idx_1 = sdvals.size(1);
      for (i7 = 0; i7 < b_unnamed_idx_1 * b_unnamed_idx_2; i7++) {
        b_tmp_data[i5 + b_tmp_size_idx_2 * i2] = dv4[i6 + 3 * loop_ub];
        loop_ub++;
        i2++;
        if (i2 > b_tmp_size_idx_1 - 1) {
          i2 = 0;
          i5++;
        }
        if (loop_ub > 29) {
          loop_ub = 0;
          i6++;
        }
      }
      for (i2 = 0; i2 < b_tmp_size_idx_1; i2++) {
        for (i5 = 0; i5 < b_tmp_size_idx_2; i5++) {
          sdvals[(i5 + sdvals.size(2) * i2) +
                 sdvals.size(2) * sdvals.size(1) * b_q] =
              b_tmp_data[i5 + b_tmp_size_idx_2 * i2];
        }
      }
    }
  } break;
  default: {
    ::coder::SizeType ub_loop;
    m2cAssert(etype == 176, "Pyramid only support up to quartic");
    //  Serial mode
    ub_loop = cs.size(0) - 1;
    for (::coder::SizeType b_q = 0; b_q <= ub_loop; b_q++) {
      pyra_55(cs[cs.size(1) * b_q], cs[cs.size(1) * b_q + 1],
              cs[cs.size(1) * b_q + 2], dv7, dv6);
      loop_ub = sfvals.size(1);
      for (i2 = 0; i2 < loop_ub; i2++) {
        sfvals[i2 + sfvals.size(1) * b_q] = dv7[i2];
      }
      b_unnamed_idx_1 = static_cast<int16_T>(sdvals.size(1));
      b_unnamed_idx_2 = static_cast<int16_T>(sdvals.size(2));
      i2 = 0;
      i5 = 0;
      loop_ub = 0;
      i6 = 0;
      b_tmp_size_idx_2 = sdvals.size(2);
      b_tmp_size_idx_1 = sdvals.size(1);
      for (i7 = 0; i7 < b_unnamed_idx_1 * b_unnamed_idx_2; i7++) {
        b_tmp_data[i5 + b_tmp_size_idx_2 * i2] = dv6[i6 + 3 * loop_ub];
        loop_ub++;
        i2++;
        if (i2 > b_tmp_size_idx_1 - 1) {
          i2 = 0;
          i5++;
        }
        if (loop_ub > 54) {
          loop_ub = 0;
          i6++;
        }
      }
      for (i2 = 0; i2 < b_tmp_size_idx_1; i2++) {
        for (i5 = 0; i5 < b_tmp_size_idx_2; i5++) {
          sdvals[(i5 + sdvals.size(2) * i2) +
                 sdvals.size(2) * sdvals.size(1) * b_q] =
              b_tmp_data[i5 + b_tmp_size_idx_2 * i2];
        }
      }
    }
  } break;
  }
}

static void sfe3_tabulate_equi_tet(::coder::SizeType etype,
                                   const ::coder::array<real_T, 2U> &cs,
                                   ::coder::array<real_T, 2U> &sfvals,
                                   ::coder::array<real_T, 3U> &sdvals)
{
  real_T b_tmp_data[1029];
  real_T tmp_data[1029];
  real_T dv10[252];
  real_T dv6[168];
  real_T dv5[105];
  real_T dv11[84];
  real_T dv4[60];
  real_T dv9[56];
  real_T dv8[35];
  real_T dv7[20];
  ::coder::SizeType b_tmp_size_idx_1;
  ::coder::SizeType b_tmp_size_idx_2;
  ::coder::SizeType i;
  ::coder::SizeType i3;
  ::coder::SizeType i5;
  ::coder::SizeType i6;
  ::coder::SizeType i7;
  ::coder::SizeType loop_ub;
  ::coder::SizeType nqp;
  int16_T b_unnamed_idx_1;
  int16_T b_unnamed_idx_2;
  //  tet
  nqp = cs.size(0) - 1;
  i = iv[etype - 1];
  sfvals.set_size(cs.size(0), i);
  sdvals.set_size(cs.size(0), i, cs.size(1));
  switch (etype) {
  case 132: {
    for (::coder::SizeType q{0}; q <= nqp; q++) {
      real_T dv2[12];
      real_T dv[4];
      ::coder::SizeType i1;
      ::coder::SizeType i2;
      ::coder::SizeType tmp_size_idx_1;
      ::coder::SizeType tmp_size_idx_2;
      ::coder::SizeType ub_loop;
      int16_T unnamed_idx_1;
      int16_T unnamed_idx_2;
      ::sfe_sfuncs::tet_4_sfunc(cs[cs.size(1) * q], cs[cs.size(1) * q + 1],
                                cs[cs.size(1) * q + 2], &dv[0], &dv2[0]);
      ub_loop = sfvals.size(1);
      for (i = 0; i < ub_loop; i++) {
        sfvals[i + sfvals.size(1) * q] = dv[i];
      }
      unnamed_idx_1 = static_cast<int16_T>(sdvals.size(1));
      unnamed_idx_2 = static_cast<int16_T>(sdvals.size(2));
      i = 0;
      i1 = 0;
      ub_loop = 0;
      i2 = 0;
      tmp_size_idx_2 = sdvals.size(2);
      tmp_size_idx_1 = sdvals.size(1);
      for (::coder::SizeType i4{0}; i4 < unnamed_idx_1 * unnamed_idx_2; i4++) {
        tmp_data[i1 + tmp_size_idx_2 * i] = dv2[i2 + 3 * ub_loop];
        ub_loop++;
        i++;
        if (i > tmp_size_idx_1 - 1) {
          i = 0;
          i1++;
        }
        if (ub_loop > 3) {
          ub_loop = 0;
          i2++;
        }
      }
      for (i = 0; i < tmp_size_idx_1; i++) {
        for (i1 = 0; i1 < tmp_size_idx_2; i1++) {
          sdvals[(i1 + sdvals.size(2) * i) +
                 sdvals.size(2) * sdvals.size(1) * q] =
              tmp_data[i1 + tmp_size_idx_2 * i];
        }
      }
    }
  } break;
  case 136: {
    for (::coder::SizeType q{0}; q <= nqp; q++) {
      real_T dv3[30];
      real_T dv1[10];
      ::coder::SizeType i1;
      ::coder::SizeType i2;
      ::coder::SizeType tmp_size_idx_1;
      ::coder::SizeType tmp_size_idx_2;
      ::coder::SizeType ub_loop;
      int16_T unnamed_idx_1;
      int16_T unnamed_idx_2;
      ::sfe_sfuncs::tet_10_sfunc(cs[cs.size(1) * q], cs[cs.size(1) * q + 1],
                                 cs[cs.size(1) * q + 2], &dv1[0], &dv3[0]);
      ub_loop = sfvals.size(1);
      for (i = 0; i < ub_loop; i++) {
        sfvals[i + sfvals.size(1) * q] = dv1[i];
      }
      unnamed_idx_1 = static_cast<int16_T>(sdvals.size(1));
      unnamed_idx_2 = static_cast<int16_T>(sdvals.size(2));
      i = 0;
      i1 = 0;
      ub_loop = 0;
      i2 = 0;
      tmp_size_idx_2 = sdvals.size(2);
      tmp_size_idx_1 = sdvals.size(1);
      for (::coder::SizeType i4{0}; i4 < unnamed_idx_1 * unnamed_idx_2; i4++) {
        tmp_data[i1 + tmp_size_idx_2 * i] = dv3[i2 + 3 * ub_loop];
        ub_loop++;
        i++;
        if (i > tmp_size_idx_1 - 1) {
          i = 0;
          i1++;
        }
        if (ub_loop > 9) {
          ub_loop = 0;
          i2++;
        }
      }
      for (i = 0; i < tmp_size_idx_1; i++) {
        for (i1 = 0; i1 < tmp_size_idx_2; i1++) {
          sdvals[(i1 + sdvals.size(2) * i) +
                 sdvals.size(2) * sdvals.size(1) * q] =
              tmp_data[i1 + tmp_size_idx_2 * i];
        }
      }
    }
  } break;
  case 140: {
    ::coder::SizeType ub_loop;
    //  Serial mode
    ub_loop = cs.size(0) - 1;
    for (::coder::SizeType b_q = 0; b_q <= ub_loop; b_q++) {
      tet_20(cs[cs.size(1) * b_q], cs[cs.size(1) * b_q + 1],
             cs[cs.size(1) * b_q + 2], dv7, dv4);
      loop_ub = sfvals.size(1);
      for (i3 = 0; i3 < loop_ub; i3++) {
        sfvals[i3 + sfvals.size(1) * b_q] = dv7[i3];
      }
      b_unnamed_idx_1 = static_cast<int16_T>(sdvals.size(1));
      b_unnamed_idx_2 = static_cast<int16_T>(sdvals.size(2));
      i3 = 0;
      i5 = 0;
      loop_ub = 0;
      i6 = 0;
      b_tmp_size_idx_2 = sdvals.size(2);
      b_tmp_size_idx_1 = sdvals.size(1);
      for (i7 = 0; i7 < b_unnamed_idx_1 * b_unnamed_idx_2; i7++) {
        b_tmp_data[i5 + b_tmp_size_idx_2 * i3] = dv4[i6 + 3 * loop_ub];
        loop_ub++;
        i3++;
        if (i3 > b_tmp_size_idx_1 - 1) {
          i3 = 0;
          i5++;
        }
        if (loop_ub > 19) {
          loop_ub = 0;
          i6++;
        }
      }
      for (i3 = 0; i3 < b_tmp_size_idx_1; i3++) {
        for (i5 = 0; i5 < b_tmp_size_idx_2; i5++) {
          sdvals[(i5 + sdvals.size(2) * i3) +
                 sdvals.size(2) * sdvals.size(1) * b_q] =
              b_tmp_data[i5 + b_tmp_size_idx_2 * i3];
        }
      }
    }
  } break;
  case 144: {
    ::coder::SizeType ub_loop;
    //  Serial mode
    ub_loop = cs.size(0) - 1;
    for (::coder::SizeType b_q = 0; b_q <= ub_loop; b_q++) {
      tet_35(cs[cs.size(1) * b_q], cs[cs.size(1) * b_q + 1],
             cs[cs.size(1) * b_q + 2], dv8, dv5);
      loop_ub = sfvals.size(1);
      for (i3 = 0; i3 < loop_ub; i3++) {
        sfvals[i3 + sfvals.size(1) * b_q] = dv8[i3];
      }
      b_unnamed_idx_1 = static_cast<int16_T>(sdvals.size(1));
      b_unnamed_idx_2 = static_cast<int16_T>(sdvals.size(2));
      i3 = 0;
      i5 = 0;
      loop_ub = 0;
      i6 = 0;
      b_tmp_size_idx_2 = sdvals.size(2);
      b_tmp_size_idx_1 = sdvals.size(1);
      for (i7 = 0; i7 < b_unnamed_idx_1 * b_unnamed_idx_2; i7++) {
        b_tmp_data[i5 + b_tmp_size_idx_2 * i3] = dv5[i6 + 3 * loop_ub];
        loop_ub++;
        i3++;
        if (i3 > b_tmp_size_idx_1 - 1) {
          i3 = 0;
          i5++;
        }
        if (loop_ub > 34) {
          loop_ub = 0;
          i6++;
        }
      }
      for (i3 = 0; i3 < b_tmp_size_idx_1; i3++) {
        for (i5 = 0; i5 < b_tmp_size_idx_2; i5++) {
          sdvals[(i5 + sdvals.size(2) * i3) +
                 sdvals.size(2) * sdvals.size(1) * b_q] =
              b_tmp_data[i5 + b_tmp_size_idx_2 * i3];
        }
      }
    }
  } break;
  case 148: {
    ::coder::SizeType ub_loop;
    //  Serial mode
    ub_loop = cs.size(0) - 1;
    for (::coder::SizeType b_q = 0; b_q <= ub_loop; b_q++) {
      tet_56(cs[cs.size(1) * b_q], cs[cs.size(1) * b_q + 1],
             cs[cs.size(1) * b_q + 2], dv9, dv6);
      loop_ub = sfvals.size(1);
      for (i3 = 0; i3 < loop_ub; i3++) {
        sfvals[i3 + sfvals.size(1) * b_q] = dv9[i3];
      }
      b_unnamed_idx_1 = static_cast<int16_T>(sdvals.size(1));
      b_unnamed_idx_2 = static_cast<int16_T>(sdvals.size(2));
      i3 = 0;
      i5 = 0;
      loop_ub = 0;
      i6 = 0;
      b_tmp_size_idx_2 = sdvals.size(2);
      b_tmp_size_idx_1 = sdvals.size(1);
      for (i7 = 0; i7 < b_unnamed_idx_1 * b_unnamed_idx_2; i7++) {
        b_tmp_data[i5 + b_tmp_size_idx_2 * i3] = dv6[i6 + 3 * loop_ub];
        loop_ub++;
        i3++;
        if (i3 > b_tmp_size_idx_1 - 1) {
          i3 = 0;
          i5++;
        }
        if (loop_ub > 55) {
          loop_ub = 0;
          i6++;
        }
      }
      for (i3 = 0; i3 < b_tmp_size_idx_1; i3++) {
        for (i5 = 0; i5 < b_tmp_size_idx_2; i5++) {
          sdvals[(i5 + sdvals.size(2) * i3) +
                 sdvals.size(2) * sdvals.size(1) * b_q] =
              b_tmp_data[i5 + b_tmp_size_idx_2 * i3];
        }
      }
    }
  } break;
  default: {
    ::coder::SizeType ub_loop;
    m2cAssert(etype == 152, "equidistant tets only supported up to sextic");
    //  Serial mode
    ub_loop = cs.size(0) - 1;
    for (::coder::SizeType b_q = 0; b_q <= ub_loop; b_q++) {
      tet_84(cs[cs.size(1) * b_q], cs[cs.size(1) * b_q + 1],
             cs[cs.size(1) * b_q + 2], dv11, dv10);
      loop_ub = sfvals.size(1);
      for (i3 = 0; i3 < loop_ub; i3++) {
        sfvals[i3 + sfvals.size(1) * b_q] = dv11[i3];
      }
      b_unnamed_idx_1 = static_cast<int16_T>(sdvals.size(1));
      b_unnamed_idx_2 = static_cast<int16_T>(sdvals.size(2));
      i3 = 0;
      i5 = 0;
      loop_ub = 0;
      i6 = 0;
      b_tmp_size_idx_2 = sdvals.size(2);
      b_tmp_size_idx_1 = sdvals.size(1);
      for (i7 = 0; i7 < b_unnamed_idx_1 * b_unnamed_idx_2; i7++) {
        b_tmp_data[i5 + b_tmp_size_idx_2 * i3] = dv10[i6 + 3 * loop_ub];
        loop_ub++;
        i3++;
        if (i3 > b_tmp_size_idx_1 - 1) {
          i3 = 0;
          i5++;
        }
        if (loop_ub > 83) {
          loop_ub = 0;
          i6++;
        }
      }
      for (i3 = 0; i3 < b_tmp_size_idx_1; i3++) {
        for (i5 = 0; i5 < b_tmp_size_idx_2; i5++) {
          sdvals[(i5 + sdvals.size(2) * i3) +
                 sdvals.size(2) * sdvals.size(1) * b_q] =
              b_tmp_data[i5 + b_tmp_size_idx_2 * i3];
        }
      }
    }
  } break;
  }
}

static void sfe3_tabulate_gl_hexa(::coder::SizeType etype,
                                  const ::coder::array<real_T, 2U> &cs,
                                  ::coder::array<real_T, 2U> &sfvals,
                                  ::coder::array<real_T, 3U> &sdvals)
{
  real_T dv6[1029];
  real_T tmp_data[1029];
  real_T dv2[648];
  real_T dv1[375];
  real_T dv7[343];
  real_T dv5[216];
  real_T dv[192];
  real_T dv4[125];
  real_T dv3[64];
  ::coder::SizeType i;
  ::coder::SizeType i1;
  ::coder::SizeType i2;
  ::coder::SizeType i3;
  ::coder::SizeType loop_ub;
  ::coder::SizeType tmp_size_idx_1;
  ::coder::SizeType tmp_size_idx_2;
  ::coder::SizeType ub_loop;
  int16_T unnamed_idx_1;
  int16_T unnamed_idx_2;
  //  hex
  ub_loop = iv[etype - 1];
  sfvals.set_size(cs.size(0), ub_loop);
  sdvals.set_size(cs.size(0), ub_loop, cs.size(1));
  switch (etype) {
  case 237:
    //  Serial mode
    ub_loop = cs.size(0) - 1;
    for (::coder::SizeType q = 0; q <= ub_loop; q++) {
      hexa_gl_64(cs[cs.size(1) * q], cs[cs.size(1) * q + 1],
                 cs[cs.size(1) * q + 2], dv3, dv);
      loop_ub = sfvals.size(1);
      for (i = 0; i < loop_ub; i++) {
        sfvals[i + sfvals.size(1) * q] = dv3[i];
      }
      unnamed_idx_1 = static_cast<int16_T>(sdvals.size(1));
      unnamed_idx_2 = static_cast<int16_T>(sdvals.size(2));
      i = 0;
      i1 = 0;
      loop_ub = 0;
      i2 = 0;
      tmp_size_idx_2 = sdvals.size(2);
      tmp_size_idx_1 = sdvals.size(1);
      for (i3 = 0; i3 < unnamed_idx_1 * unnamed_idx_2; i3++) {
        tmp_data[i1 + tmp_size_idx_2 * i] = dv[i2 + 3 * loop_ub];
        loop_ub++;
        i++;
        if (i > tmp_size_idx_1 - 1) {
          i = 0;
          i1++;
        }
        if (loop_ub > 63) {
          loop_ub = 0;
          i2++;
        }
      }
      for (i = 0; i < tmp_size_idx_1; i++) {
        for (i1 = 0; i1 < tmp_size_idx_2; i1++) {
          sdvals[(i1 + sdvals.size(2) * i) +
                 sdvals.size(2) * sdvals.size(1) * q] =
              tmp_data[i1 + tmp_size_idx_2 * i];
        }
      }
    }
    break;
  case 241:
    //  Serial mode
    ub_loop = cs.size(0) - 1;
    for (::coder::SizeType q = 0; q <= ub_loop; q++) {
      hexa_gl_125(cs[cs.size(1) * q], cs[cs.size(1) * q + 1],
                  cs[cs.size(1) * q + 2], dv4, dv1);
      loop_ub = sfvals.size(1);
      for (i = 0; i < loop_ub; i++) {
        sfvals[i + sfvals.size(1) * q] = dv4[i];
      }
      unnamed_idx_1 = static_cast<int16_T>(sdvals.size(1));
      unnamed_idx_2 = static_cast<int16_T>(sdvals.size(2));
      i = 0;
      i1 = 0;
      loop_ub = 0;
      i2 = 0;
      tmp_size_idx_2 = sdvals.size(2);
      tmp_size_idx_1 = sdvals.size(1);
      for (i3 = 0; i3 < unnamed_idx_1 * unnamed_idx_2; i3++) {
        tmp_data[i1 + tmp_size_idx_2 * i] = dv1[i2 + 3 * loop_ub];
        loop_ub++;
        i++;
        if (i > tmp_size_idx_1 - 1) {
          i = 0;
          i1++;
        }
        if (loop_ub > 124) {
          loop_ub = 0;
          i2++;
        }
      }
      for (i = 0; i < tmp_size_idx_1; i++) {
        for (i1 = 0; i1 < tmp_size_idx_2; i1++) {
          sdvals[(i1 + sdvals.size(2) * i) +
                 sdvals.size(2) * sdvals.size(1) * q] =
              tmp_data[i1 + tmp_size_idx_2 * i];
        }
      }
    }
    break;
  case 245:
    //  Serial mode
    ub_loop = cs.size(0) - 1;
    for (::coder::SizeType q = 0; q <= ub_loop; q++) {
      hexa_gl_216(cs[cs.size(1) * q], cs[cs.size(1) * q + 1],
                  cs[cs.size(1) * q + 2], dv5, dv2);
      loop_ub = sfvals.size(1);
      for (i = 0; i < loop_ub; i++) {
        sfvals[i + sfvals.size(1) * q] = dv5[i];
      }
      unnamed_idx_1 = static_cast<int16_T>(sdvals.size(1));
      unnamed_idx_2 = static_cast<int16_T>(sdvals.size(2));
      i = 0;
      i1 = 0;
      loop_ub = 0;
      i2 = 0;
      tmp_size_idx_2 = sdvals.size(2);
      tmp_size_idx_1 = sdvals.size(1);
      for (i3 = 0; i3 < unnamed_idx_1 * unnamed_idx_2; i3++) {
        tmp_data[i1 + tmp_size_idx_2 * i] = dv2[i2 + 3 * loop_ub];
        loop_ub++;
        i++;
        if (i > tmp_size_idx_1 - 1) {
          i = 0;
          i1++;
        }
        if (loop_ub > 215) {
          loop_ub = 0;
          i2++;
        }
      }
      for (i = 0; i < tmp_size_idx_1; i++) {
        for (i1 = 0; i1 < tmp_size_idx_2; i1++) {
          sdvals[(i1 + sdvals.size(2) * i) +
                 sdvals.size(2) * sdvals.size(1) * q] =
              tmp_data[i1 + tmp_size_idx_2 * i];
        }
      }
    }
    break;
  default:
    m2cAssert(etype == 249, "Gauss-Lobatto only supports up to sextic.");
    //  Serial mode
    ub_loop = cs.size(0) - 1;
    for (::coder::SizeType q = 0; q <= ub_loop; q++) {
      hexa_gl_343(cs[cs.size(1) * q], cs[cs.size(1) * q + 1],
                  cs[cs.size(1) * q + 2], dv7, dv6);
      loop_ub = sfvals.size(1);
      for (i = 0; i < loop_ub; i++) {
        sfvals[i + sfvals.size(1) * q] = dv7[i];
      }
      unnamed_idx_1 = static_cast<int16_T>(sdvals.size(1));
      unnamed_idx_2 = static_cast<int16_T>(sdvals.size(2));
      i = 0;
      i1 = 0;
      loop_ub = 0;
      i2 = 0;
      tmp_size_idx_2 = sdvals.size(2);
      tmp_size_idx_1 = sdvals.size(1);
      for (i3 = 0; i3 < unnamed_idx_1 * unnamed_idx_2; i3++) {
        tmp_data[i1 + tmp_size_idx_2 * i] = dv6[i2 + 3 * loop_ub];
        loop_ub++;
        i++;
        if (i > tmp_size_idx_1 - 1) {
          i = 0;
          i1++;
        }
        if (loop_ub > 342) {
          loop_ub = 0;
          i2++;
        }
      }
      for (i = 0; i < tmp_size_idx_1; i++) {
        for (i1 = 0; i1 < tmp_size_idx_2; i1++) {
          sdvals[(i1 + sdvals.size(2) * i) +
                 sdvals.size(2) * sdvals.size(1) * q] =
              tmp_data[i1 + tmp_size_idx_2 * i];
        }
      }
    }
    break;
  }
}

static void sfe3_tabulate_gl_prism(::coder::SizeType etype,
                                   const ::coder::array<real_T, 2U> &cs,
                                   ::coder::array<real_T, 2U> &sfvals,
                                   ::coder::array<real_T, 3U> &sdvals)
{
  real_T tmp_data[1029];
  real_T dv4[378];
  real_T dv1[225];
  real_T dv5[126];
  real_T dv[120];
  real_T dv3[75];
  real_T dv2[40];
  ::coder::SizeType i;
  ::coder::SizeType i1;
  ::coder::SizeType i2;
  ::coder::SizeType i3;
  ::coder::SizeType loop_ub;
  ::coder::SizeType tmp_size_idx_1;
  ::coder::SizeType tmp_size_idx_2;
  ::coder::SizeType ub_loop;
  int16_T unnamed_idx_1;
  int16_T unnamed_idx_2;
  //  prisms
  ub_loop = iv[etype - 1];
  sfvals.set_size(cs.size(0), ub_loop);
  sdvals.set_size(cs.size(0), ub_loop, cs.size(1));
  switch (etype) {
  case 205:
    //  Serial mode
    ub_loop = cs.size(0) - 1;
    for (::coder::SizeType q = 0; q <= ub_loop; q++) {
      prism_gl_40(cs[cs.size(1) * q], cs[cs.size(1) * q + 1],
                  cs[cs.size(1) * q + 2], dv2, dv);
      loop_ub = sfvals.size(1);
      for (i = 0; i < loop_ub; i++) {
        sfvals[i + sfvals.size(1) * q] = dv2[i];
      }
      unnamed_idx_1 = static_cast<int16_T>(sdvals.size(1));
      unnamed_idx_2 = static_cast<int16_T>(sdvals.size(2));
      i = 0;
      i1 = 0;
      loop_ub = 0;
      i2 = 0;
      tmp_size_idx_2 = sdvals.size(2);
      tmp_size_idx_1 = sdvals.size(1);
      for (i3 = 0; i3 < unnamed_idx_1 * unnamed_idx_2; i3++) {
        tmp_data[i1 + tmp_size_idx_2 * i] = dv[i2 + 3 * loop_ub];
        loop_ub++;
        i++;
        if (i > tmp_size_idx_1 - 1) {
          i = 0;
          i1++;
        }
        if (loop_ub > 39) {
          loop_ub = 0;
          i2++;
        }
      }
      for (i = 0; i < tmp_size_idx_1; i++) {
        for (i1 = 0; i1 < tmp_size_idx_2; i1++) {
          sdvals[(i1 + sdvals.size(2) * i) +
                 sdvals.size(2) * sdvals.size(1) * q] =
              tmp_data[i1 + tmp_size_idx_2 * i];
        }
      }
    }
    break;
  case 209:
    //  Serial mode
    ub_loop = cs.size(0) - 1;
    for (::coder::SizeType q = 0; q <= ub_loop; q++) {
      prism_gl_75(cs[cs.size(1) * q], cs[cs.size(1) * q + 1],
                  cs[cs.size(1) * q + 2], dv3, dv1);
      loop_ub = sfvals.size(1);
      for (i = 0; i < loop_ub; i++) {
        sfvals[i + sfvals.size(1) * q] = dv3[i];
      }
      unnamed_idx_1 = static_cast<int16_T>(sdvals.size(1));
      unnamed_idx_2 = static_cast<int16_T>(sdvals.size(2));
      i = 0;
      i1 = 0;
      loop_ub = 0;
      i2 = 0;
      tmp_size_idx_2 = sdvals.size(2);
      tmp_size_idx_1 = sdvals.size(1);
      for (i3 = 0; i3 < unnamed_idx_1 * unnamed_idx_2; i3++) {
        tmp_data[i1 + tmp_size_idx_2 * i] = dv1[i2 + 3 * loop_ub];
        loop_ub++;
        i++;
        if (i > tmp_size_idx_1 - 1) {
          i = 0;
          i1++;
        }
        if (loop_ub > 74) {
          loop_ub = 0;
          i2++;
        }
      }
      for (i = 0; i < tmp_size_idx_1; i++) {
        for (i1 = 0; i1 < tmp_size_idx_2; i1++) {
          sdvals[(i1 + sdvals.size(2) * i) +
                 sdvals.size(2) * sdvals.size(1) * q] =
              tmp_data[i1 + tmp_size_idx_2 * i];
        }
      }
    }
    break;
  default:
    m2cAssert(etype == 213, "Gauss-Lobatto only supports up to quintic.");
    //  Serial mode
    ub_loop = cs.size(0) - 1;
    for (::coder::SizeType q = 0; q <= ub_loop; q++) {
      prism_gl_126(cs[cs.size(1) * q], cs[cs.size(1) * q + 1],
                   cs[cs.size(1) * q + 2], dv5, dv4);
      loop_ub = sfvals.size(1);
      for (i = 0; i < loop_ub; i++) {
        sfvals[i + sfvals.size(1) * q] = dv5[i];
      }
      unnamed_idx_1 = static_cast<int16_T>(sdvals.size(1));
      unnamed_idx_2 = static_cast<int16_T>(sdvals.size(2));
      i = 0;
      i1 = 0;
      loop_ub = 0;
      i2 = 0;
      tmp_size_idx_2 = sdvals.size(2);
      tmp_size_idx_1 = sdvals.size(1);
      for (i3 = 0; i3 < unnamed_idx_1 * unnamed_idx_2; i3++) {
        tmp_data[i1 + tmp_size_idx_2 * i] = dv4[i2 + 3 * loop_ub];
        loop_ub++;
        i++;
        if (i > tmp_size_idx_1 - 1) {
          i = 0;
          i1++;
        }
        if (loop_ub > 125) {
          loop_ub = 0;
          i2++;
        }
      }
      for (i = 0; i < tmp_size_idx_1; i++) {
        for (i1 = 0; i1 < tmp_size_idx_2; i1++) {
          sdvals[(i1 + sdvals.size(2) * i) +
                 sdvals.size(2) * sdvals.size(1) * q] =
              tmp_data[i1 + tmp_size_idx_2 * i];
        }
      }
    }
    break;
  }
}

static void sfe3_tabulate_gl_pyra(::coder::SizeType etype,
                                  const ::coder::array<real_T, 2U> &cs,
                                  ::coder::array<real_T, 2U> &sfvals,
                                  ::coder::array<real_T, 3U> &sdvals)
{
  real_T tmp_data[1029];
  real_T dv2[165];
  real_T dv[90];
  real_T dv3[55];
  real_T dv1[30];
  ::coder::SizeType i;
  ::coder::SizeType i1;
  ::coder::SizeType i2;
  ::coder::SizeType i3;
  ::coder::SizeType loop_ub;
  ::coder::SizeType tmp_size_idx_1;
  ::coder::SizeType tmp_size_idx_2;
  ::coder::SizeType ub_loop;
  int16_T unnamed_idx_1;
  int16_T unnamed_idx_2;
  //  pyra
  ub_loop = iv[etype - 1];
  sfvals.set_size(cs.size(0), ub_loop);
  sdvals.set_size(cs.size(0), ub_loop, cs.size(1));
  if (etype == 173) {
    //  Serial mode
    ub_loop = cs.size(0) - 1;
    for (::coder::SizeType q = 0; q <= ub_loop; q++) {
      pyra_gl_30(cs[cs.size(1) * q], cs[cs.size(1) * q + 1],
                 cs[cs.size(1) * q + 2], dv1, dv);
      loop_ub = sfvals.size(1);
      for (i = 0; i < loop_ub; i++) {
        sfvals[i + sfvals.size(1) * q] = dv1[i];
      }
      unnamed_idx_1 = static_cast<int16_T>(sdvals.size(1));
      unnamed_idx_2 = static_cast<int16_T>(sdvals.size(2));
      i = 0;
      i1 = 0;
      loop_ub = 0;
      i2 = 0;
      tmp_size_idx_2 = sdvals.size(2);
      tmp_size_idx_1 = sdvals.size(1);
      for (i3 = 0; i3 < unnamed_idx_1 * unnamed_idx_2; i3++) {
        tmp_data[i1 + tmp_size_idx_2 * i] = dv[i2 + 3 * loop_ub];
        loop_ub++;
        i++;
        if (i > tmp_size_idx_1 - 1) {
          i = 0;
          i1++;
        }
        if (loop_ub > 29) {
          loop_ub = 0;
          i2++;
        }
      }
      for (i = 0; i < tmp_size_idx_1; i++) {
        for (i1 = 0; i1 < tmp_size_idx_2; i1++) {
          sdvals[(i1 + sdvals.size(2) * i) +
                 sdvals.size(2) * sdvals.size(1) * q] =
              tmp_data[i1 + tmp_size_idx_2 * i];
        }
      }
    }
  } else {
    m2cAssert(etype == 177, "Pyramid only support up to quartic");
    //  Serial mode
    ub_loop = cs.size(0) - 1;
    for (::coder::SizeType q = 0; q <= ub_loop; q++) {
      pyra_gl_55(cs[cs.size(1) * q], cs[cs.size(1) * q + 1],
                 cs[cs.size(1) * q + 2], dv3, dv2);
      loop_ub = sfvals.size(1);
      for (i = 0; i < loop_ub; i++) {
        sfvals[i + sfvals.size(1) * q] = dv3[i];
      }
      unnamed_idx_1 = static_cast<int16_T>(sdvals.size(1));
      unnamed_idx_2 = static_cast<int16_T>(sdvals.size(2));
      i = 0;
      i1 = 0;
      loop_ub = 0;
      i2 = 0;
      tmp_size_idx_2 = sdvals.size(2);
      tmp_size_idx_1 = sdvals.size(1);
      for (i3 = 0; i3 < unnamed_idx_1 * unnamed_idx_2; i3++) {
        tmp_data[i1 + tmp_size_idx_2 * i] = dv2[i2 + 3 * loop_ub];
        loop_ub++;
        i++;
        if (i > tmp_size_idx_1 - 1) {
          i = 0;
          i1++;
        }
        if (loop_ub > 54) {
          loop_ub = 0;
          i2++;
        }
      }
      for (i = 0; i < tmp_size_idx_1; i++) {
        for (i1 = 0; i1 < tmp_size_idx_2; i1++) {
          sdvals[(i1 + sdvals.size(2) * i) +
                 sdvals.size(2) * sdvals.size(1) * q] =
              tmp_data[i1 + tmp_size_idx_2 * i];
        }
      }
    }
  }
}

static void sfe3_tabulate_gl_tet(::coder::SizeType etype,
                                 const ::coder::array<real_T, 2U> &cs,
                                 ::coder::array<real_T, 2U> &sfvals,
                                 ::coder::array<real_T, 3U> &sdvals)
{
  real_T tmp_data[1029];
  real_T dv2[105];
  real_T dv[60];
  real_T dv3[35];
  real_T dv1[20];
  ::coder::SizeType i;
  ::coder::SizeType i1;
  ::coder::SizeType i2;
  ::coder::SizeType i3;
  ::coder::SizeType loop_ub;
  ::coder::SizeType tmp_size_idx_1;
  ::coder::SizeType tmp_size_idx_2;
  ::coder::SizeType ub_loop;
  int16_T unnamed_idx_1;
  int16_T unnamed_idx_2;
  //  tet
  ub_loop = iv[etype - 1];
  sfvals.set_size(cs.size(0), ub_loop);
  sdvals.set_size(cs.size(0), ub_loop, cs.size(1));
  if (etype == 141) {
    //  Serial mode
    ub_loop = cs.size(0) - 1;
    for (::coder::SizeType q = 0; q <= ub_loop; q++) {
      tet_gl_20(cs[cs.size(1) * q], cs[cs.size(1) * q + 1],
                cs[cs.size(1) * q + 2], dv1, dv);
      loop_ub = sfvals.size(1);
      for (i = 0; i < loop_ub; i++) {
        sfvals[i + sfvals.size(1) * q] = dv1[i];
      }
      unnamed_idx_1 = static_cast<int16_T>(sdvals.size(1));
      unnamed_idx_2 = static_cast<int16_T>(sdvals.size(2));
      i = 0;
      i1 = 0;
      loop_ub = 0;
      i2 = 0;
      tmp_size_idx_2 = sdvals.size(2);
      tmp_size_idx_1 = sdvals.size(1);
      for (i3 = 0; i3 < unnamed_idx_1 * unnamed_idx_2; i3++) {
        tmp_data[i1 + tmp_size_idx_2 * i] = dv[i2 + 3 * loop_ub];
        loop_ub++;
        i++;
        if (i > tmp_size_idx_1 - 1) {
          i = 0;
          i1++;
        }
        if (loop_ub > 19) {
          loop_ub = 0;
          i2++;
        }
      }
      for (i = 0; i < tmp_size_idx_1; i++) {
        for (i1 = 0; i1 < tmp_size_idx_2; i1++) {
          sdvals[(i1 + sdvals.size(2) * i) +
                 sdvals.size(2) * sdvals.size(1) * q] =
              tmp_data[i1 + tmp_size_idx_2 * i];
        }
      }
    }
  } else {
    m2cAssert(
        etype == 145,
        "Gauss-Lobatto tetrahedral elements are supported only up to quartic");
    //  Serial mode
    ub_loop = cs.size(0) - 1;
    for (::coder::SizeType q = 0; q <= ub_loop; q++) {
      tet_gl_35(cs[cs.size(1) * q], cs[cs.size(1) * q + 1],
                cs[cs.size(1) * q + 2], dv3, dv2);
      loop_ub = sfvals.size(1);
      for (i = 0; i < loop_ub; i++) {
        sfvals[i + sfvals.size(1) * q] = dv3[i];
      }
      unnamed_idx_1 = static_cast<int16_T>(sdvals.size(1));
      unnamed_idx_2 = static_cast<int16_T>(sdvals.size(2));
      i = 0;
      i1 = 0;
      loop_ub = 0;
      i2 = 0;
      tmp_size_idx_2 = sdvals.size(2);
      tmp_size_idx_1 = sdvals.size(1);
      for (i3 = 0; i3 < unnamed_idx_1 * unnamed_idx_2; i3++) {
        tmp_data[i1 + tmp_size_idx_2 * i] = dv2[i2 + 3 * loop_ub];
        loop_ub++;
        i++;
        if (i > tmp_size_idx_1 - 1) {
          i = 0;
          i1++;
        }
        if (loop_ub > 34) {
          loop_ub = 0;
          i2++;
        }
      }
      for (i = 0; i < tmp_size_idx_1; i++) {
        for (i1 = 0; i1 < tmp_size_idx_2; i1++) {
          sdvals[(i1 + sdvals.size(2) * i) +
                 sdvals.size(2) * sdvals.size(1) * q] =
              tmp_data[i1 + tmp_size_idx_2 * i];
        }
      }
    }
  }
}

// sfe3_tabulate_shapefuncs - Tabulate shape functions and sdvals at qpoints
static void sfe3_tabulate_shapefuncs(::coder::SizeType etype,
                                     const ::coder::array<real_T, 2U> &cs,
                                     ::coder::array<real_T, 2U> &sfvals,
                                     ::coder::array<real_T, 3U> &sdvals)
{
  ::coder::SizeType postype;
  postype = etype & 3;
  if (postype == 0) {
    ::coder::SizeType i;
    i = etype >> 5 & 7;
    if (i == 4) {
      sfe3_tabulate_equi_tet(etype, cs, sfvals, sdvals);
    } else if (i == 5) {
      sfe3_tabulate_equi_pyra(etype, cs, sfvals, sdvals);
    } else if (i == 6) {
      sfe3_tabulate_equi_prism(etype, cs, sfvals, sdvals);
    } else {
      sfe3_tabulate_equi_hexa(etype, cs, sfvals, sdvals);
    }
  } else {
    ::coder::SizeType i;
    m2cAssert(postype == 1,
              "Only supports Equidistant and Gauss-Lobatto points in 3D");
    i = etype >> 5 & 7;
    if (i == 4) {
      sfe3_tabulate_gl_tet(etype, cs, sfvals, sdvals);
    } else if (i == 5) {
      sfe3_tabulate_gl_pyra(etype, cs, sfvals, sdvals);
    } else if (i == 6) {
      sfe3_tabulate_gl_prism(etype, cs, sfvals, sdvals);
    } else {
      sfe3_tabulate_gl_hexa(etype, cs, sfvals, sdvals);
    }
  }
}

// sfe_init - Initialize/reinitialize an sfe object for non-boundary element
static void sfe_init(SfeObject *sfe, ::coder::SizeType etypes,
                     const ::coder::array<real_T, 2U> &xs)
{
  real_T dv[9];
  real_T v;
  ::coder::SizeType i;
  ::coder::SizeType loop_ub;
  ::coder::SizeType sfe_idx_0_tmp_tmp;
  ::coder::SizeType shape;
  ::coder::SizeType topo_dim;
  shape = etypes >> 5;
  topo_dim = ((shape > 0) + (shape > 1)) + (shape > 3);
  //  Geometric dimension
  if (xs.size(1) < topo_dim) {
    m2cErrMsgIdAndTxt("sfe_init:badDim",
                      "geometric dim cannot be smaller than topo dim");
  }
  sfe->geom_dim = xs.size(1);
  //  assign geom dimension
  sfe->topo_dim = topo_dim;
  //  assign topo dimension
  m2cAssert(iv[etypes - 1] == xs.size(0), "nnodes do not match");
  sfe->etypes[0] = etypes;
  sfe->etypes[1] = etypes;
  //  Get number of nodes per element
  i = iv[etypes - 1];
  sfe->nnodes[0] = i;
  sfe->nnodes[1] = i;
  //  Set up quadrature
  tabulate_quadratures(etypes,
                       (((etypes >> 2 & 7) << 1) + ((etypes >> 2 & 7) > 1)) +
                           (xs.size(1) > topo_dim),
                       sfe->cs, sfe->ws);
  sfe->nqp = sfe->ws.size(0);
  //  Solution space shape functions & derivs
  tabulate_shapefuncs(etypes, sfe->cs, sfe->shapes_sol, sfe->derivs_sol);
  //  Geometry space shape functions & derivs
  loop_ub = sfe->shapes_sol.size(1) * sfe->shapes_sol.size(0);
  sfe->shapes_geom.set_size(sfe->shapes_sol.size(0), sfe->shapes_sol.size(1));
  for (i = 0; i < loop_ub; i++) {
    sfe->shapes_geom[i] = sfe->shapes_sol[i];
  }
  loop_ub = sfe->derivs_sol.size(2) * sfe->derivs_sol.size(1) *
            sfe->derivs_sol.size(0);
  sfe->derivs_geom.set_size(sfe->derivs_sol.size(0), sfe->derivs_sol.size(1),
                            sfe->derivs_sol.size(2));
  for (i = 0; i < loop_ub; i++) {
    sfe->derivs_geom[i] = sfe->derivs_sol[i];
  }
  //  potentially skip re-tabulating
  sfe_idx_0_tmp_tmp = sfe->nqp;
  sfe->cs_phy.set_size(sfe_idx_0_tmp_tmp, xs.size(1));
  for (::coder::SizeType q{0}; q < sfe_idx_0_tmp_tmp; q++) {
    i = xs.size(1);
    for (::coder::SizeType k{0}; k < i; k++) {
      ::coder::SizeType m;
      m = sfe->shapes_sol.size(1);
      v = sfe->shapes_sol[sfe->shapes_sol.size(1) * q] * xs[k];
      for (::coder::SizeType b_i{2}; b_i <= m; b_i++) {
        v += sfe->shapes_sol[(b_i + sfe->shapes_sol.size(1) * q) - 1] *
             xs[k + xs.size(1) * (b_i - 1)];
      }
      sfe->cs_phy[k + sfe->cs_phy.size(1) * q] = v;
    }
  }
  //  Compute Jacobian
  sfe->wdetJ.set_size(sfe->nqp);
  if ((etypes == 68) || (etypes == 132) || (etypes == 36)) {
    real_T d;
    ::coder::SizeType geom_dim;
    ::coder::SizeType n;
    //  A single Jacobian matrix (transpose) is needed for simplex elements
    geom_dim = xs.size(1);
    topo_dim = sfe->derivs_sol.size(2);
    std::memset(&dv[0], 0, 9U * sizeof(real_T));
    n = xs.size(0);
    for (::coder::SizeType k{0}; k < n; k++) {
      for (::coder::SizeType b_i{0}; b_i < topo_dim; b_i++) {
        for (::coder::SizeType j{0}; j < geom_dim; j++) {
          i = j + 3 * b_i;
          dv[i] += xs[j + xs.size(1) * k] *
                   sfe->derivs_sol[b_i + sfe->derivs_sol.size(2) * k];
        }
      }
    }
    if (xs.size(1) == sfe->derivs_sol.size(2)) {
      if (xs.size(1) == 1) {
        d = dv[0];
      } else if (xs.size(1) == 2) {
        d = dv[0] * dv[4] - dv[1] * dv[3];
      } else {
        d = (dv[2] * (dv[3] * dv[7] - dv[4] * dv[6]) +
             dv[5] * (dv[1] * dv[6] - dv[0] * dv[7])) +
            dv[8] * (dv[0] * dv[4] - dv[1] * dv[3]);
      }
    } else if (sfe->derivs_sol.size(2) == 1) {
      d = dv[0] * dv[0] + dv[1] * dv[1];
      if (xs.size(1) == 3) {
        d += dv[2] * dv[2];
      }
      d = std::sqrt(d);
    } else {
      //  must be 2x3
      dv[6] = dv[1] * dv[5] - dv[2] * dv[4];
      dv[7] = dv[2] * dv[3] - dv[0] * dv[5];
      dv[8] = dv[0] * dv[4] - dv[1] * dv[3];
      d = std::sqrt((dv[6] * dv[6] + dv[7] * dv[7]) + dv[8] * dv[8]);
    }
    sfe->jacTs.set_size(3, 3);
    for (i = 0; i < 9; i++) {
      sfe->jacTs[i] = dv[i];
    }
    for (::coder::SizeType q{0}; q < sfe_idx_0_tmp_tmp; q++) {
      sfe->wdetJ[q] = d * sfe->ws[q];
    }
  } else {
    //  Super-parametric
    loop_ub = sfe->nqp * 3;
    sfe->jacTs.set_size(loop_ub, 3);
    for (::coder::SizeType q{0}; q < sfe_idx_0_tmp_tmp; q++) {
      ::coder::SizeType geom_dim;
      ::coder::SizeType n;
      ::coder::SizeType y;
      y = q * 3;
      geom_dim = xs.size(1);
      topo_dim = sfe->derivs_sol.size(2);
      std::memset(&dv[0], 0, 9U * sizeof(real_T));
      n = xs.size(0);
      for (::coder::SizeType k{0}; k < n; k++) {
        for (::coder::SizeType b_i{0}; b_i < topo_dim; b_i++) {
          for (::coder::SizeType j{0}; j < geom_dim; j++) {
            i = j + 3 * b_i;
            dv[i] += xs[j + xs.size(1) * k] *
                     sfe->derivs_sol[(b_i + sfe->derivs_sol.size(2) * k) +
                                     sfe->derivs_sol.size(2) *
                                         sfe->derivs_sol.size(1) * q];
          }
        }
      }
      if (xs.size(1) == sfe->derivs_sol.size(2)) {
        if (xs.size(1) == 1) {
          v = dv[0];
        } else if (xs.size(1) == 2) {
          v = dv[0] * dv[4] - dv[1] * dv[3];
        } else {
          v = (dv[2] * (dv[3] * dv[7] - dv[4] * dv[6]) +
               dv[5] * (dv[1] * dv[6] - dv[0] * dv[7])) +
              dv[8] * (dv[0] * dv[4] - dv[1] * dv[3]);
        }
      } else if (sfe->derivs_sol.size(2) == 1) {
        v = dv[0] * dv[0] + dv[1] * dv[1];
        if (xs.size(1) == 3) {
          v += dv[2] * dv[2];
        }
        v = std::sqrt(v);
      } else {
        //  must be 2x3
        dv[6] = dv[1] * dv[5] - dv[2] * dv[4];
        dv[7] = dv[2] * dv[3] - dv[0] * dv[5];
        dv[8] = dv[0] * dv[4] - dv[1] * dv[3];
        v = std::sqrt((dv[6] * dv[6] + dv[7] * dv[7]) + dv[8] * dv[8]);
      }
      for (i = 0; i < 3; i++) {
        loop_ub = i + y;
        sfe->jacTs[3 * loop_ub] = dv[3 * i];
        sfe->jacTs[3 * loop_ub + 1] = dv[3 * i + 1];
        sfe->jacTs[3 * loop_ub + 2] = dv[3 * i + 2];
      }
      sfe->wdetJ[q] = v;
      sfe->wdetJ[q] = sfe->wdetJ[q] * sfe->ws[q];
    }
  }
}

// sfe_init - Initialize/reinitialize an sfe object for non-boundary element
static void sfe_init(SfeObject *sfe, const ::coder::array<real_T, 2U> &xs)
{
  real_T dv[9];
  real_T v;
  ::coder::SizeType i;
  ::coder::SizeType sfe_idx_0_tmp_tmp;
  boolean_T cond;
  if ((sfe->etypes[0] > 0) && (iv[sfe->etypes[0] - 1] != 0)) {
    cond = true;
  } else {
    cond = false;
  }
  m2cAssert(cond, "");
  //  potentially skip re-tabulating
  sfe_idx_0_tmp_tmp = sfe->nqp;
  sfe->cs_phy.set_size(sfe_idx_0_tmp_tmp, xs.size(1));
  for (::coder::SizeType q{0}; q < sfe_idx_0_tmp_tmp; q++) {
    i = xs.size(1);
    for (::coder::SizeType k{0}; k < i; k++) {
      ::coder::SizeType m;
      m = sfe->shapes_geom.size(1);
      v = sfe->shapes_geom[sfe->shapes_geom.size(1) * q] * xs[k];
      for (::coder::SizeType b_i{2}; b_i <= m; b_i++) {
        v += sfe->shapes_geom[(b_i + sfe->shapes_geom.size(1) * q) - 1] *
             xs[k + xs.size(1) * (b_i - 1)];
      }
      sfe->cs_phy[k + sfe->cs_phy.size(1) * q] = v;
    }
  }
  //  Compute Jacobian
  sfe->wdetJ.set_size(sfe->nqp);
  if ((sfe->etypes[1] == 68) || (sfe->etypes[1] == 132) ||
      (sfe->etypes[1] == 36)) {
    real_T d;
    ::coder::SizeType geom_dim;
    ::coder::SizeType n;
    ::coder::SizeType topo_dim;
    //  A single Jacobian matrix (transpose) is needed for simplex elements
    geom_dim = xs.size(1);
    topo_dim = sfe->derivs_geom.size(2);
    std::memset(&dv[0], 0, 9U * sizeof(real_T));
    n = xs.size(0);
    for (::coder::SizeType k{0}; k < n; k++) {
      for (::coder::SizeType b_i{0}; b_i < topo_dim; b_i++) {
        for (::coder::SizeType j{0}; j < geom_dim; j++) {
          i = j + 3 * b_i;
          dv[i] += xs[j + xs.size(1) * k] *
                   sfe->derivs_geom[b_i + sfe->derivs_geom.size(2) * k];
        }
      }
    }
    if (xs.size(1) == sfe->derivs_geom.size(2)) {
      if (xs.size(1) == 1) {
        d = dv[0];
      } else if (xs.size(1) == 2) {
        d = dv[0] * dv[4] - dv[1] * dv[3];
      } else {
        d = (dv[2] * (dv[3] * dv[7] - dv[4] * dv[6]) +
             dv[5] * (dv[1] * dv[6] - dv[0] * dv[7])) +
            dv[8] * (dv[0] * dv[4] - dv[1] * dv[3]);
      }
    } else if (sfe->derivs_geom.size(2) == 1) {
      d = dv[0] * dv[0] + dv[1] * dv[1];
      if (xs.size(1) == 3) {
        d += dv[2] * dv[2];
      }
      d = std::sqrt(d);
    } else {
      //  must be 2x3
      dv[6] = dv[1] * dv[5] - dv[2] * dv[4];
      dv[7] = dv[2] * dv[3] - dv[0] * dv[5];
      dv[8] = dv[0] * dv[4] - dv[1] * dv[3];
      d = std::sqrt((dv[6] * dv[6] + dv[7] * dv[7]) + dv[8] * dv[8]);
    }
    sfe->jacTs.set_size(3, 3);
    for (i = 0; i < 9; i++) {
      sfe->jacTs[i] = dv[i];
    }
    for (::coder::SizeType q{0}; q < sfe_idx_0_tmp_tmp; q++) {
      sfe->wdetJ[q] = d * sfe->ws[q];
    }
  } else {
    ::coder::SizeType sfe_idx_0;
    //  Super-parametric
    sfe_idx_0 = sfe->nqp * 3;
    sfe->jacTs.set_size(sfe_idx_0, 3);
    for (::coder::SizeType q{0}; q < sfe_idx_0_tmp_tmp; q++) {
      ::coder::SizeType geom_dim;
      ::coder::SizeType n;
      ::coder::SizeType topo_dim;
      ::coder::SizeType y;
      y = q * 3;
      geom_dim = xs.size(1);
      topo_dim = sfe->derivs_geom.size(2);
      std::memset(&dv[0], 0, 9U * sizeof(real_T));
      n = xs.size(0);
      for (::coder::SizeType k{0}; k < n; k++) {
        for (::coder::SizeType b_i{0}; b_i < topo_dim; b_i++) {
          for (::coder::SizeType j{0}; j < geom_dim; j++) {
            i = j + 3 * b_i;
            dv[i] += xs[j + xs.size(1) * k] *
                     sfe->derivs_geom[(b_i + sfe->derivs_geom.size(2) * k) +
                                      sfe->derivs_geom.size(2) *
                                          sfe->derivs_geom.size(1) * q];
          }
        }
      }
      if (xs.size(1) == sfe->derivs_geom.size(2)) {
        if (xs.size(1) == 1) {
          v = dv[0];
        } else if (xs.size(1) == 2) {
          v = dv[0] * dv[4] - dv[1] * dv[3];
        } else {
          v = (dv[2] * (dv[3] * dv[7] - dv[4] * dv[6]) +
               dv[5] * (dv[1] * dv[6] - dv[0] * dv[7])) +
              dv[8] * (dv[0] * dv[4] - dv[1] * dv[3]);
        }
      } else if (sfe->derivs_geom.size(2) == 1) {
        v = dv[0] * dv[0] + dv[1] * dv[1];
        if (xs.size(1) == 3) {
          v += dv[2] * dv[2];
        }
        v = std::sqrt(v);
      } else {
        //  must be 2x3
        dv[6] = dv[1] * dv[5] - dv[2] * dv[4];
        dv[7] = dv[2] * dv[3] - dv[0] * dv[5];
        dv[8] = dv[0] * dv[4] - dv[1] * dv[3];
        v = std::sqrt((dv[6] * dv[6] + dv[7] * dv[7]) + dv[8] * dv[8]);
      }
      for (i = 0; i < 3; i++) {
        sfe_idx_0 = i + y;
        sfe->jacTs[3 * sfe_idx_0] = dv[3 * i];
        sfe->jacTs[3 * sfe_idx_0 + 1] = dv[3 * i + 1];
        sfe->jacTs[3 * sfe_idx_0 + 2] = dv[3 * i + 2];
      }
      sfe->wdetJ[q] = v;
      sfe->wdetJ[q] = sfe->wdetJ[q] * sfe->ws[q];
    }
  }
}

// tabulate_quadratures - Tabulate quadrature rule for given element type
static void tabulate_quadratures(::coder::SizeType etype, ::coder::SizeType qd,
                                 ::coder::array<real_T, 2U> &cs,
                                 ::coder::array<real_T, 1U> &ws)
{
  ::coder::SizeType shape;
  shape = etype >> 5 & 7;
  switch (shape) {
  case 1:
    bar_quadrules(qd, cs, ws);
    break;
  case 2:
    tri_quadrules(qd, cs, ws);
    break;
  case 3: {
    if (qd <= 1) {
      ::coder::SizeType nqp;
      nqp = ::sfe_qrules::quad_deg1_qrule();
      cs.set_size(nqp, 2);
      ws.set_size(nqp);
      ::sfe_qrules::quad_deg1_qrule(&cs[0], &(ws.data())[0]);
    } else if (qd <= 3) {
      ::coder::SizeType nqp;
      nqp = ::sfe_qrules::quad_deg3_qrule();
      cs.set_size(nqp, 2);
      ws.set_size(nqp);
      ::sfe_qrules::quad_deg3_qrule(&cs[0], &(ws.data())[0]);
    } else if (qd <= 5) {
      ::coder::SizeType nqp;
      nqp = ::sfe_qrules::quad_deg5_qrule();
      cs.set_size(nqp, 2);
      ws.set_size(nqp);
      ::sfe_qrules::quad_deg5_qrule(&cs[0], &(ws.data())[0]);
    } else if (qd <= 7) {
      ::coder::SizeType nqp;
      nqp = ::sfe_qrules::quad_deg7_qrule();
      cs.set_size(nqp, 2);
      ws.set_size(nqp);
      ::sfe_qrules::quad_deg7_qrule(&cs[0], &(ws.data())[0]);
    } else if (qd <= 9) {
      ::coder::SizeType nqp;
      nqp = ::sfe_qrules::quad_deg9_qrule();
      cs.set_size(nqp, 2);
      ws.set_size(nqp);
      ::sfe_qrules::quad_deg9_qrule(&cs[0], &(ws.data())[0]);
    } else if (qd <= 11) {
      ::coder::SizeType nqp;
      nqp = ::sfe_qrules::quad_deg11_qrule();
      cs.set_size(nqp, 2);
      ws.set_size(nqp);
      ::sfe_qrules::quad_deg11_qrule(&cs[0], &(ws.data())[0]);
    } else {
      ::coder::SizeType nqp;
      if (qd > 13) {
        m2cWarnMsgIdAndTxt("bar_quadrules:UnsupportedDegree",
                           "Only support up to degree 13");
      }
      nqp = ::sfe_qrules::quad_deg13_qrule();
      cs.set_size(nqp, 2);
      ws.set_size(nqp);
      ::sfe_qrules::quad_deg13_qrule(&cs[0], &(ws.data())[0]);
    }
  } break;
  case 4: {
    if (qd <= 1) {
      ::coder::SizeType nqp;
      nqp = ::sfe_qrules::tet_deg1_qrule();
      cs.set_size(nqp, 3);
      ws.set_size(nqp);
      ::sfe_qrules::tet_deg1_qrule(&cs[0], &(ws.data())[0]);
    } else if (qd <= 2) {
      ::coder::SizeType nqp;
      nqp = ::sfe_qrules::tet_deg2_qrule();
      cs.set_size(nqp, 3);
      ws.set_size(nqp);
      ::sfe_qrules::tet_deg2_qrule(&cs[0], &(ws.data())[0]);
    } else if (qd <= 3) {
      ::coder::SizeType nqp;
      nqp = ::sfe_qrules::tet_deg3_qrule();
      cs.set_size(nqp, 3);
      ws.set_size(nqp);
      ::sfe_qrules::tet_deg3_qrule(&cs[0], &(ws.data())[0]);
    } else if (qd <= 5) {
      ::coder::SizeType nqp;
      nqp = ::sfe_qrules::tet_deg5_qrule();
      cs.set_size(nqp, 3);
      ws.set_size(nqp);
      ::sfe_qrules::tet_deg5_qrule(&cs[0], &(ws.data())[0]);
    } else if (qd <= 6) {
      ::coder::SizeType nqp;
      nqp = ::sfe_qrules::tet_deg6_qrule();
      cs.set_size(nqp, 3);
      ws.set_size(nqp);
      ::sfe_qrules::tet_deg6_qrule(&cs[0], &(ws.data())[0]);
    } else if (qd <= 7) {
      ::coder::SizeType nqp;
      nqp = ::sfe_qrules::tet_deg7_qrule();
      cs.set_size(nqp, 3);
      ws.set_size(nqp);
      ::sfe_qrules::tet_deg7_qrule(&cs[0], &(ws.data())[0]);
    } else if (qd <= 9) {
      ::coder::SizeType nqp;
      //  The degree-8 quadrature rule KEAST9 has negative weights, so we do
      nqp = ::sfe_qrules::tet_deg9_qrule();
      cs.set_size(nqp, 3);
      ws.set_size(nqp);
      ::sfe_qrules::tet_deg9_qrule(&cs[0], &(ws.data())[0]);
    } else if (qd <= 11) {
      ::coder::SizeType nqp;
      nqp = ::sfe_qrules::tet_deg11_qrule();
      cs.set_size(nqp, 3);
      ws.set_size(nqp);
      ::sfe_qrules::tet_deg11_qrule(&cs[0], &(ws.data())[0]);
    } else {
      ::coder::SizeType nqp;
      if (qd > 13) {
        m2cWarnMsgIdAndTxt("tet_quadrules:UnsupportedDegree",
                           "Only support up to degree 13");
      }
      nqp = ::sfe_qrules::tet_deg13_qrule();
      cs.set_size(nqp, 3);
      ws.set_size(nqp);
      ::sfe_qrules::tet_deg13_qrule(&cs[0], &(ws.data())[0]);
    }
  } break;
  case 7: {
    if (qd <= 1) {
      ::coder::SizeType nqp;
      nqp = ::sfe_qrules::hexa_deg1_qrule();
      cs.set_size(nqp, 3);
      ws.set_size(nqp);
      ::sfe_qrules::hexa_deg1_qrule(&cs[0], &(ws.data())[0]);
    } else if (qd <= 3) {
      ::coder::SizeType nqp;
      nqp = ::sfe_qrules::hexa_deg3_qrule();
      cs.set_size(nqp, 3);
      ws.set_size(nqp);
      ::sfe_qrules::hexa_deg3_qrule(&cs[0], &(ws.data())[0]);
    } else if (qd <= 5) {
      ::coder::SizeType nqp;
      nqp = ::sfe_qrules::hexa_deg5_qrule();
      cs.set_size(nqp, 3);
      ws.set_size(nqp);
      ::sfe_qrules::hexa_deg5_qrule(&cs[0], &(ws.data())[0]);
    } else if (qd <= 7) {
      ::coder::SizeType nqp;
      nqp = ::sfe_qrules::hexa_deg7_qrule();
      cs.set_size(nqp, 3);
      ws.set_size(nqp);
      ::sfe_qrules::hexa_deg7_qrule(&cs[0], &(ws.data())[0]);
    } else if (qd <= 9) {
      ::coder::SizeType nqp;
      nqp = ::sfe_qrules::hexa_deg9_qrule();
      cs.set_size(nqp, 3);
      ws.set_size(nqp);
      ::sfe_qrules::hexa_deg9_qrule(&cs[0], &(ws.data())[0]);
    } else if (qd <= 11) {
      ::coder::SizeType nqp;
      nqp = ::sfe_qrules::hexa_deg11_qrule();
      cs.set_size(nqp, 3);
      ws.set_size(nqp);
      ::sfe_qrules::hexa_deg11_qrule(&cs[0], &(ws.data())[0]);
    } else {
      ::coder::SizeType nqp;
      if (qd > 13) {
        m2cWarnMsgIdAndTxt("hexa_quadrules:UnsupportedDegree",
                           "Only support up to degree 13");
      }
      nqp = ::sfe_qrules::hexa_deg13_qrule();
      cs.set_size(nqp, 3);
      ws.set_size(nqp);
      ::sfe_qrules::hexa_deg13_qrule(&cs[0], &(ws.data())[0]);
    }
  } break;
  case 6: {
    if (qd <= 1) {
      ::coder::SizeType nqp;
      nqp = ::sfe_qrules::prism_deg1_qrule();
      cs.set_size(nqp, 3);
      ws.set_size(nqp);
      ::sfe_qrules::prism_deg1_qrule(&cs[0], &(ws.data())[0]);
    } else if (qd <= 2) {
      ::coder::SizeType nqp;
      nqp = ::sfe_qrules::prism_deg2_qrule();
      cs.set_size(nqp, 3);
      ws.set_size(nqp);
      ::sfe_qrules::prism_deg2_qrule(&cs[0], &(ws.data())[0]);
    } else if (qd <= 3) {
      ::coder::SizeType nqp;
      nqp = ::sfe_qrules::prism_deg3_qrule();
      cs.set_size(nqp, 3);
      ws.set_size(nqp);
      ::sfe_qrules::prism_deg3_qrule(&cs[0], &(ws.data())[0]);
    } else if (qd <= 5) {
      ::coder::SizeType nqp;
      nqp = ::sfe_qrules::prism_deg5_qrule();
      cs.set_size(nqp, 3);
      ws.set_size(nqp);
      ::sfe_qrules::prism_deg5_qrule(&cs[0], &(ws.data())[0]);
    } else if (qd <= 7) {
      ::coder::SizeType nqp;
      nqp = ::sfe_qrules::prism_deg7_qrule();
      cs.set_size(nqp, 3);
      ws.set_size(nqp);
      ::sfe_qrules::prism_deg7_qrule(&cs[0], &(ws.data())[0]);
    } else if (qd <= 9) {
      ::coder::SizeType nqp;
      nqp = ::sfe_qrules::prism_deg9_qrule();
      cs.set_size(nqp, 3);
      ws.set_size(nqp);
      ::sfe_qrules::prism_deg9_qrule(&cs[0], &(ws.data())[0]);
    } else if (qd <= 11) {
      ::coder::SizeType nqp;
      nqp = ::sfe_qrules::prism_deg11_qrule();
      cs.set_size(nqp, 3);
      ws.set_size(nqp);
      ::sfe_qrules::prism_deg11_qrule(&cs[0], &(ws.data())[0]);
    } else {
      ::coder::SizeType nqp;
      if (qd > 13) {
        m2cWarnMsgIdAndTxt("prism_quadrules:UnsupportedDegree",
                           "Only support up to degree 13");
      }
      nqp = ::sfe_qrules::prism_deg13_qrule();
      cs.set_size(nqp, 3);
      ws.set_size(nqp);
      ::sfe_qrules::prism_deg13_qrule(&cs[0], &(ws.data())[0]);
    }
  } break;
  default:
    m2cAssert(shape == 5, "Unsupported element type");
    pyra_quadrules(qd, cs, ws);
    break;
  }
}

// tabulate_shapefuncs - kernel for tabulating shape functions
static void tabulate_shapefuncs(::coder::SizeType etype,
                                const ::coder::array<real_T, 2U> &cs,
                                ::coder::array<real_T, 2U> &sfvals,
                                ::coder::array<real_T, 3U> &sdvals)
{
  ::coder::SizeType shape;
  shape = etype >> 5 & 7;
  switch (((shape > 0) + (shape > 1)) + (shape > 3)) {
  case 3:
    sfe3_tabulate_shapefuncs(etype, cs, sfvals, sdvals);
    break;
  case 2:
    sfe2_tabulate_shapefuncs(etype, cs, sfvals, sdvals);
    break;
  default:
    sfe1_tabulate_shapefuncs(etype, cs, sfvals, sdvals);
    break;
  }
}

// tet_20 - Compute shape functions and their derivatives of tet_20
static inline
void tet_20(real_T xi, real_T eta, real_T zeta, real_T sfvals[20],
                   real_T sdvals[60])
{
  ::sfe_sfuncs::tet_20_sfunc(xi, eta, zeta, &sfvals[0], &sdvals[0]);
}

// tet_35 - Compute shape functions and their derivatives of tet_35
static inline
void tet_35(real_T xi, real_T eta, real_T zeta, real_T sfvals[35],
                   real_T sdvals[105])
{
  ::sfe_sfuncs::tet_35_sfunc(xi, eta, zeta, &sfvals[0], &sdvals[0]);
}

// tet_56 - Compute shape functions and their derivatives of tet_56
static inline
void tet_56(real_T xi, real_T eta, real_T zeta, real_T sfvals[56],
                   real_T sdvals[168])
{
  ::sfe_sfuncs::tet_56_sfunc(xi, eta, zeta, &sfvals[0], &sdvals[0]);
}

// tet_84 - Compute shape functions and their derivatives of tet_84
static inline
void tet_84(real_T xi, real_T eta, real_T zeta, real_T sfvals[84],
                   real_T sdvals[252])
{
  ::sfe_sfuncs::tet_84_sfunc(xi, eta, zeta, &sfvals[0], &sdvals[0]);
}

// tet_gl_20 - Compute shape functions and their derivatives of tet_gl_20
static inline
void tet_gl_20(real_T xi, real_T eta, real_T zeta, real_T sfvals[20],
                      real_T sdvals[60])
{
  ::sfe_sfuncs::tet_gl_20_sfunc(xi, eta, zeta, &sfvals[0], &sdvals[0]);
}

// tet_gl_35 - Compute shape functions and their derivatives of tet_gl_35
static inline
void tet_gl_35(real_T xi, real_T eta, real_T zeta, real_T sfvals[35],
                      real_T sdvals[105])
{
  ::sfe_sfuncs::tet_gl_35_sfunc(xi, eta, zeta, &sfvals[0], &sdvals[0]);
}

// tri_21 - Compute shape functions and their derivatives of tri_21
static inline
void tri_21(real_T xi, real_T eta, real_T sfvals[21], real_T sdvals[42])
{
  ::sfe_sfuncs::tri_21_sfunc(xi, eta, &sfvals[0], &sdvals[0]);
}

// tri_28 - Compute shape functions and their derivatives of tri_28
static inline
void tri_28(real_T xi, real_T eta, real_T sfvals[28], real_T sdvals[56])
{
  ::sfe_sfuncs::tri_28_sfunc(xi, eta, &sfvals[0], &sdvals[0]);
}

// tri_fek_15 - Compute shape functions and their derivatives of tri_fek_15
static inline
void tri_fek_15(real_T xi, real_T eta, real_T sfvals[15],
                       real_T sdvals[30])
{
  ::sfe_sfuncs::tri_fek_15_sfunc(xi, eta, &sfvals[0], &sdvals[0]);
}

// tri_fek_21 - Compute shape functions and their derivatives of tri_fek_21
static inline
void tri_fek_21(real_T xi, real_T eta, real_T sfvals[21],
                       real_T sdvals[42])
{
  ::sfe_sfuncs::tri_fek_21_sfunc(xi, eta, &sfvals[0], &sdvals[0]);
}

// tri_fek_28 - Compute shape functions and their derivatives of tri_fek_28
static inline
void tri_fek_28(real_T xi, real_T eta, real_T sfvals[28],
                       real_T sdvals[56])
{
  ::sfe_sfuncs::tri_fek_28_sfunc(xi, eta, &sfvals[0], &sdvals[0]);
}

// tri_gl_21 - Compute shape functions and their derivatives of tri_gl_21
static inline
void tri_gl_21(real_T xi, real_T eta, real_T sfvals[21],
                      real_T sdvals[42])
{
  ::sfe_sfuncs::tri_gl_21_sfunc(xi, eta, &sfvals[0], &sdvals[0]);
}

//  tri_quadrules - Obtain quadrature points and weights of a triangular
static void tri_quadrules(::coder::SizeType degree, ::coder::array<real_T, 2U> &cs,
                          ::coder::array<real_T, 1U> &ws)
{
  if (degree <= 1) {
    ::coder::SizeType nqp;
    nqp = ::sfe_qrules::tri_deg1_qrule();
    cs.set_size(nqp, 2);
    ws.set_size(nqp);
    ::sfe_qrules::tri_deg1_qrule(&cs[0], &(ws.data())[0]);
  } else if (degree <= 2) {
    ::coder::SizeType nqp;
    nqp = ::sfe_qrules::tri_deg2_qrule();
    cs.set_size(nqp, 2);
    ws.set_size(nqp);
    ::sfe_qrules::tri_deg2_qrule(&cs[0], &(ws.data())[0]);
  } else if (degree <= 4) {
    ::coder::SizeType nqp;
    nqp = ::sfe_qrules::tri_deg4_qrule();
    cs.set_size(nqp, 2);
    ws.set_size(nqp);
    ::sfe_qrules::tri_deg4_qrule(&cs[0], &(ws.data())[0]);
  } else if (degree <= 5) {
    ::coder::SizeType nqp;
    nqp = ::sfe_qrules::tri_deg5_qrule();
    cs.set_size(nqp, 2);
    ws.set_size(nqp);
    ::sfe_qrules::tri_deg5_qrule(&cs[0], &(ws.data())[0]);
  } else if (degree <= 7) {
    ::coder::SizeType nqp;
    nqp = ::sfe_qrules::tri_deg7_qrule();
    cs.set_size(nqp, 2);
    ws.set_size(nqp);
    ::sfe_qrules::tri_deg7_qrule(&cs[0], &(ws.data())[0]);
  } else if (degree <= 8) {
    ::coder::SizeType nqp;
    nqp = ::sfe_qrules::tri_deg8_qrule();
    cs.set_size(nqp, 2);
    ws.set_size(nqp);
    ::sfe_qrules::tri_deg8_qrule(&cs[0], &(ws.data())[0]);
  } else if (degree <= 9) {
    ::coder::SizeType nqp;
    nqp = ::sfe_qrules::tri_deg9_qrule();
    cs.set_size(nqp, 2);
    ws.set_size(nqp);
    ::sfe_qrules::tri_deg9_qrule(&cs[0], &(ws.data())[0]);
  } else if (degree <= 11) {
    ::coder::SizeType nqp;
    nqp = ::sfe_qrules::tri_deg11_qrule();
    cs.set_size(nqp, 2);
    ws.set_size(nqp);
    ::sfe_qrules::tri_deg11_qrule(&cs[0], &(ws.data())[0]);
  } else {
    ::coder::SizeType nqp;
    if (degree > 13) {
      m2cWarnMsgIdAndTxt("tri_quadrules:UnsupportedDegree",
                         "Only support up to degree 13");
    }
    nqp = ::sfe_qrules::tri_deg13_qrule();
    cs.set_size(nqp, 2);
    ws.set_size(nqp);
    ::sfe_qrules::tri_deg13_qrule(&cs[0], &(ws.data())[0]);
  }
}

// update_osusop - Update CRS given a subset stencil
static void
update_osusop(CrsMatrix *A, ::coder::array<int32_T, 1U> &nnzs,
              const ::coder::array<real_T, 2U> &mesh_coords,
              const ::coder::array<int64_T, 1U> &mesh_node2elems_row_ptr,
              const ::coder::array<int32_T, 1U> &mesh_node2elems_col_ind,
              const ::coder::array<Stencils, 1U> &mesh_stencils, ::coder::SizeType stclid)
{
  ::coder::array<boolean_T, 1U> visited_;
  ::coder::SizeType i;
  ::coder::SizeType loop_ub;
  visited_.set_size(A->ncols);
  loop_ub = A->ncols;
  for (i = 0; i < loop_ub; i++) {
    visited_[i] = false;
  }
  if (mesh_stencils[stclid - 1].vidmap.size(0) == 0) {
    i = mesh_coords.size(0);
  } else {
    i = mesh_stencils[stclid - 1].vidmap.size(0);
  }
  for (::coder::SizeType b_i{0}; b_i < i; b_i++) {
    ::coder::SizeType nid;
    if (mesh_stencils[stclid - 1].vidmap.size(0) == 0) {
      nid = b_i;
    } else {
      nid = mesh_stencils[stclid - 1].vidmap[b_i] - 1;
    }
    for (int64_T j = mesh_node2elems_row_ptr[nid];
         j < mesh_node2elems_row_ptr[nid + 1]; j++) {
      int64_T k;
      int64_T k_tmp;
      ::coder::SizeType i1;
      loop_ub = mesh_node2elems_col_ind[j - 1] - 1;
      k_tmp = A->row_ptr[loop_ub];
      k = k_tmp;
      k_tmp += nnzs[loop_ub];
      while (k <= k_tmp - 1L) {
        visited_[A->col_ind[k - 1] - 1] = true;
        k++;
      }
      i1 = (
          mesh_stencils[stclid - 1].ngbverts.row_ptr[b_i + 1] -
          mesh_stencils[stclid - 1].ngbverts.row_ptr[b_i]);
      for (::coder::SizeType b_k{0}; b_k < i1; b_k++) {
        ::coder::SizeType v;
        v = mesh_stencils[stclid - 1]
                .ngbverts
                .col_ind[static_cast<::coder::SizeType>(
                             (mesh_stencils[stclid - 1].ngbverts.row_ptr[b_i] +
                              (b_k + 1)) -
                             1L) -
                         1];
        if (!visited_[v - 1]) {
          if (nnzs[loop_ub] >= A->row_ptr[loop_ub + 1] - A->row_ptr[loop_ub]) {
            insert_mem_crs(loop_ub + 1, A->row_ptr, A->col_ind, A->val);
          }
          A->col_ind[(A->row_ptr[loop_ub] + nnzs[loop_ub]) -
                     1] = v;
          nnzs[loop_ub] = nnzs[loop_ub] + 1;
          visited_[v - 1] = true;
        }
      }
      k = A->row_ptr[loop_ub];
      k_tmp = (A->row_ptr[loop_ub] + nnzs[loop_ub]) - 1L;
      while (k <= k_tmp) {
        visited_[A->col_ind[k - 1] - 1] = false;
        k++;
      }
    }
  }
}

//  wls_buhmann_weights  Weights based on Buhmann's radial basis function
static void wls_buhmann_weights(const ::coder::array<real_T, 2U> &us,
                                ::coder::SizeType npoints, ::coder::SizeType degree,
                                const ::coder::array<real_T, 1U> &params_sh,
                                const ::coder::array<real_T, 2U> &params_pw,
                                ::coder::array<real_T, 1U> &ws)
{
  static const real_T dv[9]{2.6, 2.0, 1.6, 1.6, 1.6, 1.5, 1.4, 1.3, 1.2};
  real_T d;
  real_T dist_k;
  real_T r;
  real_T r1;
  real_T r2;
  real_T rho;
  real_T sigma;
  ::coder::SizeType abs_degree;
  ::coder::SizeType i;
  if (degree == 0) {
    degree = 2;
  }
  if (degree < 0) {
    abs_degree = 1 - degree;
  } else {
    abs_degree = degree + 1;
  }
  if ((params_sh.size(0) != 0) && (params_sh[0] != 0.0)) {
    sigma = params_sh[0];

    //  Assign default rho
  } else if (abs_degree - 1 >= 9) {
    sigma = 1.2;
  } else {
    sigma = dv[abs_degree - 2];
  }
  if (ws.size(0) == 0) {
    ws.set_size(npoints);
  } else {
    m2cAssert(ws.size(0) >= npoints,
              "length of ws cannot be smaller than npoints");
  }
  //  Compute rho to be sigma times the kth distance for k=ceil(1.5*ncoff)
  if (degree >= 0) {
    //  Compute 2-norm
    i = us.size(1);
    for (::coder::SizeType b_i{0}; b_i < npoints; b_i++) {
      d = us[us.size(1) * b_i];
      r2 = d * d;
      for (::coder::SizeType j{2}; j <= i; j++) {
        d = us[(j + us.size(1) * b_i) - 1];
        r2 += d * d;
      }
      ws[b_i] = std::sqrt(r2);
    }
  } else {
    //  Compute inf-norm for tensor-product
    i = us.size(1);
    for (::coder::SizeType b_i{0}; b_i < npoints; b_i++) {
      r = std::abs(us[us.size(1) * b_i]);
      for (::coder::SizeType j{2}; j <= i; j++) {
        r1 = std::abs(us[(j + us.size(1) * b_i) - 1]);
        if (r1 > r) {
          r = r1;
        }
      }
      ws[b_i] = r;
    }
  }
  if (us.size(1) == 1) {
    i = abs_degree;
  } else if (us.size(1) == 2) {
    if (degree < 0) {
      i = abs_degree * abs_degree;
    } else {
      i = (abs_degree + 1) * abs_degree / 2;
    }
  } else if (degree < 0) {
    i = abs_degree * abs_degree * abs_degree;
  } else {
    i = (abs_degree + 2) * (abs_degree + 1) * abs_degree / 6;
  }
  dist_k = find_kth_shortest_dist(ws, (i * 3 + 1) / 2, 1, npoints);
  rho = sigma * dist_k;
  if ((params_pw.size(0) == 0) || (params_pw.size(1) == 0)) {
    for (::coder::SizeType b_i{0}; b_i < npoints; b_i++) {
      if (degree > 0) {
        //  Compute 2-norm
        d = us[us.size(1) * b_i];
        r2 = d * d;
        i = us.size(1);
        for (::coder::SizeType j{2}; j <= i; j++) {
          d = us[(j + us.size(1) * b_i) - 1];
          r2 += d * d;
        }
        r = std::sqrt(r2);
      } else {
        //  Compute inf-norm for tensor-product
        r = std::abs(us[us.size(1) * b_i]);
        i = us.size(1);
        for (::coder::SizeType j{2}; j <= i; j++) {
          r1 = std::abs(us[(j + us.size(1) * b_i) - 1]);
          if (r1 > r) {
            r = r1;
          }
        }
      }
      if (r > rho) {
        ws[b_i] = 0.0;
      } else {
        real_T r_sqrt;
        r /= rho;
        r_sqrt = std::sqrt(r);
        ws[b_i] = r * r *
                      (r * r_sqrt *
                           (r_sqrt * (r_sqrt * 112.0 / 45.0 + -7.0) +
                            5.333333333333333) +
                       -0.93333333333333335) +
                  0.1111111111111111;
      }
    }
  } else {
    m2cAssert(params_pw.size(0) >= npoints,
              "size(params_pw,1) should be >=npoints");
    for (::coder::SizeType b_i{0}; b_i < npoints; b_i++) {
      real_T b_gamma;
      b_gamma = params_pw[params_pw.size(1) * b_i];
      if (b_gamma <= 0.0) {
        ws[b_i] = 0.0;
      } else {
        if (degree > 0) {
          //  Compute 2-norm
          d = us[us.size(1) * b_i];
          r2 = d * d;
          i = us.size(1);
          for (::coder::SizeType j{2}; j <= i; j++) {
            d = us[(j + us.size(1) * b_i) - 1];
            r2 += d * d;
          }
          r = std::sqrt(r2);
        } else {
          //  Compute inf-norm for tensor-product
          r = std::abs(us[us.size(1) * b_i]);
          i = us.size(1);
          for (::coder::SizeType j{2}; j <= i; j++) {
            r1 = std::abs(us[(j + us.size(1) * b_i) - 1]);
            if (r1 > r) {
              r = r1;
            }
          }
        }
        if (r > rho) {
          ws[b_i] = 0.0;
        } else {
          real_T r_sqrt;
          r /= rho;
          r_sqrt = std::sqrt(r);
          ws[b_i] =
              b_gamma * (r * r *
                             (r * r_sqrt *
                                  (r_sqrt * (r_sqrt * 112.0 / 45.0 + -7.0) +
                                   5.333333333333333) +
                              -0.93333333333333335) +
                         0.1111111111111111);
        }
      }
    }
  }
}

//  wls_eno_weights  WLS-ENO weights based on function values
static void wls_eno_weights(const ::coder::array<real_T, 2U> &us,
                            ::coder::SizeType npoints, ::coder::SizeType degree,
                            const ::coder::array<real_T, 2U> &us_unscaled,
                            const ::coder::array<real_T, 1U> &params_sh,
                            const ::coder::array<real_T, 2U> &params_pw,
                            ::coder::array<real_T, 1U> &ws)
{
  real_T c1dfg;
  real_T h2bar;
  real_T h2bar_tmp;
  real_T safegauard;
  m2cAssert(params_sh.size(0) >= 2, "first two shared parameters are required");
  m2cAssert(params_pw.size(0) >= npoints,
            "size(params_pw,1) should be >=npoints");
  m2cAssert(params_pw.size(1) >= 2, "size(params_pw,2) should be >=2");
  if (ws.size(0) == 0) {
    ws.set_size(npoints);
  } else {
    m2cAssert(ws.size(0) >= npoints,
              "length of ws cannot be smaller than npoints");
  }
  //  Compute hbar using ws as buffer space
  if (degree >= 0) {
    ::coder::SizeType i;
    //  Compute 2-norm
    i = us_unscaled.size(1);
    for (::coder::SizeType b_i{0}; b_i < npoints; b_i++) {
      real_T r2;
      h2bar_tmp = us_unscaled[us_unscaled.size(1) * b_i];
      r2 = h2bar_tmp * h2bar_tmp;
      for (::coder::SizeType j{2}; j <= i; j++) {
        h2bar_tmp = us_unscaled[(j + us_unscaled.size(1) * b_i) - 1];
        r2 += h2bar_tmp * h2bar_tmp;
      }
      ws[b_i] = std::sqrt(r2);
    }
  } else {
    ::coder::SizeType i;
    //  Compute inf-norm for tensor-product
    i = us_unscaled.size(1);
    for (::coder::SizeType b_i{0}; b_i < npoints; b_i++) {
      real_T r;
      r = std::abs(us_unscaled[us_unscaled.size(1) * b_i]);
      for (::coder::SizeType j{2}; j <= i; j++) {
        real_T r1;
        r1 = std::abs(us_unscaled[(j + us_unscaled.size(1) * b_i) - 1]);
        if (r1 > r) {
          r = r1;
        }
      }
      ws[b_i] = r;
    }
  }
  h2bar = ws[0] * ws[0];
  for (::coder::SizeType b_i{2}; b_i <= npoints; b_i++) {
    h2bar_tmp = ws[b_i - 1];
    h2bar += h2bar_tmp * h2bar_tmp;
  }
  h2bar /= static_cast<real_T>(npoints);
  //  Evaluate the inverse-distance weights as base
  wls_invdist_weights(us, npoints, 0.5 - static_cast<real_T>(degree < 0), ws);
  c1dfg = 0.05 * params_sh[1];
  safegauard = 0.001 * (params_sh[1] * params_sh[1]) * h2bar;
  if (params_pw.size(1) > 2) {
    for (::coder::SizeType b_i{0}; b_i < npoints; b_i++) {
      h2bar_tmp = params_pw[params_pw.size(1) * b_i] - params_sh[0];
      ws[b_i] = ws[b_i] / ((h2bar_tmp * h2bar_tmp +
                            c1dfg * params_pw[params_pw.size(1) * b_i + 1]) +
                           safegauard);
    }
  }
}

//  wls_func - Evaluate wls-fitting at one or more points.
static inline
void wls_func(e_WlsObject *wls,
                     const ::coder::array<real_T, 2U> &eval_pnts,
                     ::coder::array<real_T, 2U> &varargout_1)
{
  wls_kernel(wls, eval_pnts, varargout_1);
}

//  wls_init  Initialize WlsObject in 1D, 2D, or 3D.
static void wls_init(e_WlsObject *wls, const ::coder::array<real_T, 2U> &us,
                     const char_T weight_name_data[],
                     const ::coder::array<real_T, 1U> &weight_params_shared,
                     const ::coder::array<real_T, 2U> &weight_params_pointwise,
                     const ::coder::array<boolean_T, 1U> &weight_omit_rows,
                     ::coder::SizeType degree, ::coder::SizeType interp0, ::coder::SizeType nstpnts)
{
  static const real_T dv[7]{333.33333333333331,
                            1000.0,
                            3333.3333333333335,
                            10000.0,
                            100000.0,
                            1.0E+6,
                            1.0E+7};
  real_T timestamp;
  ::coder::SizeType dim;
  m2cAssert(us.size(1) >= 1, "");
  //  Process input arguments
  dim = us.size(1);
  wls->interp0 = interp0;
  wls->unimono = false;
  if (nstpnts <= 0) {
    nstpnts = us.size(0);
  } else {
    m2cAssert(
        nstpnts <= us.size(0),
        "Number of points cannot be greater than the first dimension of `us`.");
  }
  //  Resize buffers
  wls_resize(wls, us.size(1), nstpnts, degree);
  if (nstpnts != 0) {
    real_T maxx;
    real_T maxx_inv;
    real_T thres;
    ::coder::SizeType a;
    ::coder::SizeType b_i;
    ::coder::SizeType b_us;
    ::coder::SizeType ncols;
    ::coder::SizeType u1;
    if (wls->interp0 != 0) {
      //  Make the first node the origin in interp0 mode
      switch (us.size(1)) {
      case 1: {
        boolean_T b;
        boolean_T b1;
        wls->origin.size[1] = 1;
        wls->origin.size[0] = 1;
        wls->origin.data[0] = us[0];
        b = true;
        b1 = ((us.size(1) <= 0) || (us.size(0) <= 0));
        b_i = us.size(1) * us.size(0);
        b_us = 0;
        for (::coder::SizeType i{0}; i < nstpnts; i++) {
          if (b1 || (i >= b_i)) {
            b_us = 0;
            b = true;
          } else if (b) {
            b = false;
            b_us = i % us.size(0) * us.size(1) + i / us.size(0);
          } else {
            a = us.size(1) * us.size(0) - 1;
            if (b_us > MAX_int32_T - us.size(1)) {
              b_us = i % us.size(0) * us.size(1) + i / us.size(0);
            } else {
              b_us += us.size(1);
              if (b_us > a) {
                b_us -= a;
              }
            }
          }
          a = wls->us.size(0);
          wls->us[i % a * wls->us.size(1) + i / a] = us[b_us] - us[0];
        }
      } break;
      case 2:
        wls->origin.size[1] = 2;
        wls->origin.size[0] = 1;
        wls->origin.data[0] = us[0];
        wls->origin.data[1] = us[1];
        for (::coder::SizeType i{0}; i < nstpnts; i++) {
          wls->us[wls->us.size(1) * i] = us[us.size(1) * i] - us[0];
          wls->us[wls->us.size(1) * i + 1] = us[us.size(1) * i + 1] - us[1];
        }
        break;
      default:
        wls->origin.size[1] = 3;
        wls->origin.size[0] = 1;
        wls->origin.data[0] = us[0];
        wls->origin.data[1] = us[1];
        wls->origin.data[2] = us[2];
        for (::coder::SizeType i{0}; i < nstpnts; i++) {
          wls->us[wls->us.size(1) * i] = us[us.size(1) * i] - us[0];
          wls->us[wls->us.size(1) * i + 1] = us[us.size(1) * i + 1] - us[1];
          wls->us[wls->us.size(1) * i + 2] = us[us.size(1) * i + 2] - us[2];
        }
        break;
      }
    } else {
      switch (us.size(1)) {
      case 1: {
        boolean_T b;
        boolean_T b1;
        wls->origin.size[1] = 1;
        wls->origin.size[0] = 1;
        wls->origin.data[0] = 0.0;
        b = true;
        b1 = ((us.size(1) <= 0) || (us.size(0) <= 0));
        b_i = us.size(1) * us.size(0);
        b_us = 0;
        for (::coder::SizeType i{0}; i < nstpnts; i++) {
          if (b1 || (i >= b_i)) {
            b_us = 0;
            b = true;
          } else if (b) {
            b = false;
            b_us = i % us.size(0) * us.size(1) + i / us.size(0);
          } else {
            a = us.size(1) * us.size(0) - 1;
            if (b_us > MAX_int32_T - us.size(1)) {
              b_us = i % us.size(0) * us.size(1) + i / us.size(0);
            } else {
              b_us += us.size(1);
              if (b_us > a) {
                b_us -= a;
              }
            }
          }
          a = wls->us.size(0);
          wls->us[i % a * wls->us.size(1) + i / a] = us[b_us];
        }
      } break;
      case 2:
        wls->origin.size[1] = 2;
        wls->origin.size[0] = 1;
        wls->origin.data[0] = 0.0;
        wls->origin.data[1] = 0.0;
        for (::coder::SizeType i{0}; i < nstpnts; i++) {
          wls->us[wls->us.size(1) * i] = us[us.size(1) * i];
          wls->us[wls->us.size(1) * i + 1] = us[us.size(1) * i + 1];
        }
        break;
      default:
        wls->origin.size[1] = 3;
        wls->origin.size[0] = 1;
        wls->origin.data[0] = 0.0;
        wls->origin.data[1] = 0.0;
        wls->origin.data[2] = 0.0;
        for (::coder::SizeType i{0}; i < nstpnts; i++) {
          wls->us[wls->us.size(1) * i] = us[us.size(1) * i];
          wls->us[wls->us.size(1) * i + 1] = us[us.size(1) * i + 1];
          wls->us[wls->us.size(1) * i + 2] = us[us.size(1) * i + 2];
        }
        break;
      }
    }
    //  Scale us to be between -1 and 1
    maxx = 0.0;
    switch (us.size(1)) {
    case 1:
      b_i = wls->us.size(0);
      for (::coder::SizeType i{0}; i < nstpnts; i++) {
        maxx = std::fmax(
            maxx, std::abs(wls->us[i % b_i * wls->us.size(1) + i / b_i]));
      }
      break;
    case 2:
      for (::coder::SizeType i{0}; i < nstpnts; i++) {
        maxx = std::fmax(maxx,
                         std::fmax(std::abs(wls->us[wls->us.size(1) * i]),
                                   std::abs(wls->us[wls->us.size(1) * i + 1])));
      }
      break;
    default:
      for (::coder::SizeType i{0}; i < nstpnts; i++) {
        maxx = std::fmax(
            maxx,
            std::fmax(std::fmax(std::abs(wls->us[wls->us.size(1) * i]),
                                std::abs(wls->us[wls->us.size(1) * i + 1])),
                      std::abs(wls->us[wls->us.size(1) * i + 2])));
      }
      break;
    }
    if (maxx == 0.0) {
      maxx_inv = 1.0;
    } else {
      maxx_inv = 1.0 / maxx;
    }
    for (::coder::SizeType i{0}; i < dim; i++) {
      wls->hs_inv[i] = maxx_inv;
    }
    //  scale wls.us
    if (maxx_inv != 1.0) {
      switch (us.size(1)) {
      case 1:
        for (::coder::SizeType i{0}; i < nstpnts; i++) {
          b_i = wls->us.size(0);
          b_us = wls->us.size(0);
          wls->us[i % b_i * wls->us.size(1) + i / b_i] =
              wls->us[i % b_us * wls->us.size(1) + i / b_us] * maxx_inv;
        }
        break;
      case 2:
        for (::coder::SizeType i{0}; i < nstpnts; i++) {
          wls->us[wls->us.size(1) * i] =
              wls->us[wls->us.size(1) * i] * maxx_inv;
          wls->us[wls->us.size(1) * i + 1] =
              wls->us[wls->us.size(1) * i + 1] * maxx_inv;
        }
        break;
      default:
        for (::coder::SizeType i{0}; i < nstpnts; i++) {
          wls->us[wls->us.size(1) * i] =
              wls->us[wls->us.size(1) * i] * maxx_inv;
          wls->us[wls->us.size(1) * i + 1] =
              wls->us[wls->us.size(1) * i + 1] * maxx_inv;
          wls->us[wls->us.size(1) * i + 2] =
              wls->us[wls->us.size(1) * i + 2] * maxx_inv;
        }
        break;
      }
    }
    //  Compute point-wise weights
    if (weight_name_data[0] == 'U') {
      //  Unit weights
      wls->rweights.set_size(0);
    } else {
      wls->rweights.set_size(wls->V.size(1));
      if (weight_name_data[0] == 'U') {
        //  unit weights
        b_us = wls->rweights.size(0);
        wls->rweights.set_size(b_us);
        for (b_i = 0; b_i < b_us; b_i++) {
          wls->rweights[b_i] = 1.0;
        }
      } else if ((weight_name_data[0] == 'I') || (weight_name_data[0] == 'i')) {
        //  inverse distance
        wls_invdist_weights(wls->us, nstpnts, degree, weight_params_shared,
                            weight_params_pointwise, wls->rweights);
      } else if ((weight_name_data[0] == 'B') || (weight_name_data[0] == 'b')) {
        //  Buhmann weights. All points share same parameters
        wls_buhmann_weights(wls->us, nstpnts, degree, weight_params_shared,
                            weight_params_pointwise, wls->rweights);
      } else {
        if ((weight_name_data[0] != 'E') && (weight_name_data[0] != 'e')) {
          m2cErrMsgIdAndTxt(
              "wlslib:WrongWeightName",
              "Weighting scheme must be Unit, InvDist, Buhmann, or ENO.");
        }
        //  WLS-ENO
        wls_eno_weights(wls->us, nstpnts, degree, us, weight_params_shared,
                        weight_params_pointwise, wls->rweights);
      }
    }
    if (wls->runtimes.size[0] != 0) {
      timestamp = static_cast<std::chrono::duration<double>>(
                      std::chrono::system_clock::now().time_since_epoch())
                      .count();
    }
    //  Compute Vandermonde system
    gen_vander(wls->us, nstpnts, degree, wls->rweights, wls->V);
    a = wls->V.size(1);
    ncols = wls->V.size(0);
    //  Compact CVM if needed
    b_us = weight_omit_rows.size(0);
    u1 = wls->nrows;
    if (b_us <= u1) {
      u1 = b_us;
    }
    for (::coder::SizeType i{0}; i < u1; i++) {
      if (weight_omit_rows[i]) {
        b_us = wls->V.size(0);
        for (b_i = 0; b_i < b_us; b_i++) {
          wls->V[i + wls->V.size(1) * b_i] = 0.0;
        }
      }
    }
    if (wls->runtimes.size[0] != 0) {
      real_T timestamp1;
      timestamp1 = static_cast<std::chrono::duration<double>>(
                       std::chrono::system_clock::now().time_since_epoch())
                       .count();
      wls->runtimes.data[0] = timestamp1 - timestamp;
      timestamp = timestamp1;
    }
    wls->nrows = a / wls->stride * nstpnts;
    wls->ncols = ncols;
    //  Perform QR with column pivoting
    if ((degree > 1) && (degree < 7)) {
      thres = dv[degree - 1];
    } else {
      thres = 1.0E+8;
    }
    //  In interp0 mode, we trim off the first row and first column.
    rrqr_factor(wls->V, thres, interp0, interp0, wls->nrows - interp0,
                ncols - interp0, wls->QR, wls->jpvt, &wls->rank, wls->work);
    wls->fullrank = wls->rank == ncols - interp0;
    wls->rowmajor = true;
    if (wls->runtimes.size[0] != 0) {
      real_T t;
      t = static_cast<std::chrono::duration<double>>(
              std::chrono::system_clock::now().time_since_epoch())
              .count();
      wls->runtimes.data[1] = t - timestamp;
    }
  }
}

//  wls_invdist_weights  Weights based on inverse distance
static void wls_invdist_weights(const ::coder::array<real_T, 2U> &us,
                                ::coder::SizeType npoints, real_T degree,
                                ::coder::array<real_T, 1U> &ws)
{
  real_T alpha;
  alpha = std::abs(degree) / 2.0;
  m2cAssert(ws.size(0) >= npoints,
            "length of ws cannot be smaller than npoints");
  for (::coder::SizeType i{0}; i < npoints; i++) {
    real_T r;
    real_T r2;
    r = std::abs(us[us.size(1) * i]);
    if (us.size(1) > 1) {
      if (degree > 0.0) {
        ::coder::SizeType b_i;
        //  Compute 2-norm
        r2 = r * r;
        b_i = us.size(1);
        for (::coder::SizeType c_i{2}; c_i <= b_i; c_i++) {
          real_T d;
          d = us[(c_i + us.size(1) * i) - 1];
          r2 += d * d;
        }
      } else {
        ::coder::SizeType b_i;
        //  Compute inf-norm for tensor-product
        b_i = us.size(1);
        for (::coder::SizeType c_i{2}; c_i <= b_i; c_i++) {
          real_T r1;
          r1 = std::abs(us[(c_i + us.size(1) * i) - 1]);
          if (r1 > r) {
            r = r1;
          }
        }
        r2 = r * r;
      }
    } else {
      r2 = r * r;
    }
    //  Compute weight
    ws[i] = std::pow(std::sqrt(r2 + 0.01), -alpha);
  }
}

//  wls_invdist_weights  Weights based on inverse distance
static void wls_invdist_weights(const ::coder::array<real_T, 2U> &us,
                                ::coder::SizeType npoints, ::coder::SizeType degree,
                                const ::coder::array<real_T, 1U> &params_sh,
                                const ::coder::array<real_T, 2U> &params_pw,
                                ::coder::array<real_T, 1U> &ws)
{
  real_T alpha;
  ::coder::SizeType b_degree;
  if (degree < 0) {
    b_degree = -degree;
  } else {
    b_degree = degree;
  }
  alpha = static_cast<real_T>(b_degree) / 2.0;
  if ((params_sh.size(0) > 1) && (params_sh[1] != 0.0)) {
    alpha = params_sh[1];
  }
  if (ws.size(0) == 0) {
    ws.set_size(npoints);
  } else {
    m2cAssert(ws.size(0) >= npoints,
              "length of ws cannot be smaller than npoints");
  }
  if ((params_pw.size(0) == 0) || (params_pw.size(1) == 0)) {
    for (::coder::SizeType i{0}; i < npoints; i++) {
      real_T r;
      real_T r2;
      r = std::abs(us[us.size(1) * i]);
      if (us.size(1) > 1) {
        if (degree > 0) {
          //  Compute 2-norm
          r2 = r * r;
          b_degree = us.size(1);
          for (::coder::SizeType b_i{2}; b_i <= b_degree; b_i++) {
            real_T d;
            d = us[(b_i + us.size(1) * i) - 1];
            r2 += d * d;
          }
        } else {
          //  Compute inf-norm for tensor-product
          b_degree = us.size(1);
          for (::coder::SizeType b_i{2}; b_i <= b_degree; b_i++) {
            real_T r1;
            r1 = std::abs(us[(b_i + us.size(1) * i) - 1]);
            if (r1 > r) {
              r = r1;
            }
          }
          r2 = r * r;
        }
      } else {
        r2 = r * r;
      }
      //  Compute weight
      ws[i] = std::pow(std::sqrt(r2 + 0.01), -alpha);
    }
  } else {
    m2cAssert(params_pw.size(0) >= npoints,
              "size(params_pw,1) should be >=npoints");
    for (::coder::SizeType i{0}; i < npoints; i++) {
      real_T b_gamma;
      b_gamma = params_pw[params_pw.size(1) * i];
      if (b_gamma <= 0.0) {
        ws[i] = 0.0;
      } else {
        real_T r;
        real_T r2;
        r = std::abs(us[us.size(1) * i]);
        if (us.size(1) > 1) {
          if (degree > 0) {
            //  Compute 2-norm
            r2 = r * r;
            b_degree = us.size(1);
            for (::coder::SizeType b_i{2}; b_i <= b_degree; b_i++) {
              real_T d;
              d = us[(b_i + us.size(1) * i) - 1];
              r2 += d * d;
            }
          } else {
            //  Compute inf-norm for tensor-product
            b_degree = us.size(1);
            for (::coder::SizeType b_i{2}; b_i <= b_degree; b_i++) {
              real_T r1;
              r1 = std::abs(us[(b_i + us.size(1) * i) - 1]);
              if (r1 > r) {
                r = r1;
              }
            }
            r2 = r * r;
          }
        } else {
          r2 = r * r;
        }
        //  Compute weight
        ws[i] = b_gamma * std::pow(std::sqrt(r2 + 0.01), -alpha);
      }
    }
  }
}

//  wls_kernel - Kernel for evaluating an operator at one or more points using
static void wls_kernel(e_WlsObject *wls,
                       const ::coder::array<real_T, 2U> &eval_pnts,
                       ::coder::array<real_T, 2U> &vdops)
{
  real_T timestamp;
  ::coder::SizeType nDims;
  ::coder::SizeType nevpnts;
  nevpnts = eval_pnts.size(0);
  if (wls->runtimes.size[0] != 0) {
    timestamp = static_cast<std::chrono::duration<double>>(
                    std::chrono::system_clock::now().time_since_epoch())
                    .count();
  }
  //  Step 1: Tabulate monomial basis functions at evaluation points
  nDims = eval_pnts.size(1) - 1;
  wls->nevpnts = eval_pnts.size(0);
  //  scale the coordinates; use wls.us as buffer
  wls->us.set_size(((eval_pnts.size(0) + 3) / 4) << 2, eval_pnts.size(1));
  if (wls->interp0 != 0) {
    for (::coder::SizeType iPoint{0}; iPoint < nevpnts; iPoint++) {
      for (::coder::SizeType dim{0}; dim <= nDims; dim++) {
        wls->us[dim + wls->us.size(1) * iPoint] =
            (eval_pnts[dim + eval_pnts.size(1) * iPoint] -
             wls->origin.data[dim]) *
            wls->hs_inv[dim];
      }
    }
  } else {
    for (::coder::SizeType iPoint{0}; iPoint < nevpnts; iPoint++) {
      for (::coder::SizeType dim{0}; dim <= nDims; dim++) {
        wls->us[dim + wls->us.size(1) * iPoint] =
            eval_pnts[dim + eval_pnts.size(1) * iPoint] * wls->hs_inv[dim];
      }
    }
  }
  //  compute the weighted Vandermonde matrix
  gen_vander(wls->us, eval_pnts.size(0), wls->degree, wls->V);
  //  Step 2: Update the RHS of WLS from Vandermonde matrix
  wls_update_rhs(wls);
  if (wls->runtimes.size[0] != 0) {
    real_T timestamp1;
    timestamp1 = static_cast<std::chrono::duration<double>>(
                     std::chrono::system_clock::now().time_since_epoch())
                     .count();
    wls->runtimes.data[2] = timestamp1 - timestamp;
    timestamp = timestamp1;
  }
  //  Step 3: Solve the Vandermonde system to build the operator
  wls_solve_sys(wls, vdops);
  //  Rearrange vdops using wls.V as work space
  if (wls->runtimes.size[0] != 0) {
    real_T t;
    t = static_cast<std::chrono::duration<double>>(
            std::chrono::system_clock::now().time_since_epoch())
            .count();
    wls->runtimes.data[3] = t - timestamp;
  }
}

//  wls_resize  Reinitialize the buffers of WlsObject
static void wls_resize(e_WlsObject *wls, ::coder::SizeType dim, ::coder::SizeType nstpnts,
                       ::coder::SizeType degree)
{
  ::coder::SizeType b_degree;
  ::coder::SizeType ncols;
  ::coder::SizeType stride;
  ::coder::SizeType x;
  wls->degree = degree;
  wls->order = 0;
  //  make stride a multiple of four
  stride = ((nstpnts + 3) / 4) << 2;
  wls->stride = stride;
  wls->us.set_size(stride, dim);
  wls->rweights.set_size(stride);
  wls->nstpnts = nstpnts;
  //  determine number of columns and allocate V and QR
  switch (dim) {
  case 1:
    if (degree < 0) {
      b_degree = -degree;
    } else {
      b_degree = degree;
    }
    ncols = b_degree + 1;
    break;
  case 2:
    if (degree > 0) {
      ncols = (degree + 1) * (degree + 2) / 2;
    } else {
      ncols = (1 - degree) * (1 - degree);
    }
    break;
  default:
    m2cAssert(dim >= 1, "Dimension must be 1, 2, or 3.");
    if (degree > 0) {
      ncols = (degree + 1) * (degree + 2) * (degree + 3) / 6;
    } else {
      ncols = (1 - degree) * (1 - degree) * (1 - degree);
    }
    break;
  }
  wls->hs_inv.set_size(1, dim);
  wls->jpvt.set_size(ncols);
  wls->V.set_size(ncols, stride);
  b_degree = (ncols - wls->interp0) + 1;
  wls->QR.set_size(b_degree, stride);
  wls->rank = 0;
  //  work space
  x = ncols << 2;
  b_degree = ncols + 1;
  if (stride >= b_degree) {
    b_degree = stride;
  }
  if (x < 4160) {
    x = 4160;
  }
  wls->work.set_size((b_degree << 5) + x);
}

// wls_solve_sys - Solve Vandermonde system to obtain the operator
static void wls_solve_sys(e_WlsObject *wls, ::coder::array<real_T, 2U> &vdops)
{
  ::coder::SizeType nRhs;
  ::coder::SizeType nrows;
  ::coder::SizeType nrows_vdops;
  ::coder::SizeType u0;
  ::coder::SizeType u1;
  nRhs = wls->rhs.size(0);
  //  Multiply by generalized inverse of Vandermonde matrix
  rrqr_rtsolve(wls->QR, wls->ncols - wls->interp0, wls->rank, wls->rhs,
               wls->rhs.size(0));
  rrqr_qmulti(wls->QR, wls->nrows - wls->interp0, wls->ncols - wls->interp0,
              wls->rank, wls->rhs, nRhs, wls->work);
  u0 = wls->ncols;
  u1 = wls->nrows;
  if (u0 >= u1) {
    nrows_vdops = u0;
  } else {
    nrows_vdops = u1;
  }
  vdops.set_size(nrows_vdops, nRhs);
  //  Transpose the operator for row-major
  u0 = nrows_vdops - wls->interp0;
  for (::coder::SizeType i{0}; i < u0; i++) {
    for (::coder::SizeType j{0}; j < nRhs; j++) {
      vdops[j + vdops.size(1) * (i + wls->interp0)] =
          wls->rhs[i + wls->rhs.size(1) * j];
    }
  }
  nrows = wls->nrows;
  if (wls->rweights.size(0) != 0) {
    for (::coder::SizeType k{0}; k < nRhs; k++) {
      for (::coder::SizeType iRow{0}; iRow < nrows; iRow++) {
        vdops[k + vdops.size(1) * iRow] =
            vdops[k + vdops.size(1) * iRow] * wls->rweights[iRow];
      }
    }
  }
  if (wls->interp0 != 0) {
    //  In interp0 mode, set the first entry based on partition of unity
    for (::coder::SizeType j{0}; j < nRhs; j++) {
      real_T s;
      s = 0.0;
      u0 = wls->nstpnts;
      for (::coder::SizeType i{2}; i <= u0; i++) {
        s += vdops[j + vdops.size(1) * (i - 1)];
      }
      vdops[j] = 1.0 - s;
    }
  }
}

// wls_update_rhs - Update the RHS vectors
static void wls_update_rhs(e_WlsObject *wls)
{
  ::coder::SizeType nevpnts_tmp;
  ::coder::SizeType u0;
  ::coder::SizeType u1;
  nevpnts_tmp = wls->nevpnts;
  u0 = wls->ncols;
  u1 = wls->nrows;
  if (u0 >= u1) {
    u1 = u0;
  }
  wls->rhs.set_size(nevpnts_tmp, u1 - wls->interp0);
  //  Summing up rows in the differential operator
  u0 = wls->ncols - wls->interp0;
  for (::coder::SizeType iMonomial{0}; iMonomial < u0; iMonomial++) {
    ::coder::SizeType j;
    j = wls->jpvt[iMonomial] + wls->interp0;
    //  Do not sum if nargin<3
    for (::coder::SizeType iEval{0}; iEval < nevpnts_tmp; iEval++) {
      wls->rhs[iMonomial + wls->rhs.size(1) * iEval] =
          wls->V[iEval + wls->V.size(1) * (j - 1)];
    }
  }
}

// wlsmesh_compute_1ring - Compute 1-ring neighbor of nodes and cells
static void wlsmesh_compute_1ring(WlsMesh *mesh)
{
  real_T mykrings[3];
  append_wlsmesh_kring(mesh);
  //  Next, compute 1-ring
  if (mesh->stencils.size(0) <= 0) {
    m2cErrMsgIdAndTxt("wlsmesh_compute_stencils:badIdx",
                      "Invalid kring stencil index %d", mesh->stencils.size(0));
  } else {
  }
  //  Determine ring
  mykrings[1] = 0.0;
  mykrings[2] = 0.0;
  mykrings[0] = 1.0;
  //  default 1-ring
  if (mesh->topo_ndims == 1) {
    compute_stencils_1d(mesh, mesh->stencils.size(0), mykrings);
  } else if (mesh->topo_ndims == 2) {
    compute_stencils_2d(mesh, mesh->stencils.size(0), mykrings);
  } else {
    m2cAssert(false, "Not impl");
  }
  //  Copy to 1-ring structures in mesh
  mesh->node2nodes = mesh->stencils[mesh->stencils.size(0) - 1].ngbverts;
  mesh->node2elems = mesh->stencils[mesh->stencils.size(0) - 1].ngbelems;
}

// wlsmesh_compute_meshprop - Compute mesh properties
static void wlsmesh_compute_meshprop(WlsMesh *mesh, ::coder::SizeType nrmidx)
{
  ::coder::array<real_T, 1U> *elemh;
  ::coder::array<real_T, 2U> *mesh_coords;
  ::coder::array<ConnData, 1U> *mesh_elemtables;
  ::coder::array<int32_T, 1U> *mesh_node2nodes_col_ind;
  ::coder::array<int64_T, 1U> *mesh_node2nodes_row_ptr;
  ::coder::array<uint64_T, 1U> *mesh_teids;
  ::coder::array<real_T, 1U> *nodeh;
  ::coder::array<real_T, 2U> *nrms;
  ::coder::array<SfeObject, 1U> sfes_;
  ::coder::array<real_T, 2U> xs_;
  SfeObject sfe;
  real_T t_data[6];
  ::coder::SizeType b_loop_ub;
  ::coder::SizeType i;
  ::coder::SizeType loop_ub;
  ::coder::SizeType nelems;
  ::coder::SizeType outsize_idx_0;
  if (mesh->teids.size(0) == 0) {
    m2cErrMsgIdAndTxt("wlsmesh_compute_meshprop:missingSetup",
                      "must call setup first");
  }
  if ((mesh->coords.size(1) == mesh->topo_ndims + 1) &&
      ((nrmidx < 1) || (nrmidx > mesh->nrmstables.size(0)))) {
    m2cErrMsgIdAndTxt("wlsmesh_compute_meshprop:missingNrms",
                      "missing normals for surface");
  }
  //  Measure
  nelems = mesh->teids.size(0);
  sfe.etypes[0] = 0;
  sfe.nnodes[0] = 0;
  sfe.etypes[1] = 0;
  sfe.nnodes[1] = 0;
  sfe.geom_dim = 0;
  sfe.topo_dim = 0;
  sfe.facetid = 0;
  sfe.nqp = 0;
  sfe.ws.set_size(0);
  sfe.cs.set_size(::coder::SizeType(0), ::coder::SizeType(0));
  sfe.shapes_sol.set_size(::coder::SizeType(0), ::coder::SizeType(0));
  sfe.shapes_geom.set_size(::coder::SizeType(0), ::coder::SizeType(0));
  sfe.derivs_sol.set_size(::coder::SizeType(0), ::coder::SizeType(0), 0);
  sfe.derivs_geom.set_size(::coder::SizeType(0), ::coder::SizeType(0), 0);
  sfe.cs_phy.set_size(::coder::SizeType(0), ::coder::SizeType(0));
  sfe.grads_sol.set_size(::coder::SizeType(0), ::coder::SizeType(0));
  sfe.grads_geom.set_size(::coder::SizeType(0), ::coder::SizeType(0));
  sfe.jacTs.set_size(0, 3);
  sfe.wdetJ.set_size(0);
  sfe.dwork1.set_size(::coder::SizeType(0), ::coder::SizeType(0));
  sfe.dwork2.set_size(::coder::SizeType(0), ::coder::SizeType(0));
  sfe.xswork.set_size(::coder::SizeType(0), ::coder::SizeType(0));
  sfe.iwork.set_size(::coder::SizeType(0), ::coder::SizeType(0));
  outsize_idx_0 = mesh->elemtables.size(0);
  sfes_.set_size(outsize_idx_0);
  for (i = 0; i < outsize_idx_0; i++) {
    sfes_[i] = sfe;
  }
  //  This is local
  mesh->elemmeas.set_size(mesh->teids.size(0));
  for (::coder::SizeType b_i{0}; b_i < nelems; b_i++) {
    uint64_T c;
    ::coder::SizeType eid;
    c = mesh->teids[b_i] & 255UL;
    eid = (mesh->teids[b_i] >> 8);
    loop_ub = mesh->elemtables[c - 1].conn.size(1);
    outsize_idx_0 = mesh->coords.size(1);
    xs_.set_size(loop_ub, outsize_idx_0);
    for (i = 0; i < loop_ub; i++) {
      for (::coder::SizeType i2{0}; i2 < outsize_idx_0; i2++) {
        xs_[i2 + xs_.size(1) * i] =
            mesh->coords[i2 +
                         mesh->coords.size(1) *
                             (mesh->elemtables[c - 1].conn
                                  [i +
                                   mesh->elemtables[c - 1]
                                           .conn.size(1) *
                                       (eid - 1)] -
                              1)];
      }
    }
    if (sfes_[c - 1].etypes[0] <= 0) {
      sfe_init(&sfes_[c - 1],
               mesh->elemtables[c - 1].etype, xs_);
    } else {
      sfe_init(&sfes_[c - 1], xs_);
    }
    mesh->elemmeas[b_i] = coder::sum(sfes_[c - 1].wdetJ);
  }
  //  Sizes
  if (mesh->coords.size(1) == mesh->topo_ndims + 1) {
    real_T a_tmp;
    real_T hmax;
    ::coder::SizeType dim;
    ::coder::SizeType m;
    ::coder::SizeType n;
    //  Surface
    mesh_coords = &mesh->coords;
    mesh_elemtables = &mesh->elemtables;
    mesh_teids = &mesh->teids;
    mesh_node2nodes_row_ptr = &mesh->node2nodes.row_ptr;
    mesh_node2nodes_col_ind = &mesh->node2nodes.col_ind;
    nrms = &mesh->nrmstables[nrmidx - 1].normals;
    elemh = &mesh->elemh;
    nodeh = &mesh->nodeh;
    n = mesh_coords->size(0);
    m = mesh_teids->size(0);
    dim = mesh_coords->size(1) - 1;
    //  serial
    mesh->nodeh.set_size(mesh_coords->size(0));
    //  Nodal h
    for (::coder::SizeType b_i{0}; b_i < n; b_i++) {
      int64_T i1;
      (*nodeh)[b_i] = 0.0;
      i1 = (*mesh_node2nodes_row_ptr)[b_i + 1];
      if ((*mesh_node2nodes_row_ptr)[b_i] <= i1 - 1L) {
        b_loop_ub = mesh_coords->size(1);
      }
      for (int64_T k = (*mesh_node2nodes_row_ptr)[b_i]; k < i1; k++) {
        real_T v_data[3];
        real_T len;
        for (i = 0; i < b_loop_ub; i++) {
          v_data[i] =
              (*mesh_coords)[i + mesh_coords->size(1) *
                                     ((*mesh_node2nodes_col_ind)
                                          [k - 1] -
                                      1)] -
              (*mesh_coords)[i + mesh_coords->size(1) * b_i];
        }
        if ((nrms->size(0) != 0) && (nrms->size(1) != 0)) {
          real_T u_data[2];
          //  Surface
          outsize_idx_0 = nrms->size(1) - 1;
          if (nrms->size(1) == 2) {
            //  2D
            t_data[0] = -(*nrms)[nrms->size(1) * b_i + 1];
            t_data[outsize_idx_0] = (*nrms)[nrms->size(1) * b_i];
          } else {
            real_T a;
            a_tmp = std::abs((*nrms)[nrms->size(1) * b_i]);
            if ((a_tmp > std::abs((*nrms)[nrms->size(1) * b_i + 1])) &&
                (a_tmp > std::abs((*nrms)[nrms->size(1) * b_i + 2]))) {
              t_data[0] = -(*nrms)[nrms->size(1) * b_i] *
                          (*nrms)[nrms->size(1) * b_i + 1];
              t_data[outsize_idx_0] =
                  1.0 - (*nrms)[nrms->size(1) * b_i + 1] *
                            (*nrms)[nrms->size(1) * b_i + 1];
              t_data[outsize_idx_0 * 2] = -(*nrms)[nrms->size(1) * b_i + 1] *
                                          (*nrms)[nrms->size(1) * b_i + 2];
            } else {
              t_data[0] = 1.0 - (*nrms)[nrms->size(1) * b_i] *
                                    (*nrms)[nrms->size(1) * b_i];
              t_data[outsize_idx_0] = -(*nrms)[nrms->size(1) * b_i] *
                                      (*nrms)[nrms->size(1) * b_i + 1];
              t_data[outsize_idx_0 * 2] = -(*nrms)[nrms->size(1) * b_i] *
                                          (*nrms)[nrms->size(1) * b_i + 2];
            }
            a_tmp = t_data[outsize_idx_0 * 2];
            a = std::sqrt((t_data[0] * t_data[0] +
                           t_data[outsize_idx_0] * t_data[outsize_idx_0]) +
                          a_tmp * a_tmp);
            t_data[0] /= a;
            t_data[outsize_idx_0] /= a;
            t_data[outsize_idx_0 * 2] /= a;
            //  cross
            t_data[1] =
                t_data[outsize_idx_0 * 2] * (*nrms)[nrms->size(1) * b_i + 1] -
                t_data[outsize_idx_0] * (*nrms)[nrms->size(1) * b_i + 2];
            t_data[outsize_idx_0 + 1] =
                t_data[0] * (*nrms)[nrms->size(1) * b_i + 2] -
                (*nrms)[nrms->size(1) * b_i] * t_data[outsize_idx_0 * 2];
            t_data[outsize_idx_0 * 2 + 1] =
                (*nrms)[nrms->size(1) * b_i] * t_data[outsize_idx_0] -
                t_data[0] * (*nrms)[nrms->size(1) * b_i + 1];
          }
          //  Project
          for (::coder::SizeType ii{0}; ii < dim; ii++) {
            u_data[ii] = 0.0;
            for (::coder::SizeType jj{0}; jj <= dim; jj++) {
              u_data[ii] += v_data[jj] * t_data[ii + outsize_idx_0 * jj];
            }
          }
          if (dim + 1 == 2) {
            len = std::abs(u_data[0]);
          } else {
            len = std::sqrt(u_data[0] * u_data[0] + u_data[1] * u_data[1]);
          }
        } else {
          len = 0.0;
          for (::coder::SizeType ii{0}; ii <= dim; ii++) {
            a_tmp = v_data[ii];
            len += std::sqrt(a_tmp * a_tmp);
          }
        }
        (*nodeh)[b_i] = (*nodeh)[b_i] + len;
      }
      //  average
      (*nodeh)[b_i] = (*nodeh)[b_i] /
                      static_cast<real_T>(i1 - (*mesh_node2nodes_row_ptr)[b_i]);
    }
    elemh->set_size(mesh_teids->size(0));
    //  serial
    hmax = 0.0;
    for (::coder::SizeType b_i{0}; b_i < m; b_i++) {
      real_T h0;
      ::coder::SizeType etable;
      etable = static_cast<::coder::SizeType>((*mesh_teids)[b_i] & 255UL) - 1;
      h0 = 0.0;
      i = (*mesh_elemtables)[etable].conn.size(1) - 1;
      for (::coder::SizeType j{0}; j <= i; j++) {
        h0 += (*nodeh)
            [(*mesh_elemtables)[etable].conn
                 [j + (*mesh_elemtables)[etable].conn.size(1) *
                          (static_cast<::coder::SizeType>((*mesh_teids)[b_i] >> 8) - 1)] -
             1];
      }
      //  Average
      a_tmp = h0 / static_cast<real_T>((*mesh_elemtables)[etable].conn.size(1));
      (*elemh)[b_i] = a_tmp;
      hmax = std::fmax(hmax, a_tmp);
    }
    //  Compute global h
    mesh->globalh = hmax;
  } else {
    real_T hmax;
    ::coder::SizeType dim;
    ::coder::SizeType m;
    ::coder::SizeType n;
    n = mesh->coords.size(0);
    m = mesh->teids.size(0);
    dim = mesh->coords.size(1);
    //  serial
    mesh->nodeh.set_size(mesh->coords.size(0));
    //  Nodal h
    for (::coder::SizeType b_i{0}; b_i < n; b_i++) {
      int64_T i1;
      int64_T k;
      mesh->nodeh[b_i] = 0.0;
      k = mesh->node2nodes.row_ptr[b_i];
      ::coder::SizeType exitg1;
      do {
        exitg1 = 0;
        i1 = mesh->node2nodes.row_ptr[b_i + 1];
        if (k <= i1 - 1L) {
          real_T v_data[3];
          real_T len;
          loop_ub = mesh->coords.size(1);
          outsize_idx_0 = mesh->node2nodes.col_ind[k - 1];
          for (i = 0; i < loop_ub; i++) {
            v_data[i] =
                mesh->coords[i + mesh->coords.size(1) * (outsize_idx_0 - 1)] -
                mesh->coords[i + mesh->coords.size(1) * b_i];
          }
          len = 0.0;
          for (::coder::SizeType ii{0}; ii < dim; ii++) {
            real_T a_tmp;
            a_tmp = v_data[ii];
            len += std::sqrt(a_tmp * a_tmp);
          }
          mesh->nodeh[b_i] = mesh->nodeh[b_i] + len;
          k++;
        } else {
          exitg1 = 1;
        }
      } while (exitg1 == 0);
      //  average
      mesh->nodeh[b_i] =
          mesh->nodeh[b_i] /
          static_cast<real_T>(i1 - mesh->node2nodes.row_ptr[b_i]);
    }
    mesh->elemh.set_size(mesh->teids.size(0));
    //  serial
    hmax = 0.0;
    for (::coder::SizeType b_i{0}; b_i < m; b_i++) {
      real_T h0;
      ::coder::SizeType etable;
      etable = (mesh->teids[b_i] & 255UL) - 1;
      h0 = 0.0;
      i = mesh->elemtables[etable].conn.size(1) - 1;
      for (::coder::SizeType j{0}; j <= i; j++) {
        h0 += mesh->nodeh[mesh->elemtables[etable]
                              .conn[j + mesh->elemtables[etable].conn.size(1) *
                                            ((
                                                 mesh->teids[b_i] >> 8) -
                                             1)] -
                          1];
      }
      //  Average
      mesh->elemh[b_i] =
          h0 / static_cast<real_T>(mesh->elemtables[etable].conn.size(1));
      hmax = std::fmax(hmax, mesh->elemh[b_i]);
    }
    //  Compute global h
    mesh->globalh = hmax;
  }
}

// wlsmesh_compute_nodeparts - Compute node-based recursive partitioning
static void wlsmesh_compute_nodeparts(WlsMesh *mesh, ::coder::SizeType nparts)
{
  ::coder::array<int32_T, 2U> buf_;
  if (mesh->teids.size(0) == 0) {
    m2cErrMsgIdAndTxt("wlsmesh_compute_nodeparts:missingSetup",
                      "call setup first");
  }
  if (mesh->elemtables.size(0) == 1) {
    omp4mRecurPartMesh(mesh->coords.size(0), mesh->elemtables[0].conn,
                       mesh->topo_ndims, mesh->topo_ndims, nparts,
                       mesh->nodeparts);
  } else {
    ::coder::SizeType etype;
    ::coder::SizeType i;
    ::coder::SizeType npe;
    //  Requires buffer
    npe = 0;
    i = mesh->elemtables.size(0) - 1;
    for (::coder::SizeType etable{0}; etable <= i; etable++) {
      etype = mesh->elemtables[etable].etype - 1;
      if (iv[etype] >= npe) {
        npe = iv[etype];
      }
    }
    buf_.set_size(mesh->teids.size(0), npe);
    for (::coder::SizeType etable{0}; etable <= i; etable++) {
      ::coder::SizeType npe0;
      npe0 = iv[mesh->elemtables[etable].etype - 1] - 1;
      etype = mesh->elemtables[etable].conn.size(0);
      for (::coder::SizeType e{0}; e < etype; e++) {
        ::coder::SizeType geid;
        ::coder::SizeType i1;
        geid = (e + mesh->elemtables[etable].istart) - 1;
        for (::coder::SizeType j{0}; j <= npe0; j++) {
          buf_[j + buf_.size(1) * geid] =
              mesh->elemtables[etable]
                  .conn[j + mesh->elemtables[etable].conn.size(1) * e];
        }
        i1 = npe0 + 2;
        for (::coder::SizeType j{i1}; j <= npe; j++) {
          buf_[(j + buf_.size(1) * geid) - 1] = 0;
        }
      }
    }
    omp4mRecurPartMesh(mesh->coords.size(0), buf_, mesh->topo_ndims,
                       mesh->topo_ndims, nparts, mesh->nodeparts);
  }
}

// wlsmesh_compute_stencils - Compute stencils for WlsMesh
static void wlsmesh_compute_stencils(WlsMesh *mesh, ::coder::SizeType stclidx)
{
  real_T mykrings[3];
  if (stclidx <= 0) {
    m2cErrMsgIdAndTxt("wlsmesh_compute_stencils:badIdx",
                      "Invalid kring stencil index %d", (int)stclidx);
  } else {
    if (stclidx > mesh->stencils.size(0)) {
      m2cErrMsgIdAndTxt("wlsmesh_compute_stencils:badIdx",
                        "Invalid kring stencil index %d", (int)stclidx);
    }
  }
  //  Determine ring
  mykrings[1] = 0.0;
  mykrings[2] = 0.0;
  mykrings[0] = 1.0;
  //  default 1-ring
  if (mesh->topo_ndims == 1) {
    b_compute_stencils_1d(mesh, stclidx, mykrings);
  } else if (mesh->topo_ndims == 2) {
    b_compute_stencils_2d(mesh, stclidx, mykrings);
  } else {
    m2cAssert(false, "Not impl");
  }
}

// wlsmesh_compute_stencils - Compute stencils for WlsMesh
static void wlsmesh_compute_stencils(WlsMesh *mesh, ::coder::SizeType stclidx,
                                     real_T krings)
{
  real_T mykrings[3];
  if (stclidx <= 0) {
    m2cErrMsgIdAndTxt("wlsmesh_compute_stencils:badIdx",
                      "Invalid kring stencil index %d", (int)stclidx);
  } else {
    if (stclidx > mesh->stencils.size(0)) {
      m2cErrMsgIdAndTxt("wlsmesh_compute_stencils:badIdx",
                        "Invalid kring stencil index %d", (int)stclidx);
    }
  }
  //  Determine ring
  mykrings[1] = 0.0;
  mykrings[2] = 0.0;
  mykrings[0] = krings;
  if (krings < 1.0) {
    mykrings[0] = 1.0;
  }
  //  default 1-ring
  if (mesh->topo_ndims == 1) {
    compute_stencils_1d(mesh, stclidx, mykrings);
  } else if (mesh->topo_ndims == 2) {
    compute_stencils_2d(mesh, stclidx, mykrings);
  } else {
    m2cAssert(false, "Not impl");
  }
}

// rdi_apply - Apply RDI to obtain indicators and disc. tags
void rdi_apply(const RdiObject *rdi, const WlsMesh *mesh,
                const ::coder::array<real_T, 2U> &fs,
                const ::coder::array<real_T, 2U> &df,
                ::coder::array<real_T, 2U> &alpha,
                ::coder::array<real_T, 2U> &beta,
                ::coder::array<int8_T, 2U> &distags,
                ::coder::array<real_T, 2U> &alphanode)
{
  int64_T j;
  real_T asum_tmp;
  real_T cglobal;
  real_T clocal;
  real_T epsbeta;
  real_T epsh2;
  real_T kappa0;
  real_T kappa1;
  real_T thres;
  ::coder::SizeType atop_tmp;
  ::coder::SizeType n;
  ::coder::SizeType nrhs;
  ::coder::SizeType nrows;
  if (!rdi->fullrank) {
    m2cErrMsgIdAndTxt("rdi_compute_osusind:badA",
                      "OSUS operator must be full-rank");
  }
  if (mesh->coords.size(0) != fs.size(0)) {
    m2cErrMsgIdAndTxt("rdi_compute_osusind:badShape",
                      "Unmatched input fs size");
  }
  //  Using advanced OpenMP mode to compute cell-based alpha value
  nrows = rdi->A.row_ptr.size(0);
  nrhs = fs.size(1) - 1;
  alpha.set_size(rdi->A.row_ptr.size(0) - 1, fs.size(1));
  //  Optimized for row-major
  for (::coder::SizeType i{0}; i <= nrows - 2; i++) {
    for (::coder::SizeType k{0}; k <= nrhs; k++) {
      alpha[k + alpha.size(1) * i] = 0.0;
    }
    for (j = rdi->A.row_ptr[i]; j < rdi->A.row_ptr[i + 1]; j++) {
      for (::coder::SizeType k{0}; k <= nrhs; k++) {
        alpha[k + alpha.size(1) * i] =
            alpha[k + alpha.size(1) * i] +
            rdi->A.val[j - 1] *
                fs[k + fs.size(1) *
                           (rdi->A.col_ind[j - 1] - 1)];
      }
    }
  }
  //  Optionally, compute node based alpha values
  n = mesh->coords.size(0);
  nrhs = alpha.size(1);
  alphanode.set_size(mesh->coords.size(0), alpha.size(1));
  //  implicit barrier
  for (::coder::SizeType i{0}; i < n; i++) {
    for (::coder::SizeType k{0}; k < nrhs; k++) {
      real_T a;
      a = 0.0;
      for (j = mesh->node2elems.row_ptr[i]; j < mesh->node2elems.row_ptr[i + 1];
           j++) {
        asum_tmp = std::abs(
            alpha[k +
                  alpha.size(1) *
                      (mesh->node2elems.col_ind[j - 1] -
                       1)]);
        if (a < asum_tmp) {
          a = asum_tmp;
        }
      }
      alphanode[k + alphanode.size(1) * i] = a;
    }
  }
  //  Step 2, compute osc. indicators (beta values)
  if (alpha.size(1) != df.size(1)) {
    m2cErrMsgIdAndTxt("rdi_compute_oscind:badShape",
                      "Unmatched # of functions");
  }
  epsbeta = rdi->epsbeta;
  if (rdi->epsbeta <= 0.0) {
    epsbeta = 0.0002;
  }
  n = mesh->coords.size(0);
  nrhs = df.size(1);
  beta.set_size(mesh->coords.size(0), df.size(1));
  //  Serial
  epsh2 = epsbeta * mesh->globalh * mesh->globalh;
  //  C++
  for (::coder::SizeType i{0}; i < n; i++) {
    for (::coder::SizeType k{0}; k < nrhs; k++) {
      int64_T b_i;
      real_T abar;
      real_T asum;
      real_T atop;
      real_T wbar;
      asum = 0.0;
      wbar = 0.0;
      j = mesh->node2elems.row_ptr[i];
      ::coder::SizeType exitg1;
      do {
        exitg1 = 0;
        b_i = mesh->node2elems.row_ptr[i + 1] - 1L;
        if (j <= b_i) {
          atop_tmp = mesh->node2elems.col_ind[j - 1] - 1;
          asum_tmp = mesh->elemmeas[atop_tmp];
          asum += asum_tmp * alpha[k + alpha.size(1) * atop_tmp];
          wbar += asum_tmp;
          j++;
        } else {
          exitg1 = 1;
        }
      } while (exitg1 == 0);
      //  sum(w_e*a_e)/sum(w_e)
      abar = asum / wbar;
      atop = 0.0;
      for (j = mesh->node2elems.row_ptr[i]; j <= b_i; j++) {
        atop_tmp = mesh->node2elems.col_ind[j - 1] - 1;
        atop += mesh->elemmeas[atop_tmp] *
                std::abs(alpha[k + alpha.size(1) * atop_tmp] - abar);
      }
      beta[k + beta.size(1) * i] =
          atop /
          ((std::abs(asum) + wbar * (df[k] * epsh2)) + 2.2250738585072014E-308);
    }
  }
  //  Step 3, determine disc. nodes
  if (mesh->coords.size(0) != fs.size(0)) {
    m2cErrMsgIdAndTxt("rdi_mark_discontinuities:badShape",
                      "Unmatched input fs size");
  }
  if ((df.size(1) != fs.size(1)) || (alpha.size(1) != df.size(1)) ||
      (beta.size(1) != df.size(1))) {
    m2cErrMsgIdAndTxt("rdi_mark_discontinuities:badShape",
                      "Unmatched # of functions");
  }
  //  Get parameters
  cglobal = rdi->cglobal;
  if (rdi->cglobal < 0.0) {
    cglobal = 0.05;
  }
  clocal = rdi->clocal;
  if (rdi->clocal < 0.0) {
    clocal = 0.5;
  }
  kappa1 = rdi->kappa1;
  if (rdi->kappa1 <= 0.0) {
    kappa1 = 0.3;
  }
  //  C1 dis
  kappa0 = rdi->kappa0;
  if (rdi->kappa0 <= 0.0) {
    kappa0 = 1.0;
  }
  //  C0 dis
  if (kappa0 <= kappa1) {
    kappa0 = kappa1;
  }
  n = mesh->coords.size(0);
  nrhs = alpha.size(1);
  distags.set_size(mesh->coords.size(0), alpha.size(1));
  thres = cglobal * std::pow(mesh->globalh, 1.5);
  for (::coder::SizeType i{0}; i < n; i++) {
    for (::coder::SizeType k{0}; k < nrhs; k++) {
      real_T tauglobal;
      boolean_T discell;
      boolean_T exitg2;
      distags[k + distags.size(1) * i] = 0;
      tauglobal = thres * df[k];
      discell = false;
      j = mesh->node2elems.row_ptr[i];
      exitg2 = false;
      while ((!exitg2) && (j <= mesh->node2elems.row_ptr[i + 1] - 1L)) {
        real_T b_fmax;
        real_T b_fmin;
        uint64_T c;
        ::coder::SizeType i1;
        ::coder::SizeType leid;
        atop_tmp = mesh->node2elems.col_ind[j - 1] - 1;
        c = mesh->teids[atop_tmp] & 255UL;
        leid = (
            mesh->teids[mesh->node2elems.col_ind[j - 1] -
                        1] >>
            8);
        b_fmax = -1.7976931348623157E+308;
        b_fmin = 1.7976931348623157E+308;
        i1 =
            mesh->elemtables[(
                                 mesh->teids[mesh->node2elems.col_ind
                                                 [j - 1] -
                                             1] &
                                 255UL) -
                             1]
                .conn.size(1);
        for (::coder::SizeType ii{0}; ii < i1; ii++) {
          real_T fvalue;
          fvalue =
              fs[k +
                 fs.size(1) *
                     (mesh->elemtables[c - 1].conn
                          [ii + mesh->elemtables[c - 1]
                                        .conn.size(1) *
                                    (leid - 1)] -
                      1)];
          if (b_fmax < fvalue) {
            b_fmax = fvalue;
          }
          if (b_fmin > fvalue) {
            b_fmin = fvalue;
          }
        }
        //  Compute local df and thres
        if (std::abs(alpha[k + alpha.size(1) * atop_tmp]) >
            std::fmax(tauglobal, clocal * (b_fmax - b_fmin) *
                                     std::sqrt(mesh->elemh[atop_tmp]))) {
          //  Dis cell
          discell = true;
          exitg2 = true;
        } else {
          j++;
        }
      }
      if (discell && (beta[k + beta.size(1) * i] > kappa1)) {
        distags[k + distags.size(1) * i] = 2;
        //  C1
        if (beta[k + beta.size(1) * i] > kappa0) {
          distags[k + distags.size(1) * i] = 1;
        }
        //  C0
      }
    }
  }
}

// rdi_apply - Apply RDI to obtain indicators and disc. tags
void rdi_apply(const RdiObject *rdi, const WlsMesh *mesh,
                const ::coder::array<real_T, 2U> &fs,
                const ::coder::array<real_T, 2U> &df, ::coder::SizeType varargin_2,
                ::coder::array<real_T, 2U> &alpha,
                ::coder::array<real_T, 2U> &beta,
                ::coder::array<int8_T, 2U> &distags,
                ::coder::array<real_T, 2U> &alphanode)
{
  boolean_T m2cTryBlkErrFlag;
  ::coder::SizeType nthreads;
  if (varargin_2 <= 0) {
    nthreads = 1;
#ifdef _OPENMP
    nthreads = omp_get_max_threads();
#endif // _OPENMP
  } else {
    nthreads = varargin_2;
  }
  m2cTryBlkErrFlag = 0;
#pragma omp parallel default(shared) num_threads(nthreads)
  try { // try
    //  Step 1, compute alpha values
    rdi_compute_osusind(rdi->A.row_ptr, rdi->A.col_ind, rdi->A.val,
                        rdi->fullrank, mesh->coords, mesh->node2elems.row_ptr,
                        mesh->node2elems.col_ind, fs, alpha, alphanode);
#pragma omp barrier
    //  Step 2, compute osc. indicators (beta values)
    rdi_compute_oscind(rdi->epsbeta, mesh->coords, mesh->node2elems.row_ptr,
                       mesh->node2elems.col_ind, mesh->elemmeas, mesh->globalh,
                       df, alpha, beta);
#pragma omp barrier
    //  Step 3, determine disc. nodes
    rdi_mark_discontinuities(rdi->cglobal, rdi->clocal, rdi->kappa0,
                             rdi->kappa1, mesh->coords, mesh->elemtables,
                             mesh->teids, mesh->node2elems.row_ptr,
                             mesh->node2elems.col_ind, mesh->elemh,
                             mesh->globalh, fs, df, alpha, beta, distags);
  } catch (const std::runtime_error &m2cExc) {
    m2cTryBlkErrFlag = 1;
    m2cPrintError("runtime_error %s\n", m2cExc.what());
  } catch (const std::logic_error &m2cExc) {
    m2cTryBlkErrFlag = 1;
    m2cPrintError("logic_error %s\n", m2cExc.what());
  } catch (const std::exception &m2cExc) {
    m2cTryBlkErrFlag = 1;
    m2cPrintError("exception %s\n", m2cExc.what());
  } catch (...) {
    m2cTryBlkErrFlag = 1;
    m2cPrintError("Unknown error detected from C++ exceptions\n");
    fflush(stderr);
  } // end try
  if ((int32_T)m2cTryBlkErrFlag != 0) {
    throw std::runtime_error("omp4m:runtimeErrorInThread");
  }
}

// rdi_apply - Apply RDI to obtain indicators and disc. tags
static inline
void rdi_apply3(const RdiObject *rdi, const WlsMesh *mesh,
                const ::coder::array<real_T, 2U> &fs,
                const ::coder::array<real_T, 2U> &df,
                ::coder::array<real_T, 2U> &alpha,
                ::coder::array<real_T, 2U> &beta,
                ::coder::array<int8_T, 2U> &distags,
                ::coder::array<real_T, 2U> &alphanode)
{
  rdi_compute_osusind(rdi->A.row_ptr, rdi->A.col_ind, rdi->A.val, rdi->fullrank,
                      mesh->coords, mesh->node2elems.row_ptr,
                      mesh->node2elems.col_ind, fs, alpha, alphanode);
#pragma omp barrier
  //  Step 2, compute osc. indicators (beta values)
  rdi_compute_oscind(rdi->epsbeta, mesh->coords, mesh->node2elems.row_ptr,
                     mesh->node2elems.col_ind, mesh->elemmeas, mesh->globalh,
                     df, alpha, beta);
#pragma omp barrier
  //  Step 3, determine disc. nodes
  rdi_mark_discontinuities(rdi->cglobal, rdi->clocal, rdi->kappa0, rdi->kappa1,
                           mesh->coords, mesh->elemtables, mesh->teids,
                           mesh->node2elems.row_ptr, mesh->node2elems.col_ind,
                           mesh->elemh, mesh->globalh, fs, df, alpha, beta,
                           distags);
}

// rdi_compute_osusop - Compute the over- and under-shoot operator
void rdi_compute_osusop(RdiObject *rdi, WlsMesh *mesh)
{
  ::coder::SizeType extstclid;
  init_osusop(mesh->coords, mesh->elemtables, mesh->teids, mesh->stencils,
              rdi->stclid, rdi->A.row_ptr, rdi->A.col_ind, rdi->A.val,
              &rdi->A.ncols, rdi->nnzs);
  if (mesh->topo_ndims == mesh->coords.size(1)) {
    //  {1,2,3}-D body assembly
    assemble_body(rdi, mesh->coords, mesh->elemtables, mesh->teids,
                  mesh->node2elems.row_ptr, mesh->node2elems.col_ind,
                  mesh->stencils, rdi->stclid, true);
  } else {
    if (mesh->topo_ndims + 1 != mesh->coords.size(1)) {
      m2cErrMsgIdAndTxt("rdi_assemble_osusop:badDim",
                        "Must be either body or surface");
    }
    assemble_surf(rdi, mesh->coords, mesh->elemtables, mesh->teids,
                  mesh->nrmstables, mesh->node2elems.row_ptr,
                  mesh->node2elems.col_ind, mesh->stencils, rdi->stclid, true);
  }
  //  Determine potential RD nodes and reset rdtags. RD nodes are stored in
  extstclid = rdi->extstclid;
  //  Overcome a "bug" in coder regarding creating buffers
  determine_rdnodes(rdi->fullrank, rdi->rdtags,
                    mesh->stencils[extstclid - 1].vidmap);
  if (rdi->fullrank) {
    crsCompress(&rdi->A, rdi->nnzs);
  }
  //  Step 2, resolve potential rank-deficient entries
  while (!rdi->fullrank) {
    rdi_update_osusop(rdi, mesh, true);
  }
}

// rdi_compute_osusop - Compute the over- and under-shoot operator
void rdi_compute_osusop2(RdiObject *rdi, WlsMesh *mesh, boolean_T interp0)
{
  ::coder::SizeType extstclid;
  init_osusop(mesh->coords, mesh->elemtables, mesh->teids, mesh->stencils,
              rdi->stclid, rdi->A.row_ptr, rdi->A.col_ind, rdi->A.val,
              &rdi->A.ncols, rdi->nnzs);
  if (mesh->topo_ndims == mesh->coords.size(1)) {
    //  {1,2,3}-D body assembly
    assemble_body(rdi, mesh->coords, mesh->elemtables, mesh->teids,
                  mesh->node2elems.row_ptr, mesh->node2elems.col_ind,
                  mesh->stencils, rdi->stclid, interp0);
  } else {
    if (mesh->topo_ndims + 1 != mesh->coords.size(1)) {
      m2cErrMsgIdAndTxt("rdi_assemble_osusop:badDim",
                        "Must be either body or surface");
    }
    assemble_surf(rdi, mesh->coords, mesh->elemtables, mesh->teids,
                  mesh->nrmstables, mesh->node2elems.row_ptr,
                  mesh->node2elems.col_ind, mesh->stencils, rdi->stclid,
                  interp0);
  }
  //  Determine potential RD nodes and reset rdtags. RD nodes are stored in
  extstclid = rdi->extstclid;
  //  Overcome a "bug" in coder regarding creating buffers
  determine_rdnodes(rdi->fullrank, rdi->rdtags,
                    mesh->stencils[extstclid - 1].vidmap);
  if (rdi->fullrank) {
    crsCompress(&rdi->A, rdi->nnzs);
  }
  //  Step 2, resolve potential rank-deficient entries
  while (!rdi->fullrank) {
    rdi_update_osusop(rdi, mesh, interp0);
  }
}

// rdi_compute_osusop - Compute the over- and under-shoot operator
void rdi_compute_osusop(RdiObject *rdi, WlsMesh *mesh, boolean_T interp0,
                         ::coder::SizeType varargin_2)
{
  ::coder::SizeType extstclid;
  init_osusop(mesh->coords, mesh->elemtables, mesh->teids, mesh->stencils,
              rdi->stclid, rdi->A.row_ptr, rdi->A.col_ind, rdi->A.val,
              &rdi->A.ncols, rdi->nnzs);
  if (mesh->topo_ndims == mesh->coords.size(1)) {
    //  {1,2,3}-D body assembly
    b_assemble_body(rdi, mesh->coords, mesh->elemtables, mesh->teids,
                    mesh->node2elems.row_ptr, mesh->node2elems.col_ind,
                    mesh->stencils, mesh->nodeparts, rdi->stclid, interp0,
                    varargin_2);
  } else {
    if (mesh->topo_ndims + 1 != mesh->coords.size(1)) {
      m2cErrMsgIdAndTxt("rdi_assemble_osusop:badDim",
                        "Must be either body or surface");
    }
    b_assemble_surf(rdi, mesh->coords, mesh->elemtables, mesh->teids,
                    mesh->nrmstables, mesh->node2elems.row_ptr,
                    mesh->node2elems.col_ind, mesh->stencils, mesh->nodeparts,
                    rdi->stclid, interp0, varargin_2);
  }
  //  Determine potential RD nodes and reset rdtags. RD nodes are stored in
  extstclid = rdi->extstclid;
  //  Overcome a "bug" in coder regarding creating buffers
  determine_rdnodes(rdi->fullrank, rdi->rdtags,
                    mesh->stencils[extstclid - 1].vidmap);
  if (rdi->fullrank) {
    crsCompress(&rdi->A, rdi->nnzs);
  }
  //  Step 2, resolve potential rank-deficient entries
  while (!rdi->fullrank) {
    rdi_update_osusop(rdi, mesh, interp0);
  }
}

void rdi_compute_osusop4(RdiObject *rdi, WlsMesh *mesh, boolean_T interp0)
{
  ::coder::SizeType extstclid;
#pragma omp single
  { // single
    init_osusop(mesh->coords, mesh->elemtables, mesh->teids, mesh->stencils,
                rdi->stclid, rdi->A.row_ptr, rdi->A.col_ind, rdi->A.val,
                &rdi->A.ncols, rdi->nnzs);
  } // single
  if (mesh->topo_ndims == mesh->coords.size(1)) {
    //  {1,2,3}-D body assembly
    c_assemble_body(rdi, mesh->coords, mesh->elemtables, mesh->teids,
                    mesh->node2elems.row_ptr, mesh->node2elems.col_ind,
                    mesh->stencils, mesh->nodeparts, rdi->stclid, interp0);
  } else {
    if (mesh->topo_ndims + 1 != mesh->coords.size(1)) {
      m2cErrMsgIdAndTxt("rdi_assemble_osusop:badDim",
                        "Must be either body or surface");
    }
    c_assemble_surf(rdi, mesh->coords, mesh->elemtables, mesh->teids,
                    mesh->nrmstables, mesh->node2elems.row_ptr,
                    mesh->node2elems.col_ind, mesh->stencils, mesh->nodeparts,
                    rdi->stclid, interp0);
  }
  //  Determine potential RD nodes and reset rdtags. RD nodes are stored in
  extstclid = rdi->extstclid;
  //  Overcome a "bug" in coder regarding creating buffers
#pragma omp barrier
#pragma omp single
  { // single
    determine_rdnodes(rdi->fullrank, rdi->rdtags,
                      mesh->stencils[extstclid - 1].vidmap);
    if (rdi->fullrank) {
      crsCompress(&rdi->A, rdi->nnzs);
    }
  } // single
#pragma omp barrier
#pragma omp single
  { // single
    while (!rdi->fullrank) {
      rdi_update_osusop(rdi, mesh, interp0);
    }
  } // single
}

// rdi_create - Create RDI object
void rdi_create(WlsMesh *mesh, RdiObject *rdi)
{
  ::coder::SizeType loop_ub;
  if (mesh->teids.size(0) == 0) {
    m2cErrMsgIdAndTxt("rdi_create:emptyMesh", "Call setup on the input mesh");
  }
  if (mesh->coords.size(1) == mesh->topo_ndims + 1) {
    m2cErrMsgIdAndTxt("rdi_create:missingNormals",
                      "Must have normal ID for surface mesh");
  }
  rdi->degree = 2;
  rdi->nrmid = 0;
  rdi->ring = 2.0;
  rdi->A.row_ptr.set_size(0);
  rdi->A.ncols = 0;
  //  Construct A.row_ptr
  rdi->A.col_ind.set_size(0);
  rdi->A.val.set_size(0);
  rdi->nnzs.set_size(0);
  rdi->fullrank = true;
  rdi->epsbeta = 0.0002;
  rdi->cglobal = 0.05;
  rdi->clocal = 0.5;
  rdi->kappa0 = 1.0;
  rdi->kappa1 = 0.3;
  //  Step 1, compute 1-ring information
  if ((mesh->node2nodes.row_ptr.size(0) == 0) ||
      (mesh->node2elems.row_ptr.size(0) == 0)) {
    wlsmesh_compute_1ring(mesh);
  }
  //  Step 2, compute mesh properties (sizes, measures, etc)
  if (mesh->elemmeas.size(0) == 0) {
    wlsmesh_compute_meshprop(mesh, 0);
  }
  //  Step 3, compute nodal partitioning
  rdi->stclid = b_append_wlsmesh_kring(mesh);
  rdi->extstclid = c_append_wlsmesh_kring(mesh);
  //  Step 5, compute primary stencil object
  wlsmesh_compute_stencils(mesh, rdi->stclid, 2.0);
  //  Step 6, initialize rank-deficient tags
  rdi->rdtags.set_size(mesh->coords.size(0));
  loop_ub = mesh->coords.size(0);
  for (::coder::SizeType i{0}; i < loop_ub; i++) {
    rdi->rdtags[i] = false;
  }
}

// rdi_create - Create RDI object
void rdi_create(WlsMesh *mesh, ::coder::SizeType degree, RdiObject *rdi)
{
  ::coder::SizeType loop_ub;
  if (mesh->teids.size(0) == 0) {
    m2cErrMsgIdAndTxt("rdi_create:emptyMesh", "Call setup on the input mesh");
  }
  rdi->ring = degree;
  if (mesh->coords.size(1) == mesh->topo_ndims + 1) {
    m2cErrMsgIdAndTxt("rdi_create:missingNormals",
                      "Must have normal ID for surface mesh");
  }
  rdi->degree = degree;
  rdi->nrmid = 0;
  rdi->A.row_ptr.set_size(0);
  rdi->A.ncols = 0;
  //  Construct A.row_ptr
  rdi->A.col_ind.set_size(0);
  rdi->A.val.set_size(0);
  rdi->nnzs.set_size(0);
  rdi->fullrank = true;
  rdi->epsbeta = 0.0002;
  rdi->cglobal = 0.05;
  rdi->clocal = 0.5;
  rdi->kappa0 = 1.0;
  rdi->kappa1 = 0.3;
  //  Step 1, compute 1-ring information
  if ((mesh->node2nodes.row_ptr.size(0) == 0) ||
      (mesh->node2elems.row_ptr.size(0) == 0)) {
    wlsmesh_compute_1ring(mesh);
  }
  //  Step 2, compute mesh properties (sizes, measures, etc)
  if (mesh->elemmeas.size(0) == 0) {
    wlsmesh_compute_meshprop(mesh, 0);
  }
  //  Step 3, compute nodal partitioning
  rdi->stclid = b_append_wlsmesh_kring(mesh);
  rdi->extstclid = c_append_wlsmesh_kring(mesh);
  //  Step 5, compute primary stencil object
  wlsmesh_compute_stencils(mesh, rdi->stclid, static_cast<real_T>(degree));
  //  Step 6, initialize rank-deficient tags
  rdi->rdtags.set_size(mesh->coords.size(0));
  loop_ub = mesh->coords.size(0);
  for (::coder::SizeType i{0}; i < loop_ub; i++) {
    rdi->rdtags[i] = false;
  }
}

// rdi_create - Create RDI object
void rdi_create(WlsMesh *mesh, ::coder::SizeType degree, ::coder::SizeType nparts, RdiObject *rdi)
{
  ::coder::SizeType loop_ub;
  if (mesh->teids.size(0) == 0) {
    m2cErrMsgIdAndTxt("rdi_create:emptyMesh", "Call setup on the input mesh");
  }
  rdi->ring = degree;
  if (mesh->coords.size(1) == mesh->topo_ndims + 1) {
    m2cErrMsgIdAndTxt("rdi_create:missingNormals",
                      "Must have normal ID for surface mesh");
  }
  rdi->degree = degree;
  rdi->nrmid = 0;
  rdi->A.row_ptr.set_size(0);
  rdi->A.ncols = 0;
  //  Construct A.row_ptr
  rdi->A.col_ind.set_size(0);
  rdi->A.val.set_size(0);
  rdi->nnzs.set_size(0);
  rdi->fullrank = true;
  rdi->epsbeta = 0.0002;
  rdi->cglobal = 0.05;
  rdi->clocal = 0.5;
  rdi->kappa0 = 1.0;
  rdi->kappa1 = 0.3;
  //  Step 1, compute 1-ring information
  if ((mesh->node2nodes.row_ptr.size(0) == 0) ||
      (mesh->node2elems.row_ptr.size(0) == 0)) {
    wlsmesh_compute_1ring(mesh);
  }
  //  Step 2, compute mesh properties (sizes, measures, etc)
  if (mesh->elemmeas.size(0) == 0) {
    wlsmesh_compute_meshprop(mesh, 0);
  }
  //  Step 3, compute nodal partitioning
  if (nparts == 0) {
    nparts = 1;
#ifdef _OPENMP
    nparts = omp_get_max_threads();
#endif // _OPENMP
  }
  if (nparts > 0) {
    wlsmesh_compute_nodeparts(mesh, nparts);
  }
  //  Step 4, initialize stencils
  rdi->stclid = b_append_wlsmesh_kring(mesh);
  rdi->extstclid = c_append_wlsmesh_kring(mesh);
  //  Step 5, compute primary stencil object
  wlsmesh_compute_stencils(mesh, rdi->stclid, static_cast<real_T>(degree));
  //  Step 6, initialize rank-deficient tags
  rdi->rdtags.set_size(mesh->coords.size(0));
  loop_ub = mesh->coords.size(0);
  for (::coder::SizeType i{0}; i < loop_ub; i++) {
    rdi->rdtags[i] = false;
  }
}

// rdi_create - Create RDI object
void rdi_create(WlsMesh *mesh, ::coder::SizeType degree, ::coder::SizeType nparts, ::coder::SizeType nrmid,
                 RdiObject *rdi)
{
  ::coder::SizeType loop_ub;
  if (mesh->teids.size(0) == 0) {
    m2cErrMsgIdAndTxt("rdi_create:emptyMesh", "Call setup on the input mesh");
  }
  rdi->ring = degree;
  if ((mesh->coords.size(1) == mesh->topo_ndims + 1) &&
      ((nrmid < 1) || (nrmid > mesh->nrmstables.size(0)))) {
    m2cErrMsgIdAndTxt("rdi_create:missingNormals",
                      "Must have normal ID for surface mesh");
  }
  rdi->degree = degree;
  rdi->nrmid = nrmid;
  rdi->A.row_ptr.set_size(0);
  rdi->A.ncols = 0;
  //  Construct A.row_ptr
  rdi->A.col_ind.set_size(0);
  rdi->A.val.set_size(0);
  rdi->nnzs.set_size(0);
  rdi->fullrank = true;
  rdi->epsbeta = 0.0002;
  rdi->cglobal = 0.05;
  rdi->clocal = 0.5;
  rdi->kappa0 = 1.0;
  rdi->kappa1 = 0.3;
  //  Step 1, compute 1-ring information
  if ((mesh->node2nodes.row_ptr.size(0) == 0) ||
      (mesh->node2elems.row_ptr.size(0) == 0)) {
    wlsmesh_compute_1ring(mesh);
  }
  //  Step 2, compute mesh properties (sizes, measures, etc)
  if (mesh->elemmeas.size(0) == 0) {
    wlsmesh_compute_meshprop(mesh, nrmid);
  }
  //  Step 3, compute nodal partitioning
  if (nparts == 0) {
    nparts = 1;
#ifdef _OPENMP
    nparts = omp_get_max_threads();
#endif // _OPENMP
  }
  if (nparts > 0) {
    wlsmesh_compute_nodeparts(mesh, nparts);
  }
  //  Step 4, initialize stencils
  rdi->stclid = b_append_wlsmesh_kring(mesh);
  rdi->extstclid = c_append_wlsmesh_kring(mesh);
  //  Step 5, compute primary stencil object
  wlsmesh_compute_stencils(mesh, rdi->stclid, static_cast<real_T>(degree));
  //  Step 6, initialize rank-deficient tags
  rdi->rdtags.set_size(mesh->coords.size(0));
  loop_ub = mesh->coords.size(0);
  for (::coder::SizeType i{0}; i < loop_ub; i++) {
    rdi->rdtags[i] = false;
  }
}

// rdi_create - Create RDI object
void rdi_create5(WlsMesh *mesh, ::coder::SizeType degree, ::coder::SizeType nparts, ::coder::SizeType nrmid,
                 real_T ring, RdiObject *rdi)
{
  ::coder::SizeType loop_ub;
  if (mesh->teids.size(0) == 0) {
    m2cErrMsgIdAndTxt("rdi_create:emptyMesh", "Call setup on the input mesh");
  }
  if ((mesh->coords.size(1) == mesh->topo_ndims + 1) &&
      ((nrmid < 1) || (nrmid > mesh->nrmstables.size(0)))) {
    m2cErrMsgIdAndTxt("rdi_create:missingNormals",
                      "Must have normal ID for surface mesh");
  }
  rdi->degree = degree;
  rdi->nrmid = nrmid;
  rdi->ring = ring;
  rdi->A.row_ptr.set_size(0);
  rdi->A.ncols = 0;
  //  Construct A.row_ptr
  rdi->A.col_ind.set_size(0);
  rdi->A.val.set_size(0);
  rdi->nnzs.set_size(0);
  rdi->fullrank = true;
  rdi->epsbeta = 0.0002;
  rdi->cglobal = 0.05;
  rdi->clocal = 0.5;
  rdi->kappa0 = 1.0;
  rdi->kappa1 = 0.3;
  //  Step 1, compute 1-ring information
  if ((mesh->node2nodes.row_ptr.size(0) == 0) ||
      (mesh->node2elems.row_ptr.size(0) == 0)) {
    wlsmesh_compute_1ring(mesh);
  }
  //  Step 2, compute mesh properties (sizes, measures, etc)
  if (mesh->elemmeas.size(0) == 0) {
    wlsmesh_compute_meshprop(mesh, nrmid);
  }
  //  Step 3, compute nodal partitioning
  if (nparts == 0) {
    nparts = 1;
#ifdef _OPENMP
    nparts = omp_get_max_threads();
#endif // _OPENMP
  }
  if (nparts > 0) {
    wlsmesh_compute_nodeparts(mesh, nparts);
  }
  //  Step 4, initialize stencils
  rdi->stclid = b_append_wlsmesh_kring(mesh);
  rdi->extstclid = c_append_wlsmesh_kring(mesh);
  //  Step 5, compute primary stencil object
  wlsmesh_compute_stencils(mesh, rdi->stclid, ring);
  //  Step 6, initialize rank-deficient tags
  rdi->rdtags.set_size(mesh->coords.size(0));
  loop_ub = mesh->coords.size(0);
  for (::coder::SizeType i{0}; i < loop_ub; i++) {
    rdi->rdtags[i] = false;
  }
}

// rdi_create - Create RDI object
void rdi_create(WlsMesh *mesh, ::coder::SizeType degree, ::coder::SizeType nparts, ::coder::SizeType nrmid,
                 real_T ring, ::coder::SizeType varargin_2, RdiObject *rdi)
{
  boolean_T m2cTryBlkErrFlag;
  ::coder::SizeType loop_ub;
  ::coder::SizeType nthreads;
  if (mesh->teids.size(0) == 0) {
    m2cErrMsgIdAndTxt("rdi_create:emptyMesh", "Call setup on the input mesh");
  }
  if ((mesh->coords.size(1) == mesh->topo_ndims + 1) &&
      ((nrmid < 1) || (nrmid > mesh->nrmstables.size(0)))) {
    m2cErrMsgIdAndTxt("rdi_create:missingNormals",
                      "Must have normal ID for surface mesh");
  }
  rdi->degree = degree;
  rdi->nrmid = nrmid;
  rdi->ring = ring;
  rdi->A.row_ptr.set_size(0);
  rdi->A.ncols = 0;
  //  Construct A.row_ptr
  rdi->A.col_ind.set_size(0);
  rdi->A.val.set_size(0);
  rdi->nnzs.set_size(0);
  rdi->fullrank = true;
  rdi->epsbeta = 0.0002;
  rdi->cglobal = 0.05;
  rdi->clocal = 0.5;
  rdi->kappa0 = 1.0;
  rdi->kappa1 = 0.3;
  if (varargin_2 <= 0) {
    nthreads = 1;
#ifdef _OPENMP
    nthreads = omp_get_max_threads();
#endif // _OPENMP
  } else {
    nthreads = varargin_2;
  }
  m2cTryBlkErrFlag = 0;
#pragma omp parallel default(shared) num_threads(nthreads)
  try { // try
    //  Step 1, compute 1-ring information
    if ((mesh->node2nodes.row_ptr.size(0) == 0) ||
        (mesh->node2elems.row_ptr.size(0) == 0)) {
#pragma omp barrier
      b_wlsmesh_compute_1ring(mesh);
    }
    //  Step 2, compute mesh properties (sizes, measures, etc)
    if (mesh->elemmeas.size(0) == 0) {
#pragma omp barrier
      b_wlsmesh_compute_meshprop(mesh, nrmid);
    }
    //  Step 3, compute nodal partitioning
    if (nparts == 0) {
      nparts = 1;
#ifdef _OPENMP
      nparts = omp_get_max_threads();
#endif // _OPENMP
    }
#pragma omp barrier
#pragma omp single
    { // single
      if (nparts > 0) {
        wlsmesh_compute_nodeparts(mesh, nparts);
      }
    } // single
#pragma omp barrier
#pragma omp single
    { // single
      rdi->stclid = b_append_wlsmesh_kring(mesh);
      rdi->extstclid = c_append_wlsmesh_kring(mesh);
    } // single
    //  Step 5, compute primary stencil object
    b_wlsmesh_compute_stencils(mesh, rdi->stclid, ring);
#pragma omp single
    { // single
      rdi->rdtags.set_size(mesh->coords.size(0));
      loop_ub = mesh->coords.size(0);
      for (::coder::SizeType i{0}; i < loop_ub; i++) {
        rdi->rdtags[i] = false;
      }
    } // single
  } catch (const std::runtime_error &m2cExc) {
    m2cTryBlkErrFlag = 1;
    m2cPrintError("runtime_error %s\n", m2cExc.what());
  } catch (const std::logic_error &m2cExc) {
    m2cTryBlkErrFlag = 1;
    m2cPrintError("logic_error %s\n", m2cExc.what());
  } catch (const std::exception &m2cExc) {
    m2cTryBlkErrFlag = 1;
    m2cPrintError("exception %s\n", m2cExc.what());
  } catch (...) {
    m2cTryBlkErrFlag = 1;
    m2cPrintError("Unknown error detected from C++ exceptions\n");
    fflush(stderr);
  } // end try
  if ((int32_T)m2cTryBlkErrFlag != 0) {
    throw std::runtime_error("omp4m:runtimeErrorInThread");
  }
}

// rdi_create - Create RDI object
void rdi_create7(WlsMesh *mesh, ::coder::SizeType degree, ::coder::SizeType nparts, ::coder::SizeType nrmid,
                 real_T ring, RdiObject *rdi)
{
  ::coder::SizeType loop_ub;
  if (mesh->teids.size(0) == 0) {
    m2cErrMsgIdAndTxt("rdi_create:emptyMesh", "Call setup on the input mesh");
  }
  if ((mesh->coords.size(1) == mesh->topo_ndims + 1) &&
      ((nrmid < 1) || (nrmid > mesh->nrmstables.size(0)))) {
    m2cErrMsgIdAndTxt("rdi_create:missingNormals",
                      "Must have normal ID for surface mesh");
  }
  rdi->degree = degree;
  rdi->nrmid = nrmid;
  rdi->ring = ring;
  rdi->A.row_ptr.set_size(0);
  rdi->A.ncols = 0;
  //  Construct A.row_ptr
  rdi->A.col_ind.set_size(0);
  rdi->A.val.set_size(0);
  rdi->nnzs.set_size(0);
  rdi->fullrank = true;
  rdi->epsbeta = 0.0002;
  rdi->cglobal = 0.05;
  rdi->clocal = 0.5;
  rdi->kappa0 = 1.0;
  rdi->kappa1 = 0.3;
  //  Step 1, compute 1-ring information
  if ((mesh->node2nodes.row_ptr.size(0) == 0) ||
      (mesh->node2elems.row_ptr.size(0) == 0)) {
#pragma omp barrier
    b_wlsmesh_compute_1ring(mesh);
  }
  //  Step 2, compute mesh properties (sizes, measures, etc)
  if (mesh->elemmeas.size(0) == 0) {
#pragma omp barrier
    b_wlsmesh_compute_meshprop(mesh, nrmid);
  }
  //  Step 3, compute nodal partitioning
  if (nparts == 0) {
    nparts = 1;
#ifdef _OPENMP
    nparts = omp_get_max_threads();
#endif // _OPENMP
  }
#pragma omp barrier
#pragma omp single
  { // single
    if (nparts > 0) {
      wlsmesh_compute_nodeparts(mesh, nparts);
    }
  } // single
#pragma omp barrier
#pragma omp single
  { // single
    rdi->stclid = b_append_wlsmesh_kring(mesh);
    rdi->extstclid = c_append_wlsmesh_kring(mesh);
  } // single
  //  Step 5, compute primary stencil object
  b_wlsmesh_compute_stencils(mesh, rdi->stclid, ring);
#pragma omp single
  { // single
    rdi->rdtags.set_size(mesh->coords.size(0));
    loop_ub = mesh->coords.size(0);
    for (::coder::SizeType i{0}; i < loop_ub; i++) {
      rdi->rdtags[i] = false;
    }
  } // single
}

// rdi_postproc_markers - Improve disc. markers via post-processing
static inline
void rdi_postproc_markers(::coder::array<int8_T, 2U> &distags,
                           const WlsMesh *mesh)
{
  rdi_expand_markers(distags, mesh->node2nodes.row_ptr,
                     mesh->node2nodes.col_ind, 4);
  //  Step 2, shrink/contract markers
  rdi_contract_markers(distags, mesh->coords, mesh->elemtables, mesh->teids,
                       mesh->node2nodes.row_ptr, mesh->node2nodes.col_ind,
                       mesh->node2elems.row_ptr, mesh->node2elems.col_ind, 4);
}

// rdi_postproc_markers - Improve disc. markers via post-processing
void rdi_postproc_markers(::coder::array<int8_T, 2U> &distags,
                           const WlsMesh *mesh, ::coder::SizeType nlayers)
{
  if (nlayers < 0) {
    nlayers = 4;
  }
  //  Step 1, expand markers
  rdi_expand_markers(distags, mesh->node2nodes.row_ptr,
                     mesh->node2nodes.col_ind, nlayers);
  //  Step 2, shrink/contract markers
  rdi_contract_markers(distags, mesh->coords, mesh->elemtables, mesh->teids,
                       mesh->node2nodes.row_ptr, mesh->node2nodes.col_ind,
                       mesh->node2elems.row_ptr, mesh->node2elems.col_ind,
                       nlayers);
}

} // namespace rdi_kernel

// End of code generation (librdi.cpp)
