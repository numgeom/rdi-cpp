/*

Copyright (c) 2022, NumGeom Group at Stony Brook University
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this
   list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

*/

/*!
 * @file rdi.cpp
 * @brief Implementations of RDI
 * @authors Qiao Chen
 */

#include <algorithm>
#include <cmath>
#include <cstdarg>
#include <cstdio>
#include <stdexcept>
#include <vector>

#ifdef RDI_USE_MPI
#  include <mpi.h>
#endif

#ifdef RDI_USE_METIS
#  define OMP4M_HAS_METIS
#endif

#include "coder_array.h"

#include "rdi.hpp"

// internal generated sources

// parameters
#include "codegen_src/librdi.h"

// stencils
#include "codegen_src/rdi_compute_stencils.h"

namespace rdi {

using rdi_kernel::RdiMesh;  // defined in codegen_src/librdi_types.h

namespace internal {

// error handling, consistent with rdi_error
inline static void rdi_error(const char *msg, ...) {
  va_list args1;
  va_start(args1, msg);
  va_list args2;
  va_copy(args2, args1);
  std::vector<char> buf(1 + std::vsnprintf(nullptr, 0, msg, args1));
  va_end(args1);
  std::vsnprintf(buf.data(), buf.size(), msg, args2);
  va_end(args2);
#ifndef RDI_USE_MPI
  throw std::runtime_error(buf.data());
#else
  std::fprintf(stderr, buf.data());
  std::fflush(stderr);
  MPI_Abort(MPI_COMM_WORLD, 1);
#endif
}

// assertion
inline static void rdi_assert(const bool flag, const char *msg = nullptr) {
#ifndef NDEBUG
  if (!flag) rdi_error(msg ? msg : "assertion failed");
#else
  (void)flag;
  (void)msg;
#endif
}

// transpose data
template <class T>
inline static void transpose_c2f(const T *v, const int m, const int n,
                                 T *v_out) {
  for (int j = 0; j < n; ++j) {
    v += j;
    v_out += j * m;
    for (int i = 0; i < m; ++i) v_out[i] = v[i * n];
  }
}

}  // namespace internal

// default parameters
void set_default_params(const int topo_dim, RdiParams &params,
                        const int nthreads) {
  rdi_kernel::rdi_default_params(topo_dim, nthreads, &params);
}

// default constructor
RDI::RDI() {
  // set all array attributes
  _mesh.reset(new RdiMesh);
  _mesh->hGlobal  = 0.0;
  _mesh->surfType = 0;
  _stencils.reset(new ::coder::array<int, 2>());
  _M.ncols = 0;
  _M.idx.indptr.reset(new ::coder::array<int, 1>());
  _M.idx.indices.reset(new ::coder::array<int, 1>());
  _M.vals.reset(new ::coder::array<double, 1>());
}

// constructor with mesh input
RDI::RDI(const double *xs, const int *conn, const int attrs[],
         const RdiParams &params, const double *dirs)
    : RDI() {
  _init(xs, conn, attrs, params, dirs);
}

RDI::~RDI() {}

// check empty
bool RDI::empty() const { return _mesh->xs.size(0) == 0; }

// get geo dimension
int RDI::dim() const { return _mesh->xs.size(1); }

// get topo dimension
int RDI::topo_dim() const { return _mesh->topoDim; }

// get number of nodes
int RDI::nnodes() const { return _mesh->xs.size(0); }

// get number of cells
int RDI::ncells() const { return _mesh->conn.size(0); }

// get number of nodes per cell
int RDI::nnpcell() const { return _mesh->conn.size(1); }

// get global h
double RDI::global_h() const { return _mesh->hGlobal; }

// set global h
double &RDI::global_h() { return _mesh->hGlobal; }

// get data to cell weights
const double *RDI::cell_weights() const { return _mesh->cellWeights.data(); }

// get data to cell weights
double *RDI::cell_weights() { return _mesh->cellWeights.data(); }

// get data to cell sizes
const double *RDI::cell_sizes() const { return _mesh->cellSizes.data(); }

// get data to cell sizes
double *RDI::cell_sizes() { return _mesh->cellSizes.data(); }

// compute the OSUS operator
void RDI::compute_osusop(RdiParams &params) {
  // first, compute the operator assuming no adaptivity is needed
  if (params.verbose > 0) std::printf("Start computing the OSUS operator...\n");
  if (_M.vals->size(0)) {
    if (params.verbose > 1)
      std::printf(" Clear previous OSUS operator first...\n");
    _M.ncols = 0;
    // maybe use set_size is better for reusing same RDI?
    // _M.idx.indptr->clear();
    // _M.idx.indices->clear();
    // _M.vals->clear();
    _M.idx.indptr->set_size(0);
    _M.idx.indices->set_size(0);
    _M.vals->set_size(0);
  }
  ::coder::array<int, 1> rd_nodes, nnz_pr;  // rank deficient nodes and nnz pr
  _compute_stencils(params, rd_nodes);      // build stencils
  _assemble_osusop(params, nnz_pr, rd_nodes);  // assemble operator
  // then, if we have rd nodes, continue with adaptivity
  int adap_steps(0);
  while (rd_nodes.size(0) != 0) {
    params.ring += 1.0 / params.dim;
    ++adap_steps;
    if (params.verbose > 0)
      std::printf(
          "%d rd nodes were detected, updating the op with ring %g...\n",
          rd_nodes.size(0), params.ring);
    _compute_stencils(params, rd_nodes);         // build stencils
    _assemble_osusop(params, nnz_pr, rd_nodes);  // assemble operator
  }
  if (params.verbose > 0) {
    if (!adap_steps)
      std::printf("Finish computing the OSUS operator...\n");
    else
      std::printf("Finish computing the OSUS operator with %d adap. steps...\n",
                  adap_steps);
  }
  _M.ncols = nnodes();
}

// compute OSUS indicators
void RDI::compute_osusind(const double *fs, const RdiParams &params,
                          double *alpha_cell, double *alpha_node) const {
  if (empty() || _M.vals->size(0) == 0)
    internal::rdi_error("an empty RDI or empty OSUS operator!");
  ::coder::array<double, 2> fs_wrap, alpha_cell_wrap, alpha_node_wrap;
  fs_wrap.set((double *)fs, nnodes(), 1);
  alpha_cell_wrap.set(alpha_cell, ncells(), 1);
  alpha_node_wrap.set(alpha_node, nnodes(), 1);

  int nthreads = params.nThreads;
#ifdef _OPENMP
  if (nthreads < 1) nthreads = omp_get_max_threads();
#endif
#pragma omp parallel default(shared) num_threads(nthreads)
  {
    rdi_kernel::rdi_compute_osusind(*_M.idx.indptr, *_M.idx.indices, *_M.vals,
                                    fs_wrap, _mesh.get(), &params,
                                    alpha_cell_wrap, alpha_node_wrap);
  }
}

// compute osc. indicators
void RDI::compute_oscind(const double df_global, const double *alpha_cell,
                         const RdiParams &params, double *beta) const {
  if (empty()) internal::rdi_error("an empty RDI!");
  ::coder::array<double, 2> df_wrap, alpha_cell_wrap, beta_wrap;
  df_wrap.set((double *)&df_global, 1, 1);
  alpha_cell_wrap.set((double *)alpha_cell, ncells(), 1);
  beta_wrap.set(beta, nnodes(), 1);

  int nthreads = params.nThreads;
#ifdef _OPENMP
  if (nthreads < 1) nthreads = omp_get_max_threads();
#endif
#pragma omp parallel default(shared) num_threads(nthreads)
  {
    rdi_kernel::rdi_compute_oscind(df_wrap, alpha_cell_wrap, _mesh.get(),
                                   &params, beta_wrap);
  }
}

void RDI::mark_discontinuities(const double *fs, const double df_global,
                               const double *alpha_cell, const double *beta,
                               const RdiParams &params,
                               signed char *    dis_tags) const {
  if (empty()) internal::rdi_error("an empty RDI!");
  ::coder::array<double, 2>      df_wrap, fs_wrap, alpha_cell_wrap, beta_wrap;
  ::coder::array<signed char, 2> dis_tags_wrap;
  df_wrap.set((double *)&df_global, 1, 1);
  fs_wrap.set((double *)fs, nnodes(), 1);
  alpha_cell_wrap.set((double *)alpha_cell, ncells(), 1);
  beta_wrap.set((double *)beta, nnodes(), 1);
  dis_tags_wrap.set(dis_tags, nnodes(), 1);
  rdi_kernel::rdi_mark_discontinuities(fs_wrap, alpha_cell_wrap, beta_wrap,
                                       df_wrap, _mesh.get(), &params,
                                       dis_tags_wrap);
}

void RDI::compute_inds(const double *fs, const double df_global,
                       const RdiParams &params, double *alpha_cell,
                       double *alpha_node, double *beta) const {
  if (empty() || _M.vals->size(0) == 0)
    internal::rdi_error("an empty RDI or empty OSUS operator!");

  ::coder::array<double, 2> fs_wrap, df_wrap, alpha_cell_wrap, alpha_node_wrap,
      beta_wrap;
  df_wrap.set((double *)&df_global, 1, 1);
  fs_wrap.set((double *)fs, nnodes(), 1);
  alpha_cell_wrap.set(alpha_cell, ncells(), 1);
  alpha_node_wrap.set(alpha_node, nnodes(), 1);

  int nthreads = params.nThreads;
#ifdef _OPENMP
  if (nthreads < 1) nthreads = omp_get_max_threads();
#endif
#pragma omp parallel default(shared) num_threads(nthreads)
  {
    rdi_kernel::rdi_compute_inds(*_M.idx.indptr, *_M.idx.indices, *_M.vals,
                                 fs_wrap, df_wrap, _mesh.get(), &params,
                                 alpha_cell_wrap, alpha_node_wrap, beta_wrap);
  }
}

void RDI::_init(const double *xs, const int *conn, const int attrs[],
                const RdiParams &params, const double *dirs) {
  if (params.verbose > 0) std::printf("Begin setting up the mesh...\n");

  // get attributes of the mesh
  const int n = attrs[0], dim = attrs[1], m = attrs[2], npc = attrs[3];
  // check dimension
  if (dim < 1 || dim > 3) internal::rdi_error("bad dimension value %d", dim);
  if (params.dim > dim)
    internal::rdi_error("bad topological dimension %d", params.dim);

  if (n <= 0 || m <= 0 || npc <= 0) {
    if (params.verbose > 0)
      std::printf(
          "Invalid nnodes (%d), ncells (%d), and/or nnpcell (%d) detected, do "
          "nothing...\n",
          n, m, npc);
    return;
  }

  if (params.dim < dim && !dirs && !params.surfType)
    internal::rdi_error("missing normal or tangent vector fields");
  _mesh->topoDim  = params.dim;
  _mesh->surfType = params.surfType;

  // setup the mesh
  _mesh->xs.set((double *)xs, n, dim);
  internal::rdi_assert(!_mesh->xs.is_owner(), "cannot be owner");
  _mesh->conn.set((int *)conn, m, npc);
  internal::rdi_assert(!_mesh->conn.is_owner(), "cannot be owner");
  if (dirs) {
    _mesh->dirs.set((double *)dirs, n, dim);
    internal::rdi_assert(!_mesh->dirs.is_owner(), "cannot be owner");
  }

  // compute node to cell
  if (params.verbose > 1)
    std::printf(" Build node-to-cell adj. structure...\n");
  _build_node2cell();

  // compute node to node
  if (params.verbose > 1)
    std::printf(" Build node-to-node adj. structure...\n");
  _build_node2node();

  // compute cell weights
  _compute_cellweights(params);

  // compute cell sizes
  _compute_cellsizes(params);
  _mesh->hGlobal =
      *std::max_element(_mesh->cellSizes.data(), _mesh->cellSizes.data() + m);

  // compute the partitioning
  if (params.verbose > 1)
    std::printf(" Compute recursive nodal partitioning...\n");
  rdi_kernel::rdi_partition(&params, _mesh.get());

  if (params.verbose > 0) std::printf("Finish setting up the mesh...\n");
}

// build node-to-cell
void RDI::_build_node2cell() {
  rdi_kernel::rdi_build_node2cell(nnodes(), _mesh->conn, _mesh->n2cPtr,
                                  _mesh->n2cList);
}

// build node-to-node
void RDI::_build_node2node() {
  rdi_kernel::rdi_build_node2node(nnodes(), _mesh->conn, _mesh->topoDim,
                                  _mesh->n2cPtr, _mesh->n2cList, _mesh->n2nPtr,
                                  _mesh->n2nList);
}

// compute cell weights
void RDI::_compute_cellweights(const RdiParams &params) {
  rdi_kernel::rdi_compute_cellweights(_mesh->xs, _mesh->conn, &params,
                                      _mesh->cellWeights);
}

// compute cell sizes
void RDI::_compute_cellsizes(const RdiParams &params) {
  rdi_kernel::rdi_compute_cellsizes(_mesh->xs, _mesh->conn, _mesh->n2cPtr,
                                    _mesh->n2cList, &params, _mesh->dirs,
                                    _mesh->cellSizes);
}

// compute stencils
// this function is overridable
void RDI::_compute_stencils(const RdiParams &             params,
                            const ::coder::array<int, 1> &rd_nodes) {
  // use our default impl with AHF, this is not feature-aware for surface
  rdi_stencils::rdi_compute_stencils(nnodes(), _mesh->conn, &params, rd_nodes,
                                   *_stencils);
}

// assemble stencils
void RDI::_assemble_osusop(const RdiParams &       params,
                           ::coder::array<int, 1> &nnz_pr,
                           ::coder::array<int, 1> &rd_nodes) {
  rdi_kernel::rdi_assemble_osusop(_mesh.get(), *_stencils, &params,
                                  *_M.idx.indptr, *_M.idx.indices, *_M.vals,
                                  nnz_pr, rd_nodes);
}

}  // namespace rdi
