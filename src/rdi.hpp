/*

Copyright (c) 2022, The NumGeom Group at Stony Brook University
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this
   list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.

3. Neither the name of the copyright holder nor the names of its
   contributors may be used to endorse or promote products derived from
   this software without specific prior written permission.

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
 * @file rdi.hpp
 * @brief The interface of Robust Discontinuity Indicators (RDI) library
 * @authors Qiao Chen
 */

#ifndef RDI_HPP_
#define RDI_HPP_

#include <memory>

#include "rdi_params.hpp"

// declaration of coder::array
namespace coder {

template <typename T, int N>
class array;

}  // namespace coder


// definition of parameter structure
namespace rdi_kernel {

// declaration of mesh object
struct RdiMesh;

}

namespace rdi {

using rdi_kernel::RdiParams;

/*!
 * @brief Set a parameter structure to tis default values
 *
 * @param[in] topo_dim Topological dimension, ranging from 1 to 3
 * @param[out] params Upon output, @a params will be set to default values
 * @param[in] nthreads (optional) Number of parallel threads
 *
 * Notice that the topological dimension (i.e., @a topo_dim ) represents the
 * dimension of the cells, and it can be smaller than the actual geometric
 * dimension. For instance, for a surface in 3D, its topological dimension is
 * just 2. In RDI, for surface computations in 2D and 3D, set the corresponding
 * topological dimensions to be 1D and 2D, respectively.
 */
void set_default_params(const int topo_dim, RdiParams &params,
                        const int nthreads = 0);

/*!
 * @brief Set parameters for spherical surface (sphsurf) for 3D
 *
 * @sa See set_set_default_params() for more details
 */
inline void set_sphsurf_params(RdiParams &params, const int nthreads = 0) {
  set_default_params(2, params, nthreads);
  params.surfType = 1;  // spherical surface type
}

/*!
 * @struct AdjStruct
 * @brief Adjacent structure for indices
 *
 */
struct AdjStruct {
  std::unique_ptr<::coder::array<int, 1>> indptr;   ///< index pointer array
  std::unique_ptr<::coder::array<int, 1>> indices;  ///< list of indices
};

/*!
 * @struct CRS
 * @brief CRS (Compressed Row Structure) formart for sparse matrix
 *
 */
struct CRS {
  int                                        ncols;  ///< number of columns
  AdjStruct                                  idx;    ///< index adjacencies
  std::unique_ptr<::coder::array<double, 1>> vals;   ///< numerical values
};

/*!
 * @class RDI
 * @brief RDI class definition
 *
 * The core algorithms are implemented as member functions in the @a RDI class
 * in which each instance holds an unstructured mesh, i.e., coordinates and
 * connectivity. Given an mesh, the RDI will automatically computes the
 * corresponding essential mesh-related fields, such as cell weights, cell
 * sizes, etc.
 *
 * @code{.cpp}
 *    #include "rdi.hpp"
 *
 *    using namespace rdi;
 *
 *    int main() {
 *      RdiParams params;
 *      int dim = 2;
 *      get_default_params(dim, params); // set up parameters
 *      // pionts of the unit square 4-by-2
 *      const double xs[] = {0.0, 0.0,
 *                           1.0, 0.0,
 *                           1.0, 1.0,
 *                           0.0, 1.0};
 *      // 2 triangular cells
 *      // 4---3
 *      // | / |
 *      // 1---2
 *      const int conn[] = {1, 2, 3,
 *                          1, 3, 4};
 *      // mesh attribues
 *      const int attrs[] = {4,  // four nodes
 *                           2,  // geometric dimension is 2
 *                           2,  // two cells
 *                           3}; // (maximum) number of nodes per cell
 *      RDI rdi_obj(xs, conn, attrs, params); // constructor RDI instance
 *    }
 * @endcode
 */
class RDI {
 public:
  /*!
   * @brief Construct an empty RDI instance
   *
   */
  RDI();

  /*!
   * @brief Construct a RDI object with a given mesh
   *
   * @param[in] xs Coordinate array
   * @param[in] conn Connectivity array
   * @param[in] attrs Mesh attributes, at least length of six
   * @param[in] params Control parameters
   * @param[in] dirs (optional) Directional vector field for surface/curve
   *
   * @sa For more details regarding the input parameters, see set_mesh()
   */
  RDI(const double *xs, const int *conn, const int attrs[],
      const RdiParams &params, const double *dirs = nullptr);

  RDI(const RDI &other) = delete;  // disable copy constructor

  virtual ~RDI();

  /*!
   * @brief Set up a mesh to a RDI instance
   *
   * @param[in] xs Coordinate array
   * @param[in] conn Connectivity array
   * @param[in] attrs Mesh attributes, at least length of six
   * @param[in] params Control parameters
   * @param[in] dirs (optional) Directional vector field for surface/curve
   *
   * For the mesh attributes (i.e., @a attrs ), it must contain six essential
   * information fields:
   *  - attrs[0]: Number of nodes in @a xs
   *  - attrs[1]: Geometric dimension of the mesh
   *  - attrs[2]: Number of cells in @a conn
   *  - attrs[3]: Maximum number of nodes per cell
   *
   * As we can see, @a attrs[0]*attrs[1] gives us the total length of @a xs
   * array, and @a attrs[2]*attrs[3] gives us the total length of @a conn array.
   * An important note is the ordering of @a xs and @a conn arrays: In RDI, we
   * assume row-major arrays for the meshes, i.e., @a xs is stored as
   * @a {x1,y1,x2,y2,...} and @a conn is stored as
   * @a {n11,n12,n13,...,nx1,nx2,nx3,...}, where @a nxy means the @a y-th node
   * in the @a x-th cell.
   *
   * In addition, we assume one-based indexing system for @a conn, i.e., the
   * first node is labelled as 1. For mixed mesh, simply pad unused entries in
   * @a conn with values <=0.
   *
   * The optional input argument is the directional vector field, which is not
   * needed for body computation (i.e., topological and geometric dimensions are
   * the same). For surface computations, (i.e., @a params.dim=attrs[1]-1 ),
   * then @a dirs must be the unit normal vectors of each node. For curve
   * computation in 3D (i.e., @a params.dim==1 and @a attrs[1]==3 ), @a dirs
   * must be the unit tangent vectors of each node. Note that if @a dirs is
   * given, then it must be the same size as the @a xs field does. One
   * exception for surface computation is that @a dirs can be omitted if the
   * surface is unit sphere in 3D or unit circle in 2D; for more details,
   * see @ref RdiParams::surfType
   */
  inline void set_mesh(const double *xs, const int *conn, const int attrs[],
                       const RdiParams &params, const double *dirs = nullptr) {
    _init(xs, conn, attrs, params, dirs);
  }

  // utilities

  /*!
   * @brief Check if a RDI instance is empty (mesh hasn't been set up)
   *
   */
  bool empty() const;

  /*!
   * @brief Return the geometric dimension
   *
   */
  int dim() const;

  /*!
   * @brief Return the topological dimension
   *
   */
  int topo_dim() const;

  /*!
   * @brief Return the total number of nodes
   *
   */
  int nnodes() const;

  /*!
   * @brief Return the total number of cells
   *
   */
  int ncells() const;

  /*!
   * @brief Return the maximum number of nodes per cell
   *
   */
  int nnpcell() const;

  /*!
   * @brief Return the maximum global mesh size @f$h_{\text{global}}@f$
   *
   */
  double global_h() const;

  /*!
   * @brief Update the global mesh size
   *
   */
  double &global_h();

  /*!
   * @brief Return the cell weights array _cell_wgts (length is @ref ncells())
   *
   * @note If @ref empty() is @a true, then @a nullptr is returned
   */
  const double *cell_weights() const;

  /*!
   * @brief Return the cell weights array _cell_wgts (length is @ref ncells())
   *
   * @note If @ref empty() is @a true, then @a nullptr is returned
   */
  double *cell_weights();

  /*!
   * @brief Return the local cell sizes (length is @ref ncells())
   *
   * @note If @ref empty() is @a true, then @a nullptr is returned
   */
  const double *cell_sizes() const;

  /*!
   * @brief Return the local cell sizes (length is @ref ncells())
   *
   * @note If @ref empty() is @a true, then @a nullptr is returned
   */
  double *cell_sizes();

  /*!
   * @brief Compute the over- and under-shoot (OSUS) operator
   *
   * @param[in,out] params Control parameters
   *
   * Upon input, @a params.ring is the initial ring size. Upon output,
   * @a params.ring is the final ring size if adaptivity occurred in computing
   * the USOS operator.
   *
   * @sa compute_osusind() uses the OSUS operator to build OSUS indicators.
   */
  void compute_osusop(RdiParams &params);

  /*!
   * @brief Computing the cell-based and nodal OSUS indicators
   *
   * @param[in] fs User input discrete function value at nodes
   * @param[in] params Control parameters
   * @param[out] alpha_cell Cell-based OSUS indicator values
   * @param[out] alpha_node Nodal OSUS indicator values
   *
   * This function computes the OSUS indicators, which can be used to determine
   * OSUS cells and used to compute WLS-ENO weights. Regarding input parameters
   * of this function, @a fs must be at least length @ref nnodes(), and similar
   * for the output nodal OSUS indicators @a alpha_node. For output cell-based
   * OSUS indicators @a alpha_cell, it must be at least @ref ncells()
   *
   * Note that this function is designed for scalar function fields in which
   * each function field is stored contiguous in memory.
   *
   * @sa This function must be called after compute_osusop(). Also,
   *     @a alpha_cell is used in mark_discontinuities() to mark
   *     discontinuity nodes.
   */
  void compute_osusind(const double *fs, const RdiParams &params,
                       double *alpha_cell, double *alpha_node) const;

  /*!
   * @brief Compute the oscillation indicators
   *
   * @param[in] df_global Global range of a discrete functions
   * @param[in] alpha_cell Cell-based OSUS indicators
   * @param[in] params Control parameters
   * @param[out] beta Oscillation indicators
   *
   * This function compute the oscillation indicators (@f$\beta@f$). Regarding
   * the input parameters of this function, @a df_global approxmates the global
   * range of each target functions, i.e., @f$\max f - \min f@f$. The output
   * oscillation indicators @a beta has size of @ref nnodes()
   *
   * Note that this function is designed for scalar function fields in which
   * each function field is stored contiguous in memory.
   *
   * @sa See compute_osusind() for computing @a alpha_cell, and also see
   *     mark_discontinuities() for using @a beta in marking discontinuity
   *     nodes.
   */
  void compute_oscind(const double df_global, const double *alpha_cell,
                      const RdiParams &params, double *beta) const;

  /*!
   * @brief Mark the discontinuity nodes
   *
   * @param[in] fs User input discrete function value at nodes
   * @param[in] df_global Global range of all @a fs
   * @param[in] alpha_cell Cell-based OSUS indicators
   * @param[in] beta Nodal oscillation indicators
   * @param[in] params Control parameters
   * @param[out] dis_tags Discontinuity tags
   *
   * This function marks and categorizes the discontinuity nodes. Regarding the
   * input parameters of this function, @a fs and @a beta must be at least
   * length @ref nnodes(). @a df_global approxmates the global
   * range of the function field @a fs. @a alpha_cell must be at least
   * length @ref ncells(). The output @a dis_tags are nodal value
   * thus its length is at least @ref nnodes()
   *
   * As to the output @a dis_tags, smooth nodes are marked as 0, @f$C_0@f$
   * discontinuity (shocks) nodes are marked as 1, and @f$C_1@f$ discontinuity
   * (whose direvatives are shocks) are marked as 2.
   *
   * Note that this function is designed for scalar function fields in which
   * each function field is stored contiguous in memory.
   *
   * @sa See compute_osusind() and compute_oscind() for computations of OSUS
   *     indicators @a alpha_cell and oscillation indicators @a beta.
   */
  void mark_discontinuities(const double *fs, const double df_global,
                            const double *alpha_cell, const double *beta,
                            const RdiParams &params,
                            signed char *    dis_tags) const;

  /*!
   * @brief Compute all indicators at once
   */
  void compute_inds(const double *fs, const double df_global,
                    const RdiParams &params, double *alpha_cell,
                    double *alpha_node, double *beta) const;

 protected:
  /*!
   * @brief Initialized the mesh database in a RDI instance
   *
   * This function is called inside both constructor and @ref set_mesh()
   */
  void _init(const double *xs, const int *conn, const int attrs[],
             const RdiParams &params, const double *dirs);

  /*!
   * @brief Compute the node-to-cell adjacency structure
   *
   * This function builds an adj. structure that stores the 1-ring neighbor
   * cells of each node.
   *
   * @note This function is called inside @ref _init()
   */
  void _build_node2cell();

  /*!
   * @brief Compute the node-to-node adjacency structure (1-ring)
   *
   * This function builds and adj. structure that stores the 1-ring neighbor
   * nodes of each node.
   *
   * @note This function is called inside @ref _init()
   */
  void _build_node2node();

  /*!
   * @brief Compute the cell sizes (@f$h@f$)
   *
   * This function computes averaged local @f$h@f$ for each cell
   *
   * @note This function is called inside @ref _init()
   */
  void _compute_cellsizes(const RdiParams &params);

  /*!
   * @brief Compute the cell weights (e.g., areas in 2D)
   *
   * @note This function is called inside @ref _init()
   */
  void _compute_cellweights(const RdiParams &params);

  /*!
   * @brief Compute or update the k-ring stencil of each node
   *
   * This function is overridable so that one can build the stencil in either
   * way they want. Note that after calling this function, _stencils must be
   * set up properly. For the i-th node, its stencil is a compressed list of
   * near-by nodal IDs like @a {i,n1,n2,n3,...}, and node i must be the first
   * one in its stencil
   *
   * If @a rd_nodes is not empty, then stencils corresponding to nodal IDs in
   * @a rd_nodes need to be updated (or recomputed).
   *
   * @note This function is called inside @ref compute_osusop()
   */
  virtual void _compute_stencils(const RdiParams &             params,
                                 const ::coder::array<int, 1> &rd_nodes);

  /*!
   * @brief Assemble or update the OSUS operator given nodal stencils
   *
   * This function essentially builds the CRS OSUS operator @ref _M
   *
   * @note This function is called inside @ref compute_osusop()
   */
  void _assemble_osusop(const RdiParams &params, ::coder::array<int, 1> &nnz_pr,
                        ::coder::array<int, 1> &rd_nodes);

 protected:
  std::unique_ptr<rdi_kernel::RdiMesh> _mesh;  ///< mesh object
  ///< node-to-cell adj. structure
  std::unique_ptr<::coder::array<int, 2>> _stencils;
  ///< stencils structure
  CRS _M;  ///< CRS for OSUS operator
};

}  // namespace rdi

#endif  // RDI_HPP_
