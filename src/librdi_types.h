// Copyright 2022 The NumGeom Group, Stony Brook University
//
// librdi_types.h
//
// Code generation for function 'rdi_create'
//

#ifndef LIBRDI_TYPES_H
#define LIBRDI_TYPES_H

// Include files
#include "rtwtypes.h"
#ifndef CODER_ARRAY_SIZE_TYPE_DEFINED
#define CODER_ARRAY_SIZE_TYPE_DEFINED
namespace coder {
#ifdef M2C_SIZETYPE64
typedef int64_T SizeType;
#else
typedef int32_T SizeType;
#endif // M2C_SIZETYPE64
} // namespace coder
#endif // CODER_ARRAY_SIZE_TYPE_DEFINED

#include "coder_array.h"
#define MAX_THREADS omp_get_max_threads()

// Type Definitions
namespace rdi_kernel {
struct ConnData {
  int32_T etype;
  ::coder::array<char_T, 2U> name;
  int32_T istart;
  ::coder::array<int32_T, 2U> conn;
};

struct Omp4mPart {
  int32_T nparts;
  ::coder::array<int32_T, 1U> part_ptr;
  ::coder::array<int32_T, 1U> part_list;
  ::coder::array<int32_T, 1U> shared_ents;
};

struct CrsGraph {
  ::coder::array<int64_T, 1U> row_ptr;
  ::coder::array<int32_T, 1U> col_ind;
  int32_T ncols;
};

struct Stencils {
  ::coder::array<char_T, 2U> name;
  ::coder::array<boolean_T, 1U> reflected;
  ::coder::array<int32_T, 1U> vidmap;
  CrsGraph ngbverts;
  CrsGraph ngbelems;
};

struct CrsMatrix {
  ::coder::array<int64_T, 1U> row_ptr;
  ::coder::array<int32_T, 1U> col_ind;
  ::coder::array<real_T, 1U> val;
  int32_T ncols;
};

struct RdiObject {
  int32_T degree;
  int32_T nrmid;
  int32_T stclid;
  int32_T extstclid;
  real_T ring;
  CrsMatrix A;
  ::coder::array<int32_T, 1U> nnzs;
  ::coder::array<boolean_T, 1U> rdtags;
  boolean_T fullrank;
  real_T epsbeta;
  real_T cglobal;
  real_T clocal;
  real_T kappa0;
  real_T kappa1;
};

struct NodeSet {
  ::coder::array<char_T, 2U> name;
  ::coder::array<int32_T, 1U> nids;
};

struct ElementSet {
  ::coder::array<char_T, 2U> name;
  ::coder::array<int32_T, 1U> eids;
};

struct FacetSet {
  ::coder::array<char_T, 2U> name;
  ::coder::array<int32_T, 1U> eids;
  ::coder::array<int8_T, 1U> facets;
};

struct EdgeSet {
  ::coder::array<char_T, 2U> name;
  ::coder::array<int32_T, 1U> eids;
  ::coder::array<int8_T, 1U> edges;
};

struct NormalsData {
  int8_T nrmtype;
  int32_T esetidx;
  ::coder::array<char_T, 2U> name;
  ::coder::array<real_T, 2U> normals;
  ::coder::array<int32_T, 2U> indices;
};

struct WlsMesh {
  ::coder::array<real_T, 2U> coords;
  ::coder::array<ConnData, 1U> elemtables;
  ::coder::array<uint64_T, 1U> teids;
  int32_T topo_ndims;
  ::coder::array<NodeSet, 1U> nodesets;
  ::coder::array<ElementSet, 1U> elemsets;
  ::coder::array<FacetSet, 1U> facetsets;
  ::coder::array<EdgeSet, 1U> edgesets;
  ::coder::array<NormalsData, 1U> nrmstables;
  ::coder::array<Omp4mPart, 1U> elemparts;
  ::coder::array<uint64_T, 2U> sibhfs;
  ::coder::array<uint64_T, 1U> v2hfid;
  ::coder::array<boolean_T, 1U> bwork1;
  ::coder::array<boolean_T, 1U> bwork2;
  ::coder::array<int32_T, 1U> iwork1;
  ::coder::array<int32_T, 1U> iwork2;
  CrsGraph node2nodes;
  CrsGraph node2elems;
  ::coder::array<Stencils, 1U> stencils;
  ::coder::array<int32_T, 1U> bridges;
  ::coder::array<Omp4mPart, 1U> nodeparts;
  ::coder::array<real_T, 1U> elemmeas;
  ::coder::array<real_T, 1U> elemh;
  ::coder::array<real_T, 1U> nodeh;
  real_T globalh;
};

} // namespace rdi_kernel

#endif
// End of code generation (librdi_types.h)
