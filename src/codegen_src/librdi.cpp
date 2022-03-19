// Copyright 2022 The NumGeom Group, Stony Brook University
// Main developers:
//     rdilib: Qiao Chen
//     momp2cpp: Xiangmin Jiao, Qiao Chen
//     wlslib: Xiangmin Jiao, Qiao Chen, Jacob Jones
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
#ifdef OMP4M_HAS_METIS
#include "metis.h"
#endif
#include "wls_lapack.hpp"
#include <cmath>
#include <cstdio>
#include <cstring>
#ifdef _OPENMP
#include <omp.h>
#endif
#include <stdexcept>
#include <stdio.h>

// Type Definitions
namespace rdi_kernel
{
  struct emxArray_real_T_3x1
  {
    double data[3];
    int size[2];
  };

  struct emxArray_char_T_7x1
  {
    char data[7];
    int size[2];
  };

  struct emxArray_real_T_1x512
  {
    double data[512];
    int size[2];
  };

  struct Omp4mPartContext
  {
    ::coder::array<int, 1U> atomic_counters;
  };

  struct emxArray_boolean_T_0
  {
    int size[1];
  };

  struct emxArray_real_T_0
  {
    int size[1];
  };

  struct WlsWeight
  {
    emxArray_char_T_7x1 name;
    emxArray_real_T_0 params_shared;
    emxArray_real_T_1x512 params_pointwise;
    emxArray_boolean_T_0 omit_rows;
  };

  struct emxArray_real_T_0x0
  {
    int size[2];
  };

  struct WlsObject
  {
    int npoints;
    int degree;
    int order;
    int interp0;
    boolean_T use_dag;
    int stride;
    ::coder::array<double, 2U> us;
    emxArray_real_T_3x1 origin;
    ::coder::array<double, 1U> rweights;
    emxArray_real_T_3x1 hs_inv;
    ::coder::array<unsigned char, 2U> dag;
    ::coder::array<double, 2U> V;
    ::coder::array<double, 2U> QR;
    ::coder::array<double, 2U> vdops;
    int nrows;
    int ncols;
    int rank;
    ::coder::array<int, 1U> jpvt;
    ::coder::array<double, 1U> work;
    boolean_T rowmajor;
    emxArray_real_T_0x0 QRt;
  };

  struct WlsDataStruct
  {
    WlsObject wlsObj;
    WlsWeight wlsWgts;
    ::coder::array<double, 2U> coeffs;
    boolean_T interp0;
    boolean_T useDag;
  };
}

// Variable Definitions
namespace rdi_kernel
{
  static const signed char iv[12]{ 5, 10, 15, 20, 15, 35, 55, 75, 30, 60, 80,
    100 };

  static const double dv[7]{ 333.33333333333331, 1000.0, 3333.3333333333335,
    10000.0, 100000.0, 1.0E+6, 1.0E+7 };
}

// Function Declarations
namespace rdi_kernel
{
  static inline
  void assemble_body(const ::coder::array<double, 2U> &mesh_xs, const ::
    coder::array<int, 2U> &mesh_conn, const ::coder::array<int, 1U> &mesh_n2cPtr,
    const ::coder::array<int, 1U> &mesh_n2cList, const ::coder::array<Omp4mPart,
    1U> &mesh_parts, const ::coder::array<int, 2U> &stcls, const ::coder::array<
    int, 1U> &rowPtr, const ::coder::array<int, 1U> &colInd, ::coder::array<
    double, 1U> &vals, const ::coder::array<int, 1U> &nnzPr, ::coder::array<int,
    1U> &rdNodes, const RdiParams *params);
  static inline
  int assemble_body_kernel(const ::coder::array<double, 2U> &xs, const ::
    coder::array<int, 2U> &stcls, int degree, int iend, WlsDataStruct *wls);
  static inline
  void assemble_body_par(const ::coder::array<double, 2U> &xs, const ::
    coder::array<int, 2U> &conn, const ::coder::array<int, 2U> &stcls, const ::
    coder::array<int, 1U> &n2cPtr, const ::coder::array<int, 1U> &n2cList, int
    degree, boolean_T interp0, boolean_T useDag, const ::coder::array<int, 1U>
    &rowPtr, const ::coder::array<int, 1U> &colInd, ::coder::array<double, 1U>
    &vals, const ::coder::array<int, 1U> &nnzPr, const ::coder::array<Omp4mPart,
    1U> &parts, ::coder::array<WlsDataStruct, 1U> &wls, ::coder::array<boolean_T,
    1U> &rdTags, ::coder::array<int, 1U> &rdCounts, Omp4mPartContext
    *partContext);
  static inline
  void assemble_body_task(const ::coder::array<double, 2U> &xs, const ::
    coder::array<int, 2U> &conn, const ::coder::array<int, 2U> &stcls, const ::
    coder::array<int, 1U> &n2cPtr, const ::coder::array<int, 1U> &n2cList, int
    degree, boolean_T interp0, boolean_T useDag, const ::coder::array<int, 1U>
    &rowPtr, const ::coder::array<int, 1U> &colInd, ::coder::array<double, 1U>
    &vals, const ::coder::array<int, 1U> &nnzPr, const ::coder::array<Omp4mPart,
    1U> &parts, ::coder::array<WlsDataStruct, 1U> &wls, ::coder::array<boolean_T,
    1U> &rdTags, ::coder::array<int, 1U> &rdCounts);
  static inline
  void assemble_body_task_kernel(const ::coder::array<double, 2U> &xs,
    const ::coder::array<int, 2U> &conn, const ::coder::array<int, 2U> &stcls,
    const ::coder::array<int, 1U> &n2cPtr, const ::coder::array<int, 1U>
    &n2cList, int degree, const ::coder::array<int, 1U> &rowPtr, const ::coder::
    array<int, 1U> &colInd, ::coder::array<double, 1U> &vals, const ::coder::
    array<int, 1U> &nnzPr, const ::coder::array<Omp4mPart, 1U> &parts, int
    mypart, int lvl, ::coder::array<WlsDataStruct, 1U> &wls, ::coder::array<
    boolean_T, 1U> &rdTags, ::coder::array<int, 1U> &rdCounts);
  static inline
  void assemble_surf(const ::coder::array<double, 2U> &mesh_xs, const ::
    coder::array<int, 2U> &mesh_conn, const ::coder::array<double, 2U>
    &mesh_dirs, const ::coder::array<int, 1U> &mesh_n2cPtr, const ::coder::array<
    int, 1U> &mesh_n2cList, const ::coder::array<Omp4mPart, 1U> &mesh_parts,
    const ::coder::array<int, 2U> &stcls, const ::coder::array<int, 1U> &rowPtr,
    const ::coder::array<int, 1U> &colInd, ::coder::array<double, 1U> &vals,
    const ::coder::array<int, 1U> &nnzPr, ::coder::array<int, 1U> &rdNodes,
    const RdiParams *params);
  static inline
  int assemble_surf_kernel(const ::coder::array<double, 2U> &xs, const ::
    coder::array<int, 2U> &conn, const ::coder::array<int, 2U> &stcls, const ::
    coder::array<int, 1U> &n2cPtr, const ::coder::array<int, 1U> &n2cList, int
    surfType, int degree, const ::coder::array<double, 2U> &nrms, const ::coder::
    array<int, 1U> &rowPtr, const ::coder::array<int, 1U> &colInd, ::coder::
    array<double, 1U> &vals, const ::coder::array<int, 1U> &nnzPr, const ::coder::
    array<int, 1U> &nrange, int iend, WlsDataStruct *wls, ::coder::array<
    boolean_T, 1U> &rdTags);
  static inline
  void assemble_surf_par(const ::coder::array<double, 2U> &xs, const ::
    coder::array<int, 2U> &conn, const ::coder::array<int, 2U> &stcls, const ::
    coder::array<int, 1U> &n2cPtr, const ::coder::array<int, 1U> &n2cList, int
    surfType, int degree, boolean_T interp0, boolean_T useDag, const ::coder::
    array<double, 2U> &nrms, const ::coder::array<int, 1U> &rowPtr, const ::
    coder::array<int, 1U> &colInd, ::coder::array<double, 1U> &vals, const ::
    coder::array<int, 1U> &nnzPr, const ::coder::array<Omp4mPart, 1U> &parts, ::
    coder::array<WlsDataStruct, 1U> &wls, ::coder::array<boolean_T, 1U> &rdTags,
    ::coder::array<int, 1U> &rdCounts, Omp4mPartContext *partContext);
  static inline
  void assemble_surf_task(const ::coder::array<double, 2U> &xs, const ::
    coder::array<int, 2U> &conn, const ::coder::array<int, 2U> &stcls, const ::
    coder::array<int, 1U> &n2cPtr, const ::coder::array<int, 1U> &n2cList, int
    surfType, int degree, boolean_T interp0, boolean_T useDag, const ::coder::
    array<double, 2U> &nrms, const ::coder::array<int, 1U> &rowPtr, const ::
    coder::array<int, 1U> &colInd, ::coder::array<double, 1U> &vals, const ::
    coder::array<int, 1U> &nnzPr, const ::coder::array<Omp4mPart, 1U> &parts, ::
    coder::array<WlsDataStruct, 1U> &wls, ::coder::array<boolean_T, 1U> &rdTags,
    ::coder::array<int, 1U> &rdCounts);
  static inline
  void assemble_surf_task_kernel(const ::coder::array<double, 2U> &xs,
    const ::coder::array<int, 2U> &conn, const ::coder::array<int, 2U> &stcls,
    const ::coder::array<int, 1U> &n2cPtr, const ::coder::array<int, 1U>
    &n2cList, int surfType, int degree, const ::coder::array<double, 2U> &nrms,
    const ::coder::array<int, 1U> &rowPtr, const ::coder::array<int, 1U> &colInd,
    ::coder::array<double, 1U> &vals, const ::coder::array<int, 1U> &nnzPr,
    const ::coder::array<Omp4mPart, 1U> &parts, int mypart, int lvl, ::coder::
    array<WlsDataStruct, 1U> &wls, ::coder::array<boolean_T, 1U> &rdTags, ::
    coder::array<int, 1U> &rdCounts);
  static inline
  void b_WlsDataStruct(int degree, boolean_T interp0, boolean_T useDag,
    WlsObject *wls_wlsObj, WlsWeight *wls_wlsWgts, ::coder::array<double, 2U>
    &wls_coeffs, boolean_T *wls_interp0, boolean_T *wls_useDag);
  static inline
  int b_assemble_body_kernel(const ::coder::array<double, 2U> &xs, const ::
    coder::array<int, 2U> &stcls, int degree, int iend, WlsDataStruct *wls);
  static inline
  void b_assemble_body_par(const ::coder::array<double, 2U> &xs, const ::
    coder::array<int, 2U> &conn, const ::coder::array<int, 2U> &stcls, const ::
    coder::array<int, 1U> &n2cPtr, const ::coder::array<int, 1U> &n2cList, int
    degree, boolean_T interp0, boolean_T useDag, const ::coder::array<int, 1U>
    &rowPtr, const ::coder::array<int, 1U> &colInd, ::coder::array<double, 1U>
    &vals, const ::coder::array<int, 1U> &nnzPr, const ::coder::array<Omp4mPart,
    1U> &parts, ::coder::array<WlsDataStruct, 1U> &wls, ::coder::array<boolean_T,
    1U> &rdTags, ::coder::array<int, 1U> &rdCounts, Omp4mPartContext
    *partContext);
  static inline
  void b_assemble_body_task(const ::coder::array<double, 2U> &xs, const ::
    coder::array<int, 2U> &conn, const ::coder::array<int, 2U> &stcls, const ::
    coder::array<int, 1U> &n2cPtr, const ::coder::array<int, 1U> &n2cList, int
    degree, boolean_T interp0, boolean_T useDag, const ::coder::array<int, 1U>
    &rowPtr, const ::coder::array<int, 1U> &colInd, ::coder::array<double, 1U>
    &vals, const ::coder::array<int, 1U> &nnzPr, const ::coder::array<Omp4mPart,
    1U> &parts, ::coder::array<WlsDataStruct, 1U> &wls, ::coder::array<boolean_T,
    1U> &rdTags, ::coder::array<int, 1U> &rdCounts);
  static inline
  void b_assemble_body_task_kernel(const ::coder::array<double, 2U> &xs,
    const ::coder::array<int, 2U> &conn, const ::coder::array<int, 2U> &stcls,
    const ::coder::array<int, 1U> &n2cPtr, const ::coder::array<int, 1U>
    &n2cList, int degree, const ::coder::array<int, 1U> &rowPtr, const ::coder::
    array<int, 1U> &colInd, ::coder::array<double, 1U> &vals, const ::coder::
    array<int, 1U> &nnzPr, const ::coder::array<Omp4mPart, 1U> &parts, int
    mypart, int lvl, ::coder::array<WlsDataStruct, 1U> &wls, ::coder::array<
    boolean_T, 1U> &rdTags, ::coder::array<int, 1U> &rdCounts);
  static inline
  void b_assemble_surf_kernel(const ::coder::array<double, 2U> &xs, const
    ::coder::array<int, 2U> &conn, const ::coder::array<int, 2U> &stcls, const ::
    coder::array<int, 1U> &n2cPtr, const ::coder::array<int, 1U> &n2cList, int
    surfType, int degree, const ::coder::array<double, 2U> &nrms, const ::coder::
    array<int, 1U> &rowPtr, const ::coder::array<int, 1U> &colInd, ::coder::
    array<double, 1U> &vals, const ::coder::array<int, 1U> &nnzPr, const ::coder::
    array<int, 1U> &nrange, int istart, int iend, WlsDataStruct *wls, ::coder::
    array<boolean_T, 1U> &rdTags, int *rdCounts);
  static inline
  void b_assemble_surf_par(const ::coder::array<double, 2U> &xs, const ::
    coder::array<int, 2U> &conn, const ::coder::array<int, 2U> &stcls, const ::
    coder::array<int, 1U> &n2cPtr, const ::coder::array<int, 1U> &n2cList, int
    surfType, int degree, boolean_T interp0, boolean_T useDag, const ::coder::
    array<double, 2U> &nrms, const ::coder::array<int, 1U> &rowPtr, const ::
    coder::array<int, 1U> &colInd, ::coder::array<double, 1U> &vals, const ::
    coder::array<int, 1U> &nnzPr, const ::coder::array<Omp4mPart, 1U> &parts, ::
    coder::array<WlsDataStruct, 1U> &wls, ::coder::array<boolean_T, 1U> &rdTags,
    ::coder::array<int, 1U> &rdCounts, Omp4mPartContext *partContext);
  static inline
  void b_assemble_surf_task(const ::coder::array<double, 2U> &xs, const ::
    coder::array<int, 2U> &conn, const ::coder::array<int, 2U> &stcls, const ::
    coder::array<int, 1U> &n2cPtr, const ::coder::array<int, 1U> &n2cList, int
    surfType, int degree, boolean_T interp0, boolean_T useDag, const ::coder::
    array<double, 2U> &nrms, const ::coder::array<int, 1U> &rowPtr, const ::
    coder::array<int, 1U> &colInd, ::coder::array<double, 1U> &vals, const ::
    coder::array<int, 1U> &nnzPr, const ::coder::array<Omp4mPart, 1U> &parts, ::
    coder::array<WlsDataStruct, 1U> &wls, ::coder::array<boolean_T, 1U> &rdTags,
    ::coder::array<int, 1U> &rdCounts);
  static inline
  void b_assemble_surf_task_kernel(const ::coder::array<double, 2U> &xs,
    const ::coder::array<int, 2U> &conn, const ::coder::array<int, 2U> &stcls,
    const ::coder::array<int, 1U> &n2cPtr, const ::coder::array<int, 1U>
    &n2cList, int surfType, int degree, const ::coder::array<double, 2U> &nrms,
    const ::coder::array<int, 1U> &rowPtr, const ::coder::array<int, 1U> &colInd,
    ::coder::array<double, 1U> &vals, const ::coder::array<int, 1U> &nnzPr,
    const ::coder::array<Omp4mPart, 1U> &parts, int mypart, int lvl, ::coder::
    array<WlsDataStruct, 1U> &wls, ::coder::array<boolean_T, 1U> &rdTags, ::
    coder::array<int, 1U> &rdCounts);
  static inline
  void b_wls_init(WlsObject *wls, const double us_data[], int us_size,
    const char weight_name_data[], const double weight_params_pointwise_data[],
    const int weight_params_pointwise_size[2], int degree, boolean_T interp0,
    boolean_T use_dag, int npoints);
  static inline
  void b_wls_resize(WlsObject *wls, int npoints, int degree, boolean_T
    use_dag);
  static inline
  void build_part(int nParts, const ::coder::array<int, 1U> &nparts,
    const ::coder::array<int, 1U> &cparts, const ::coder::array<int, 1U> &eptr,
    const ::coder::array<int, 1U> &eind, ::coder::array<boolean_T, 1U> &ctags, ::
    coder::array<int, 1U> &iwork, ::coder::array<int, 1U> &partptr, ::coder::
    array<int, 1U> &partlist, ::coder::array<int, 1U> &sharedents);
  static inline
  int c_assemble_body_kernel(const ::coder::array<double, 2U> &xs, const ::
    coder::array<int, 2U> &stcls, int degree, int iend, WlsDataStruct *wls);
  static inline
  void c_assemble_body_par(const ::coder::array<double, 2U> &xs, const ::
    coder::array<int, 2U> &conn, const ::coder::array<int, 2U> &stcls, const ::
    coder::array<int, 1U> &n2cPtr, const ::coder::array<int, 1U> &n2cList, int
    degree, boolean_T interp0, boolean_T useDag, const ::coder::array<int, 1U>
    &rowPtr, const ::coder::array<int, 1U> &colInd, ::coder::array<double, 1U>
    &vals, const ::coder::array<int, 1U> &nnzPr, const ::coder::array<Omp4mPart,
    1U> &parts, ::coder::array<WlsDataStruct, 1U> &wls, ::coder::array<boolean_T,
    1U> &rdTags, ::coder::array<int, 1U> &rdCounts, Omp4mPartContext
    *partContext);
  static inline
  void c_assemble_body_task(const ::coder::array<double, 2U> &xs, const ::
    coder::array<int, 2U> &conn, const ::coder::array<int, 2U> &stcls, const ::
    coder::array<int, 1U> &n2cPtr, const ::coder::array<int, 1U> &n2cList, int
    degree, boolean_T interp0, boolean_T useDag, const ::coder::array<int, 1U>
    &rowPtr, const ::coder::array<int, 1U> &colInd, ::coder::array<double, 1U>
    &vals, const ::coder::array<int, 1U> &nnzPr, const ::coder::array<Omp4mPart,
    1U> &parts, ::coder::array<WlsDataStruct, 1U> &wls, ::coder::array<boolean_T,
    1U> &rdTags, ::coder::array<int, 1U> &rdCounts);
  static inline
  void c_assemble_body_task_kernel(const ::coder::array<double, 2U> &xs,
    const ::coder::array<int, 2U> &conn, const ::coder::array<int, 2U> &stcls,
    const ::coder::array<int, 1U> &n2cPtr, const ::coder::array<int, 1U>
    &n2cList, int degree, const ::coder::array<int, 1U> &rowPtr, const ::coder::
    array<int, 1U> &colInd, ::coder::array<double, 1U> &vals, const ::coder::
    array<int, 1U> &nnzPr, const ::coder::array<Omp4mPart, 1U> &parts, int
    mypart, int lvl, ::coder::array<WlsDataStruct, 1U> &wls, ::coder::array<
    boolean_T, 1U> &rdTags, ::coder::array<int, 1U> &rdCounts);
  static inline
  int c_assemble_surf_kernel(const ::coder::array<double, 2U> &xs, const ::
    coder::array<int, 2U> &conn, const ::coder::array<int, 2U> &stcls, const ::
    coder::array<int, 1U> &n2cPtr, const ::coder::array<int, 1U> &n2cList, int
    surfType, int degree, const ::coder::array<double, 2U> &nrms, const ::coder::
    array<int, 1U> &rowPtr, const ::coder::array<int, 1U> &colInd, ::coder::
    array<double, 1U> &vals, const ::coder::array<int, 1U> &nnzPr, const ::coder::
    array<int, 1U> &nrange, int iend, WlsDataStruct *wls, ::coder::array<
    boolean_T, 1U> &rdTags);
  static inline
  void call_metis_mesh(int n, ::coder::array<int, 1U> &eptr, ::coder::
    array<int, 1U> &eind, int nParts, ::coder::array<int, 1U> &nparts, ::coder::
    array<int, 1U> &cparts);
  static inline
  void compute_area(const ::coder::array<double, 2U> &xs, const ::coder::
    array<int, 2U> &conn, ::coder::array<double, 1U> &A, int nThreads);
  static inline
  void compute_area_kernel(int m, int dim, const ::coder::array<double,
    2U> &xs, const ::coder::array<int, 2U> &conn, ::coder::array<double, 1U> &A);
  static inline
  void compute_area_kernel(int m, int dim, const ::coder::array<double,
    2U> &xs, const ::coder::array<int, 2U> &conn, ::coder::array<double, 1U> &A,
    int nThreads);
  static inline
  void compute_beta_kernel(int n, double epsBeta, double hGlobal, const ::
    coder::array<double, 2U> &dfGlobal, const ::coder::array<double, 2U>
    &alphaCell, const ::coder::array<double, 1U> &mesh_cellWeights, const ::
    coder::array<int, 1U> &mesh_n2cPtr, const ::coder::array<int, 1U>
    &mesh_n2cList, int nRhs, ::coder::array<double, 2U> &beta);
  static inline
  void compute_body_h(const ::coder::array<double, 2U> &xs, const ::coder::
    array<int, 2U> &conn, ::coder::array<double, 1U> &h, int nThreads);
  static inline
  void compute_body_h_kernel(int m, int dim, const ::coder::array<double,
    2U> &xs, const ::coder::array<int, 2U> &conn, ::coder::array<double, 1U> &h,
    int nThreads);
  static inline
  void compute_body_h_kernel(int m, int dim, const ::coder::array<double,
    2U> &xs, const ::coder::array<int, 2U> &conn, ::coder::array<double, 1U> &h);
  static inline
  void compute_nodal_alpha(int n, const ::coder::array<double, 2U>
    &alphaCell, const ::coder::array<int, 1U> &mesh_n2cPtr, const ::coder::array<
    int, 1U> &mesh_n2cList, int nRhs, ::coder::array<double, 2U> &alphaNode);
  static inline
  void compute_surf_h(const ::coder::array<double, 2U> &xs, const ::coder::
    array<int, 2U> &conn, const ::coder::array<int, 1U> &n2nPtr, const ::coder::
    array<int, 1U> &n2nList, int surfType, ::coder::array<double, 1U> &h, int
    nThreads);
  static inline
  void compute_surf_h(const ::coder::array<double, 2U> &xs, const ::coder::
    array<int, 2U> &conn, const ::coder::array<int, 1U> &n2nPtr, const ::coder::
    array<int, 1U> &n2nList, const ::coder::array<double, 2U> &nrms, int
    surfType, ::coder::array<double, 1U> &h, int nThreads);
  static inline
  void compute_surf_h_kernel(int n, int m, int dim, const ::coder::array<
    double, 2U> &xs, const ::coder::array<int, 2U> &conn, const ::coder::array<
    int, 1U> &n2nPtr, const ::coder::array<int, 1U> &n2nList, int surfType, ::
    coder::array<double, 1U> &h, ::coder::array<double, 1U> &buf);
  static inline
  void compute_surf_h_kernel(int n, int m, int dim, const ::coder::array<
    double, 2U> &xs, const ::coder::array<int, 2U> &conn, const ::coder::array<
    int, 1U> &n2nPtr, const ::coder::array<int, 1U> &n2nList, int surfType, ::
    coder::array<double, 1U> &h, ::coder::array<double, 1U> &buf, int nThreads);
  static inline
  void compute_surf_h_kernel(int n, int m, int dim, const ::coder::array<
    double, 2U> &xs, const ::coder::array<int, 2U> &conn, const ::coder::array<
    int, 1U> &n2nPtr, const ::coder::array<int, 1U> &n2nList, const ::coder::
    array<double, 2U> &nrms, int surfType, ::coder::array<double, 1U> &h, ::
    coder::array<double, 1U> &buf);
  static inline
  void compute_surf_h_kernel(int n, int m, int dim, const ::coder::array<
    double, 2U> &xs, const ::coder::array<int, 2U> &conn, const ::coder::array<
    int, 1U> &n2nPtr, const ::coder::array<int, 1U> &n2nList, const ::coder::
    array<double, 2U> &nrms, int surfType, ::coder::array<double, 1U> &h, ::
    coder::array<double, 1U> &buf, int nThreads);
  static inline
  void compute_volume_tet(const ::coder::array<double, 2U> &xs, const ::
    coder::array<int, 2U> &conn, ::coder::array<double, 1U> &V, int nThreads);
  static inline
  void compute_volume_tet_kernel(int m, const ::coder::array<double, 2U>
    &xs, const ::coder::array<int, 2U> &conn, ::coder::array<double, 1U> &V, int
    nThreads);
  static inline
  void compute_volume_tet_kernel(int m, const ::coder::array<double, 2U>
    &xs, const ::coder::array<int, 2U> &conn, ::coder::array<double, 1U> &V);
  static inline
  void crsAx_kernel(const ::coder::array<int, 1U> &row_ptr, const ::coder::
    array<int, 1U> &col_ind, const ::coder::array<double, 1U> &val, int nrows,
    const ::coder::array<double, 2U> &x, int nrhs, ::coder::array<double, 2U> &b);
  static inline
  void crs_prod_mat_vec(const ::coder::array<int, 1U> &A_rowptr, const ::
    coder::array<int, 1U> &A_colind, const ::coder::array<double, 1U> &A_val,
    const ::coder::array<double, 2U> &x, ::coder::array<double, 2U> &b);
  static inline
  void d_assemble_body_kernel(const ::coder::array<double, 2U> &xs, const
    ::coder::array<int, 2U> &conn, const ::coder::array<int, 2U> &stcls, const ::
    coder::array<int, 1U> &n2cPtr, const ::coder::array<int, 1U> &n2cList, int
    degree, const ::coder::array<int, 1U> &rowPtr, const ::coder::array<int, 1U>
    &colInd, ::coder::array<double, 1U> &vals, const ::coder::array<int, 1U>
    &nnzPr, const ::coder::array<int, 1U> &nrange, int istart, int iend,
    WlsDataStruct *wls, ::coder::array<boolean_T, 1U> &rdTags, int *rdCounts);
  static inline
  void d_assemble_surf_kernel(const ::coder::array<double, 2U> &xs, const
    ::coder::array<int, 2U> &conn, const ::coder::array<int, 2U> &stcls, const ::
    coder::array<int, 1U> &n2cPtr, const ::coder::array<int, 1U> &n2cList, int
    surfType, int degree, const ::coder::array<double, 2U> &nrms, const ::coder::
    array<int, 1U> &rowPtr, const ::coder::array<int, 1U> &colInd, ::coder::
    array<double, 1U> &vals, const ::coder::array<int, 1U> &nnzPr, const ::coder::
    array<int, 1U> &nrange, int istart, int iend, WlsDataStruct *wls, ::coder::
    array<boolean_T, 1U> &rdTags, int *rdCounts);
  static inline
  void e_assemble_body_kernel(const ::coder::array<double, 2U> &xs, const
    ::coder::array<int, 2U> &conn, const ::coder::array<int, 2U> &stcls, const ::
    coder::array<int, 1U> &n2cPtr, const ::coder::array<int, 1U> &n2cList, int
    degree, const ::coder::array<int, 1U> &rowPtr, const ::coder::array<int, 1U>
    &colInd, ::coder::array<double, 1U> &vals, const ::coder::array<int, 1U>
    &nnzPr, const ::coder::array<int, 1U> &nrange, int istart, int iend,
    WlsDataStruct *wls, ::coder::array<boolean_T, 1U> &rdTags, int *rdCounts);
  static inline
  void extract_sub(int n, const ::coder::array<int, 1U> &crange, const ::
    coder::array<int, 1U> &eptr, const ::coder::array<int, 1U> &eind, ::coder::
    array<int, 1U> &iwork, ::coder::array<boolean_T, 1U> &ntags, ::coder::array<
    int, 1U> &eptrloc, ::coder::array<int, 1U> &eindloc, int *nnodes);
  static inline
  void f_assemble_body_kernel(const ::coder::array<double, 2U> &xs, const
    ::coder::array<int, 2U> &conn, const ::coder::array<int, 2U> &stcls, const ::
    coder::array<int, 1U> &n2cPtr, const ::coder::array<int, 1U> &n2cList, int
    degree, const ::coder::array<int, 1U> &rowPtr, const ::coder::array<int, 1U>
    &colInd, ::coder::array<double, 1U> &vals, const ::coder::array<int, 1U>
    &nnzPr, const ::coder::array<int, 1U> &nrange, int istart, int iend,
    WlsDataStruct *wls, ::coder::array<boolean_T, 1U> &rdTags, int *rdCounts);
  static inline
  double find_kth_shortest_dist(::coder::array<double, 1U> &arr, int k,
    int l, int r);
  static inline
  int g_assemble_body_kernel(const ::coder::array<double, 2U> &xs, const ::
    coder::array<int, 2U> &conn, const ::coder::array<int, 2U> &stcls, const ::
    coder::array<int, 1U> &n2cPtr, const ::coder::array<int, 1U> &n2cList, int
    degree, const ::coder::array<int, 1U> &rowPtr, const ::coder::array<int, 1U>
    &colInd, ::coder::array<double, 1U> &vals, const ::coder::array<int, 1U>
    &nnzPr, const ::coder::array<int, 1U> &nrange, int istart, int iend,
    WlsDataStruct *wls, ::coder::array<boolean_T, 1U> &rdTags);
  static inline
  void gen_vander(const ::coder::array<double, 2U> &us, int npoints, int
    degree, ::coder::array<double, 2U> &V);
  static inline
  void gen_vander(const ::coder::array<double, 2U> &us, int npoints, int
    degree, const ::coder::array<double, 1U> &weights, ::coder::array<double, 2U>
    &V);
  static inline
  void gen_vander_1d_dag(int degree, ::coder::array<unsigned char, 1U>
    &dag);
  static inline
  void gen_vander_2d(const ::coder::array<double, 2U> &us, int npoints,
    int degree, ::coder::array<double, 2U> &V);
  static inline
  void gen_vander_2d(const ::coder::array<double, 2U> &us, int npoints,
    int degree, const ::coder::array<double, 1U> &weights, ::coder::array<double,
    2U> &V);
  static inline
  void gen_vander_2d_dag(int degree, ::coder::array<unsigned char, 2U>
    &dag);
  static inline
  void gen_vander_3d(const ::coder::array<double, 2U> &us, int npoints,
    int degree, const ::coder::array<double, 1U> &weights, ::coder::array<double,
    2U> &V);
  static inline
  void gen_vander_3d(const ::coder::array<double, 2U> &us, int npoints,
    int degree, ::coder::array<double, 2U> &V);
  static inline
  void gen_vander_3d_dag(int degree, ::coder::array<unsigned char, 2U>
    &dag);
  static inline
  int h_assemble_body_kernel(const ::coder::array<double, 2U> &xs, const ::
    coder::array<int, 2U> &conn, const ::coder::array<int, 2U> &stcls, const ::
    coder::array<int, 1U> &n2cPtr, const ::coder::array<int, 1U> &n2cList, int
    degree, const ::coder::array<int, 1U> &rowPtr, const ::coder::array<int, 1U>
    &colInd, ::coder::array<double, 1U> &vals, const ::coder::array<int, 1U>
    &nnzPr, const ::coder::array<int, 1U> &nrange, int istart, int iend,
    WlsDataStruct *wls, ::coder::array<boolean_T, 1U> &rdTags);
  static inline
  int i_assemble_body_kernel(const ::coder::array<double, 2U> &xs, const ::
    coder::array<int, 2U> &conn, const ::coder::array<int, 2U> &stcls, const ::
    coder::array<int, 1U> &n2cPtr, const ::coder::array<int, 1U> &n2cList, int
    degree, const ::coder::array<int, 1U> &rowPtr, const ::coder::array<int, 1U>
    &colInd, ::coder::array<double, 1U> &vals, const ::coder::array<int, 1U>
    &nnzPr, const ::coder::array<int, 1U> &nrange, int istart, int iend,
    WlsDataStruct *wls, ::coder::array<boolean_T, 1U> &rdTags);
  static inline
  void init_osusop(const ::coder::array<int, 2U> &conn, const ::coder::
    array<int, 2U> &stcls, int maxNnzPr, ::coder::array<int, 1U> &rowPtr, ::
    coder::array<int, 1U> &colInd, ::coder::array<double, 1U> &vals, ::coder::
    array<int, 1U> &nnzPr);
  static inline
  void mark_kernel(int n, double hGlobal, double cGlobal, double cLocal,
    double kappa1, double kappa0, const ::coder::array<double, 2U> &fs, const ::
    coder::array<double, 2U> &alphaCell, const ::coder::array<double, 2U> &beta,
    const ::coder::array<double, 2U> &dfGlobal, const ::coder::array<int, 2U>
    &mesh_conn, const ::coder::array<double, 1U> &mesh_cellSizes, const ::coder::
    array<int, 1U> &mesh_n2cPtr, const ::coder::array<int, 1U> &mesh_n2cList,
    int nRhs, ::coder::array<signed char, 2U> &disTags);
  static inline
  void mark_kernel(int n, double hGlobal, double cGlobal, double cLocal,
    double kappa1, double kappa0, const ::coder::array<double, 2U> &fs, const ::
    coder::array<double, 2U> &alphaCell, const ::coder::array<double, 2U> &beta,
    const ::coder::array<double, 2U> &dfGlobal, const ::coder::array<int, 2U>
    &mesh_conn, const ::coder::array<double, 1U> &mesh_cellSizes, const ::coder::
    array<int, 1U> &mesh_n2cPtr, const ::coder::array<int, 1U> &mesh_n2cList,
    const ::coder::array<int, 1U> &mesh_n2nPtr, const ::coder::array<int, 1U>
    &mesh_n2nList, int nRhs, ::coder::array<signed char, 2U> &disTags);
  static inline
  void omp4mRecurPartMesh(int n, const ::coder::array<int, 2U> &cells,
    int dim, int nLevels, int nParts, ::coder::array<Omp4mPart, 1U> &parts);
  static inline
  void rdi_compute_oscind(const ::coder::array<double, 2U> &dfGlobal,
    const ::coder::array<double, 2U> &alphaCell, const ::coder::array<double, 2U>
    &mesh_xs, double mesh_hGlobal, const ::coder::array<double, 1U>
    &mesh_cellWeights, const ::coder::array<int, 1U> &mesh_n2cPtr, const ::coder::
    array<int, 1U> &mesh_n2cList, int params_dim, double params_epsBeta, ::coder::
    array<double, 2U> &beta);
  static inline
  void rdi_compute_osusind(const ::coder::array<int, 1U> &rowPtr, const ::
    coder::array<int, 1U> &colInd, const ::coder::array<double, 1U> &vals, const
    ::coder::array<double, 2U> &fs, const ::coder::array<int, 1U> &mesh_n2cPtr,
    const ::coder::array<int, 1U> &mesh_n2cList, ::coder::array<double, 2U>
    &alphaCell, ::coder::array<double, 2U> &alphaNode);
  static inline
  void rrqr_factor(const ::coder::array<double, 2U> &A, double thres, int
    rowoffset, int coloffset, int m, int n, ::coder::array<double, 2U> &QR, ::
    coder::array<int, 1U> &p, int *rank, ::coder::array<double, 1U> &work);
  static inline
  void rrqr_qmulti(const ::coder::array<double, 2U> &QR, int m, int n,
    int rank, ::coder::array<double, 2U> &bs, int nrhs, ::coder::array<double,
    1U> &work);
  static inline
  void rrqr_rtsolve(const ::coder::array<double, 2U> &QR, int n, int rank,
    ::coder::array<double, 2U> &bs, int nrhs);
  static inline
  void update_osusop(int n, const ::coder::array<int, 2U> &stcls, const ::
    coder::array<int, 1U> &nrange, const ::coder::array<int, 1U> &n2cPtr, const ::
    coder::array<int, 1U> &n2cList, ::coder::array<int, 1U> &rowPtr, ::coder::
    array<int, 1U> &colInd, ::coder::array<double, 1U> &vals, ::coder::array<int,
    1U> &nnzPr);
  static inline
  void wls_buhmann_weights(const ::coder::array<double, 2U> &us, int
    npoints, int degree, const double params_pw_data[], const int
    params_pw_size[2], ::coder::array<double, 1U> &ws);
  static inline
  void wls_func(WlsObject *wls, const double pnts_data[], const int
                       pnts_size[2], int npoints, ::coder::array<double, 2U>
                       &vdops);
  static inline
  void wls_func(WlsObject *wls, const double pnts_data[], int npoints, ::
                       coder::array<double, 2U> &vdops);
  static inline
  void wls_init(WlsObject *wls, const double us_data[], const int
                       us_size[2], const char weight_name_data[], const double
                       weight_params_pointwise_data[], const int
                       weight_params_pointwise_size[2], int degree, boolean_T
                       interp0, boolean_T use_dag, int npoints);
  static inline
  void wls_invdist_weights(const ::coder::array<double, 2U> &us, int
    npoints, int degree, const double params_pw_data[], const int
    params_pw_size[2], ::coder::array<double, 1U> &ws);
  static inline
  void wls_invdist_weights(const ::coder::array<double, 2U> &us, int
    npoints, double degree, ::coder::array<double, 1U> &ws);
  static inline
  void wls_resize(WlsObject *wls, int dim, int npoints, int degree,
    boolean_T use_dag);
}

// Function Definitions
namespace rdi_kernel
{
  static void assemble_body(const ::coder::array<double, 2U> &mesh_xs, const ::
    coder::array<int, 2U> &mesh_conn, const ::coder::array<int, 1U> &mesh_n2cPtr,
    const ::coder::array<int, 1U> &mesh_n2cList, const ::coder::array<Omp4mPart,
    1U> &mesh_parts, const ::coder::array<int, 2U> &stcls, const ::coder::array<
    int, 1U> &rowPtr, const ::coder::array<int, 1U> &colInd, ::coder::array<
    double, 1U> &vals, const ::coder::array<int, 1U> &nnzPr, ::coder::array<int,
    1U> &rdNodes, const RdiParams *params)
  {
    ::coder::array<WlsDataStruct, 1U> wls_;
    ::coder::array<int, 1U> rdCounts_;
    ::coder::array<boolean_T, 1U> rdTags_;
    Omp4mPartContext b_mesh_xs;
    WlsDataStruct b_wls_;
    int i;
    int loop_ub;
    int nrd;

    //  Body {1,2,3}-D assmebly
    rdTags_.set_size(mesh_xs.size(0));
    loop_ub = mesh_xs.size(0);
    for (i = 0; i < loop_ub; i++) {
      rdTags_[i] = false;
    }

    if (rdNodes.size(0) != 0) {
      //  adaptivity is needed for resolving RD nodes
      b_WlsDataStruct(params->degree, params->wlsInterp0, params->wlsUseDag,
                      &b_wls_.wlsObj, &b_wls_.wlsWgts, b_wls_.coeffs,
                      &b_wls_.interp0, &b_wls_.useDag);
      if (mesh_xs.size(1) == 3) {
        loop_ub = assemble_body_kernel(mesh_xs, stcls, params->degree,
          rdNodes.size(0), &b_wls_);
        rdCounts_.set_size(1);
        rdCounts_[0] = loop_ub;
      } else if (mesh_xs.size(1) == 2) {
        loop_ub = b_assemble_body_kernel(mesh_xs, stcls, params->degree,
          rdNodes.size(0), &b_wls_);
        rdCounts_.set_size(1);
        rdCounts_[0] = loop_ub;
      } else {
        loop_ub = c_assemble_body_kernel(mesh_xs, stcls, params->degree,
          rdNodes.size(0), &b_wls_);
        rdCounts_.set_size(1);
        rdCounts_[0] = loop_ub;
      }
    } else {
      wls_.set_size(params->nThreads);
      rdCounts_.set_size(params->nThreads);
      if (!params->parTask) {
        int m2cTryBlkErrCode;
        loop_ub = params->nThreads;

        m2cTryBlkErrCode = (0);
        if (params->nThreads <= 0) {
          loop_ub = 1;

#ifdef _OPENMP

          loop_ub = omp_get_max_threads();

#endif //_OPENMP

        }

#pragma omp parallel default(shared) num_threads(loop_ub)

        try { //try
          if (mesh_xs.size(1) == 3) {
            assemble_body_par(mesh_xs, mesh_conn, stcls, mesh_n2cPtr,
                              mesh_n2cList, params->degree, params->wlsInterp0,
                              params->wlsUseDag, rowPtr, colInd, vals, nnzPr,
                              mesh_parts, wls_, rdTags_, rdCounts_, &b_mesh_xs);
          } else if (mesh_xs.size(1) == 2) {
            b_assemble_body_par(mesh_xs, mesh_conn, stcls, mesh_n2cPtr,
                                mesh_n2cList, params->degree, params->wlsInterp0,
                                params->wlsUseDag, rowPtr, colInd, vals, nnzPr,
                                mesh_parts, wls_, rdTags_, rdCounts_, &b_mesh_xs);
          } else {
            c_assemble_body_par(mesh_xs, mesh_conn, stcls, mesh_n2cPtr,
                                mesh_n2cList, params->degree, params->wlsInterp0,
                                params->wlsUseDag, rowPtr, colInd, vals, nnzPr,
                                mesh_parts, wls_, rdTags_, rdCounts_, &b_mesh_xs);
          }

        } catch (const std::runtime_error& m2cExc) {
          m2cTryBlkErrCode = (1);
          m2cPrintError("runtime_error %s\n", m2cExc.what());

        }

        catch (const std::logic_error& m2cExc)
        {
          m2cTryBlkErrCode = (1);
          m2cPrintError("logic_error %s\n", m2cExc.what());

        }

        catch (const std::exception& m2cExc)
        {
          m2cTryBlkErrCode = (1);
          m2cPrintError("exception %s\n", m2cExc.what());

        }

        catch (...)
        {
          m2cTryBlkErrCode = (1);
          m2cPrintError("Unknown error detected from C++ exceptions\n");

        } //try

        if ((int)m2cTryBlkErrCode != 0) {
          throw std::runtime_error("omp4m:runtimeErrorInThread");
        }
      } else {
        int m2cTryBlkErrCode;

        //  task
        loop_ub = params->nThreads;

        m2cTryBlkErrCode = (0);
        if (params->nThreads <= 0) {
          loop_ub = 1;

#ifdef _OPENMP

          loop_ub = omp_get_max_threads();

#endif //_OPENMP

        }

#pragma omp parallel default(shared) num_threads(loop_ub)

        try { //try
          if (mesh_xs.size(1) == 3) {
            assemble_body_task(mesh_xs, mesh_conn, stcls, mesh_n2cPtr,
                               mesh_n2cList, params->degree, params->wlsInterp0,
                               params->wlsUseDag, rowPtr, colInd, vals, nnzPr,
                               mesh_parts, wls_, rdTags_, rdCounts_);
          } else if (mesh_xs.size(1) == 2) {
            b_assemble_body_task(mesh_xs, mesh_conn, stcls, mesh_n2cPtr,
                                 mesh_n2cList, params->degree,
                                 params->wlsInterp0, params->wlsUseDag, rowPtr,
                                 colInd, vals, nnzPr, mesh_parts, wls_, rdTags_,
                                 rdCounts_);
          } else {
            c_assemble_body_task(mesh_xs, mesh_conn, stcls, mesh_n2cPtr,
                                 mesh_n2cList, params->degree,
                                 params->wlsInterp0, params->wlsUseDag, rowPtr,
                                 colInd, vals, nnzPr, mesh_parts, wls_, rdTags_,
                                 rdCounts_);
          }

        } catch (const std::runtime_error& m2cExc) {
          m2cTryBlkErrCode = (1);
          m2cPrintError("runtime_error %s\n", m2cExc.what());

        }

        catch (const std::logic_error& m2cExc)
        {
          m2cTryBlkErrCode = (1);
          m2cPrintError("logic_error %s\n", m2cExc.what());

        }

        catch (const std::exception& m2cExc)
        {
          m2cTryBlkErrCode = (1);
          m2cPrintError("exception %s\n", m2cExc.what());

        }

        catch (...)
        {
          m2cTryBlkErrCode = (1);
          m2cPrintError("Unknown error detected from C++ exceptions\n");

        } //try

        if ((int)m2cTryBlkErrCode != 0) {
          throw std::runtime_error("omp4m:runtimeErrorInThread");
        }
      }
    }

    nrd = 0;
    i = rdCounts_.size(0);
    for (int b_i{0}; b_i < i; b_i++) {
      if (rdCounts_[b_i] < 0) {
        m2cErrMsgIdAndTxt("assemble_surf:badGEQP3",
                          "GEQP3 returned negative info at node %d in thread %d",
                          -rdCounts_[b_i], b_i);
      }

      nrd += rdCounts_[b_i];
    }

    rdNodes.set_size(nrd);
    if (nrd != 0) {
      int j;
      j = 0;
      i = rdTags_.size(0);
      for (int b_i{0}; b_i < i; b_i++) {
        if (rdTags_[b_i]) {
          rdNodes[j] = b_i + 1;
          j++;
        }
      }

      if (params->verbose > 1) {
        for (int b_i{0}; b_i < nrd; b_i++) {
          m2cPrintf(" Node %d is rank-deficient...\n", rdNodes[b_i]);
          fflush(stdout);
        }
      }
    }
  }

  static int assemble_body_kernel(const ::coder::array<double, 2U> &xs, const ::
    coder::array<int, 2U> &stcls, int degree, int iend, WlsDataStruct *wls)
  {
    double xsLocal[3072];
    int gIDs[1024];
    int b_iv[2];
    int rdCounts;
    boolean_T exitg1;
    rdCounts = 0;

    //  kernel
    exitg1 = false;
    while ((!exitg1) && (iend >= 1)) {
      double ori_idx_0;
      double ori_idx_1;
      double ori_idx_2;
      int nPoints;
      int xsLocal_tmp;

      //  get nodal ID
      nPoints = stcls[stcls.size(1) - 1] - 1;

      //  fetch coords, global IDs, and then localize coordinates
      for (int j{0}; j <= nPoints; j++) {
        gIDs[j] = stcls[j];
      }

      //  fetch coordinates
      for (int k{0}; k <= nPoints; k++) {
        xsLocal_tmp = gIDs[k];
        xsLocal[3 * k] = xs[xs.size(1) * (xsLocal_tmp - 1)];
        xsLocal[3 * k + 1] = xs[xs.size(1) * (xsLocal_tmp - 1) + 1];
        xsLocal[3 * k + 2] = xs[xs.size(1) * (xsLocal_tmp - 1) + 2];
      }

      //  origin
      ori_idx_0 = xsLocal[0];
      ori_idx_1 = xsLocal[1];
      ori_idx_2 = xsLocal[2];

      //  localize
      for (int k{0}; k <= nPoints; k++) {
        xsLocal[3 * k] -= ori_idx_0;
        xsLocal_tmp = 3 * k + 1;
        xsLocal[xsLocal_tmp] -= ori_idx_1;
        xsLocal_tmp = 3 * k + 2;
        xsLocal[xsLocal_tmp] -= ori_idx_2;
      }

      //  compute wls
      b_iv[0] = 1024;
      b_iv[1] = 3;
      wls_init(&wls->wlsObj, xsLocal, b_iv, wls->wlsWgts.name.data,
               wls->wlsWgts.params_pointwise.data,
               wls->wlsWgts.params_pointwise.size, degree, wls->interp0,
               wls->useDag, nPoints + 1);
      if (wls->wlsObj.rank < 0) {
        //  record error node and return
        exitg1 = true;
      }
    }

    return rdCounts;
  }

  static void assemble_body_par(const ::coder::array<double, 2U> &xs, const ::
    coder::array<int, 2U> &conn, const ::coder::array<int, 2U> &stcls, const ::
    coder::array<int, 1U> &n2cPtr, const ::coder::array<int, 1U> &n2cList, int
    degree, boolean_T interp0, boolean_T useDag, const ::coder::array<int, 1U>
    &rowPtr, const ::coder::array<int, 1U> &colInd, ::coder::array<double, 1U>
    &vals, const ::coder::array<int, 1U> &nnzPr, const ::coder::array<Omp4mPart,
    1U> &parts, ::coder::array<WlsDataStruct, 1U> &wls, ::coder::array<boolean_T,
    1U> &rdTags, ::coder::array<int, 1U> &rdCounts, Omp4mPartContext
    *partContext)
  {
    int lvl;
    int n;
    int partid;
    int rdCountsLoc;
    boolean_T exitg1;

    //  parallel kernel
    n = 0;

#ifdef _OPENMP

    n = omp_get_thread_num();

#endif //_OPENMP

    rdCounts[n] = 0;
    rdCountsLoc = 0;
    b_WlsDataStruct(degree, interp0, useDag, &wls[n].wlsObj, &wls[n].wlsWgts,
                    wls[n].coeffs, &wls[n].interp0, &wls[n].useDag);

#pragma omp single

    { //single
      partContext->atomic_counters.set_size(parts.size(0));
      partid = parts.size(0);
      for (int i{0}; i < partid; i++) {
        partContext->atomic_counters[i] = 1;
      }

    } //single

    //  safe and shared across all threads
    lvl = 0;
    exitg1 = false;
    while ((!exitg1) && (lvl <= parts.size(0) - 1)) {
      //  assemble for this part
      int exitg2;
      do {
        exitg2 = 0;

#pragma omp atomic capture

        {                              // atomic capture();
          //  to prevent Coder putting this line before pragma
          partid = (partContext->atomic_counters[lvl]);
          partContext->atomic_counters[lvl] = partContext->atomic_counters[lvl]
            + 1;
        }                              // atomic capture();

        //  to prevent Coder detecting v is for index and automatically subtract
        if ((lvl + 1 <= parts.size(0)) && (partid <= parts[lvl].nparts)) {
          d_assemble_body_kernel(xs, conn, stcls, n2cPtr, n2cList, degree,
            rowPtr, colInd, vals, nnzPr, parts[lvl].part_list, parts[lvl].
            part_ptr[partid - 1], parts[lvl].part_ptr[partid] - 1, &wls[n],
            rdTags, &rdCountsLoc);
          if (rdCountsLoc < 0) {
            exitg2 = 1;
          } else {
            //  LAPACK error
            rdCounts[n] = rdCounts[n] + rdCountsLoc;
          }
        } else {
          exitg2 = 1;
        }
      } while (exitg2 == 0);

      if (rdCountsLoc < 0) {
        rdCounts[n] = rdCountsLoc;

        //  error
        exitg1 = true;
      } else {

#pragma omp barrier

        //  barrier here
        if (parts[lvl].shared_ents.size(0) != 0) {

#pragma omp single nowait

          { //single
            d_assemble_body_kernel(xs, conn, stcls, n2cPtr, n2cList, degree,
              rowPtr, colInd, vals, nnzPr, parts[lvl].shared_ents, 1, parts[lvl]
              .shared_ents.size(0), &wls[n], rdTags, &rdCountsLoc);
            if (rdCountsLoc < 0) {
              rdCounts[n] = rdCountsLoc;

              //  error
            } else {
              rdCounts[n] = rdCounts[n] + rdCountsLoc;
            }

          } //single

          exitg1 = true;
        } else {
          lvl++;
        }
      }
    }
  }

  static void assemble_body_task(const ::coder::array<double, 2U> &xs, const ::
    coder::array<int, 2U> &conn, const ::coder::array<int, 2U> &stcls, const ::
    coder::array<int, 1U> &n2cPtr, const ::coder::array<int, 1U> &n2cList, int
    degree, boolean_T interp0, boolean_T useDag, const ::coder::array<int, 1U>
    &rowPtr, const ::coder::array<int, 1U> &colInd, ::coder::array<double, 1U>
    &vals, const ::coder::array<int, 1U> &nnzPr, const ::coder::array<Omp4mPart,
    1U> &parts, ::coder::array<WlsDataStruct, 1U> &wls, ::coder::array<boolean_T,
    1U> &rdTags, ::coder::array<int, 1U> &rdCounts)
  {
    int lvl;
    int n;
    boolean_T exitg1;

    //  kernel for running in task-based parallelism
    n = 0;

#ifdef _OPENMP

    n = omp_get_thread_num();

#endif //_OPENMP

    rdCounts[n] = 0;
    b_WlsDataStruct(degree, interp0, useDag, &wls[n].wlsObj, &wls[n].wlsWgts,
                    wls[n].coeffs, &wls[n].interp0, &wls[n].useDag);

#pragma omp barrier

    lvl = 0;
    exitg1 = false;
    while ((!exitg1) && (lvl <= parts.size(0) - 1)) {
      int i;

#pragma omp single

      { //single
        i = parts[lvl].nparts;
        for (int mypart{0}; mypart < i; mypart++) {

#pragma omp task default(shared)

          { //task
            assemble_body_task_kernel(xs, conn, stcls, n2cPtr, n2cList, degree,
              rowPtr, colInd, vals, nnzPr, parts, mypart + 1, lvl + 1, wls,
              rdTags, rdCounts);

          } //task
        }

      } //single

      //  implicit barrier
      if (parts[lvl].shared_ents.size(0) != 0) {
        int rdCountsLoc;

#pragma omp single nowait

        { //single
          rdCountsLoc = g_assemble_body_kernel(xs, conn, stcls, n2cPtr, n2cList,
            degree, rowPtr, colInd, vals, nnzPr, parts[lvl].shared_ents, 1,
            parts[lvl].shared_ents.size(0), &wls[n], rdTags);
          if (rdCountsLoc < 0) {
            rdCounts[n] = rdCountsLoc;

            //  error
          } else {
            rdCounts[n] = rdCounts[n] + rdCountsLoc;
          }

        } //single

        exitg1 = true;
      } else {
        lvl++;
      }
    }
  }

  static void assemble_body_task_kernel(const ::coder::array<double, 2U> &xs,
    const ::coder::array<int, 2U> &conn, const ::coder::array<int, 2U> &stcls,
    const ::coder::array<int, 1U> &n2cPtr, const ::coder::array<int, 1U>
    &n2cList, int degree, const ::coder::array<int, 1U> &rowPtr, const ::coder::
    array<int, 1U> &colInd, ::coder::array<double, 1U> &vals, const ::coder::
    array<int, 1U> &nnzPr, const ::coder::array<Omp4mPart, 1U> &parts, int
    mypart, int lvl, ::coder::array<WlsDataStruct, 1U> &wls, ::coder::array<
    boolean_T, 1U> &rdTags, ::coder::array<int, 1U> &rdCounts)
  {
    int n;

    //  Gets the thread number of the thread, within the team, making this call.
    n = 0;

#ifdef _OPENMP

    n = omp_get_thread_num();

#endif //_OPENMP

    //  if a previously assembled operator had a LAPACK error, then stop
    if (rdCounts[n] >= 0) {
      int rdCountsLoc;

      rdCountsLoc = g_assemble_body_kernel(xs, conn, stcls, n2cPtr, n2cList,
        degree, rowPtr, colInd, vals, nnzPr, parts[lvl - 1].part_list, parts[lvl
        - 1].part_ptr[mypart - 1], parts[lvl - 1].part_ptr[mypart] - 1, &wls[n],
        rdTags);
      if (rdCountsLoc < 0) {
        rdCounts[n] = rdCountsLoc;

        //  error
      } else {
        rdCounts[n] = rdCounts[n] + rdCountsLoc;
      }
    }
  }

  static void assemble_surf(const ::coder::array<double, 2U> &mesh_xs, const ::
    coder::array<int, 2U> &mesh_conn, const ::coder::array<double, 2U>
    &mesh_dirs, const ::coder::array<int, 1U> &mesh_n2cPtr, const ::coder::array<
    int, 1U> &mesh_n2cList, const ::coder::array<Omp4mPart, 1U> &mesh_parts,
    const ::coder::array<int, 2U> &stcls, const ::coder::array<int, 1U> &rowPtr,
    const ::coder::array<int, 1U> &colInd, ::coder::array<double, 1U> &vals,
    const ::coder::array<int, 1U> &nnzPr, ::coder::array<int, 1U> &rdNodes,
    const RdiParams *params)
  {
    ::coder::array<WlsDataStruct, 1U> wls_;
    ::coder::array<int, 1U> rdCounts_;
    ::coder::array<boolean_T, 1U> rdTags_;
    Omp4mPartContext b_mesh_xs;
    WlsDataStruct a;
    int i;
    int loop_ub;
    int nrd;

    //  Surface assembly
    rdTags_.set_size(mesh_xs.size(0));
    loop_ub = mesh_xs.size(0);
    for (i = 0; i < loop_ub; i++) {
      rdTags_[i] = false;
    }

    if (mesh_xs.size(1) == 3) {
      if (rdNodes.size(0) != 0) {
        b_WlsDataStruct(params->degree, params->wlsInterp0, params->wlsUseDag,
                        &a.wlsObj, &a.wlsWgts, a.coeffs, &a.interp0, &a.useDag);
        loop_ub = assemble_surf_kernel(mesh_xs, mesh_conn, stcls, mesh_n2cPtr,
          mesh_n2cList, params->surfType, params->degree, mesh_dirs, rowPtr,
          colInd, vals, nnzPr, rdNodes, rdNodes.size(0), &a, rdTags_);
        rdCounts_.set_size(1);
        rdCounts_[0] = loop_ub;
      } else {
        wls_.set_size(params->nThreads);
        rdCounts_.set_size(params->nThreads);
        if (params->parTask) {
          int m2cTryBlkErrCode;
          loop_ub = params->nThreads;

          m2cTryBlkErrCode = (0);
          if (params->nThreads <= 0) {
            loop_ub = 1;

#ifdef _OPENMP

            loop_ub = omp_get_max_threads();

#endif //_OPENMP

          }

#pragma omp parallel default(shared) num_threads(loop_ub)

          try { //try
            assemble_surf_task(mesh_xs, mesh_conn, stcls, mesh_n2cPtr,
                               mesh_n2cList, params->surfType, params->degree,
                               params->wlsInterp0, params->wlsUseDag, mesh_dirs,
                               rowPtr, colInd, vals, nnzPr, mesh_parts, wls_,
                               rdTags_, rdCounts_);

          } catch (const std::runtime_error& m2cExc) {
            m2cTryBlkErrCode = (1);
            m2cPrintError("runtime_error %s\n", m2cExc.what());

          }

          catch (const std::logic_error& m2cExc)
          {
            m2cTryBlkErrCode = (1);
            m2cPrintError("logic_error %s\n", m2cExc.what());

          }

          catch (const std::exception& m2cExc)
          {
            m2cTryBlkErrCode = (1);
            m2cPrintError("exception %s\n", m2cExc.what());

          }

          catch (...)
          {
            m2cTryBlkErrCode = (1);
            m2cPrintError("Unknown error detected from C++ exceptions\n");

          } //try

          if ((int)m2cTryBlkErrCode != 0) {
            throw std::runtime_error("omp4m:runtimeErrorInThread");
          }
        } else {
          int m2cTryBlkErrCode;
          loop_ub = params->nThreads;

          m2cTryBlkErrCode = (0);
          if (params->nThreads <= 0) {
            loop_ub = 1;

#ifdef _OPENMP

            loop_ub = omp_get_max_threads();

#endif //_OPENMP

          }

#pragma omp parallel default(shared) num_threads(loop_ub)

          try { //try
            assemble_surf_par(mesh_xs, mesh_conn, stcls, mesh_n2cPtr,
                              mesh_n2cList, params->surfType, params->degree,
                              params->wlsInterp0, params->wlsUseDag, mesh_dirs,
                              rowPtr, colInd, vals, nnzPr, mesh_parts, wls_,
                              rdTags_, rdCounts_, &b_mesh_xs);

          } catch (const std::runtime_error& m2cExc) {
            m2cTryBlkErrCode = (1);
            m2cPrintError("runtime_error %s\n", m2cExc.what());

          }

          catch (const std::logic_error& m2cExc)
          {
            m2cTryBlkErrCode = (1);
            m2cPrintError("logic_error %s\n", m2cExc.what());

          }

          catch (const std::exception& m2cExc)
          {
            m2cTryBlkErrCode = (1);
            m2cPrintError("exception %s\n", m2cExc.what());

          }

          catch (...)
          {
            m2cTryBlkErrCode = (1);
            m2cPrintError("Unknown error detected from C++ exceptions\n");

          } //try

          if ((int)m2cTryBlkErrCode != 0) {
            throw std::runtime_error("omp4m:runtimeErrorInThread");
          }
        }
      }
    } else if (rdNodes.size(0) != 0) {
      b_WlsDataStruct(params->degree, params->wlsInterp0, params->wlsUseDag,
                      &a.wlsObj, &a.wlsWgts, a.coeffs, &a.interp0, &a.useDag);
      loop_ub = c_assemble_surf_kernel(mesh_xs, mesh_conn, stcls, mesh_n2cPtr,
        mesh_n2cList, params->surfType, params->degree, mesh_dirs, rowPtr,
        colInd, vals, nnzPr, rdNodes, rdNodes.size(0), &a, rdTags_);
      rdCounts_.set_size(1);
      rdCounts_[0] = loop_ub;
    } else {
      wls_.set_size(params->nThreads);
      rdCounts_.set_size(params->nThreads);
      if (params->parTask) {
        int m2cTryBlkErrCode;
        loop_ub = params->nThreads;

        m2cTryBlkErrCode = (0);
        if (params->nThreads <= 0) {
          loop_ub = 1;

#ifdef _OPENMP

          loop_ub = omp_get_max_threads();

#endif //_OPENMP

        }

#pragma omp parallel default(shared) num_threads(loop_ub)

        try { //try
          b_assemble_surf_task(mesh_xs, mesh_conn, stcls, mesh_n2cPtr,
                               mesh_n2cList, params->surfType, params->degree,
                               params->wlsInterp0, params->wlsUseDag, mesh_dirs,
                               rowPtr, colInd, vals, nnzPr, mesh_parts, wls_,
                               rdTags_, rdCounts_);

        } catch (const std::runtime_error& m2cExc) {
          m2cTryBlkErrCode = (1);
          m2cPrintError("runtime_error %s\n", m2cExc.what());

        }

        catch (const std::logic_error& m2cExc)
        {
          m2cTryBlkErrCode = (1);
          m2cPrintError("logic_error %s\n", m2cExc.what());

        }

        catch (const std::exception& m2cExc)
        {
          m2cTryBlkErrCode = (1);
          m2cPrintError("exception %s\n", m2cExc.what());

        }

        catch (...)
        {
          m2cTryBlkErrCode = (1);
          m2cPrintError("Unknown error detected from C++ exceptions\n");

        } //try

        if ((int)m2cTryBlkErrCode != 0) {
          throw std::runtime_error("omp4m:runtimeErrorInThread");
        }
      } else {
        int m2cTryBlkErrCode;
        loop_ub = params->nThreads;

        m2cTryBlkErrCode = (0);
        if (params->nThreads <= 0) {
          loop_ub = 1;

#ifdef _OPENMP

          loop_ub = omp_get_max_threads();

#endif //_OPENMP

        }

#pragma omp parallel default(shared) num_threads(loop_ub)

        try { //try
          b_assemble_surf_par(mesh_xs, mesh_conn, stcls, mesh_n2cPtr,
                              mesh_n2cList, params->surfType, params->degree,
                              params->wlsInterp0, params->wlsUseDag, mesh_dirs,
                              rowPtr, colInd, vals, nnzPr, mesh_parts, wls_,
                              rdTags_, rdCounts_, &b_mesh_xs);

        } catch (const std::runtime_error& m2cExc) {
          m2cTryBlkErrCode = (1);
          m2cPrintError("runtime_error %s\n", m2cExc.what());

        }

        catch (const std::logic_error& m2cExc)
        {
          m2cTryBlkErrCode = (1);
          m2cPrintError("logic_error %s\n", m2cExc.what());

        }

        catch (const std::exception& m2cExc)
        {
          m2cTryBlkErrCode = (1);
          m2cPrintError("exception %s\n", m2cExc.what());

        }

        catch (...)
        {
          m2cTryBlkErrCode = (1);
          m2cPrintError("Unknown error detected from C++ exceptions\n");

        } //try

        if ((int)m2cTryBlkErrCode != 0) {
          throw std::runtime_error("omp4m:runtimeErrorInThread");
        }
      }
    }

    nrd = 0;
    i = rdCounts_.size(0);
    for (int b_i{0}; b_i < i; b_i++) {
      if (rdCounts_[b_i] < 0) {
        m2cErrMsgIdAndTxt("assemble_surf:badGEQP3",
                          "GEQP3 returned negative info at node %d in thread %d",
                          -rdCounts_[b_i], b_i);
      }

      nrd += rdCounts_[b_i];
    }

    rdNodes.set_size(nrd);
    if (nrd != 0) {
      int j;
      j = 0;
      i = rdTags_.size(0);
      for (int b_i{0}; b_i < i; b_i++) {
        if (rdTags_[b_i]) {
          rdNodes[j] = b_i + 1;
          j++;
        }
      }

      if (params->verbose > 1) {
        for (int b_i{0}; b_i < nrd; b_i++) {
          m2cPrintf(" Node %d is rank-deficient...\n", rdNodes[b_i]);
          fflush(stdout);
        }
      }
    }
  }

  static int assemble_surf_kernel(const ::coder::array<double, 2U> &xs, const ::
    coder::array<int, 2U> &conn, const ::coder::array<int, 2U> &stcls, const ::
    coder::array<int, 1U> &n2cPtr, const ::coder::array<int, 1U> &n2cList, int
    surfType, int degree, const ::coder::array<double, 2U> &nrms, const ::coder::
    array<int, 1U> &rowPtr, const ::coder::array<int, 1U> &colInd, ::coder::
    array<double, 1U> &vals, const ::coder::array<int, 1U> &nnzPr, const ::coder::
    array<int, 1U> &nrange, int iend, WlsDataStruct *wls, ::coder::array<
    boolean_T, 1U> &rdTags)
  {
    double xsLocal[1536];
    double us[1024];
    double Ns[512];
    double nrm0_idx_0;
    double nrm0_idx_1;
    double nrm0_idx_2;
    int eids[512];
    int gIDs[512];
    int b_iv[2];
    int iv1[2];
    int i;
    int rdCounts;
    boolean_T exitg1;
    boolean_T isSphSurf;
    rdCounts = 0;

    //  kernel for surface computation
    wls->wlsWgts.params_pointwise.size[1] = 1;
    wls->wlsWgts.params_pointwise.size[0] = 512;

    //  check if spherical surface
    isSphSurf = (surfType == 1);

    //  loop begins
    i = 0;
    exitg1 = false;
    while ((!exitg1) && (i + 1 <= iend)) {
      double a;
      double d;
      double d1;
      double d2;
      double d3;
      double ori_idx_0;
      double ori_idx_1;
      double ori_idx_2;
      double t_idx_0;
      double t_idx_1;
      double t_idx_2;
      double t_idx_3;
      double t_idx_4;
      double t_idx_5;
      int k;
      int nPoints;
      int us_tmp;
      int xsLocal_tmp;

      //  get nodal ID
      nPoints = stcls[(stcls.size(1) + stcls.size(1) * i) - 1] - 1;

      //  fetch data to local buffer and localize the coordinates
      for (int j{0}; j <= nPoints; j++) {
        gIDs[j] = stcls[j + stcls.size(1) * i];
      }

      //  fetch coordinates
      for (k = 0; k <= nPoints; k++) {
        us_tmp = gIDs[k];
        xsLocal[3 * k] = xs[xs.size(1) * (us_tmp - 1)];
        xsLocal[3 * k + 1] = xs[xs.size(1) * (us_tmp - 1) + 1];
        xsLocal[3 * k + 2] = xs[xs.size(1) * (us_tmp - 1) + 2];
      }

      //  origin
      ori_idx_0 = xsLocal[0];
      ori_idx_1 = xsLocal[1];
      ori_idx_2 = xsLocal[2];

      //  localize
      for (k = 0; k <= nPoints; k++) {
        xsLocal[3 * k] -= ori_idx_0;
        xsLocal_tmp = 3 * k + 1;
        xsLocal[xsLocal_tmp] -= ori_idx_1;
        xsLocal_tmp = 3 * k + 2;
        xsLocal[xsLocal_tmp] -= ori_idx_2;
      }

      //  get normal direction
      if ((nrms.size(0) != 0) && (nrms.size(1) != 0)) {
        nrm0_idx_0 = nrms[nrms.size(1) * (nrange[i] - 1)];
        nrm0_idx_1 = nrms[nrms.size(1) * (nrange[i] - 1) + 1];
        nrm0_idx_2 = nrms[nrms.size(1) * (nrange[i] - 1) + 2];
      } else if (isSphSurf) {
        //  spherical or circle
        nrm0_idx_0 = xs[xs.size(1) * (nrange[i] - 1)];
        nrm0_idx_1 = xs[xs.size(1) * (nrange[i] - 1) + 1];
        nrm0_idx_2 = xs[xs.size(1) * (nrange[i] - 1) + 2];
      }

      //  compute tagent
      d = std::abs(nrm0_idx_0);
      if ((d > std::abs(nrm0_idx_1)) && (d > std::abs(nrm0_idx_2))) {
        t_idx_0 = -nrm0_idx_0 * nrm0_idx_1;
        t_idx_2 = 1.0 - nrm0_idx_1 * nrm0_idx_1;
        t_idx_4 = -nrm0_idx_1 * nrm0_idx_2;
      } else {
        t_idx_0 = 1.0 - nrm0_idx_0 * nrm0_idx_0;
        t_idx_2 = -nrm0_idx_0 * nrm0_idx_1;
        t_idx_4 = -nrm0_idx_0 * nrm0_idx_2;
      }

      a = std::sqrt((t_idx_0 * t_idx_0 + t_idx_2 * t_idx_2) + t_idx_4 * t_idx_4);
      t_idx_0 /= a;
      t_idx_2 /= a;
      t_idx_4 /= a;

      //  cross
      t_idx_1 = nrm0_idx_1 * t_idx_4 - nrm0_idx_2 * t_idx_2;
      t_idx_3 = t_idx_0 * nrm0_idx_2 - nrm0_idx_0 * t_idx_4;
      t_idx_5 = nrm0_idx_0 * t_idx_2 - t_idx_0 * nrm0_idx_1;

      //  project onto the tangent plan
      for (k = 0; k <= nPoints; k++) {
        us_tmp = k << 1;
        us[us_tmp] = 0.0;
        us[us_tmp + 1] = 0.0;
      }

      //  matrix-matrix, using mem-efficient loop
      for (int ii{0}; ii <= nPoints; ii++) {
        d = xsLocal[3 * ii];
        d1 = xsLocal[3 * ii + 1];
        us_tmp = ii << 1;
        d2 = (us[us_tmp + 1] + d * t_idx_1) + d1 * t_idx_3;
        d3 = xsLocal[3 * ii + 2];
        us[us_tmp] = ((us[us_tmp] + d * t_idx_0) + t_idx_2 * d1) + t_idx_4 * d3;
        d2 += t_idx_5 * d3;
        us[us_tmp + 1] = d2;
      }

      //  compute scaling to the WLS weights by inner product of normals
      if ((nrms.size(0) != 0) && (nrms.size(1) != 0)) {
        for (int j{2}; j <= nPoints + 1; j++) {
          us_tmp = stcls[(j + stcls.size(1) * i) - 1] - 1;
          a = (nrm0_idx_0 * nrms[nrms.size(1) * us_tmp] + nrm0_idx_1 *
               nrms[nrms.size(1) * us_tmp + 1]) + nrm0_idx_2 * nrms[nrms.size(1)
            * us_tmp + 2];
          if (a < 0.0) {
            a = 0.0;
          }

          wls->wlsWgts.params_pointwise.data[j - 1] = a;
        }
      } else if (isSphSurf) {
        //  spherical or circle
        for (int j{2}; j <= nPoints + 1; j++) {
          int idx;
          idx = stcls[(j + stcls.size(1) * i) - 1] - 1;
          a = (nrm0_idx_0 * xs[xs.size(1) * idx] + nrm0_idx_1 * xs[xs.size(1) *
               idx + 1]) + nrm0_idx_2 * xs[xs.size(1) * idx + 2];
          if (a < 0.0) {
            a = 0.0;
          }

          wls->wlsWgts.params_pointwise.data[j - 1] = a;
        }
      }

      wls->wlsWgts.params_pointwise.data[0] = 1.0;

      //  compute wls
      b_iv[0] = 512;
      b_iv[1] = 2;
      wls_init(&wls->wlsObj, us, b_iv, wls->wlsWgts.name.data,
               wls->wlsWgts.params_pointwise.data,
               wls->wlsWgts.params_pointwise.size, degree, wls->interp0,
               wls->useDag, nPoints + 1);
      if (wls->wlsObj.rank < 0) {
        rdCounts = -nrange[i];

        //  record error node and return
        exitg1 = true;
      } else {
        if (wls->wlsObj.rank != wls->wlsObj.ncols) {
          rdCounts++;

          //  increment counter for RD nodes
          rdTags[nrange[i] - 1] = true;
        } else {
          int b_xsLocal_tmp;
          int nCells;

          //  get centers and shape function Ns at centers
          us_tmp = n2cPtr[nrange[i] - 1];
          nCells = (n2cPtr[nrange[i]] - us_tmp) - 1;
          for (int j{0}; j <= nCells; j++) {
            int eid;
            int npts;
            eid = n2cList[(us_tmp + j) - 1] - 1;

            eids[j] = eid + 1;
            for (npts = conn.size(1); conn[(npts + conn.size(1) * eid) - 1] <= 0;
                 npts--) {
            }

            Ns[j] = 1.0 / static_cast<double>(npts);

            //  compute center
            xsLocal[3 * j] = xs[xs.size(1) * (conn[conn.size(1) * eid] - 1)];
            xsLocal_tmp = 3 * j + 1;
            xsLocal[xsLocal_tmp] = xs[xs.size(1) * (conn[conn.size(1) * eid] - 1)
              + 1];
            b_xsLocal_tmp = 3 * j + 2;
            xsLocal[b_xsLocal_tmp] = xs[xs.size(1) * (conn[conn.size(1) * eid] -
              1) + 2];

            //  C++
            for (int ii{2}; ii <= npts; ii++) {
              int c_xsLocal_tmp;
              c_xsLocal_tmp = conn[(ii + conn.size(1) * eid) - 1] - 1;
              xsLocal[3 * j] += xs[xs.size(1) * c_xsLocal_tmp];
              xsLocal[xsLocal_tmp] += xs[xs.size(1) * c_xsLocal_tmp + 1];
              xsLocal[b_xsLocal_tmp] += xs[xs.size(1) * c_xsLocal_tmp + 2];
            }

            d = Ns[j];
            d1 = xsLocal[3 * j] * d;
            d2 = xsLocal[xsLocal_tmp] * d;
            d *= xsLocal[b_xsLocal_tmp];
            if (isSphSurf) {
              a = std::sqrt((d1 * d1 + d2 * d2) + d * d);
              d1 /= a;
              d2 /= a;
              d /= a;
            }

            //  localize
            d1 -= ori_idx_0;
            xsLocal[3 * j] = d1;
            d2 -= ori_idx_1;
            xsLocal[xsLocal_tmp] = d2;
            d -= ori_idx_2;
            xsLocal[b_xsLocal_tmp] = d;
          }

          for (k = 0; k <= nCells; k++) {
            us_tmp = k << 1;
            us[us_tmp] = 0.0;
            us[us_tmp + 1] = 0.0;
          }

          //  matrix-matrix, using mem-efficient loop
          for (int ii{0}; ii <= nCells; ii++) {
            d = xsLocal[3 * ii];
            d1 = xsLocal[3 * ii + 1];
            us_tmp = ii << 1;
            d2 = (us[us_tmp + 1] + d * t_idx_1) + d1 * t_idx_3;
            d3 = xsLocal[3 * ii + 2];
            us[us_tmp] = ((us[us_tmp] + d * t_idx_0) + t_idx_2 * d1) + t_idx_4 *
              d3;
            d2 += t_idx_5 * d3;
            us[us_tmp + 1] = d2;
          }

          iv1[0] = 512;
          iv1[1] = 2;
          wls_func(&wls->wlsObj, us, iv1, nCells + 1, wls->coeffs);

          //  Compute local OSUS coefficients
          for (k = 0; k <= nPoints; k++) {
            for (int j{0}; j <= nCells; j++) {
              wls->coeffs[j + wls->coeffs.size(1) * k] = wls->coeffs[j +
                wls->coeffs.size(1) * k] * Ns[j];
            }
          }

          //  substract linear interp
          for (int j{0}; j <= nCells; j++) {
            wls->coeffs[j] = wls->coeffs[j] - Ns[j];
          }

          //  kernel for updating crs matrix
          for (int j{0}; j <= nCells; j++) {
            int rowIdx;
            rowIdx = eids[j] - 1;
            for (int ii{0}; ii <= nPoints; ii++) {
              int b_i;
              b_xsLocal_tmp = gIDs[ii];
              b_i = rowPtr[rowIdx] - 1;
              us_tmp = (rowPtr[rowIdx] + nnzPr[eids[j] - 1]) - 1;
              int exitg2;
              do {
                exitg2 = 0;
                if (b_i + 1 <= us_tmp) {
                  if (colInd[b_i] == b_xsLocal_tmp) {
                    k = b_i;
                    exitg2 = 1;
                  } else {
                    b_i++;
                  }
                } else {
                  k = -1;

                  m2cErrMsgIdAndTxt("add_crs:missingIndex",
                                    "could not find column index %d",
                                    b_xsLocal_tmp);
                  exitg2 = 1;
                }
              } while (exitg2 == 0);

              vals[k] = vals[k] + wls->coeffs[j + wls->coeffs.size(1) * ii];
            }
          }
        }

        i++;
      }
    }

    return rdCounts;
  }

  static void assemble_surf_par(const ::coder::array<double, 2U> &xs, const ::
    coder::array<int, 2U> &conn, const ::coder::array<int, 2U> &stcls, const ::
    coder::array<int, 1U> &n2cPtr, const ::coder::array<int, 1U> &n2cList, int
    surfType, int degree, boolean_T interp0, boolean_T useDag, const ::coder::
    array<double, 2U> &nrms, const ::coder::array<int, 1U> &rowPtr, const ::
    coder::array<int, 1U> &colInd, ::coder::array<double, 1U> &vals, const ::
    coder::array<int, 1U> &nnzPr, const ::coder::array<Omp4mPart, 1U> &parts, ::
    coder::array<WlsDataStruct, 1U> &wls, ::coder::array<boolean_T, 1U> &rdTags,
    ::coder::array<int, 1U> &rdCounts, Omp4mPartContext *partContext)
  {
    int lvl;
    int n;
    int partid;
    int rdCountsLoc;
    boolean_T exitg1;

    //  kernel for parallel
    n = 0;

#ifdef _OPENMP

    n = omp_get_thread_num();

#endif //_OPENMP

    rdCounts[n] = 0;
    rdCountsLoc = 0;
    b_WlsDataStruct(degree, interp0, useDag, &wls[n].wlsObj, &wls[n].wlsWgts,
                    wls[n].coeffs, &wls[n].interp0, &wls[n].useDag);

#pragma omp single

    { //single
      partContext->atomic_counters.set_size(parts.size(0));
      partid = parts.size(0);
      for (int i{0}; i < partid; i++) {
        partContext->atomic_counters[i] = 1;
      }

    } //single

    //  safe and shared across threads
    lvl = 0;
    exitg1 = false;
    while ((!exitg1) && (lvl <= parts.size(0) - 1)) {
      //  Assemble for this part
      int exitg2;
      do {
        exitg2 = 0;

#pragma omp atomic capture

        {                              // atomic capture();
          //  to prevent Coder putting this line before pragma
          partid = (partContext->atomic_counters[lvl]);
          partContext->atomic_counters[lvl] = partContext->atomic_counters[lvl]
            + 1;
        }                              // atomic capture();

        //  to prevent Coder detecting v is for index and automatically subtract
        if ((lvl + 1 <= parts.size(0)) && (partid <= parts[lvl].nparts)) {
          b_assemble_surf_kernel(xs, conn, stcls, n2cPtr, n2cList, surfType,
            degree, nrms, rowPtr, colInd, vals, nnzPr, parts[lvl].part_list,
            parts[lvl].part_ptr[partid - 1], parts[lvl].part_ptr[partid] - 1,
            &wls[n], rdTags, &rdCountsLoc);
          if (rdCountsLoc < 0) {
            exitg2 = 1;
          } else {
            rdCounts[n] = rdCounts[n] + rdCountsLoc;
          }
        } else {
          exitg2 = 1;
        }
      } while (exitg2 == 0);

      if (rdCountsLoc < 0) {
        rdCounts[n] = rdCountsLoc;

        //  error
        exitg1 = true;
      } else {

#pragma omp barrier

        //  barrier here to avoid reset token before all break
        if (parts[lvl].shared_ents.size(0) != 0) {

#pragma omp single nowait

          { //single
            b_assemble_surf_kernel(xs, conn, stcls, n2cPtr, n2cList, surfType,
              degree, nrms, rowPtr, colInd, vals, nnzPr, parts[lvl].shared_ents,
              1, parts[lvl].shared_ents.size(0), &wls[n], rdTags, &rdCountsLoc);
            if (rdCountsLoc < 0) {
              rdCounts[n] = rdCountsLoc;

              //  error
            } else {
              rdCounts[n] = rdCounts[n] + rdCountsLoc;
            }

          } //single

          exitg1 = true;
        } else {
          lvl++;
        }
      }
    }
  }

  static void assemble_surf_task(const ::coder::array<double, 2U> &xs, const ::
    coder::array<int, 2U> &conn, const ::coder::array<int, 2U> &stcls, const ::
    coder::array<int, 1U> &n2cPtr, const ::coder::array<int, 1U> &n2cList, int
    surfType, int degree, boolean_T interp0, boolean_T useDag, const ::coder::
    array<double, 2U> &nrms, const ::coder::array<int, 1U> &rowPtr, const ::
    coder::array<int, 1U> &colInd, ::coder::array<double, 1U> &vals, const ::
    coder::array<int, 1U> &nnzPr, const ::coder::array<Omp4mPart, 1U> &parts, ::
    coder::array<WlsDataStruct, 1U> &wls, ::coder::array<boolean_T, 1U> &rdTags,
    ::coder::array<int, 1U> &rdCounts)
  {
    int lvl;
    int n;
    int rdCountsLoc;
    boolean_T exitg1;

    n = 0;

#ifdef _OPENMP

    n = omp_get_thread_num();

#endif //_OPENMP

    rdCounts[n] = 0;
    b_WlsDataStruct(degree, interp0, useDag, &wls[n].wlsObj, &wls[n].wlsWgts,
                    wls[n].coeffs, &wls[n].interp0, &wls[n].useDag);

#pragma omp barrier

    lvl = 0;
    exitg1 = false;
    while ((!exitg1) && (lvl <= parts.size(0) - 1)) {
      int i;

#pragma omp single

      { //single
        i = parts[lvl].nparts;
        for (int mypart{0}; mypart < i; mypart++) {

#pragma omp task default(shared)

          { //task
            assemble_surf_task_kernel(xs, conn, stcls, n2cPtr, n2cList, surfType,
              degree, nrms, rowPtr, colInd, vals, nnzPr, parts, mypart + 1, lvl
              + 1, wls, rdTags, rdCounts);

          } //task
        }

      } //single

      //  implicit barrier
      if (parts[lvl].shared_ents.size(0) != 0) {

#pragma omp single nowait

        { //single
          rdCountsLoc = 0;
          b_assemble_surf_kernel(xs, conn, stcls, n2cPtr, n2cList, surfType,
            degree, nrms, rowPtr, colInd, vals, nnzPr, parts[lvl].shared_ents, 1,
            parts[lvl].shared_ents.size(0), &wls[n], rdTags, &rdCountsLoc);
          if (rdCountsLoc < 0) {
            rdCounts[n] = rdCountsLoc;

            //  error
          } else {
            rdCounts[n] = rdCounts[n] + rdCountsLoc;
          }

        } //single

        exitg1 = true;
      } else {
        lvl++;
      }
    }
  }

  static void assemble_surf_task_kernel(const ::coder::array<double, 2U> &xs,
    const ::coder::array<int, 2U> &conn, const ::coder::array<int, 2U> &stcls,
    const ::coder::array<int, 1U> &n2cPtr, const ::coder::array<int, 1U>
    &n2cList, int surfType, int degree, const ::coder::array<double, 2U> &nrms,
    const ::coder::array<int, 1U> &rowPtr, const ::coder::array<int, 1U> &colInd,
    ::coder::array<double, 1U> &vals, const ::coder::array<int, 1U> &nnzPr,
    const ::coder::array<Omp4mPart, 1U> &parts, int mypart, int lvl, ::coder::
    array<WlsDataStruct, 1U> &wls, ::coder::array<boolean_T, 1U> &rdTags, ::
    coder::array<int, 1U> &rdCounts)
  {
    int n;
    int rdCountsLoc;

    //  Gets the thread number of the thread, within the team, making this call.
    n = 0;

#ifdef _OPENMP

    n = omp_get_thread_num();

#endif //_OPENMP

    //  if a previously assembled operator had a LAPACK error, then stop
    if (rdCounts[n] >= 0) {
      rdCountsLoc = 0;
      b_assemble_surf_kernel(xs, conn, stcls, n2cPtr, n2cList, surfType, degree,
        nrms, rowPtr, colInd, vals, nnzPr, parts[lvl - 1].part_list, parts[lvl -
        1].part_ptr[mypart - 1], parts[lvl - 1].part_ptr[mypart] - 1, &wls[n],
        rdTags, &rdCountsLoc);
      if (rdCountsLoc < 0) {
        rdCounts[n] = rdCountsLoc;

        //  error
      } else {
        rdCounts[n] = rdCounts[n] + rdCountsLoc;
      }
    }
  }

  static void b_WlsDataStruct(int degree, boolean_T interp0, boolean_T useDag,
    WlsObject *wls_wlsObj, WlsWeight *wls_wlsWgts, ::coder::array<double, 2U>
    &wls_coeffs, boolean_T *wls_interp0, boolean_T *wls_useDag)
  {
    static const char name[7]{ 'B', 'u', 'h', 'm', 'a', 'n', 'n' };

    //  Create a buffer struct for WLS
    wls_coeffs.set_size(0, 1);
    wls_wlsWgts->omit_rows.size[0] = 0;
    wls_wlsWgts->params_pointwise.size[1] = 0;
    wls_wlsWgts->params_pointwise.size[0] = 0;
    wls_wlsWgts->params_shared.size[0] = 0;
    wls_wlsWgts->name.size[1] = 7;
    wls_wlsWgts->name.size[0] = 1;
    for (int i{0}; i < 7; i++) {
      wls_wlsWgts->name.data[i] = name[i];
    }

    wls_wlsObj->QRt.size[1] = 0;
    wls_wlsObj->QRt.size[0] = 0;
    wls_wlsObj->rowmajor = true;
    wls_wlsObj->work.set_size(0);
    wls_wlsObj->jpvt.set_size(0);
    wls_wlsObj->rank = 0;
    wls_wlsObj->ncols = 0;
    wls_wlsObj->nrows = 0;
    wls_wlsObj->vdops.set_size(0, 0);
    wls_wlsObj->QR.set_size(0, 0);
    wls_wlsObj->V.set_size(0, 0);
    wls_wlsObj->dag.set_size(0, 0);
    wls_wlsObj->hs_inv.size[1] = 0;
    wls_wlsObj->hs_inv.size[0] = 1;
    wls_wlsObj->rweights.set_size(0);
    wls_wlsObj->origin.size[1] = 3;
    wls_wlsObj->origin.size[0] = 1;
    wls_wlsObj->origin.data[0] = 0.0;
    wls_wlsObj->origin.data[1] = 0.0;
    wls_wlsObj->origin.data[2] = 0.0;
    wls_wlsObj->us.set_size(0, 3);
    wls_wlsObj->stride = 0;
    wls_wlsObj->use_dag = true;
    wls_wlsObj->interp0 = 0;
    wls_wlsObj->order = 0;
    wls_wlsObj->degree = degree;
    wls_wlsObj->npoints = 0;
    *wls_interp0 = interp0;
    *wls_useDag = useDag;
  }

  static int b_assemble_body_kernel(const ::coder::array<double, 2U> &xs, const ::
    coder::array<int, 2U> &stcls, int degree, int iend, WlsDataStruct *wls)
  {
    double xsLocal[1024];
    int gIDs[512];
    int b_iv[2];
    int rdCounts;
    boolean_T exitg1;
    rdCounts = 0;

    //  kernel
    exitg1 = false;
    while ((!exitg1) && (iend >= 1)) {
      double ori_idx_0;
      double ori_idx_1;
      int nPoints;
      int xsLocal_tmp;

      //  get nodal ID
      nPoints = stcls[stcls.size(1) - 1] - 1;

      //  fetch coords, global IDs, and then localize coordinates
      for (int j{0}; j <= nPoints; j++) {
        gIDs[j] = stcls[j];
      }

      //  fetch coordinates
      for (int k{0}; k <= nPoints; k++) {
        int i;
        xsLocal_tmp = k << 1;
        i = gIDs[k];
        xsLocal[xsLocal_tmp] = xs[xs.size(1) * (i - 1)];
        xsLocal[xsLocal_tmp + 1] = xs[xs.size(1) * (i - 1) + 1];
      }

      //  origin
      ori_idx_0 = xsLocal[0];
      ori_idx_1 = xsLocal[1];

      //  localize
      for (int k{0}; k <= nPoints; k++) {
        xsLocal_tmp = k << 1;
        xsLocal[xsLocal_tmp] -= ori_idx_0;
        xsLocal[xsLocal_tmp + 1] -= ori_idx_1;
      }

      //  compute wls
      b_iv[0] = 512;
      b_iv[1] = 2;
      wls_init(&wls->wlsObj, xsLocal, b_iv, wls->wlsWgts.name.data,
               wls->wlsWgts.params_pointwise.data,
               wls->wlsWgts.params_pointwise.size, degree, wls->interp0,
               wls->useDag, nPoints + 1);
      if (wls->wlsObj.rank < 0) {
        //  record error node and return
        exitg1 = true;
      }
    }

    return rdCounts;
  }

  static void b_assemble_body_par(const ::coder::array<double, 2U> &xs, const ::
    coder::array<int, 2U> &conn, const ::coder::array<int, 2U> &stcls, const ::
    coder::array<int, 1U> &n2cPtr, const ::coder::array<int, 1U> &n2cList, int
    degree, boolean_T interp0, boolean_T useDag, const ::coder::array<int, 1U>
    &rowPtr, const ::coder::array<int, 1U> &colInd, ::coder::array<double, 1U>
    &vals, const ::coder::array<int, 1U> &nnzPr, const ::coder::array<Omp4mPart,
    1U> &parts, ::coder::array<WlsDataStruct, 1U> &wls, ::coder::array<boolean_T,
    1U> &rdTags, ::coder::array<int, 1U> &rdCounts, Omp4mPartContext
    *partContext)
  {
    int lvl;
    int n;
    int partid;
    int rdCountsLoc;
    boolean_T exitg1;

    //  parallel kernel
    n = 0;

#ifdef _OPENMP

    n = omp_get_thread_num();

#endif //_OPENMP

    rdCounts[n] = 0;
    rdCountsLoc = 0;
    b_WlsDataStruct(degree, interp0, useDag, &wls[n].wlsObj, &wls[n].wlsWgts,
                    wls[n].coeffs, &wls[n].interp0, &wls[n].useDag);

#pragma omp single

    { //single
      partContext->atomic_counters.set_size(parts.size(0));
      partid = parts.size(0);
      for (int i{0}; i < partid; i++) {
        partContext->atomic_counters[i] = 1;
      }

    } //single

    //  safe and shared across all threads
    lvl = 0;
    exitg1 = false;
    while ((!exitg1) && (lvl <= parts.size(0) - 1)) {
      //  assemble for this part
      int exitg2;
      do {
        exitg2 = 0;

#pragma omp atomic capture

        {                              // atomic capture();
          //  to prevent Coder putting this line before pragma
          partid = (partContext->atomic_counters[lvl]);
          partContext->atomic_counters[lvl] = partContext->atomic_counters[lvl]
            + 1;
        }                              // atomic capture();

        //  to prevent Coder detecting v is for index and automatically subtract
        if ((lvl + 1 <= parts.size(0)) && (partid <= parts[lvl].nparts)) {
          e_assemble_body_kernel(xs, conn, stcls, n2cPtr, n2cList, degree,
            rowPtr, colInd, vals, nnzPr, parts[lvl].part_list, parts[lvl].
            part_ptr[partid - 1], parts[lvl].part_ptr[partid] - 1, &wls[n],
            rdTags, &rdCountsLoc);
          if (rdCountsLoc < 0) {
            exitg2 = 1;
          } else {
            //  LAPACK error
            rdCounts[n] = rdCounts[n] + rdCountsLoc;
          }
        } else {
          exitg2 = 1;
        }
      } while (exitg2 == 0);

      if (rdCountsLoc < 0) {
        rdCounts[n] = rdCountsLoc;

        //  error
        exitg1 = true;
      } else {

#pragma omp barrier

        //  barrier here
        if (parts[lvl].shared_ents.size(0) != 0) {

#pragma omp single nowait

          { //single
            e_assemble_body_kernel(xs, conn, stcls, n2cPtr, n2cList, degree,
              rowPtr, colInd, vals, nnzPr, parts[lvl].shared_ents, 1, parts[lvl]
              .shared_ents.size(0), &wls[n], rdTags, &rdCountsLoc);
            if (rdCountsLoc < 0) {
              rdCounts[n] = rdCountsLoc;

              //  error
            } else {
              rdCounts[n] = rdCounts[n] + rdCountsLoc;
            }

          } //single

          exitg1 = true;
        } else {
          lvl++;
        }
      }
    }
  }

  static void b_assemble_body_task(const ::coder::array<double, 2U> &xs, const ::
    coder::array<int, 2U> &conn, const ::coder::array<int, 2U> &stcls, const ::
    coder::array<int, 1U> &n2cPtr, const ::coder::array<int, 1U> &n2cList, int
    degree, boolean_T interp0, boolean_T useDag, const ::coder::array<int, 1U>
    &rowPtr, const ::coder::array<int, 1U> &colInd, ::coder::array<double, 1U>
    &vals, const ::coder::array<int, 1U> &nnzPr, const ::coder::array<Omp4mPart,
    1U> &parts, ::coder::array<WlsDataStruct, 1U> &wls, ::coder::array<boolean_T,
    1U> &rdTags, ::coder::array<int, 1U> &rdCounts)
  {
    int lvl;
    int n;
    boolean_T exitg1;

    //  kernel for running in task-based parallelism
    n = 0;

#ifdef _OPENMP

    n = omp_get_thread_num();

#endif //_OPENMP

    rdCounts[n] = 0;
    b_WlsDataStruct(degree, interp0, useDag, &wls[n].wlsObj, &wls[n].wlsWgts,
                    wls[n].coeffs, &wls[n].interp0, &wls[n].useDag);

#pragma omp barrier

    lvl = 0;
    exitg1 = false;
    while ((!exitg1) && (lvl <= parts.size(0) - 1)) {
      int i;

#pragma omp single

      { //single
        i = parts[lvl].nparts;
        for (int mypart{0}; mypart < i; mypart++) {

#pragma omp task default(shared)

          { //task
            b_assemble_body_task_kernel(xs, conn, stcls, n2cPtr, n2cList, degree,
              rowPtr, colInd, vals, nnzPr, parts, mypart + 1, lvl + 1, wls,
              rdTags, rdCounts);

          } //task
        }

      } //single

      //  implicit barrier
      if (parts[lvl].shared_ents.size(0) != 0) {
        int rdCountsLoc;

#pragma omp single nowait

        { //single
          rdCountsLoc = h_assemble_body_kernel(xs, conn, stcls, n2cPtr, n2cList,
            degree, rowPtr, colInd, vals, nnzPr, parts[lvl].shared_ents, 1,
            parts[lvl].shared_ents.size(0), &wls[n], rdTags);
          if (rdCountsLoc < 0) {
            rdCounts[n] = rdCountsLoc;

            //  error
          } else {
            rdCounts[n] = rdCounts[n] + rdCountsLoc;
          }

        } //single

        exitg1 = true;
      } else {
        lvl++;
      }
    }
  }

  static void b_assemble_body_task_kernel(const ::coder::array<double, 2U> &xs,
    const ::coder::array<int, 2U> &conn, const ::coder::array<int, 2U> &stcls,
    const ::coder::array<int, 1U> &n2cPtr, const ::coder::array<int, 1U>
    &n2cList, int degree, const ::coder::array<int, 1U> &rowPtr, const ::coder::
    array<int, 1U> &colInd, ::coder::array<double, 1U> &vals, const ::coder::
    array<int, 1U> &nnzPr, const ::coder::array<Omp4mPart, 1U> &parts, int
    mypart, int lvl, ::coder::array<WlsDataStruct, 1U> &wls, ::coder::array<
    boolean_T, 1U> &rdTags, ::coder::array<int, 1U> &rdCounts)
  {
    int n;

    //  Gets the thread number of the thread, within the team, making this call.
    n = 0;

#ifdef _OPENMP

    n = omp_get_thread_num();

#endif //_OPENMP

    //  if a previously assembled operator had a LAPACK error, then stop
    if (rdCounts[n] >= 0) {
      int rdCountsLoc;

      rdCountsLoc = h_assemble_body_kernel(xs, conn, stcls, n2cPtr, n2cList,
        degree, rowPtr, colInd, vals, nnzPr, parts[lvl - 1].part_list, parts[lvl
        - 1].part_ptr[mypart - 1], parts[lvl - 1].part_ptr[mypart] - 1, &wls[n],
        rdTags);
      if (rdCountsLoc < 0) {
        rdCounts[n] = rdCountsLoc;

        //  error
      } else {
        rdCounts[n] = rdCounts[n] + rdCountsLoc;
      }
    }
  }

  static void b_assemble_surf_kernel(const ::coder::array<double, 2U> &xs, const
    ::coder::array<int, 2U> &conn, const ::coder::array<int, 2U> &stcls, const ::
    coder::array<int, 1U> &n2cPtr, const ::coder::array<int, 1U> &n2cList, int
    surfType, int degree, const ::coder::array<double, 2U> &nrms, const ::coder::
    array<int, 1U> &rowPtr, const ::coder::array<int, 1U> &colInd, ::coder::
    array<double, 1U> &vals, const ::coder::array<int, 1U> &nnzPr, const ::coder::
    array<int, 1U> &nrange, int istart, int iend, WlsDataStruct *wls, ::coder::
    array<boolean_T, 1U> &rdTags, int *rdCounts)
  {
    double xsLocal[1536];
    double us[1024];
    double Ns[512];
    double nrm0_idx_0;
    double nrm0_idx_1;
    double nrm0_idx_2;
    int eids[512];
    int gIDs[512];
    int b_iv[2];
    int iv1[2];
    int i;
    boolean_T exitg1;
    boolean_T isSphSurf;

    //  kernel for surface computation
    wls->wlsWgts.params_pointwise.size[1] = 1;
    wls->wlsWgts.params_pointwise.size[0] = 512;

    //  check if spherical surface
    isSphSurf = (surfType == 1);

    //  loop begins
    i = istart - 1;
    exitg1 = false;
    while ((!exitg1) && (i + 1 <= iend)) {
      double a;
      double d;
      double d1;
      double d2;
      double d3;
      double ori_idx_0;
      double ori_idx_1;
      double ori_idx_2;
      double t_idx_0;
      double t_idx_1;
      double t_idx_2;
      double t_idx_3;
      double t_idx_4;
      double t_idx_5;
      int b_xsLocal_tmp;
      int k;
      int nPoints;
      int nid;
      int us_tmp;
      int xsLocal_tmp;

      //  get nodal ID
      nid = nrange[i] - 1;
      nPoints = stcls[(stcls.size(1) + stcls.size(1) * (nrange[i] - 1)) - 1] - 1;

      //  fetch data to local buffer and localize the coordinates
      for (int j{0}; j <= nPoints; j++) {
        gIDs[j] = stcls[j + stcls.size(1) * nid];
      }

      //  fetch coordinates
      for (k = 0; k <= nPoints; k++) {
        xsLocal_tmp = gIDs[k];
        xsLocal[3 * k] = xs[xs.size(1) * (xsLocal_tmp - 1)];
        xsLocal[3 * k + 1] = xs[xs.size(1) * (xsLocal_tmp - 1) + 1];
        xsLocal[3 * k + 2] = xs[xs.size(1) * (xsLocal_tmp - 1) + 2];
      }

      //  origin
      ori_idx_0 = xsLocal[0];
      ori_idx_1 = xsLocal[1];
      ori_idx_2 = xsLocal[2];

      //  localize
      for (k = 0; k <= nPoints; k++) {
        xsLocal[3 * k] -= ori_idx_0;
        b_xsLocal_tmp = 3 * k + 1;
        xsLocal[b_xsLocal_tmp] -= ori_idx_1;
        b_xsLocal_tmp = 3 * k + 2;
        xsLocal[b_xsLocal_tmp] -= ori_idx_2;
      }

      //  get normal direction
      if ((nrms.size(0) != 0) && (nrms.size(1) != 0)) {
        nrm0_idx_0 = nrms[nrms.size(1) * (nrange[i] - 1)];
        nrm0_idx_1 = nrms[nrms.size(1) * (nrange[i] - 1) + 1];
        nrm0_idx_2 = nrms[nrms.size(1) * (nrange[i] - 1) + 2];
      } else if (isSphSurf) {
        //  spherical or circle
        nrm0_idx_0 = xs[xs.size(1) * (nrange[i] - 1)];
        nrm0_idx_1 = xs[xs.size(1) * (nrange[i] - 1) + 1];
        nrm0_idx_2 = xs[xs.size(1) * (nrange[i] - 1) + 2];
      }

      //  compute tagent
      d = std::abs(nrm0_idx_0);
      if ((d > std::abs(nrm0_idx_1)) && (d > std::abs(nrm0_idx_2))) {
        t_idx_0 = -nrm0_idx_0 * nrm0_idx_1;
        t_idx_2 = 1.0 - nrm0_idx_1 * nrm0_idx_1;
        t_idx_4 = -nrm0_idx_1 * nrm0_idx_2;
      } else {
        t_idx_0 = 1.0 - nrm0_idx_0 * nrm0_idx_0;
        t_idx_2 = -nrm0_idx_0 * nrm0_idx_1;
        t_idx_4 = -nrm0_idx_0 * nrm0_idx_2;
      }

      a = std::sqrt((t_idx_0 * t_idx_0 + t_idx_2 * t_idx_2) + t_idx_4 * t_idx_4);
      t_idx_0 /= a;
      t_idx_2 /= a;
      t_idx_4 /= a;

      //  cross
      t_idx_1 = nrm0_idx_1 * t_idx_4 - nrm0_idx_2 * t_idx_2;
      t_idx_3 = t_idx_0 * nrm0_idx_2 - nrm0_idx_0 * t_idx_4;
      t_idx_5 = nrm0_idx_0 * t_idx_2 - t_idx_0 * nrm0_idx_1;

      //  project onto the tangent plan
      for (k = 0; k <= nPoints; k++) {
        us_tmp = k << 1;
        us[us_tmp] = 0.0;
        us[us_tmp + 1] = 0.0;
      }

      //  matrix-matrix, using mem-efficient loop
      for (int ii{0}; ii <= nPoints; ii++) {
        d = xsLocal[3 * ii];
        d1 = xsLocal[3 * ii + 1];
        xsLocal_tmp = ii << 1;
        d2 = (us[xsLocal_tmp + 1] + d * t_idx_1) + d1 * t_idx_3;
        d3 = xsLocal[3 * ii + 2];
        us[xsLocal_tmp] = ((us[xsLocal_tmp] + d * t_idx_0) + t_idx_2 * d1) +
          t_idx_4 * d3;
        d2 += t_idx_5 * d3;
        us[xsLocal_tmp + 1] = d2;
      }

      //  compute scaling to the WLS weights by inner product of normals
      if ((nrms.size(0) != 0) && (nrms.size(1) != 0)) {
        for (int j{2}; j <= nPoints + 1; j++) {
          us_tmp = stcls[(j + stcls.size(1) * nid) - 1] - 1;
          a = (nrm0_idx_0 * nrms[nrms.size(1) * us_tmp] + nrm0_idx_1 *
               nrms[nrms.size(1) * us_tmp + 1]) + nrm0_idx_2 * nrms[nrms.size(1)
            * us_tmp + 2];
          if (a < 0.0) {
            a = 0.0;
          }

          wls->wlsWgts.params_pointwise.data[j - 1] = a;
        }
      } else if (isSphSurf) {
        //  spherical or circle
        for (int j{2}; j <= nPoints + 1; j++) {
          int idx;
          idx = stcls[(j + stcls.size(1) * nid) - 1] - 1;
          a = (nrm0_idx_0 * xs[xs.size(1) * idx] + nrm0_idx_1 * xs[xs.size(1) *
               idx + 1]) + nrm0_idx_2 * xs[xs.size(1) * idx + 2];
          if (a < 0.0) {
            a = 0.0;
          }

          wls->wlsWgts.params_pointwise.data[j - 1] = a;
        }
      }

      wls->wlsWgts.params_pointwise.data[0] = 1.0;

      //  compute wls
      b_iv[0] = 512;
      b_iv[1] = 2;
      wls_init(&wls->wlsObj, us, b_iv, wls->wlsWgts.name.data,
               wls->wlsWgts.params_pointwise.data,
               wls->wlsWgts.params_pointwise.size, degree, wls->interp0,
               wls->useDag, nPoints + 1);
      if (wls->wlsObj.rank < 0) {
        *rdCounts = -nrange[i];

        //  record error node and return
        exitg1 = true;
      } else {
        if (wls->wlsObj.rank != wls->wlsObj.ncols) {
          (*rdCounts)++;

          //  increment counter for RD nodes
          rdTags[nrange[i] - 1] = true;
        } else {
          int nCells;

          //  get centers and shape function Ns at centers
          nCells = (n2cPtr[nrange[i]] - n2cPtr[nrange[i] - 1]) - 1;
          for (int j{0}; j <= nCells; j++) {
            int eid;
            int npts;
            eid = n2cList[(n2cPtr[nid] + j) - 1] - 1;

            eids[j] = eid + 1;
            for (npts = conn.size(1); conn[(npts + conn.size(1) * eid) - 1] <= 0;
                 npts--) {
            }

            Ns[j] = 1.0 / static_cast<double>(npts);

            //  compute center
            xsLocal[3 * j] = xs[xs.size(1) * (conn[conn.size(1) * eid] - 1)];
            b_xsLocal_tmp = 3 * j + 1;
            xsLocal[b_xsLocal_tmp] = xs[xs.size(1) * (conn[conn.size(1) * eid] -
              1) + 1];
            us_tmp = 3 * j + 2;
            xsLocal[us_tmp] = xs[xs.size(1) * (conn[conn.size(1) * eid] - 1) + 2];

            //  C++
            for (int ii{2}; ii <= npts; ii++) {
              xsLocal_tmp = conn[(ii + conn.size(1) * eid) - 1] - 1;
              xsLocal[3 * j] += xs[xs.size(1) * xsLocal_tmp];
              xsLocal[b_xsLocal_tmp] += xs[xs.size(1) * xsLocal_tmp + 1];
              xsLocal[us_tmp] += xs[xs.size(1) * xsLocal_tmp + 2];
            }

            d = Ns[j];
            d1 = xsLocal[3 * j] * d;
            d2 = xsLocal[b_xsLocal_tmp] * d;
            d *= xsLocal[us_tmp];
            if (isSphSurf) {
              a = std::sqrt((d1 * d1 + d2 * d2) + d * d);
              d1 /= a;
              d2 /= a;
              d /= a;
            }

            //  localize
            d1 -= ori_idx_0;
            xsLocal[3 * j] = d1;
            d2 -= ori_idx_1;
            xsLocal[b_xsLocal_tmp] = d2;
            d -= ori_idx_2;
            xsLocal[us_tmp] = d;
          }

          for (k = 0; k <= nCells; k++) {
            us_tmp = k << 1;
            us[us_tmp] = 0.0;
            us[us_tmp + 1] = 0.0;
          }

          //  matrix-matrix, using mem-efficient loop
          for (int ii{0}; ii <= nCells; ii++) {
            d = xsLocal[3 * ii];
            d1 = xsLocal[3 * ii + 1];
            xsLocal_tmp = ii << 1;
            d2 = (us[xsLocal_tmp + 1] + d * t_idx_1) + d1 * t_idx_3;
            d3 = xsLocal[3 * ii + 2];
            us[xsLocal_tmp] = ((us[xsLocal_tmp] + d * t_idx_0) + t_idx_2 * d1) +
              t_idx_4 * d3;
            d2 += t_idx_5 * d3;
            us[xsLocal_tmp + 1] = d2;
          }

          iv1[0] = 512;
          iv1[1] = 2;
          wls_func(&wls->wlsObj, us, iv1, nCells + 1, wls->coeffs);

          //  Compute local OSUS coefficients
          for (k = 0; k <= nPoints; k++) {
            for (int j{0}; j <= nCells; j++) {
              wls->coeffs[j + wls->coeffs.size(1) * k] = wls->coeffs[j +
                wls->coeffs.size(1) * k] * Ns[j];
            }
          }

          //  substract linear interp
          for (int j{0}; j <= nCells; j++) {
            wls->coeffs[j] = wls->coeffs[j] - Ns[j];
          }

          //  kernel for updating crs matrix
          for (int j{0}; j <= nCells; j++) {
            int rowIdx;
            rowIdx = eids[j] - 1;
            for (int ii{0}; ii <= nPoints; ii++) {
              int b_i;
              us_tmp = gIDs[ii];
              b_i = rowPtr[rowIdx] - 1;
              xsLocal_tmp = (rowPtr[rowIdx] + nnzPr[eids[j] - 1]) - 1;
              int exitg2;
              do {
                exitg2 = 0;
                if (b_i + 1 <= xsLocal_tmp) {
                  if (colInd[b_i] == us_tmp) {
                    k = b_i;
                    exitg2 = 1;
                  } else {
                    b_i++;
                  }
                } else {
                  k = -1;

                  m2cErrMsgIdAndTxt("add_crs:missingIndex",
                                    "could not find column index %d", us_tmp);
                  exitg2 = 1;
                }
              } while (exitg2 == 0);

              vals[k] = vals[k] + wls->coeffs[j + wls->coeffs.size(1) * ii];
            }
          }
        }

        i++;
      }
    }
  }

  static void b_assemble_surf_par(const ::coder::array<double, 2U> &xs, const ::
    coder::array<int, 2U> &conn, const ::coder::array<int, 2U> &stcls, const ::
    coder::array<int, 1U> &n2cPtr, const ::coder::array<int, 1U> &n2cList, int
    surfType, int degree, boolean_T interp0, boolean_T useDag, const ::coder::
    array<double, 2U> &nrms, const ::coder::array<int, 1U> &rowPtr, const ::
    coder::array<int, 1U> &colInd, ::coder::array<double, 1U> &vals, const ::
    coder::array<int, 1U> &nnzPr, const ::coder::array<Omp4mPart, 1U> &parts, ::
    coder::array<WlsDataStruct, 1U> &wls, ::coder::array<boolean_T, 1U> &rdTags,
    ::coder::array<int, 1U> &rdCounts, Omp4mPartContext *partContext)
  {
    int lvl;
    int n;
    int partid;
    int rdCountsLoc;
    boolean_T exitg1;

    //  kernel for parallel
    n = 0;

#ifdef _OPENMP

    n = omp_get_thread_num();

#endif //_OPENMP

    rdCounts[n] = 0;
    rdCountsLoc = 0;
    b_WlsDataStruct(degree, interp0, useDag, &wls[n].wlsObj, &wls[n].wlsWgts,
                    wls[n].coeffs, &wls[n].interp0, &wls[n].useDag);

#pragma omp single

    { //single
      partContext->atomic_counters.set_size(parts.size(0));
      partid = parts.size(0);
      for (int i{0}; i < partid; i++) {
        partContext->atomic_counters[i] = 1;
      }

    } //single

    //  safe and shared across threads
    lvl = 0;
    exitg1 = false;
    while ((!exitg1) && (lvl <= parts.size(0) - 1)) {
      //  Assemble for this part
      int exitg2;
      do {
        exitg2 = 0;

#pragma omp atomic capture

        {                              // atomic capture();
          //  to prevent Coder putting this line before pragma
          partid = (partContext->atomic_counters[lvl]);
          partContext->atomic_counters[lvl] = partContext->atomic_counters[lvl]
            + 1;
        }                              // atomic capture();

        //  to prevent Coder detecting v is for index and automatically subtract
        if ((lvl + 1 <= parts.size(0)) && (partid <= parts[lvl].nparts)) {
          d_assemble_surf_kernel(xs, conn, stcls, n2cPtr, n2cList, surfType,
            degree, nrms, rowPtr, colInd, vals, nnzPr, parts[lvl].part_list,
            parts[lvl].part_ptr[partid - 1], parts[lvl].part_ptr[partid] - 1,
            &wls[n], rdTags, &rdCountsLoc);
          if (rdCountsLoc < 0) {
            exitg2 = 1;
          } else {
            rdCounts[n] = rdCounts[n] + rdCountsLoc;
          }
        } else {
          exitg2 = 1;
        }
      } while (exitg2 == 0);

      if (rdCountsLoc < 0) {
        rdCounts[n] = rdCountsLoc;

        //  error
        exitg1 = true;
      } else {

#pragma omp barrier

        //  barrier here to avoid reset token before all break
        if (parts[lvl].shared_ents.size(0) != 0) {

#pragma omp single nowait

          { //single
            d_assemble_surf_kernel(xs, conn, stcls, n2cPtr, n2cList, surfType,
              degree, nrms, rowPtr, colInd, vals, nnzPr, parts[lvl].shared_ents,
              1, parts[lvl].shared_ents.size(0), &wls[n], rdTags, &rdCountsLoc);
            if (rdCountsLoc < 0) {
              rdCounts[n] = rdCountsLoc;

              //  error
            } else {
              rdCounts[n] = rdCounts[n] + rdCountsLoc;
            }

          } //single

          exitg1 = true;
        } else {
          lvl++;
        }
      }
    }
  }

  static void b_assemble_surf_task(const ::coder::array<double, 2U> &xs, const ::
    coder::array<int, 2U> &conn, const ::coder::array<int, 2U> &stcls, const ::
    coder::array<int, 1U> &n2cPtr, const ::coder::array<int, 1U> &n2cList, int
    surfType, int degree, boolean_T interp0, boolean_T useDag, const ::coder::
    array<double, 2U> &nrms, const ::coder::array<int, 1U> &rowPtr, const ::
    coder::array<int, 1U> &colInd, ::coder::array<double, 1U> &vals, const ::
    coder::array<int, 1U> &nnzPr, const ::coder::array<Omp4mPart, 1U> &parts, ::
    coder::array<WlsDataStruct, 1U> &wls, ::coder::array<boolean_T, 1U> &rdTags,
    ::coder::array<int, 1U> &rdCounts)
  {
    int lvl;
    int n;
    int rdCountsLoc;
    boolean_T exitg1;

    n = 0;

#ifdef _OPENMP

    n = omp_get_thread_num();

#endif //_OPENMP

    rdCounts[n] = 0;
    b_WlsDataStruct(degree, interp0, useDag, &wls[n].wlsObj, &wls[n].wlsWgts,
                    wls[n].coeffs, &wls[n].interp0, &wls[n].useDag);

#pragma omp barrier

    lvl = 0;
    exitg1 = false;
    while ((!exitg1) && (lvl <= parts.size(0) - 1)) {
      int i;

#pragma omp single

      { //single
        i = parts[lvl].nparts;
        for (int mypart{0}; mypart < i; mypart++) {

#pragma omp task default(shared)

          { //task
            b_assemble_surf_task_kernel(xs, conn, stcls, n2cPtr, n2cList,
              surfType, degree, nrms, rowPtr, colInd, vals, nnzPr, parts, mypart
              + 1, lvl + 1, wls, rdTags, rdCounts);

          } //task
        }

      } //single

      //  implicit barrier
      if (parts[lvl].shared_ents.size(0) != 0) {

#pragma omp single nowait

        { //single
          rdCountsLoc = 0;
          d_assemble_surf_kernel(xs, conn, stcls, n2cPtr, n2cList, surfType,
            degree, nrms, rowPtr, colInd, vals, nnzPr, parts[lvl].shared_ents, 1,
            parts[lvl].shared_ents.size(0), &wls[n], rdTags, &rdCountsLoc);
          if (rdCountsLoc < 0) {
            rdCounts[n] = rdCountsLoc;

            //  error
          } else {
            rdCounts[n] = rdCounts[n] + rdCountsLoc;
          }

        } //single

        exitg1 = true;
      } else {
        lvl++;
      }
    }
  }

  static void b_assemble_surf_task_kernel(const ::coder::array<double, 2U> &xs,
    const ::coder::array<int, 2U> &conn, const ::coder::array<int, 2U> &stcls,
    const ::coder::array<int, 1U> &n2cPtr, const ::coder::array<int, 1U>
    &n2cList, int surfType, int degree, const ::coder::array<double, 2U> &nrms,
    const ::coder::array<int, 1U> &rowPtr, const ::coder::array<int, 1U> &colInd,
    ::coder::array<double, 1U> &vals, const ::coder::array<int, 1U> &nnzPr,
    const ::coder::array<Omp4mPart, 1U> &parts, int mypart, int lvl, ::coder::
    array<WlsDataStruct, 1U> &wls, ::coder::array<boolean_T, 1U> &rdTags, ::
    coder::array<int, 1U> &rdCounts)
  {
    int n;
    int rdCountsLoc;

    //  Gets the thread number of the thread, within the team, making this call.
    n = 0;

#ifdef _OPENMP

    n = omp_get_thread_num();

#endif //_OPENMP

    //  if a previously assembled operator had a LAPACK error, then stop
    if (rdCounts[n] >= 0) {
      rdCountsLoc = 0;
      d_assemble_surf_kernel(xs, conn, stcls, n2cPtr, n2cList, surfType, degree,
        nrms, rowPtr, colInd, vals, nnzPr, parts[lvl - 1].part_list, parts[lvl -
        1].part_ptr[mypart - 1], parts[lvl - 1].part_ptr[mypart] - 1, &wls[n],
        rdTags, &rdCountsLoc);
      if (rdCountsLoc < 0) {
        rdCounts[n] = rdCountsLoc;

        //  error
      } else {
        rdCounts[n] = rdCounts[n] + rdCountsLoc;
      }
    }
  }

  static void b_wls_init(WlsObject *wls, const double us_data[], int us_size,
    const char weight_name_data[], const double weight_params_pointwise_data[],
    const int weight_params_pointwise_size[2], int degree, boolean_T interp0,
    boolean_T use_dag, int npoints)
  {
    ::coder::array<unsigned char, 1U> r;
    double maxx;
    double maxx_inv;
    double thres;
    int b_interp0;
    int i;
    int loop_ub;

    //  wls_init  Initialize WlsObject in 1D, 2D, or 3D.
    m2cAssert(true, "");

    //  Process input arguments
    wls->interp0 = interp0;
    b_interp0 = wls->interp0;
    wls->use_dag = use_dag;
    if (npoints <= 0) {
      npoints = us_size;
    } else {
      m2cAssert(npoints <= us_size,
                "Number of points cannot be greater than the first dimension of `us`.");
    }

    //  Resize buffers
    b_wls_resize(wls, npoints, degree, use_dag);

    //  Recompute DAG if use_dag and its signature does not match
    if (use_dag) {
      i = wls->dag.size(1) * wls->dag.size(0) - 1;
      loop_ub = wls->dag.size(0);
      if (wls->dag[i % loop_ub * wls->dag.size(1) + i / loop_ub] != degree + 127)
      {
        //  Wrapper function for building DAG in nD.
        gen_vander_1d_dag(degree, r);
        wls->dag.set_size(r.size(0), 1);
        loop_ub = r.size(0);
        for (i = 0; i < loop_ub; i++) {
          wls->dag[i] = r[i];
        }
      }
    }

    if (wls->interp0 != 0) {
      //  Make the first node the origin in interp0 mode
      wls->origin.size[1] = 1;
      wls->origin.size[0] = 1;
      wls->origin.data[0] = us_data[0];
      for (int b_i{0}; b_i < npoints; b_i++) {
        i = wls->us.size(0);
        wls->us[b_i % i * wls->us.size(1) + b_i / i] = us_data[b_i] - us_data[0];
      }
    } else {
      wls->origin.size[1] = 1;
      wls->origin.size[0] = 1;
      wls->origin.data[0] = 0.0;
      for (int b_i{0}; b_i < npoints; b_i++) {
        i = wls->us.size(0);
        wls->us[b_i % i * wls->us.size(1) + b_i / i] = us_data[b_i];
      }
    }

    //  Scale us to be between -1 and 1
    maxx = 0.0;
    for (int b_i{0}; b_i < npoints; b_i++) {
      i = wls->us.size(0);
      maxx = std::fmax(maxx, std::abs(wls->us[b_i % i * wls->us.size(1) + b_i /
        i]));
    }

    if (maxx == 0.0) {
      maxx_inv = 1.0;
    } else {
      maxx_inv = 1.0 / maxx;
    }

    wls->hs_inv.data[0] = maxx_inv;

    //  scale wls.us
    if (maxx_inv != 1.0) {
      for (int b_i{0}; b_i < npoints; b_i++) {
        i = wls->us.size(0);
        loop_ub = wls->us.size(0);
        wls->us[b_i % i * wls->us.size(1) + b_i / i] = wls->us[b_i % loop_ub *
          wls->us.size(1) + b_i / loop_ub] * maxx_inv;
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
        loop_ub = wls->rweights.size(0);
        wls->rweights.set_size(loop_ub);
        for (i = 0; i < loop_ub; i++) {
          wls->rweights[i] = 1.0;
        }
      } else if ((weight_name_data[0] == 'I') || (weight_name_data[0] == 'i')) {
        //  inverse distance
        wls_invdist_weights(wls->us, npoints, degree,
                            weight_params_pointwise_data,
                            weight_params_pointwise_size, wls->rweights);
      } else if ((weight_name_data[0] == 'B') || (weight_name_data[0] == 'b')) {
        //  Buhmann weights. All points share same parameters
        wls_buhmann_weights(wls->us, npoints, degree,
                            weight_params_pointwise_data,
                            weight_params_pointwise_size, wls->rweights);
      } else {
        if ((weight_name_data[0] != 'E') && (weight_name_data[0] != 'e')) {
          m2cErrMsgIdAndTxt("wlslib:WrongWeightName",
                            "Weighting scheme must be unit, inverse, Buhmann, or WLS-ENO.");
        }

        //  WLS-ENO
        m2cAssert(false, "first two shared parameters are required");

        m2cAssert(weight_params_pointwise_size[0] >= npoints,
                  "size(params_pw,1) should be >=npoints");

        m2cAssert(false, "size(params_pw,2) should be >=2");
        wls->rweights.set_size(npoints);

        //  Compute hbar using ws as buffer space
        if (degree >= 0) {
          //  Compute 2-norm
          for (int b_i{0}; b_i < npoints; b_i++) {
            wls->rweights[b_i] = std::sqrt(us_data[b_i] * us_data[b_i]);
          }
        } else {
          //  Compute inf-norm for tensor-product
          for (int b_i{0}; b_i < npoints; b_i++) {
            wls->rweights[b_i] = std::abs(us_data[b_i]);
          }
        }

        //  Evaluate the inverse-distance weights as base
        wls_invdist_weights(wls->us, npoints, 0.5 - static_cast<double>(degree <
          0), wls->rweights);

        // A check that is always false is detected at compile-time. Eliminating code that follows.
      }
    }

    //  Compute Vandermonde system and recompute DAG if needed
    gen_vander(wls->us, npoints, degree, wls->rweights, wls->V);
    wls->ncols = wls->V.size(0);

    //  Compact CVM if needed
    wls->nrows = wls->V.size(1) / wls->stride * npoints;

    //  Omit rows in CVM if needed
    if ((degree > 1) && (degree < 7)) {
      thres = dv[degree - 1];
    } else {
      thres = 1.0E+8;
    }

    //  In interp0 mode, we trim off the first row and first column.
    rrqr_factor(wls->V, thres, b_interp0, b_interp0, wls->nrows - b_interp0,
                wls->ncols - b_interp0, wls->QR, wls->jpvt, &wls->rank,
                wls->work);
  }

  static void b_wls_resize(WlsObject *wls, int npoints, int degree, boolean_T
    use_dag)
  {
    int stride;
    int stride_idx_0;
    int y;

    //  wls_resize  Reinitialize the buffers of WlsObject
    wls->degree = degree;
    wls->order = 0;
    wls->use_dag = use_dag;

    //  make stride a multiple of four
    stride = ((npoints + 3) / 4) << 2;
    wls->stride = stride;
    wls->us.set_size(stride, 1);
    wls->rweights.set_size(stride);
    wls->npoints = npoints;

    //  determine number of columns and allocate V and QR
    if (degree < 0) {
      y = 2 - degree;
    } else {
      y = degree + 2;
    }

    wls->hs_inv.size[1] = 1;
    wls->hs_inv.size[0] = 1;
    if (use_dag) {
      if ((wls->dag.size(1) != 1) || (y != wls->dag.size(0))) {
        //  Reset DAG if dimension or degree has changed.
        wls->dag.set_size(y, 1);
        stride_idx_0 = wls->dag.size(0);
        wls->dag[(y - 1) % stride_idx_0 * wls->dag.size(1) + (y - 1) /
          stride_idx_0] = MAX_uint8_T;
      }
    } else {
      wls->dag.set_size(0, 1);
    }

    wls->jpvt.set_size(y - 1);

    //  V is always full, but QR has one fewer row and column in interp0 mode
    wls->V.set_size(y - 1, stride);
    stride_idx_0 = y - wls->interp0;
    wls->QR.set_size(stride_idx_0, stride);
    wls->rank = 0;

    //  work space
    stride_idx_0 = (y - 1) << 2;
    if (stride >= y) {
      y = stride;
    }

    if (stride_idx_0 < 4160) {
      stride_idx_0 = 4160;
    }

    wls->work.set_size((y << 5) + stride_idx_0);
  }

  static void build_part(int nParts, const ::coder::array<int, 1U> &nparts,
    const ::coder::array<int, 1U> &cparts, const ::coder::array<int, 1U> &eptr,
    const ::coder::array<int, 1U> &eind, ::coder::array<boolean_T, 1U> &ctags, ::
    coder::array<int, 1U> &iwork, ::coder::array<int, 1U> &partptr, ::coder::
    array<int, 1U> &partlist, ::coder::array<int, 1U> &sharedents)
  {
    int i;
    int j;
    int k;
    int m;
    int n0;

    //  local function to compute the partitioning
    m = cparts.size(0) - 1;
    partptr.set_size(nParts + 1);
    for (i = 0; i <= nParts; i++) {
      partptr[i] = 0;
    }

    partptr[0] = 1;
    for (int b_i{0}; b_i <= m; b_i++) {
      partptr[cparts[b_i]] = partptr[cparts[b_i]] + 1;
    }

    //  accumulate
    for (int b_i{0}; b_i < nParts; b_i++) {
      partptr[b_i + 1] = partptr[b_i + 1] + partptr[b_i];
    }

    if (iwork.size(0) < partptr[nParts] - 1) {
      iwork.set_size(partptr[nParts] - 1);
    }

    for (int b_i{0}; b_i <= m; b_i++) {
      i = partptr[cparts[b_i] - 1];
      iwork[i - 1] = b_i + 1;
      partptr[cparts[b_i] - 1] = i + 1;
    }

    //  reset
    for (int b_i{nParts}; b_i >= 1; b_i--) {
      partptr[b_i] = partptr[b_i - 1];
    }

    partptr[0] = 1;

    //  determine shared region
    for (int part{0}; part < nParts; part++) {
      int i1;
      i = partptr[part];
      i1 = partptr[part + 1] - 1;
      for (int b_i{i}; b_i <= i1; b_i++) {
        int j_tmp;
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
    for (int b_i{0}; b_i < nParts; b_i++) {
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
    for (int b_i{0}; b_i <= m; b_i++) {
      if (!ctags[iwork[b_i] - 1]) {
        k++;
        partlist[k] = iwork[b_i];
      }
    }

    //  get shared cells
    sharedents.set_size(n0);
    k = -1;
    for (int b_i{0}; b_i <= m; b_i++) {
      if (ctags[b_i]) {
        k++;
        sharedents[k] = b_i + 1;
        ctags[b_i] = false;

        //  reset
      }
    }
  }

  static int c_assemble_body_kernel(const ::coder::array<double, 2U> &xs, const ::
    coder::array<int, 2U> &stcls, int degree, int iend, WlsDataStruct *wls)
  {
    double xsLocal[64];
    int gIDs[64];
    int rdCounts;
    boolean_T exitg1;
    rdCounts = 0;

    //  kernel
    exitg1 = false;
    while ((!exitg1) && (iend >= 1)) {
      double ori;
      int nPoints;

      //  get nodal ID
      nPoints = stcls[stcls.size(1) - 1] - 1;

      //  fetch coords, global IDs, and then localize coordinates
      for (int j{0}; j <= nPoints; j++) {
        gIDs[j] = stcls[j];
      }

      //  fetch coordinates
      for (int k{0}; k <= nPoints; k++) {
        xsLocal[k] = xs[xs.size(1) * (gIDs[k] - 1)];
      }

      //  origin
      ori = xsLocal[0];

      //  localize
      for (int k{0}; k <= nPoints; k++) {
        xsLocal[k] -= ori;
      }

      //  compute wls
      b_wls_init(&wls->wlsObj, xsLocal, 64, wls->wlsWgts.name.data,
                 wls->wlsWgts.params_pointwise.data,
                 wls->wlsWgts.params_pointwise.size, degree, wls->interp0,
                 wls->useDag, nPoints + 1);
      if (wls->wlsObj.rank < 0) {
        //  record error node and return
        exitg1 = true;
      }
    }

    return rdCounts;
  }

  static void c_assemble_body_par(const ::coder::array<double, 2U> &xs, const ::
    coder::array<int, 2U> &conn, const ::coder::array<int, 2U> &stcls, const ::
    coder::array<int, 1U> &n2cPtr, const ::coder::array<int, 1U> &n2cList, int
    degree, boolean_T interp0, boolean_T useDag, const ::coder::array<int, 1U>
    &rowPtr, const ::coder::array<int, 1U> &colInd, ::coder::array<double, 1U>
    &vals, const ::coder::array<int, 1U> &nnzPr, const ::coder::array<Omp4mPart,
    1U> &parts, ::coder::array<WlsDataStruct, 1U> &wls, ::coder::array<boolean_T,
    1U> &rdTags, ::coder::array<int, 1U> &rdCounts, Omp4mPartContext
    *partContext)
  {
    int lvl;
    int n;
    int partid;
    int rdCountsLoc;
    boolean_T exitg1;

    //  parallel kernel
    n = 0;

#ifdef _OPENMP

    n = omp_get_thread_num();

#endif //_OPENMP

    rdCounts[n] = 0;
    rdCountsLoc = 0;
    b_WlsDataStruct(degree, interp0, useDag, &wls[n].wlsObj, &wls[n].wlsWgts,
                    wls[n].coeffs, &wls[n].interp0, &wls[n].useDag);

#pragma omp single

    { //single
      partContext->atomic_counters.set_size(parts.size(0));
      partid = parts.size(0);
      for (int i{0}; i < partid; i++) {
        partContext->atomic_counters[i] = 1;
      }

    } //single

    //  safe and shared across all threads
    lvl = 0;
    exitg1 = false;
    while ((!exitg1) && (lvl <= parts.size(0) - 1)) {
      //  assemble for this part
      int exitg2;
      do {
        exitg2 = 0;

#pragma omp atomic capture

        {                              // atomic capture();
          //  to prevent Coder putting this line before pragma
          partid = (partContext->atomic_counters[lvl]);
          partContext->atomic_counters[lvl] = partContext->atomic_counters[lvl]
            + 1;
        }                              // atomic capture();

        //  to prevent Coder detecting v is for index and automatically subtract
        if ((lvl + 1 <= parts.size(0)) && (partid <= parts[lvl].nparts)) {
          f_assemble_body_kernel(xs, conn, stcls, n2cPtr, n2cList, degree,
            rowPtr, colInd, vals, nnzPr, parts[lvl].part_list, parts[lvl].
            part_ptr[partid - 1], parts[lvl].part_ptr[partid] - 1, &wls[n],
            rdTags, &rdCountsLoc);
          if (rdCountsLoc < 0) {
            exitg2 = 1;
          } else {
            //  LAPACK error
            rdCounts[n] = rdCounts[n] + rdCountsLoc;
          }
        } else {
          exitg2 = 1;
        }
      } while (exitg2 == 0);

      if (rdCountsLoc < 0) {
        rdCounts[n] = rdCountsLoc;

        //  error
        exitg1 = true;
      } else {

#pragma omp barrier

        //  barrier here
        if (parts[lvl].shared_ents.size(0) != 0) {

#pragma omp single nowait

          { //single
            f_assemble_body_kernel(xs, conn, stcls, n2cPtr, n2cList, degree,
              rowPtr, colInd, vals, nnzPr, parts[lvl].shared_ents, 1, parts[lvl]
              .shared_ents.size(0), &wls[n], rdTags, &rdCountsLoc);
            if (rdCountsLoc < 0) {
              rdCounts[n] = rdCountsLoc;

              //  error
            } else {
              rdCounts[n] = rdCounts[n] + rdCountsLoc;
            }

          } //single

          exitg1 = true;
        } else {
          lvl++;
        }
      }
    }
  }

  static void c_assemble_body_task(const ::coder::array<double, 2U> &xs, const ::
    coder::array<int, 2U> &conn, const ::coder::array<int, 2U> &stcls, const ::
    coder::array<int, 1U> &n2cPtr, const ::coder::array<int, 1U> &n2cList, int
    degree, boolean_T interp0, boolean_T useDag, const ::coder::array<int, 1U>
    &rowPtr, const ::coder::array<int, 1U> &colInd, ::coder::array<double, 1U>
    &vals, const ::coder::array<int, 1U> &nnzPr, const ::coder::array<Omp4mPart,
    1U> &parts, ::coder::array<WlsDataStruct, 1U> &wls, ::coder::array<boolean_T,
    1U> &rdTags, ::coder::array<int, 1U> &rdCounts)
  {
    int lvl;
    int n;
    boolean_T exitg1;

    //  kernel for running in task-based parallelism
    n = 0;

#ifdef _OPENMP

    n = omp_get_thread_num();

#endif //_OPENMP

    rdCounts[n] = 0;
    b_WlsDataStruct(degree, interp0, useDag, &wls[n].wlsObj, &wls[n].wlsWgts,
                    wls[n].coeffs, &wls[n].interp0, &wls[n].useDag);

#pragma omp barrier

    lvl = 0;
    exitg1 = false;
    while ((!exitg1) && (lvl <= parts.size(0) - 1)) {
      int i;

#pragma omp single

      { //single
        i = parts[lvl].nparts;
        for (int mypart{0}; mypart < i; mypart++) {

#pragma omp task default(shared)

          { //task
            c_assemble_body_task_kernel(xs, conn, stcls, n2cPtr, n2cList, degree,
              rowPtr, colInd, vals, nnzPr, parts, mypart + 1, lvl + 1, wls,
              rdTags, rdCounts);

          } //task
        }

      } //single

      //  implicit barrier
      if (parts[lvl].shared_ents.size(0) != 0) {
        int rdCountsLoc;

#pragma omp single nowait

        { //single
          rdCountsLoc = i_assemble_body_kernel(xs, conn, stcls, n2cPtr, n2cList,
            degree, rowPtr, colInd, vals, nnzPr, parts[lvl].shared_ents, 1,
            parts[lvl].shared_ents.size(0), &wls[n], rdTags);
          if (rdCountsLoc < 0) {
            rdCounts[n] = rdCountsLoc;

            //  error
          } else {
            rdCounts[n] = rdCounts[n] + rdCountsLoc;
          }

        } //single

        exitg1 = true;
      } else {
        lvl++;
      }
    }
  }

  static void c_assemble_body_task_kernel(const ::coder::array<double, 2U> &xs,
    const ::coder::array<int, 2U> &conn, const ::coder::array<int, 2U> &stcls,
    const ::coder::array<int, 1U> &n2cPtr, const ::coder::array<int, 1U>
    &n2cList, int degree, const ::coder::array<int, 1U> &rowPtr, const ::coder::
    array<int, 1U> &colInd, ::coder::array<double, 1U> &vals, const ::coder::
    array<int, 1U> &nnzPr, const ::coder::array<Omp4mPart, 1U> &parts, int
    mypart, int lvl, ::coder::array<WlsDataStruct, 1U> &wls, ::coder::array<
    boolean_T, 1U> &rdTags, ::coder::array<int, 1U> &rdCounts)
  {
    int n;

    //  Gets the thread number of the thread, within the team, making this call.
    n = 0;

#ifdef _OPENMP

    n = omp_get_thread_num();

#endif //_OPENMP

    //  if a previously assembled operator had a LAPACK error, then stop
    if (rdCounts[n] >= 0) {
      int rdCountsLoc;

      rdCountsLoc = i_assemble_body_kernel(xs, conn, stcls, n2cPtr, n2cList,
        degree, rowPtr, colInd, vals, nnzPr, parts[lvl - 1].part_list, parts[lvl
        - 1].part_ptr[mypart - 1], parts[lvl - 1].part_ptr[mypart] - 1, &wls[n],
        rdTags);
      if (rdCountsLoc < 0) {
        rdCounts[n] = rdCountsLoc;

        //  error
      } else {
        rdCounts[n] = rdCounts[n] + rdCountsLoc;
      }
    }
  }

  static int c_assemble_surf_kernel(const ::coder::array<double, 2U> &xs, const ::
    coder::array<int, 2U> &conn, const ::coder::array<int, 2U> &stcls, const ::
    coder::array<int, 1U> &n2cPtr, const ::coder::array<int, 1U> &n2cList, int
    surfType, int degree, const ::coder::array<double, 2U> &nrms, const ::coder::
    array<int, 1U> &rowPtr, const ::coder::array<int, 1U> &colInd, ::coder::
    array<double, 1U> &vals, const ::coder::array<int, 1U> &nnzPr, const ::coder::
    array<int, 1U> &nrange, int iend, WlsDataStruct *wls, ::coder::array<
    boolean_T, 1U> &rdTags)
  {
    double xsLocal[1024];
    double Ns[512];
    double us[512];
    double nrm0_idx_0;
    double nrm0_idx_1;
    int eids[512];
    int gIDs[512];
    int i;
    int rdCounts;
    boolean_T exitg1;
    boolean_T isSphSurf;
    rdCounts = 0;

    //  kernel for surface computation
    wls->wlsWgts.params_pointwise.size[1] = 1;
    wls->wlsWgts.params_pointwise.size[0] = 512;

    //  check if spherical surface
    isSphSurf = (surfType == 1);

    //  loop begins
    i = 0;
    exitg1 = false;
    while ((!exitg1) && (i + 1 <= iend)) {
      double a;
      double ori_idx_0;
      double ori_idx_1;
      int b_xsLocal_tmp;
      int k;
      int nPoints;
      int us_tmp;
      int xsLocal_tmp;

      //  get nodal ID
      nPoints = stcls[(stcls.size(1) + stcls.size(1) * i) - 1] - 1;

      //  fetch data to local buffer and localize the coordinates
      for (int j{0}; j <= nPoints; j++) {
        gIDs[j] = stcls[j + stcls.size(1) * i];
      }

      //  fetch coordinates
      for (k = 0; k <= nPoints; k++) {
        xsLocal_tmp = k << 1;
        b_xsLocal_tmp = gIDs[k];
        xsLocal[xsLocal_tmp] = xs[xs.size(1) * (b_xsLocal_tmp - 1)];
        xsLocal[xsLocal_tmp + 1] = xs[xs.size(1) * (b_xsLocal_tmp - 1) + 1];
      }

      //  origin
      ori_idx_0 = xsLocal[0];
      ori_idx_1 = xsLocal[1];

      //  localize
      for (k = 0; k <= nPoints; k++) {
        xsLocal_tmp = k << 1;
        xsLocal[xsLocal_tmp] -= ori_idx_0;
        xsLocal[xsLocal_tmp + 1] -= ori_idx_1;
      }

      //  get normal direction
      if ((nrms.size(0) != 0) && (nrms.size(1) != 0)) {
        nrm0_idx_0 = nrms[nrms.size(1) * (nrange[i] - 1)];
        nrm0_idx_1 = nrms[nrms.size(1) * (nrange[i] - 1) + 1];
      } else if (isSphSurf) {
        //  spherical or circle
        nrm0_idx_0 = xs[xs.size(1) * (nrange[i] - 1)];
        nrm0_idx_1 = xs[xs.size(1) * (nrange[i] - 1) + 1];
      }

      //  compute tagent
      if (nPoints >= 0) {
        std::memset(&us[0], 0, (nPoints + 1) * sizeof(double));
      }

      //  matrix-matrix, using mem-efficient loop
      for (int ii{0}; ii <= nPoints; ii++) {
        us_tmp = ii << 1;
        us[ii] = (us[ii] + xsLocal[us_tmp] * -nrm0_idx_1) + nrm0_idx_0 *
          xsLocal[us_tmp + 1];
      }

      //  compute scaling to the WLS weights by inner product of normals
      if ((nrms.size(0) != 0) && (nrms.size(1) != 0)) {
        for (int j{2}; j <= nPoints + 1; j++) {
          us_tmp = stcls[(j + stcls.size(1) * i) - 1] - 1;
          a = nrm0_idx_0 * nrms[nrms.size(1) * us_tmp] + nrm0_idx_1 *
            nrms[nrms.size(1) * us_tmp + 1];
          if (a < 0.0) {
            a = 0.0;
          }

          wls->wlsWgts.params_pointwise.data[j - 1] = a;
        }
      } else if (isSphSurf) {
        //  spherical or circle
        for (int j{2}; j <= nPoints + 1; j++) {
          int idx;
          idx = stcls[(j + stcls.size(1) * i) - 1] - 1;
          a = nrm0_idx_0 * xs[xs.size(1) * idx] + nrm0_idx_1 * xs[xs.size(1) *
            idx + 1];
          if (a < 0.0) {
            a = 0.0;
          }

          wls->wlsWgts.params_pointwise.data[j - 1] = a;
        }
      }

      wls->wlsWgts.params_pointwise.data[0] = 1.0;

      //  compute wls
      b_wls_init(&wls->wlsObj, us, 512, wls->wlsWgts.name.data,
                 wls->wlsWgts.params_pointwise.data,
                 wls->wlsWgts.params_pointwise.size, degree, wls->interp0,
                 wls->useDag, nPoints + 1);
      if (wls->wlsObj.rank < 0) {
        rdCounts = -nrange[i];

        //  record error node and return
        exitg1 = true;
      } else {
        if (wls->wlsObj.rank != wls->wlsObj.ncols) {
          rdCounts++;

          //  increment counter for RD nodes
          rdTags[nrange[i] - 1] = true;
        } else {
          int nCells;

          //  get centers and shape function Ns at centers
          us_tmp = n2cPtr[nrange[i] - 1];
          nCells = (n2cPtr[nrange[i]] - us_tmp) - 1;
          for (int j{0}; j <= nCells; j++) {
            double d;
            double d1;
            int eid;
            int npts;
            eid = n2cList[(us_tmp + j) - 1] - 1;

            eids[j] = eid + 1;
            for (npts = conn.size(1); conn[(npts + conn.size(1) * eid) - 1] <= 0;
                 npts--) {
            }

            Ns[j] = 1.0 / static_cast<double>(npts);

            //  compute center
            xsLocal_tmp = j << 1;
            xsLocal[xsLocal_tmp] = xs[xs.size(1) * (conn[conn.size(1) * eid] - 1)];
            xsLocal[xsLocal_tmp + 1] = xs[xs.size(1) * (conn[conn.size(1) * eid]
              - 1) + 1];

            //  C++
            for (int ii{2}; ii <= npts; ii++) {
              b_xsLocal_tmp = conn[(ii + conn.size(1) * eid) - 1] - 1;
              xsLocal[xsLocal_tmp] += xs[xs.size(1) * b_xsLocal_tmp];
              xsLocal[xsLocal_tmp + 1] += xs[xs.size(1) * b_xsLocal_tmp + 1];
            }

            d = Ns[j];
            d1 = xsLocal[xsLocal_tmp] * d;
            d *= xsLocal[xsLocal_tmp + 1];
            if (isSphSurf) {
              a = std::sqrt(d1 * d1 + d * d);
              d1 /= a;
              d /= a;
            }

            //  localize
            d1 -= ori_idx_0;
            xsLocal[xsLocal_tmp] = d1;
            d -= ori_idx_1;
            xsLocal[xsLocal_tmp + 1] = d;
          }

          if (nCells >= 0) {
            std::memset(&us[0], 0, (nCells + 1) * sizeof(double));
          }

          //  matrix-matrix, using mem-efficient loop
          for (int ii{0}; ii <= nCells; ii++) {
            us_tmp = ii << 1;
            us[ii] = (us[ii] + xsLocal[us_tmp] * -nrm0_idx_1) + nrm0_idx_0 *
              xsLocal[us_tmp + 1];
          }

          wls_func(&wls->wlsObj, us, nCells + 1, wls->coeffs);

          //  Compute local OSUS coefficients
          for (k = 0; k <= nPoints; k++) {
            for (int j{0}; j <= nCells; j++) {
              wls->coeffs[j + wls->coeffs.size(1) * k] = wls->coeffs[j +
                wls->coeffs.size(1) * k] * Ns[j];
            }
          }

          //  substract linear interp
          for (int j{0}; j <= nCells; j++) {
            wls->coeffs[j] = wls->coeffs[j] - Ns[j];
          }

          //  kernel for updating crs matrix
          for (int j{0}; j <= nCells; j++) {
            int rowIdx;
            rowIdx = eids[j] - 1;
            for (int ii{0}; ii <= nPoints; ii++) {
              int b_i;
              us_tmp = gIDs[ii];
              b_i = rowPtr[rowIdx] - 1;
              b_xsLocal_tmp = (rowPtr[rowIdx] + nnzPr[eids[j] - 1]) - 1;
              int exitg2;
              do {
                exitg2 = 0;
                if (b_i + 1 <= b_xsLocal_tmp) {
                  if (colInd[b_i] == us_tmp) {
                    k = b_i;
                    exitg2 = 1;
                  } else {
                    b_i++;
                  }
                } else {
                  k = -1;

                  m2cErrMsgIdAndTxt("add_crs:missingIndex",
                                    "could not find column index %d", us_tmp);
                  exitg2 = 1;
                }
              } while (exitg2 == 0);

              vals[k] = vals[k] + wls->coeffs[j + wls->coeffs.size(1) * ii];
            }
          }
        }

        i++;
      }
    }

    return rdCounts;
  }

  static void call_metis_mesh(int n, ::coder::array<int, 1U> &eptr, ::coder::
    array<int, 1U> &eind, int nParts, ::coder::array<int, 1U> &nparts, ::coder::
    array<int, 1U> &cparts)
  {
    int loop_ub;
    int m;
    int status;

    //  this funtion calls METIS
    m = eptr.size(0) - 1;

    //  num. of cells
    nparts.set_size(n);
    for (int i{0}; i < n; i++) {
      nparts[i] = 1;
    }

    cparts.set_size(eptr.size(0) - 1);
    loop_ub = eptr.size(0);
    for (int i{0}; i <= loop_ub - 2; i++) {
      cparts[i] = 1;
    }

    int opts[40];
    int ncuts;

#ifdef OMP4M_HAS_METIS

    //  if we have METIS
    METIS_SetDefaultOptions(&opts[0]);
    opts[17] = 1;

    //  1-based
    status = METIS_PartMeshNodal(&m, &n, &(eptr.data())[0], &(eind.data())[0],
      NULL, NULL, &nParts, NULL, &opts[0], &ncuts, &(cparts.data())[0],
      &(nparts.data())[0]);
    if (status != 1) {
      m2cErrMsgIdAndTxt("metis:MetisError", "METIS returned error %d.", status);
    }

#else

    //  if we do NOT have METIS
    m2cWarnMsgIdAndTxt("metis:missingMetis",
                       "METIS is not available. All nodes are assigned to partition 1.");

#endif

  }

  static inline
  void compute_area(const ::coder::array<double, 2U> &xs, const ::coder::
    array<int, 2U> &conn, ::coder::array<double, 1U> &A, int nThreads)
  {
    //  kernel for computing areas
    if (nThreads <= 1) {
      compute_area_kernel(conn.size(0), xs.size(1), xs, conn, A);
    } else {
      int m2cTryBlkErrCode;

      m2cTryBlkErrCode = (0);

#pragma omp parallel default(shared) num_threads(nThreads)

      try { //try
        compute_area_kernel(conn.size(0), xs.size(1), xs, conn, A, nThreads);

      } catch (const std::runtime_error& m2cExc) {
        m2cTryBlkErrCode = (1);
        m2cPrintError("runtime_error %s\n", m2cExc.what());

      }

      catch (const std::logic_error& m2cExc)
      {
        m2cTryBlkErrCode = (1);
        m2cPrintError("logic_error %s\n", m2cExc.what());

      }

      catch (const std::exception& m2cExc)
      {
        m2cTryBlkErrCode = (1);
        m2cPrintError("exception %s\n", m2cExc.what());

      }

      catch (...)
      {
        m2cTryBlkErrCode = (1);
        m2cPrintError("Unknown error detected from C++ exceptions\n");

      } //try

      if ((int)m2cTryBlkErrCode != 0) {
        throw std::runtime_error("omp4m:runtimeErrorInThread");
      }
    }
  }

  static void compute_area_kernel(int m, int dim, const ::coder::array<double,
    2U> &xs, const ::coder::array<int, 2U> &conn, ::coder::array<double, 1U> &A)
  {
    //  serial
    if (dim == 2) {
      for (int e{0}; e < m; e++) {
        double a;
        double a1_idx_0_tmp;
        double a1_idx_1_tmp;
        double a_tmp;
        double b_a_tmp;
        int npts;
        for (npts = conn.size(1); conn[(npts + conn.size(1) * e) - 1] <= 0; npts
             --) {
        }

        a_tmp = xs[xs.size(1) * (conn[conn.size(1) * e] - 1) + 1];
        a1_idx_0_tmp = xs[xs.size(1) * (conn[conn.size(1) * e] - 1)];
        a1_idx_1_tmp = xs[xs.size(1) * (conn[conn.size(1) * e + 2] - 1)] -
          a1_idx_0_tmp;
        b_a_tmp = xs[xs.size(1) * (conn[conn.size(1) * e + 2] - 1) + 1] - a_tmp;
        a = 0.5 * std::abs((xs[xs.size(1) * (conn[conn.size(1) * e + 1] - 1)] -
                            a1_idx_0_tmp) * b_a_tmp - (xs[xs.size(1) *
          (conn[conn.size(1) * e + 1] - 1) + 1] - a_tmp) * a1_idx_1_tmp);
        if (npts > 3) {
          a += 0.5 * std::abs(a1_idx_1_tmp * (xs[xs.size(1) * (conn[conn.size(1)
            * e + 3] - 1) + 1] - a_tmp) - b_a_tmp * (xs[xs.size(1) *
            (conn[conn.size(1) * e + 3] - 1)] - a1_idx_0_tmp));
        }

        A[e] = a;
      }
    } else {
      for (int e{0}; e < m; e++) {
        double a;
        double a1_idx_0_tmp;
        double a1_idx_1_tmp;
        double a1_idx_2_tmp;
        double a2_idx_0_tmp;
        double a2_idx_1_tmp;
        double a2_idx_2_tmp;
        double a_tmp;
        double b_a_tmp;
        double d;
        int npts;
        for (npts = conn.size(1); conn[(npts + conn.size(1) * e) - 1] <= 0; npts
             --) {
        }

        double a1_idx_0;
        double a1_idx_1;
        double a1_idx_2;
        a1_idx_0_tmp = xs[xs.size(1) * (conn[conn.size(1) * e] - 1)];
        a1_idx_0 = xs[xs.size(1) * (conn[conn.size(1) * e + 1] - 1)] -
          a1_idx_0_tmp;
        a1_idx_1_tmp = xs[xs.size(1) * (conn[conn.size(1) * e] - 1) + 1];
        a1_idx_1 = xs[xs.size(1) * (conn[conn.size(1) * e + 1] - 1) + 1] -
          a1_idx_1_tmp;
        a1_idx_2_tmp = xs[xs.size(1) * (conn[conn.size(1) * e] - 1) + 2];
        a1_idx_2 = xs[xs.size(1) * (conn[conn.size(1) * e + 1] - 1) + 2] -
          a1_idx_2_tmp;
        a2_idx_0_tmp = xs[xs.size(1) * (conn[conn.size(1) * e + 2] - 1)] -
          a1_idx_0_tmp;
        a2_idx_1_tmp = xs[xs.size(1) * (conn[conn.size(1) * e + 2] - 1) + 1] -
          a1_idx_1_tmp;
        a2_idx_2_tmp = xs[xs.size(1) * (conn[conn.size(1) * e + 2] - 1) + 2] -
          a1_idx_2_tmp;

        //  cross
        d = a1_idx_1 * a2_idx_2_tmp - a2_idx_1_tmp * a1_idx_2;
        b_a_tmp = a2_idx_0_tmp * a1_idx_2 - a1_idx_0 * a2_idx_2_tmp;
        a_tmp = a1_idx_0 * a2_idx_1_tmp - a2_idx_0_tmp * a1_idx_1;
        a = 0.5 * std::sqrt((d * d + b_a_tmp * b_a_tmp) + a_tmp * a_tmp);
        if (npts > 3) {
          double a2_idx_0;
          double a2_idx_1;
          double a2_idx_2;
          a2_idx_0 = xs[xs.size(1) * (conn[conn.size(1) * e + 3] - 1)] -
            a1_idx_0_tmp;
          a2_idx_1 = xs[xs.size(1) * (conn[conn.size(1) * e + 3] - 1) + 1] -
            a1_idx_1_tmp;
          a2_idx_2 = xs[xs.size(1) * (conn[conn.size(1) * e + 3] - 1) + 2] -
            a1_idx_2_tmp;

          //  cross
          d = a2_idx_1_tmp * a2_idx_2 - a2_idx_1 * a2_idx_2_tmp;
          b_a_tmp = a2_idx_0 * a2_idx_2_tmp - a2_idx_0_tmp * a2_idx_2;
          a_tmp = a2_idx_0_tmp * a2_idx_1 - a2_idx_0 * a2_idx_1_tmp;
          a += 0.5 * std::sqrt((d * d + b_a_tmp * b_a_tmp) + a_tmp * a_tmp);
        }

        A[e] = a;
      }
    }
  }

  static void compute_area_kernel(int m, int dim, const ::coder::array<double,
    2U> &xs, const ::coder::array<int, 2U> &conn, ::coder::array<double, 1U> &A,
    int nThreads)
  {
    int iend;
    int istart;
    int threadID;
    int u1;

    //  Obtains starting and ending indices of local chunk for current thread
    threadID = 0;

#ifdef _OPENMP

    threadID = omp_get_thread_num();

#endif //_OPENMP

    if (threadID >= nThreads) {
      //  The excess threads would not be assigned any task
      istart = 1;
      iend = 0;
    } else {
      int b_remainder;
      int chunk;
      chunk = m / nThreads;
      b_remainder = m - nThreads * chunk;
      u1 = threadID;
      if (b_remainder <= threadID) {
        u1 = b_remainder;
      }

      istart = (threadID * chunk + u1) + 1;
      iend = ((istart + chunk) + (threadID < b_remainder)) - 1;
    }

    //  Note that we do a very rough estimation
    if (dim == 2) {
      for (int e{istart}; e <= iend; e++) {
        double a;
        double a1_idx_0_tmp;
        double a1_idx_1_tmp;
        double a_tmp;
        double c_a_tmp;
        int npts;
        int u0;
        for (npts = conn.size(1); conn[(npts + conn.size(1) * (e - 1)) - 1] <= 0;
             npts--) {
        }

        int b_a_tmp;
        u0 = conn[conn.size(1) * (e - 1)] - 1;
        u1 = conn[conn.size(1) * (e - 1) + 1] - 1;
        a_tmp = xs[xs.size(1) * u0 + 1];
        b_a_tmp = conn[conn.size(1) * (e - 1) + 2] - 1;
        a1_idx_0_tmp = xs[xs.size(1) * u0];
        a1_idx_1_tmp = xs[xs.size(1) * b_a_tmp] - a1_idx_0_tmp;
        c_a_tmp = xs[xs.size(1) * b_a_tmp + 1] - a_tmp;
        a = 0.5 * std::abs((xs[xs.size(1) * u1] - a1_idx_0_tmp) * c_a_tmp -
                           (xs[xs.size(1) * u1 + 1] - a_tmp) * a1_idx_1_tmp);
        if (npts > 3) {
          u0 = conn[conn.size(1) * (e - 1) + 3] - 1;
          a += 0.5 * std::abs(a1_idx_1_tmp * (xs[xs.size(1) * u0 + 1] - a_tmp) -
                              c_a_tmp * (xs[xs.size(1) * u0] - a1_idx_0_tmp));
        }

        A[e - 1] = a;
      }
    } else {
      for (int e{istart}; e <= iend; e++) {
        double a;
        double a1_idx_0_tmp;
        double a1_idx_1_tmp;
        double a1_idx_2_tmp;
        double a2_idx_0_tmp;
        double a2_idx_1_tmp;
        double a2_idx_2_tmp;
        double a_tmp;
        double c_a_tmp;
        double d;
        int npts;
        int u0;
        for (npts = conn.size(1); conn[(npts + conn.size(1) * (e - 1)) - 1] <= 0;
             npts--) {
        }

        double a1_idx_0;
        double a1_idx_1;
        double a1_idx_2;
        u0 = conn[conn.size(1) * (e - 1) + 1] - 1;
        u1 = conn[conn.size(1) * (e - 1)] - 1;
        a1_idx_0_tmp = xs[xs.size(1) * u1];
        a1_idx_0 = xs[xs.size(1) * u0] - a1_idx_0_tmp;
        a1_idx_1_tmp = xs[xs.size(1) * u1 + 1];
        a1_idx_1 = xs[xs.size(1) * u0 + 1] - a1_idx_1_tmp;
        a1_idx_2_tmp = xs[xs.size(1) * u1 + 2];
        a1_idx_2 = xs[xs.size(1) * u0 + 2] - a1_idx_2_tmp;
        u0 = conn[conn.size(1) * (e - 1) + 2] - 1;
        a2_idx_0_tmp = xs[xs.size(1) * u0] - a1_idx_0_tmp;
        a2_idx_1_tmp = xs[xs.size(1) * u0 + 1] - a1_idx_1_tmp;
        a2_idx_2_tmp = xs[xs.size(1) * u0 + 2] - a1_idx_2_tmp;

        //  cross
        d = a1_idx_1 * a2_idx_2_tmp - a2_idx_1_tmp * a1_idx_2;
        c_a_tmp = a2_idx_0_tmp * a1_idx_2 - a1_idx_0 * a2_idx_2_tmp;
        a_tmp = a1_idx_0 * a2_idx_1_tmp - a2_idx_0_tmp * a1_idx_1;
        a = 0.5 * std::sqrt((d * d + c_a_tmp * c_a_tmp) + a_tmp * a_tmp);
        if (npts > 3) {
          double a2_idx_0;
          double a2_idx_1;
          double a2_idx_2;
          u0 = conn[conn.size(1) * (e - 1) + 3] - 1;
          a2_idx_0 = xs[xs.size(1) * u0] - a1_idx_0_tmp;
          a2_idx_1 = xs[xs.size(1) * u0 + 1] - a1_idx_1_tmp;
          a2_idx_2 = xs[xs.size(1) * u0 + 2] - a1_idx_2_tmp;

          //  cross
          d = a2_idx_1_tmp * a2_idx_2 - a2_idx_1 * a2_idx_2_tmp;
          c_a_tmp = a2_idx_0 * a2_idx_2_tmp - a2_idx_0_tmp * a2_idx_2;
          a_tmp = a2_idx_0_tmp * a2_idx_1 - a2_idx_0 * a2_idx_1_tmp;
          a += 0.5 * std::sqrt((d * d + c_a_tmp * c_a_tmp) + a_tmp * a_tmp);
        }

        A[e - 1] = a;
      }
    }
  }

  static void compute_beta_kernel(int n, double epsBeta, double hGlobal, const ::
    coder::array<double, 2U> &dfGlobal, const ::coder::array<double, 2U>
    &alphaCell, const ::coder::array<double, 1U> &mesh_cellWeights, const ::
    coder::array<int, 1U> &mesh_n2cPtr, const ::coder::array<int, 1U>
    &mesh_n2cList, int nRhs, ::coder::array<double, 2U> &beta)
  {
    double epsh2;
    int iend;
    int istart;
    int nthreads;
    int u1;

#pragma omp single

    { //single
      beta.set_size(n, nRhs);

    } //single

    //  Obtains starting and ending indices of local chunk for current thread
    nthreads = 1;

#ifdef _OPENMP

    nthreads = omp_get_num_threads();

#endif //_OPENMP

    if (nthreads == 1) {
      istart = 0;
      iend = n;
    } else {
      int b_remainder;
      int chunk;
      int threadID;

      //  Gets the thread number of the thread, within the team, making this call.
      threadID = 0;

#ifdef _OPENMP

      threadID = omp_get_thread_num();

#endif //_OPENMP

      chunk = n / nthreads;
      b_remainder = n - nthreads * chunk;
      u1 = threadID;
      if (b_remainder <= threadID) {
        u1 = b_remainder;
      }

      istart = threadID * chunk + u1;
      iend = (istart + chunk) + (threadID < b_remainder);
    }

    epsh2 = epsBeta * hGlobal * hGlobal;

    //  C++
    for (int i{istart + 1}; i <= iend; i++) {
      for (int k{0}; k < nRhs; k++) {
        double aBar;
        double aSum;
        double aTop;
        double wBar;
        int b_i;
        int u0;
        aSum = 0.0;
        wBar = 0.0;
        u1 = mesh_n2cPtr[i - 1];
        b_i = mesh_n2cPtr[i] - 1;
        for (int j{u1}; j <= b_i; j++) {
          double aSum_tmp;
          u0 = mesh_n2cList[j - 1] - 1;
          aSum_tmp = mesh_cellWeights[u0];
          aSum += aSum_tmp * alphaCell[k + alphaCell.size(1) * u0];
          wBar += aSum_tmp;
        }

        //  sum(w_e*a_e)/sum(w_e)
        aBar = aSum / wBar;
        aTop = 0.0;
        for (int j{u1}; j <= b_i; j++) {
          u0 = mesh_n2cList[j - 1] - 1;
          aTop += mesh_cellWeights[u0] * std::abs(alphaCell[k + alphaCell.size(1)
            * u0] - aBar);
        }

        beta[k + beta.size(1) * (i - 1)] = aTop / ((std::abs(aSum) + wBar *
          (dfGlobal[k] * epsh2)) + 2.2250738585072014E-308);
      }
    }
  }

  static inline
  void compute_body_h(const ::coder::array<double, 2U> &xs, const ::coder::
    array<int, 2U> &conn, ::coder::array<double, 1U> &h, int nThreads)
  {
    //  kernel for computing body cell sizes
    if (nThreads <= 1) {
      compute_body_h_kernel(conn.size(0), xs.size(1), xs, conn, h);
    } else {
      int m2cTryBlkErrCode;

      m2cTryBlkErrCode = (0);

#pragma omp parallel default(shared) num_threads(nThreads)

      try { //try
        compute_body_h_kernel(conn.size(0), xs.size(1), xs, conn, h, nThreads);

      } catch (const std::runtime_error& m2cExc) {
        m2cTryBlkErrCode = (1);
        m2cPrintError("runtime_error %s\n", m2cExc.what());

      }

      catch (const std::logic_error& m2cExc)
      {
        m2cTryBlkErrCode = (1);
        m2cPrintError("logic_error %s\n", m2cExc.what());

      }

      catch (const std::exception& m2cExc)
      {
        m2cTryBlkErrCode = (1);
        m2cPrintError("exception %s\n", m2cExc.what());

      }

      catch (...)
      {
        m2cTryBlkErrCode = (1);
        m2cPrintError("Unknown error detected from C++ exceptions\n");

      } //try

      if ((int)m2cTryBlkErrCode != 0) {
        throw std::runtime_error("omp4m:runtimeErrorInThread");
      }
    }
  }

  static void compute_body_h_kernel(int m, int dim, const ::coder::array<double,
    2U> &xs, const ::coder::array<int, 2U> &conn, ::coder::array<double, 1U> &h,
    int nThreads)
  {
    int iend;
    int istart;
    int threadID;
    int u1;

    //  Obtains starting and ending indices of local chunk for current thread
    threadID = 0;

#ifdef _OPENMP

    threadID = omp_get_thread_num();

#endif //_OPENMP

    if (threadID >= nThreads) {
      //  The excess threads would not be assigned any task
      istart = 1;
      iend = 0;
    } else {
      int b_remainder;
      int chunk;
      chunk = m / nThreads;
      b_remainder = m - nThreads * chunk;
      u1 = threadID;
      if (b_remainder <= threadID) {
        u1 = b_remainder;
      }

      istart = (threadID * chunk + u1) + 1;
      iend = ((istart + chunk) + (threadID < b_remainder)) - 1;
    }

    switch (dim) {
     case 1:
      {
        for (int e{istart}; e <= iend; e++) {
          int u0;
          u0 = conn[conn.size(1) * (e - 1) + 1] - 1;
          u1 = conn[conn.size(1) * (e - 1)] - 1;
          h[e - 1] = std::abs(xs[u0 % xs.size(0) * xs.size(1) + u0 / xs.size(0)]
                              - xs[u1 % xs.size(0) * xs.size(1) + u1 / xs.size(0)]);
        }
      }
      break;

     case 2:
      {
        for (int e{istart}; e <= iend; e++) {
          int npts;
          for (npts = conn.size(1); conn[(npts + conn.size(1) * (e - 1)) - 1] <=
               0; npts--) {
          }

          if (npts == 3) {
            double v_idx_0;
            double v_idx_1;
            double v_idx_2;
            double v_idx_3;
            double v_idx_4;
            double v_idx_5;
            int b_v_idx_2_tmp;
            int u0;
            u0 = conn[conn.size(1) * (e - 1) + 1] - 1;
            u1 = conn[conn.size(1) * (e - 1)] - 1;
            v_idx_0 = xs[xs.size(1) * u0] - xs[xs.size(1) * u1];
            v_idx_1 = xs[xs.size(1) * u0 + 1] - xs[xs.size(1) * u1 + 1];
            b_v_idx_2_tmp = conn[conn.size(1) * (e - 1) + 2] - 1;
            v_idx_2 = xs[xs.size(1) * b_v_idx_2_tmp] - xs[xs.size(1) * u1];
            v_idx_3 = xs[xs.size(1) * b_v_idx_2_tmp + 1] - xs[xs.size(1) * u1 +
              1];
            v_idx_4 = xs[xs.size(1) * b_v_idx_2_tmp] - xs[xs.size(1) * u0];
            v_idx_5 = xs[xs.size(1) * b_v_idx_2_tmp + 1] - xs[xs.size(1) * u0 +
              1];
            h[e - 1] = ((std::sqrt(v_idx_0 * v_idx_0 + v_idx_1 * v_idx_1) + std::
                         sqrt(v_idx_2 * v_idx_2 + v_idx_3 * v_idx_3)) + std::
                        sqrt(v_idx_4 * v_idx_4 + v_idx_5 * v_idx_5)) / 3.0;
          } else {
            double b_v_idx_0_tmp;
            double b_v_idx_1_tmp;
            double v_idx_0;
            double v_idx_0_tmp;
            double v_idx_1;
            double v_idx_1_tmp;
            double v_idx_2;
            double v_idx_2_tmp;
            double v_idx_3;
            double v_idx_3_tmp;
            double v_idx_4;
            double v_idx_5;
            double v_idx_6;
            double v_idx_7;
            int b_v_idx_2_tmp;
            int u0;
            u0 = conn[conn.size(1) * (e - 1) + 1] - 1;
            u1 = conn[conn.size(1) * (e - 1)] - 1;
            v_idx_0_tmp = xs[xs.size(1) * u1];
            b_v_idx_0_tmp = xs[xs.size(1) * u0];
            v_idx_0 = b_v_idx_0_tmp - v_idx_0_tmp;
            b_v_idx_1_tmp = xs[xs.size(1) * u1 + 1];
            v_idx_1_tmp = xs[xs.size(1) * u0 + 1];
            v_idx_1 = v_idx_1_tmp - b_v_idx_1_tmp;
            b_v_idx_2_tmp = conn[conn.size(1) * (e - 1) + 3] - 1;
            v_idx_2_tmp = xs[xs.size(1) * b_v_idx_2_tmp];
            v_idx_2 = v_idx_2_tmp - v_idx_0_tmp;
            v_idx_3_tmp = xs[xs.size(1) * b_v_idx_2_tmp + 1];
            v_idx_3 = v_idx_3_tmp - b_v_idx_1_tmp;
            u0 = conn[conn.size(1) * (e - 1) + 2] - 1;
            v_idx_0_tmp = xs[xs.size(1) * u0];
            v_idx_4 = v_idx_0_tmp - b_v_idx_0_tmp;
            b_v_idx_1_tmp = xs[xs.size(1) * u0 + 1];
            v_idx_5 = b_v_idx_1_tmp - v_idx_1_tmp;
            v_idx_6 = v_idx_2_tmp - v_idx_0_tmp;
            v_idx_7 = v_idx_3_tmp - b_v_idx_1_tmp;
            h[e - 1] = (((std::sqrt(v_idx_0 * v_idx_0 + v_idx_1 * v_idx_1) + std::
                          sqrt(v_idx_2 * v_idx_2 + v_idx_3 * v_idx_3)) + std::
                         sqrt(v_idx_4 * v_idx_4 + v_idx_5 * v_idx_5)) + std::
                        sqrt(v_idx_6 * v_idx_6 + v_idx_7 * v_idx_7)) * 0.25;
          }
        }
      }
      break;

     case 3:
      {
        //  only tet meshes are supported
        for (int e{istart}; e <= iend; e++) {
          double a;
          double a_tmp;
          double b;
          double b_a_tmp;
          double b_c_tmp;
          double b_tmp;
          double b_v_idx_0_tmp;
          double b_v_idx_1_tmp;
          double c;
          double c_tmp;
          double h0;
          double v_idx_0_tmp;
          double v_idx_1_tmp;
          double v_idx_2_tmp;
          double v_idx_3_tmp;
          int u0;
          u0 = conn[conn.size(1) * (e - 1)];
          u1 = conn[conn.size(1) * (e - 1) + 1];
          v_idx_1_tmp = xs[xs.size(1) * (u1 - 1)];
          v_idx_0_tmp = xs[xs.size(1) * (u0 - 1)];
          a = v_idx_0_tmp - v_idx_1_tmp;
          v_idx_2_tmp = xs[xs.size(1) * (u1 - 1) + 1];
          b_v_idx_1_tmp = xs[xs.size(1) * (u0 - 1) + 1];
          b = b_v_idx_1_tmp - v_idx_2_tmp;
          v_idx_3_tmp = xs[xs.size(1) * (u1 - 1) + 2];
          c_tmp = xs[xs.size(1) * (u0 - 1) + 2];
          c = c_tmp - v_idx_3_tmp;
          h0 = std::sqrt((a * a + b * b) + c * c);
          u0 = conn[conn.size(1) * (e - 1) + 2];
          a_tmp = xs[xs.size(1) * (u0 - 1)];
          a = v_idx_1_tmp - a_tmp;
          b_tmp = xs[xs.size(1) * (u0 - 1) + 1];
          b = v_idx_2_tmp - b_tmp;
          b_c_tmp = xs[xs.size(1) * (u0 - 1) + 2];
          c = v_idx_3_tmp - b_c_tmp;
          h0 += std::sqrt((a * a + b * b) + c * c);
          a = a_tmp - v_idx_0_tmp;
          b = b_tmp - b_v_idx_1_tmp;
          c = b_c_tmp - c_tmp;
          h0 += std::sqrt((a * a + b * b) + c * c);
          u0 = conn[conn.size(1) * (e - 1) + 3];
          b_a_tmp = xs[xs.size(1) * (u0 - 1)];
          a = v_idx_0_tmp - b_a_tmp;
          b_v_idx_0_tmp = xs[xs.size(1) * (u0 - 1) + 1];
          b = b_v_idx_1_tmp - b_v_idx_0_tmp;
          v_idx_0_tmp = xs[xs.size(1) * (u0 - 1) + 2];
          c = c_tmp - v_idx_0_tmp;
          h0 += std::sqrt((a * a + b * b) + c * c);
          a = v_idx_1_tmp - b_a_tmp;
          b = v_idx_2_tmp - b_v_idx_0_tmp;
          c = v_idx_3_tmp - v_idx_0_tmp;
          h0 += std::sqrt((a * a + b * b) + c * c);
          a = a_tmp - b_a_tmp;
          b = b_tmp - b_v_idx_0_tmp;
          c = b_c_tmp - v_idx_0_tmp;
          h0 += std::sqrt((a * a + b * b) + c * c);
          h[e - 1] = h0 / 6.0;
        }
      }
      break;
    }
  }

  static void compute_body_h_kernel(int m, int dim, const ::coder::array<double,
    2U> &xs, const ::coder::array<int, 2U> &conn, ::coder::array<double, 1U> &h)
  {
    //  serial
    switch (dim) {
     case 1:
      {
        for (int e{0}; e < m; e++) {
          int i;
          int i1;
          i = conn[conn.size(1) * e + 1] - 1;
          i1 = conn[conn.size(1) * e] - 1;
          h[e] = std::abs(xs[i % xs.size(0) * xs.size(1) + i / xs.size(0)] -
                          xs[i1 % xs.size(0) * xs.size(1) + i1 / xs.size(0)]);
        }
      }
      break;

     case 2:
      {
        for (int e{0}; e < m; e++) {
          int npts;
          for (npts = conn.size(1); conn[(npts + conn.size(1) * e) - 1] <= 0;
               npts--) {
          }

          if (npts == 3) {
            double v_idx_0;
            double v_idx_1;
            double v_idx_2;
            double v_idx_3;
            double v_idx_4;
            double v_idx_5;
            v_idx_0 = xs[xs.size(1) * (conn[conn.size(1) * e + 1] - 1)] -
              xs[xs.size(1) * (conn[conn.size(1) * e] - 1)];
            v_idx_1 = xs[xs.size(1) * (conn[conn.size(1) * e + 1] - 1) + 1] -
              xs[xs.size(1) * (conn[conn.size(1) * e] - 1) + 1];
            v_idx_2 = xs[xs.size(1) * (conn[conn.size(1) * e + 2] - 1)] -
              xs[xs.size(1) * (conn[conn.size(1) * e] - 1)];
            v_idx_3 = xs[xs.size(1) * (conn[conn.size(1) * e + 2] - 1) + 1] -
              xs[xs.size(1) * (conn[conn.size(1) * e] - 1) + 1];
            v_idx_4 = xs[xs.size(1) * (conn[conn.size(1) * e + 2] - 1)] -
              xs[xs.size(1) * (conn[conn.size(1) * e + 1] - 1)];
            v_idx_5 = xs[xs.size(1) * (conn[conn.size(1) * e + 2] - 1) + 1] -
              xs[xs.size(1) * (conn[conn.size(1) * e + 1] - 1) + 1];
            h[e] = ((std::sqrt(v_idx_0 * v_idx_0 + v_idx_1 * v_idx_1) + std::
                     sqrt(v_idx_2 * v_idx_2 + v_idx_3 * v_idx_3)) + std::sqrt
                    (v_idx_4 * v_idx_4 + v_idx_5 * v_idx_5)) / 3.0;
          } else {
            double b_v_idx_0_tmp;
            double b_v_idx_1_tmp;
            double v_idx_0;
            double v_idx_0_tmp;
            double v_idx_1;
            double v_idx_1_tmp;
            double v_idx_2;
            double v_idx_2_tmp;
            double v_idx_3;
            double v_idx_3_tmp;
            double v_idx_4;
            double v_idx_5;
            double v_idx_6;
            double v_idx_7;
            v_idx_0_tmp = xs[xs.size(1) * (conn[conn.size(1) * e] - 1)];
            b_v_idx_0_tmp = xs[xs.size(1) * (conn[conn.size(1) * e + 1] - 1)];
            v_idx_0 = b_v_idx_0_tmp - v_idx_0_tmp;
            b_v_idx_1_tmp = xs[xs.size(1) * (conn[conn.size(1) * e] - 1) + 1];
            v_idx_1_tmp = xs[xs.size(1) * (conn[conn.size(1) * e + 1] - 1) + 1];
            v_idx_1 = v_idx_1_tmp - b_v_idx_1_tmp;
            v_idx_2_tmp = xs[xs.size(1) * (conn[conn.size(1) * e + 3] - 1)];
            v_idx_2 = v_idx_2_tmp - v_idx_0_tmp;
            v_idx_3_tmp = xs[xs.size(1) * (conn[conn.size(1) * e + 3] - 1) + 1];
            v_idx_3 = v_idx_3_tmp - b_v_idx_1_tmp;
            v_idx_0_tmp = xs[xs.size(1) * (conn[conn.size(1) * e + 2] - 1)];
            v_idx_4 = v_idx_0_tmp - b_v_idx_0_tmp;
            b_v_idx_1_tmp = xs[xs.size(1) * (conn[conn.size(1) * e + 2] - 1) + 1];
            v_idx_5 = b_v_idx_1_tmp - v_idx_1_tmp;
            v_idx_6 = v_idx_2_tmp - v_idx_0_tmp;
            v_idx_7 = v_idx_3_tmp - b_v_idx_1_tmp;
            h[e] = (((std::sqrt(v_idx_0 * v_idx_0 + v_idx_1 * v_idx_1) + std::
                      sqrt(v_idx_2 * v_idx_2 + v_idx_3 * v_idx_3)) + std::sqrt
                     (v_idx_4 * v_idx_4 + v_idx_5 * v_idx_5)) + std::sqrt
                    (v_idx_6 * v_idx_6 + v_idx_7 * v_idx_7)) * 0.25;
          }
        }
      }
      break;

     case 3:
      {
        //  only tet meshes are supported
        for (int e{0}; e < m; e++) {
          double a;
          double a_tmp;
          double b;
          double b_a_tmp;
          double b_c_tmp;
          double b_tmp;
          double b_v_idx_0_tmp;
          double b_v_idx_1_tmp;
          double c;
          double c_tmp;
          double h0;
          double v_idx_0_tmp;
          double v_idx_1_tmp;
          double v_idx_2_tmp;
          double v_idx_3_tmp;
          int i;
          int i1;
          i = conn[conn.size(1) * e];
          i1 = conn[conn.size(1) * e + 1];
          v_idx_1_tmp = xs[xs.size(1) * (i1 - 1)];
          v_idx_0_tmp = xs[xs.size(1) * (i - 1)];
          a = v_idx_0_tmp - v_idx_1_tmp;
          v_idx_2_tmp = xs[xs.size(1) * (i1 - 1) + 1];
          b_v_idx_1_tmp = xs[xs.size(1) * (i - 1) + 1];
          b = b_v_idx_1_tmp - v_idx_2_tmp;
          v_idx_3_tmp = xs[xs.size(1) * (i1 - 1) + 2];
          c_tmp = xs[xs.size(1) * (i - 1) + 2];
          c = c_tmp - v_idx_3_tmp;
          h0 = std::sqrt((a * a + b * b) + c * c);
          i = conn[conn.size(1) * e + 2];
          a_tmp = xs[xs.size(1) * (i - 1)];
          a = v_idx_1_tmp - a_tmp;
          b_tmp = xs[xs.size(1) * (i - 1) + 1];
          b = v_idx_2_tmp - b_tmp;
          b_c_tmp = xs[xs.size(1) * (i - 1) + 2];
          c = v_idx_3_tmp - b_c_tmp;
          h0 += std::sqrt((a * a + b * b) + c * c);
          a = a_tmp - v_idx_0_tmp;
          b = b_tmp - b_v_idx_1_tmp;
          c = b_c_tmp - c_tmp;
          h0 += std::sqrt((a * a + b * b) + c * c);
          i = conn[conn.size(1) * e + 3];
          b_a_tmp = xs[xs.size(1) * (i - 1)];
          a = v_idx_0_tmp - b_a_tmp;
          b_v_idx_0_tmp = xs[xs.size(1) * (i - 1) + 1];
          b = b_v_idx_1_tmp - b_v_idx_0_tmp;
          v_idx_0_tmp = xs[xs.size(1) * (i - 1) + 2];
          c = c_tmp - v_idx_0_tmp;
          h0 += std::sqrt((a * a + b * b) + c * c);
          a = v_idx_1_tmp - b_a_tmp;
          b = v_idx_2_tmp - b_v_idx_0_tmp;
          c = v_idx_3_tmp - v_idx_0_tmp;
          h0 += std::sqrt((a * a + b * b) + c * c);
          a = a_tmp - b_a_tmp;
          b = b_tmp - b_v_idx_0_tmp;
          c = b_c_tmp - v_idx_0_tmp;
          h0 += std::sqrt((a * a + b * b) + c * c);
          h[e] = h0 / 6.0;
        }
      }
      break;
    }
  }

  static void compute_nodal_alpha(int n, const ::coder::array<double, 2U>
    &alphaCell, const ::coder::array<int, 1U> &mesh_n2cPtr, const ::coder::array<
    int, 1U> &mesh_n2cList, int nRhs, ::coder::array<double, 2U> &alphaNode)
  {
    int iend;
    int istart;
    int nthreads;
    int u1;

#pragma omp single

    { //single
      alphaNode.set_size(n, nRhs);

    } //single

    //  implicit barrier
    nthreads = 1;

#ifdef _OPENMP

    nthreads = omp_get_num_threads();

#endif //_OPENMP

    if (nthreads == 1) {
      istart = 0;
      iend = n;
    } else {
      int b_remainder;
      int chunk;
      int threadID;

      //  Gets the thread number of the thread, within the team, making this call.
      threadID = 0;

#ifdef _OPENMP

      threadID = omp_get_thread_num();

#endif //_OPENMP

      chunk = n / nthreads;
      b_remainder = n - nthreads * chunk;
      u1 = threadID;
      if (b_remainder <= threadID) {
        u1 = b_remainder;
      }

      istart = threadID * chunk + u1;
      iend = (istart + chunk) + (threadID < b_remainder);
    }

    //  compute nodal alpha values
    for (int i{istart + 1}; i <= iend; i++) {
      for (int k{0}; k < nRhs; k++) {
        double a;
        int u0;
        a = 0.0;
        u0 = mesh_n2cPtr[i - 1];
        u1 = mesh_n2cPtr[i] - 1;
        for (int j{u0}; j <= u1; j++) {
          double d;
          d = std::abs(alphaCell[k + alphaCell.size(1) * (mesh_n2cList[j - 1] -
            1)]);
          if (a < d) {
            a = d;
          }
        }

        alphaNode[k + alphaNode.size(1) * (i - 1)] = a;
      }
    }
  }

  static void compute_surf_h(const ::coder::array<double, 2U> &xs, const ::coder::
    array<int, 2U> &conn, const ::coder::array<int, 1U> &n2nPtr, const ::coder::
    array<int, 1U> &n2nList, const ::coder::array<double, 2U> &nrms, int
    surfType, ::coder::array<double, 1U> &h, int nThreads)
  {
    ::coder::array<double, 1U> buf_;
    int loop_ub;

    //  kernel for computing the cell sizes on surfaces in 2D and 3D
    buf_.set_size(xs.size(0));
    loop_ub = xs.size(0);
    for (int i{0}; i < loop_ub; i++) {
      buf_[i] = 0.0;
    }

    if (nThreads <= 1) {
      compute_surf_h_kernel(xs.size(0), conn.size(0), xs.size(1), xs, conn,
                            n2nPtr, n2nList, nrms, surfType, h, buf_);
    } else {
      int m2cTryBlkErrCode;

      //  parallel
      m2cTryBlkErrCode = (0);

#pragma omp parallel default(shared) num_threads(nThreads)

      try { //try
        compute_surf_h_kernel(xs.size(0), conn.size(0), xs.size(1), xs, conn,
                              n2nPtr, n2nList, nrms, surfType, h, buf_, nThreads);

      } catch (const std::runtime_error& m2cExc) {
        m2cTryBlkErrCode = (1);
        m2cPrintError("runtime_error %s\n", m2cExc.what());

      }

      catch (const std::logic_error& m2cExc)
      {
        m2cTryBlkErrCode = (1);
        m2cPrintError("logic_error %s\n", m2cExc.what());

      }

      catch (const std::exception& m2cExc)
      {
        m2cTryBlkErrCode = (1);
        m2cPrintError("exception %s\n", m2cExc.what());

      }

      catch (...)
      {
        m2cTryBlkErrCode = (1);
        m2cPrintError("Unknown error detected from C++ exceptions\n");

      } //try

      if ((int)m2cTryBlkErrCode != 0) {
        throw std::runtime_error("omp4m:runtimeErrorInThread");
      }
    }

  }

  static void compute_surf_h(const ::coder::array<double, 2U> &xs, const ::coder::
    array<int, 2U> &conn, const ::coder::array<int, 1U> &n2nPtr, const ::coder::
    array<int, 1U> &n2nList, int surfType, ::coder::array<double, 1U> &h, int
    nThreads)
  {
    ::coder::array<double, 1U> buf_;
    int loop_ub;

    //  kernel for computing the cell sizes on surfaces in 2D and 3D
    buf_.set_size(xs.size(0));
    loop_ub = xs.size(0);
    for (int i{0}; i < loop_ub; i++) {
      buf_[i] = 0.0;
    }

    if (nThreads <= 1) {
      compute_surf_h_kernel(xs.size(0), conn.size(0), xs.size(1), xs, conn,
                            n2nPtr, n2nList, surfType, h, buf_);
    } else {
      int m2cTryBlkErrCode;

      //  parallel
      m2cTryBlkErrCode = (0);

#pragma omp parallel default(shared) num_threads(nThreads)

      try { //try
        compute_surf_h_kernel(xs.size(0), conn.size(0), xs.size(1), xs, conn,
                              n2nPtr, n2nList, surfType, h, buf_, nThreads);

      } catch (const std::runtime_error& m2cExc) {
        m2cTryBlkErrCode = (1);
        m2cPrintError("runtime_error %s\n", m2cExc.what());

      }

      catch (const std::logic_error& m2cExc)
      {
        m2cTryBlkErrCode = (1);
        m2cPrintError("logic_error %s\n", m2cExc.what());

      }

      catch (const std::exception& m2cExc)
      {
        m2cTryBlkErrCode = (1);
        m2cPrintError("exception %s\n", m2cExc.what());

      }

      catch (...)
      {
        m2cTryBlkErrCode = (1);
        m2cPrintError("Unknown error detected from C++ exceptions\n");

      } //try

      if ((int)m2cTryBlkErrCode != 0) {
        throw std::runtime_error("omp4m:runtimeErrorInThread");
      }
    }

  }

  static void compute_surf_h_kernel(int n, int m, int dim, const ::coder::array<
    double, 2U> &xs, const ::coder::array<int, 2U> &conn, const ::coder::array<
    int, 1U> &n2nPtr, const ::coder::array<int, 1U> &n2nList, const ::coder::
    array<double, 2U> &nrms, int surfType, ::coder::array<double, 1U> &h, ::
    coder::array<double, 1U> &buf, int nThreads)
  {
    double d1;
    double h0;
    double nrm02_idx_0;
    double nrm02_idx_1;
    double nrm03_idx_0;
    double nrm03_idx_1;
    double nrm03_idx_2;
    int b_remainder;
    int chunk;
    int iendM;
    int iendN;
    int istartM;
    int istartN;
    int threadID;
    int u1;

    //  Obtains starting and ending indices of local chunk for current thread
    threadID = 0;

#ifdef _OPENMP

    threadID = omp_get_thread_num();

#endif //_OPENMP

    if (threadID >= nThreads) {
      //  The excess threads would not be assigned any task
      istartN = 1;
      iendN = 0;
    } else {
      chunk = n / nThreads;
      b_remainder = n - nThreads * chunk;
      u1 = threadID;
      if (b_remainder <= threadID) {
        u1 = b_remainder;
      }

      istartN = (threadID * chunk + u1) + 1;
      iendN = ((istartN + chunk) + (threadID < b_remainder)) - 1;
    }

    if (dim == 3) {
      double xs3[384];
      double us3[256];

      //  nodal h
      for (int i{istartN}; i <= iendN; i++) {
        double a;
        double b_xs2_tmp;
        double d;
        double t3_idx_1;
        double t3_idx_3;
        double t3_idx_5;
        double xs2_tmp;
        double xs3_tmp;
        int b_i;
        int b_xs3_tmp;
        int nPoints;
        int u0;
        boolean_T guard1{ false };

        u0 = n2nPtr[i - 1];
        nPoints = n2nPtr[i] - u0;

        //  i is excluded
        xs2_tmp = xs[xs.size(1) * (i - 1)];
        xs3[0] = xs2_tmp;
        b_xs2_tmp = xs[xs.size(1) * (i - 1) + 1];
        xs3[1] = b_xs2_tmp;
        xs3_tmp = xs[xs.size(1) * (i - 1) + 2];
        xs3[2] = xs3_tmp;
        for (int j{0}; j < nPoints; j++) {
          b_xs3_tmp = n2nList[(u0 + j) - 1] - 1;
          u1 = 3 * (j + 1);
          xs3[u1] = xs[xs.size(1) * b_xs3_tmp];
          xs3[u1 + 1] = xs[xs.size(1) * b_xs3_tmp + 1];
          xs3[u1 + 2] = xs[xs.size(1) * b_xs3_tmp + 2];
        }

        b_i = nPoints + 1;
        for (int j{2}; j <= b_i; j++) {
          b_xs3_tmp = 3 * (j - 1);
          xs3[b_xs3_tmp] -= xs3[0];
          xs3[b_xs3_tmp + 1] -= xs3[1];
          xs3[b_xs3_tmp + 2] -= xs3[2];
        }

        if ((nrms.size(0) != 0) && (nrms.size(1) != 0)) {
          nrm03_idx_0 = nrms[nrms.size(1) * (i - 1)];
          nrm03_idx_1 = nrms[nrms.size(1) * (i - 1) + 1];
          nrm03_idx_2 = nrms[nrms.size(1) * (i - 1) + 2];
        } else if (surfType == 1) {
          //  sphere
          a = std::sqrt((xs2_tmp * xs2_tmp + b_xs2_tmp * b_xs2_tmp) + xs3_tmp *
                        xs3_tmp);
          nrm03_idx_0 = xs2_tmp / a;
          nrm03_idx_1 = b_xs2_tmp / a;
          nrm03_idx_2 = xs3_tmp / a;
        }

        d = std::abs(nrm03_idx_0);
        guard1 = false;
        if (d > std::abs(nrm03_idx_1)) {
          d1 = nrm03_idx_2;
          if (d > std::abs(nrm03_idx_2)) {
            d = -nrm03_idx_0 * nrm03_idx_1;
            b_xs2_tmp = 1.0 - nrm03_idx_1 * nrm03_idx_1;
            xs3_tmp = -nrm03_idx_1 * nrm03_idx_2;
          } else {
            guard1 = true;
          }
        } else {
          guard1 = true;
        }

        if (guard1) {
          d = 1.0 - nrm03_idx_0 * nrm03_idx_0;
          b_xs2_tmp = -nrm03_idx_0 * nrm03_idx_1;
          d1 = nrm03_idx_2;
          xs3_tmp = -nrm03_idx_0 * nrm03_idx_2;
        }

        a = std::sqrt((d * d + b_xs2_tmp * b_xs2_tmp) + xs3_tmp * xs3_tmp);
        d /= a;
        b_xs2_tmp /= a;
        xs3_tmp /= a;

        //  cross
        t3_idx_1 = nrm03_idx_1 * xs3_tmp - d1 * b_xs2_tmp;
        t3_idx_3 = d1 * d - nrm03_idx_0 * xs3_tmp;
        t3_idx_5 = nrm03_idx_0 * b_xs2_tmp - nrm03_idx_1 * d;

        //  compute tangent
        for (int j{0}; j <= nPoints; j++) {
          d1 = xs3[3 * j + 1];
          xs2_tmp = xs3[3 * j + 2];
          u0 = j << 1;
          us3[u0] = (xs3[3 * j] * d + b_xs2_tmp * d1) + xs3_tmp * xs2_tmp;
          us3[u0 + 1] = (d1 * t3_idx_1 + d1 * t3_idx_3) + t3_idx_5 * xs2_tmp;
        }

        h0 = 0.0;
        for (int j{2}; j <= b_i; j++) {
          u0 = (j - 1) << 1;
          d = us3[u0];
          d1 = us3[u0 + 1];
          h0 += std::sqrt(d * d + d1 * d1);
        }

        buf[i - 1] = h0 / static_cast<double>(nPoints);
      }
    } else {
      double xs2[256];
      double us2[128];

      //  2D
      for (int i{istartN}; i <= iendN; i++) {
        double b_xs2_tmp;
        double xs2_tmp;
        int b_i;
        int nPoints;
        int u0;
        nPoints = n2nPtr[i] - n2nPtr[i - 1];

        //  i is excluded
        xs2_tmp = xs[xs.size(1) * (i - 1)];
        xs2[0] = xs2_tmp;
        b_xs2_tmp = xs[xs.size(1) * (i - 1) + 1];
        xs2[1] = b_xs2_tmp;
        for (int j{0}; j < nPoints; j++) {
          u1 = (n2nList[i - 1] + j) - 1;
          u0 = (j + 1) << 1;
          xs2[u0] = xs[xs.size(1) * u1];
          xs2[u0 + 1] = xs[xs.size(1) * u1 + 1];
        }

        b_i = nPoints + 1;
        for (int j{2}; j <= b_i; j++) {
          u1 = (j - 1) << 1;
          xs2[u1] -= xs[0];
          xs2[u1 + 1] -= xs[1];
        }

        if ((nrms.size(0) != 0) && (nrms.size(1) != 0)) {
          nrm02_idx_0 = nrms[nrms.size(1) * (i - 1)];
          nrm02_idx_1 = nrms[nrms.size(1) * (i - 1) + 1];
        } else if (surfType == 1) {
          double a;

          //  circle
          a = std::sqrt(xs2_tmp * xs2_tmp + b_xs2_tmp * b_xs2_tmp);
          nrm02_idx_0 = xs2_tmp / a;
          nrm02_idx_1 = b_xs2_tmp / a;
        }

        //  compute tangent
        for (int j{0}; j <= nPoints; j++) {
          u0 = j << 1;
          us2[j] = xs2[u0] * -nrm02_idx_1 + xs2[u0 + 1] * nrm02_idx_0;
        }

        h0 = 0.0;
        for (int j{2}; j <= b_i; j++) {
          h0 += std::abs(us2[j - 1]);
        }

        buf[i - 1] = h0 / static_cast<double>(nPoints);
      }
    }

#pragma omp barrier

    //  Obtains starting and ending indices of local chunk for current thread
    threadID = 0;

#ifdef _OPENMP

    threadID = omp_get_thread_num();

#endif //_OPENMP

    if (threadID >= nThreads) {
      //  The excess threads would not be assigned any task
      istartM = 0;
      iendM = 0;
    } else {
      chunk = m / nThreads;
      b_remainder = m - nThreads * chunk;
      u1 = threadID;
      if (b_remainder <= threadID) {
        u1 = b_remainder;
      }

      istartM = threadID * chunk + u1;
      iendM = (istartM + chunk) + (threadID < b_remainder);
    }

    //  average to cell-based average h
    for (int e{istartM + 1}; e <= iendM; e++) {
      int npts;
      for (npts = conn.size(1); conn[(npts + conn.size(1) * (e - 1)) - 1] <= 0;
           npts--) {
      }

      h0 = 0.0;
      for (int i{0}; i < npts; i++) {
        h0 += buf[conn[i + conn.size(1) * (e - 1)] - 1];
      }

      h[e - 1] = h0 / static_cast<double>(npts);
    }
  }

  static void compute_surf_h_kernel(int n, int m, int dim, const ::coder::array<
    double, 2U> &xs, const ::coder::array<int, 2U> &conn, const ::coder::array<
    int, 1U> &n2nPtr, const ::coder::array<int, 1U> &n2nList, const ::coder::
    array<double, 2U> &nrms, int surfType, ::coder::array<double, 1U> &h, ::
    coder::array<double, 1U> &buf)
  {
    double d1;
    double d2;
    double d3;
    double h0;
    double nrm02_idx_0;
    double nrm02_idx_1;
    double nrm03_idx_0;
    double nrm03_idx_1;
    double nrm03_idx_2;

    //  serial
    if (dim == 3) {
      double xs3[384];
      double us3[256];

      //  nodal h
      for (int i{0}; i < n; i++) {
        double a;
        double d;
        double t3_idx_1;
        double t3_idx_3;
        double t3_idx_5;
        int b_i;
        int b_xs2_tmp;
        int nPoints;
        int xs2_tmp;
        boolean_T guard1{ false };

        nPoints = n2nPtr[i + 1] - n2nPtr[i];

        //  i is excluded
        xs3[0] = xs[xs.size(1) * i];
        xs3[1] = xs[xs.size(1) * i + 1];
        xs3[2] = xs[xs.size(1) * i + 2];
        for (int j{0}; j < nPoints; j++) {
          xs2_tmp = n2nList[(n2nPtr[i] + j) - 1] - 1;
          b_xs2_tmp = 3 * (j + 1);
          xs3[b_xs2_tmp] = xs[xs.size(1) * xs2_tmp];
          xs3[b_xs2_tmp + 1] = xs[xs.size(1) * xs2_tmp + 1];
          xs3[b_xs2_tmp + 2] = xs[xs.size(1) * xs2_tmp + 2];
        }

        b_i = nPoints + 1;
        for (int j{2}; j <= b_i; j++) {
          xs2_tmp = 3 * (j - 1);
          xs3[xs2_tmp] -= xs3[0];
          xs3[xs2_tmp + 1] -= xs3[1];
          xs3[xs2_tmp + 2] -= xs3[2];
        }

        if ((nrms.size(0) != 0) && (nrms.size(1) != 0)) {
          nrm03_idx_0 = nrms[nrms.size(1) * i];
          nrm03_idx_1 = nrms[nrms.size(1) * i + 1];
          nrm03_idx_2 = nrms[nrms.size(1) * i + 2];
        } else if (surfType == 1) {
          //  sphere
          nrm03_idx_0 = xs[xs.size(1) * i];
          nrm03_idx_1 = xs[xs.size(1) * i + 1];
          nrm03_idx_2 = xs[xs.size(1) * i + 2];
          a = std::sqrt((nrm03_idx_0 * nrm03_idx_0 + nrm03_idx_1 * nrm03_idx_1)
                        + nrm03_idx_2 * nrm03_idx_2);
          nrm03_idx_0 /= a;
          nrm03_idx_1 /= a;
          nrm03_idx_2 /= a;
        }

        d = std::abs(nrm03_idx_0);
        guard1 = false;
        if (d > std::abs(nrm03_idx_1)) {
          d1 = nrm03_idx_2;
          if (d > std::abs(nrm03_idx_2)) {
            d = -nrm03_idx_0 * nrm03_idx_1;
            d2 = 1.0 - nrm03_idx_1 * nrm03_idx_1;
            d3 = -nrm03_idx_1 * nrm03_idx_2;
          } else {
            guard1 = true;
          }
        } else {
          guard1 = true;
        }

        if (guard1) {
          d = 1.0 - nrm03_idx_0 * nrm03_idx_0;
          d2 = -nrm03_idx_0 * nrm03_idx_1;
          d1 = nrm03_idx_2;
          d3 = -nrm03_idx_0 * nrm03_idx_2;
        }

        a = std::sqrt((d * d + d2 * d2) + d3 * d3);
        d /= a;
        d2 /= a;
        d3 /= a;

        //  cross
        t3_idx_1 = nrm03_idx_1 * d3 - d1 * d2;
        t3_idx_3 = d1 * d - nrm03_idx_0 * d3;
        t3_idx_5 = nrm03_idx_0 * d2 - nrm03_idx_1 * d;

        //  compute tangent
        for (int j{0}; j <= nPoints; j++) {
          double d4;
          d1 = xs3[3 * j + 1];
          d4 = xs3[3 * j + 2];
          b_xs2_tmp = j << 1;
          us3[b_xs2_tmp] = (xs3[3 * j] * d + d2 * d1) + d3 * d4;
          us3[b_xs2_tmp + 1] = (d1 * t3_idx_1 + d1 * t3_idx_3) + t3_idx_5 * d4;
        }

        h0 = 0.0;
        for (int j{2}; j <= b_i; j++) {
          b_xs2_tmp = (j - 1) << 1;
          d = us3[b_xs2_tmp];
          d1 = us3[b_xs2_tmp + 1];
          h0 += std::sqrt(d * d + d1 * d1);
        }

        buf[i] = h0 / static_cast<double>(nPoints);
      }
    } else {
      double xs2[256];
      double us2[128];

      //  2D
      for (int i{0}; i < n; i++) {
        int b_i;
        int b_xs2_tmp;
        int nPoints;
        int xs2_tmp;
        nPoints = n2nPtr[i + 1] - n2nPtr[i];

        //  i is excluded
        xs2[0] = xs[xs.size(1) * i];
        xs2[1] = xs[xs.size(1) * i + 1];
        for (int j{0}; j < nPoints; j++) {
          xs2_tmp = (n2nList[i] + j) - 1;
          b_xs2_tmp = (j + 1) << 1;
          xs2[b_xs2_tmp] = xs[xs.size(1) * xs2_tmp];
          xs2[b_xs2_tmp + 1] = xs[xs.size(1) * xs2_tmp + 1];
        }

        b_i = nPoints + 1;
        for (int j{2}; j <= b_i; j++) {
          xs2_tmp = (j - 1) << 1;
          xs2[xs2_tmp] -= xs[0];
          xs2[xs2_tmp + 1] -= xs[1];
        }

        if ((nrms.size(0) != 0) && (nrms.size(1) != 0)) {
          nrm02_idx_0 = nrms[nrms.size(1) * i];
          nrm02_idx_1 = nrms[nrms.size(1) * i + 1];
        } else if (surfType == 1) {
          double a;

          //  circle
          nrm02_idx_0 = xs[xs.size(1) * i];
          nrm02_idx_1 = xs[xs.size(1) * i + 1];
          a = std::sqrt(nrm02_idx_0 * nrm02_idx_0 + nrm02_idx_1 * nrm02_idx_1);
          nrm02_idx_0 /= a;
          nrm02_idx_1 /= a;
        }

        //  compute tangent
        for (int j{0}; j <= nPoints; j++) {
          b_xs2_tmp = j << 1;
          us2[j] = xs2[b_xs2_tmp] * -nrm02_idx_1 + xs2[b_xs2_tmp + 1] *
            nrm02_idx_0;
        }

        h0 = 0.0;
        for (int j{2}; j <= b_i; j++) {
          h0 += std::abs(us2[j - 1]);
        }

        buf[i] = h0 / static_cast<double>(nPoints);
      }
    }

    //  serial
    for (int e{0}; e < m; e++) {
      int npts;
      for (npts = conn.size(1); conn[(npts + conn.size(1) * e) - 1] <= 0; npts--)
      {
      }

      h0 = 0.0;
      for (int i{0}; i < npts; i++) {
        h0 += buf[conn[i + conn.size(1) * e] - 1];
      }

      h[e] = h0 / static_cast<double>(npts);
    }
  }

  static void compute_surf_h_kernel(int n, int m, int dim, const ::coder::array<
    double, 2U> &xs, const ::coder::array<int, 2U> &conn, const ::coder::array<
    int, 1U> &n2nPtr, const ::coder::array<int, 1U> &n2nList, int surfType, ::
    coder::array<double, 1U> &h, ::coder::array<double, 1U> &buf)
  {
    double d1;
    double d2;
    double d3;
    double h0;
    double nrm02_idx_0;
    double nrm02_idx_1;
    double nrm03_idx_0;
    double nrm03_idx_1;
    double nrm03_idx_2;

    //  serial
    if (dim == 3) {
      double xs3[384];
      double us3[256];

      //  nodal h
      for (int i{0}; i < n; i++) {
        double a;
        double d;
        double t3_idx_1;
        double t3_idx_3;
        double t3_idx_5;
        int b_i;
        int b_xs2_tmp;
        int nPoints;
        int xs2_tmp;
        boolean_T guard1{ false };

        nPoints = n2nPtr[i + 1] - n2nPtr[i];

        //  i is excluded
        xs3[0] = xs[xs.size(1) * i];
        xs3[1] = xs[xs.size(1) * i + 1];
        xs3[2] = xs[xs.size(1) * i + 2];
        for (int j{0}; j < nPoints; j++) {
          xs2_tmp = n2nList[(n2nPtr[i] + j) - 1] - 1;
          b_xs2_tmp = 3 * (j + 1);
          xs3[b_xs2_tmp] = xs[xs.size(1) * xs2_tmp];
          xs3[b_xs2_tmp + 1] = xs[xs.size(1) * xs2_tmp + 1];
          xs3[b_xs2_tmp + 2] = xs[xs.size(1) * xs2_tmp + 2];
        }

        b_i = nPoints + 1;
        for (int j{2}; j <= b_i; j++) {
          xs2_tmp = 3 * (j - 1);
          xs3[xs2_tmp] -= xs3[0];
          xs3[xs2_tmp + 1] -= xs3[1];
          xs3[xs2_tmp + 2] -= xs3[2];
        }

        if (surfType == 1) {
          //  sphere
          nrm03_idx_0 = xs[xs.size(1) * i];
          nrm03_idx_1 = xs[xs.size(1) * i + 1];
          nrm03_idx_2 = xs[xs.size(1) * i + 2];
          a = std::sqrt((nrm03_idx_0 * nrm03_idx_0 + nrm03_idx_1 * nrm03_idx_1)
                        + nrm03_idx_2 * nrm03_idx_2);
          nrm03_idx_0 /= a;
          nrm03_idx_1 /= a;
          nrm03_idx_2 /= a;
        }

        d = std::abs(nrm03_idx_0);
        guard1 = false;
        if (d > std::abs(nrm03_idx_1)) {
          d1 = nrm03_idx_2;
          if (d > std::abs(nrm03_idx_2)) {
            d = -nrm03_idx_0 * nrm03_idx_1;
            d2 = 1.0 - nrm03_idx_1 * nrm03_idx_1;
            d3 = -nrm03_idx_1 * nrm03_idx_2;
          } else {
            guard1 = true;
          }
        } else {
          guard1 = true;
        }

        if (guard1) {
          d = 1.0 - nrm03_idx_0 * nrm03_idx_0;
          d2 = -nrm03_idx_0 * nrm03_idx_1;
          d1 = nrm03_idx_2;
          d3 = -nrm03_idx_0 * nrm03_idx_2;
        }

        a = std::sqrt((d * d + d2 * d2) + d3 * d3);
        d /= a;
        d2 /= a;
        d3 /= a;

        //  cross
        t3_idx_1 = nrm03_idx_1 * d3 - d1 * d2;
        t3_idx_3 = d1 * d - nrm03_idx_0 * d3;
        t3_idx_5 = nrm03_idx_0 * d2 - nrm03_idx_1 * d;

        //  compute tangent
        for (int j{0}; j <= nPoints; j++) {
          double d4;
          d1 = xs3[3 * j + 1];
          d4 = xs3[3 * j + 2];
          b_xs2_tmp = j << 1;
          us3[b_xs2_tmp] = (xs3[3 * j] * d + d2 * d1) + d3 * d4;
          us3[b_xs2_tmp + 1] = (d1 * t3_idx_1 + d1 * t3_idx_3) + t3_idx_5 * d4;
        }

        h0 = 0.0;
        for (int j{2}; j <= b_i; j++) {
          b_xs2_tmp = (j - 1) << 1;
          d = us3[b_xs2_tmp];
          d1 = us3[b_xs2_tmp + 1];
          h0 += std::sqrt(d * d + d1 * d1);
        }

        buf[i] = h0 / static_cast<double>(nPoints);
      }
    } else {
      double xs2[256];
      double us2[128];

      //  2D
      for (int i{0}; i < n; i++) {
        int b_i;
        int b_xs2_tmp;
        int nPoints;
        int xs2_tmp;
        nPoints = n2nPtr[i + 1] - n2nPtr[i];

        //  i is excluded
        xs2[0] = xs[xs.size(1) * i];
        xs2[1] = xs[xs.size(1) * i + 1];
        for (int j{0}; j < nPoints; j++) {
          xs2_tmp = (n2nList[i] + j) - 1;
          b_xs2_tmp = (j + 1) << 1;
          xs2[b_xs2_tmp] = xs[xs.size(1) * xs2_tmp];
          xs2[b_xs2_tmp + 1] = xs[xs.size(1) * xs2_tmp + 1];
        }

        b_i = nPoints + 1;
        for (int j{2}; j <= b_i; j++) {
          xs2_tmp = (j - 1) << 1;
          xs2[xs2_tmp] -= xs[0];
          xs2[xs2_tmp + 1] -= xs[1];
        }

        if (surfType == 1) {
          double a;

          //  circle
          nrm02_idx_0 = xs[xs.size(1) * i];
          nrm02_idx_1 = xs[xs.size(1) * i + 1];
          a = std::sqrt(nrm02_idx_0 * nrm02_idx_0 + nrm02_idx_1 * nrm02_idx_1);
          nrm02_idx_0 /= a;
          nrm02_idx_1 /= a;
        }

        //  compute tangent
        for (int j{0}; j <= nPoints; j++) {
          b_xs2_tmp = j << 1;
          us2[j] = xs2[b_xs2_tmp] * -nrm02_idx_1 + xs2[b_xs2_tmp + 1] *
            nrm02_idx_0;
        }

        h0 = 0.0;
        for (int j{2}; j <= b_i; j++) {
          h0 += std::abs(us2[j - 1]);
        }

        buf[i] = h0 / static_cast<double>(nPoints);
      }
    }

    //  serial
    for (int e{0}; e < m; e++) {
      int npts;
      for (npts = conn.size(1); conn[(npts + conn.size(1) * e) - 1] <= 0; npts--)
      {
      }

      h0 = 0.0;
      for (int i{0}; i < npts; i++) {
        h0 += buf[conn[i + conn.size(1) * e] - 1];
      }

      h[e] = h0 / static_cast<double>(npts);
    }
  }

  static void compute_surf_h_kernel(int n, int m, int dim, const ::coder::array<
    double, 2U> &xs, const ::coder::array<int, 2U> &conn, const ::coder::array<
    int, 1U> &n2nPtr, const ::coder::array<int, 1U> &n2nList, int surfType, ::
    coder::array<double, 1U> &h, ::coder::array<double, 1U> &buf, int nThreads)
  {
    double d1;
    double h0;
    double nrm02_idx_0;
    double nrm02_idx_1;
    double nrm03_idx_0;
    double nrm03_idx_1;
    double nrm03_idx_2;
    int b_remainder;
    int chunk;
    int iendM;
    int iendN;
    int istartM;
    int istartN;
    int threadID;
    int u1;

    //  Obtains starting and ending indices of local chunk for current thread
    threadID = 0;

#ifdef _OPENMP

    threadID = omp_get_thread_num();

#endif //_OPENMP

    if (threadID >= nThreads) {
      //  The excess threads would not be assigned any task
      istartN = 1;
      iendN = 0;
    } else {
      chunk = n / nThreads;
      b_remainder = n - nThreads * chunk;
      u1 = threadID;
      if (b_remainder <= threadID) {
        u1 = b_remainder;
      }

      istartN = (threadID * chunk + u1) + 1;
      iendN = ((istartN + chunk) + (threadID < b_remainder)) - 1;
    }

    if (dim == 3) {
      double xs3[384];
      double us3[256];

      //  nodal h
      for (int i{istartN}; i <= iendN; i++) {
        double a;
        double b_xs2_tmp;
        double d;
        double t3_idx_1;
        double t3_idx_3;
        double t3_idx_5;
        double xs2_tmp;
        double xs3_tmp;
        int b_i;
        int b_xs3_tmp;
        int nPoints;
        int u0;
        boolean_T guard1{ false };

        u0 = n2nPtr[i - 1];
        nPoints = n2nPtr[i] - u0;

        //  i is excluded
        xs2_tmp = xs[xs.size(1) * (i - 1)];
        xs3[0] = xs2_tmp;
        b_xs2_tmp = xs[xs.size(1) * (i - 1) + 1];
        xs3[1] = b_xs2_tmp;
        xs3_tmp = xs[xs.size(1) * (i - 1) + 2];
        xs3[2] = xs3_tmp;
        for (int j{0}; j < nPoints; j++) {
          b_xs3_tmp = n2nList[(u0 + j) - 1] - 1;
          u1 = 3 * (j + 1);
          xs3[u1] = xs[xs.size(1) * b_xs3_tmp];
          xs3[u1 + 1] = xs[xs.size(1) * b_xs3_tmp + 1];
          xs3[u1 + 2] = xs[xs.size(1) * b_xs3_tmp + 2];
        }

        b_i = nPoints + 1;
        for (int j{2}; j <= b_i; j++) {
          b_xs3_tmp = 3 * (j - 1);
          xs3[b_xs3_tmp] -= xs3[0];
          xs3[b_xs3_tmp + 1] -= xs3[1];
          xs3[b_xs3_tmp + 2] -= xs3[2];
        }

        if (surfType == 1) {
          //  sphere
          a = std::sqrt((xs2_tmp * xs2_tmp + b_xs2_tmp * b_xs2_tmp) + xs3_tmp *
                        xs3_tmp);
          nrm03_idx_0 = xs2_tmp / a;
          nrm03_idx_1 = b_xs2_tmp / a;
          nrm03_idx_2 = xs3_tmp / a;
        }

        d = std::abs(nrm03_idx_0);
        guard1 = false;
        if (d > std::abs(nrm03_idx_1)) {
          d1 = nrm03_idx_2;
          if (d > std::abs(nrm03_idx_2)) {
            d = -nrm03_idx_0 * nrm03_idx_1;
            b_xs2_tmp = 1.0 - nrm03_idx_1 * nrm03_idx_1;
            xs3_tmp = -nrm03_idx_1 * nrm03_idx_2;
          } else {
            guard1 = true;
          }
        } else {
          guard1 = true;
        }

        if (guard1) {
          d = 1.0 - nrm03_idx_0 * nrm03_idx_0;
          b_xs2_tmp = -nrm03_idx_0 * nrm03_idx_1;
          d1 = nrm03_idx_2;
          xs3_tmp = -nrm03_idx_0 * nrm03_idx_2;
        }

        a = std::sqrt((d * d + b_xs2_tmp * b_xs2_tmp) + xs3_tmp * xs3_tmp);
        d /= a;
        b_xs2_tmp /= a;
        xs3_tmp /= a;

        //  cross
        t3_idx_1 = nrm03_idx_1 * xs3_tmp - d1 * b_xs2_tmp;
        t3_idx_3 = d1 * d - nrm03_idx_0 * xs3_tmp;
        t3_idx_5 = nrm03_idx_0 * b_xs2_tmp - nrm03_idx_1 * d;

        //  compute tangent
        for (int j{0}; j <= nPoints; j++) {
          d1 = xs3[3 * j + 1];
          xs2_tmp = xs3[3 * j + 2];
          u0 = j << 1;
          us3[u0] = (xs3[3 * j] * d + b_xs2_tmp * d1) + xs3_tmp * xs2_tmp;
          us3[u0 + 1] = (d1 * t3_idx_1 + d1 * t3_idx_3) + t3_idx_5 * xs2_tmp;
        }

        h0 = 0.0;
        for (int j{2}; j <= b_i; j++) {
          u0 = (j - 1) << 1;
          d = us3[u0];
          d1 = us3[u0 + 1];
          h0 += std::sqrt(d * d + d1 * d1);
        }

        buf[i - 1] = h0 / static_cast<double>(nPoints);
      }
    } else {
      double xs2[256];
      double us2[128];

      //  2D
      for (int i{istartN}; i <= iendN; i++) {
        double b_xs2_tmp;
        double xs2_tmp;
        int b_i;
        int nPoints;
        int u0;
        nPoints = n2nPtr[i] - n2nPtr[i - 1];

        //  i is excluded
        xs2_tmp = xs[xs.size(1) * (i - 1)];
        xs2[0] = xs2_tmp;
        b_xs2_tmp = xs[xs.size(1) * (i - 1) + 1];
        xs2[1] = b_xs2_tmp;
        for (int j{0}; j < nPoints; j++) {
          u1 = (n2nList[i - 1] + j) - 1;
          u0 = (j + 1) << 1;
          xs2[u0] = xs[xs.size(1) * u1];
          xs2[u0 + 1] = xs[xs.size(1) * u1 + 1];
        }

        b_i = nPoints + 1;
        for (int j{2}; j <= b_i; j++) {
          u1 = (j - 1) << 1;
          xs2[u1] -= xs[0];
          xs2[u1 + 1] -= xs[1];
        }

        if (surfType == 1) {
          double a;

          //  circle
          a = std::sqrt(xs2_tmp * xs2_tmp + b_xs2_tmp * b_xs2_tmp);
          nrm02_idx_0 = xs2_tmp / a;
          nrm02_idx_1 = b_xs2_tmp / a;
        }

        //  compute tangent
        for (int j{0}; j <= nPoints; j++) {
          u0 = j << 1;
          us2[j] = xs2[u0] * -nrm02_idx_1 + xs2[u0 + 1] * nrm02_idx_0;
        }

        h0 = 0.0;
        for (int j{2}; j <= b_i; j++) {
          h0 += std::abs(us2[j - 1]);
        }

        buf[i - 1] = h0 / static_cast<double>(nPoints);
      }
    }

#pragma omp barrier

    //  Obtains starting and ending indices of local chunk for current thread
    threadID = 0;

#ifdef _OPENMP

    threadID = omp_get_thread_num();

#endif //_OPENMP

    if (threadID >= nThreads) {
      //  The excess threads would not be assigned any task
      istartM = 0;
      iendM = 0;
    } else {
      chunk = m / nThreads;
      b_remainder = m - nThreads * chunk;
      u1 = threadID;
      if (b_remainder <= threadID) {
        u1 = b_remainder;
      }

      istartM = threadID * chunk + u1;
      iendM = (istartM + chunk) + (threadID < b_remainder);
    }

    //  average to cell-based average h
    for (int e{istartM + 1}; e <= iendM; e++) {
      int npts;
      for (npts = conn.size(1); conn[(npts + conn.size(1) * (e - 1)) - 1] <= 0;
           npts--) {
      }

      h0 = 0.0;
      for (int i{0}; i < npts; i++) {
        h0 += buf[conn[i + conn.size(1) * (e - 1)] - 1];
      }

      h[e - 1] = h0 / static_cast<double>(npts);
    }
  }

  static inline
  void compute_volume_tet(const ::coder::array<double, 2U> &xs, const ::
    coder::array<int, 2U> &conn, ::coder::array<double, 1U> &V, int nThreads)
  {
    //  kernel for estimating the volumes for a tet mesh
    if (nThreads <= 1) {
      compute_volume_tet_kernel(conn.size(0), xs, conn, V);
    } else {
      int m2cTryBlkErrCode;

      m2cTryBlkErrCode = (0);

#pragma omp parallel default(shared) num_threads(nThreads)

      try { //try
        compute_volume_tet_kernel(conn.size(0), xs, conn, V, nThreads);

      } catch (const std::runtime_error& m2cExc) {
        m2cTryBlkErrCode = (1);
        m2cPrintError("runtime_error %s\n", m2cExc.what());

      }

      catch (const std::logic_error& m2cExc)
      {
        m2cTryBlkErrCode = (1);
        m2cPrintError("logic_error %s\n", m2cExc.what());

      }

      catch (const std::exception& m2cExc)
      {
        m2cTryBlkErrCode = (1);
        m2cPrintError("exception %s\n", m2cExc.what());

      }

      catch (...)
      {
        m2cTryBlkErrCode = (1);
        m2cPrintError("Unknown error detected from C++ exceptions\n");

      } //try

      if ((int)m2cTryBlkErrCode != 0) {
        throw std::runtime_error("omp4m:runtimeErrorInThread");
      }
    }
  }

  static void compute_volume_tet_kernel(int m, const ::coder::array<double, 2U>
    &xs, const ::coder::array<int, 2U> &conn, ::coder::array<double, 1U> &V, int
    nThreads)
  {
    int iend;
    int istart;
    int threadID;
    int u1;

    //  Obtains starting and ending indices of local chunk for current thread
    threadID = 0;

#ifdef _OPENMP

    threadID = omp_get_thread_num();

#endif //_OPENMP

    if (threadID >= nThreads) {
      //  The excess threads would not be assigned any task
      istart = 0;
      iend = 0;
    } else {
      int b_remainder;
      int chunk;
      chunk = m / nThreads;
      b_remainder = m - nThreads * chunk;
      u1 = threadID;
      if (b_remainder <= threadID) {
        u1 = b_remainder;
      }

      istart = threadID * chunk + u1;
      iend = (istart + chunk) + (threadID < b_remainder);
    }

    for (int e{istart + 1}; e <= iend; e++) {
      double J_idx_0;
      double J_idx_0_tmp;
      double J_idx_1;
      double J_idx_3;
      double J_idx_3_tmp;
      double J_idx_4;
      double J_idx_6;
      double J_idx_6_tmp;
      double J_idx_7;
      int i;
      int i1;
      int u0;
      u0 = conn[conn.size(1) * (e - 1) + 1];
      u1 = conn[conn.size(1) * (e - 1)];
      J_idx_0_tmp = xs[xs.size(1) * (u1 - 1)];
      J_idx_0 = xs[xs.size(1) * (u0 - 1)] - J_idx_0_tmp;
      i = conn[conn.size(1) * (e - 1) + 2];
      J_idx_1 = xs[xs.size(1) * (i - 1)] - J_idx_0_tmp;
      i1 = conn[conn.size(1) * (e - 1) + 3];
      J_idx_3_tmp = xs[xs.size(1) * (u1 - 1) + 1];
      J_idx_3 = xs[xs.size(1) * (u0 - 1) + 1] - J_idx_3_tmp;
      J_idx_4 = xs[xs.size(1) * (i - 1) + 1] - J_idx_3_tmp;
      J_idx_6_tmp = xs[xs.size(1) * (u1 - 1) + 2];
      J_idx_6 = xs[xs.size(1) * (u0 - 1) + 2] - J_idx_6_tmp;
      J_idx_7 = xs[xs.size(1) * (i - 1) + 2] - J_idx_6_tmp;
      V[e - 1] = std::abs(((xs[xs.size(1) * (i1 - 1)] - J_idx_0_tmp) * (J_idx_3 *
        J_idx_7 - J_idx_4 * J_idx_6) + (xs[xs.size(1) * (i1 - 1) + 1] -
        J_idx_3_tmp) * (J_idx_1 * J_idx_6 - J_idx_0 * J_idx_7)) + (xs[xs.size(1)
        * (i1 - 1) + 2] - J_idx_6_tmp) * (J_idx_0 * J_idx_4 - J_idx_1 * J_idx_3))
        / 6.0;
    }
  }

  static void compute_volume_tet_kernel(int m, const ::coder::array<double, 2U>
    &xs, const ::coder::array<int, 2U> &conn, ::coder::array<double, 1U> &V)
  {
    for (int e{0}; e < m; e++) {
      double J_idx_0;
      double J_idx_0_tmp;
      double J_idx_1;
      double J_idx_3;
      double J_idx_3_tmp;
      double J_idx_4;
      double J_idx_6;
      double J_idx_6_tmp;
      double J_idx_7;
      int i;
      int i1;
      int i2;
      int i3;
      i = conn[conn.size(1) * e + 1];
      i1 = conn[conn.size(1) * e];
      J_idx_0_tmp = xs[xs.size(1) * (i1 - 1)];
      J_idx_0 = xs[xs.size(1) * (i - 1)] - J_idx_0_tmp;
      i2 = conn[conn.size(1) * e + 2];
      J_idx_1 = xs[xs.size(1) * (i2 - 1)] - J_idx_0_tmp;
      i3 = conn[conn.size(1) * e + 3];
      J_idx_3_tmp = xs[xs.size(1) * (i1 - 1) + 1];
      J_idx_3 = xs[xs.size(1) * (i - 1) + 1] - J_idx_3_tmp;
      J_idx_4 = xs[xs.size(1) * (i2 - 1) + 1] - J_idx_3_tmp;
      J_idx_6_tmp = xs[xs.size(1) * (i1 - 1) + 2];
      J_idx_6 = xs[xs.size(1) * (i - 1) + 2] - J_idx_6_tmp;
      J_idx_7 = xs[xs.size(1) * (i2 - 1) + 2] - J_idx_6_tmp;
      V[e] = std::abs(((xs[xs.size(1) * (i3 - 1)] - J_idx_0_tmp) * (J_idx_3 *
        J_idx_7 - J_idx_4 * J_idx_6) + (xs[xs.size(1) * (i3 - 1) + 1] -
        J_idx_3_tmp) * (J_idx_1 * J_idx_6 - J_idx_0 * J_idx_7)) + (xs[xs.size(1)
        * (i3 - 1) + 2] - J_idx_6_tmp) * (J_idx_0 * J_idx_4 - J_idx_1 * J_idx_3))
        / 6.0;
    }
  }

  static void crsAx_kernel(const ::coder::array<int, 1U> &row_ptr, const ::coder::
    array<int, 1U> &col_ind, const ::coder::array<double, 1U> &val, int nrows,
    const ::coder::array<double, 2U> &x, int nrhs, ::coder::array<double, 2U> &b)
  {
    int iend;
    int istart;
    int nthreads;
    int u1;

    //  Never inline the function, so that all variables are private
    nthreads = 1;

#ifdef _OPENMP

    nthreads = omp_get_num_threads();

#endif //_OPENMP

    if (nthreads == 1) {
      istart = 0;
      iend = nrows;
    } else {
      int b_remainder;
      int chunk;
      int threadID;

      //  Gets the thread number of the thread, within the team, making this call.
      threadID = 0;

#ifdef _OPENMP

      threadID = omp_get_thread_num();

#endif //_OPENMP

      chunk = nrows / nthreads;
      b_remainder = nrows - nthreads * chunk;
      u1 = threadID;
      if (b_remainder <= threadID) {
        u1 = b_remainder;
      }

      istart = threadID * chunk + u1;
      iend = (istart + chunk) + (threadID < b_remainder);
    }

    //  Optimized for row-major
    for (int i{istart + 1}; i <= iend; i++) {
      int u0;
      for (int k{0}; k < nrhs; k++) {
        b[k + b.size(1) * (i - 1)] = 0.0;
      }

      u0 = row_ptr[i - 1];
      u1 = row_ptr[i] - 1;
      for (int j{u0}; j <= u1; j++) {
        for (int k{0}; k < nrhs; k++) {
          b[k + b.size(1) * (i - 1)] = b[k + b.size(1) * (i - 1)] + val[j - 1] *
            x[k + x.size(1) * (col_ind[j - 1] - 1)];
        }
      }
    }
  }

  static inline
  void crs_prod_mat_vec(const ::coder::array<int, 1U> &A_rowptr, const ::
    coder::array<int, 1U> &A_colind, const ::coder::array<double, 1U> &A_val,
    const ::coder::array<double, 2U> &x, ::coder::array<double, 2U> &b)
  {
    //  Compute b=A*x for an m-by-n CRS matrix A and column vector(s) x.

#pragma omp single

    { //single
      b.set_size(A_rowptr.size(0) - 1, x.size(1));

    } //single
 
    //      Compute b=A*x in parallel
    crsAx_kernel(A_rowptr, A_colind, A_val, A_rowptr.size(0) - 1, x, x.size(1),
                 b);
  }

  static void d_assemble_body_kernel(const ::coder::array<double, 2U> &xs, const
    ::coder::array<int, 2U> &conn, const ::coder::array<int, 2U> &stcls, const ::
    coder::array<int, 1U> &n2cPtr, const ::coder::array<int, 1U> &n2cList, int
    degree, const ::coder::array<int, 1U> &rowPtr, const ::coder::array<int, 1U>
    &colInd, ::coder::array<double, 1U> &vals, const ::coder::array<int, 1U>
    &nnzPr, const ::coder::array<int, 1U> &nrange, int istart, int iend,
    WlsDataStruct *wls, ::coder::array<boolean_T, 1U> &rdTags, int *rdCounts)
  {
    double xsLocal[3072];
    double Ns[1024];
    int eids[1024];
    int gIDs[1024];
    int b_iv[2];
    int iv1[2];
    int i;
    boolean_T exitg1;

    //  kernel
    i = istart - 1;
    exitg1 = false;
    while ((!exitg1) && (i + 1 <= iend)) {
      double ori_idx_0;
      double ori_idx_1;
      double ori_idx_2;
      int b_xsLocal_tmp;
      int k;
      int nPoints;
      int nid;
      int xsLocal_tmp;

      //  get nodal ID
      nid = nrange[i] - 1;
      nPoints = stcls[(stcls.size(1) + stcls.size(1) * (nrange[i] - 1)) - 1] - 1;

      //  fetch coords, global IDs, and then localize coordinates
      for (int j{0}; j <= nPoints; j++) {
        gIDs[j] = stcls[j + stcls.size(1) * nid];
      }

      //  fetch coordinates
      for (k = 0; k <= nPoints; k++) {
        xsLocal_tmp = gIDs[k];
        xsLocal[3 * k] = xs[xs.size(1) * (xsLocal_tmp - 1)];
        xsLocal[3 * k + 1] = xs[xs.size(1) * (xsLocal_tmp - 1) + 1];
        xsLocal[3 * k + 2] = xs[xs.size(1) * (xsLocal_tmp - 1) + 2];
      }

      //  origin
      ori_idx_0 = xsLocal[0];
      ori_idx_1 = xsLocal[1];
      ori_idx_2 = xsLocal[2];

      //  localize
      for (k = 0; k <= nPoints; k++) {
        xsLocal[3 * k] -= ori_idx_0;
        b_xsLocal_tmp = 3 * k + 1;
        xsLocal[b_xsLocal_tmp] -= ori_idx_1;
        b_xsLocal_tmp = 3 * k + 2;
        xsLocal[b_xsLocal_tmp] -= ori_idx_2;
      }

      //  compute wls
      b_iv[0] = 1024;
      b_iv[1] = 3;
      wls_init(&wls->wlsObj, xsLocal, b_iv, wls->wlsWgts.name.data,
               wls->wlsWgts.params_pointwise.data,
               wls->wlsWgts.params_pointwise.size, degree, wls->interp0,
               wls->useDag, nPoints + 1);
      if (wls->wlsObj.rank < 0) {
        *rdCounts = -nrange[i];

        //  record error node and return
        exitg1 = true;
      } else {
        if (wls->wlsObj.rank != wls->wlsObj.ncols) {
          (*rdCounts)++;

          //  increment counter for RD nodes
          rdTags[nrange[i] - 1] = true;
        } else {
          int c_xsLocal_tmp;
          int nCells;

          //  get centers and shape function Ns at centers
          nCells = (n2cPtr[nrange[i]] - n2cPtr[nrange[i] - 1]) - 1;
          for (int j{0}; j <= nCells; j++) {
            int eid;
            int npts;
            eid = n2cList[(n2cPtr[nid] + j) - 1] - 1;

            eids[j] = eid + 1;
            for (npts = conn.size(1); conn[(npts + conn.size(1) * eid) - 1] <= 0;
                 npts--) {
            }

            Ns[j] = 1.0 / static_cast<double>(npts);

            //  compute center
            xsLocal[3 * j] = xs[xs.size(1) * (conn[conn.size(1) * eid] - 1)];
            b_xsLocal_tmp = 3 * j + 1;
            xsLocal[b_xsLocal_tmp] = xs[xs.size(1) * (conn[conn.size(1) * eid] -
              1) + 1];
            c_xsLocal_tmp = 3 * j + 2;
            xsLocal[c_xsLocal_tmp] = xs[xs.size(1) * (conn[conn.size(1) * eid] -
              1) + 2];

            //  C++
            for (int ii{2}; ii <= npts; ii++) {
              xsLocal_tmp = conn[(ii + conn.size(1) * eid) - 1] - 1;
              xsLocal[3 * j] += xs[xs.size(1) * xsLocal_tmp];
              xsLocal[b_xsLocal_tmp] += xs[xs.size(1) * xsLocal_tmp + 1];
              xsLocal[c_xsLocal_tmp] += xs[xs.size(1) * xsLocal_tmp + 2];
            }

            double d;
            double d1;
            double d2;
            d = Ns[j];
            d1 = xsLocal[b_xsLocal_tmp] * d;
            d2 = xsLocal[c_xsLocal_tmp] * d;

            //  localize
            xsLocal[3 * j] = xsLocal[3 * j] * d - ori_idx_0;
            d1 -= ori_idx_1;
            xsLocal[b_xsLocal_tmp] = d1;
            d2 -= ori_idx_2;
            xsLocal[c_xsLocal_tmp] = d2;
          }

          iv1[0] = 1024;
          iv1[1] = 3;
          wls_func(&wls->wlsObj, xsLocal, iv1, nCells + 1, wls->coeffs);

          //  Compute local OSUS coefficients
          for (k = 0; k <= nPoints; k++) {
            for (int j{0}; j <= nCells; j++) {
              wls->coeffs[j + wls->coeffs.size(1) * k] = wls->coeffs[j +
                wls->coeffs.size(1) * k] * Ns[j];
            }
          }

          //  substract linear interp
          for (int j{0}; j <= nCells; j++) {
            wls->coeffs[j] = wls->coeffs[j] - Ns[j];
          }

          //  kernel for updating crs matrix
          for (int j{0}; j <= nCells; j++) {
            int rowIdx;
            rowIdx = eids[j] - 1;
            for (int ii{0}; ii <= nPoints; ii++) {
              int b_i;
              c_xsLocal_tmp = gIDs[ii];
              b_i = rowPtr[rowIdx] - 1;
              xsLocal_tmp = (rowPtr[rowIdx] + nnzPr[eids[j] - 1]) - 1;
              int exitg2;
              do {
                exitg2 = 0;
                if (b_i + 1 <= xsLocal_tmp) {
                  if (colInd[b_i] == c_xsLocal_tmp) {
                    k = b_i;
                    exitg2 = 1;
                  } else {
                    b_i++;
                  }
                } else {
                  k = -1;

                  m2cErrMsgIdAndTxt("add_crs:missingIndex",
                                    "could not find column index %d",
                                    c_xsLocal_tmp);
                  exitg2 = 1;
                }
              } while (exitg2 == 0);

              vals[k] = vals[k] + wls->coeffs[j + wls->coeffs.size(1) * ii];
            }
          }
        }

        i++;
      }
    }
  }

  static void d_assemble_surf_kernel(const ::coder::array<double, 2U> &xs, const
    ::coder::array<int, 2U> &conn, const ::coder::array<int, 2U> &stcls, const ::
    coder::array<int, 1U> &n2cPtr, const ::coder::array<int, 1U> &n2cList, int
    surfType, int degree, const ::coder::array<double, 2U> &nrms, const ::coder::
    array<int, 1U> &rowPtr, const ::coder::array<int, 1U> &colInd, ::coder::
    array<double, 1U> &vals, const ::coder::array<int, 1U> &nnzPr, const ::coder::
    array<int, 1U> &nrange, int istart, int iend, WlsDataStruct *wls, ::coder::
    array<boolean_T, 1U> &rdTags, int *rdCounts)
  {
    double xsLocal[1024];
    double Ns[512];
    double us[512];
    double nrm0_idx_0;
    double nrm0_idx_1;
    int eids[512];
    int gIDs[512];
    int i;
    boolean_T exitg1;
    boolean_T isSphSurf;

    //  kernel for surface computation
    wls->wlsWgts.params_pointwise.size[1] = 1;
    wls->wlsWgts.params_pointwise.size[0] = 512;

    //  check if spherical surface
    isSphSurf = (surfType == 1);

    //  loop begins
    i = istart - 1;
    exitg1 = false;
    while ((!exitg1) && (i + 1 <= iend)) {
      double a;
      double ori_idx_0;
      double ori_idx_1;
      int b_i;
      int k;
      int nPoints;
      int nid;
      int us_tmp;
      int xsLocal_tmp;

      //  get nodal ID
      nid = nrange[i] - 1;
      nPoints = stcls[(stcls.size(1) + stcls.size(1) * (nrange[i] - 1)) - 1] - 1;

      //  fetch data to local buffer and localize the coordinates
      for (int j{0}; j <= nPoints; j++) {
        gIDs[j] = stcls[j + stcls.size(1) * nid];
      }

      //  fetch coordinates
      for (k = 0; k <= nPoints; k++) {
        xsLocal_tmp = k << 1;
        b_i = gIDs[k];
        xsLocal[xsLocal_tmp] = xs[xs.size(1) * (b_i - 1)];
        xsLocal[xsLocal_tmp + 1] = xs[xs.size(1) * (b_i - 1) + 1];
      }

      //  origin
      ori_idx_0 = xsLocal[0];
      ori_idx_1 = xsLocal[1];

      //  localize
      for (k = 0; k <= nPoints; k++) {
        xsLocal_tmp = k << 1;
        xsLocal[xsLocal_tmp] -= ori_idx_0;
        xsLocal[xsLocal_tmp + 1] -= ori_idx_1;
      }

      //  get normal direction
      if ((nrms.size(0) != 0) && (nrms.size(1) != 0)) {
        nrm0_idx_0 = nrms[nrms.size(1) * (nrange[i] - 1)];
        nrm0_idx_1 = nrms[nrms.size(1) * (nrange[i] - 1) + 1];
      } else if (isSphSurf) {
        //  spherical or circle
        nrm0_idx_0 = xs[xs.size(1) * (nrange[i] - 1)];
        nrm0_idx_1 = xs[xs.size(1) * (nrange[i] - 1) + 1];
      }

      //  compute tagent
      if (nPoints >= 0) {
        std::memset(&us[0], 0, (nPoints + 1) * sizeof(double));
      }

      //  matrix-matrix, using mem-efficient loop
      for (int ii{0}; ii <= nPoints; ii++) {
        us_tmp = ii << 1;
        us[ii] = (us[ii] + xsLocal[us_tmp] * -nrm0_idx_1) + nrm0_idx_0 *
          xsLocal[us_tmp + 1];
      }

      //  compute scaling to the WLS weights by inner product of normals
      if ((nrms.size(0) != 0) && (nrms.size(1) != 0)) {
        for (int j{2}; j <= nPoints + 1; j++) {
          us_tmp = stcls[(j + stcls.size(1) * nid) - 1] - 1;
          a = nrm0_idx_0 * nrms[nrms.size(1) * us_tmp] + nrm0_idx_1 *
            nrms[nrms.size(1) * us_tmp + 1];
          if (a < 0.0) {
            a = 0.0;
          }

          wls->wlsWgts.params_pointwise.data[j - 1] = a;
        }
      } else if (isSphSurf) {
        //  spherical or circle
        for (int j{2}; j <= nPoints + 1; j++) {
          int idx;
          idx = stcls[(j + stcls.size(1) * nid) - 1] - 1;
          a = nrm0_idx_0 * xs[xs.size(1) * idx] + nrm0_idx_1 * xs[xs.size(1) *
            idx + 1];
          if (a < 0.0) {
            a = 0.0;
          }

          wls->wlsWgts.params_pointwise.data[j - 1] = a;
        }
      }

      wls->wlsWgts.params_pointwise.data[0] = 1.0;

      //  compute wls
      b_wls_init(&wls->wlsObj, us, 512, wls->wlsWgts.name.data,
                 wls->wlsWgts.params_pointwise.data,
                 wls->wlsWgts.params_pointwise.size, degree, wls->interp0,
                 wls->useDag, nPoints + 1);
      if (wls->wlsObj.rank < 0) {
        *rdCounts = -nrange[i];

        //  record error node and return
        exitg1 = true;
      } else {
        if (wls->wlsObj.rank != wls->wlsObj.ncols) {
          (*rdCounts)++;

          //  increment counter for RD nodes
          rdTags[nrange[i] - 1] = true;
        } else {
          int nCells;

          //  get centers and shape function Ns at centers
          nCells = (n2cPtr[nrange[i]] - n2cPtr[nrange[i] - 1]) - 1;
          for (int j{0}; j <= nCells; j++) {
            double d;
            double d1;
            int eid;
            int npts;
            eid = n2cList[(n2cPtr[nid] + j) - 1] - 1;

            eids[j] = eid + 1;
            for (npts = conn.size(1); conn[(npts + conn.size(1) * eid) - 1] <= 0;
                 npts--) {
            }

            Ns[j] = 1.0 / static_cast<double>(npts);

            //  compute center
            xsLocal_tmp = j << 1;
            xsLocal[xsLocal_tmp] = xs[xs.size(1) * (conn[conn.size(1) * eid] - 1)];
            xsLocal[xsLocal_tmp + 1] = xs[xs.size(1) * (conn[conn.size(1) * eid]
              - 1) + 1];

            //  C++
            for (int ii{2}; ii <= npts; ii++) {
              us_tmp = conn[(ii + conn.size(1) * eid) - 1] - 1;
              xsLocal[xsLocal_tmp] += xs[xs.size(1) * us_tmp];
              xsLocal[xsLocal_tmp + 1] += xs[xs.size(1) * us_tmp + 1];
            }

            d = Ns[j];
            d1 = xsLocal[xsLocal_tmp] * d;
            d *= xsLocal[xsLocal_tmp + 1];
            if (isSphSurf) {
              a = std::sqrt(d1 * d1 + d * d);
              d1 /= a;
              d /= a;
            }

            //  localize
            d1 -= ori_idx_0;
            xsLocal[xsLocal_tmp] = d1;
            d -= ori_idx_1;
            xsLocal[xsLocal_tmp + 1] = d;
          }

          if (nCells >= 0) {
            std::memset(&us[0], 0, (nCells + 1) * sizeof(double));
          }

          //  matrix-matrix, using mem-efficient loop
          for (int ii{0}; ii <= nCells; ii++) {
            us_tmp = ii << 1;
            us[ii] = (us[ii] + xsLocal[us_tmp] * -nrm0_idx_1) + nrm0_idx_0 *
              xsLocal[us_tmp + 1];
          }

          wls_func(&wls->wlsObj, us, nCells + 1, wls->coeffs);

          //  Compute local OSUS coefficients
          for (k = 0; k <= nPoints; k++) {
            for (int j{0}; j <= nCells; j++) {
              wls->coeffs[j + wls->coeffs.size(1) * k] = wls->coeffs[j +
                wls->coeffs.size(1) * k] * Ns[j];
            }
          }

          //  substract linear interp
          for (int j{0}; j <= nCells; j++) {
            wls->coeffs[j] = wls->coeffs[j] - Ns[j];
          }

          //  kernel for updating crs matrix
          for (int j{0}; j <= nCells; j++) {
            int rowIdx;
            rowIdx = eids[j] - 1;
            for (int ii{0}; ii <= nPoints; ii++) {
              int c_i;
              us_tmp = gIDs[ii];
              c_i = rowPtr[rowIdx] - 1;
              b_i = (rowPtr[rowIdx] + nnzPr[eids[j] - 1]) - 1;
              int exitg2;
              do {
                exitg2 = 0;
                if (c_i + 1 <= b_i) {
                  if (colInd[c_i] == us_tmp) {
                    k = c_i;
                    exitg2 = 1;
                  } else {
                    c_i++;
                  }
                } else {
                  k = -1;

                  m2cErrMsgIdAndTxt("add_crs:missingIndex",
                                    "could not find column index %d", us_tmp);
                  exitg2 = 1;
                }
              } while (exitg2 == 0);

              vals[k] = vals[k] + wls->coeffs[j + wls->coeffs.size(1) * ii];
            }
          }
        }

        i++;
      }
    }
  }

  static void e_assemble_body_kernel(const ::coder::array<double, 2U> &xs, const
    ::coder::array<int, 2U> &conn, const ::coder::array<int, 2U> &stcls, const ::
    coder::array<int, 1U> &n2cPtr, const ::coder::array<int, 1U> &n2cList, int
    degree, const ::coder::array<int, 1U> &rowPtr, const ::coder::array<int, 1U>
    &colInd, ::coder::array<double, 1U> &vals, const ::coder::array<int, 1U>
    &nnzPr, const ::coder::array<int, 1U> &nrange, int istart, int iend,
    WlsDataStruct *wls, ::coder::array<boolean_T, 1U> &rdTags, int *rdCounts)
  {
    double xsLocal[1024];
    double Ns[512];
    int eids[512];
    int gIDs[512];
    int b_iv[2];
    int iv1[2];
    int i;
    boolean_T exitg1;

    //  kernel
    i = istart - 1;
    exitg1 = false;
    while ((!exitg1) && (i + 1 <= iend)) {
      double ori_idx_0;
      double ori_idx_1;
      int b_i;
      int k;
      int nPoints;
      int nid;
      int xsLocal_tmp;

      //  get nodal ID
      nid = nrange[i] - 1;
      nPoints = stcls[(stcls.size(1) + stcls.size(1) * (nrange[i] - 1)) - 1] - 1;

      //  fetch coords, global IDs, and then localize coordinates
      for (int j{0}; j <= nPoints; j++) {
        gIDs[j] = stcls[j + stcls.size(1) * nid];
      }

      //  fetch coordinates
      for (k = 0; k <= nPoints; k++) {
        xsLocal_tmp = k << 1;
        b_i = gIDs[k];
        xsLocal[xsLocal_tmp] = xs[xs.size(1) * (b_i - 1)];
        xsLocal[xsLocal_tmp + 1] = xs[xs.size(1) * (b_i - 1) + 1];
      }

      //  origin
      ori_idx_0 = xsLocal[0];
      ori_idx_1 = xsLocal[1];

      //  localize
      for (k = 0; k <= nPoints; k++) {
        xsLocal_tmp = k << 1;
        xsLocal[xsLocal_tmp] -= ori_idx_0;
        xsLocal[xsLocal_tmp + 1] -= ori_idx_1;
      }

      //  compute wls
      b_iv[0] = 512;
      b_iv[1] = 2;
      wls_init(&wls->wlsObj, xsLocal, b_iv, wls->wlsWgts.name.data,
               wls->wlsWgts.params_pointwise.data,
               wls->wlsWgts.params_pointwise.size, degree, wls->interp0,
               wls->useDag, nPoints + 1);
      if (wls->wlsObj.rank < 0) {
        *rdCounts = -nrange[i];

        //  record error node and return
        exitg1 = true;
      } else {
        if (wls->wlsObj.rank != wls->wlsObj.ncols) {
          (*rdCounts)++;

          //  increment counter for RD nodes
          rdTags[nrange[i] - 1] = true;
        } else {
          int b_xsLocal_tmp;
          int nCells;

          //  get centers and shape function Ns at centers
          nCells = (n2cPtr[nrange[i]] - n2cPtr[nrange[i] - 1]) - 1;
          for (int j{0}; j <= nCells; j++) {
            int eid;
            int npts;
            eid = n2cList[(n2cPtr[nid] + j) - 1] - 1;

            eids[j] = eid + 1;
            for (npts = conn.size(1); conn[(npts + conn.size(1) * eid) - 1] <= 0;
                 npts--) {
            }

            Ns[j] = 1.0 / static_cast<double>(npts);

            //  compute center
            xsLocal_tmp = j << 1;
            xsLocal[xsLocal_tmp] = xs[xs.size(1) * (conn[conn.size(1) * eid] - 1)];
            xsLocal[xsLocal_tmp + 1] = xs[xs.size(1) * (conn[conn.size(1) * eid]
              - 1) + 1];

            //  C++
            for (int ii{2}; ii <= npts; ii++) {
              b_xsLocal_tmp = conn[(ii + conn.size(1) * eid) - 1] - 1;
              xsLocal[xsLocal_tmp] += xs[xs.size(1) * b_xsLocal_tmp];
              xsLocal[xsLocal_tmp + 1] += xs[xs.size(1) * b_xsLocal_tmp + 1];
            }

            double d;
            double d1;
            d = Ns[j];
            d1 = xsLocal[xsLocal_tmp + 1] * d;

            //  localize
            xsLocal[xsLocal_tmp] = xsLocal[xsLocal_tmp] * d - ori_idx_0;
            d1 -= ori_idx_1;
            xsLocal[xsLocal_tmp + 1] = d1;
          }

          iv1[0] = 512;
          iv1[1] = 2;
          wls_func(&wls->wlsObj, xsLocal, iv1, nCells + 1, wls->coeffs);

          //  Compute local OSUS coefficients
          for (k = 0; k <= nPoints; k++) {
            for (int j{0}; j <= nCells; j++) {
              wls->coeffs[j + wls->coeffs.size(1) * k] = wls->coeffs[j +
                wls->coeffs.size(1) * k] * Ns[j];
            }
          }

          //  substract linear interp
          for (int j{0}; j <= nCells; j++) {
            wls->coeffs[j] = wls->coeffs[j] - Ns[j];
          }

          //  kernel for updating crs matrix
          for (int j{0}; j <= nCells; j++) {
            int rowIdx;
            rowIdx = eids[j] - 1;
            for (int ii{0}; ii <= nPoints; ii++) {
              int c_i;
              b_xsLocal_tmp = gIDs[ii];
              c_i = rowPtr[rowIdx] - 1;
              b_i = (rowPtr[rowIdx] + nnzPr[eids[j] - 1]) - 1;
              int exitg2;
              do {
                exitg2 = 0;
                if (c_i + 1 <= b_i) {
                  if (colInd[c_i] == b_xsLocal_tmp) {
                    k = c_i;
                    exitg2 = 1;
                  } else {
                    c_i++;
                  }
                } else {
                  k = -1;

                  m2cErrMsgIdAndTxt("add_crs:missingIndex",
                                    "could not find column index %d",
                                    b_xsLocal_tmp);
                  exitg2 = 1;
                }
              } while (exitg2 == 0);

              vals[k] = vals[k] + wls->coeffs[j + wls->coeffs.size(1) * ii];
            }
          }
        }

        i++;
      }
    }
  }

  static void extract_sub(int n, const ::coder::array<int, 1U> &crange, const ::
    coder::array<int, 1U> &eptr, const ::coder::array<int, 1U> &eind, ::coder::
    array<int, 1U> &iwork, ::coder::array<boolean_T, 1U> &ntags, ::coder::array<
    int, 1U> &eptrloc, ::coder::array<int, 1U> &eindloc, int *nnodes)
  {
    ::coder::array<int, 1U> l2g_;
    int j;
    int k;
    int m;
    int u0;
    int y;

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
    for (int i{0}; i <= m; i++) {
      eptrloc[i + 1] = (eptrloc[i] + eptr[crange[i]]) - eptr[crange[i] - 1];
    }

    eindloc.set_size(eptrloc[crange.size(0)] - 1);

    //  build local to global cell ID map
    for (int i{0}; i <= m; i++) {
      int n0;
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
    for (int i{0}; i < n; i++) {
      if (ntags[i]) {
        j++;
        l2g_[j] = i + 1;
        ntags[i] = false;

        //  reset
      }
    }

    //  global to local
    for (int i{0}; i < y; i++) {
      iwork[l2g_[i] - 1] = i + 1;
    }

    k = eindloc.size(0);
    for (int i{0}; i < k; i++) {
      eindloc[i] = iwork[eindloc[i] - 1];
    }

    *nnodes = y;
  }

  static void f_assemble_body_kernel(const ::coder::array<double, 2U> &xs, const
    ::coder::array<int, 2U> &conn, const ::coder::array<int, 2U> &stcls, const ::
    coder::array<int, 1U> &n2cPtr, const ::coder::array<int, 1U> &n2cList, int
    degree, const ::coder::array<int, 1U> &rowPtr, const ::coder::array<int, 1U>
    &colInd, ::coder::array<double, 1U> &vals, const ::coder::array<int, 1U>
    &nnzPr, const ::coder::array<int, 1U> &nrange, int istart, int iend,
    WlsDataStruct *wls, ::coder::array<boolean_T, 1U> &rdTags, int *rdCounts)
  {
    double Ns[64];
    double xsLocal[64];
    int eids[64];
    int gIDs[64];
    int i;
    boolean_T exitg1;

    //  kernel
    i = istart - 1;
    exitg1 = false;
    while ((!exitg1) && (i + 1 <= iend)) {
      double ori;
      int k;
      int nPoints;
      int nid;

      //  get nodal ID
      nid = nrange[i] - 1;
      nPoints = stcls[(stcls.size(1) + stcls.size(1) * (nrange[i] - 1)) - 1] - 1;

      //  fetch coords, global IDs, and then localize coordinates
      for (int j{0}; j <= nPoints; j++) {
        gIDs[j] = stcls[j + stcls.size(1) * nid];
      }

      //  fetch coordinates
      for (k = 0; k <= nPoints; k++) {
        xsLocal[k] = xs[xs.size(1) * (gIDs[k] - 1)];
      }

      //  origin
      ori = xsLocal[0];

      //  localize
      for (k = 0; k <= nPoints; k++) {
        xsLocal[k] -= ori;
      }

      //  compute wls
      b_wls_init(&wls->wlsObj, xsLocal, 64, wls->wlsWgts.name.data,
                 wls->wlsWgts.params_pointwise.data,
                 wls->wlsWgts.params_pointwise.size, degree, wls->interp0,
                 wls->useDag, nPoints + 1);
      if (wls->wlsObj.rank < 0) {
        *rdCounts = -nrange[i];

        //  record error node and return
        exitg1 = true;
      } else {
        if (wls->wlsObj.rank != wls->wlsObj.ncols) {
          (*rdCounts)++;

          //  increment counter for RD nodes
          rdTags[nrange[i] - 1] = true;
        } else {
          int nCells;

          //  get centers and shape function Ns at centers
          nCells = (n2cPtr[nrange[i]] - n2cPtr[nrange[i] - 1]) - 1;
          for (int j{0}; j <= nCells; j++) {
            int eid;
            int npts;
            eid = n2cList[(n2cPtr[nid] + j) - 1] - 1;

            eids[j] = eid + 1;
            for (npts = conn.size(1); conn[(npts + conn.size(1) * eid) - 1] <= 0;
                 npts--) {
            }

            Ns[j] = 1.0 / static_cast<double>(npts);

            //  compute center
            xsLocal[j] = xs[xs.size(1) * (conn[conn.size(1) * eid] - 1)];

            //  C++
            for (int ii{2}; ii <= npts; ii++) {
              xsLocal[j] += xs[xs.size(1) * (conn[(ii + conn.size(1) * eid) - 1]
                - 1)];
            }

            //  localize
            xsLocal[j] = xsLocal[j] * Ns[j] - ori;
          }

          wls_func(&wls->wlsObj, xsLocal, nCells + 1, wls->coeffs);

          //  Compute local OSUS coefficients
          for (k = 0; k <= nPoints; k++) {
            for (int j{0}; j <= nCells; j++) {
              wls->coeffs[j + wls->coeffs.size(1) * k] = wls->coeffs[j +
                wls->coeffs.size(1) * k] * Ns[j];
            }
          }

          //  substract linear interp
          for (int j{0}; j <= nCells; j++) {
            wls->coeffs[j] = wls->coeffs[j] - Ns[j];
          }

          //  kernel for updating crs matrix
          for (int j{0}; j <= nCells; j++) {
            int rowIdx;
            rowIdx = eids[j] - 1;
            for (int ii{0}; ii <= nPoints; ii++) {
              int b_i;
              int c_i;
              int cidx;
              cidx = gIDs[ii];
              b_i = rowPtr[rowIdx] - 1;
              c_i = (rowPtr[rowIdx] + nnzPr[eids[j] - 1]) - 1;
              int exitg2;
              do {
                exitg2 = 0;
                if (b_i + 1 <= c_i) {
                  if (colInd[b_i] == cidx) {
                    k = b_i;
                    exitg2 = 1;
                  } else {
                    b_i++;
                  }
                } else {
                  k = -1;

                  m2cErrMsgIdAndTxt("add_crs:missingIndex",
                                    "could not find column index %d", cidx);
                  exitg2 = 1;
                }
              } while (exitg2 == 0);

              vals[k] = vals[k] + wls->coeffs[j + wls->coeffs.size(1) * ii];
            }
          }
        }

        i++;
      }
    }
  }

  static double find_kth_shortest_dist(::coder::array<double, 1U> &arr, int k,
    int l, int r)
  {
    double dist;
    double val;
    int i;
    int j;

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
      double d;
      double d1;
      int exitg1;
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

  static int g_assemble_body_kernel(const ::coder::array<double, 2U> &xs, const ::
    coder::array<int, 2U> &conn, const ::coder::array<int, 2U> &stcls, const ::
    coder::array<int, 1U> &n2cPtr, const ::coder::array<int, 1U> &n2cList, int
    degree, const ::coder::array<int, 1U> &rowPtr, const ::coder::array<int, 1U>
    &colInd, ::coder::array<double, 1U> &vals, const ::coder::array<int, 1U>
    &nnzPr, const ::coder::array<int, 1U> &nrange, int istart, int iend,
    WlsDataStruct *wls, ::coder::array<boolean_T, 1U> &rdTags)
  {
    double xsLocal[3072];
    double Ns[1024];
    int eids[1024];
    int gIDs[1024];
    int b_iv[2];
    int iv1[2];
    int i;
    int rdCounts;
    boolean_T exitg1;
    rdCounts = 0;

    //  kernel
    i = istart - 1;
    exitg1 = false;
    while ((!exitg1) && (i + 1 <= iend)) {
      double ori_idx_0;
      double ori_idx_1;
      double ori_idx_2;
      int b_xsLocal_tmp;
      int k;
      int nPoints;
      int nid;
      int xsLocal_tmp;

      //  get nodal ID
      nid = nrange[i] - 1;
      nPoints = stcls[(stcls.size(1) + stcls.size(1) * (nrange[i] - 1)) - 1] - 1;

      //  fetch coords, global IDs, and then localize coordinates
      for (int j{0}; j <= nPoints; j++) {
        gIDs[j] = stcls[j + stcls.size(1) * nid];
      }

      //  fetch coordinates
      for (k = 0; k <= nPoints; k++) {
        xsLocal_tmp = gIDs[k];
        xsLocal[3 * k] = xs[xs.size(1) * (xsLocal_tmp - 1)];
        xsLocal[3 * k + 1] = xs[xs.size(1) * (xsLocal_tmp - 1) + 1];
        xsLocal[3 * k + 2] = xs[xs.size(1) * (xsLocal_tmp - 1) + 2];
      }

      //  origin
      ori_idx_0 = xsLocal[0];
      ori_idx_1 = xsLocal[1];
      ori_idx_2 = xsLocal[2];

      //  localize
      for (k = 0; k <= nPoints; k++) {
        xsLocal[3 * k] -= ori_idx_0;
        b_xsLocal_tmp = 3 * k + 1;
        xsLocal[b_xsLocal_tmp] -= ori_idx_1;
        b_xsLocal_tmp = 3 * k + 2;
        xsLocal[b_xsLocal_tmp] -= ori_idx_2;
      }

      //  compute wls
      b_iv[0] = 1024;
      b_iv[1] = 3;
      wls_init(&wls->wlsObj, xsLocal, b_iv, wls->wlsWgts.name.data,
               wls->wlsWgts.params_pointwise.data,
               wls->wlsWgts.params_pointwise.size, degree, wls->interp0,
               wls->useDag, nPoints + 1);
      if (wls->wlsObj.rank < 0) {
        rdCounts = -nrange[i];

        //  record error node and return
        exitg1 = true;
      } else {
        if (wls->wlsObj.rank != wls->wlsObj.ncols) {
          rdCounts++;

          //  increment counter for RD nodes
          rdTags[nrange[i] - 1] = true;
        } else {
          int c_xsLocal_tmp;
          int nCells;

          //  get centers and shape function Ns at centers
          nCells = (n2cPtr[nrange[i]] - n2cPtr[nrange[i] - 1]) - 1;
          for (int j{0}; j <= nCells; j++) {
            int eid;
            int npts;
            eid = n2cList[(n2cPtr[nid] + j) - 1] - 1;

            eids[j] = eid + 1;
            for (npts = conn.size(1); conn[(npts + conn.size(1) * eid) - 1] <= 0;
                 npts--) {
            }

            Ns[j] = 1.0 / static_cast<double>(npts);

            //  compute center
            xsLocal[3 * j] = xs[xs.size(1) * (conn[conn.size(1) * eid] - 1)];
            b_xsLocal_tmp = 3 * j + 1;
            xsLocal[b_xsLocal_tmp] = xs[xs.size(1) * (conn[conn.size(1) * eid] -
              1) + 1];
            c_xsLocal_tmp = 3 * j + 2;
            xsLocal[c_xsLocal_tmp] = xs[xs.size(1) * (conn[conn.size(1) * eid] -
              1) + 2];

            //  C++
            for (int ii{2}; ii <= npts; ii++) {
              xsLocal_tmp = conn[(ii + conn.size(1) * eid) - 1] - 1;
              xsLocal[3 * j] += xs[xs.size(1) * xsLocal_tmp];
              xsLocal[b_xsLocal_tmp] += xs[xs.size(1) * xsLocal_tmp + 1];
              xsLocal[c_xsLocal_tmp] += xs[xs.size(1) * xsLocal_tmp + 2];
            }

            double d;
            double d1;
            double d2;
            d = Ns[j];
            d1 = xsLocal[b_xsLocal_tmp] * d;
            d2 = xsLocal[c_xsLocal_tmp] * d;

            //  localize
            xsLocal[3 * j] = xsLocal[3 * j] * d - ori_idx_0;
            d1 -= ori_idx_1;
            xsLocal[b_xsLocal_tmp] = d1;
            d2 -= ori_idx_2;
            xsLocal[c_xsLocal_tmp] = d2;
          }

          iv1[0] = 1024;
          iv1[1] = 3;
          wls_func(&wls->wlsObj, xsLocal, iv1, nCells + 1, wls->coeffs);

          //  Compute local OSUS coefficients
          for (k = 0; k <= nPoints; k++) {
            for (int j{0}; j <= nCells; j++) {
              wls->coeffs[j + wls->coeffs.size(1) * k] = wls->coeffs[j +
                wls->coeffs.size(1) * k] * Ns[j];
            }
          }

          //  substract linear interp
          for (int j{0}; j <= nCells; j++) {
            wls->coeffs[j] = wls->coeffs[j] - Ns[j];
          }

          //  kernel for updating crs matrix
          for (int j{0}; j <= nCells; j++) {
            int rowIdx;
            rowIdx = eids[j] - 1;
            for (int ii{0}; ii <= nPoints; ii++) {
              int b_i;
              c_xsLocal_tmp = gIDs[ii];
              b_i = rowPtr[rowIdx] - 1;
              xsLocal_tmp = (rowPtr[rowIdx] + nnzPr[eids[j] - 1]) - 1;
              int exitg2;
              do {
                exitg2 = 0;
                if (b_i + 1 <= xsLocal_tmp) {
                  if (colInd[b_i] == c_xsLocal_tmp) {
                    k = b_i;
                    exitg2 = 1;
                  } else {
                    b_i++;
                  }
                } else {
                  k = -1;

                  m2cErrMsgIdAndTxt("add_crs:missingIndex",
                                    "could not find column index %d",
                                    c_xsLocal_tmp);
                  exitg2 = 1;
                }
              } while (exitg2 == 0);

              vals[k] = vals[k] + wls->coeffs[j + wls->coeffs.size(1) * ii];
            }
          }
        }

        i++;
      }
    }

    return rdCounts;
  }

  static void gen_vander(const ::coder::array<double, 2U> &us, int npoints, int
    degree, ::coder::array<double, 2U> &V)
  {
    //  Wrapper function for computing confluent Vandermonde matrix in 1D, 2D, or 3D.
    switch (us.size(1)) {
     case 1:
      {
        int b_n;
        int b_npoints;
        int i;
        int n;
        boolean_T b;
        boolean_T b1;
        b_npoints = npoints - 1;

        //  Generate (confluent) Vandermonde matrix in 1D.
        m2cAssert(us.size(1) == 1, "");

        //  Handle input arguments
        if (npoints == 0) {
          b_npoints = us.size(0) - 1;
        } else {
          m2cAssert(npoints <= us.size(0), "Input us is too small.");
        }

        m2cAssert(degree >= 0, "Degree must be nonnegative");

        //  Number of row blocks
        n = (degree + 1);

        b_n = (us.size(0));
        V.set_size(n, b_n);

        //  Compute rows corresponding to function values
        if (degree != 0) {
          b = true;
          b1 = ((us.size(1) <= 0) || (us.size(0) <= 0));
          n = us.size(1) * us.size(0);
          b_n = 0;
          for (int iPnt{0}; iPnt <= b_npoints; iPnt++) {
            if (b1 || (iPnt >= n)) {
              b_n = 0;
              b = true;
            } else if (b) {
              b = false;
              b_n = iPnt % us.size(0) * us.size(1) + iPnt / us.size(0);
            } else {
              i = us.size(1) * us.size(0) - 1;
              if (b_n > MAX_int32_T - us.size(1)) {
                b_n = iPnt % us.size(0) * us.size(1) + iPnt / us.size(0);
              } else {
                b_n += us.size(1);
                if (b_n > i) {
                  b_n -= i;
                }
              }
            }

            V[iPnt] = 1.0;
            V[iPnt + V.size(1)] = us[b_n];
          }
        } else {
          for (int iPnt{0}; iPnt <= b_npoints; iPnt++) {
            V[iPnt] = 1.0;
          }
        }

        n = degree + 1;
        for (int ii{2}; ii <= n; ii++) {
          b = true;
          b1 = ((us.size(1) <= 0) || (us.size(0) <= 0));
          b_n = us.size(1) * us.size(0);
          i = 0;
          for (int iPnt{0}; iPnt <= b_npoints; iPnt++) {
            if (b1 || (iPnt >= b_n)) {
              i = 0;
              b = true;
            } else if (b) {
              b = false;
              i = iPnt % us.size(0) * us.size(1) + iPnt / us.size(0);
            } else {
              int i1;
              i1 = us.size(1) * us.size(0) - 1;
              if (i > MAX_int32_T - us.size(1)) {
                i = iPnt % us.size(0) * us.size(1) + iPnt / us.size(0);
              } else {
                i += us.size(1);
                if (i > i1) {
                  i -= i1;
                }
              }
            }

            V[iPnt + V.size(1) * (ii - 1)] = V[iPnt + V.size(1) * (ii - 2)] *
              us[i];
          }
        }

        //  Add row blocks corresponding to kth derivatives
      }
      break;

     case 2:
      gen_vander_2d(us, npoints, degree, V);
      break;

     default:
      gen_vander_3d(us, npoints, degree, V);
      break;
    }
  }

  static void gen_vander(const ::coder::array<double, 2U> &us, int npoints, int
    degree, const ::coder::array<double, 1U> &weights, ::coder::array<double, 2U>
    &V)
  {
    //  Wrapper function for computing confluent Vandermonde matrix in 1D, 2D, or 3D.
    switch (us.size(1)) {
     case 1:
      {
        int b_n;
        int i;
        int n;
        boolean_T b;
        boolean_T b1;

        //  Generate (confluent) Vandermonde matrix in 1D.
        m2cAssert(us.size(1) == 1, "");

        //  Handle input arguments
        m2cAssert(npoints <= us.size(0), "Input us is too small.");

        m2cAssert(degree >= 0, "Degree must be nonnegative");

        //  Number of row blocks
        n = (degree + 1);

        b_n = (us.size(0));
        V.set_size(n, b_n);

        //  Compute rows corresponding to function values
        if (weights.size(0) == 0) {
          if (degree != 0) {
            b = true;
            b1 = ((us.size(1) <= 0) || (us.size(0) <= 0));
            n = us.size(1) * us.size(0);
            b_n = 0;
            for (int iPnt{0}; iPnt < npoints; iPnt++) {
              if (b1 || (iPnt >= n)) {
                b_n = 0;
                b = true;
              } else if (b) {
                b = false;
                b_n = iPnt % us.size(0) * us.size(1) + iPnt / us.size(0);
              } else {
                i = us.size(1) * us.size(0) - 1;
                if (b_n > MAX_int32_T - us.size(1)) {
                  b_n = iPnt % us.size(0) * us.size(1) + iPnt / us.size(0);
                } else {
                  b_n += us.size(1);
                  if (b_n > i) {
                    b_n -= i;
                  }
                }
              }

              V[iPnt] = 1.0;
              V[iPnt + V.size(1)] = us[b_n];
            }
          } else {
            for (int iPnt{0}; iPnt < npoints; iPnt++) {
              V[iPnt] = 1.0;
            }
          }
        } else if (degree != 0) {
          b = true;
          b1 = ((us.size(1) <= 0) || (us.size(0) <= 0));
          n = us.size(1) * us.size(0);
          b_n = 0;
          for (int iPnt{0}; iPnt < npoints; iPnt++) {
            if (b1 || (iPnt >= n)) {
              b_n = 0;
              b = true;
            } else if (b) {
              b = false;
              b_n = iPnt % us.size(0) * us.size(1) + iPnt / us.size(0);
            } else {
              i = us.size(1) * us.size(0) - 1;
              if (b_n > MAX_int32_T - us.size(1)) {
                b_n = iPnt % us.size(0) * us.size(1) + iPnt / us.size(0);
              } else {
                b_n += us.size(1);
                if (b_n > i) {
                  b_n -= i;
                }
              }
            }

            V[iPnt] = weights[iPnt];
            V[iPnt + V.size(1)] = us[b_n] * weights[iPnt];
          }
        } else {
          for (int iPnt{0}; iPnt < npoints; iPnt++) {
            V[iPnt] = weights[iPnt];
          }
        }

        n = degree + 1;
        for (int ii{2}; ii <= n; ii++) {
          b = true;
          b1 = ((us.size(1) <= 0) || (us.size(0) <= 0));
          b_n = us.size(1) * us.size(0);
          i = 0;
          for (int iPnt{0}; iPnt < npoints; iPnt++) {
            if (b1 || (iPnt >= b_n)) {
              i = 0;
              b = true;
            } else if (b) {
              b = false;
              i = iPnt % us.size(0) * us.size(1) + iPnt / us.size(0);
            } else {
              int i1;
              i1 = us.size(1) * us.size(0) - 1;
              if (i > MAX_int32_T - us.size(1)) {
                i = iPnt % us.size(0) * us.size(1) + iPnt / us.size(0);
              } else {
                i += us.size(1);
                if (i > i1) {
                  i -= i1;
                }
              }
            }

            V[iPnt + V.size(1) * (ii - 1)] = V[iPnt + V.size(1) * (ii - 2)] *
              us[i];
          }
        }

        //  Add row blocks corresponding to kth derivatives
      }
      break;

     case 2:
      gen_vander_2d(us, npoints, degree, weights, V);
      break;

     default:
      gen_vander_3d(us, npoints, degree, weights, V);
      break;
    }
  }

  static inline
  void gen_vander_1d_dag(int degree, ::coder::array<unsigned char, 1U>
    &dag)
  {
    //  Build a dag for Vandermonde matrix in 1D.
    dag.set_size(degree + 2);
    for (int i{0}; i < degree; i++) {
      dag[i] = 1U;
    }

    dag[degree] = 0U;

    //  a leaf has no child
    dag[dag.size(0) - 1] = static_cast<unsigned char>(degree + 127);
  }

  static void gen_vander_2d(const ::coder::array<double, 2U> &us, int npoints,
    int degree, ::coder::array<double, 2U> &V)
  {
    int b_degree;
    int c;
    int n;

    //  Generate generalized/confluent Vandermonde matrix in 2D.
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

    n = (b_degree);

    b_degree = (us.size(0));
    V.set_size(n, b_degree);

    //  compute 0th order generalized Vandermonde matrix
    if (degree != 0) {
      for (int iPnt{0}; iPnt < npoints; iPnt++) {
        V[iPnt] = 1.0;
        V[iPnt + V.size(1)] = us[us.size(1) * iPnt];
        V[iPnt + V.size(1) * 2] = us[us.size(1) * iPnt + 1];
      }
    } else {
      for (int iPnt{0}; iPnt < npoints; iPnt++) {
        V[iPnt] = 1.0;
      }
    }

    c = 3;
    if (degree < 0) {
      n = -degree;
    } else {
      n = degree;
    }

    for (int deg{2}; deg <= n; deg++) {
      for (int j{0}; j < deg; j++) {
        for (int iPnt{0}; iPnt < npoints; iPnt++) {
          V[iPnt + V.size(1) * c] = V[iPnt + V.size(1) * (c - deg)] * us[us.size
            (1) * iPnt];
        }

        c++;
      }

      for (int iPnt{0}; iPnt < npoints; iPnt++) {
        V[iPnt + V.size(1) * c] = V[iPnt + V.size(1) * ((c - deg) - 1)] *
          us[us.size(1) * iPnt + 1];
      }

      c++;
    }

    //  Compute the bi-degree terms if degree<0
    n = -degree;
    for (int deg{n}; deg >= 1; deg--) {
      for (int k{0}; k < deg; k++) {
        for (int iPnt{0}; iPnt < npoints; iPnt++) {
          V[iPnt + V.size(1) * c] = V[iPnt + V.size(1) * (c - deg)] * us[us.size
            (1) * iPnt];
        }

        c++;
      }
    }

    //  compute higher order confluent Vandermonde matrix blocks incrementally
  }

  static void gen_vander_2d(const ::coder::array<double, 2U> &us, int npoints,
    int degree, const ::coder::array<double, 1U> &weights, ::coder::array<double,
    2U> &V)
  {
    int b_degree;
    int c;
    int n;

    //  Generate generalized/confluent Vandermonde matrix in 2D.
    if (npoints > us.size(0)) {
      m2cErrMsgIdAndTxt("wlslib:BufferTooSmall", "Input us is too small.");
    }

    //  Number of row blocks
    if (degree >= 0) {
      b_degree = (degree + 1) * (degree + 2) / 2;
    } else {
      b_degree = (1 - degree) * (1 - degree);
    }

    n = (b_degree);

    b_degree = (us.size(0));
    V.set_size(n, b_degree);

    //  compute 0th order generalized Vandermonde matrix
    if (weights.size(0) == 0) {
      if (degree != 0) {
        for (int iPnt{0}; iPnt < npoints; iPnt++) {
          V[iPnt] = 1.0;
          V[iPnt + V.size(1)] = us[us.size(1) * iPnt];
          V[iPnt + V.size(1) * 2] = us[us.size(1) * iPnt + 1];
        }
      } else {
        for (int iPnt{0}; iPnt < npoints; iPnt++) {
          V[iPnt] = 1.0;
        }
      }
    } else if (degree != 0) {
      for (int iPnt{0}; iPnt < npoints; iPnt++) {
        V[iPnt] = weights[iPnt];
        V[iPnt + V.size(1)] = us[us.size(1) * iPnt] * weights[iPnt];
        V[iPnt + V.size(1) * 2] = us[us.size(1) * iPnt + 1] * weights[iPnt];
      }
    } else {
      for (int iPnt{0}; iPnt < npoints; iPnt++) {
        V[iPnt] = weights[iPnt];
      }
    }

    c = 3;
    if (degree < 0) {
      n = -degree;
    } else {
      n = degree;
    }

    for (int deg{2}; deg <= n; deg++) {
      for (int j{0}; j < deg; j++) {
        for (int iPnt{0}; iPnt < npoints; iPnt++) {
          V[iPnt + V.size(1) * c] = V[iPnt + V.size(1) * (c - deg)] * us[us.size
            (1) * iPnt];
        }

        c++;
      }

      for (int iPnt{0}; iPnt < npoints; iPnt++) {
        V[iPnt + V.size(1) * c] = V[iPnt + V.size(1) * ((c - deg) - 1)] *
          us[us.size(1) * iPnt + 1];
      }

      c++;
    }

    //  Compute the bi-degree terms if degree<0
    n = -degree;
    for (int deg{n}; deg >= 1; deg--) {
      for (int k{0}; k < deg; k++) {
        for (int iPnt{0}; iPnt < npoints; iPnt++) {
          V[iPnt + V.size(1) * c] = V[iPnt + V.size(1) * (c - deg)] * us[us.size
            (1) * iPnt];
        }

        c++;
      }
    }

    //  compute higher order confluent Vandermonde matrix blocks incrementally
  }

  static void gen_vander_2d_dag(int degree, ::coder::array<unsigned char, 2U>
    &dag)
  {
    int b_degree;
    int c;
    int j;

    //  Build a dag for Vandermonde matrix in 2D.
    if (degree >= 0) {
      b_degree = (degree + 1) * (degree + 2) / 2;
    } else {
      b_degree = (1 - degree) * (1 - degree);
    }

    dag.set_size(b_degree + 1, 2);
    if (degree != 0) {
      dag[0] = 1U;

      //  x-child
      dag[1] = 2U;

      //  y-child
    } else {
      dag[0] = 0U;
      dag[1] = 0U;

      //  No children
    }

    c = 2;
    if (degree < 0) {
      b_degree = -degree;
    } else {
      b_degree = degree;
    }

    for (int deg{2}; deg <= b_degree; deg++) {
      dag[2 * ((c - deg) + 1)] = static_cast<unsigned char>(deg);

      //  x-child
      c += 2;
      for (j = 2; j <= deg; j++) {
        dag[2 * (((c + j) - deg) - 2)] = static_cast<unsigned char>(deg);

        //  x-child
        dag[2 * (((c + j) - deg) - 3) + 1] = static_cast<unsigned char>(deg + 1);

        //  y-child
      }

      c = (c + deg) - 1;
      dag[2 * ((c - deg) - 1) + 1] = static_cast<unsigned char>(deg + 1);

      //  y-child
    }

    //  Set the children of last row to zero
    if (degree > 0) {
      c -= degree;
      for (j = 0; j <= degree; j++) {
        dag[2 * c] = 0U;
        dag[2 * c + 1] = 0U;

        //  no children
        c++;
      }
    } else if (degree < 0) {
      //  Compute the bi-degree terms if degree<0
      b_degree = -degree;
      for (int deg{b_degree}; deg >= 1; deg--) {
        j = c - deg;
        dag[2 * j] = 0U;

        //  no x-child
        dag[2 * j + 1] = static_cast<unsigned char>(deg + 1);

        //  y-child
        dag[2 * (j + 1)] = static_cast<unsigned char>(deg);

        //  x-child
        c++;
        for (int k{2}; k <= deg; k++) {
          j = (c + k) - deg;
          dag[2 * (j - 1)] = static_cast<unsigned char>(deg);

          //  x-child
          dag[2 * (j - 2) + 1] = static_cast<unsigned char>(deg + 1);

          //  y-child
        }

        c = (c + deg) - 1;
        dag[2 * (c - deg) + 1] = 0U;

        //  no y-child
      }

      dag[2 * c] = 0U;
      dag[2 * c + 1] = 0U;
    }

    //  Use last entry as signature
    b_degree = (dag.size(0) << 1) - 1;
    dag[((b_degree % dag.size(0)) << 1) + b_degree / dag.size(0)] = static_cast<
      unsigned char>(degree + 127);
  }

  static void gen_vander_3d(const ::coder::array<double, 2U> &us, int npoints,
    int degree, ::coder::array<double, 2U> &V)
  {
    int b_degree;
    int c;
    int d;
    int deg;
    int n;

    //  Generate generalized/confluent Vandermonde matrix in 3D.
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

    b_degree = (b_degree);

    n = (us.size(0));
    V.set_size(b_degree, n);

    //  compute 0th order generalized Vandermonde matrix
    if (degree != 0) {
      for (int iPnt{0}; iPnt < npoints; iPnt++) {
        V[iPnt] = 1.0;
        V[iPnt + V.size(1)] = us[us.size(1) * iPnt];
        V[iPnt + V.size(1) * 2] = us[us.size(1) * iPnt + 1];
        V[iPnt + V.size(1) * 3] = us[us.size(1) * iPnt + 2];
      }
    } else {
      for (int iPnt{0}; iPnt < npoints; iPnt++) {
        V[iPnt] = 1.0;
      }
    }

    c = 4;
    d = 3;
    if (degree < 0) {
      n = -degree;
    } else {
      n = degree;
    }

    for (deg = 2; deg <= n; deg++) {
      //  Within each level, use convention of Pascal triangle with x^deg at peak
      for (int j{0}; j < deg; j++) {
        for (int iPnt{0}; iPnt < npoints; iPnt++) {
          V[iPnt + V.size(1) * c] = V[iPnt + V.size(1) * (c - d)] * us[us.size(1)
            * iPnt];
        }

        c++;
      }

      for (int iPnt{0}; iPnt < npoints; iPnt++) {
        V[iPnt + V.size(1) * c] = V[iPnt + V.size(1) * ((c - d) - 1)] *
          us[us.size(1) * iPnt + 1];
      }

      c++;
      for (int j{0}; j < d; j++) {
        for (int iPnt{0}; iPnt < npoints; iPnt++) {
          V[iPnt + V.size(1) * c] = V[iPnt + V.size(1) * (((c - d) - deg) - 1)] *
            us[us.size(1) * iPnt + 2];
        }

        c++;
      }

      d = (d + deg) + 1;
    }

    //  Compute the tri-degree terms if degree<0
    if (degree < 0) {
      int cornerTriangle;
      int excess;
      int maxLayers;
      int nTermsInLayer;
      deg = -degree;
      maxLayers = -degree * 3;

      // max number of layers needed in the Pascal tetrahedron
      cornerTriangle = 0;

      // number of elements subtracted in each corner Pascal triangle
      nTermsInLayer = d;

      // initializing number of elements in layer
      excess = 0;

      // excess based on overlapping of growing Pascal triangles
      n = 1 - degree;
      for (int p{n}; p <= maxLayers; p++) {
        int counterBottomRow;
        int gap;
        int nTermsInPrevLayer;

        //  Within each level, x^deg is at the peak of Pascal triangle
        cornerTriangle = (cornerTriangle + p) + degree;
        counterBottomRow = 1;

        // counter for the bottom row to be subtracted later
        for (int k{0}; k < deg; k++) {
          for (int iPnt{0}; iPnt < npoints; iPnt++) {
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
        for (int j{0}; j <= b_degree; j++) {
          for (int iPnt{0}; iPnt < npoints; iPnt++) {
            V[iPnt + V.size(1) * c] = V[iPnt + V.size(1) * (c - gap)] *
              us[us.size(1) * iPnt + 2];
          }

          c++;
        }
      }
    }

    m2cAssert(true, "");
  }

  static void gen_vander_3d(const ::coder::array<double, 2U> &us, int npoints,
    int degree, const ::coder::array<double, 1U> &weights, ::coder::array<double,
    2U> &V)
  {
    int b_degree;
    int c;
    int d;
    int deg;
    int n;

    //  Generate generalized/confluent Vandermonde matrix in 3D.
    if (npoints > us.size(0)) {
      m2cErrMsgIdAndTxt("wlslib:BufferTooSmall", "Input us is too small.");
    }

    //  Allocate storage for V
    if (degree >= 0) {
      b_degree = (degree + 1) * (degree + 2) * (degree + 3) / 6;
    } else {
      b_degree = (1 - degree) * (1 - degree) * (1 - degree);
    }

    b_degree = (b_degree);

    n = (us.size(0));
    V.set_size(b_degree, n);

    //  compute 0th order generalized Vandermonde matrix
    if (weights.size(0) == 0) {
      if (degree != 0) {
        for (int iPnt{0}; iPnt < npoints; iPnt++) {
          V[iPnt] = 1.0;
          V[iPnt + V.size(1)] = us[us.size(1) * iPnt];
          V[iPnt + V.size(1) * 2] = us[us.size(1) * iPnt + 1];
          V[iPnt + V.size(1) * 3] = us[us.size(1) * iPnt + 2];
        }
      } else {
        for (int iPnt{0}; iPnt < npoints; iPnt++) {
          V[iPnt] = 1.0;
        }
      }
    } else if (degree != 0) {
      for (int iPnt{0}; iPnt < npoints; iPnt++) {
        V[iPnt] = weights[iPnt];
        V[iPnt + V.size(1)] = us[us.size(1) * iPnt] * weights[iPnt];
        V[iPnt + V.size(1) * 2] = us[us.size(1) * iPnt + 1] * weights[iPnt];
        V[iPnt + V.size(1) * 3] = us[us.size(1) * iPnt + 2] * weights[iPnt];
      }
    } else {
      for (int iPnt{0}; iPnt < npoints; iPnt++) {
        V[iPnt] = weights[iPnt];
      }
    }

    c = 4;
    d = 3;
    if (degree < 0) {
      n = -degree;
    } else {
      n = degree;
    }

    for (deg = 2; deg <= n; deg++) {
      //  Within each level, use convention of Pascal triangle with x^deg at peak
      for (int j{0}; j < deg; j++) {
        for (int iPnt{0}; iPnt < npoints; iPnt++) {
          V[iPnt + V.size(1) * c] = V[iPnt + V.size(1) * (c - d)] * us[us.size(1)
            * iPnt];
        }

        c++;
      }

      for (int iPnt{0}; iPnt < npoints; iPnt++) {
        V[iPnt + V.size(1) * c] = V[iPnt + V.size(1) * ((c - d) - 1)] *
          us[us.size(1) * iPnt + 1];
      }

      c++;
      for (int j{0}; j < d; j++) {
        for (int iPnt{0}; iPnt < npoints; iPnt++) {
          V[iPnt + V.size(1) * c] = V[iPnt + V.size(1) * (((c - d) - deg) - 1)] *
            us[us.size(1) * iPnt + 2];
        }

        c++;
      }

      d = (d + deg) + 1;
    }

    //  Compute the tri-degree terms if degree<0
    if (degree < 0) {
      int cornerTriangle;
      int excess;
      int maxLayers;
      int nTermsInLayer;
      deg = -degree;
      maxLayers = -degree * 3;

      // max number of layers needed in the Pascal tetrahedron
      cornerTriangle = 0;

      // number of elements subtracted in each corner Pascal triangle
      nTermsInLayer = d;

      // initializing number of elements in layer
      excess = 0;

      // excess based on overlapping of growing Pascal triangles
      n = 1 - degree;
      for (int p{n}; p <= maxLayers; p++) {
        int counterBottomRow;
        int gap;
        int nTermsInPrevLayer;

        //  Within each level, x^deg is at the peak of Pascal triangle
        cornerTriangle = (cornerTriangle + p) + degree;
        counterBottomRow = 1;

        // counter for the bottom row to be subtracted later
        for (int k{0}; k < deg; k++) {
          for (int iPnt{0}; iPnt < npoints; iPnt++) {
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
        for (int j{0}; j <= b_degree; j++) {
          for (int iPnt{0}; iPnt < npoints; iPnt++) {
            V[iPnt + V.size(1) * c] = V[iPnt + V.size(1) * (c - gap)] *
              us[us.size(1) * iPnt + 2];
          }

          c++;
        }
      }
    }

    m2cAssert(true, "");
  }

  static void gen_vander_3d_dag(int degree, ::coder::array<unsigned char, 2U>
    &dag)
  {
    int b_degree;
    int c;
    int d;
    int i;
    int maxterms;

    //  Build a DAG for Vandermonde matrix in 3D.
    if (degree >= 0) {
      b_degree = (degree + 1) * (degree + 2) * (degree + 3) / 6;
    } else {
      b_degree = (1 - degree) * (1 - degree) * (1 - degree);
    }

    dag.set_size(b_degree + 1, 3);
    if (degree != 0) {
      dag[0] = 1U;

      //  x-child
      dag[1] = 2U;

      //  y-child
      dag[2] = 3U;

      //  z-child
    } else {
      dag[0] = 0U;
      dag[1] = 0U;
      dag[2] = 0U;

      //  No children
    }

    c = 1;
    d = 3;
    if (degree < 0) {
      i = -degree;
    } else {
      i = degree;
    }

    for (int deg{2}; deg <= i; deg++) {
      maxterms = (d + deg) + 1;
      for (int j{deg}; j >= 1; j--) {
        for (int b_i{0}; b_i < j; b_i++) {
          b_degree = c + b_i;
          dag[3 * b_degree] = static_cast<unsigned char>(d);
          dag[3 * b_degree + 1] = static_cast<unsigned char>(d + 1);
          dag[3 * b_degree + 2] = static_cast<unsigned char>(maxterms);
        }

        c += j;
        d++;
      }

      d = maxterms;
    }

    if (degree > 0) {
      i = dag.size(0);
      for (int b_i{c + 1}; b_i <= i; b_i++) {
        dag[3 * (b_i - 1)] = 0U;
        dag[3 * (b_i - 1) + 1] = 0U;
        dag[3 * (b_i - 1) + 2] = 0U;
      }
    } else if (degree < 0) {
      int cornerTriangle;
      int excess;
      int maxlayers;
      int num_elem_group;
      maxlayers = -3 * degree + 1;
      cornerTriangle = 0;

      // number of elements subtracted in each corner Pascal triangle
      excess = 0;

      // excess based on overlapping of growing Pascal triangles
      num_elem_group = -1;
      i = -degree;
      for (int p{i}; p <= maxlayers; p++) {
        int ntermsinlayer;
        int x_tmp;
        int y;
        cornerTriangle = (cornerTriangle + p) + degree;
        y = excess << 1;
        x_tmp = p + degree;
        b_degree = (x_tmp << 1) - p;
        if (b_degree < 0) {
          b_degree = 0;
        }

        excess += b_degree;
        maxterms = ((((d + (-degree << 1)) - 3 * cornerTriangle) - p) + y) +
          excess;
        ntermsinlayer = d + 3 * (excess - cornerTriangle);
        for (int group{0}; group <= i; group++) {
          y = x_tmp - group;
          if (y < 0) {
            b_degree = -y;
          } else {
            b_degree = y;
          }

          num_elem_group = -degree - b_degree;
          if (num_elem_group + 1 < 1) {
            ntermsinlayer -= 3;
          } else if (((-degree - p) + group) - 1 < 0) {
            dag[3 * c] = 0U;
            for (int b_i{0}; b_i < num_elem_group; b_i++) {
              b_degree = c + b_i;
              dag[3 * (b_degree + 1)] = static_cast<unsigned char>(ntermsinlayer
                - 1);
              dag[3 * b_degree + 1] = static_cast<unsigned char>(ntermsinlayer);
              dag[3 * b_degree + 2] = static_cast<unsigned char>(maxterms);
            }

            c += num_elem_group;
            dag[3 * c + 1] = 0U;
            dag[3 * c + 2] = static_cast<unsigned char>(maxterms);
            c++;
            if (y > 1) {
              y = 1;
            }

            ntermsinlayer -= y;
          } else {
            for (int b_i{0}; b_i <= num_elem_group; b_i++) {
              b_degree = c + b_i;
              dag[3 * b_degree] = static_cast<unsigned char>(ntermsinlayer - 1);
              dag[3 * b_degree + 1] = static_cast<unsigned char>(ntermsinlayer);
              dag[3 * b_degree + 2] = static_cast<unsigned char>(maxterms);
            }

            c = (c + num_elem_group) + 1;
            ntermsinlayer++;
          }
        }

        for (int j{0}; j <= num_elem_group; j++) {
          dag[3 * (((c + j) - num_elem_group) - 1) + 2] = 0U;
        }

        d = (d + p) + 2;
      }
    }

    //  Use last entry as signature
    i = dag.size(0) * 3 - 1;
    dag[i % dag.size(0) * 3 + i / dag.size(0)] = static_cast<unsigned char>
      (degree + 127);
  }

  static int h_assemble_body_kernel(const ::coder::array<double, 2U> &xs, const ::
    coder::array<int, 2U> &conn, const ::coder::array<int, 2U> &stcls, const ::
    coder::array<int, 1U> &n2cPtr, const ::coder::array<int, 1U> &n2cList, int
    degree, const ::coder::array<int, 1U> &rowPtr, const ::coder::array<int, 1U>
    &colInd, ::coder::array<double, 1U> &vals, const ::coder::array<int, 1U>
    &nnzPr, const ::coder::array<int, 1U> &nrange, int istart, int iend,
    WlsDataStruct *wls, ::coder::array<boolean_T, 1U> &rdTags)
  {
    double xsLocal[1024];
    double Ns[512];
    int eids[512];
    int gIDs[512];
    int b_iv[2];
    int iv1[2];
    int i;
    int rdCounts;
    boolean_T exitg1;
    rdCounts = 0;

    //  kernel
    i = istart - 1;
    exitg1 = false;
    while ((!exitg1) && (i + 1 <= iend)) {
      double ori_idx_0;
      double ori_idx_1;
      int b_i;
      int k;
      int nPoints;
      int nid;
      int xsLocal_tmp;

      //  get nodal ID
      nid = nrange[i] - 1;
      nPoints = stcls[(stcls.size(1) + stcls.size(1) * (nrange[i] - 1)) - 1] - 1;

      //  fetch coords, global IDs, and then localize coordinates
      for (int j{0}; j <= nPoints; j++) {
        gIDs[j] = stcls[j + stcls.size(1) * nid];
      }

      //  fetch coordinates
      for (k = 0; k <= nPoints; k++) {
        xsLocal_tmp = k << 1;
        b_i = gIDs[k];
        xsLocal[xsLocal_tmp] = xs[xs.size(1) * (b_i - 1)];
        xsLocal[xsLocal_tmp + 1] = xs[xs.size(1) * (b_i - 1) + 1];
      }

      //  origin
      ori_idx_0 = xsLocal[0];
      ori_idx_1 = xsLocal[1];

      //  localize
      for (k = 0; k <= nPoints; k++) {
        xsLocal_tmp = k << 1;
        xsLocal[xsLocal_tmp] -= ori_idx_0;
        xsLocal[xsLocal_tmp + 1] -= ori_idx_1;
      }

      //  compute wls
      b_iv[0] = 512;
      b_iv[1] = 2;
      wls_init(&wls->wlsObj, xsLocal, b_iv, wls->wlsWgts.name.data,
               wls->wlsWgts.params_pointwise.data,
               wls->wlsWgts.params_pointwise.size, degree, wls->interp0,
               wls->useDag, nPoints + 1);
      if (wls->wlsObj.rank < 0) {
        rdCounts = -nrange[i];

        //  record error node and return
        exitg1 = true;
      } else {
        if (wls->wlsObj.rank != wls->wlsObj.ncols) {
          rdCounts++;

          //  increment counter for RD nodes
          rdTags[nrange[i] - 1] = true;
        } else {
          int b_xsLocal_tmp;
          int nCells;

          //  get centers and shape function Ns at centers
          nCells = (n2cPtr[nrange[i]] - n2cPtr[nrange[i] - 1]) - 1;
          for (int j{0}; j <= nCells; j++) {
            int eid;
            int npts;
            eid = n2cList[(n2cPtr[nid] + j) - 1] - 1;

            eids[j] = eid + 1;
            for (npts = conn.size(1); conn[(npts + conn.size(1) * eid) - 1] <= 0;
                 npts--) {
            }

            Ns[j] = 1.0 / static_cast<double>(npts);

            //  compute center
            xsLocal_tmp = j << 1;
            xsLocal[xsLocal_tmp] = xs[xs.size(1) * (conn[conn.size(1) * eid] - 1)];
            xsLocal[xsLocal_tmp + 1] = xs[xs.size(1) * (conn[conn.size(1) * eid]
              - 1) + 1];

            //  C++
            for (int ii{2}; ii <= npts; ii++) {
              b_xsLocal_tmp = conn[(ii + conn.size(1) * eid) - 1] - 1;
              xsLocal[xsLocal_tmp] += xs[xs.size(1) * b_xsLocal_tmp];
              xsLocal[xsLocal_tmp + 1] += xs[xs.size(1) * b_xsLocal_tmp + 1];
            }

            double d;
            double d1;
            d = Ns[j];
            d1 = xsLocal[xsLocal_tmp + 1] * d;

            //  localize
            xsLocal[xsLocal_tmp] = xsLocal[xsLocal_tmp] * d - ori_idx_0;
            d1 -= ori_idx_1;
            xsLocal[xsLocal_tmp + 1] = d1;
          }

          iv1[0] = 512;
          iv1[1] = 2;
          wls_func(&wls->wlsObj, xsLocal, iv1, nCells + 1, wls->coeffs);

          //  Compute local OSUS coefficients
          for (k = 0; k <= nPoints; k++) {
            for (int j{0}; j <= nCells; j++) {
              wls->coeffs[j + wls->coeffs.size(1) * k] = wls->coeffs[j +
                wls->coeffs.size(1) * k] * Ns[j];
            }
          }

          //  substract linear interp
          for (int j{0}; j <= nCells; j++) {
            wls->coeffs[j] = wls->coeffs[j] - Ns[j];
          }

          //  kernel for updating crs matrix
          for (int j{0}; j <= nCells; j++) {
            int rowIdx;
            rowIdx = eids[j] - 1;
            for (int ii{0}; ii <= nPoints; ii++) {
              int c_i;
              b_xsLocal_tmp = gIDs[ii];
              c_i = rowPtr[rowIdx] - 1;
              b_i = (rowPtr[rowIdx] + nnzPr[eids[j] - 1]) - 1;
              int exitg2;
              do {
                exitg2 = 0;
                if (c_i + 1 <= b_i) {
                  if (colInd[c_i] == b_xsLocal_tmp) {
                    k = c_i;
                    exitg2 = 1;
                  } else {
                    c_i++;
                  }
                } else {
                  k = -1;

                  m2cErrMsgIdAndTxt("add_crs:missingIndex",
                                    "could not find column index %d",
                                    b_xsLocal_tmp);
                  exitg2 = 1;
                }
              } while (exitg2 == 0);

              vals[k] = vals[k] + wls->coeffs[j + wls->coeffs.size(1) * ii];
            }
          }
        }

        i++;
      }
    }

    return rdCounts;
  }

  static int i_assemble_body_kernel(const ::coder::array<double, 2U> &xs, const ::
    coder::array<int, 2U> &conn, const ::coder::array<int, 2U> &stcls, const ::
    coder::array<int, 1U> &n2cPtr, const ::coder::array<int, 1U> &n2cList, int
    degree, const ::coder::array<int, 1U> &rowPtr, const ::coder::array<int, 1U>
    &colInd, ::coder::array<double, 1U> &vals, const ::coder::array<int, 1U>
    &nnzPr, const ::coder::array<int, 1U> &nrange, int istart, int iend,
    WlsDataStruct *wls, ::coder::array<boolean_T, 1U> &rdTags)
  {
    double Ns[64];
    double xsLocal[64];
    int eids[64];
    int gIDs[64];
    int i;
    int rdCounts;
    boolean_T exitg1;
    rdCounts = 0;

    //  kernel
    i = istart - 1;
    exitg1 = false;
    while ((!exitg1) && (i + 1 <= iend)) {
      double ori;
      int k;
      int nPoints;
      int nid;

      //  get nodal ID
      nid = nrange[i] - 1;
      nPoints = stcls[(stcls.size(1) + stcls.size(1) * (nrange[i] - 1)) - 1] - 1;

      //  fetch coords, global IDs, and then localize coordinates
      for (int j{0}; j <= nPoints; j++) {
        gIDs[j] = stcls[j + stcls.size(1) * nid];
      }

      //  fetch coordinates
      for (k = 0; k <= nPoints; k++) {
        xsLocal[k] = xs[xs.size(1) * (gIDs[k] - 1)];
      }

      //  origin
      ori = xsLocal[0];

      //  localize
      for (k = 0; k <= nPoints; k++) {
        xsLocal[k] -= ori;
      }

      //  compute wls
      b_wls_init(&wls->wlsObj, xsLocal, 64, wls->wlsWgts.name.data,
                 wls->wlsWgts.params_pointwise.data,
                 wls->wlsWgts.params_pointwise.size, degree, wls->interp0,
                 wls->useDag, nPoints + 1);
      if (wls->wlsObj.rank < 0) {
        rdCounts = -nrange[i];

        //  record error node and return
        exitg1 = true;
      } else {
        if (wls->wlsObj.rank != wls->wlsObj.ncols) {
          rdCounts++;

          //  increment counter for RD nodes
          rdTags[nrange[i] - 1] = true;
        } else {
          int nCells;

          //  get centers and shape function Ns at centers
          nCells = (n2cPtr[nrange[i]] - n2cPtr[nrange[i] - 1]) - 1;
          for (int j{0}; j <= nCells; j++) {
            int eid;
            int npts;
            eid = n2cList[(n2cPtr[nid] + j) - 1] - 1;

            eids[j] = eid + 1;
            for (npts = conn.size(1); conn[(npts + conn.size(1) * eid) - 1] <= 0;
                 npts--) {
            }

            Ns[j] = 1.0 / static_cast<double>(npts);

            //  compute center
            xsLocal[j] = xs[xs.size(1) * (conn[conn.size(1) * eid] - 1)];

            //  C++
            for (int ii{2}; ii <= npts; ii++) {
              xsLocal[j] += xs[xs.size(1) * (conn[(ii + conn.size(1) * eid) - 1]
                - 1)];
            }

            //  localize
            xsLocal[j] = xsLocal[j] * Ns[j] - ori;
          }

          wls_func(&wls->wlsObj, xsLocal, nCells + 1, wls->coeffs);

          //  Compute local OSUS coefficients
          for (k = 0; k <= nPoints; k++) {
            for (int j{0}; j <= nCells; j++) {
              wls->coeffs[j + wls->coeffs.size(1) * k] = wls->coeffs[j +
                wls->coeffs.size(1) * k] * Ns[j];
            }
          }

          //  substract linear interp
          for (int j{0}; j <= nCells; j++) {
            wls->coeffs[j] = wls->coeffs[j] - Ns[j];
          }

          //  kernel for updating crs matrix
          for (int j{0}; j <= nCells; j++) {
            int rowIdx;
            rowIdx = eids[j] - 1;
            for (int ii{0}; ii <= nPoints; ii++) {
              int b_i;
              int c_i;
              int cidx;
              cidx = gIDs[ii];
              b_i = rowPtr[rowIdx] - 1;
              c_i = (rowPtr[rowIdx] + nnzPr[eids[j] - 1]) - 1;
              int exitg2;
              do {
                exitg2 = 0;
                if (b_i + 1 <= c_i) {
                  if (colInd[b_i] == cidx) {
                    k = b_i;
                    exitg2 = 1;
                  } else {
                    b_i++;
                  }
                } else {
                  k = -1;

                  m2cErrMsgIdAndTxt("add_crs:missingIndex",
                                    "could not find column index %d", cidx);
                  exitg2 = 1;
                }
              } while (exitg2 == 0);

              vals[k] = vals[k] + wls->coeffs[j + wls->coeffs.size(1) * ii];
            }
          }
        }

        i++;
      }
    }

    return rdCounts;
  }

  static void init_osusop(const ::coder::array<int, 2U> &conn, const ::coder::
    array<int, 2U> &stcls, int maxNnzPr, ::coder::array<int, 1U> &rowPtr, ::
    coder::array<int, 1U> &colInd, ::coder::array<double, 1U> &vals, ::coder::
    array<int, 1U> &nnzPr)
  {
    ::coder::array<boolean_T, 1U> visited_;
    int nodes[512];
    int i;
    int loop_ub;
    int m;

    //  function for allocating uncompressed CRS matrix
    m = conn.size(0);

    rowPtr.set_size(conn.size(0) + 1);
    rowPtr[0] = 1;
    colInd.set_size(static_cast<int>(std::round(static_cast<double>(conn.size(0)
      * maxNnzPr) * 1.25)));
    vals.set_size(colInd.size(0));

    //  allocate for nnz per row
    nnzPr.set_size(conn.size(0));

    //  buffers
    visited_.set_size(stcls.size(0));
    loop_ub = stcls.size(0);
    for (i = 0; i < loop_ub; i++) {
      visited_[i] = false;
    }

    for (int e{0}; e < m; e++) {
      int j;
      int npts;
      for (npts = conn.size(1) - 1; conn[npts + conn.size(1) * e] <= 0; npts--)
      {
      }

      j = -1;

      //  number of unique points of all stencils associated e
      for (int b_i{0}; b_i <= npts; b_i++) {
        int n0;
        int nid;
        nid = conn[b_i + conn.size(1) * e] - 1;

        //  get stencil size
        n0 = stcls[(stcls.size(1) + stcls.size(1) * nid) - 1];
        for (int k{0}; k < n0; k++) {
          if (!visited_[stcls[k + stcls.size(1) * nid] - 1]) {
            j++;
            nodes[j] = stcls[k + stcls.size(1) * nid];
            visited_[stcls[k + stcls.size(1) * nid] - 1] = true;
          }
        }
      }

      //  reset visited nodes
      for (int b_i{0}; b_i <= j; b_i++) {
        visited_[nodes[b_i] - 1] = false;
      }

      //  actual nnz
      nnzPr[e] = j + 1;

      //  with extra reserved space
      i = rowPtr[e] + static_cast<int>(std::ceil(1.25 * static_cast<double>(j +
        1)));
      rowPtr[e + 1] = i;
      if (i - 1 > colInd.size(0)) {
        int exSpace;

        //  helper function to enlarge an array
        exSpace = static_cast<int>(std::round(static_cast<double>(colInd.size(0))
          * 0.2));
        if (colInd.size(0) + exSpace < i - 1) {
          exSpace = (i - colInd.size(0)) - 1;
        }

        i = colInd.size(0);
        colInd.set_size(colInd.size(0) + exSpace);
        for (loop_ub = 0; loop_ub < exSpace; loop_ub++) {
          colInd[i + loop_ub] = 0;
        }

        //  helper function to enlarge an array
        exSpace = static_cast<int>(std::round(static_cast<double>(vals.size(0)) *
          0.2));
        i = rowPtr[e + 1];
        if (vals.size(0) + exSpace < i - 1) {
          exSpace = (i - vals.size(0)) - 1;
        }

        i = vals.size(0);
        vals.set_size(vals.size(0) + exSpace);
        for (loop_ub = 0; loop_ub < exSpace; loop_ub++) {
          vals[i + loop_ub] = 0.0;
        }
      }

      for (int k{0}; k <= j; k++) {
        colInd[(rowPtr[e] + k) - 1] = nodes[k];
      }
    }

    //  zeros out vals
    loop_ub = vals.size(0);
    vals.set_size(loop_ub);
    for (i = 0; i < loop_ub; i++) {
      vals[i] = 0.0;
    }
  }

  static void mark_kernel(int n, double hGlobal, double cGlobal, double cLocal,
    double kappa1, double kappa0, const ::coder::array<double, 2U> &fs, const ::
    coder::array<double, 2U> &alphaCell, const ::coder::array<double, 2U> &beta,
    const ::coder::array<double, 2U> &dfGlobal, const ::coder::array<int, 2U>
    &mesh_conn, const ::coder::array<double, 1U> &mesh_cellSizes, const ::coder::
    array<int, 1U> &mesh_n2cPtr, const ::coder::array<int, 1U> &mesh_n2cList,
    int nRhs, ::coder::array<signed char, 2U> &disTags)
  {
    double thres;
    int iend;
    int istart;
    int nthreads;

#pragma omp single

    { //single
      disTags.set_size(n, nRhs);

    } //single

    //  Obtains starting and ending indices of local chunk for current thread
    nthreads = 1;

#ifdef _OPENMP

    nthreads = omp_get_num_threads();

#endif //_OPENMP

    if (nthreads == 1) {
      istart = 1;
      iend = n;
    } else {
      int b_remainder;
      int chunk;
      int threadID;
      int u1;

      //  Gets the thread number of the thread, within the team, making this call.
      threadID = 0;

#ifdef _OPENMP

      threadID = omp_get_thread_num();

#endif //_OPENMP

      chunk = n / nthreads;
      b_remainder = n - nthreads * chunk;
      u1 = threadID;
      if (b_remainder <= threadID) {
        u1 = b_remainder;
      }

      istart = (threadID * chunk + u1) + 1;
      iend = ((istart + chunk) + (threadID < b_remainder)) - 1;
    }

    thres = cGlobal * std::pow(hGlobal, 1.5);

    //  C++
    for (int i{istart}; i <= iend; i++) {
      for (int k{0}; k < nRhs; k++) {
        double tauGlobal;
        int j;
        boolean_T disCell;
        boolean_T exitg1;
        disTags[k + disTags.size(1) * (i - 1)] = 0;
        tauGlobal = thres * dfGlobal[k];
        disCell = false;
        j = mesh_n2cPtr[i - 1] - 1;
        exitg1 = false;
        while ((!exitg1) && (j + 1 <= mesh_n2cPtr[i] - 1)) {
          double fMax;
          double fMin;
          int eid;
          int npts;
          eid = mesh_n2cList[j] - 1;
          for (npts = mesh_conn.size(1) - 1; mesh_conn[npts + mesh_conn.size(1) *
               eid] <= 0; npts--) {
          }

          fMax = -1.7976931348623157E+308;
          fMin = 1.7976931348623157E+308;
          for (int ii{0}; ii <= npts; ii++) {
            double fValue;
            fValue = fs[k + fs.size(1) * (mesh_conn[ii + mesh_conn.size(1) * eid]
              - 1)];
            if (fMax < fValue) {
              fMax = fValue;
            }

            if (fMin > fValue) {
              fMin = fValue;
            }
          }

          //  compute local df and thres
          if (std::abs(alphaCell[k + alphaCell.size(1) * (mesh_n2cList[j] - 1)])
              > std::fmax(tauGlobal, cLocal * (fMax - fMin) * std::sqrt
                          (mesh_cellSizes[mesh_n2cList[j] - 1]))) {
            //  dis cell
            disCell = true;
            exitg1 = true;
          } else {
            j++;
          }
        }

        if (disCell) {
          double d;
          d = beta[k + beta.size(1) * (i - 1)];
          if (d > kappa1) {
            disTags[k + disTags.size(1) * (i - 1)] = 2;

            //  C1
            if (d > kappa0) {
              disTags[k + disTags.size(1) * (i - 1)] = 1;
            }

            //  C0
          }
        }
      }
    }
  }

  static void mark_kernel(int n, double hGlobal, double cGlobal, double cLocal,
    double kappa1, double kappa0, const ::coder::array<double, 2U> &fs, const ::
    coder::array<double, 2U> &alphaCell, const ::coder::array<double, 2U> &beta,
    const ::coder::array<double, 2U> &dfGlobal, const ::coder::array<int, 2U>
    &mesh_conn, const ::coder::array<double, 1U> &mesh_cellSizes, const ::coder::
    array<int, 1U> &mesh_n2cPtr, const ::coder::array<int, 1U> &mesh_n2cList,
    const ::coder::array<int, 1U> &mesh_n2nPtr, const ::coder::array<int, 1U>
    &mesh_n2nList, int nRhs, ::coder::array<signed char, 2U> &disTags)
  {
    double thres;
    disTags.set_size(n, nRhs);
    thres = cGlobal * std::pow(hGlobal, 1.5);

    //  C++
    for (int i{0}; i < n; i++) {
      for (int k{0}; k < nRhs; k++) {
        disTags[k + disTags.size(1) * i] = 0;
      }
    }

    for (int i{0}; i < n; i++) {
      for (int k{0}; k < nRhs; k++) {
        double tauGlobal;
        int j;
        boolean_T disCell;
        boolean_T exitg1;
        tauGlobal = thres * dfGlobal[k];
        disCell = false;
        j = mesh_n2cPtr[i] - 1;
        exitg1 = false;
        while ((!exitg1) && (j + 1 <= mesh_n2cPtr[i + 1] - 1)) {
          double fMax;
          double fMin;
          int eid;
          int npts;
          eid = mesh_n2cList[j] - 1;
          for (npts = mesh_conn.size(1) - 1; mesh_conn[npts + mesh_conn.size(1) *
               eid] <= 0; npts--) {
          }

          fMax = -1.7976931348623157E+308;
          fMin = 1.7976931348623157E+308;
          for (int ii{0}; ii <= npts; ii++) {
            double fValue;
            fValue = fs[k + fs.size(1) * (mesh_conn[ii + mesh_conn.size(1) * eid]
              - 1)];
            if (fMax < fValue) {
              fMax = fValue;
            }

            if (fMin > fValue) {
              fMin = fValue;
            }
          }

          //  compute local df and thres
          if (std::abs(alphaCell[k + alphaCell.size(1) * (mesh_n2cList[j] - 1)])
              > std::fmax(tauGlobal, cLocal * (fMax - fMin) * std::sqrt
                          (mesh_cellSizes[mesh_n2cList[j] - 1]))) {
            //  dis cell
            disCell = true;
            exitg1 = true;
          } else {
            j++;
          }
        }

        if (disCell) {
          double d;
          d = beta[k + beta.size(1) * i];
          if (d > kappa1) {
            int b_i;
            int i1;
            signed char tag;
            tag = 2;

            //  C1
            if (d > kappa0) {
              tag = 1;
            }

            //  C0
            disTags[k + disTags.size(1) * i] = tag;

            //  if we extend to 1-ring neighborhood
            b_i = mesh_n2nPtr[i];
            i1 = mesh_n2nPtr[i + 1] - 1;
            for (j = b_i; j <= i1; j++) {
              int i2;
              i2 = mesh_n2nList[j - 1] - 1;
              if (disTags[k + disTags.size(1) * i2] != 1) {
                disTags[k + disTags.size(1) * i2] = tag;
              }
            }
          }
        }
      }
    }
  }

  static void omp4mRecurPartMesh(int n, const ::coder::array<int, 2U> &cells,
    int dim, int nLevels, int nParts, ::coder::array<Omp4mPart, 1U> &parts)
  {
    ::coder::array<int, 1U> cparts_;
    ::coder::array<int, 1U> eind_;
    ::coder::array<int, 1U> eindloc_;
    ::coder::array<int, 1U> eptr_;
    ::coder::array<int, 1U> eptrloc_;
    ::coder::array<int, 1U> iwork_;
    ::coder::array<int, 1U> l2gmap_;
    ::coder::array<int, 1U> r;
    ::coder::array<int, 1U> r1;
    ::coder::array<int, 1U> vparts_;
    ::coder::array<boolean_T, 1U> visited_;
    Omp4mPart b_parts;
    int nnodes;

    if ((n < 1) || ((cells.size(0) == 0) || (cells.size(1) == 0)) || (dim < 1) ||
        (dim > 3)) {
      parts.set_size(0);
    } else {
      //  check inputs
      if (nLevels < 1) {
        nLevels = dim;
      }

      if (nParts < 1) {
        nParts = 1;

#ifdef _OPENMP

        nParts = omp_get_max_threads();

#endif //_OPENMP

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
        for (int i{0}; i < n; i++) {
          parts[0].part_list[i] = i + 1;
        }
      } else {
        int b_i;
        int idx;
        int m;
        int npts;
        int u0;
        int u1;

        parts.set_size(nLevels);

        //  buffer
        u0 = cells.size(0);
        u1 = n;
        if (u0 >= n) {
          u1 = u0;
        }

        iwork_.set_size(0);

        //  dual
        eptr_.set_size(n + 1);
        for (b_i = 0; b_i <= n; b_i++) {
          eptr_[b_i] = 0;
        }

        eptr_[0] = 1;

        //  determine number of incident cells
        m = cells.size(0) - 1;

        //  number of cells
        for (int cid{0}; cid <= m; cid++) {
          //  local function to get number of points per cell
          for (npts = cells.size(1) - 1; cells[npts + cells.size(1) * cid] <= 0;
               npts--) {
          }

          for (int i{0}; i <= npts; i++) {
            idx = cells[i + cells.size(1) * cid];
            eptr_[idx] = eptr_[idx] + 1;
          }
        }

        for (int i{0}; i < n; i++) {
          eptr_[i + 1] = eptr_[i + 1] + eptr_[i];
        }

        //  allocate n2cList
        eind_.set_size(eptr_[n] - 1);
        for (int cid{0}; cid <= m; cid++) {
          //  local function to get number of points per cell
          for (npts = cells.size(1) - 1; cells[npts + cells.size(1) * cid] <= 0;
               npts--) {
          }

          for (int i{0}; i <= npts; i++) {
            idx = cells[i + cells.size(1) * cid] - 1;
            eind_[eptr_[idx] - 1] = cid + 1;
            eptr_[idx] = eptr_[idx] + 1;
          }
        }

        for (int i{n}; i >= 1; i--) {
          eptr_[i] = eptr_[i - 1];
        }

        eptr_[0] = 1;

        //  call metis
        call_metis_mesh(cells.size(0), eptr_, eind_, nParts, vparts_, cparts_);

        //  build partition, note if nodal-based, then nparts is for cell
        visited_.set_size(u1);
        for (b_i = 0; b_i < u1; b_i++) {
          visited_[b_i] = false;
        }

        build_part(nParts, vparts_, cparts_, eptr_, eind_, visited_, iwork_,
                   parts[0].part_ptr, parts[0].part_list, parts[0].shared_ents);
        parts[0].nparts = nParts;
        if ((nLevels == 1) || (parts[0].shared_ents.size(0) == 0)) {
          if (nLevels > 1) {
            b_parts = parts[0];
            parts.set_size(1);
            parts[0] = b_parts;
          }
        } else {
          int lvl;
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
            extract_sub(cells.size(0), l2gmap_, eptr_, eind_, iwork_, visited_,
                        eptrloc_, eindloc_, &nnodes);

            //  call metis
            call_metis_mesh(nnodes, eptrloc_, eindloc_, nParts, vparts_, cparts_);

            //  build partition, note if nodal-based, then nparts is for cell
            build_part(nParts, vparts_, cparts_, eptrloc_, eindloc_, visited_,
                       iwork_, parts[lvl].part_ptr, r, r1);
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
            for (int i{0}; i < b_i; i++) {
              parts[lvl].part_list[i] = l2gmap_[parts[lvl].part_list[i] - 1];
            }

            b_i = parts[lvl].shared_ents.size(0);
            for (int i{0}; i < b_i; i++) {
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
        }
      }
    }
  }

  static void rdi_compute_oscind(const ::coder::array<double, 2U> &dfGlobal,
    const ::coder::array<double, 2U> &alphaCell, const ::coder::array<double, 2U>
    &mesh_xs, double mesh_hGlobal, const ::coder::array<double, 1U>
    &mesh_cellWeights, const ::coder::array<int, 1U> &mesh_n2cPtr, const ::coder::
    array<int, 1U> &mesh_n2cList, int params_dim, double params_epsBeta, ::coder::
    array<double, 2U> &beta)
  {
    double epsBeta;
    double hGlobal;

    // rdi_compute_oscind - Compute oscillation indicators (beta values)
    epsBeta = params_epsBeta;
    if (params_epsBeta <= 0.0) {
      epsBeta = 0.0002;
    }

    hGlobal = mesh_hGlobal;
    if (mesh_hGlobal <= 0.0) {
      double ex;
      int last;
      last = mesh_cellWeights.size(0);
      if (mesh_cellWeights.size(0) <= 2) {
        if (mesh_cellWeights.size(0) == 1) {
          ex = mesh_cellWeights[0];
        } else if (mesh_cellWeights[0] < mesh_cellWeights[mesh_cellWeights.size
                   (0) - 1]) {
          ex = mesh_cellWeights[mesh_cellWeights.size(0) - 1];
        } else {
          ex = mesh_cellWeights[0];
        }
      } else {
        ex = mesh_cellWeights[0];
        for (int k{2}; k <= last; k++) {
          double d;
          d = mesh_cellWeights[k - 1];
          if (ex < d) {
            ex = d;
          }
        }
      }

      hGlobal = std::pow(ex, 1.0 / static_cast<double>(params_dim));
    }

    compute_beta_kernel(mesh_xs.size(0), epsBeta, hGlobal, dfGlobal, alphaCell,
                        mesh_cellWeights, mesh_n2cPtr, mesh_n2cList,
                        alphaCell.size(1), beta);
  }

  static inline
  void rdi_compute_osusind(const ::coder::array<int, 1U> &rowPtr, const ::
    coder::array<int, 1U> &colInd, const ::coder::array<double, 1U> &vals, const
    ::coder::array<double, 2U> &fs, const ::coder::array<int, 1U> &mesh_n2cPtr,
    const ::coder::array<int, 1U> &mesh_n2cList, ::coder::array<double, 2U>
    &alphaCell, ::coder::array<double, 2U> &alphaNode)
  {
    // rdi_compute_osusind - Compute over-/under-shoot indicators
    crs_prod_mat_vec(rowPtr, colInd, vals, fs, alphaCell);
    compute_nodal_alpha(fs.size(0), alphaCell, mesh_n2cPtr, mesh_n2cList,
                        fs.size(1), alphaNode);
  }

  static void rrqr_factor(const ::coder::array<double, 2U> &A, double thres, int
    rowoffset, int coloffset, int m, int n, ::coder::array<double, 2U> &QR, ::
    coder::array<int, 1U> &p, int *rank, ::coder::array<double, 1U> &work)
  {
    int i;
    int wsize;

    //  rrqr_factor  Compute rank-revealing QR with column pivoting
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

    m2cAssert(p.size(0) >= n,
              "Length of permutation vector must be no smaller than the number of columns.");

    //  Allocate work space if needed
    wsize = wls::query_work_size(m, n);
    work.set_size(wsize);

    //  Invoke C++ function
    p[0] = 0;

    //  Note: A and Q are always stored in column major
    i = coloffset * A.size(1) + rowoffset;
    *rank = wls::rrqr_factor_nodag(&A[i % A.size(0) * A.size(1) + i / A.size(0)],
      thres, m, n, &QR[0], &(p.data())[0], &(work.data())[0], wsize, A.size(1));
  }

  static void rrqr_qmulti(const ::coder::array<double, 2U> &QR, int m, int n,
    int rank, ::coder::array<double, 2U> &bs, int nrhs, ::coder::array<double,
    1U> &work)
  {
    int stride_bs;
    int u1;
    int wsize;

    //  Perform Q*bs, where Q is stored implicitly in QR
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
      m2cErrMsgIdAndTxt("wlslib:WrongRank",
                        "Rank %d must be a positive value no greater than min(%d, %d).",
                        rank, m, n);
    }

    if (nrhs == 0) {
      nrhs = bs.size(0);
    }

    //  Resize work space if needed
    wsize = wls::query_work_size(m, n);
    work.set_size(wsize);

    //  zero out extra rows in bs to avoid errors in LAPACK
    u1 = n + 1;
    for (int i{u1}; i <= m; i++) {
      for (int j{0}; j < nrhs; j++) {
        bs[(i + bs.size(1) * j) - 1] = 0.0;
      }
    }

    //  Invoke C++ function
    wls::rrqr_qmulti(&QR[0], m, n, rank, QR.size(1), nrhs, &bs[0], stride_bs,
                     &(work.data())[0], wsize);
  }

  static void rrqr_rtsolve(const ::coder::array<double, 2U> &QR, int n, int rank,
    ::coder::array<double, 2U> &bs, int nrhs)
  {
    int i;

    //  Perform forward substitution to compute bs=R'\bs, where R is stored in QR
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
      m2cErrMsgIdAndTxt("wlslib:WrongRank",
                        "Rank %d must be a positive value no greater than min(%d, %d).",
                        rank, QR.size(1), n);
    }

    if (nrhs == 0) {
      nrhs = bs.size(0);
    }

    //  Obtain stride
    wls::rrqr_rtsolve(&QR[0], n, rank, QR.size(1), nrhs, &bs[0], bs.size(1));
  }

  static void update_osusop(int n, const ::coder::array<int, 2U> &stcls, const ::
    coder::array<int, 1U> &nrange, const ::coder::array<int, 1U> &n2cPtr, const ::
    coder::array<int, 1U> &n2cList, ::coder::array<int, 1U> &rowPtr, ::coder::
    array<int, 1U> &colInd, ::coder::array<double, 1U> &vals, ::coder::array<int,
    1U> &nnzPr)
  {
    ::coder::array<double, 1U> b_vals;
    ::coder::array<int, 1U> b_colInd;
    ::coder::array<boolean_T, 1U> visited_;
    int N;
    int i;

    //  function for updating uncompressed CRS
    N = nrange.size(0);
    visited_.set_size(n);
    for (i = 0; i < n; i++) {
      visited_[i] = false;
    }

    for (int b_i{0}; b_i < N; b_i++) {
      int i1;
      i = n2cPtr[nrange[b_i] - 1];
      i1 = n2cPtr[nrange[b_i]] - 1;
      for (int j{i}; j <= i1; j++) {
        int i2;
        int i3;
        int i4;
        int n0;
        i2 = n2cList[j - 1];
        i3 = rowPtr[i2 - 1];
        i4 = (i3 + nnzPr[i2 - 1]) - 1;
        for (int k{i3}; k <= i4; k++) {
          visited_[colInd[k - 1] - 1] = true;
        }

        n0 = stcls[(stcls.size(1) + stcls.size(1) * b_i) - 1];
        for (int k{0}; k < n0; k++) {
          if (!visited_[stcls[k + stcls.size(1) * b_i] - 1]) {
            i3 = rowPtr[i2] - rowPtr[i2 - 1];
            if (nnzPr[i2 - 1] >= i3) {
              int b_n;
              int dif;
              int u0;
              int u1;
              u0 = static_cast<int>(std::round(static_cast<double>(i3) * 1.2));
              u1 = i3 + 1;
              if (u0 >= u1) {
                u1 = u0;
              }

              dif = u1 - i3;
              if (rowPtr[i2] - 1 < 1) {
                u1 = 0;
              } else {
                u1 = rowPtr[i2] - 1;
              }

              if (rowPtr[i2] > colInd.size(0)) {
                i3 = 0;
                i4 = 0;
              } else {
                i3 = rowPtr[i2] - 1;
                i4 = colInd.size(0);
              }

              b_colInd.set_size(((u1 + dif) + i4) - i3);
              for (u0 = 0; u0 < u1; u0++) {
                b_colInd[u0] = colInd[u0];
              }

              for (u0 = 0; u0 < dif; u0++) {
                b_colInd[u0 + u1] = 0;
              }

              u0 = i4 - i3;
              for (i4 = 0; i4 < u0; i4++) {
                b_colInd[(i4 + u1) + dif] = colInd[i3 + i4];
              }

              colInd.set_size(b_colInd.size(0));
              u1 = b_colInd.size(0);
              for (i3 = 0; i3 < u1; i3++) {
                colInd[i3] = b_colInd[i3];
              }

              if (rowPtr[i2] - 1 < 1) {
                u1 = 0;
              } else {
                u1 = rowPtr[i2] - 1;
              }

              if (rowPtr[i2] > vals.size(0)) {
                i3 = 0;
                i4 = 0;
              } else {
                i3 = rowPtr[i2] - 1;
                i4 = vals.size(0);
              }

              b_vals.set_size(((u1 + dif) + i4) - i3);
              for (u0 = 0; u0 < u1; u0++) {
                b_vals[u0] = vals[u0];
              }

              for (u0 = 0; u0 < dif; u0++) {
                b_vals[u0 + u1] = 0.0;
              }

              u0 = i4 - i3;
              for (i4 = 0; i4 < u0; i4++) {
                b_vals[(i4 + u1) + dif] = vals[i3 + i4];
              }

              vals.set_size(b_vals.size(0));
              u1 = b_vals.size(0);
              for (i3 = 0; i3 < u1; i3++) {
                vals[i3] = b_vals[i3];
              }

              b_n = rowPtr.size(0);
              i3 = i2 + 1;
              for (int b_j{i3}; b_j <= b_n; b_j++) {
                rowPtr[b_j - 1] = rowPtr[b_j - 1] + dif;
              }
            }

            i3 = stcls[k + stcls.size(1) * b_i];
            colInd[(rowPtr[i2 - 1] + nnzPr[i2 - 1]) - 1] = i3;
            nnzPr[i2 - 1] = nnzPr[i2 - 1] + 1;
            visited_[i3 - 1] = true;
          }
        }

        i3 = rowPtr[i2 - 1];
        i2 = (rowPtr[i2 - 1] + nnzPr[i2 - 1]) - 1;
        for (int k{i3}; k <= i2; k++) {
          visited_[colInd[k - 1] - 1] = false;
        }
      }
    }
  }

  static void wls_buhmann_weights(const ::coder::array<double, 2U> &us, int
    npoints, int degree, const double params_pw_data[], const int
    params_pw_size[2], ::coder::array<double, 1U> &ws)
  {
    static const double b_dv[7]{ 2.6, 2.0, 1.6, 1.6, 1.6, 1.5, 1.4 };

    double d;
    double dist_k;
    double r;
    double r1;
    double r2;
    double rho;
    double sigma;
    int abs_degree;
    int i;

    //  Weights based on Buhmann's radial basis function
    if (degree == 0) {
      degree = 2;
    }

    if (degree < 0) {
      abs_degree = 1 - degree;
    } else {
      abs_degree = degree + 1;
    }

    //  Assign default rho
    if (abs_degree - 1 >= 7) {
      sigma = 1.4;
    } else {
      sigma = b_dv[abs_degree - 2];
    }

    ws.set_size(npoints);

    //  Compute rho to be sigma times the kth distance for k=ceil(1.5*ncoff)
    if (degree >= 0) {
      //  Compute 2-norm
      i = us.size(1);
      for (int b_i{0}; b_i < npoints; b_i++) {
        d = us[us.size(1) * b_i];
        r2 = d * d;
        for (int j{2}; j <= i; j++) {
          d = us[(j + us.size(1) * b_i) - 1];
          r2 += d * d;
        }

        ws[b_i] = std::sqrt(r2);
      }
    } else {
      //  Compute inf-norm for tensor-product
      i = us.size(1);
      for (int b_i{0}; b_i < npoints; b_i++) {
        r = std::abs(us[us.size(1) * b_i]);
        for (int j{2}; j <= i; j++) {
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
    if ((params_pw_size[0] == 0) || (params_pw_size[1] == 0)) {
      for (int b_i{0}; b_i < npoints; b_i++) {
        if (degree > 0) {
          //  Compute 2-norm
          d = us[us.size(1) * b_i];
          r2 = d * d;
          i = us.size(1);
          for (int j{2}; j <= i; j++) {
            d = us[(j + us.size(1) * b_i) - 1];
            r2 += d * d;
          }

          r = std::sqrt(r2);
        } else {
          //  Compute inf-norm for tensor-product
          r = std::abs(us[us.size(1) * b_i]);
          i = us.size(1);
          for (int j{2}; j <= i; j++) {
            r1 = std::abs(us[(j + us.size(1) * b_i) - 1]);
            if (r1 > r) {
              r = r1;
            }
          }
        }

        if (r > rho) {
          ws[b_i] = 0.0;
        } else {
          double r_sqrt;
          r /= rho;
          r_sqrt = std::sqrt(r);
          ws[b_i] = r * r * (r * r_sqrt * (r_sqrt * (r_sqrt * 112.0 / 45.0 +
            -7.0) + 5.333333333333333) + -0.93333333333333335) +
            0.1111111111111111;
        }
      }
    } else {
      m2cAssert(params_pw_size[0] >= npoints,
                "size(params_pw,1) should be >=npoints");
      for (int b_i{0}; b_i < npoints; b_i++) {
        double b_gamma;
        b_gamma = params_pw_data[params_pw_size[1] * b_i];
        if (b_gamma <= 0.0) {
          ws[b_i] = 0.0;
        } else {
          if (degree > 0) {
            //  Compute 2-norm
            d = us[us.size(1) * b_i];
            r2 = d * d;
            i = us.size(1);
            for (int j{2}; j <= i; j++) {
              d = us[(j + us.size(1) * b_i) - 1];
              r2 += d * d;
            }

            r = std::sqrt(r2);
          } else {
            //  Compute inf-norm for tensor-product
            r = std::abs(us[us.size(1) * b_i]);
            i = us.size(1);
            for (int j{2}; j <= i; j++) {
              r1 = std::abs(us[(j + us.size(1) * b_i) - 1]);
              if (r1 > r) {
                r = r1;
              }
            }
          }

          if (r > rho) {
            ws[b_i] = 0.0;
          } else {
            double r_sqrt;
            r /= rho;
            r_sqrt = std::sqrt(r);
            ws[b_i] = b_gamma * (r * r * (r * r_sqrt * (r_sqrt * (r_sqrt * 112.0
              / 45.0 + -7.0) + 5.333333333333333) + -0.93333333333333335) +
                                 0.1111111111111111);
          }
        }
      }
    }
  }

  static void wls_func(WlsObject *wls, const double pnts_data[], const int
                       pnts_size[2], int npoints, ::coder::array<double, 2U>
                       &vdops)
  {
    int j;
    int nDims;
    int nrows;
    int nrows_vdops;
    int u0;
    int u1;

    //  wls_func  Compute wls--fitting at one or more points.
    nDims = pnts_size[1] - 1;

    //  scale the coordinates; use wls.us as buffer
    wls->us.set_size(((npoints + 3) / 4) << 2, pnts_size[1]);
    if (wls->interp0 != 0) {
      //  coordinate system is centered at first node in interp0 mode
      for (int dim{0}; dim <= nDims; dim++) {
        for (int iPoint{0}; iPoint < npoints; iPoint++) {
          wls->us[dim + wls->us.size(1) * iPoint] = (pnts_data[dim + pnts_size[1]
            * iPoint] - wls->origin.data[dim]) * wls->hs_inv.data[dim];
        }
      }
    } else {
      for (int dim{0}; dim <= nDims; dim++) {
        for (int iPoint{0}; iPoint < npoints; iPoint++) {
          wls->us[dim + wls->us.size(1) * iPoint] = pnts_data[dim + pnts_size[1]
            * iPoint] * wls->hs_inv.data[dim];
        }
      }
    }

    //  compute the generalized Vandermonde matrix and right-hand side
    gen_vander(wls->us, npoints, wls->degree, wls->V);
    u0 = wls->ncols;
    u1 = wls->nrows;
    if (u0 >= u1) {
      nrows_vdops = u0;
    } else {
      nrows_vdops = u1;
    }

    //  force each operator (rhs) to be stored contiguously
    u1 = nrows_vdops - wls->interp0;
    wls->vdops.set_size(npoints, u1);

    //  Extract vopts from Vandermonde matrix
    u0 = wls->ncols - wls->interp0;
    for (int iMonomial{0}; iMonomial < u0; iMonomial++) {
      j = wls->jpvt[iMonomial] + wls->interp0;
      for (int iPoint{0}; iPoint < npoints; iPoint++) {
        wls->vdops[iMonomial + wls->vdops.size(1) * iPoint] = wls->V[iPoint
          + wls->V.size(1) * (j - 1)];
      }
    }

    //  Multiply by generalized inverse of Vandermonde matrix
    rrqr_rtsolve(wls->QR, wls->ncols - wls->interp0, wls->rank, wls->vdops,
                 npoints);
    rrqr_qmulti(wls->QR, wls->nrows - wls->interp0, wls->ncols - wls->interp0,
                wls->rank, wls->vdops, npoints, wls->work);
    vdops.set_size(nrows_vdops, npoints);

    //  Transpose the operator for row-major
    for (int i{0}; i < u1; i++) {
      for (j = 0; j < npoints; j++) {
        vdops[j + vdops.size(1) * (i + wls->interp0)] = wls->vdops[i +
          wls->vdops.size(1) * j];
      }
    }

    nrows = wls->nrows;
    if (wls->rweights.size(0) != 0) {
      for (int k{0}; k < npoints; k++) {
        for (int iRow{0}; iRow < nrows; iRow++) {
          vdops[k + vdops.size(1) * iRow] = vdops[k + vdops.size(1) * iRow] *
            wls->rweights[iRow];
        }
      }
    }

    if (wls->interp0 != 0) {
      //  In interp0 mode, we set the first entry based on partition of unity
      for (j = 0; j < npoints; j++) {
        double s;
        s = 0.0;
        u1 = wls->npoints;
        for (int i{2}; i <= u1; i++) {
          s += vdops[j + vdops.size(1) * (i - 1)];
        }

        vdops[j] = 1.0 - s;
      }
    }
  }

  static void wls_func(WlsObject *wls, const double pnts_data[], int npoints, ::
                       coder::array<double, 2U> &vdops)
  {
    int j;
    int nrows;
    int nrows_vdops;
    int u0;
    int u1;

    //  wls_func  Compute wls--fitting at one or more points.
    wls->us.set_size(((npoints + 3) / 4) << 2, 1);
    if (wls->interp0 != 0) {
      //  coordinate system is centered at first node in interp0 mode
      for (int iPoint{0}; iPoint < npoints; iPoint++) {
        wls->us[wls->us.size(1) * iPoint] = (pnts_data[iPoint] -
          wls->origin.data[0]) * wls->hs_inv.data[0];
      }
    } else {
      for (int iPoint{0}; iPoint < npoints; iPoint++) {
        wls->us[wls->us.size(1) * iPoint] = pnts_data[iPoint] * wls->
          hs_inv.data[0];
      }
    }

    //  compute the generalized Vandermonde matrix and right-hand side
    gen_vander(wls->us, npoints, wls->degree, wls->V);
    u0 = wls->ncols;
    u1 = wls->nrows;
    if (u0 >= u1) {
      nrows_vdops = u0;
    } else {
      nrows_vdops = u1;
    }

    //  force each operator (rhs) to be stored contiguously
    u1 = nrows_vdops - wls->interp0;
    wls->vdops.set_size(npoints, u1);

    //  Extract vopts from Vandermonde matrix
    u0 = wls->ncols - wls->interp0;
    for (int iMonomial{0}; iMonomial < u0; iMonomial++) {
      j = wls->jpvt[iMonomial] + wls->interp0;
      for (int iPoint{0}; iPoint < npoints; iPoint++) {
        wls->vdops[iMonomial + wls->vdops.size(1) * iPoint] = wls->V[iPoint
          + wls->V.size(1) * (j - 1)];
      }
    }

    //  Multiply by generalized inverse of Vandermonde matrix
    rrqr_rtsolve(wls->QR, wls->ncols - wls->interp0, wls->rank, wls->vdops,
                 npoints);
    rrqr_qmulti(wls->QR, wls->nrows - wls->interp0, wls->ncols - wls->interp0,
                wls->rank, wls->vdops, npoints, wls->work);
    vdops.set_size(nrows_vdops, npoints);

    //  Transpose the operator for row-major
    for (int i{0}; i < u1; i++) {
      for (j = 0; j < npoints; j++) {
        vdops[j + vdops.size(1) * (i + wls->interp0)] = wls->vdops[i +
          wls->vdops.size(1) * j];
      }
    }

    nrows = wls->nrows;
    if (wls->rweights.size(0) != 0) {
      for (int k{0}; k < npoints; k++) {
        for (int iRow{0}; iRow < nrows; iRow++) {
          vdops[k + vdops.size(1) * iRow] = vdops[k + vdops.size(1) * iRow] *
            wls->rweights[iRow];
        }
      }
    }

    if (wls->interp0 != 0) {
      //  In interp0 mode, we set the first entry based on partition of unity
      for (j = 0; j < npoints; j++) {
        double s;
        s = 0.0;
        u1 = wls->npoints;
        for (int i{2}; i <= u1; i++) {
          s += vdops[j + vdops.size(1) * (i - 1)];
        }

        vdops[j] = 1.0 - s;
      }
    }
  }

  static void wls_init(WlsObject *wls, const double us_data[], const int
                       us_size[2], const char weight_name_data[], const double
                       weight_params_pointwise_data[], const int
                       weight_params_pointwise_size[2], int degree, boolean_T
                       interp0, boolean_T use_dag, int npoints)
  {
    double maxx;
    double maxx_inv;
    double thres;
    int b_interp0;
    int dim;
    int i;
    int r_tmp;

    //  wls_init  Initialize WlsObject in 1D, 2D, or 3D.
    m2cAssert(true, "");

    //  Process input arguments
    dim = us_size[1];
    wls->interp0 = interp0;
    b_interp0 = wls->interp0;
    wls->use_dag = use_dag;
    if (npoints <= 0) {
      npoints = us_size[0];
    } else {
      m2cAssert(npoints <= us_size[0],
                "Number of points cannot be greater than the first dimension of `us`.");
    }

    //  Resize buffers
    wls_resize(wls, us_size[1], npoints, degree, use_dag);

    //  Recompute DAG if use_dag and its signature does not match
    if (use_dag) {
      i = wls->dag.size(1) * wls->dag.size(0) - 1;
      r_tmp = wls->dag.size(0);
      if (wls->dag[i % r_tmp * wls->dag.size(1) + i / r_tmp] != degree + 127) {
        //  Wrapper function for building DAG in nD.
        if (us_size[1] == 2) {
          gen_vander_2d_dag(degree, wls->dag);
        } else {
          gen_vander_3d_dag(degree, wls->dag);
        }
      }
    }

    if (wls->interp0 != 0) {
      //  Make the first node the origin in interp0 mode
      if (us_size[1] == 2) {
        wls->origin.size[1] = 2;
        wls->origin.size[0] = 1;
        wls->origin.data[0] = us_data[0];
        wls->origin.data[1] = us_data[1];
        for (int b_i{0}; b_i < npoints; b_i++) {
          wls->us[wls->us.size(1) * b_i] = us_data[2 * b_i] - us_data[0];
          wls->us[wls->us.size(1) * b_i + 1] = us_data[2 * b_i + 1] - us_data[1];
        }
      } else {
        wls->origin.size[1] = 3;
        wls->origin.size[0] = 1;
        wls->origin.data[0] = us_data[0];
        wls->origin.data[1] = us_data[1];
        wls->origin.data[2] = us_data[2];
        for (int b_i{0}; b_i < npoints; b_i++) {
          i = us_size[1] * b_i;
          wls->us[wls->us.size(1) * b_i] = us_data[i] - us_data[0];
          wls->us[wls->us.size(1) * b_i + 1] = us_data[i + 1] - us_data[1];
          wls->us[wls->us.size(1) * b_i + 2] = us_data[i + 2] - us_data[2];
        }
      }
    } else if (us_size[1] == 2) {
      wls->origin.size[1] = 2;
      wls->origin.size[0] = 1;
      wls->origin.data[0] = 0.0;
      wls->origin.data[1] = 0.0;
      for (int b_i{0}; b_i < npoints; b_i++) {
        wls->us[wls->us.size(1) * b_i] = us_data[2 * b_i];
        wls->us[wls->us.size(1) * b_i + 1] = us_data[2 * b_i + 1];
      }
    } else {
      wls->origin.size[1] = 3;
      wls->origin.size[0] = 1;
      wls->origin.data[0] = 0.0;
      wls->origin.data[1] = 0.0;
      wls->origin.data[2] = 0.0;
      for (int b_i{0}; b_i < npoints; b_i++) {
        i = us_size[1] * b_i;
        wls->us[wls->us.size(1) * b_i] = us_data[i];
        wls->us[wls->us.size(1) * b_i + 1] = us_data[i + 1];
        wls->us[wls->us.size(1) * b_i + 2] = us_data[i + 2];
      }
    }

    //  Scale us to be between -1 and 1
    maxx = 0.0;
    if (us_size[1] == 2) {
      for (int b_i{0}; b_i < npoints; b_i++) {
        maxx = std::fmax(maxx, std::fmax(std::abs(wls->us[wls->us.size(1) * b_i]),
          std::abs(wls->us[wls->us.size(1) * b_i + 1])));
      }
    } else {
      for (int b_i{0}; b_i < npoints; b_i++) {
        maxx = std::fmax(maxx, std::fmax(std::fmax(std::abs(wls->us[wls->us.size
          (1) * b_i]), std::abs(wls->us[wls->us.size(1) * b_i + 1])), std::abs
          (wls->us[wls->us.size(1) * b_i + 2])));
      }
    }

    if (maxx == 0.0) {
      maxx_inv = 1.0;
    } else {
      maxx_inv = 1.0 / maxx;
    }

    for (int b_i{0}; b_i < dim; b_i++) {
      wls->hs_inv.data[b_i] = maxx_inv;
    }

    //  scale wls.us
    if (maxx_inv != 1.0) {
      if (us_size[1] == 2) {
        for (int b_i{0}; b_i < npoints; b_i++) {
          wls->us[wls->us.size(1) * b_i] = wls->us[wls->us.size(1) * b_i] *
            maxx_inv;
          wls->us[wls->us.size(1) * b_i + 1] = wls->us[wls->us.size(1) * b_i + 1]
            * maxx_inv;
        }
      } else {
        for (int b_i{0}; b_i < npoints; b_i++) {
          wls->us[wls->us.size(1) * b_i] = wls->us[wls->us.size(1) * b_i] *
            maxx_inv;
          wls->us[wls->us.size(1) * b_i + 1] = wls->us[wls->us.size(1) * b_i + 1]
            * maxx_inv;
          wls->us[wls->us.size(1) * b_i + 2] = wls->us[wls->us.size(1) * b_i + 2]
            * maxx_inv;
        }
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
        r_tmp = wls->rweights.size(0);
        wls->rweights.set_size(r_tmp);
        for (i = 0; i < r_tmp; i++) {
          wls->rweights[i] = 1.0;
        }
      } else if ((weight_name_data[0] == 'I') || (weight_name_data[0] == 'i')) {
        //  inverse distance
        wls_invdist_weights(wls->us, npoints, degree,
                            weight_params_pointwise_data,
                            weight_params_pointwise_size, wls->rweights);
      } else if ((weight_name_data[0] == 'B') || (weight_name_data[0] == 'b')) {
        //  Buhmann weights. All points share same parameters
        wls_buhmann_weights(wls->us, npoints, degree,
                            weight_params_pointwise_data,
                            weight_params_pointwise_size, wls->rweights);
      } else {
        if ((weight_name_data[0] != 'E') && (weight_name_data[0] != 'e')) {
          m2cErrMsgIdAndTxt("wlslib:WrongWeightName",
                            "Weighting scheme must be unit, inverse, Buhmann, or WLS-ENO.");
        }

        //  WLS-ENO
        m2cAssert(false, "first two shared parameters are required");

        m2cAssert(weight_params_pointwise_size[0] >= npoints,
                  "size(params_pw,1) should be >=npoints");

        m2cAssert(false, "size(params_pw,2) should be >=2");
        wls->rweights.set_size(npoints);

        //  Compute hbar using ws as buffer space
        if (degree >= 0) {
          //  Compute 2-norm
          i = us_size[1];
          for (int b_i{0}; b_i < npoints; b_i++) {
            double d;
            double r2;
            r_tmp = us_size[1] * b_i;
            d = us_data[r_tmp];
            r2 = d * d;
            for (int j{2}; j <= i; j++) {
              d = us_data[(j + r_tmp) - 1];
              r2 += d * d;
            }

            wls->rweights[b_i] = std::sqrt(r2);
          }
        } else {
          //  Compute inf-norm for tensor-product
          i = us_size[1];
          for (int b_i{0}; b_i < npoints; b_i++) {
            double r;
            r_tmp = us_size[1] * b_i;
            r = std::abs(us_data[r_tmp]);
            for (int j{2}; j <= i; j++) {
              double r1;
              r1 = std::abs(us_data[(j + r_tmp) - 1]);
              if (r1 > r) {
                r = r1;
              }
            }

            wls->rweights[b_i] = r;
          }
        }

        //  Evaluate the inverse-distance weights as base
        wls_invdist_weights(wls->us, npoints, 0.5 - static_cast<double>(degree <
          0), wls->rweights);

        // A check that is always false is detected at compile-time. Eliminating code that follows.
      }
    }

    //  Compute Vandermonde system and recompute DAG if needed
    gen_vander(wls->us, npoints, degree, wls->rweights, wls->V);
    wls->ncols = wls->V.size(0);

    //  Compact CVM if needed
    wls->nrows = wls->V.size(1) / wls->stride * npoints;

    //  Omit rows in CVM if needed
    if ((degree > 1) && (degree < 7)) {
      thres = dv[degree - 1];
    } else {
      thres = 1.0E+8;
    }

    //  In interp0 mode, we trim off the first row and first column.
    rrqr_factor(wls->V, thres, b_interp0, b_interp0, wls->nrows - b_interp0,
                wls->ncols - b_interp0, wls->QR, wls->jpvt, &wls->rank,
                wls->work);
  }

  static void wls_invdist_weights(const ::coder::array<double, 2U> &us, int
    npoints, int degree, const double params_pw_data[], const int
    params_pw_size[2], ::coder::array<double, 1U> &ws)
  {
    double alpha;
    int b_degree;

    //  Weights based on inverse distance
    if (degree < 0) {
      b_degree = -degree;
    } else {
      b_degree = degree;
    }

    alpha = static_cast<double>(b_degree) / 2.0;
    ws.set_size(npoints);
    if ((params_pw_size[0] == 0) || (params_pw_size[1] == 0)) {
      for (int i{0}; i < npoints; i++) {
        double r;
        double r2;
        r = std::abs(us[us.size(1) * i]);
        if (us.size(1) > 1) {
          if (degree > 0) {
            //  Compute 2-norm
            r2 = r * r;
            b_degree = us.size(1);
            for (int b_i{2}; b_i <= b_degree; b_i++) {
              double d;
              d = us[(b_i + us.size(1) * i) - 1];
              r2 += d * d;
            }
          } else {
            //  Compute inf-norm for tensor-product
            b_degree = us.size(1);
            for (int b_i{2}; b_i <= b_degree; b_i++) {
              double r1;
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
        ws[i] = std::pow(r2 + 0.01, -alpha);
      }
    } else {
      m2cAssert(params_pw_size[0] >= npoints,
                "size(params_pw,1) should be >=npoints");
      for (int i{0}; i < npoints; i++) {
        double b_gamma;
        b_gamma = params_pw_data[params_pw_size[1] * i];
        if (b_gamma <= 0.0) {
          ws[i] = 0.0;
        } else {
          double r;
          double r2;
          r = std::abs(us[us.size(1) * i]);
          if (us.size(1) > 1) {
            if (degree > 0) {
              //  Compute 2-norm
              r2 = r * r;
              b_degree = us.size(1);
              for (int b_i{2}; b_i <= b_degree; b_i++) {
                double d;
                d = us[(b_i + us.size(1) * i) - 1];
                r2 += d * d;
              }
            } else {
              //  Compute inf-norm for tensor-product
              b_degree = us.size(1);
              for (int b_i{2}; b_i <= b_degree; b_i++) {
                double r1;
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
          ws[i] = b_gamma * std::pow(r2 + 0.01, -alpha);
        }
      }
    }
  }

  static void wls_invdist_weights(const ::coder::array<double, 2U> &us, int
    npoints, double degree, ::coder::array<double, 1U> &ws)
  {
    double alpha;

    //  Weights based on inverse distance
    alpha = std::abs(degree) / 2.0;
    ws.set_size(npoints);
    for (int i{0}; i < npoints; i++) {
      double r;
      double r2;
      r = std::abs(us[us.size(1) * i]);
      if (us.size(1) > 1) {
        if (degree > 0.0) {
          int b_i;

          //  Compute 2-norm
          r2 = r * r;
          b_i = us.size(1);
          for (int c_i{2}; c_i <= b_i; c_i++) {
            double d;
            d = us[(c_i + us.size(1) * i) - 1];
            r2 += d * d;
          }
        } else {
          int b_i;

          //  Compute inf-norm for tensor-product
          b_i = us.size(1);
          for (int c_i{2}; c_i <= b_i; c_i++) {
            double r1;
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
      ws[i] = std::pow(r2 + 0.01, -alpha);
    }
  }

  static void wls_resize(WlsObject *wls, int dim, int npoints, int degree,
    boolean_T use_dag)
  {
    int ncols;
    int stride;
    int stride_idx_0;
    int u0;

    //  wls_resize  Reinitialize the buffers of WlsObject
    wls->degree = degree;
    wls->order = 0;
    wls->use_dag = use_dag;

    //  make stride a multiple of four
    stride = ((npoints + 3) / 4) << 2;
    wls->stride = stride;
    wls->us.set_size(stride, dim);
    wls->rweights.set_size(stride);
    wls->npoints = npoints;

    //  determine number of columns and allocate V and QR
    if (dim == 2) {
      if (degree > 0) {
        ncols = (degree + 1) * (degree + 2) / 2;
      } else {
        ncols = (1 - degree) * (1 - degree);
      }
    } else {
      if (degree > 0) {
        ncols = (degree + 1) * (degree + 2) * (degree + 3) / 6;
      } else {
        ncols = (1 - degree) * (1 - degree) * (1 - degree);
      }
    }

    wls->hs_inv.size[1] = dim;
    wls->hs_inv.size[0] = 1;
    if (use_dag) {
      if ((dim != wls->dag.size(1)) || (ncols + 1 != wls->dag.size(0))) {
        //  Reset DAG if dimension or degree has changed.
        wls->dag.set_size(ncols + 1, dim);
        stride_idx_0 = wls->dag.size(1) * wls->dag.size(0) - 1;
        u0 = wls->dag.size(0);
        wls->dag[stride_idx_0 % u0 * wls->dag.size(1) + stride_idx_0 / u0] =
          MAX_uint8_T;
      }
    } else {
      wls->dag.set_size(0, dim);
    }

    wls->jpvt.set_size(ncols);

    //  V is always full, but QR has one fewer row and column in interp0 mode
    wls->V.set_size(ncols, stride);
    stride_idx_0 = (ncols - wls->interp0) + 1;
    wls->QR.set_size(stride_idx_0, stride);
    wls->rank = 0;

    //  work space
    u0 = ncols << 2;
    stride_idx_0 = ncols + 1;
    if (stride >= stride_idx_0) {
      stride_idx_0 = stride;
    }

    if (u0 < 4160) {
      u0 = 4160;
    }

    wls->work.set_size((stride_idx_0 << 5) + u0);
  }

  void librdi_initialize()
  {
  }

  void librdi_terminate()
  {
  }

  void rdi_assemble_osusop(const RdiMesh *mesh, const ::coder::array<int, 2U>
    &stcls, const RdiParams *params, ::coder::array<int, 1U> &rowPtr, ::coder::
    array<int, 1U> &colInd, ::coder::array<double, 1U> &vals, ::coder::array<int,
    1U> &nnzPr, ::coder::array<int, 1U> &rdNodes)
  {
    ::coder::array<double, 1U> b_vals;
    ::coder::array<int, 1U> b_colInd;
    ::coder::array<int, 1U> b_rowPtr;
    double tEnd;
    double tStart;
    int maxStclPr;
    boolean_T isSurf;

    // rdi_assemble_osusop - Assemble the WALF-based OSUS operator
    isSurf = (params->dim == mesh->xs.size(1) - 1);

    //  surface computation tag
    maxStclPr = params->maxStclSize;
    if (params->maxStclSize <= 0) {
      if (params->degree < 5) {
        maxStclPr = iv[(params->degree + ((params->dim - 1) << 2)) - 1];
      } else {
        maxStclPr = iv[((params->dim - 1) << 2) + 3];
      }
    }

    if (isSurf && (mesh->xs.size(1) < 2)) {
      m2cErrMsgIdAndTxt("rdi_assemble_rdiop:badDim",
                        "surface computation only supports 2D or 3D");
    }

    // Provides a portable wall clock timing routine.
    tStart = 0.0;

#ifdef _OPENMP

    tStart = omp_get_wtime();

#endif //_OPENMP

    if (rowPtr.size(0) == 0) {
      //  call on no rank-deficient nodes
      init_osusop(mesh->conn, stcls, maxStclPr, rowPtr, colInd, vals, nnzPr);
      rdNodes.set_size(0);
    } else {
      //  if we get here, it means we need to resolve rd nodes
      update_osusop(mesh->xs.size(0), stcls, rdNodes, mesh->n2cPtr,
                    mesh->n2cList, rowPtr, colInd, vals, nnzPr);
    }

    // Provides a portable wall clock timing routine.
    tEnd = 0.0;

#ifdef _OPENMP

    tEnd = omp_get_wtime();

#endif //_OPENMP

    if (params->verbose > 1) {
      m2cPrintf(" Init or updated OSUS operator in %gs...\n", tEnd - tStart);
      fflush(stdout);
    }

    // Provides a portable wall clock timing routine.
    tStart = 0.0;

#ifdef _OPENMP

    tStart = omp_get_wtime();

#endif //_OPENMP

    if (isSurf) {
      //  surface assembly
      assemble_surf(mesh->xs, mesh->conn, mesh->dirs, mesh->n2cPtr,
                    mesh->n2cList, mesh->parts, stcls, rowPtr, colInd, vals,
                    nnzPr, rdNodes, params);
    } else {
      //  body assembly
      assemble_body(mesh->xs, mesh->conn, mesh->n2cPtr, mesh->n2cList,
                    mesh->parts, stcls, rowPtr, colInd, vals, nnzPr, rdNodes,
                    params);
    }

    // Provides a portable wall clock timing routine.
    tEnd = 0.0;

#ifdef _OPENMP

    tEnd = omp_get_wtime();

#endif //_OPENMP

    if (params->verbose > 1) {
      m2cPrintf(" Assembly OSUS operator finished in %gs...\n", tEnd - tStart);
      fflush(stdout);
    }

    if (rdNodes.size(0) == 0) {
      int b_i;
      int loop_ub;
      int n;

      n = rowPtr.size(0) - 2;
      b_rowPtr.set_size(rowPtr.size(0));
      b_rowPtr[0] = 1;
      for (int i{0}; i <= n; i++) {
        b_rowPtr[i + 1] = b_rowPtr[i] + nnzPr[i];
      }

      b_colInd.set_size(b_rowPtr[rowPtr.size(0) - 1] - 1);
      b_vals.set_size(b_rowPtr[rowPtr.size(0) - 1] - 1);
      for (int i{0}; i <= n; i++) {
        int j;
        j = rowPtr[i] - 1;
        b_i = b_rowPtr[i];
        loop_ub = b_rowPtr[i + 1] - 1;
        for (int k{b_i}; k <= loop_ub; k++) {
          int colInd_tmp;
          colInd_tmp = (j + k) - b_i;
          b_colInd[k - 1] = colInd[colInd_tmp];
          b_vals[k - 1] = vals[colInd_tmp];
        }
      }

      rowPtr.set_size(b_rowPtr.size(0));
      loop_ub = b_rowPtr.size(0);
      for (b_i = 0; b_i < loop_ub; b_i++) {
        rowPtr[b_i] = b_rowPtr[b_i];
      }

      colInd.set_size(b_colInd.size(0));
      loop_ub = b_colInd.size(0);
      for (b_i = 0; b_i < loop_ub; b_i++) {
        colInd[b_i] = b_colInd[b_i];
      }

      vals.set_size(b_vals.size(0));
      loop_ub = b_vals.size(0);
      for (b_i = 0; b_i < loop_ub; b_i++) {
        vals[b_i] = b_vals[b_i];
      }
    }
  }

  void rdi_assemble_osusop2(const RdiMesh *mesh, const ::coder::array<int, 2U>
    &stcls, const RdiParams *params, ::coder::array<int, 1U> &rowPtr, ::coder::
    array<int, 1U> &colInd, ::coder::array<double, 1U> &vals, ::coder::array<int,
    1U> &nnzPr, ::coder::array<int, 1U> &rdNodes)
  {
    ::coder::array<double, 1U> b_vals;
    ::coder::array<int, 1U> b_colInd;
    ::coder::array<int, 1U> b_rowPtr;
    double tEnd;
    double tStart;
    int maxStclPr;
    boolean_T isSurf;

    // rdi_assemble_osusop - Assemble the WALF-based OSUS operator
    isSurf = (params->dim == mesh->xs.size(1) - 1);

    //  surface computation tag
    maxStclPr = params->maxStclSize;
    if (params->maxStclSize <= 0) {
      if (params->degree < 5) {
        maxStclPr = iv[(params->degree + ((params->dim - 1) << 2)) - 1];
      } else {
        maxStclPr = iv[((params->dim - 1) << 2) + 3];
      }
    }

    if (isSurf && (mesh->xs.size(1) < 2)) {
      m2cErrMsgIdAndTxt("rdi_assemble_rdiop:badDim",
                        "surface computation only supports 2D or 3D");
    }

    // Provides a portable wall clock timing routine.
    tStart = 0.0;

#ifdef _OPENMP

    tStart = omp_get_wtime();

#endif //_OPENMP

    //  call on no rank-deficient nodes
    init_osusop(mesh->conn, stcls, maxStclPr, rowPtr, colInd, vals, nnzPr);
    rdNodes.set_size(0);

    // Provides a portable wall clock timing routine.
    tEnd = 0.0;

#ifdef _OPENMP

    tEnd = omp_get_wtime();

#endif //_OPENMP

    if (params->verbose > 1) {
      m2cPrintf(" Init or updated OSUS operator in %gs...\n", tEnd - tStart);
      fflush(stdout);
    }

    // Provides a portable wall clock timing routine.
    tStart = 0.0;

#ifdef _OPENMP

    tStart = omp_get_wtime();

#endif //_OPENMP

    if (isSurf) {
      //  surface assembly
      assemble_surf(mesh->xs, mesh->conn, mesh->dirs, mesh->n2cPtr,
                    mesh->n2cList, mesh->parts, stcls, rowPtr, colInd, vals,
                    nnzPr, rdNodes, params);
    } else {
      //  body assembly
      assemble_body(mesh->xs, mesh->conn, mesh->n2cPtr, mesh->n2cList,
                    mesh->parts, stcls, rowPtr, colInd, vals, nnzPr, rdNodes,
                    params);
    }

    // Provides a portable wall clock timing routine.
    tEnd = 0.0;

#ifdef _OPENMP

    tEnd = omp_get_wtime();

#endif //_OPENMP

    if (params->verbose > 1) {
      m2cPrintf(" Assembly OSUS operator finished in %gs...\n", tEnd - tStart);
      fflush(stdout);
    }

    if (rdNodes.size(0) == 0) {
      int b_i;
      int loop_ub;
      int n;

      n = rowPtr.size(0) - 2;
      b_rowPtr.set_size(rowPtr.size(0));
      b_rowPtr[0] = 1;
      for (int i{0}; i <= n; i++) {
        b_rowPtr[i + 1] = b_rowPtr[i] + nnzPr[i];
      }

      b_colInd.set_size(b_rowPtr[rowPtr.size(0) - 1] - 1);
      b_vals.set_size(b_rowPtr[rowPtr.size(0) - 1] - 1);
      for (int i{0}; i <= n; i++) {
        int j;
        j = rowPtr[i] - 1;
        b_i = b_rowPtr[i];
        loop_ub = b_rowPtr[i + 1] - 1;
        for (int k{b_i}; k <= loop_ub; k++) {
          int colInd_tmp;
          colInd_tmp = (j + k) - b_i;
          b_colInd[k - 1] = colInd[colInd_tmp];
          b_vals[k - 1] = vals[colInd_tmp];
        }
      }

      rowPtr.set_size(b_rowPtr.size(0));
      loop_ub = b_rowPtr.size(0);
      for (b_i = 0; b_i < loop_ub; b_i++) {
        rowPtr[b_i] = b_rowPtr[b_i];
      }

      colInd.set_size(b_colInd.size(0));
      loop_ub = b_colInd.size(0);
      for (b_i = 0; b_i < loop_ub; b_i++) {
        colInd[b_i] = b_colInd[b_i];
      }

      vals.set_size(b_vals.size(0));
      loop_ub = b_vals.size(0);
      for (b_i = 0; b_i < loop_ub; b_i++) {
        vals[b_i] = b_vals[b_i];
      }
    }
  }

  void rdi_build_node2cell(int n, const ::coder::array<int, 2U> &conn, ::coder::
    array<int, 1U> &n2cPtr, ::coder::array<int, 1U> &n2cList)
  {
    int idx;
    int m;
    int npts;

    // rdi_build_node2cell - Initilize an unstructured mesh
    n2cPtr.set_size(n + 1);
    for (int i{0}; i <= n; i++) {
      n2cPtr[i] = 0;
    }

    n2cPtr[0] = 1;

    //  determine number of incident cells
    m = conn.size(0) - 1;

    //  number of cells
    for (int e{0}; e <= m; e++) {
      for (npts = conn.size(1) - 1; conn[npts + conn.size(1) * e] <= 0; npts--)
      {
      }

      for (int b_i{0}; b_i <= npts; b_i++) {
        idx = conn[b_i + conn.size(1) * e];
        n2cPtr[idx] = n2cPtr[idx] + 1;
      }
    }

    for (int b_i{0}; b_i < n; b_i++) {
      n2cPtr[b_i + 1] = n2cPtr[b_i + 1] + n2cPtr[b_i];
    }

    //  allocate n2cList
    n2cList.set_size(n2cPtr[n] - 1);
    for (int e{0}; e <= m; e++) {
      for (npts = conn.size(1) - 1; conn[npts + conn.size(1) * e] <= 0; npts--)
      {
      }

      for (int b_i{0}; b_i <= npts; b_i++) {
        idx = conn[b_i + conn.size(1) * e] - 1;
        n2cList[n2cPtr[idx] - 1] = e + 1;
        n2cPtr[idx] = n2cPtr[idx] + 1;
      }
    }

    for (int b_i{n}; b_i >= 1; b_i--) {
      n2cPtr[b_i] = n2cPtr[b_i - 1];
    }

    n2cPtr[0] = 1;
  }

  void rdi_build_node2node(int n, const ::coder::array<int, 2U> &conn, int dim,
    const ::coder::array<int, 1U> &n2cPtr, const ::coder::array<int, 1U>
    &n2cList, ::coder::array<int, 1U> &n2nPtr, ::coder::array<int, 1U> &n2nList)
  {
    static const signed char next[8]{ 2, 3, 1, 0, 2, 3, 4, 1 };

    static const signed char prev[8]{ 3, 1, 2, 0, 4, 1, 2, 3 };

    ::coder::array<boolean_T, 1U> visited;
    int nodes[256];
    int i;

    // rdi_build_node2node - Compute node-to-node adjacency (1-ring)
    n2nPtr.set_size(n + 1);
    n2nPtr[0] = 1;
    visited.set_size(n);
    for (i = 0; i < n; i++) {
      visited[i] = false;
    }

    if ((dim == 1) || (dim == 3) || ((dim == 2) && (conn.size(1) == 3))) {
      int i1;
      int nPc;
      int nPoints;
      int nid;

      //  for 1D, 3D (tet only), and 2D triangles, we can simply determine
      nPc = conn.size(1) - 1;
      for (int b_i{0}; b_i < n; b_i++) {
        nPoints = -1;
        i = n2cPtr[b_i];
        i1 = n2cPtr[b_i + 1] - 1;
        for (int j{i}; j <= i1; j++) {
          for (int k{0}; k <= nPc; k++) {
            nid = conn[k + conn.size(1) * (n2cList[j - 1] - 1)];
            if ((nid != b_i + 1) && (!visited[nid - 1])) {
              nPoints++;
              nodes[nPoints] = nid;
              visited[nid - 1] = true;
            }
          }
        }

        n2nPtr[b_i + 1] = (n2nPtr[b_i] + nPoints) + 1;

        //  reset tags
        for (int j{0}; j <= nPoints; j++) {
          visited[nodes[j] - 1] = false;
        }
      }

      n2nList.set_size(n2nPtr[n] - 1);
      for (int b_i{0}; b_i < n; b_i++) {
        nPoints = -1;
        i = n2cPtr[b_i];
        i1 = n2cPtr[b_i + 1] - 1;
        for (int j{i}; j <= i1; j++) {
          for (int k{0}; k <= nPc; k++) {
            nid = conn[k + conn.size(1) * (n2cList[j - 1] - 1)];
            if ((nid != b_i + 1) && (!visited[nid - 1])) {
              n2nList[n2nPtr[b_i] + nPoints] = nid;
              nPoints++;
              visited[nid - 1] = true;
            }
          }
        }

        //  reset tags
        i = n2nPtr[b_i];
        i1 = n2nPtr[b_i + 1] - 1;
        for (int j{i}; j <= i1; j++) {
          visited[n2nList[j - 1] - 1] = false;
        }
      }
    } else {
      int i1;
      int i2;
      int k;
      int nPc;
      int nPoints;
      int nidNext;
      int nidNext_tmp;
      int nidPrev;

      //  2D quad or mixed
      for (int b_i{0}; b_i < n; b_i++) {
        nPoints = -1;
        i = n2cPtr[b_i];
        i1 = n2cPtr[b_i + 1] - 1;
        for (int j{i}; j <= i1; j++) {
          nPc = 1;
          i2 = n2cList[j - 1] - 1;
          if (conn[conn.size(1) * i2 + 3] <= 0) {
            nPc = 0;
          }

          for (k = 0; conn[k + conn.size(1) * i2] != b_i + 1; k++) {
          }

          nidNext_tmp = k + (nPc << 2);
          nidNext = conn[(next[nidNext_tmp] + conn.size(1) * i2) - 1];
          nidPrev = conn[(prev[nidNext_tmp] + conn.size(1) * i2) - 1];
          if (!visited[nidNext - 1]) {
            nPoints++;
            nodes[nPoints] = nidNext;
            visited[nidNext - 1] = true;
          }

          if (!visited[nidPrev - 1]) {
            nPoints++;
            nodes[nPoints] = nidPrev;
            visited[nidPrev - 1] = true;
          }
        }

        n2nPtr[b_i + 1] = (n2nPtr[b_i] + nPoints) + 1;

        //  reset tags
        for (int j{0}; j <= nPoints; j++) {
          visited[nodes[j] - 1] = false;
        }
      }

      n2nList.set_size(n2nPtr[n] - 1);
      for (int b_i{0}; b_i < n; b_i++) {
        nPoints = -1;
        i = n2cPtr[b_i];
        i1 = n2cPtr[b_i + 1] - 1;
        for (int j{i}; j <= i1; j++) {
          nPc = 1;
          i2 = n2cList[j - 1] - 1;
          if (conn[conn.size(1) * i2 + 3] <= 0) {
            nPc = 0;
          }

          for (k = 0; conn[k + conn.size(1) * i2] != b_i + 1; k++) {
          }

          nidNext_tmp = k + (nPc << 2);
          nidNext = conn[(next[nidNext_tmp] + conn.size(1) * i2) - 1];
          nidPrev = conn[(prev[nidNext_tmp] + conn.size(1) * i2) - 1];
          if (!visited[nidNext - 1]) {
            n2nList[n2nPtr[b_i] + nPoints] = nidNext;
            nPoints++;
            visited[nidNext - 1] = true;
          }

          if (!visited[nidPrev - 1]) {
            n2nList[n2nPtr[b_i] + nPoints] = nidPrev;
            nPoints++;
            visited[nidPrev - 1] = true;
          }
        }

        //  reset tags
        i = n2nPtr[b_i];
        i1 = n2nPtr[b_i + 1] - 1;
        for (int j{i}; j <= i1; j++) {
          visited[n2nList[j - 1] - 1] = false;
        }
      }
    }
  }

  void rdi_compute_cellsizes(const ::coder::array<double, 2U> &xs, const ::
    coder::array<int, 2U> &conn, const ::coder::array<int, 1U> &n2nPtr, const ::
    coder::array<int, 1U> &n2nList, const RdiParams *params, const ::coder::
    array<double, 2U> &nrms, ::coder::array<double, 1U> &h)
  {
    ::coder::array<double, 1U> buf;

    // rdi_compute_cellsizes - Compute cell sizes (h)
    if (h.size(0) == 0) {

#pragma omp single

      { //single
        h.set_size(conn.size(0));

      } //single
    }

    //  rdi_assert(size(h, 1) >= size(conn, 1), 'insufficient memory');
    if (xs.size(1) == params->dim) {
      //  body
      if (params->verbose > 1) {
        m2cPrintf(" Compute cell sizes for body mesh...\n");
        fflush(stdout);
      }

      compute_body_h(xs, conn, h, params->nThreads);
    } else if (xs.size(1) == params->dim + 1) {
      //  surface
      if ((params->surfType == 0) && ((nrms.size(0) == 0) || (nrms.size(1) == 0)))
      {
        m2cErrMsgIdAndTxt("rdi_compute_cellsizes:missingNorms",
                          "normal vectors must be passed in");
      }

      if (params->verbose > 1) {
        m2cPrintf(" Compute cell sizes for surface mesh (type=%d)...\n",
               params->surfType);
        fflush(stdout);
      }

      compute_surf_h(xs, conn, n2nPtr, n2nList, nrms, params->surfType, h,
                     params->nThreads);
    } else if (xs.size(1) == params->dim + 2) {
      int i;
      int loop_ub;
      int m;
      int n;

      //  3D curve
      if ((nrms.size(0) == 0) || (nrms.size(1) == 0)) {
        m2cErrMsgIdAndTxt("rdi_compute_cellsizes:missingTangent",
                          "for 3D curve, the tangent vectors must be passed in");
      }

      if (params->verbose > 1) {
        m2cPrintf(" Compute cell sizes for curve mesh in 3D...\n");
        fflush(stdout);
      }

      //  use nrms to store tangent vectors at each node
      m = conn.size(0);
      n = xs.size(0);
      buf.set_size(xs.size(0));
      loop_ub = xs.size(0);
      for (i = 0; i < loop_ub; i++) {
        buf[i] = 0.0;
      }

      for (int b_i{0}; b_i < n; b_i++) {
        double t1;
        double t2;
        double t3;
        int i1;
        t1 = nrms[nrms.size(1) * b_i];
        t2 = nrms[nrms.size(1) * b_i + 1];
        t3 = nrms[nrms.size(1) * b_i + 2];
        i = n2nPtr[b_i];
        loop_ub = n2nPtr[b_i + 1];
        i1 = loop_ub - 1;
        for (int j{i}; j <= i1; j++) {
          int buf_tmp;
          buf_tmp = n2nList[j - 1] - 1;
          buf[b_i] = buf[b_i] + std::abs(((xs[xs.size(1) * buf_tmp] - xs[xs.size
            (1) * b_i]) * t1 + (xs[xs.size(1) * buf_tmp + 1] - xs[xs.size(1) *
                                b_i + 1]) * t2) + (xs[xs.size(1) * buf_tmp + 2]
            - xs[xs.size(1) * b_i + 2]) * t3);
        }

        buf[b_i] = buf[b_i] / static_cast<double>(loop_ub - n2nPtr[b_i]);
      }

      for (int e{0}; e < m; e++) {
        h[e] = 0.5 * (buf[conn[conn.size(1) * e + 1] - 1] + buf[conn[conn.size(1)
                      * e] - 1]);
      }
    } else {
      m2cErrMsgIdAndTxt("rdi_compute_cellsizes:badDims",
                        "invalid dim(%d) and topoDim(%d)", xs.size(1),
                        params->dim);
    }
  }

  void rdi_compute_cellsizes2(const ::coder::array<double, 2U> &xs, const ::
    coder::array<int, 2U> &conn, const ::coder::array<int, 1U> &n2nPtr, const ::
    coder::array<int, 1U> &n2nList, const RdiParams *params, const ::coder::
    array<double, 2U> &nrms, ::coder::array<double, 1U> &h)
  {
    ::coder::array<double, 1U> buf;

    // rdi_compute_cellsizes - Compute cell sizes (h)
#pragma omp single

    { //single
      h.set_size(conn.size(0));

    } //single

    //  rdi_assert(size(h, 1) >= size(conn, 1), 'insufficient memory');
    if (xs.size(1) == params->dim) {
      //  body
      if (params->verbose > 1) {
        m2cPrintf(" Compute cell sizes for body mesh...\n");
        fflush(stdout);
      }

      compute_body_h(xs, conn, h, params->nThreads);
    } else if (xs.size(1) == params->dim + 1) {
      //  surface
      if ((params->surfType == 0) && ((nrms.size(0) == 0) || (nrms.size(1) == 0)))
      {
        m2cErrMsgIdAndTxt("rdi_compute_cellsizes:missingNorms",
                          "normal vectors must be passed in");
      }

      if (params->verbose > 1) {
        m2cPrintf(" Compute cell sizes for surface mesh (type=%d)...\n",
               params->surfType);
        fflush(stdout);
      }

      compute_surf_h(xs, conn, n2nPtr, n2nList, nrms, params->surfType, h,
                     params->nThreads);
    } else if (xs.size(1) == params->dim + 2) {
      int i;
      int loop_ub;
      int m;
      int n;

      //  3D curve
      if ((nrms.size(0) == 0) || (nrms.size(1) == 0)) {
        m2cErrMsgIdAndTxt("rdi_compute_cellsizes:missingTangent",
                          "for 3D curve, the tangent vectors must be passed in");
      }

      if (params->verbose > 1) {
        m2cPrintf(" Compute cell sizes for curve mesh in 3D...\n");
        fflush(stdout);
      }

      //  use nrms to store tangent vectors at each node
      m = conn.size(0);
      n = xs.size(0);
      buf.set_size(xs.size(0));
      loop_ub = xs.size(0);
      for (i = 0; i < loop_ub; i++) {
        buf[i] = 0.0;
      }

      for (int b_i{0}; b_i < n; b_i++) {
        double t1;
        double t2;
        double t3;
        int i1;
        t1 = nrms[nrms.size(1) * b_i];
        t2 = nrms[nrms.size(1) * b_i + 1];
        t3 = nrms[nrms.size(1) * b_i + 2];
        i = n2nPtr[b_i];
        loop_ub = n2nPtr[b_i + 1];
        i1 = loop_ub - 1;
        for (int j{i}; j <= i1; j++) {
          int buf_tmp;
          buf_tmp = n2nList[j - 1] - 1;
          buf[b_i] = buf[b_i] + std::abs(((xs[xs.size(1) * buf_tmp] - xs[xs.size
            (1) * b_i]) * t1 + (xs[xs.size(1) * buf_tmp + 1] - xs[xs.size(1) *
                                b_i + 1]) * t2) + (xs[xs.size(1) * buf_tmp + 2]
            - xs[xs.size(1) * b_i + 2]) * t3);
        }

        buf[b_i] = buf[b_i] / static_cast<double>(loop_ub - n2nPtr[b_i]);
      }

      for (int e{0}; e < m; e++) {
        h[e] = 0.5 * (buf[conn[conn.size(1) * e + 1] - 1] + buf[conn[conn.size(1)
                      * e] - 1]);
      }
    } else {
      m2cErrMsgIdAndTxt("rdi_compute_cellsizes:badDims",
                        "invalid dim(%d) and topoDim(%d)", xs.size(1),
                        params->dim);
    }
  }

  static inline
  void rdi_compute_cellsizes(const ::coder::array<double, 2U> &xs, const ::
    coder::array<int, 2U> &conn, const ::coder::array<int, 1U> &n2nPtr, const ::
    coder::array<int, 1U> &n2nList, const RdiParams *params, ::coder::array<
    double, 1U> &h)
  {
    // rdi_compute_cellsizes - Compute cell sizes (h)

#pragma omp single

    { //single
      h.set_size(conn.size(0));

    } //single

    //  rdi_assert(size(h, 1) >= size(conn, 1), 'insufficient memory');
    if (xs.size(1) == params->dim) {
      //  body
      if (params->verbose > 1) {
        m2cPrintf(" Compute cell sizes for body mesh...\n");
        fflush(stdout);
      }

      compute_body_h(xs, conn, h, params->nThreads);
    } else if (xs.size(1) == params->dim + 1) {
      //  surface
      if (params->surfType == 0) {
        m2cErrMsgIdAndTxt("rdi_compute_cellsizes:missingNorms",
                          "normal vectors must be passed in");
      }

      if (params->verbose > 1) {
        m2cPrintf(" Compute cell sizes for surface mesh (type=%d)...\n",
               params->surfType);
        fflush(stdout);
      }

      compute_surf_h(xs, conn, n2nPtr, n2nList, params->surfType, h,
                     params->nThreads);
    } else if (xs.size(1) == params->dim + 2) {
      int m;

      //  3D curve
      m2cErrMsgIdAndTxt("rdi_compute_cellsizes:missingTangent",
                        "for 3D curve, the tangent vectors must be passed in");
      if (params->verbose > 1) {
        m2cPrintf(" Compute cell sizes for curve mesh in 3D...\n");
        fflush(stdout);
      }

      //  use nrms to store tangent vectors at each node
      m = conn.size(0);
      for (int e{0}; e < m; e++) {
        h[e] = 0.0;
      }
    } else {
      m2cErrMsgIdAndTxt("rdi_compute_cellsizes:badDims",
                        "invalid dim(%d) and topoDim(%d)", xs.size(1),
                        params->dim);
    }
  }

  static inline
  void rdi_compute_cellweights(const ::coder::array<double, 2U> &xs, const ::
    coder::array<int, 2U> &conn, const RdiParams *params, ::coder::array<double,
    1U> &w)
  {
    // rdi_compute_cellweights - Compute the cell weights (area, volume, etc.)
    if (w.size(0) == 0) {

#pragma omp single

      { //single
        w.set_size(conn.size(0));

      } //single
    }

    //  rdi_assert(size(w, 1) >= size(conn, 1), 'insufficient size of w');
    if (params->dim == 3) {
      if (conn.size(1) != 4) {
        m2cErrMsgIdAndTxt("rdi_compute_cellweights:unsupportedCell",
                          "only linear tet meshes are supported in 3D");
      }

      if (params->verbose > 1) {
        m2cPrintf(" Compute cell weights for 3D...\n");
        fflush(stdout);
      }

      compute_volume_tet(xs, conn, w, params->nThreads);
    } else if (params->dim == 2) {
      if ((conn.size(1) != 3) && (conn.size(1) != 4)) {
        m2cErrMsgIdAndTxt("rdi_compute_cellweights:unsupportedCell",
                          "must be linear tri or bilinear quad");
      }

      if (params->verbose > 1) {
        m2cPrintf(" Compute cell weights for (topo) 2D...\n");
        fflush(stdout);
      }

      compute_area(xs, conn, w, params->nThreads);
    } else {
      int m;
      if (conn.size(1) != 2) {
        m2cErrMsgIdAndTxt("rdi_compute_cellweights:unsupportedCell",
                          "must be linear bars");
      }

      if (params->verbose > 1) {
        m2cPrintf(" Compute cell weights for (topo) 1D...\n");
        fflush(stdout);
      }

      //  kernel for computing length of bars
      m = conn.size(0) - 1;
      if (xs.size(1) == 1) {
        for (int e{0}; e <= m; e++) {
          int i;
          int i1;
          i = conn[conn.size(1) * e + 1] - 1;
          i1 = conn[conn.size(1) * e] - 1;
          w[e] = std::abs(xs[i % xs.size(0) + i / xs.size(0)] - xs[i1 % xs.size
                          (0) + i1 / xs.size(0)]);
        }
      } else if (xs.size(1) == 2) {
        for (int e{0}; e <= m; e++) {
          double a_idx_0;
          double a_idx_1;
          int i;
          int i1;
          i = conn[conn.size(1) * e + 1];
          i1 = conn[conn.size(1) * e];
          a_idx_0 = xs[xs.size(1) * (i - 1)] - xs[xs.size(1) * (i1 - 1)];
          a_idx_1 = xs[xs.size(1) * (i - 1) + 1] - xs[xs.size(1) * (i1 - 1) + 1];
          w[e] = std::sqrt(a_idx_0 * a_idx_0 + a_idx_1 * a_idx_1);
        }
      } else {
        for (int e{0}; e <= m; e++) {
          double a_idx_0;
          double a_idx_1;
          double a_idx_2;
          int i;
          int i1;
          i = conn[conn.size(1) * e + 1];
          i1 = conn[conn.size(1) * e];
          a_idx_0 = xs[xs.size(1) * (i - 1)] - xs[xs.size(1) * (i1 - 1)];
          a_idx_1 = xs[xs.size(1) * (i - 1) + 1] - xs[xs.size(1) * (i1 - 1) + 1];
          a_idx_2 = xs[xs.size(1) * (i - 1) + 2] - xs[xs.size(1) * (i1 - 1) + 2];
          w[e] = std::sqrt((a_idx_0 * a_idx_0 + a_idx_1 * a_idx_1) + a_idx_2 *
                           a_idx_2);
        }
      }
    }
  }

  static inline
  void rdi_compute_cellweights2(const ::coder::array<double, 2U> &xs, const ::
    coder::array<int, 2U> &conn, const RdiParams *params, ::coder::array<double,
    1U> &w)
  {
    // rdi_compute_cellweights - Compute the cell weights (area, volume, etc.)

#pragma omp single

    { //single
      w.set_size(conn.size(0));

    } //single

    //  rdi_assert(size(w, 1) >= size(conn, 1), 'insufficient size of w');
    if (params->dim == 3) {
      if (conn.size(1) != 4) {
        m2cErrMsgIdAndTxt("rdi_compute_cellweights:unsupportedCell",
                          "only linear tet meshes are supported in 3D");
      }

      if (params->verbose > 1) {
        m2cPrintf(" Compute cell weights for 3D...\n");
        fflush(stdout);
      }

      compute_volume_tet(xs, conn, w, params->nThreads);
    } else if (params->dim == 2) {
      if ((conn.size(1) != 3) && (conn.size(1) != 4)) {
        m2cErrMsgIdAndTxt("rdi_compute_cellweights:unsupportedCell",
                          "must be linear tri or bilinear quad");
      }

      if (params->verbose > 1) {
        m2cPrintf(" Compute cell weights for (topo) 2D...\n");
        fflush(stdout);
      }

      compute_area(xs, conn, w, params->nThreads);
    } else {
      int m;
      if (conn.size(1) != 2) {
        m2cErrMsgIdAndTxt("rdi_compute_cellweights:unsupportedCell",
                          "must be linear bars");
      }

      if (params->verbose > 1) {
        m2cPrintf(" Compute cell weights for (topo) 1D...\n");
        fflush(stdout);
      }

      //  kernel for computing length of bars
      m = conn.size(0) - 1;
      if (xs.size(1) == 1) {
        for (int e{0}; e <= m; e++) {
          int i;
          int i1;
          i = conn[conn.size(1) * e + 1] - 1;
          i1 = conn[conn.size(1) * e] - 1;
          w[e] = std::abs(xs[i % xs.size(0) + i / xs.size(0)] - xs[i1 % xs.size
                          (0) + i1 / xs.size(0)]);
        }
      } else if (xs.size(1) == 2) {
        for (int e{0}; e <= m; e++) {
          double a_idx_0;
          double a_idx_1;
          int i;
          int i1;
          i = conn[conn.size(1) * e + 1];
          i1 = conn[conn.size(1) * e];
          a_idx_0 = xs[xs.size(1) * (i - 1)] - xs[xs.size(1) * (i1 - 1)];
          a_idx_1 = xs[xs.size(1) * (i - 1) + 1] - xs[xs.size(1) * (i1 - 1) + 1];
          w[e] = std::sqrt(a_idx_0 * a_idx_0 + a_idx_1 * a_idx_1);
        }
      } else {
        for (int e{0}; e <= m; e++) {
          double a_idx_0;
          double a_idx_1;
          double a_idx_2;
          int i;
          int i1;
          i = conn[conn.size(1) * e + 1];
          i1 = conn[conn.size(1) * e];
          a_idx_0 = xs[xs.size(1) * (i - 1)] - xs[xs.size(1) * (i1 - 1)];
          a_idx_1 = xs[xs.size(1) * (i - 1) + 1] - xs[xs.size(1) * (i1 - 1) + 1];
          a_idx_2 = xs[xs.size(1) * (i - 1) + 2] - xs[xs.size(1) * (i1 - 1) + 2];
          w[e] = std::sqrt((a_idx_0 * a_idx_0 + a_idx_1 * a_idx_1) + a_idx_2 *
                           a_idx_2);
        }
      }
    }
  }

  static inline
  void rdi_compute_inds(const ::coder::array<int, 1U> &rowPtr, const ::coder::
    array<int, 1U> &colInd, const ::coder::array<double, 1U> &vals, const ::
    coder::array<double, 2U> &fs, const ::coder::array<double, 2U> &dfGlobal,
    const RdiMesh *mesh, const RdiParams *params, ::coder::array<double, 2U>
    &alphaCell, ::coder::array<double, 2U> &alphaNode, ::coder::array<double, 2U>
    &beta)
  {
    // rdi_compute_inds - Compute all indicators at once
    rdi_compute_osusind(rowPtr, colInd, vals, fs, mesh->n2cPtr, mesh->n2cList,
                        alphaCell, alphaNode);

    //  ompBarrier; % No need to have barrier here
    rdi_compute_oscind(dfGlobal, alphaCell, mesh->xs, mesh->hGlobal,
                       mesh->cellWeights, mesh->n2cPtr, mesh->n2cList,
                       params->dim, params->epsBeta, beta);
  }

  void rdi_compute_inds2(const ::coder::array<int, 1U> &rowPtr, const ::coder::
    array<int, 1U> &colInd, const ::coder::array<double, 1U> &vals, const ::
    coder::array<double, 2U> &fs, const ::coder::array<double, 2U> &dfGlobal,
    const RdiMesh *mesh, const RdiParams *params, ::coder::array<double, 2U>
    &alphaCell, ::coder::array<double, 2U> &alphaNode, ::coder::array<double, 2U>
    &beta)
  {
    int m2cTryBlkErrCode;
    int nthreads;

    // rdi_compute_inds - Compute all indicators at once
    nthreads = params->nThreads;

    m2cTryBlkErrCode = (0);
    if (params->nThreads <= 0) {
      nthreads = 1;

#ifdef _OPENMP

      nthreads = omp_get_max_threads();

#endif //_OPENMP

    }

#pragma omp parallel default(shared) num_threads(nthreads)

    try { //try
      rdi_compute_osusind(rowPtr, colInd, vals, fs, mesh->n2cPtr, mesh->n2cList,
                          alphaCell, alphaNode);

      //  ompBarrier; % No need to have barrier here
      rdi_compute_oscind(dfGlobal, alphaCell, mesh->xs, mesh->hGlobal,
                         mesh->cellWeights, mesh->n2cPtr, mesh->n2cList,
                         params->dim, params->epsBeta, beta);

    } catch (const std::runtime_error& m2cExc) {
      m2cTryBlkErrCode = (1);
      m2cPrintError("runtime_error %s\n", m2cExc.what());

    }

    catch (const std::logic_error& m2cExc)
    {
      m2cTryBlkErrCode = (1);
      m2cPrintError("logic_error %s\n", m2cExc.what());

    }

    catch (const std::exception& m2cExc)
    {
      m2cTryBlkErrCode = (1);
      m2cPrintError("exception %s\n", m2cExc.what());

    }

    catch (...)
    {
      m2cTryBlkErrCode = (1);
      m2cPrintError("Unknown error detected from C++ exceptions\n");

    } //try

    if ((int)m2cTryBlkErrCode != 0) {
      throw std::runtime_error("omp4m:runtimeErrorInThread");
    }

  }

  void rdi_compute_oscind(const ::coder::array<double, 2U> &dfGlobal, const ::
    coder::array<double, 2U> &alphaCell, const RdiMesh *mesh, const RdiParams
    *params, ::coder::array<double, 2U> &beta)
  {
    double epsBeta;
    double hGlobal;

    // rdi_compute_oscind - Compute oscillation indicators (beta values)
    epsBeta = params->epsBeta;
    if (params->epsBeta <= 0.0) {
      epsBeta = 0.0002;
    }

    hGlobal = mesh->hGlobal;
    if (mesh->hGlobal <= 0.0) {
      double ex;
      int last;
      last = mesh->cellWeights.size(0);
      if (mesh->cellWeights.size(0) <= 2) {
        if (mesh->cellWeights.size(0) == 1) {
          ex = mesh->cellWeights[0];
        } else if (mesh->cellWeights[0] < mesh->cellWeights
                   [mesh->cellWeights.size(0) - 1]) {
          ex = mesh->cellWeights[mesh->cellWeights.size(0) - 1];
        } else {
          ex = mesh->cellWeights[0];
        }
      } else {
        ex = mesh->cellWeights[0];
        for (int k{2}; k <= last; k++) {
          double d;
          d = mesh->cellWeights[k - 1];
          if (ex < d) {
            ex = d;
          }
        }
      }

      hGlobal = std::pow(ex, 1.0 / static_cast<double>(params->dim));
    }

    compute_beta_kernel(mesh->xs.size(0), epsBeta, hGlobal, dfGlobal, alphaCell,
                        mesh->cellWeights, mesh->n2cPtr, mesh->n2cList,
                        alphaCell.size(1), beta);
  }

  void rdi_compute_oscind2(const ::coder::array<double, 2U> &dfGlobal, const ::
    coder::array<double, 2U> &alphaCell, const RdiMesh *mesh, const RdiParams
    *params, ::coder::array<double, 2U> &beta)
  {
    int m2cTryBlkErrCode;
    double epsBeta;
    double hGlobal;
    int last;

    // rdi_compute_oscind - Compute oscillation indicators (beta values)
    epsBeta = params->epsBeta;
    if (params->epsBeta <= 0.0) {
      epsBeta = 0.0002;
    }

    hGlobal = mesh->hGlobal;
    if (mesh->hGlobal <= 0.0) {
      double ex;
      last = mesh->cellWeights.size(0);
      if (mesh->cellWeights.size(0) <= 2) {
        if (mesh->cellWeights.size(0) == 1) {
          ex = mesh->cellWeights[0];
        } else if (mesh->cellWeights[0] < mesh->cellWeights
                   [mesh->cellWeights.size(0) - 1]) {
          ex = mesh->cellWeights[mesh->cellWeights.size(0) - 1];
        } else {
          ex = mesh->cellWeights[0];
        }
      } else {
        ex = mesh->cellWeights[0];
        for (int k{2}; k <= last; k++) {
          double d;
          d = mesh->cellWeights[k - 1];
          if (ex < d) {
            ex = d;
          }
        }
      }

      hGlobal = std::pow(ex, 1.0 / static_cast<double>(params->dim));
    }

    last = params->nThreads;

    m2cTryBlkErrCode = (0);
    if (params->nThreads <= 0) {
      last = 1;

#ifdef _OPENMP

      last = omp_get_max_threads();

#endif //_OPENMP

    }

#pragma omp parallel default(shared) num_threads(last)

    try { //try
      compute_beta_kernel(mesh->xs.size(0), epsBeta, hGlobal, dfGlobal,
                          alphaCell, mesh->cellWeights, mesh->n2cPtr,
                          mesh->n2cList, alphaCell.size(1), beta);

    } catch (const std::runtime_error& m2cExc) {
      m2cTryBlkErrCode = (1);
      m2cPrintError("runtime_error %s\n", m2cExc.what());

    }

    catch (const std::logic_error& m2cExc)
    {
      m2cTryBlkErrCode = (1);
      m2cPrintError("logic_error %s\n", m2cExc.what());

    }

    catch (const std::exception& m2cExc)
    {
      m2cTryBlkErrCode = (1);
      m2cPrintError("exception %s\n", m2cExc.what());

    }

    catch (...)
    {
      m2cTryBlkErrCode = (1);
      m2cPrintError("Unknown error detected from C++ exceptions\n");

    } //try

    if ((int)m2cTryBlkErrCode != 0) {
      throw std::runtime_error("omp4m:runtimeErrorInThread");
    }

  }

  static inline
  void rdi_compute_osusind(const ::coder::array<int, 1U> &rowPtr, const ::coder::
    array<int, 1U> &colInd, const ::coder::array<double, 1U> &vals, const ::
    coder::array<double, 2U> &fs, const RdiMesh *mesh, const RdiParams *, ::
    coder::array<double, 2U> &alphaCell, ::coder::array<double, 2U> &alphaNode)
  {
    // rdi_compute_osusind - Compute over-/under-shoot indicators
    crs_prod_mat_vec(rowPtr, colInd, vals, fs, alphaCell);
    compute_nodal_alpha(fs.size(0), alphaCell, mesh->n2cPtr, mesh->n2cList,
                        fs.size(1), alphaNode);
  }

  void rdi_compute_osusind2(const ::coder::array<int, 1U> &rowPtr, const ::coder::
    array<int, 1U> &colInd, const ::coder::array<double, 1U> &vals, const ::
    coder::array<double, 2U> &fs, const RdiMesh *mesh, const RdiParams *params, ::
    coder::array<double, 2U> &alphaCell, ::coder::array<double, 2U> &alphaNode)
  {
    int m2cTryBlkErrCode;
    int nthreads;

    // rdi_compute_osusind - Compute over-/under-shoot indicators
    nthreads = params->nThreads;

    m2cTryBlkErrCode = (0);
    if (params->nThreads <= 0) {
      nthreads = 1;

#ifdef _OPENMP

      nthreads = omp_get_max_threads();

#endif //_OPENMP

    }

#pragma omp parallel default(shared) num_threads(nthreads)

    try { //try
      //  using advanced OpenMP mode
      crs_prod_mat_vec(rowPtr, colInd, vals, fs, alphaCell);
      compute_nodal_alpha(fs.size(0), alphaCell, mesh->n2cPtr, mesh->n2cList,
                          fs.size(1), alphaNode);

    } catch (const std::runtime_error& m2cExc) {
      m2cTryBlkErrCode = (1);
      m2cPrintError("runtime_error %s\n", m2cExc.what());

    }

    catch (const std::logic_error& m2cExc)
    {
      m2cTryBlkErrCode = (1);
      m2cPrintError("logic_error %s\n", m2cExc.what());

    }

    catch (const std::exception& m2cExc)
    {
      m2cTryBlkErrCode = (1);
      m2cPrintError("exception %s\n", m2cExc.what());

    }

    catch (...)
    {
      m2cTryBlkErrCode = (1);
      m2cPrintError("Unknown error detected from C++ exceptions\n");

    } //try

    if ((int)m2cTryBlkErrCode != 0) {
      throw std::runtime_error("omp4m:runtimeErrorInThread");
    }

  }

  static inline
  void rdi_default_params(int dim, int nThreads, RdiParams *params)
  {
    // rdi_default_params - Create a structure of default parameters
    if ((dim < 1) || (dim > 3)) {
      m2cErrMsgIdAndTxt("rdi_default_params:badDim",
                        "wrong (topological) dimension %d", dim);
    }

    params->dim = dim;
    params->ring = 2.0;
    params->maxStclSize = 0;
    params->verbose = 1;
    params->degree = 2;
    params->surfType = 0;
    params->epsBeta = 0.0002;
    params->cGlobal = 0.08;
    params->cLocal = 0.5;
    params->kappa1 = 0.3;
    params->kappa0 = 1.0;
    params->markNearDis = true;
    params->wlsInterp0 = false;
    params->wlsUseDag = false;
    if (nThreads < 1) {
      nThreads = 1;

#ifdef _OPENMP

      nThreads = omp_get_max_threads();

#endif //_OPENMP

    }

    params->nThreads = nThreads;
    params->parTask = true;
  }

  static inline
  void rdi_default_params2(int dim, int nThreads, RdiParams *params)
  {
    // rdi_default_params - Create a structure of default parameters
    if ((dim < 1) || (dim > 3)) {
      m2cErrMsgIdAndTxt("rdi_default_params:badDim",
                        "wrong (topological) dimension %d", dim);
    }

    params->dim = dim;
    params->ring = 2.0;
    params->maxStclSize = 0;
    params->verbose = 1;
    params->degree = 2;
    params->surfType = 0;
    params->epsBeta = 0.0002;
    params->cGlobal = 0.08;
    params->cLocal = 0.5;
    params->kappa1 = 0.3;
    params->kappa0 = 1.0;
    params->markNearDis = true;
    params->wlsInterp0 = false;
    params->wlsUseDag = false;
    if (nThreads < 1) {
      nThreads = 1;

#ifdef _OPENMP

      nThreads = omp_get_max_threads();

#endif //_OPENMP

    }

    params->nThreads = nThreads;
    params->parTask = true;
  }

  void rdi_default_params(int dim, RdiParams *params)
  {
    int nThreads;

    // rdi_default_params - Create a structure of default parameters
    if ((dim < 1) || (dim > 3)) {
      m2cErrMsgIdAndTxt("rdi_default_params:badDim",
                        "wrong (topological) dimension %d", dim);
    }

    params->dim = dim;
    params->ring = 2.0;
    params->maxStclSize = 0;
    params->verbose = 1;
    params->degree = 2;
    params->surfType = 0;
    params->epsBeta = 0.0002;
    params->cGlobal = 0.08;
    params->cLocal = 0.5;
    params->kappa1 = 0.3;
    params->kappa0 = 1.0;
    params->markNearDis = true;
    params->wlsInterp0 = false;
    params->wlsUseDag = false;

    nThreads = 1;

#ifdef _OPENMP

    params->nThreads = omp_get_max_threads();

#endif //_OPENMP

    params->parTask = true;
  }

  void rdi_mark_discontinuities(const ::coder::array<double, 2U> &fs, const ::
    coder::array<double, 2U> &alphaCell, const ::coder::array<double, 2U> &beta,
    const ::coder::array<double, 2U> &dfGlobal, const RdiMesh *mesh, const
    RdiParams *params, ::coder::array<signed char, 2U> &disTags)
  {
    double cGlobal;
    double cLocal;
    double hGlobal;
    double kappa0;
    double kappa1;
    int last;

    // rdi_mark_discontinuities - Determine discontinuous nodes
    hGlobal = mesh->hGlobal;
    if (mesh->hGlobal <= 0.0) {
      last = mesh->cellSizes.size(0);
      if (mesh->cellSizes.size(0) <= 2) {
        if (mesh->cellSizes.size(0) == 1) {
          hGlobal = mesh->cellSizes[0];
        } else if (mesh->cellSizes[0] < mesh->cellSizes[mesh->cellSizes.size(0)
                   - 1]) {
          hGlobal = mesh->cellSizes[mesh->cellSizes.size(0) - 1];
        } else {
          hGlobal = mesh->cellSizes[0];
        }
      } else {
        double ex;
        ex = mesh->cellSizes[0];
        for (int k{2}; k <= last; k++) {
          double d;
          d = mesh->cellSizes[k - 1];
          if (ex < d) {
            ex = d;
          }
        }

        hGlobal = ex;
      }
    }

    cGlobal = params->cGlobal;
    if (params->cGlobal < 0.0) {
      cGlobal = 0.05;
    }

    cLocal = params->cLocal;
    if (params->cLocal < 0.0) {
      cLocal = 0.5;
    }

    kappa1 = params->kappa1;
    if (params->kappa1 <= 0.0) {
      kappa1 = 0.3;
    }

    //  C1 dis
    kappa0 = params->kappa0;
    if (params->kappa0 <= 0.0) {
      kappa0 = 1.0;
    }

    //  C0 dis
    if (kappa0 <= kappa1) {
      kappa0 = kappa1;
    }

    if (params->markNearDis) {
      //  cannot run in parallel
      mark_kernel(fs.size(0), hGlobal, cGlobal, cLocal, kappa1, kappa0, fs,
                  alphaCell, beta, dfGlobal, mesh->conn, mesh->cellSizes,
                  mesh->n2cPtr, mesh->n2cList, mesh->n2nPtr, mesh->n2nList,
                  fs.size(1), disTags);
    } else {
      int m2cTryBlkErrCode;
      last = params->nThreads;

      m2cTryBlkErrCode = (0);
      if (params->nThreads <= 0) {
        last = 1;

#ifdef _OPENMP

        last = omp_get_max_threads();

#endif //_OPENMP

      }

#pragma omp parallel default(shared) num_threads(last)

      try { //try
        mark_kernel(fs.size(0), hGlobal, cGlobal, cLocal, kappa1, kappa0, fs,
                    alphaCell, beta, dfGlobal, mesh->conn, mesh->cellSizes,
                    mesh->n2cPtr, mesh->n2cList, fs.size(1), disTags);

      } catch (const std::runtime_error& m2cExc) {
        m2cTryBlkErrCode = (1);
        m2cPrintError("runtime_error %s\n", m2cExc.what());

      }

      catch (const std::logic_error& m2cExc)
      {
        m2cTryBlkErrCode = (1);
        m2cPrintError("logic_error %s\n", m2cExc.what());

      }

      catch (const std::exception& m2cExc)
      {
        m2cTryBlkErrCode = (1);
        m2cPrintError("exception %s\n", m2cExc.what());

      }

      catch (...)
      {
        m2cTryBlkErrCode = (1);
        m2cPrintError("Unknown error detected from C++ exceptions\n");

      } //try

      if ((int)m2cTryBlkErrCode != 0) {
        throw std::runtime_error("omp4m:runtimeErrorInThread");
      }

    }
  }

  static inline
  void rdi_partition(const RdiParams *params, RdiMesh *mesh)
  {
    // rdi_partition - Compute a recursive nodal partition of the mesh
    omp4mRecurPartMesh(mesh->xs.size(0), mesh->conn, params->dim, params->dim,
                       params->nThreads, mesh->parts);
  }
}

// End of code generation (librdi.cpp)
