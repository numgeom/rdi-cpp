// Copyright 2022 The NumGeom Group, Stony Brook University
// Main developers:
//     rdilib: Qiao Chen
//     momp2cpp: Xiangmin Jiao, Qiao Chen
//     wlslib: Xiangmin Jiao, Qiao Chen, Jacob Jones
//
// rdi_compute_stencils.cpp
//
// Code generation for function 'rdi_compute_stencils'
//

// Include files
#include "rdi_compute_stencils.h"
#include "m2c_lib.h"
#include "coder_array.h"
#include "rdi_params.hpp"
#include <algorithm>
#include <cmath>
#include <cstdio>
#ifdef _OPENMP
#include <omp.h>
#endif
#include <stdexcept>
#include <stdio.h>

// Variable Definitions
namespace rdi_stencils
{
  static const signed char iv[12]{ 5, 10, 15, 20, 12, 20, 35, 55, 30, 60, 80,
    100 };

  static const signed char iv1[20]{ 1, 2, 5, 4, 2, 3, 6, 5, 3, 1, 4, 6, 1, 3, 2,
    0, 4, 5, 6, 0 };

  static const signed char iv2[12]{ 2, 4, 3, 1, 3, 4, 4, 2, 1, 3, 1, 2 };

  static const signed char iv3[12]{ 2, 4, 1, 1, 3, 2, 3, 1, 4, 4, 2, 3 };

  static const signed char iv4[12]{ 1, 3, 2, 1, 2, 4, 2, 3, 4, 3, 1, 4 };

  static const signed char iv5[20]{ 1, 4, 3, 2, 1, 2, 5, 0, 2, 3, 5, 0, 3, 4, 5,
    0, 4, 1, 5, 0 };

  static const signed char iv6[24]{ 2, 5, 1, 3, 2, 1, 4, 3, 1, 5, 4, 1, 6, 5, 2,
    6, 2, 3, 6, 3, 4, 6, 4, 5 };

  static const signed char iv7[18]{ 1, 3, 4, 2, 1, 4, 3, 2, 4, 3, 1, 5, 1, 2, 5,
    2, 3, 5 };

  static const signed char iv8[20]{ 2, 5, 1, 0, 3, 2, 1, 0, 4, 3, 1, 0, 5, 4, 1,
    0, 2, 3, 4, 5 };

  static const signed char iv9[12]{ 1, 2, 4, 1, 2, 3, 1, 3, 4, 2, 3, 4 };
}

// Function Declarations
namespace rdi_stencils
{
  static inline
  void append_one_ring(int vid, const ::coder::array<int, 2U> &tets,
    const ::coder::array<int, 2U> &sibhfs, const ::coder::array<int, 1U> &v2hf,
    int ngbvs[2048], ::coder::array<boolean_T, 1U> &vtags, ::coder::array<
    boolean_T, 1U> &etags, int ngbes[2048], int *nverts, int *nelems);
  static inline
  int b_append_one_ring(int vid, const ::coder::array<int, 2U> &tets,
    const ::coder::array<int, 2U> &sibhfs, const ::coder::array<int, 1U> &v2hf,
    int ngbvs[2048], int *nverts, ::coder::array<boolean_T, 1U> &vtags, ::coder::
    array<boolean_T, 1U> &etags, int ngbes[128]);
  static inline
  void compute_stcl_kernel1(int n, const ::coder::array<int, 2U> &conn,
    double ring, ::coder::array<int, 2U> &stcls, const ::coder::array<int, 1U>
    &nrange);
  static inline
  void compute_stcl_kernel1(int n, const ::coder::array<int, 2U> &conn,
    double ring, ::coder::array<int, 2U> &stcls);
  static inline
  void compute_stcl_kernel2(int n, const ::coder::array<int, 2U> &conn,
    double ring, ::coder::array<int, 2U> &stcls, const ::coder::array<int, 1U>
    &nrange);
  static inline
  void compute_stcl_kernel2(int n, const ::coder::array<int, 2U> &conn,
    double ring, ::coder::array<int, 2U> &stcls);
  static inline
  void compute_stcl_kernel_tet(int n, const ::coder::array<int, 2U> &conn,
    double ring, ::coder::array<int, 2U> &stcls, const ::coder::array<int, 1U>
    &nrange);
  static inline
  void compute_stcl_kernel_tet(int n, const ::coder::array<int, 2U> &conn,
    double ring, ::coder::array<int, 2U> &stcls);
  static inline
  void determine_border_vertices_vol(int nv, const ::coder::array<int, 2U>
    &elems, ::coder::array<int, 2U> &sibhfs, ::coder::array<boolean_T, 1U>
    &isborder);
  static inline
  void determine_incident_halfedges(int nv, const ::coder::array<int, 2U>
    &elems, const ::coder::array<int, 2U> &sibhes, ::coder::array<int, 1U> &v2he);
  static inline
  void determine_incident_halffaces(int nv, const ::coder::array<int, 2U>
    &elems, const ::coder::array<int, 2U> &sibhfs, ::coder::array<int, 1U> &v2hf);
  static inline
  void determine_sibling_halfedges(int nv, const ::coder::array<int, 2U>
    &elems, ::coder::array<int, 2U> &sibhes);
  static inline
  void determine_sibling_halffaces(int nv, const ::coder::array<int, 2U>
    &elems, ::coder::array<int, 2U> &sibhfs);
  static inline
  void determine_sibling_halffaces_mixed(int nv, const ::coder::array<int,
    2U> &elems, ::coder::array<int, 1U> &sibhfs);
  static inline
  void determine_sibling_halffaces_pyramid(int nv, const ::coder::array<
    int, 2U> &elems, ::coder::array<int, 2U> &sibhfs);
  static inline
  void determine_sibling_halffaces_tet(int nv, const ::coder::array<int,
    2U> &elems, ::coder::array<int, 2U> &sibhfs, boolean_T *manifold, boolean_T *
    oriented);
  static inline
  void determine_sibling_halfverts(int nv, const ::coder::array<int, 2U>
    &edges, ::coder::array<int, 2U> &sibhvs);
  static inline
  int elem_one_ring(int vid, const ::coder::array<int, 2U> &tets, const ::
    coder::array<int, 2U> &sibhfs, const ::coder::array<int, 1U> &v2hf, ::coder::
    array<boolean_T, 1U> &etags, int ngbes[128]);
  static inline
  void obtain_nring_curv(int vid, double ring, const ::coder::array<int,
    2U> &edgs, const ::coder::array<int, 2U> &sibhvs, const ::coder::array<int,
    1U> &v2hv, int ngbvs[128], ::coder::array<boolean_T, 1U> &vtags, ::coder::
    array<boolean_T, 1U> &etags, int ngbes[128], int *nverts, int *nedges);
  static inline
  void obtain_nring_quad(int vid, double ring, const ::coder::array<int,
    2U> &elems, const ::coder::array<int, 2U> &sibhes, const ::coder::array<int,
    1U> &v2he, int ngbvs[1024], ::coder::array<boolean_T, 1U> &vtags, ::coder::
    array<boolean_T, 1U> &ftags, int ngbfs[1024], int *nverts, int *nfaces);
  static inline
  void obtain_nring_vol(int vid, double ring, const ::coder::array<int,
    2U> &tets, const ::coder::array<int, 2U> &sibhfs, const ::coder::array<int,
    1U> &v2hf, int ngbvs[2048], int ngbes[2048], ::coder::array<boolean_T, 1U>
    &vtags, ::coder::array<boolean_T, 1U> &etags, ::coder::array<boolean_T, 1U>
    &etags_elem, int *nverts, int *nelems);
}

// Function Definitions
namespace rdi_stencils
{
  static void append_one_ring(int vid, const ::coder::array<int, 2U> &tets,
    const ::coder::array<int, 2U> &sibhfs, const ::coder::array<int, 1U> &v2hf,
    int ngbvs[2048], ::coder::array<boolean_T, 1U> &vtags, ::coder::array<
    boolean_T, 1U> &etags, int ngbes[2048], int *nverts, int *nelems)
  {
    int stack[1024];
    *nelems = 0;
    *nverts = 0;

    //  This determines the 1-ring of vertices aroung the element VID
    if (v2hf[vid - 1] != 0) {
      coder::SizeType size_stack;
      boolean_T exitg1;
      boolean_T overflow;

      stack[0] = static_cast<int>(static_cast<coder::SizeType>(v2hf[vid - 1]) >> 3);

      //  Element (region) ID
      overflow = false;

      //  Create a stack for storing tets
      size_stack = 1;

      //  sibhfs_tet(lvid, :) gives the faces that border on local vertex lvid.
      exitg1 = false;
      while ((!exitg1) && (size_stack > 0)) {
        coder::SizeType rid;

        //  Pop the element from top of stack
        rid = stack[size_stack - 1] - 1;
        size_stack--;

        //  Append element
        if (*nelems >= 1024) {
          m2cPrintf("Overflow in elements in append_one_ring.m \n");
          fflush(stdout);
          exitg1 = true;
        } else {
          if (!etags[rid]) {
            coder::SizeType ii;
            coder::SizeType lvid;
            boolean_T exitg2;
            etags[rid] = true;
            (*nelems)++;
            ngbes[*nelems - 1] = rid + 1;
            lvid = 0;

            //  Stores the local vertex of vid
            ii = 0;
            exitg2 = false;
            while ((!exitg2) && (ii < 4)) {
              coder::SizeType v;
              v = tets[ii + tets.size(1) * rid];
              if (v == vid) {
                lvid = ii + 1;
              }

              if (!vtags[v - 1]) {
                if (*nverts < 1024) {
                  vtags[v - 1] = true;
                  (*nverts)++;
                  ngbvs[*nverts - 1] = v;
                  ii++;
                } else {
                  overflow = true;
                  exitg2 = true;
                }
              } else {
                ii++;
              }
            }

            //  Push unvisited neighbor tets onto stack
            ii = 0;
            while ((ii < 3) && (lvid != 0)) {
              unsigned int c;

              c = static_cast<coder::SizeType>(sibhfs[(iv9[ii + 3 * (lvid - 1)] +
                sibhfs.size(1) * rid) - 1]) >> 3;
              if ((static_cast<int>(c) != 0) && (!etags[static_cast<coder::SizeType>(c) - 1]))
              {
                if (size_stack >= 1024) {
                  m2cPrintf("Overflow in stack in append_one_ring.m \n");
                  fflush(stdout);
                  overflow = true;
                } else {
                  size_stack++;
                  stack[size_stack - 1] = static_cast<int>(c);
                }
              }

              ii++;
            }
          }

          if (overflow) {
            exitg1 = true;
          }
        }
      }
    }
  }

  static int b_append_one_ring(int vid, const ::coder::array<int, 2U> &tets,
    const ::coder::array<int, 2U> &sibhfs, const ::coder::array<int, 1U> &v2hf,
    int ngbvs[2048], int *nverts, ::coder::array<boolean_T, 1U> &vtags, ::coder::
    array<boolean_T, 1U> &etags, int ngbes[128])
  {
    int stack[128];
    coder::SizeType nelems;
    nelems = 0;

    //  This determines the 1-ring of vertices aroung the element VID
    if (v2hf[vid - 1] != 0) {
      coder::SizeType size_stack;
      boolean_T exitg1;
      boolean_T overflow;

      stack[0] = static_cast<int>(static_cast<coder::SizeType>(v2hf[vid - 1]) >> 3);

      //  Element (region) ID
      overflow = false;

      //  Create a stack for storing tets
      size_stack = 1;

      //  sibhfs_tet(lvid, :) gives the faces that border on local vertex lvid.
      exitg1 = false;
      while ((!exitg1) && (size_stack > 0)) {
        coder::SizeType rid;

        //  Pop the element from top of stack
        rid = stack[size_stack - 1] - 1;
        size_stack--;

        //  Append element
        if (nelems >= 128) {
          m2cPrintf("Overflow in elements in append_one_ring.m \n");
          fflush(stdout);
          exitg1 = true;
        } else {
          if (!etags[rid]) {
            coder::SizeType ii;
            coder::SizeType lvid;
            boolean_T exitg2;
            etags[rid] = true;
            nelems++;
            ngbes[nelems - 1] = rid + 1;
            lvid = 0;

            //  Stores the local vertex of vid
            ii = 0;
            exitg2 = false;
            while ((!exitg2) && (ii < 4)) {
              coder::SizeType v;
              v = tets[ii + tets.size(1) * rid];
              if (v == vid) {
                lvid = ii + 1;
              }

              if (!vtags[v - 1]) {
                if (*nverts < 1024) {
                  vtags[v - 1] = true;
                  (*nverts)++;
                  ngbvs[*nverts - 1] = v;
                  ii++;
                } else {
                  overflow = true;
                  exitg2 = true;
                }
              } else {
                ii++;
              }
            }

            //  Push unvisited neighbor tets onto stack
            ii = 0;
            while ((ii < 3) && (lvid != 0)) {
              unsigned int c;

              c = static_cast<coder::SizeType>(sibhfs[(iv9[ii + 3 * (lvid - 1)] +
                sibhfs.size(1) * rid) - 1]) >> 3;
              if ((static_cast<int>(c) != 0) && (!etags[static_cast<coder::SizeType>(c) - 1]))
              {
                if (size_stack >= 128) {
                  m2cPrintf("Overflow in stack in append_one_ring.m \n");
                  fflush(stdout);
                  overflow = true;
                } else {
                  size_stack++;
                  stack[size_stack - 1] = static_cast<int>(c);
                }
              }

              ii++;
            }
          }

          if (overflow) {
            exitg1 = true;
          }
        }
      }
    }

    return nelems;
  }

  static void compute_stcl_kernel1(int n, const ::coder::array<int, 2U> &conn,
    double ring, ::coder::array<int, 2U> &stcls, const ::coder::array<int, 1U>
    &nrange)
  {
    ::coder::array<int, 2U> sibhvs_;
    ::coder::array<int, 1U> v2hv_;
    ::coder::array<boolean_T, 1U> etags_;
    ::coder::array<boolean_T, 1U> vtags_;
    int ngbes[128];
    int ngbvs[128];
    coder::SizeType N;
    int a__1;
    coder::SizeType loop_ub;
    coder::SizeType maxStclSize;
    int nVerts;
    coder::SizeType nedgs;

    //  kernel for 1D
    determine_sibling_halfverts(n, conn, sibhvs_);

    //  DETERMINE_INCIDENT_HALFVERTS Determine an incident half-vertex.
    nedgs = conn.size(0);

    //  Construct a vertex to halfedge mapping.
    v2hv_.set_size(n);
    for (coder::SizeType i{0}; i < n; i++) {
      v2hv_[i] = 0;
    }

    //  Compute v2hv.
    for (coder::SizeType kk{0}; kk < nedgs; kk++) {
      coder::SizeType v;
      v = conn[conn.size(1) * kk];
      if ((v > 0) && (v2hv_[v - 1] == 0)) {
        //  Encode <eid,lvid> pair into a hvid.
        v2hv_[v - 1] = (kk + 1) << 1;
      }

      v = conn[conn.size(1) * kk + 1];
      if ((v > 0) && (v2hv_[v - 1] == 0)) {
        //  Encode <eid,lvid> pair into a hvid.
        v2hv_[v - 1] = ((kk + 1) << 1) + 1;
      }
    }

    //  buffers
    vtags_.set_size(n);
    for (coder::SizeType i{0}; i < n; i++) {
      vtags_[i] = false;
    }

    etags_.set_size(conn.size(0));
    loop_ub = conn.size(0);
    for (coder::SizeType i{0}; i < loop_ub; i++) {
      etags_[i] = false;
    }

    maxStclSize = stcls.size(1) - 1;
    N = n;
    if (nrange.size(0) != 0) {
      N = nrange.size(0);
    }

    for (coder::SizeType b_i{0}; b_i < N; b_i++) {
      coder::SizeType n0;
      coder::SizeType nid;
      if (nrange.size(0) != 0) {
        nid = nrange[b_i];
      } else {
        nid = b_i + 1;
      }

      obtain_nring_curv(nid, ring, conn, sibhvs_, v2hv_, ngbvs, vtags_, etags_,
                        ngbes, &nVerts, &a__1);

      //  NOTE: ngbvs doesn't count i
      n0 = nVerts + 1;
      if (nVerts + 1 > maxStclSize) {
        n0 = maxStclSize;
      }

      stcls[stcls.size(1) * b_i] = nid;
      stcls[(stcls.size(1) + stcls.size(1) * b_i) - 1] = n0;
      for (coder::SizeType j{2}; j <= n0; j++) {
        stcls[(j + stcls.size(1) * b_i) - 1] = ngbvs[j - 2];
      }
    }
  }

  static void compute_stcl_kernel1(int n, const ::coder::array<int, 2U> &conn,
    double ring, ::coder::array<int, 2U> &stcls)
  {
    ::coder::array<int, 2U> sibhvs_;
    ::coder::array<int, 1U> v2hv_;
    ::coder::array<boolean_T, 1U> etags_;
    ::coder::array<boolean_T, 1U> vtags_;
    int ngbes[128];
    int ngbvs[128];
    int a__1;
    coder::SizeType loop_ub;
    coder::SizeType maxStclSize;
    int nVerts;
    coder::SizeType nedgs;

    //  kernel for 1D
    determine_sibling_halfverts(n, conn, sibhvs_);

    //  DETERMINE_INCIDENT_HALFVERTS Determine an incident half-vertex.
    nedgs = conn.size(0);

    //  Construct a vertex to halfedge mapping.
    v2hv_.set_size(n);
    for (coder::SizeType i{0}; i < n; i++) {
      v2hv_[i] = 0;
    }

    //  Compute v2hv.
    for (coder::SizeType kk{0}; kk < nedgs; kk++) {
      coder::SizeType v;
      v = conn[conn.size(1) * kk];
      if ((v > 0) && (v2hv_[v - 1] == 0)) {
        //  Encode <eid,lvid> pair into a hvid.
        v2hv_[v - 1] = (kk + 1) << 1;
      }

      v = conn[conn.size(1) * kk + 1];
      if ((v > 0) && (v2hv_[v - 1] == 0)) {
        //  Encode <eid,lvid> pair into a hvid.
        v2hv_[v - 1] = ((kk + 1) << 1) + 1;
      }
    }

    //  buffers
    vtags_.set_size(n);
    for (coder::SizeType i{0}; i < n; i++) {
      vtags_[i] = false;
    }

    etags_.set_size(conn.size(0));
    loop_ub = conn.size(0);
    for (coder::SizeType i{0}; i < loop_ub; i++) {
      etags_[i] = false;
    }

    maxStclSize = stcls.size(1) - 1;
    for (coder::SizeType b_i{0}; b_i < n; b_i++) {
      coder::SizeType n0;
      obtain_nring_curv(b_i + 1, ring, conn, sibhvs_, v2hv_, ngbvs, vtags_,
                        etags_, ngbes, &nVerts, &a__1);

      //  NOTE: ngbvs doesn't count i
      n0 = nVerts + 1;
      if (nVerts + 1 > maxStclSize) {
        n0 = maxStclSize;
      }

      stcls[stcls.size(1) * b_i] = b_i + 1;
      stcls[(stcls.size(1) + stcls.size(1) * b_i) - 1] = n0;
      for (coder::SizeType j{2}; j <= n0; j++) {
        stcls[(j + stcls.size(1) * b_i) - 1] = ngbvs[j - 2];
      }
    }
  }

  static void compute_stcl_kernel2(int n, const ::coder::array<int, 2U> &conn,
    double ring, ::coder::array<int, 2U> &stcls, const ::coder::array<int, 1U>
    &nrange)
  {
    ::coder::array<int, 2U> opphes;
    ::coder::array<int, 1U> v2he;
    ::coder::array<boolean_T, 1U> etags_;
    ::coder::array<boolean_T, 1U> vtags_;
    int ngbfs[1024];
    int ngbvs[1024];
    coder::SizeType N;
    int a__1;
    coder::SizeType loop_ub;
    coder::SizeType maxStclSize;
    int numNbs;

    //  kernel for 2D
    determine_sibling_halfedges(n, conn, opphes);
    determine_incident_halfedges(n, conn, opphes, v2he);
    N = n;
    if (nrange.size(0) != 0) {
      N = nrange.size(0);
    }

    maxStclSize = stcls.size(1) - 1;

    //  buffers
    vtags_.set_size(n);
    for (coder::SizeType i{0}; i < n; i++) {
      vtags_[i] = false;
    }

    etags_.set_size(conn.size(0));
    loop_ub = conn.size(0);
    for (coder::SizeType i{0}; i < loop_ub; i++) {
      etags_[i] = false;
    }

    for (coder::SizeType b_i{0}; b_i < N; b_i++) {
      coder::SizeType n0;
      coder::SizeType nid;
      if (nrange.size(0) != 0) {
        nid = nrange[b_i];
      } else {
        nid = b_i + 1;
      }

      obtain_nring_quad(nid, ring, conn, opphes, v2he, ngbvs, vtags_, etags_,
                        ngbfs, &numNbs, &a__1);

      //  NOTE: num_nbs doesn't count nid
      n0 = numNbs + 1;
      if (numNbs + 1 > maxStclSize) {
        n0 = maxStclSize;
      }

      stcls[stcls.size(1) * b_i] = nid;
      stcls[(stcls.size(1) + stcls.size(1) * b_i) - 1] = n0;
      for (coder::SizeType j{2}; j <= n0; j++) {
        stcls[(j + stcls.size(1) * b_i) - 1] = ngbvs[j - 2];
      }
    }
  }

  static void compute_stcl_kernel2(int n, const ::coder::array<int, 2U> &conn,
    double ring, ::coder::array<int, 2U> &stcls)
  {
    ::coder::array<int, 2U> opphes;
    ::coder::array<int, 1U> v2he;
    ::coder::array<boolean_T, 1U> etags_;
    ::coder::array<boolean_T, 1U> vtags_;
    int ngbfs[1024];
    int ngbvs[1024];
    int a__1;
    coder::SizeType loop_ub;
    coder::SizeType maxStclSize;
    int numNbs;

    //  kernel for 2D
    determine_sibling_halfedges(n, conn, opphes);
    determine_incident_halfedges(n, conn, opphes, v2he);
    maxStclSize = stcls.size(1) - 1;

    //  buffers
    vtags_.set_size(n);
    for (coder::SizeType i{0}; i < n; i++) {
      vtags_[i] = false;
    }

    etags_.set_size(conn.size(0));
    loop_ub = conn.size(0);
    for (coder::SizeType i{0}; i < loop_ub; i++) {
      etags_[i] = false;
    }

    for (coder::SizeType b_i{0}; b_i < n; b_i++) {
      coder::SizeType n0;
      obtain_nring_quad(b_i + 1, ring, conn, opphes, v2he, ngbvs, vtags_, etags_,
                        ngbfs, &numNbs, &a__1);

      //  NOTE: num_nbs doesn't count nid
      n0 = numNbs + 1;
      if (numNbs + 1 > maxStclSize) {
        n0 = maxStclSize;
      }

      stcls[stcls.size(1) * b_i] = b_i + 1;
      stcls[(stcls.size(1) + stcls.size(1) * b_i) - 1] = n0;
      for (coder::SizeType j{2}; j <= n0; j++) {
        stcls[(j + stcls.size(1) * b_i) - 1] = ngbvs[j - 2];
      }
    }
  }

  static void compute_stcl_kernel_tet(int n, const ::coder::array<int, 2U> &conn,
    double ring, ::coder::array<int, 2U> &stcls, const ::coder::array<int, 1U>
    &nrange)
  {
    ::coder::array<int, 2U> sibhfs_;
    ::coder::array<int, 1U> v2hf_;
    ::coder::array<boolean_T, 1U> etagsElem_;
    ::coder::array<boolean_T, 1U> etags_;
    ::coder::array<boolean_T, 1U> vtags_;
    int ngbes[2048];
    int ngbvs[2048];
    coder::SizeType N;
    int a__1;
    coder::SizeType loop_ub;
    coder::SizeType maxStclSize;
    int nVerts;

    //  kernel for 3D
    determine_sibling_halffaces(n, conn, sibhfs_);
    determine_incident_halffaces(n, conn, sibhfs_, v2hf_);

    //  buffer for AHF
    vtags_.set_size(n);
    for (coder::SizeType i{0}; i < n; i++) {
      vtags_[i] = false;
    }

    etags_.set_size(conn.size(0));
    loop_ub = conn.size(0);
    for (coder::SizeType i{0}; i < loop_ub; i++) {
      etags_[i] = false;
    }

    etagsElem_.set_size(conn.size(0));
    loop_ub = conn.size(0);
    for (coder::SizeType i{0}; i < loop_ub; i++) {
      etagsElem_[i] = false;
    }

    N = n;
    if (nrange.size(0) != 0) {
      N = nrange.size(0);
    }

    maxStclSize = stcls.size(1) - 1;
    for (coder::SizeType b_i{0}; b_i < N; b_i++) {
      coder::SizeType n0;
      coder::SizeType nid;
      if (nrange.size(0) != 0) {
        nid = nrange[b_i];
      } else {
        nid = b_i + 1;
      }

      obtain_nring_vol(nid, ring, conn, sibhfs_, v2hf_, ngbvs, ngbes, vtags_,
                       etags_, etagsElem_, &nVerts, &a__1);

      //  NOTE: ngbvs doesn't count i
      n0 = nVerts + 1;
      if (nVerts + 1 > maxStclSize) {
        n0 = maxStclSize;
      }

      stcls[stcls.size(1) * b_i] = nid;
      stcls[(stcls.size(1) + stcls.size(1) * b_i) - 1] = n0;
      for (coder::SizeType j{2}; j <= n0; j++) {
        stcls[(j + stcls.size(1) * b_i) - 1] = ngbvs[j - 2];
      }
    }
  }

  static void compute_stcl_kernel_tet(int n, const ::coder::array<int, 2U> &conn,
    double ring, ::coder::array<int, 2U> &stcls)
  {
    ::coder::array<int, 2U> sibhfs_;
    ::coder::array<int, 1U> v2hf_;
    ::coder::array<boolean_T, 1U> etagsElem_;
    ::coder::array<boolean_T, 1U> etags_;
    ::coder::array<boolean_T, 1U> vtags_;
    int ngbes[2048];
    int ngbvs[2048];
    int a__1;
    coder::SizeType loop_ub;
    coder::SizeType maxStclSize;
    int nVerts;

    //  kernel for 3D
    determine_sibling_halffaces(n, conn, sibhfs_);
    determine_incident_halffaces(n, conn, sibhfs_, v2hf_);

    //  buffer for AHF
    vtags_.set_size(n);
    for (coder::SizeType i{0}; i < n; i++) {
      vtags_[i] = false;
    }

    etags_.set_size(conn.size(0));
    loop_ub = conn.size(0);
    for (coder::SizeType i{0}; i < loop_ub; i++) {
      etags_[i] = false;
    }

    etagsElem_.set_size(conn.size(0));
    loop_ub = conn.size(0);
    for (coder::SizeType i{0}; i < loop_ub; i++) {
      etagsElem_[i] = false;
    }

    maxStclSize = stcls.size(1) - 1;
    for (coder::SizeType b_i{0}; b_i < n; b_i++) {
      coder::SizeType n0;
      obtain_nring_vol(b_i + 1, ring, conn, sibhfs_, v2hf_, ngbvs, ngbes, vtags_,
                       etags_, etagsElem_, &nVerts, &a__1);

      //  NOTE: ngbvs doesn't count i
      n0 = nVerts + 1;
      if (nVerts + 1 > maxStclSize) {
        n0 = maxStclSize;
      }

      stcls[stcls.size(1) * b_i] = b_i + 1;
      stcls[(stcls.size(1) + stcls.size(1) * b_i) - 1] = n0;
      for (coder::SizeType j{2}; j <= n0; j++) {
        stcls[(j + stcls.size(1) * b_i) - 1] = ngbvs[j - 2];
      }
    }
  }

  static void determine_border_vertices_vol(int nv, const ::coder::array<int, 2U>
    &elems, ::coder::array<int, 2U> &sibhfs, ::coder::array<boolean_T, 1U>
    &isborder)
  {
    static const signed char hf_hex[54]{ 1, 4, 3, 2, 12, 11, 10, 9, 21, 1, 2, 6,
      5, 9, 14, 17, 13, 22, 2, 3, 7, 6, 10, 15, 18, 14, 23, 3, 4, 8, 7, 11, 16,
      19, 15, 24, 4, 1, 5, 8, 13, 20, 16, 12, 25, 5, 6, 7, 8, 17, 18, 19, 20, 26
    };

    static const signed char hf_pri[45]{ 1, 2, 5, 4, 7, 11, 13, 10, 16, 2, 3, 6,
      5, 8, 12, 14, 11, 17, 3, 1, 4, 6, 9, 10, 15, 12, 18, 1, 3, 2, 9, 8, 7, 0,
      0, 0, 4, 5, 6, 13, 14, 15, 0, 0, 0 };

    static const signed char hf_pyr[45]{ 1, 4, 3, 2, 9, 8, 7, 6, 14, 1, 2, 5, 6,
      11, 10, 0, 0, 0, 2, 3, 5, 7, 12, 11, 0, 0, 0, 3, 4, 5, 8, 13, 12, 0, 0, 0,
      4, 1, 5, 9, 10, 13, 0, 0, 0 };

    static const signed char hf_tet[24]{ 1, 3, 2, 7, 5, 6, 1, 2, 4, 5, 9, 8, 2,
      3, 4, 6, 10, 9, 3, 1, 4, 7, 8, 10 };

    coder::SizeType hf_size_idx_0;
    coder::SizeType hf_size_idx_1;
    coder::SizeType i;
    coder::SizeType nvpf;
    signed char hf_data[54];

    //  DETERMINE_BORDER_VERTICES_VOL Determine border vertices of a volume mesh.
    if ((sibhfs.size(0) == 0) || (sibhfs.size(1) == 0)) {
      determine_sibling_halffaces(nv, elems, sibhfs);
    }

    isborder.set_size(nv);
    for (i = 0; i < nv; i++) {
      isborder[i] = false;
    }

    //  List vertices in counter-clockwise order, so that faces are outwards.
    if (elems.size(1) == 1) {
      coder::SizeType offset;
      coder::SizeType offset_o;

      //  Mixed elements
      offset = 0;
      offset_o = 1;
      while (offset + 1 < elems.size(0)) {
        coder::SizeType i1;
        boolean_T b;
        boolean_T b1;
        boolean_T guard1;
        boolean_T guard2;
        boolean_T guard3;
        boolean_T guard4;
        guard1 = false;
        guard2 = false;
        guard3 = false;
        guard4 = false;
        switch (elems[offset % elems.size(0) + offset / elems.size(0)]) {
         case 4:
          guard1 = true;
          break;

         case 10:
          guard1 = true;
          break;

         case 5:
          guard2 = true;
          break;

         case 14:
          guard2 = true;
          break;

         case 6:
          guard3 = true;
          break;

         case 15:
          guard3 = true;
          break;

         case 18:
          guard3 = true;
          break;

         case 8:
          guard4 = true;
          break;

         case 20:
          guard4 = true;
          break;

         case 27:
          guard4 = true;
          break;

         default:
          m2cErrMsgIdAndTxt("momp2cpp:runtimeError",
                            "Unrecognized element type.");
          break;
        }

        if (guard4) {
          //  Hexahedral
          nvpf = (elems[offset % elems.size(0) + offset / elems.size(0)] == 27)
            + 3;
          b = true;
          b1 = ((sibhfs.size(1) <= 0) || (sibhfs.size(0) <= 0));
          i = sibhfs.size(1) * sibhfs.size(0);
          hf_size_idx_0 = 0;
          for (coder::SizeType jj{0}; jj < 6; jj++) {
            i1 = offset_o + jj;
            if (b1 || (i1 < 0) || (i1 >= i)) {
              hf_size_idx_0 = 0;
              b = true;
            } else if (b) {
              b = false;
              hf_size_idx_0 = i1 % sibhfs.size(0) * sibhfs.size(1) + i1 /
                sibhfs.size(0);
            } else {
              hf_size_idx_1 = sibhfs.size(1) * sibhfs.size(0) - 1;
              if (hf_size_idx_0 > MAX_int32_T - sibhfs.size(1)) {
                hf_size_idx_0 = i1 % sibhfs.size(0) * sibhfs.size(1) + i1 /
                  sibhfs.size(0);
              } else {
                hf_size_idx_0 += sibhfs.size(1);
                if (hf_size_idx_0 > hf_size_idx_1) {
                  hf_size_idx_0 -= hf_size_idx_1;
                }
              }
            }

            if (sibhfs[hf_size_idx_0] == 0) {
              for (coder::SizeType kk{0}; kk <= nvpf; kk++) {
                i1 = offset + hf_hex[kk + 9 * jj];
                isborder[elems[i1 % elems.size(0) + i1 / elems.size(0)] - 1] =
                  true;
              }
            }
          }
        }

        if (guard3) {
          //  Prism
          b = true;
          b1 = ((sibhfs.size(1) <= 0) || (sibhfs.size(0) <= 0));
          i = sibhfs.size(1) * sibhfs.size(0);
          hf_size_idx_0 = 0;
          for (coder::SizeType jj{0}; jj < 5; jj++) {
            i1 = offset_o + jj;
            if (b1 || (i1 < 0) || (i1 >= i)) {
              hf_size_idx_0 = 0;
              b = true;
            } else if (b) {
              b = false;
              hf_size_idx_0 = i1 % sibhfs.size(0) * sibhfs.size(1) + i1 /
                sibhfs.size(0);
            } else {
              hf_size_idx_1 = sibhfs.size(1) * sibhfs.size(0) - 1;
              if (hf_size_idx_0 > MAX_int32_T - sibhfs.size(1)) {
                hf_size_idx_0 = i1 % sibhfs.size(0) * sibhfs.size(1) + i1 /
                  sibhfs.size(0);
              } else {
                hf_size_idx_0 += sibhfs.size(1);
                if (hf_size_idx_0 > hf_size_idx_1) {
                  hf_size_idx_0 -= hf_size_idx_1;
                }
              }
            }

            if (sibhfs[hf_size_idx_0] == 0) {
              boolean_T b2;
              if ((elems[offset % elems.size(0) + offset / elems.size(0)] == 18)
                  && (jj + 1 < 4)) {
                b2 = true;
              } else {
                b2 = false;
              }

              nvpf = ((jj + 1 < 4) + b2) + 2;
              for (coder::SizeType kk{0}; kk <= nvpf; kk++) {
                i1 = offset + hf_pri[kk + 9 * jj];
                isborder[elems[i1 % elems.size(0) + i1 / elems.size(0)] - 1] =
                  true;
              }
            }
          }
        }

        if (guard2) {
          //  Pyramid
          b = true;
          b1 = ((sibhfs.size(1) <= 0) || (sibhfs.size(0) <= 0));
          i = sibhfs.size(1) * sibhfs.size(0);
          hf_size_idx_0 = 0;
          for (coder::SizeType jj{0}; jj < 5; jj++) {
            i1 = offset_o + jj;
            if (b1 || (i1 < 0) || (i1 >= i)) {
              hf_size_idx_0 = 0;
              b = true;
            } else if (b) {
              b = false;
              hf_size_idx_0 = i1 % sibhfs.size(0) * sibhfs.size(1) + i1 /
                sibhfs.size(0);
            } else {
              hf_size_idx_1 = sibhfs.size(1) * sibhfs.size(0) - 1;
              if (hf_size_idx_0 > MAX_int32_T - sibhfs.size(1)) {
                hf_size_idx_0 = i1 % sibhfs.size(0) * sibhfs.size(1) + i1 /
                  sibhfs.size(0);
              } else {
                hf_size_idx_0 += sibhfs.size(1);
                if (hf_size_idx_0 > hf_size_idx_1) {
                  hf_size_idx_0 -= hf_size_idx_1;
                }
              }
            }

            if (sibhfs[hf_size_idx_0] == 0) {
              nvpf = (jj + 1 == 1) + 2;
              for (coder::SizeType kk{0}; kk <= nvpf; kk++) {
                i1 = offset + hf_pyr[kk + 9 * jj];
                isborder[elems[i1 % elems.size(0) + i1 / elems.size(0)] - 1] =
                  true;
              }
            }
          }
        }

        if (guard1) {
          //  Tetrahedral
          b = true;
          b1 = ((sibhfs.size(1) <= 0) || (sibhfs.size(0) <= 0));
          i = sibhfs.size(1) * sibhfs.size(0);
          hf_size_idx_0 = 0;
          for (coder::SizeType jj{0}; jj < 4; jj++) {
            i1 = offset_o + jj;
            if (b1 || (i1 < 0) || (i1 >= i)) {
              hf_size_idx_0 = 0;
              b = true;
            } else if (b) {
              b = false;
              hf_size_idx_0 = i1 % sibhfs.size(0) * sibhfs.size(1) + i1 /
                sibhfs.size(0);
            } else {
              hf_size_idx_1 = sibhfs.size(1) * sibhfs.size(0) - 1;
              if (hf_size_idx_0 > MAX_int32_T - sibhfs.size(1)) {
                hf_size_idx_0 = i1 % sibhfs.size(0) * sibhfs.size(1) + i1 /
                  sibhfs.size(0);
              } else {
                hf_size_idx_0 += sibhfs.size(1);
                if (hf_size_idx_0 > hf_size_idx_1) {
                  hf_size_idx_0 -= hf_size_idx_1;
                }
              }
            }

            if (sibhfs[hf_size_idx_0] == 0) {
              i1 = offset + hf_tet[6 * jj];
              isborder[elems[i1 % elems.size(0) + i1 / elems.size(0)] - 1] =
                true;
              i1 = offset + hf_tet[6 * jj + 1];
              isborder[elems[i1 % elems.size(0) + i1 / elems.size(0)] - 1] =
                true;
              i1 = offset + hf_tet[6 * jj + 2];
              isborder[elems[i1 % elems.size(0) + i1 / elems.size(0)] - 1] =
                true;
            }
          }
        }

        offset = (offset + elems[offset % elems.size(0) + offset / elems.size(0)])
          + 1;
        offset_o = (offset_o + sibhfs[(offset_o - 1) % sibhfs.size(0) *
                    sibhfs.size(1) + (offset_o - 1) / sibhfs.size(0)]) + 1;
      }
    } else {
      coder::SizeType ii;
      boolean_T guard1{ false };

      boolean_T guard2{ false };

      boolean_T guard3{ false };

      boolean_T guard4{ false };

      //  Table for local IDs of incident faces of each vertex.
      guard1 = false;
      guard2 = false;
      guard3 = false;
      guard4 = false;
      switch (elems.size(1)) {
       case 4:
        guard1 = true;
        break;

       case 10:
        guard1 = true;
        break;

       case 5:
        guard2 = true;
        break;

       case 14:
        guard2 = true;
        break;

       case 6:
        guard3 = true;
        break;

       case 15:
        guard3 = true;
        break;

       case 18:
        guard3 = true;
        break;

       case 8:
        guard4 = true;
        break;

       case 20:
        guard4 = true;
        break;

       case 27:
        guard4 = true;
        break;

       default:
        hf_size_idx_1 = 6;
        hf_size_idx_0 = 4;
        for (i = 0; i < 24; i++) {
          hf_data[i] = hf_tet[i];
        }

        nvpf = -4;

        m2cErrMsgIdAndTxt("momp2cpp:runtimeError", "Unsupported element type");
        break;
      }

      if (guard4) {
        hf_size_idx_1 = 9;
        hf_size_idx_0 = 6;
        for (i = 0; i < 54; i++) {
          hf_data[i] = hf_hex[i];
        }

        nvpf = (elems.size(1) == 27);
      }

      if (guard3) {
        hf_size_idx_1 = 9;
        hf_size_idx_0 = 5;
        for (i = 0; i < 45; i++) {
          hf_data[i] = hf_pri[i];
        }

        nvpf = (elems.size(1) == 18);
      }

      if (guard2) {
        hf_size_idx_1 = 9;
        hf_size_idx_0 = 5;
        for (i = 0; i < 45; i++) {
          hf_data[i] = hf_pyr[i];
        }

        nvpf = 1;
      }

      if (guard1) {
        hf_size_idx_1 = 6;
        hf_size_idx_0 = 4;
        for (i = 0; i < 24; i++) {
          hf_data[i] = hf_tet[i];
        }

        nvpf = -1;
      }

      ii = 0;
      while ((ii + 1 <= elems.size(0)) && (elems[elems.size(1) * ii] != 0)) {
        for (coder::SizeType jj{0}; jj < hf_size_idx_0; jj++) {
          if (sibhfs[jj + sibhfs.size(1) * ii] == 0) {
            for (coder::SizeType kk{0}; kk <= nvpf + 3; kk++) {
              isborder[elems[(hf_data[kk + hf_size_idx_1 * jj] + elems.size(1) *
                              ii) - 1] - 1] = true;
            }
          }
        }

        ii++;
      }
    }
  }

  static void determine_incident_halfedges(int nv, const ::coder::array<int, 2U>
    &elems, const ::coder::array<int, 2U> &sibhes, ::coder::array<int, 1U> &v2he)
  {
    coder::SizeType i;
    coder::SizeType kk;

    // DETERMINE_INCIDENT_HALFEDGES Determine an incident halfedges.
    v2he.set_size(nv);
    for (i = 0; i < nv; i++) {
      v2he[i] = 0;
    }

    kk = 0;
    while ((kk <= elems.size(0) - 1) && (elems[elems.size(1) * kk] != 0)) {
      i = elems.size(1);
      for (coder::SizeType lid{0}; lid < i; lid++) {
        coder::SizeType v;
        v = elems[lid + elems.size(1) * kk];
        if ((v > 0) && ((v2he[v - 1] == 0) || (sibhes[lid + sibhes.size(1) * kk]
              == 0))) {
          //  Encode <fid,leid> pair into a heid.
          v2he[v - 1] = static_cast<int>(static_cast<unsigned int>(kk + 1) << 4)
            + lid;
        } else {
        }
      }

      kk++;
    }
  }

  static void determine_incident_halffaces(int nv, const ::coder::array<int, 2U>
    &elems, const ::coder::array<int, 2U> &sibhfs, ::coder::array<int, 1U> &v2hf)
  {
    ::coder::array<int, 2U> b_sibhfs;
    ::coder::array<boolean_T, 1U> isborder;
    coder::SizeType i;
    coder::SizeType ncvpE;
    coder::SizeType v2f_size_idx_1;
    signed char v2f_data[32];

    // DETERMINE_INCIDENT_HALFFACES Determine an incident half-faces.
    v2hf.set_size(nv);
    for (i = 0; i < nv; i++) {
      v2hf[i] = 0;
    }

    //  We use three bits for local-face ID.
    if (elems.size(1) == 1) {
      coder::SizeType ii;
      coder::SizeType loop_ub;
      coder::SizeType offset;
      coder::SizeType offset_o;

      //  Mixed elements
      b_sibhfs.set_size(sibhfs.size(0), sibhfs.size(1));
      loop_ub = sibhfs.size(1) * sibhfs.size(0) - 1;
      for (i = 0; i <= loop_ub; i++) {
        b_sibhfs[i] = sibhfs[i];
      }

      determine_border_vertices_vol(nv, elems, b_sibhfs, isborder);
      offset = 1;
      offset_o = 0;
      ii = 1;
      while (offset < elems.size(0)) {
        coder::SizeType v;
        boolean_T b;
        boolean_T b1;
        boolean_T guard1;
        boolean_T guard2;
        boolean_T guard3;
        boolean_T guard4;
        guard1 = false;
        guard2 = false;
        guard3 = false;
        guard4 = false;
        switch (elems[(offset - 1) % elems.size(0) + (offset - 1) / elems.size(0)])
        {
         case 4:
          guard1 = true;
          break;

         case 10:
          guard1 = true;
          break;

         case 5:
          guard2 = true;
          break;

         case 14:
          guard2 = true;
          break;

         case 6:
          guard3 = true;
          break;

         case 15:
          guard3 = true;
          break;

         case 18:
          guard3 = true;
          break;

         case 8:
          guard4 = true;
          break;

         case 20:
          guard4 = true;
          break;

         case 27:
          guard4 = true;
          break;

         default:
          m2cErrMsgIdAndTxt("momp2cpp:runtimeError", "unsupported element type");
          break;
        }

        if (guard4) {
          //  Hexahedral
          b = true;
          b1 = (elems.size(0) <= 0);
          i = 0;
          for (coder::SizeType jj{0}; jj < 8; jj++) {
            loop_ub = offset + jj;
            if (b1 || (loop_ub < 0) || (loop_ub >= elems.size(0))) {
              i = 0;
              b = true;
            } else if (b) {
              b = false;
              i = loop_ub % elems.size(0) + loop_ub / elems.size(0);
            } else if (i > 2147483646) {
              i = loop_ub % elems.size(0) + loop_ub / elems.size(0);
            } else {
              i++;
              if (i > elems.size(0) - 1) {
                i = (i - elems.size(0)) + 1;
              }
            }

            v = elems[i] - 1;
            if (v2hf[elems[i] - 1] == 0) {
              if (!isborder[v]) {
                //  Encode <cid,lfid> pair into a hfid.
                v2hf[v] = ((ii << 3) + iv6[3 * jj]) - 1;
              } else {
                loop_ub = offset_o + iv6[3 * jj];
                if (sibhfs[loop_ub % sibhfs.size(0) + loop_ub / sibhfs.size(0)] ==
                    0) {
                  //  Encode <cid,lfid> pair into a hfid.
                  v2hf[v] = ((ii << 3) + iv6[3 * jj]) - 1;
                }
              }

              if (!isborder[v]) {
                //  Encode <cid,lfid> pair into a hfid.
                v2hf[v] = ((ii << 3) + iv6[3 * jj + 1]) - 1;
              } else {
                loop_ub = offset_o + iv6[3 * jj + 1];
                if (sibhfs[loop_ub % sibhfs.size(0) + loop_ub / sibhfs.size(0)] ==
                    0) {
                  //  Encode <cid,lfid> pair into a hfid.
                  v2hf[v] = ((ii << 3) + iv6[3 * jj + 1]) - 1;
                }
              }

              if (!isborder[v]) {
                //  Encode <cid,lfid> pair into a hfid.
                v2hf[v] = ((ii << 3) + iv6[3 * jj + 2]) - 1;
              } else {
                loop_ub = offset_o + iv6[3 * jj + 2];
                if (sibhfs[loop_ub % sibhfs.size(0) + loop_ub / sibhfs.size(0)] ==
                    0) {
                  //  Encode <cid,lfid> pair into a hfid.
                  v2hf[v] = ((ii << 3) + iv6[3 * jj + 2]) - 1;
                }
              }

              if (!isborder[v]) {
                //  Encode <cid,lfid> pair into a hfid.
                v2hf[v] = ((ii << 3) + iv6[3 * jj + 3]) - 1;
              } else {
                loop_ub = offset_o + iv6[3 * jj + 3];
                if (sibhfs[loop_ub % sibhfs.size(0) + loop_ub / sibhfs.size(0)] ==
                    0) {
                  //  Encode <cid,lfid> pair into a hfid.
                  v2hf[v] = ((ii << 3) + iv6[3 * jj + 3]) - 1;
                }
              }
            }
          }
        }

        if (guard3) {
          //  Prism
          b = true;
          b1 = (elems.size(0) <= 0);
          i = 0;
          for (coder::SizeType jj{0}; jj < 6; jj++) {
            loop_ub = offset + jj;
            if (b1 || (loop_ub < 0) || (loop_ub >= elems.size(0))) {
              i = 0;
              b = true;
            } else if (b) {
              b = false;
              i = loop_ub % elems.size(0) + loop_ub / elems.size(0);
            } else if (i > 2147483646) {
              i = loop_ub % elems.size(0) + loop_ub / elems.size(0);
            } else {
              i++;
              if (i > elems.size(0) - 1) {
                i = (i - elems.size(0)) + 1;
              }
            }

            v = elems[i] - 1;
            if (v2hf[elems[i] - 1] == 0) {
              if (!isborder[v]) {
                //  Encode <cid,lfid> pair into a hfid.
                v2hf[v] = ((ii << 3) + iv7[3 * jj]) - 1;
              } else {
                loop_ub = offset_o + iv7[3 * jj];
                if (sibhfs[loop_ub % sibhfs.size(0) + loop_ub / sibhfs.size(0)] ==
                    0) {
                  //  Encode <cid,lfid> pair into a hfid.
                  v2hf[v] = ((ii << 3) + iv7[3 * jj]) - 1;
                }
              }

              if (!isborder[v]) {
                //  Encode <cid,lfid> pair into a hfid.
                v2hf[v] = ((ii << 3) + iv7[3 * jj + 1]) - 1;
              } else {
                loop_ub = offset_o + iv7[3 * jj + 1];
                if (sibhfs[loop_ub % sibhfs.size(0) + loop_ub / sibhfs.size(0)] ==
                    0) {
                  //  Encode <cid,lfid> pair into a hfid.
                  v2hf[v] = ((ii << 3) + iv7[3 * jj + 1]) - 1;
                }
              }

              if (!isborder[v]) {
                //  Encode <cid,lfid> pair into a hfid.
                v2hf[v] = ((ii << 3) + iv7[3 * jj + 2]) - 1;
              } else {
                loop_ub = offset_o + iv7[3 * jj + 2];
                if (sibhfs[loop_ub % sibhfs.size(0) + loop_ub / sibhfs.size(0)] ==
                    0) {
                  //  Encode <cid,lfid> pair into a hfid.
                  v2hf[v] = ((ii << 3) + iv7[3 * jj + 2]) - 1;
                }
              }
            }
          }
        }

        if (guard2) {
          //  Pyramid
          b = true;
          b1 = (elems.size(0) <= 0);
          i = 0;
          for (coder::SizeType jj{0}; jj < 5; jj++) {
            loop_ub = offset + jj;
            if (b1 || (loop_ub < 0) || (loop_ub >= elems.size(0))) {
              i = 0;
              b = true;
            } else if (b) {
              b = false;
              i = loop_ub % elems.size(0) + loop_ub / elems.size(0);
            } else if (i > 2147483646) {
              i = loop_ub % elems.size(0) + loop_ub / elems.size(0);
            } else {
              i++;
              if (i > elems.size(0) - 1) {
                i = (i - elems.size(0)) + 1;
              }
            }

            v = elems[i] - 1;
            if (v2hf[elems[i] - 1] == 0) {
              loop_ub = (jj + 1 == 4) + 2;
              for (coder::SizeType kk{0}; kk <= loop_ub; kk++) {
                if (!isborder[v]) {
                  //  Encode <cid,lfid> pair into a hfid.
                  v2hf[v] = ((ii << 3) + iv8[kk + (jj << 2)]) - 1;
                } else {
                  v2f_size_idx_1 = offset_o + iv8[kk + (jj << 2)];
                  if (sibhfs[v2f_size_idx_1 % sibhfs.size(0) + v2f_size_idx_1 /
                      sibhfs.size(0)] == 0) {
                    //  Encode <cid,lfid> pair into a hfid.
                    v2hf[v] = ((ii << 3) + iv8[kk + (jj << 2)]) - 1;
                  }
                }
              }
            }
          }
        }

        if (guard1) {
          //  Tetrahedral
          b = true;
          b1 = (elems.size(0) <= 0);
          i = 0;
          for (coder::SizeType jj{0}; jj < 4; jj++) {
            loop_ub = offset + jj;
            if (b1 || (loop_ub < 0) || (loop_ub >= elems.size(0))) {
              i = 0;
              b = true;
            } else if (b) {
              b = false;
              i = loop_ub % elems.size(0) + loop_ub / elems.size(0);
            } else if (i > 2147483646) {
              i = loop_ub % elems.size(0) + loop_ub / elems.size(0);
            } else {
              i++;
              if (i > elems.size(0) - 1) {
                i = (i - elems.size(0)) + 1;
              }
            }

            v = elems[i] - 1;
            if (v2hf[elems[i] - 1] == 0) {
              if (!isborder[v]) {
                //  Encode <cid,lfid> pair into a hfid.
                v2hf[v] = ((ii << 3) + iv3[3 * jj]) - 1;
              } else {
                loop_ub = offset_o + iv3[3 * jj];
                if (sibhfs[loop_ub % sibhfs.size(0) + loop_ub / sibhfs.size(0)] ==
                    0) {
                  //  Encode <cid,lfid> pair into a hfid.
                  v2hf[v] = ((ii << 3) + iv3[3 * jj]) - 1;
                }
              }

              if (!isborder[v]) {
                //  Encode <cid,lfid> pair into a hfid.
                v2hf[v] = ((ii << 3) + iv3[3 * jj + 1]) - 1;
              } else {
                loop_ub = offset_o + iv3[3 * jj + 1];
                if (sibhfs[loop_ub % sibhfs.size(0) + loop_ub / sibhfs.size(0)] ==
                    0) {
                  //  Encode <cid,lfid> pair into a hfid.
                  v2hf[v] = ((ii << 3) + iv3[3 * jj + 1]) - 1;
                }
              }

              if (!isborder[v]) {
                //  Encode <cid,lfid> pair into a hfid.
                v2hf[v] = ((ii << 3) + iv3[3 * jj + 2]) - 1;
              } else {
                loop_ub = offset_o + iv3[3 * jj + 2];
                if (sibhfs[loop_ub % sibhfs.size(0) + loop_ub / sibhfs.size(0)] ==
                    0) {
                  //  Encode <cid,lfid> pair into a hfid.
                  v2hf[v] = ((ii << 3) + iv3[3 * jj + 2]) - 1;
                }
              }
            }
          }
        }

        ii++;
        offset = (offset + elems[(offset - 1) % elems.size(0) + (offset - 1) /
                  elems.size(0)]) + 1;
        offset_o = (offset_o + sibhfs[offset_o % sibhfs.size(0) + offset_o /
                    sibhfs.size(0)]) + 1;
      }
    } else {
      coder::SizeType ii;
      coder::SizeType loop_ub;
      boolean_T guard1{ false };

      boolean_T guard2{ false };

      boolean_T guard3{ false };

      boolean_T guard4{ false };

      //  Table for local IDs of incident faces of each vertex.
      guard1 = false;
      guard2 = false;
      guard3 = false;
      guard4 = false;
      switch (elems.size(1)) {
       case 4:
        guard1 = true;
        break;

       case 10:
        guard1 = true;
        break;

       case 5:
        guard2 = true;
        break;

       case 14:
        guard2 = true;
        break;

       case 6:
        guard3 = true;
        break;

       case 15:
        guard3 = true;
        break;

       case 18:
        guard3 = true;
        break;

       case 8:
        guard4 = true;
        break;

       case 20:
        guard4 = true;
        break;

       case 27:
        guard4 = true;
        break;

       default:
        m2cErrMsgIdAndTxt("momp2cpp:runtimeError", "Unsupported element type");
        ncvpE = -1;
        v2f_size_idx_1 = 3;
        for (i = 0; i < 24; i++) {
          v2f_data[i] = iv6[i];
        }
        break;
      }

      if (guard4) {
        ncvpE = 7;
        v2f_size_idx_1 = 3;
        for (i = 0; i < 24; i++) {
          v2f_data[i] = iv6[i];
        }
      }

      if (guard3) {
        ncvpE = 5;
        v2f_size_idx_1 = 3;
        for (i = 0; i < 18; i++) {
          v2f_data[i] = iv7[i];
        }
      }

      if (guard2) {
        ncvpE = 4;
        v2f_size_idx_1 = 4;
        for (i = 0; i < 20; i++) {
          v2f_data[i] = iv8[i];
        }
      }

      if (guard1) {
        ncvpE = 3;
        v2f_size_idx_1 = 3;
        for (i = 0; i < 12; i++) {
          v2f_data[i] = iv3[i];
        }
      }

      b_sibhfs.set_size(sibhfs.size(0), sibhfs.size(1));
      loop_ub = sibhfs.size(1) * sibhfs.size(0) - 1;
      for (i = 0; i <= loop_ub; i++) {
        b_sibhfs[i] = sibhfs[i];
      }

      determine_border_vertices_vol(nv, elems, b_sibhfs, isborder);

      //  Construct a vertex-to-halfface mapping.
      ii = 0;
      while ((ii + 1 <= elems.size(0)) && (elems[elems.size(1) * ii] != 0)) {
        coder::SizeType jj;
        boolean_T exitg1;
        jj = 0;
        exitg1 = false;
        while ((!exitg1) && (jj <= ncvpE)) {
          coder::SizeType v;
          v = elems[jj + elems.size(1) * ii] - 1;
          if (v + 1 == 0) {
            exitg1 = true;
          } else {
            if (v2hf[v] == 0) {
              boolean_T b2;
              if ((v2f_size_idx_1 > 3) && (v2f_data[4 * jj + 3] != 0)) {
                b2 = true;
              } else {
                b2 = false;
              }

              i = b2 + 2;
              for (coder::SizeType kk{0}; kk <= i; kk++) {
                if ((!isborder[v]) || (sibhfs[(v2f_data[kk + v2f_size_idx_1 * jj]
                      + sibhfs.size(1) * ii) - 1] == 0)) {
                  //  Encode <cid,lfid> pair into a hfid.
                  v2hf[v] = (((ii + 1) << 3) + v2f_data[kk + v2f_size_idx_1 * jj])
                    - 1;
                }
              }
            }

            jj++;
          }
        }

        ii++;
      }
    }
  }

  static void determine_sibling_halfedges(int nv, const ::coder::array<int, 2U>
    &elems, ::coder::array<int, 2U> &sibhes)
  {
    static const signed char b_iv[4]{ 2, 3, 4, 1 };

    ::coder::array<int, 1U> is_index;
    ::coder::array<int, 1U> v2he_fid;
    ::coder::array<int, 1U> v2nv;
    ::coder::array<signed char, 1U> v2he_leid;
    coder::SizeType i;
    coder::SizeType ii;
    coder::SizeType k;
    coder::SizeType nelems;
    coder::SizeType nepE;
    coder::SizeType unnamed_idx_0_tmp;
    boolean_T exitg1;
    boolean_T hasthree;

    // DETERMINE_SIBLING_HALFEDGES constructs an extended half-edge data structure
    if ((elems.size(1) == 4) || (elems.size(1) == 8) || (elems.size(1) == 9)) {
      nepE = 4;

      //  Number of edges per element
    } else {
      nepE = 3;

      //  Number of edges per element
    }

    //  First, build is_index to store starting position for each vertex.
    is_index.set_size(nv + 1);
    for (i = 0; i <= nv; i++) {
      is_index[i] = 0;
    }

    nelems = elems.size(0);
    ii = 0;
    exitg1 = false;
    while ((!exitg1) && (ii <= nelems - 1)) {
      if (elems[elems.size(1) * ii] == 0) {
        nelems = ii;
        exitg1 = true;
      } else {
        if ((nepE != 4) || (elems[elems.size(1) * ii + 3] == 0)) {
          hasthree = true;
        } else {
          hasthree = false;
        }

        i = 3 - hasthree;
        for (coder::SizeType j{0}; j <= i; j++) {
          k = elems[j + elems.size(1) * ii];
          is_index[k] = is_index[k] + 1;
        }

        ii++;
      }
    }

    is_index[0] = 1;
    for (ii = 0; ii < nv; ii++) {
      is_index[ii + 1] = is_index[ii] + is_index[ii + 1];
    }

    unnamed_idx_0_tmp = nelems * nepE;
    v2nv.set_size(unnamed_idx_0_tmp);

    //  Vertex to next vertex in each halfedge.
    v2he_fid.set_size(unnamed_idx_0_tmp);

    //  Vertex to halfedge.
    v2he_leid.set_size(unnamed_idx_0_tmp);

    //  Vertex to halfedge.
    for (ii = 0; ii < nelems; ii++) {
      if ((nepE != 4) || (elems[elems.size(1) * ii + 3] == 0)) {
        hasthree = true;
      } else {
        hasthree = false;
      }

      i = 3 - hasthree;
      for (coder::SizeType j{0}; j <= i; j++) {
        k = elems[j + elems.size(1) * ii] - 1;
        v2nv[is_index[k] - 1] = elems[(b_iv[j + (hasthree && (j + 1 == 3))] +
          elems.size(1) * ii) - 1];
        v2he_fid[is_index[k] - 1] = ii + 1;

        //  Vertex to halfedge.
        v2he_leid[is_index[k] - 1] = static_cast<signed char>(j + 1);

        //  Vertex to halfedge.
        is_index[k] = is_index[k] + 1;
      }
    }

    i = nv - 1;
    for (ii = i; ii >= 1; ii--) {
      is_index[ii] = is_index[ii - 1];
    }

    is_index[0] = 1;

    //  Set sibhes
    sibhes.set_size(nelems, nepE);
    for (i = 0; i < unnamed_idx_0_tmp; i++) {
      sibhes[i] = 0;
    }

    for (ii = 0; ii < nelems; ii++) {
      if ((nepE != 4) || (elems[elems.size(1) * ii + 3] == 0)) {
        hasthree = true;
      } else {
        hasthree = false;
      }

      i = 3 - hasthree;
      for (coder::SizeType jj{0}; jj <= i; jj++) {
        //  Process each edge only once
        if (sibhes[jj + sibhes.size(1) * ii] == 0) {
          coder::SizeType i1;
          coder::SizeType i2;
          coder::SizeType prev_heid_fid;
          coder::SizeType prev_heid_leid;
          coder::SizeType sibhes_tmp;
          coder::SizeType v;
          coder::SizeType vn;
          v = elems[jj + elems.size(1) * ii];
          vn = elems[(b_iv[jj + (hasthree && (jj + 1 == 3))] + elems.size(1) *
                      ii) - 1];
          prev_heid_fid = ii;
          prev_heid_leid = jj;

          //  LOCATE: Locate index in v2nv(first:last)
          i1 = is_index[vn - 1];
          i2 = is_index[vn] - 1;
          for (coder::SizeType b_index{i1}; b_index <= i2; b_index++) {
            if (v2nv[b_index - 1] == v) {
              //  Encode <fid,leid> pair into a heid.
              sibhes_tmp = v2he_fid[b_index - 1];
              unnamed_idx_0_tmp = v2he_leid[b_index - 1];
              sibhes[prev_heid_leid + sibhes.size(1) * prev_heid_fid] = (
                static_cast<int>(static_cast<unsigned int>(sibhes_tmp) << 4) +
                unnamed_idx_0_tmp) - 1;
              prev_heid_fid = sibhes_tmp - 1;
              prev_heid_leid = unnamed_idx_0_tmp - 1;
            }
          }

          //  Check for halfedges in the same orientation
          i1 = is_index[v - 1];
          i2 = is_index[v] - 1;
          for (coder::SizeType b_index{i1}; b_index <= i2; b_index++) {
            if (v2nv[b_index - 1] == vn) {
              unnamed_idx_0_tmp = v2he_fid[b_index - 1];
              if (unnamed_idx_0_tmp != ii + 1) {
                //  Encode <fid,leid> pair into a heid.
                sibhes_tmp = v2he_leid[b_index - 1];
                sibhes[prev_heid_leid + sibhes.size(1) * prev_heid_fid] = (
                  static_cast<int>(static_cast<unsigned int>(unnamed_idx_0_tmp) <<
                                   4) + sibhes_tmp) - 1;
                prev_heid_fid = unnamed_idx_0_tmp - 1;
                prev_heid_leid = sibhes_tmp - 1;
              }
            }
          }

          if (prev_heid_fid != ii) {
            //  Close up the cycle
            sibhes[prev_heid_leid + sibhes.size(1) * prev_heid_fid] =
              static_cast<int>(static_cast<unsigned int>(ii + 1) << 4) + jj;
          }
        }
      }
    }
  }

  static void determine_sibling_halffaces(int nv, const ::coder::array<int, 2U>
    &elems, ::coder::array<int, 2U> &sibhfs)
  {
    static const signed char hf_hex[24]{ 1, 4, 3, 2, 1, 2, 6, 5, 2, 3, 7, 6, 3,
      4, 8, 7, 1, 5, 8, 4, 5, 6, 7, 8 };

    static const signed char next[8]{ 2, 3, 1, 0, 2, 3, 4, 1 };

    static const signed char prev[8]{ 3, 1, 2, 0, 4, 1, 2, 3 };

    static const signed char b_iv[4]{ 2, 3, 4, 1 };

    static const signed char b_iv1[4]{ 4, 1, 2, 3 };

    ::coder::array<int, 1U> b_sibhfs;
    ::coder::array<int, 1U> is_index;
    ::coder::array<int, 1U> v2hf_cid;
    ::coder::array<int, 1U> v2oe_v1;
    ::coder::array<int, 1U> v2oe_v2;
    ::coder::array<signed char, 1U> v2hf_lfid;
    coder::SizeType b_index;
    coder::SizeType ex;
    coder::SizeType i;
    coder::SizeType ii;
    coder::SizeType iindx;
    coder::SizeType jj;
    coder::SizeType loop_ub;
    coder::SizeType nelems;
    coder::SizeType v;
    coder::SizeType v1;
    coder::SizeType v2;
    signed char tmp_data[4];
    boolean_T exitg7;
    boolean_T found;
    boolean_T guard1{ false };

    boolean_T guard2{ false };

    boolean_T manifold;
    boolean_T oriented;

    //  DETERMINE_SIBLING_HALFFACES determines the sibling half-face of each
    guard1 = false;
    guard2 = false;
    switch (elems.size(1)) {
     case 4:
      //  tet
      determine_sibling_halffaces_tet(nv, elems, sibhfs, &manifold, &oriented);
      break;

     case 10:
      //  tet
      determine_sibling_halffaces_tet(nv, elems, sibhfs, &manifold, &oriented);
      break;

     case 5:
      //  pyramid
      determine_sibling_halffaces_pyramid(nv, elems, sibhfs);
      break;

     case 14:
      //  pyramid
      determine_sibling_halffaces_pyramid(nv, elems, sibhfs);
      break;

     case 6:
      guard1 = true;
      break;

     case 15:
      guard1 = true;
      break;

     case 18:
      guard1 = true;
      break;

     case 8:
      guard2 = true;
      break;

     case 20:
      guard2 = true;
      break;

     case 27:
      guard2 = true;
      break;

     default:
      if (elems.size(1) != 1) {
        m2cErrMsgIdAndTxt("momp2cpp:runtimeError", "Unsupported element type.");
      }

      determine_sibling_halffaces_mixed(nv, elems, b_sibhfs);
      sibhfs.set_size(b_sibhfs.size(0), 1);
      loop_ub = b_sibhfs.size(0);
      for (i = 0; i < loop_ub; i++) {
        sibhfs[i] = b_sibhfs[i];
      }
      break;
    }

    if (guard2) {
      boolean_T exitg6;

      //  hex
      is_index.set_size(nv + 1);
      for (i = 0; i <= nv; i++) {
        is_index[i] = 0;
      }

      nelems = elems.size(0);
      ii = 0;
      exitg7 = false;
      while ((!exitg7) && (ii <= nelems - 1)) {
        if (elems[elems.size(1) * ii] == 0) {
          nelems = ii;
          exitg7 = true;
        } else {
          for (jj = 0; jj < 6; jj++) {
            loop_ub = jj << 2;
            v = elems[(hf_hex[loop_ub] + elems.size(1) * ii) - 1];
            i = elems[(hf_hex[loop_ub + 1] + elems.size(1) * ii) - 1];
            if (v < i) {
              v = i;
            }

            i = elems[(hf_hex[loop_ub + 2] + elems.size(1) * ii) - 1];
            if (v < i) {
              v = i;
            }

            i = elems[(hf_hex[loop_ub + 3] + elems.size(1) * ii) - 1];
            if (v < i) {
              v = i;
            }

            is_index[v] = is_index[v] + 1;
          }

          ii++;
        }
      }

      is_index[0] = 1;
      for (ii = 0; ii < nv; ii++) {
        is_index[ii + 1] = is_index[ii] + is_index[ii + 1];
      }

      //  v2hf stores mapping from each vertex to half-face ID.
      v2hf_cid.set_size(is_index[nv]);
      v2hf_lfid.set_size(is_index[nv]);
      v2oe_v1.set_size(is_index[nv]);
      v2oe_v2.set_size(is_index[nv]);
      for (ii = 0; ii < nelems; ii++) {
        for (jj = 0; jj < 6; jj++) {
          iindx = 0;
          loop_ub = jj << 2;
          ex = elems[(hf_hex[loop_ub] + elems.size(1) * ii) - 1] - 1;
          i = elems[(hf_hex[loop_ub + 1] + elems.size(1) * ii) - 1];
          if (ex + 1 < i) {
            ex = i - 1;
            iindx = 1;
          }

          i = elems[(hf_hex[loop_ub + 2] + elems.size(1) * ii) - 1];
          if (ex + 1 < i) {
            ex = i - 1;
            iindx = 2;
          }

          i = elems[(hf_hex[loop_ub + 3] + elems.size(1) * ii) - 1];
          if (ex + 1 < i) {
            ex = i - 1;
            iindx = 3;
          }

          v2oe_v1[is_index[ex] - 1] = elems[(hf_hex[(b_iv[iindx] + loop_ub) - 1]
            + elems.size(1) * ii) - 1];
          v2oe_v2[is_index[ex] - 1] = elems[(hf_hex[(b_iv1[iindx] + loop_ub) - 1]
            + elems.size(1) * ii) - 1];

          // v2hf(is_index(v)) = clfids2hfid(ii,jj);
          v2hf_cid[is_index[ex] - 1] = ii + 1;
          v2hf_lfid[is_index[ex] - 1] = static_cast<signed char>(jj + 1);
          is_index[ex] = is_index[ex] + 1;
        }
      }

      i = nv - 1;
      for (ii = i; ii >= 1; ii--) {
        is_index[ii] = is_index[ii - 1];
      }

      is_index[0] = 1;

      //  Fill in sibhfs for each half-face.
      sibhfs.set_size(nelems, 6);
      loop_ub = 6 * nelems;
      for (i = 0; i < loop_ub; i++) {
        sibhfs[i] = 0;
      }

      ii = 0;
      exitg6 = false;
      while ((!exitg6) && (ii <= nelems - 1)) {
        coder::SizeType exitg5;
        jj = 0;
        do {
          exitg5 = 0;
          if (jj < 6) {
            //  local face ID
            if (sibhfs[jj + sibhfs.size(1) * ii] != 0) {
              jj++;
            } else {
              //  list of vertices of face
              iindx = 0;
              loop_ub = jj << 2;
              ex = elems[(hf_hex[loop_ub] + elems.size(1) * ii) - 1];
              i = elems[(hf_hex[loop_ub + 1] + elems.size(1) * ii) - 1];
              if (ex < i) {
                ex = i;
                iindx = 1;
              }

              i = elems[(hf_hex[loop_ub + 2] + elems.size(1) * ii) - 1];
              if (ex < i) {
                ex = i;
                iindx = 2;
              }

              i = elems[(hf_hex[loop_ub + 3] + elems.size(1) * ii) - 1];
              if (ex < i) {
                ex = i;
                iindx = 3;
              }

              found = false;
              v1 = elems[(hf_hex[(b_iv1[iindx] + loop_ub) - 1] + elems.size(1) *
                          ii) - 1];
              v2 = elems[(hf_hex[(b_iv[iindx] + loop_ub) - 1] + elems.size(1) *
                          ii) - 1];

              //  Search for sibling half-face.
              b_index = is_index[ex - 1] - 1;
              exitg7 = false;
              while ((!exitg7) && (b_index + 1 <= is_index[ex] - 1)) {
                if ((v2oe_v1[b_index] == v1) && (v2oe_v2[b_index] == v2)) {
                  //  Encode <cid,lfid> pair into a hfid.
                  sibhfs[jj + 6 * ii] = ((v2hf_cid[b_index] << 3) +
                    v2hf_lfid[b_index]) - 1;

                  //  Encode <cid,lfid> pair into a hfid.
                  sibhfs[(v2hf_lfid[b_index] + 6 * (v2hf_cid[b_index] - 1)) - 1]
                    = ((ii + 1) << 3) + jj;
                  found = true;
                  exitg7 = true;
                } else {
                  b_index++;
                }
              }

              if (!found) {
                coder::SizeType exitg4;
                b_index = is_index[ex - 1] - 1;
                do {
                  exitg4 = 0;
                  if (b_index + 1 <= is_index[ex] - 1) {
                    if ((v2oe_v1[b_index] == v2) && (v2oe_v2[b_index] == v1) &&
                        (v2hf_cid[b_index] != ii + 1)) {
                      sibhfs.set_size(0, 6);
                      exitg4 = 1;
                    } else {
                      b_index++;
                    }
                  } else {
                    exitg4 = 2;
                  }
                } while (exitg4 == 0);

                if (exitg4 == 1) {
                  exitg5 = 2;
                } else {
                  jj++;
                }
              } else {
                jj++;
              }
            }
          } else {
            ii++;
            exitg5 = 1;
          }
        } while (exitg5 == 0);

        if (exitg5 != 1) {
          exitg6 = true;
        }
      }
    }

    if (guard1) {
      coder::SizeType nvpf;
      boolean_T exitg3;

      //  prism
      is_index.set_size(nv + 1);
      for (i = 0; i <= nv; i++) {
        is_index[i] = 0;
      }

      nelems = elems.size(0);
      ii = 0;
      exitg7 = false;
      while ((!exitg7) && (ii <= nelems - 1)) {
        if (elems[elems.size(1) * ii] == 0) {
          nelems = ii;
          exitg7 = true;
        } else {
          for (jj = 0; jj < 5; jj++) {
            loop_ub = (jj + 1 < 4) + 3;
            std::copy(&iv1[jj * 4], &iv1[jj * 4 + loop_ub], &tmp_data[0]);
            v = elems[(tmp_data[0] + elems.size(1) * ii) - 1];
            for (coder::SizeType k{2}; k <= loop_ub; k++) {
              i = elems[(tmp_data[k - 1] + elems.size(1) * ii) - 1];
              if (v < i) {
                v = i;
              }
            }

            is_index[v] = is_index[v] + 1;
          }

          ii++;
        }
      }

      is_index[0] = 1;
      for (ii = 0; ii < nv; ii++) {
        is_index[ii + 1] = is_index[ii] + is_index[ii + 1];
      }

      //  v2hf stores mapping from each vertex to half-face ID.
      v2hf_cid.set_size(is_index[nv]);
      v2hf_lfid.set_size(is_index[nv]);
      v2oe_v1.set_size(is_index[nv]);
      v2oe_v2.set_size(is_index[nv]);
      for (ii = 0; ii < nelems; ii++) {
        for (jj = 0; jj < 5; jj++) {
          nvpf = (jj + 1 < 4);
          loop_ub = nvpf + 3;
          std::copy(&iv1[jj * 4], &iv1[jj * 4 + loop_ub], &tmp_data[0]);
          iindx = 0;
          ex = elems[(tmp_data[0] + elems.size(1) * ii) - 1] - 1;
          for (coder::SizeType k{2}; k <= loop_ub; k++) {
            i = elems[(tmp_data[k - 1] + elems.size(1) * ii) - 1];
            if (ex + 1 < i) {
              ex = i - 1;
              iindx = k - 1;
            }
          }

          loop_ub = iindx + (nvpf << 2);
          v2oe_v1[is_index[ex] - 1] = elems[(tmp_data[next[loop_ub] - 1] +
            elems.size(1) * ii) - 1];
          v2oe_v2[is_index[ex] - 1] = elems[(tmp_data[prev[loop_ub] - 1] +
            elems.size(1) * ii) - 1];
          v2hf_cid[is_index[ex] - 1] = ii + 1;
          v2hf_lfid[is_index[ex] - 1] = static_cast<signed char>(jj + 1);
          is_index[ex] = is_index[ex] + 1;
        }
      }

      i = nv - 1;
      for (ii = i; ii >= 1; ii--) {
        is_index[ii] = is_index[ii - 1];
      }

      is_index[0] = 1;

      //  Fill in sibhfs for each half-face.
      sibhfs.set_size(nelems, 5);
      loop_ub = 5 * nelems;
      for (i = 0; i < loop_ub; i++) {
        sibhfs[i] = 0;
      }

      ii = 0;
      exitg3 = false;
      while ((!exitg3) && (ii <= nelems - 1)) {
        coder::SizeType exitg2;
        jj = 0;
        do {
          exitg2 = 0;
          if (jj < 5) {
            //  local face ID
            if (sibhfs[jj + sibhfs.size(1) * ii] != 0) {
              jj++;
            } else {
              nvpf = (jj + 1 < 4);
              loop_ub = nvpf + 3;
              for (i = 0; i < loop_ub; i++) {
                tmp_data[i] = iv1[i + (jj << 2)];
              }

              //  list of vertices of face
              iindx = 0;
              ex = elems[(tmp_data[0] + elems.size(1) * ii) - 1];
              for (coder::SizeType k{2}; k <= loop_ub; k++) {
                i = elems[(tmp_data[k - 1] + elems.size(1) * ii) - 1];
                if (ex < i) {
                  ex = i;
                  iindx = k - 1;
                }
              }

              found = false;
              loop_ub = iindx + (nvpf << 2);
              v1 = elems[(tmp_data[prev[loop_ub] - 1] + elems.size(1) * ii) - 1];
              v2 = elems[(tmp_data[next[loop_ub] - 1] + elems.size(1) * ii) - 1];

              //  Search for sibling half-face.
              b_index = is_index[ex - 1] - 1;
              exitg7 = false;
              while ((!exitg7) && (b_index + 1 <= is_index[ex] - 1)) {
                if ((v2oe_v1[b_index] == v1) && (v2oe_v2[b_index] == v2)) {
                  //  Encode <cid,lfid> pair into a hfid.
                  sibhfs[jj + 5 * ii] = ((v2hf_cid[b_index] << 3) +
                    v2hf_lfid[b_index]) - 1;

                  //  Encode <cid,lfid> pair into a hfid.
                  sibhfs[(v2hf_lfid[b_index] + 5 * (v2hf_cid[b_index] - 1)) - 1]
                    = ((ii + 1) << 3) + jj;
                  found = true;
                  exitg7 = true;
                } else {
                  b_index++;
                }
              }

              if (!found) {
                coder::SizeType exitg1;
                b_index = is_index[ex - 1] - 1;
                do {
                  exitg1 = 0;
                  if (b_index + 1 <= is_index[ex] - 1) {
                    if ((v2oe_v1[b_index] == v2) && (v2oe_v2[b_index] == v1) &&
                        (v2hf_cid[b_index] != ii + 1)) {
                      sibhfs.set_size(0, 5);
                      exitg1 = 1;
                    } else {
                      b_index++;
                    }
                  } else {
                    exitg1 = 2;
                  }
                } while (exitg1 == 0);

                if (exitg1 == 1) {
                  exitg2 = 2;
                } else {
                  jj++;
                }
              } else {
                jj++;
              }
            }
          } else {
            ii++;
            exitg2 = 1;
          }
        } while (exitg2 == 0);

        if (exitg2 != 1) {
          exitg3 = true;
        }
      }
    }
  }

  static void determine_sibling_halffaces_mixed(int nv, const ::coder::array<int,
    2U> &elems, ::coder::array<int, 1U> &sibhfs)
  {
    static const signed char v2av_hex[24]{ 2, 5, 4, 3, 6, 1, 4, 7, 2, 1, 8, 3, 6,
      8, 1, 7, 5, 2, 8, 6, 3, 5, 7, 4 };

    static const signed char v2av_pyr[20]{ 2, 5, 4, 0, 3, 5, 1, 0, 4, 5, 2, 0, 1,
      5, 3, 0, 1, 2, 3, 4 };

    static const signed char v2av_pri[18]{ 2, 4, 3, 3, 5, 1, 1, 6, 2, 6, 1, 5, 4,
      2, 6, 5, 3, 4 };

    ::coder::array<int, 2U> r;
    ::coder::array<int, 2U> r1;
    ::coder::array<int, 2U> vs_elem;
    ::coder::array<int, 1U> b_is_index_v;
    ::coder::array<int, 1U> is_index_o;
    ::coder::array<int, 1U> is_index_v;
    ::coder::array<int, 1U> v2hf;
    ::coder::array<int, 1U> v2oe_v1;
    ::coder::array<int, 1U> v2oe_v2;
    int b_vs_elem[8];
    int c_vs_elem[6];
    int d_vs_elem[5];
    int e_vs_elem[4];
    coder::SizeType i;
    coder::SizeType i1;
    coder::SizeType i2;
    coder::SizeType ii;
    coder::SizeType jj;
    coder::SizeType k;
    coder::SizeType nelems;
    coder::SizeType nf;
    coder::SizeType nfpE;
    coder::SizeType nopphf;
    coder::SizeType nvpE;
    coder::SizeType nvpf;
    coder::SizeType offset;
    coder::SizeType offset_ohf;
    coder::SizeType start_index_tmp;
    coder::SizeType v2oe_v1_tmp;
    boolean_T exitg13;
    boolean_T found;

    //  Determine the sibling half-faces of a mixed mesh.
    nf = 0;
    is_index_v.set_size(nv + 1);
    for (i = 0; i <= nv; i++) {
      is_index_v[i] = 0;
    }

    offset = 0;
    nelems = -1;
    if (elems.size(0) > 1) {
      i1 = elems.size(1);
      i2 = elems.size(0);
    }

    while (offset + 1 < elems.size(0)) {
      start_index_tmp = elems[offset % elems.size(0) * elems.size(1) + offset /
        elems.size(0)];
      b_is_index_v.set_size(elems[offset % elems.size(0) * elems.size(1) +
                            offset / elems.size(0)]);
      for (i = 0; i < start_index_tmp; i++) {
        b_is_index_v[i] = (i + offset) + 2;
      }

      r.set_size(1, b_is_index_v.size(0));
      start_index_tmp = b_is_index_v.size(0);
      for (i = 0; i < start_index_tmp; i++) {
        r[i] = b_is_index_v[i] - 1;
      }

      vs_elem.set_size(1, r.size(1));
      start_index_tmp = r.size(1);
      for (i = 0; i < start_index_tmp; i++) {
        vs_elem[i] = elems[r[i] % i2 * i1 + r[i] / i2];
      }

      r.set_size(1, vs_elem.size(1));
      start_index_tmp = vs_elem.size(1);
      for (i = 0; i < start_index_tmp; i++) {
        r[i] = vs_elem[i] + 1;
      }

      r1.set_size(1, vs_elem.size(1));
      start_index_tmp = vs_elem.size(1);
      for (i = 0; i < start_index_tmp; i++) {
        r1[i] = is_index_v[vs_elem[i]];
      }

      start_index_tmp = r.size(1);
      for (i = 0; i < start_index_tmp; i++) {
        is_index_v[r[i] - 1] = r1[i] + 3;
      }

      nf += 3 * elems[offset % elems.size(0) * elems.size(1) + offset /
        elems.size(0)];
      if ((elems[offset % elems.size(0) * elems.size(1) + offset / elems.size(0)]
           == 5) || (elems[offset % elems.size(0) * elems.size(1) + offset /
                     elems.size(0)] == 14)) {
        // pyramid
        is_index_v[vs_elem[4]] = is_index_v[vs_elem[4]] + 1;
        nf++;
      }

      nelems++;
      offset = (offset + elems[offset % elems.size(0) * elems.size(1) + offset /
                elems.size(0)]) + 1;
    }

    is_index_v[0] = 1;
    for (ii = 0; ii < nv; ii++) {
      is_index_v[ii + 1] = is_index_v[ii] + is_index_v[ii + 1];
    }

    //  v2hf stores mapping from each vertex to half-face ID.
    v2hf.set_size(nf);
    v2oe_v1.set_size(nf);
    v2oe_v2.set_size(nf);
    is_index_o.set_size(nelems + 2);
    for (i = 0; i <= nelems + 1; i++) {
      is_index_o[i] = 0;
    }

    is_index_o[0] = 1;
    offset = 1;
    nopphf = 1;
    for (ii = 0; ii <= nelems; ii++) {
      int f_is_index_v[4];
      coder::SizeType b_v2oe_v1_tmp;
      coder::SizeType c_v2oe_v1_tmp;
      boolean_T guard1{ false };

      boolean_T guard2{ false };

      boolean_T guard3{ false };

      boolean_T guard4{ false };

      nvpE = elems[(offset - 1) % elems.size(0) * elems.size(1) + (offset - 1) /
        elems.size(0)];
      guard1 = false;
      guard2 = false;
      guard3 = false;
      guard4 = false;
      switch (nvpE) {
       case 4:
        guard1 = true;
        break;

       case 10:
        guard1 = true;
        break;

       case 5:
        guard2 = true;
        break;

       case 14:
        guard2 = true;
        break;

       case 6:
        guard3 = true;
        break;

       case 15:
        guard3 = true;
        break;

       case 18:
        guard3 = true;
        break;

       case 8:
        guard4 = true;
        break;

       case 20:
        guard4 = true;
        break;

       case 27:
        guard4 = true;
        break;

       default:
        nfpE = 0;

        m2cErrMsgIdAndTxt("momp2cpp:runtimeError", "Unrecognized element type.");
        break;
      }

      if (guard4) {
        int c_is_index_v[8];
        int start_index[8];
        i = elems.size(0);
        for (i1 = 0; i1 < 8; i1++) {
          i2 = i1 + offset;
          i2 = elems[i2 % i * elems.size(1) + i2 / i];
          b_vs_elem[i1] = i2;
          start_index_tmp = is_index_v[i2 - 1];
          start_index[i1] = start_index_tmp;
          c_is_index_v[i1] = start_index_tmp + 3;
        }

        start_index_tmp = (ii + 1) << 3;
        for (jj = 0; jj < 8; jj++) {
          is_index_v[b_vs_elem[jj] - 1] = c_is_index_v[jj];
          v2oe_v1_tmp = b_vs_elem[v2av_hex[3 * jj] - 1];
          i = start_index[jj];
          v2oe_v1[i - 1] = v2oe_v1_tmp;
          i1 = 3 * jj + 1;
          b_v2oe_v1_tmp = b_vs_elem[v2av_hex[i1] - 1];
          v2oe_v1[i] = b_v2oe_v1_tmp;
          i2 = 3 * jj + 2;
          c_v2oe_v1_tmp = b_vs_elem[v2av_hex[i2] - 1];
          v2oe_v1[i + 1] = c_v2oe_v1_tmp;
          v2oe_v2[i - 1] = b_v2oe_v1_tmp;
          v2oe_v2[i] = c_v2oe_v1_tmp;
          v2oe_v2[i + 1] = v2oe_v1_tmp;

          //  Encode <cid,lfid> pair into a hfid.
          v2hf[i - 1] = (iv6[3 * jj] + start_index_tmp) - 1;
          v2hf[i] = (iv6[i1] + start_index_tmp) - 1;
          v2hf[i + 1] = (iv6[i2] + start_index_tmp) - 1;
        }

        nfpE = 6;
      }

      if (guard3) {
        int b_start_index[6];
        int d_is_index_v[6];
        i = elems.size(0);
        for (i1 = 0; i1 < 6; i1++) {
          i2 = i1 + offset;
          i2 = elems[i2 % i * elems.size(1) + i2 / i];
          c_vs_elem[i1] = i2;
          start_index_tmp = is_index_v[i2 - 1];
          b_start_index[i1] = start_index_tmp;
          d_is_index_v[i1] = start_index_tmp + 3;
        }

        start_index_tmp = (ii + 1) << 3;
        for (jj = 0; jj < 6; jj++) {
          is_index_v[c_vs_elem[jj] - 1] = d_is_index_v[jj];
          v2oe_v1_tmp = c_vs_elem[v2av_pri[3 * jj] - 1];
          i = b_start_index[jj];
          v2oe_v1[i - 1] = v2oe_v1_tmp;
          i1 = 3 * jj + 1;
          b_v2oe_v1_tmp = c_vs_elem[v2av_pri[i1] - 1];
          v2oe_v1[i] = b_v2oe_v1_tmp;
          i2 = 3 * jj + 2;
          c_v2oe_v1_tmp = c_vs_elem[v2av_pri[i2] - 1];
          v2oe_v1[i + 1] = c_v2oe_v1_tmp;
          v2oe_v2[i - 1] = b_v2oe_v1_tmp;
          v2oe_v2[i] = c_v2oe_v1_tmp;
          v2oe_v2[i + 1] = v2oe_v1_tmp;

          //  Encode <cid,lfid> pair into a hfid.
          v2hf[i - 1] = (iv7[3 * jj] + start_index_tmp) - 1;
          v2hf[i] = (iv7[i1] + start_index_tmp) - 1;
          v2hf[i + 1] = (iv7[i2] + start_index_tmp) - 1;
        }

        nfpE = 5;
      }

      if (guard2) {
        int c_start_index[5];
        int e_is_index_v[5];
        i = elems.size(0);
        for (i1 = 0; i1 < 5; i1++) {
          i2 = i1 + offset;
          i2 = elems[i2 % i * elems.size(1) + i2 / i];
          d_vs_elem[i1] = i2;
          start_index_tmp = is_index_v[i2 - 1];
          c_start_index[i1] = start_index_tmp;
          e_is_index_v[i1] = start_index_tmp + 3;
        }

        for (i = 0; i < 5; i++) {
          is_index_v[d_vs_elem[i] - 1] = e_is_index_v[i];
        }

        is_index_v[d_vs_elem[4] - 1] = is_index_v[d_vs_elem[4] - 1] + 1;
        start_index_tmp = (ii + 1) << 3;
        for (jj = 0; jj < 4; jj++) {
          i = jj << 2;
          v2oe_v1_tmp = d_vs_elem[v2av_pyr[i] - 1];
          i1 = c_start_index[jj];
          v2oe_v1[i1 - 1] = v2oe_v1_tmp;
          b_v2oe_v1_tmp = d_vs_elem[v2av_pyr[i + 1] - 1];
          v2oe_v1[i1] = b_v2oe_v1_tmp;
          c_v2oe_v1_tmp = d_vs_elem[v2av_pyr[i + 2] - 1];
          v2oe_v1[i1 + 1] = c_v2oe_v1_tmp;
          v2oe_v2[i1 - 1] = b_v2oe_v1_tmp;
          v2oe_v2[i1] = c_v2oe_v1_tmp;
          v2oe_v2[i1 + 1] = v2oe_v1_tmp;

          //  Encode <cid,lfid> pair into a hfid.
          v2hf[i1 - 1] = (iv8[i] + start_index_tmp) - 1;
          v2hf[i1] = (iv8[i + 1] + start_index_tmp) - 1;
          v2hf[i1 + 1] = (iv8[i + 2] + start_index_tmp) - 1;
          f_is_index_v[jj] = jj + c_start_index[4];
        }

        v2oe_v1[f_is_index_v[0] - 1] = d_vs_elem[0];
        v2oe_v1[f_is_index_v[1] - 1] = d_vs_elem[1];
        v2oe_v1[f_is_index_v[2] - 1] = d_vs_elem[2];
        v2oe_v1[f_is_index_v[3] - 1] = d_vs_elem[3];
        v2oe_v2[f_is_index_v[0] - 1] = d_vs_elem[1];
        v2oe_v2[f_is_index_v[1] - 1] = d_vs_elem[2];
        v2oe_v2[f_is_index_v[2] - 1] = d_vs_elem[3];
        v2oe_v2[f_is_index_v[3] - 1] = d_vs_elem[0];

        //  Encode <cid,lfid> pair into a hfid.
        v2hf[f_is_index_v[0] - 1] = start_index_tmp + 1;
        v2hf[f_is_index_v[1] - 1] = start_index_tmp + 2;
        v2hf[f_is_index_v[2] - 1] = start_index_tmp + 3;
        v2hf[f_is_index_v[3] - 1] = start_index_tmp + 4;
        nfpE = 5;
      }

      if (guard1) {
        int d_start_index[4];
        e_vs_elem[0] = elems[offset % elems.size(0) * elems.size(1) + offset /
          elems.size(0)];
        e_vs_elem[1] = elems[(offset + 1) % elems.size(0) * elems.size(1) +
          (offset + 1) / elems.size(0)];
        e_vs_elem[2] = elems[(offset + 2) % elems.size(0) * elems.size(1) +
          (offset + 2) / elems.size(0)];
        e_vs_elem[3] = elems[(offset + 3) % elems.size(0) * elems.size(1) +
          (offset + 3) / elems.size(0)];
        start_index_tmp = is_index_v[e_vs_elem[0] - 1];
        d_start_index[0] = start_index_tmp;
        f_is_index_v[0] = start_index_tmp + 3;
        start_index_tmp = is_index_v[e_vs_elem[1] - 1];
        d_start_index[1] = start_index_tmp;
        f_is_index_v[1] = start_index_tmp + 3;
        start_index_tmp = is_index_v[e_vs_elem[2] - 1];
        d_start_index[2] = start_index_tmp;
        f_is_index_v[2] = start_index_tmp + 3;
        start_index_tmp = is_index_v[e_vs_elem[3] - 1];
        d_start_index[3] = start_index_tmp;
        f_is_index_v[3] = start_index_tmp + 3;
        start_index_tmp = (ii + 1) << 3;
        for (jj = 0; jj < 4; jj++) {
          is_index_v[e_vs_elem[jj] - 1] = f_is_index_v[jj];
          v2oe_v1_tmp = e_vs_elem[iv2[3 * jj] - 1];
          i = d_start_index[jj];
          v2oe_v1[i - 1] = v2oe_v1_tmp;
          i1 = 3 * jj + 1;
          b_v2oe_v1_tmp = e_vs_elem[iv2[i1] - 1];
          v2oe_v1[i] = b_v2oe_v1_tmp;
          i2 = 3 * jj + 2;
          c_v2oe_v1_tmp = e_vs_elem[iv2[i2] - 1];
          v2oe_v1[i + 1] = c_v2oe_v1_tmp;
          v2oe_v2[i - 1] = b_v2oe_v1_tmp;
          v2oe_v2[i] = c_v2oe_v1_tmp;
          v2oe_v2[i + 1] = v2oe_v1_tmp;

          //  Encode <cid,lfid> pair into a hfid.
          v2hf[i - 1] = (iv3[3 * jj] + start_index_tmp) - 1;
          v2hf[i] = (iv3[i1] + start_index_tmp) - 1;
          v2hf[i + 1] = (iv3[i2] + start_index_tmp) - 1;
        }

        nfpE = 4;
      }

      offset = (offset + nvpE) + 1;
      nopphf += nfpE;
      is_index_o[ii + 1] = (is_index_o[ii] + nfpE) + 1;
    }

    if (nv - 1 < 1) {
      start_index_tmp = 0;
    } else {
      start_index_tmp = nv - 1;
    }

    i = (nv >= 2);
    b_is_index_v.set_size(start_index_tmp);
    for (i1 = 0; i1 < start_index_tmp; i1++) {
      b_is_index_v[i1] = is_index_v[i1];
    }

    start_index_tmp = b_is_index_v.size(0);
    for (i1 = 0; i1 < start_index_tmp; i1++) {
      is_index_v[i + i1] = b_is_index_v[i1];
    }

    is_index_v[0] = 1;

    //  Fill in sibhfs for each half-face.
    if (elems.size(0) == 0) {
      start_index_tmp = nopphf + nelems;
      sibhfs.set_size(start_index_tmp);
      for (i = 0; i < start_index_tmp; i++) {
        sibhfs[i] = 0;
      }
    } else {
      sibhfs.set_size(elems.size(0));
      start_index_tmp = elems.size(0);
      for (i = 0; i < start_index_tmp; i++) {
        sibhfs[i] = 0;
      }
    }

    offset = 1;
    offset_ohf = 1;
    ii = 0;
    exitg13 = false;
    while ((!exitg13) && (ii <= nelems)) {
      coder::SizeType b_index;
      boolean_T b_guard1{ false };

      boolean_T b_guard2{ false };

      boolean_T b_guard3{ false };

      boolean_T b_guard4{ false };

      boolean_T guard5{ false };

      boolean_T guard6{ false };

      nvpE = elems[(offset - 1) % elems.size(0) * elems.size(1) + (offset - 1) /
        elems.size(0)];
      b_guard1 = false;
      b_guard2 = false;
      b_guard3 = false;
      b_guard4 = false;
      guard5 = false;
      switch (nvpE) {
       case 4:
        b_guard2 = true;
        break;

       case 10:
        b_guard2 = true;
        break;

       case 5:
        b_guard3 = true;
        break;

       case 14:
        b_guard3 = true;
        break;

       case 6:
        b_guard4 = true;
        break;

       case 15:
        b_guard4 = true;
        break;

       case 18:
        b_guard4 = true;
        break;

       case 8:
        guard5 = true;
        break;

       case 20:
        guard5 = true;
        break;

       case 27:
        guard5 = true;
        break;

       default:
        nfpE = 0;

        m2cErrMsgIdAndTxt("momp2cpp:runtimeError", "Unrecognized element type.");
        b_guard1 = true;
        break;
      }

      if (guard5) {
        coder::SizeType exitg12;
        nfpE = 6;
        i = elems.size(0);
        for (i1 = 0; i1 < 8; i1++) {
          i2 = i1 + offset;
          b_vs_elem[i1] = elems[i2 % i * elems.size(1) + i2 / i];
        }

        jj = 0;
        do {
          exitg12 = 0;
          if (jj < 6) {
            //  local face ID
            guard6 = false;
            if (sibhfs[ii] != 0) {
              jj++;
            } else {
              coder::SizeType exitg10;

              //  list of vertices of face
              found = false;

              //  Search for sibling half-face.
              b_index = is_index_v[b_vs_elem[0] - 1] - 1;
              do {
                exitg10 = 0;
                i = is_index_v[b_vs_elem[0]] - 1;
                if (b_index + 1 <= i) {
                  if ((v2oe_v1[b_index] == b_vs_elem[1]) && (v2oe_v2[b_index] ==
                       b_vs_elem[3])) {
                    sibhfs[offset_ohf] = v2hf[b_index];

                    k = is_index_o[static_cast<int>(static_cast<unsigned int>
                      (v2hf[b_index]) >> 3) - 1] + static_cast<coder::SizeType>(v2hf[b_index]
                      & 7U);
                    exitg10 = 1;
                  } else {
                    b_index++;
                  }
                } else {
                  exitg10 = 2;
                }
              } while (exitg10 == 0);

              if (exitg10 == 1) {
                if (sibhfs[k] != 0) {
                  sibhfs.set_size(0);
                  exitg12 = 2;
                } else {
                  //  Encode <cid,lfid> pair into a hfid.
                  sibhfs[k] = (ii + 1) << 3;
                  found = true;
                  guard6 = true;
                }
              } else {
                guard6 = true;
              }
            }

            if (guard6) {
              if (!found) {
                coder::SizeType exitg11;
                b_index = is_index_v[b_vs_elem[0] - 1] - 1;
                do {
                  exitg11 = 0;
                  if (b_index + 1 <= i) {
                    if ((v2oe_v1[b_index] == b_vs_elem[3]) && (v2oe_v2[b_index] ==
                         b_vs_elem[1]) && (static_cast<coder::SizeType>(static_cast<unsigned
                          int>(v2hf[b_index]) >> 3) != ii + 1)) {
                      sibhfs.set_size(0);
                      exitg11 = 1;
                    } else {
                      b_index++;
                    }
                  } else {
                    exitg11 = 2;
                  }
                } while (exitg11 == 0);

                if (exitg11 == 1) {
                  exitg12 = 2;
                } else {
                  jj++;
                }
              } else {
                jj++;
              }
            }
          } else {
            exitg12 = 1;
          }
        } while (exitg12 == 0);

        if (exitg12 == 1) {
          b_guard1 = true;
        } else {
          exitg13 = true;
        }
      }

      if (b_guard4) {
        coder::SizeType exitg9;
        nfpE = 5;
        i = elems.size(0);
        for (i1 = 0; i1 < 6; i1++) {
          i2 = i1 + offset;
          c_vs_elem[i1] = elems[i2 % i * elems.size(1) + i2 / i];
        }

        jj = 0;
        do {
          exitg9 = 0;
          if (jj < 5) {
            //  local face ID
            i = offset_ohf + jj;
            guard6 = false;
            if (sibhfs[i] != 0) {
              jj++;
            } else {
              coder::SizeType exitg7;
              nvpf = (jj + 1 < 4) + 3;

              //  list of vertices of face
              found = false;

              //  Search for sibling half-face.
              start_index_tmp = jj << 2;
              v2oe_v1_tmp = c_vs_elem[iv1[start_index_tmp] - 1];
              b_index = is_index_v[v2oe_v1_tmp - 1] - 1;
              do {
                exitg7 = 0;
                i1 = is_index_v[v2oe_v1_tmp] - 1;
                if (b_index + 1 <= i1) {
                  if ((v2oe_v1[b_index] == c_vs_elem[iv1[((nvpf - 1) % nvpf +
                        (nvpf - 1) / nvpf * nvpf) + start_index_tmp] - 1]) &&
                      (v2oe_v2[b_index] == c_vs_elem[iv1[1 % nvpf +
                       start_index_tmp] - 1])) {
                    sibhfs[i] = v2hf[b_index];

                    k = is_index_o[static_cast<int>(static_cast<unsigned int>
                      (v2hf[b_index]) >> 3) - 1] + static_cast<coder::SizeType>(v2hf[b_index]
                      & 7U);
                    exitg7 = 1;
                  } else {
                    b_index++;
                  }
                } else {
                  exitg7 = 2;
                }
              } while (exitg7 == 0);

              if (exitg7 == 1) {
                if (sibhfs[k] != 0) {
                  sibhfs.set_size(0);
                  exitg9 = 2;
                } else {
                  //  Encode <cid,lfid> pair into a hfid.
                  sibhfs[k] = ((ii + 1) << 3) + jj;
                  found = true;
                  guard6 = true;
                }
              } else {
                guard6 = true;
              }
            }

            if (guard6) {
              if (!found) {
                coder::SizeType exitg8;
                b_index = is_index_v[v2oe_v1_tmp - 1] - 1;
                do {
                  exitg8 = 0;
                  if (b_index + 1 <= i1) {
                    if ((v2oe_v1[b_index] == c_vs_elem[iv1[1 % nvpf + (jj << 2)]
                         - 1]) && (v2oe_v2[b_index] == c_vs_elem[iv1[((nvpf - 1)
                          % nvpf + (nvpf - 1) / nvpf * nvpf) + (jj << 2)] - 1]) &&
                        (static_cast<int>(static_cast<coder::SizeType>(v2hf[b_index])
                                          >> 3) != ii + 1)) {
                      sibhfs.set_size(0);
                      exitg8 = 1;
                    } else {
                      b_index++;
                    }
                  } else {
                    exitg8 = 2;
                  }
                } while (exitg8 == 0);

                if (exitg8 == 1) {
                  exitg9 = 2;
                } else {
                  jj++;
                }
              } else {
                jj++;
              }
            }
          } else {
            exitg9 = 1;
          }
        } while (exitg9 == 0);

        if (exitg9 == 1) {
          b_guard1 = true;
        } else {
          exitg13 = true;
        }
      }

      if (b_guard3) {
        coder::SizeType exitg6;
        i = elems.size(0);
        for (i1 = 0; i1 < 5; i1++) {
          i2 = i1 + offset;
          d_vs_elem[i1] = elems[i2 % i * elems.size(1) + i2 / i];
        }

        nfpE = 5;
        jj = 0;
        do {
          exitg6 = 0;
          if (jj < 5) {
            //  local face ID
            i = offset_ohf + jj;
            guard6 = false;
            if (sibhfs[i] != 0) {
              jj++;
            } else {
              coder::SizeType exitg4;
              nvpf = (jj + 1 == 1) + 3;

              //  list of vertices of face
              found = false;

              //  Search for sibling half-face.
              start_index_tmp = jj << 2;
              v2oe_v1_tmp = d_vs_elem[iv5[start_index_tmp] - 1];
              b_index = is_index_v[v2oe_v1_tmp - 1] - 1;
              do {
                exitg4 = 0;
                i1 = is_index_v[v2oe_v1_tmp] - 1;
                if (b_index + 1 <= i1) {
                  if ((v2oe_v1[b_index] == d_vs_elem[iv5[((nvpf - 1) % nvpf +
                        (nvpf - 1) / nvpf * nvpf) + start_index_tmp] - 1]) &&
                      (v2oe_v2[b_index] == d_vs_elem[iv5[1 % nvpf +
                       start_index_tmp] - 1])) {
                    sibhfs[i] = v2hf[b_index];

                    k = is_index_o[static_cast<int>(static_cast<unsigned int>
                      (v2hf[b_index]) >> 3) - 1] + static_cast<coder::SizeType>(v2hf[b_index]
                      & 7U);
                    exitg4 = 1;
                  } else {
                    b_index++;
                  }
                } else {
                  exitg4 = 2;
                }
              } while (exitg4 == 0);

              if (exitg4 == 1) {
                if (sibhfs[k] != 0) {
                  sibhfs.set_size(0);
                  exitg6 = 2;
                } else {
                  //  Encode <cid,lfid> pair into a hfid.
                  sibhfs[k] = ((ii + 1) << 3) + jj;
                  found = true;
                  guard6 = true;
                }
              } else {
                guard6 = true;
              }
            }

            if (guard6) {
              //  check whether the surface is oriented
              if (!found) {
                coder::SizeType exitg5;
                b_index = is_index_v[v2oe_v1_tmp - 1] - 1;
                do {
                  exitg5 = 0;
                  if (b_index + 1 <= i1) {
                    if ((v2oe_v1[b_index] == d_vs_elem[iv5[1 % nvpf + (jj << 2)]
                         - 1]) && (v2oe_v2[b_index] == d_vs_elem[iv5[((nvpf - 1)
                          % nvpf + (nvpf - 1) / nvpf * nvpf) + (jj << 2)] - 1]) &&
                        (static_cast<int>(static_cast<coder::SizeType>(v2hf[b_index])
                                          >> 3) != ii + 1)) {
                      sibhfs.set_size(0);
                      exitg5 = 1;
                    } else {
                      b_index++;
                    }
                  } else {
                    exitg5 = 2;
                  }
                } while (exitg5 == 0);

                if (exitg5 == 1) {
                  exitg6 = 2;
                } else {
                  jj++;
                }
              } else {
                jj++;
              }
            }
          } else {
            exitg6 = 1;
          }
        } while (exitg6 == 0);

        if (exitg6 == 1) {
          b_guard1 = true;
        } else {
          exitg13 = true;
        }
      }

      if (b_guard2) {
        coder::SizeType exitg3;
        e_vs_elem[0] = elems[offset % elems.size(0) * elems.size(1) + offset /
          elems.size(0)];
        e_vs_elem[1] = elems[(offset + 1) % elems.size(0) * elems.size(1) +
          (offset + 1) / elems.size(0)];
        e_vs_elem[2] = elems[(offset + 2) % elems.size(0) * elems.size(1) +
          (offset + 2) / elems.size(0)];
        e_vs_elem[3] = elems[(offset + 3) % elems.size(0) * elems.size(1) +
          (offset + 3) / elems.size(0)];
        nfpE = 4;
        jj = 0;
        do {
          exitg3 = 0;
          if (jj < 4) {
            //  local face ID
            i = offset_ohf + jj;
            guard6 = false;
            if (sibhfs[i] != 0) {
              jj++;
            } else {
              coder::SizeType exitg1;

              //  list of vertices of face
              found = false;

              //  Search for sibling half-face.
              start_index_tmp = e_vs_elem[iv4[3 * jj] - 1];
              b_index = is_index_v[start_index_tmp - 1] - 1;
              do {
                exitg1 = 0;
                i1 = is_index_v[start_index_tmp] - 1;
                if (b_index + 1 <= i1) {
                  if ((v2oe_v1[b_index] == e_vs_elem[iv4[3 * jj + 2] - 1]) &&
                      (v2oe_v2[b_index] == e_vs_elem[iv4[3 * jj + 1] - 1])) {
                    sibhfs[i] = v2hf[b_index];

                    k = is_index_o[static_cast<int>(static_cast<unsigned int>
                      (v2hf[b_index]) >> 3) - 1] + static_cast<coder::SizeType>(v2hf[b_index]
                      & 7U);
                    exitg1 = 1;
                  } else {
                    b_index++;
                  }
                } else {
                  exitg1 = 2;
                }
              } while (exitg1 == 0);

              if (exitg1 == 1) {
                if (sibhfs[k] != 0) {
                  sibhfs.set_size(0);
                  exitg3 = 2;
                } else {
                  //  Encode <cid,lfid> pair into a hfid.
                  sibhfs[k] = ((ii + 1) << 3) + jj;
                  found = true;
                  guard6 = true;
                }
              } else {
                guard6 = true;
              }
            }

            if (guard6) {
              //  check whether the surface is oriented
              if (!found) {
                coder::SizeType exitg2;
                b_index = is_index_v[start_index_tmp - 1] - 1;
                do {
                  exitg2 = 0;
                  if (b_index + 1 <= i1) {
                    if ((v2oe_v1[b_index] == e_vs_elem[iv4[3 * jj + 1] - 1]) &&
                        (v2oe_v2[b_index] == e_vs_elem[iv4[3 * jj + 2] - 1]) &&
                        (static_cast<int>(static_cast<coder::SizeType>(v2hf[b_index])
                                          >> 3) != ii + 1)) {
                      sibhfs.set_size(0);
                      exitg2 = 1;
                    } else {
                      b_index++;
                    }
                  } else {
                    exitg2 = 2;
                  }
                } while (exitg2 == 0);

                if (exitg2 == 1) {
                  exitg3 = 2;
                } else {
                  jj++;
                }
              } else {
                jj++;
              }
            }
          } else {
            exitg3 = 1;
          }
        } while (exitg3 == 0);

        if (exitg3 == 1) {
          b_guard1 = true;
        } else {
          exitg13 = true;
        }
      }

      if (b_guard1) {
        offset = (offset + nvpE) + 1;
        sibhfs[offset_ohf - 1] = nfpE;
        offset_ohf = (offset_ohf + nfpE) + 1;
        ii++;
      }
    }
  }

  static void determine_sibling_halffaces_pyramid(int nv, const ::coder::array<
    int, 2U> &elems, ::coder::array<int, 2U> &sibhfs)
  {
    static const signed char next[8]{ 2, 3, 1, 0, 2, 3, 4, 1 };

    static const signed char prev[8]{ 3, 1, 2, 0, 4, 1, 2, 3 };

    ::coder::array<int, 1U> is_index;
    ::coder::array<int, 1U> v2hf_cid;
    ::coder::array<int, 1U> v2oe_v1;
    ::coder::array<int, 1U> v2oe_v2;
    ::coder::array<signed char, 1U> v2hf_lfid;
    coder::SizeType ex;
    coder::SizeType i;
    coder::SizeType ii;
    coder::SizeType iindx;
    coder::SizeType jj;
    coder::SizeType jj_idx_0;
    coder::SizeType nelems;
    coder::SizeType nvpf;
    signed char tmp_data[4];
    boolean_T exitg1;

    // DETERMINE_SIBLING_HALFFACE_PYRAMID Determine the sibling half-face.
    is_index.set_size(nv + 1);
    for (i = 0; i <= nv; i++) {
      is_index[i] = 0;
    }

    nelems = elems.size(0) - 1;
    ii = 0;
    exitg1 = false;
    while ((!exitg1) && (ii <= nelems)) {
      if (elems[elems.size(1) * ii] == 0) {
        nelems = ii - 1;
        exitg1 = true;
      } else {
        for (jj = 0; jj < 5; jj++) {
          coder::SizeType v;
          jj_idx_0 = (jj + 1 == 1) + 3;
          std::copy(&iv5[jj * 4], &iv5[jj * 4 + jj_idx_0], &tmp_data[0]);
          v = elems[(tmp_data[0] + elems.size(1) * ii) - 1];
          for (coder::SizeType k{2}; k <= jj_idx_0; k++) {
            i = elems[(tmp_data[k - 1] + elems.size(1) * ii) - 1];
            if (v < i) {
              v = i;
            }
          }

          is_index[v] = is_index[v] + 1;
        }

        ii++;
      }
    }

    is_index[0] = 1;
    for (ii = 0; ii < nv; ii++) {
      is_index[ii + 1] = is_index[ii] + is_index[ii + 1];
    }

    //  v2hf stores mapping from each vertex to half-face ID.
    v2hf_cid.set_size(is_index[nv]);
    v2hf_lfid.set_size(is_index[nv]);
    v2oe_v1.set_size(is_index[nv]);
    v2oe_v2.set_size(is_index[nv]);
    for (ii = 0; ii <= nelems; ii++) {
      for (jj = 0; jj < 5; jj++) {
        nvpf = (jj + 1 < 4);
        jj_idx_0 = nvpf + 3;
        std::copy(&iv5[jj * 4], &iv5[jj * 4 + jj_idx_0], &tmp_data[0]);
        iindx = 0;
        ex = elems[(tmp_data[0] + elems.size(1) * ii) - 1] - 1;
        for (coder::SizeType k{2}; k <= jj_idx_0; k++) {
          i = elems[(tmp_data[k - 1] + elems.size(1) * ii) - 1];
          if (ex + 1 < i) {
            ex = i - 1;
            iindx = k - 1;
          }
        }

        jj_idx_0 = iindx + (nvpf << 2);
        v2oe_v1[is_index[ex] - 1] = elems[(tmp_data[next[jj_idx_0] - 1] +
          elems.size(1) * ii) - 1];
        v2oe_v2[is_index[ex] - 1] = elems[(tmp_data[prev[jj_idx_0] - 1] +
          elems.size(1) * ii) - 1];

        // v2hf(is_index(v)) = clfids2hfid(ii,jj);
        v2hf_cid[is_index[ex] - 1] = ii + 1;
        v2hf_lfid[is_index[ex] - 1] = static_cast<signed char>(jj + 1);
        is_index[ex] = is_index[ex] + 1;
      }
    }

    i = nv - 1;
    for (ii = i; ii >= 1; ii--) {
      is_index[ii] = is_index[ii - 1];
    }

    is_index[0] = 1;

    //  Fill in sibhfs for each half-face.
    sibhfs.set_size(elems.size(0), elems.size(1));
    jj_idx_0 = elems.size(1) * elems.size(0);
    for (i = 0; i < jj_idx_0; i++) {
      sibhfs[i] = 0;
    }

    ii = 0;
    exitg1 = false;
    while ((!exitg1) && (ii <= nelems)) {
      coder::SizeType exitg3;
      jj = 0;
      do {
        exitg3 = 0;
        if (jj < 5) {
          //  local face ID
          if (sibhfs[jj + sibhfs.size(1) * ii] != 0) {
            jj++;
          } else {
            coder::SizeType b_index;
            coder::SizeType v1;
            coder::SizeType v2;
            boolean_T exitg4;
            boolean_T found;
            nvpf = (jj + 1 == 1);
            jj_idx_0 = nvpf + 3;
            for (i = 0; i < jj_idx_0; i++) {
              tmp_data[i] = iv5[i + (jj << 2)];
            }

            //  list of vertices of face
            iindx = 0;
            ex = elems[(tmp_data[0] + elems.size(1) * ii) - 1];
            for (coder::SizeType k{2}; k <= jj_idx_0; k++) {
              i = elems[(tmp_data[k - 1] + elems.size(1) * ii) - 1];
              if (ex < i) {
                ex = i;
                iindx = k - 1;
              }
            }

            found = false;
            jj_idx_0 = iindx + (nvpf << 2);
            v1 = elems[(tmp_data[prev[jj_idx_0] - 1] + elems.size(1) * ii) - 1];
            v2 = elems[(tmp_data[next[jj_idx_0] - 1] + elems.size(1) * ii) - 1];

            //  Search for sibling half-face.
            b_index = is_index[ex - 1] - 1;
            exitg4 = false;
            while ((!exitg4) && (b_index + 1 <= is_index[ex] - 1)) {
              if ((v2oe_v1[b_index] == v1) && (v2oe_v2[b_index] == v2)) {
                //  Encode <cid,lfid> pair into a hfid.
                sibhfs[jj + sibhfs.size(1) * ii] = ((v2hf_cid[b_index] << 3) +
                  v2hf_lfid[b_index]) - 1;

                //  Encode <cid,lfid> pair into a hfid.
                sibhfs[(v2hf_lfid[b_index] + sibhfs.size(1) * (v2hf_cid[b_index]
                         - 1)) - 1] = ((ii + 1) << 3) + jj;
                found = true;
                exitg4 = true;
              } else {
                b_index++;
              }
            }

            if (!found) {
              coder::SizeType exitg2;
              b_index = is_index[ex - 1] - 1;
              do {
                exitg2 = 0;
                if (b_index + 1 <= is_index[ex] - 1) {
                  if ((v2oe_v1[b_index] == v2) && (v2oe_v2[b_index] == v1) &&
                      (v2hf_cid[b_index] != ii + 1)) {
                    sibhfs.set_size(0, 5);
                    exitg2 = 1;
                  } else {
                    b_index++;
                  }
                } else {
                  exitg2 = 2;
                }
              } while (exitg2 == 0);

              if (exitg2 == 1) {
                exitg3 = 2;
              } else {
                jj++;
              }
            } else {
              jj++;
            }
          }
        } else {
          ii++;
          exitg3 = 1;
        }
      } while (exitg3 == 0);

      if (exitg3 != 1) {
        exitg1 = true;
      }
    }
  }

  static void determine_sibling_halffaces_tet(int nv, const ::coder::array<int,
    2U> &elems, ::coder::array<int, 2U> &sibhfs, boolean_T *manifold, boolean_T *
    oriented)
  {
    static const signed char b_iv[3]{ 3, 1, 2 };

    static const signed char b_iv1[3]{ 2, 3, 1 };

    ::coder::array<int, 1U> is_index;
    ::coder::array<int, 1U> v2hf_cid;
    ::coder::array<int, 1U> v2oe_v1;
    ::coder::array<int, 1U> v2oe_v2;
    ::coder::array<signed char, 1U> v2hf_lfid;
    coder::SizeType i;
    coder::SizeType i1;
    coder::SizeType i2;
    coder::SizeType ii;
    coder::SizeType nelems;
    coder::SizeType unnamed_idx_0;
    coder::SizeType v;
    boolean_T exitg1;

    // DETERMINE_SIBLING_HALFFACE_TET Determine the sibling half-faces.
    *manifold = true;
    *oriented = true;

    //  First, build is_index to store starting position for each vertex.
    is_index.set_size(nv + 1);
    for (i = 0; i <= nv; i++) {
      is_index[i] = 0;
    }

    nelems = elems.size(0);
    ii = 0;
    exitg1 = false;
    while ((!exitg1) && (ii <= nelems - 1)) {
      if (elems[elems.size(1) * ii] == 0) {
        nelems = ii;
        exitg1 = true;
      } else {
        for (coder::SizeType jj{0}; jj < 4; jj++) {
          v = elems[(iv4[3 * jj] + elems.size(1) * ii) - 1];
          i = elems[(iv4[3 * jj + 1] + elems.size(1) * ii) - 1];
          if (v < i) {
            v = i;
          }

          i = elems[(iv4[3 * jj + 2] + elems.size(1) * ii) - 1];
          if (v < i) {
            v = i;
          }

          is_index[v] = is_index[v] + 1;
        }

        ii++;
      }
    }

    is_index[0] = 1;
    for (ii = 0; ii < nv; ii++) {
      is_index[ii + 1] = is_index[ii] + is_index[ii + 1];
    }

    //  Store dimensions of objects.
    unnamed_idx_0 = nelems * 12;
    v2hf_cid.set_size(unnamed_idx_0);
    v2hf_lfid.set_size(unnamed_idx_0);
    v2oe_v1.set_size(unnamed_idx_0);
    v2oe_v2.set_size(unnamed_idx_0);
    for (ii = 0; ii < nelems; ii++) {
      for (coder::SizeType jj{0}; jj < 4; jj++) {
        v = elems[jj + elems.size(1) * ii] - 1;
        i = elems[(iv2[3 * jj] + elems.size(1) * ii) - 1];
        if (v + 1 > i) {
          i1 = elems[(iv2[3 * jj + 1] + elems.size(1) * ii) - 1];
          if (v + 1 > i1) {
            v2oe_v1[is_index[v] - 1] = i;
            v2oe_v2[is_index[v] - 1] = i1;
            v2hf_cid[is_index[v] - 1] = ii + 1;
            v2hf_lfid[is_index[v] - 1] = iv3[3 * jj];
            is_index[v] = is_index[v] + 1;
          }
        }

        i1 = 3 * jj + 1;
        i2 = elems[(iv2[i1] + elems.size(1) * ii) - 1];
        if (v + 1 > i2) {
          unnamed_idx_0 = elems[(iv2[3 * jj + 2] + elems.size(1) * ii) - 1];
          if (v + 1 > unnamed_idx_0) {
            v2oe_v1[is_index[v] - 1] = i2;
            v2oe_v2[is_index[v] - 1] = unnamed_idx_0;
            v2hf_cid[is_index[v] - 1] = ii + 1;
            v2hf_lfid[is_index[v] - 1] = iv3[i1];
            is_index[v] = is_index[v] + 1;
          }
        }

        i1 = 3 * jj + 2;
        i2 = elems[(iv2[i1] + elems.size(1) * ii) - 1];
        if ((v + 1 > i2) && (v + 1 > i)) {
          v2oe_v1[is_index[v] - 1] = i2;
          v2oe_v2[is_index[v] - 1] = i;
          v2hf_cid[is_index[v] - 1] = ii + 1;
          v2hf_lfid[is_index[v] - 1] = iv3[i1];
          is_index[v] = is_index[v] + 1;
        }
      }
    }

    i = nv - 1;
    for (ii = i; ii >= 1; ii--) {
      is_index[ii] = is_index[ii - 1];
    }

    is_index[0] = 1;

    //  Fill in sibhfs for each half-face.
    sibhfs.set_size(elems.size(0), elems.size(1));
    unnamed_idx_0 = elems.size(1) * elems.size(0);
    for (i = 0; i < unnamed_idx_0; i++) {
      sibhfs[i] = 0;
    }

    for (ii = 0; ii < nelems; ii++) {
      for (coder::SizeType jj{0}; jj < 4; jj++) {
        //  local face ID
        if (sibhfs[jj + sibhfs.size(1) * ii] == 0) {
          coder::SizeType iindx;
          coder::SizeType nhfs;
          coder::SizeType prev_cid;
          signed char prev_lfid;
          signed char sibhfs_tmp;

          //  list of vertices of face
          iindx = 0;
          unnamed_idx_0 = elems[(iv4[3 * jj] + elems.size(1) * ii) - 1];
          i = elems[(iv4[3 * jj + 1] + elems.size(1) * ii) - 1];
          if (unnamed_idx_0 < i) {
            unnamed_idx_0 = i;
            iindx = 1;
          }

          i = elems[(iv4[3 * jj + 2] + elems.size(1) * ii) - 1];
          if (unnamed_idx_0 < i) {
            unnamed_idx_0 = i;
            iindx = 2;
          }

          prev_cid = ii;
          prev_lfid = static_cast<signed char>(jj + 1);
          nhfs = 0;

          //  Search for half-face in the opposite orientation
          i = is_index[unnamed_idx_0 - 1];
          i1 = is_index[unnamed_idx_0] - 1;
          for (coder::SizeType b_index{i}; b_index <= i1; b_index++) {
            if ((v2oe_v1[b_index - 1] == elems[(iv4[(b_iv[iindx] + 3 * jj) - 1]
                  + elems.size(1) * ii) - 1]) && (v2oe_v2[b_index - 1] == elems
                 [(iv4[(b_iv1[iindx] + 3 * jj) - 1] + elems.size(1) * ii) - 1]))
            {
              //  Encode <cid,lfid> pair into a hfid.
              unnamed_idx_0 = v2hf_cid[b_index - 1];
              sibhfs_tmp = v2hf_lfid[b_index - 1];
              sibhfs[(prev_lfid + sibhfs.size(1) * prev_cid) - 1] =
                ((unnamed_idx_0 << 3) + sibhfs_tmp) - 1;
              prev_cid = unnamed_idx_0 - 1;
              prev_lfid = sibhfs_tmp;
              nhfs++;
            }
          }

          //  Check for halfface in the same orientation
          for (coder::SizeType b_index{i}; b_index <= i1; b_index++) {
            if ((v2oe_v1[b_index - 1] == elems[(iv4[(b_iv1[iindx] + 3 * jj) - 1]
                  + elems.size(1) * ii) - 1]) && (v2oe_v2[b_index - 1] == elems
                 [(iv4[(b_iv[iindx] + 3 * jj) - 1] + elems.size(1) * ii) - 1]))
            {
              i2 = v2hf_cid[b_index - 1];
              if (i2 != ii + 1) {
                //  Encode <cid,lfid> pair into a hfid.
                sibhfs_tmp = v2hf_lfid[b_index - 1];
                sibhfs[(prev_lfid + sibhfs.size(1) * prev_cid) - 1] = ((i2 << 3)
                  + sibhfs_tmp) - 1;
                prev_cid = i2 - 1;
                prev_lfid = sibhfs_tmp;
                nhfs++;
                *oriented = false;
              }
            }
          }

          if (prev_cid != ii) {
            //  Close up the cycle
            sibhfs[(prev_lfid + sibhfs.size(1) * prev_cid) - 1] = ((ii + 1) << 3)
              + jj;
            nhfs++;
          }

          if ((*manifold) && (nhfs > 2)) {
            *manifold = false;
            *oriented = false;
          }
        }
      }
    }
  }

  static void determine_sibling_halfverts(int nv, const ::coder::array<int, 2U>
    &edges, ::coder::array<int, 2U> &sibhvs)
  {
    ::coder::array<int, 1U> is_index;
    ::coder::array<int, 1U> v2hv;
    coder::SizeType b_v2hv_tmp;
    coder::SizeType i;
    coder::SizeType ii;
    coder::SizeType nedgs;
    coder::SizeType v2hv_tmp;
    boolean_T exitg1;

    //  DETERMINE_SIBLING_HALFVERTS determines the sibling half-vertices for each vertex.
    is_index.set_size(nv + 1);
    for (i = 0; i <= nv; i++) {
      is_index[i] = 0;
    }

    nedgs = edges.size(0);
    ii = 0;
    exitg1 = false;
    while ((!exitg1) && (ii <= nedgs - 1)) {
      if (edges[edges.size(1) * ii] == 0) {
        nedgs = ii;
        exitg1 = true;
      } else {
        is_index[edges[edges.size(1) * ii]] = is_index[edges[edges.size(1) * ii]]
          + 1;
        is_index[edges[edges.size(1) * ii + 1]] = is_index[edges[edges.size(1) *
          ii + 1]] + 1;
        ii++;
      }
    }

    is_index[0] = 1;
    for (ii = 0; ii < nv; ii++) {
      is_index[ii + 1] = is_index[ii] + is_index[ii + 1];
    }

    v2hv.set_size(nedgs << 1);

    //  Vertex to half-vertex.
    for (ii = 0; ii < nedgs; ii++) {
      //  Encode <eid,lvid> pair into a hvid.
      i = edges[edges.size(1) * ii];
      v2hv_tmp = is_index[i - 1];
      b_v2hv_tmp = (ii + 1) << 1;
      v2hv[v2hv_tmp - 1] = b_v2hv_tmp;
      is_index[i - 1] = v2hv_tmp + 1;

      //  Encode <eid,lvid> pair into a hvid.
      i = edges[edges.size(1) * ii + 1];
      v2hv_tmp = is_index[i - 1];
      v2hv[v2hv_tmp - 1] = b_v2hv_tmp + 1;
      is_index[i - 1] = v2hv_tmp + 1;
    }

    i = nv - 1;
    for (ii = i; ii >= 1; ii--) {
      is_index[ii] = is_index[ii - 1];
    }

    is_index[0] = 1;

    //  Set sibhvs
    sibhvs.set_size(edges.size(0), 2);
    v2hv_tmp = edges.size(0) << 1;
    for (i = 0; i < v2hv_tmp; i++) {
      sibhvs[i] = 0;
    }

    for (coder::SizeType v{0}; v < nv; v++) {
      v2hv_tmp = is_index[v + 1];
      b_v2hv_tmp = v2hv_tmp - 1;
      i = is_index[v];
      if (v2hv_tmp - 1 > is_index[v]) {
        coder::SizeType hvid_prev;

        //  The vertex has two or more incident halfedges
        hvid_prev = v2hv[v2hv_tmp - 2];
        for (ii = i; ii <= b_v2hv_tmp; ii++) {
          v2hv_tmp = v2hv[ii - 1];
          sibhvs[static_cast<int>(hvid_prev & 1U) + 2 * (static_cast<int>(
            static_cast<coder::SizeType>(hvid_prev) >> 1) - 1)] = v2hv_tmp;
          hvid_prev = v2hv_tmp;
        }
      }
    }
  }

  static int elem_one_ring(int vid, const ::coder::array<int, 2U> &tets, const ::
    coder::array<int, 2U> &sibhfs, const ::coder::array<int, 1U> &v2hf, ::coder::
    array<boolean_T, 1U> &etags, int ngbes[128])
  {
    int stack[128];
    coder::SizeType nelems;
    nelems = 0;

    //  This is similar to append_one_ring.m but this function only finds the 1-ring of elements.
    if (v2hf[vid - 1] != 0) {
      coder::SizeType size_stack;
      boolean_T exitg1;
      boolean_T overflow;

      stack[0] = static_cast<int>(static_cast<coder::SizeType>(v2hf[vid - 1]) >> 3);

      //  Element (region) ID
      overflow = false;

      //  Create a stack for storing tets
      size_stack = 1;

      //  sibhfs_tet(lvid, :) gives the faces that border on local vertex lvid.
      exitg1 = false;
      while ((!exitg1) && (size_stack > 0)) {
        coder::SizeType rid;

        //  Pop the element from top of stack
        rid = stack[size_stack - 1] - 1;
        size_stack--;

        //  Append element
        if (nelems >= 128) {
          m2cPrintf("Overflow in elements in append_one_ring.m \n");
          fflush(stdout);
          overflow = true;
        } else if (!etags[rid]) {
          coder::SizeType ii;
          coder::SizeType lvid;
          etags[rid] = true;
          nelems++;
          ngbes[nelems - 1] = rid + 1;
          lvid = -1;

          //  Stores which vertex vid is within the tetrahedron.
          if (tets[tets.size(1) * rid] == vid) {
            lvid = 0;
          }

          if (tets[tets.size(1) * rid + 1] == vid) {
            lvid = 1;
          }

          if (tets[tets.size(1) * rid + 2] == vid) {
            lvid = 2;
          }

          if (tets[tets.size(1) * rid + 3] == vid) {
            lvid = 3;
          }

          //  Push unvisited neighbor tets onto stack
          ii = 0;
          while ((ii < 3) && (lvid + 1 != 0)) {
            unsigned int c;

            c = static_cast<coder::SizeType>(sibhfs[(iv9[ii + 3 * lvid] +
              sibhfs.size(1) * rid) - 1]) >> 3;
            if ((static_cast<int>(c) != 0) && (!etags[static_cast<coder::SizeType>(c) - 1]))
            {
              if (size_stack >= 128) {
                m2cPrintf("Overflow in stack in append_one_ring.m \n");
                fflush(stdout);
                overflow = true;
              } else {
                size_stack++;
                stack[size_stack - 1] = static_cast<int>(c);
              }
            }

            ii++;
          }
        }

        if (overflow) {
          exitg1 = true;
        }
      }
    }

    return nelems;
  }

  static void obtain_nring_curv(int vid, double ring, const ::coder::array<int,
    2U> &edgs, const ::coder::array<int, 2U> &sibhvs, const ::coder::array<int,
    1U> &v2hv, int ngbvs[128], ::coder::array<boolean_T, 1U> &vtags, ::coder::
    array<boolean_T, 1U> &etags, int ngbes[128], int *nverts, int *nedges)
  {
    int ebuf[13];
    coder::SizeType b_c;
    unsigned int c;

    // OBTAIN_NRING_CURV Collect n-ring vertices and edges.
    c = static_cast<unsigned int>(v2hv[vid - 1]) >> 1;

    b_c = static_cast<int>(v2hv[vid - 1] & 1U);
    *nverts = 0;
    *nedges = 0;
    if (static_cast<int>(c) != 0) {
      int hvbuf[13];
      coder::SizeType opp;

      //  Collect one-ring vertices and edges
      ngbvs[0] = edgs[(edgs.size(1) * (static_cast<coder::SizeType>(c) - 1) - b_c) + 1];
      *nverts = 1;
      hvbuf[0] = sibhvs[(2 * (static_cast<coder::SizeType>(c) - 1) - b_c) + 1];
      ebuf[0] = static_cast<int>(c);
      *nedges = 1;
      opp = sibhvs[b_c + 2 * (static_cast<coder::SizeType>(c) - 1)];
      if (opp != 0) {
        c = static_cast<unsigned int>(opp) >> 1;

        b_c = static_cast<int>(opp & 1U);
        ngbvs[1] = edgs[(edgs.size(1) * (static_cast<coder::SizeType>(c) - 1) - b_c) + 1];
        *nverts = 2;
        hvbuf[1] = sibhvs[(2 * (static_cast<coder::SizeType>(c) - 1) - b_c) + 1];
        ebuf[1] = static_cast<int>(c);
        *nedges = 2;
      }

      if (ring <= 1.0) {
        std::copy(&ebuf[0], &ebuf[*nedges], &ngbes[0]);
      } else {
        coder::SizeType cur_ring;
        coder::SizeType i;
        coder::SizeType nverts_pre;
        vtags[vid - 1] = true;
        for (i = 0; i < *nverts; i++) {
          vtags[ngbvs[i] - 1] = true;
        }

        for (i = 0; i < *nedges; i++) {
          etags[ebuf[i] - 1] = true;
        }

        //  Define buffers and prepare tags for further processing
        nverts_pre = 0;

        //  Second, build full-size ring
        cur_ring = 1;
        ring = std::fmin(6.0, ring);
        coder::SizeType exitg1;
        coder::SizeType nedges_pre;
        coder::SizeType nverts_last;
        do {
          exitg1 = 0;

          //  Collect next level of ring
          nverts_last = *nverts;
          nedges_pre = *nedges;
          i = nverts_pre + 1;
          for (coder::SizeType ii{i}; ii <= nverts_last; ii++) {
            c = static_cast<unsigned int>(hvbuf[ii - 1]) >> 1;

            //  If the edge has already been inserted, then the vertex must be
            if ((static_cast<int>(c) != 0) && (!etags[static_cast<coder::SizeType>(c) - 1]))
            {
              coder::SizeType v;

              //  Insert edge into list
              (*nedges)++;
              ebuf[*nedges - 1] = static_cast<int>(c);
              etags[static_cast<coder::SizeType>(c) - 1] = true;

              b_c = static_cast<int>(hvbuf[ii - 1] & 1U);
              v = edgs[(edgs.size(1) * (static_cast<coder::SizeType>(c) - 1) - b_c) + 1];
              if (!vtags[v - 1]) {
                (*nverts)++;
                ngbvs[*nverts - 1] = v;
                vtags[v - 1] = true;

                //  Insert opposite halfvertex
                hvbuf[*nverts - 1] = sibhvs[(2 * (static_cast<coder::SizeType>(c) - 1) - b_c)
                  + 1];
              }
            }
          }

          cur_ring++;
          if ((cur_ring >= ring) || (*nedges == nedges_pre)) {
            exitg1 = 1;
          } else {
            nverts_pre = nverts_last;
          }
        } while (exitg1 == 0);

        //  Reset flags
        if (*nedges - 1 >= 0) {
          std::copy(&ebuf[0], &ebuf[*nedges], &ngbes[0]);
        }

        vtags[vid - 1] = false;
        for (i = 0; i < *nverts; i++) {
          vtags[ngbvs[i] - 1] = false;
        }

        for (i = 0; i < *nedges; i++) {
          etags[ebuf[i] - 1] = false;
        }
      }
    }
  }

  static void obtain_nring_quad(int vid, double ring, const ::coder::array<int,
    2U> &elems, const ::coder::array<int, 2U> &sibhes, const ::coder::array<int,
    1U> &v2he, int ngbvs[1024], ::coder::array<boolean_T, 1U> &vtags, ::coder::
    array<boolean_T, 1U> &ftags, int ngbfs[1024], int *nverts, int *nfaces)
  {
    static const signed char nxt[8]{ 2, 3, 1, 0, 2, 3, 4, 1 };

    static const signed char prv[8]{ 3, 1, 2, 0, 4, 1, 2, 3 };

    unsigned int c;
    coder::SizeType c_tmp;
    coder::SizeType fid;
    coder::SizeType lid;
    boolean_T overflow;

    //  OBTAIN_NRING_QUAD Collect n-ring vertices of a quad or mixed mesh.
    c = static_cast<unsigned int>(v2he[vid - 1]) >> 4;
    fid = static_cast<int>(static_cast<unsigned int>(v2he[vid - 1]) >> 4);

    c_tmp = static_cast<int>(v2he[vid - 1] & 15U);
    lid = c_tmp;
    *nverts = 0;
    *nfaces = 0;
    overflow = false;
    if (static_cast<int>(c) != 0) {
      int hebuf[1024];
      coder::SizeType exitg1;
      coder::SizeType fid_in;
      coder::SizeType lid_prv;
      coder::SizeType opp;
      boolean_T b;

      //  Optimized version for collecting one-ring vertices
      if (sibhes[c_tmp + sibhes.size(1) * (static_cast<coder::SizeType>(c) - 1)] != 0) {
        fid_in = static_cast<int>(c);
      } else {
        fid_in = 0;

        //  If vertex is border edge, insert its incident border vertex.
        if ((elems.size(1) == 4) && (elems[elems.size(1) * (static_cast<coder::SizeType>(c)
              - 1) + 3] != 0)) {
          b = true;
        } else {
          b = false;
        }

        ngbvs[0] = elems[(nxt[c_tmp + (b << 2)] + elems.size(1) * (static_cast<
          coder::SizeType>(c) - 1)) - 1];
        *nverts = 1;
        hebuf[0] = 0;
      }

      //  Rotate counterclockwise order around vertex and insert vertices
      do {
        exitg1 = 0;
        if ((elems.size(1) == 4) && (elems[elems.size(1) * (fid - 1) + 3] != 0))
        {
          b = true;
        } else {
          b = false;
        }

        lid_prv = prv[lid + (b << 2)] - 1;
        if ((*nverts < 1024) && (*nfaces < 1024)) {
          (*nverts)++;
          ngbvs[*nverts - 1] = elems[lid_prv + elems.size(1) * (fid - 1)];

          //  Save a starting edge for newly inserted vertex to allow early
          hebuf[*nverts - 1] = static_cast<int>(static_cast<unsigned int>(fid) <<
            4) + lid_prv;
          (*nfaces)++;
          ngbfs[*nfaces - 1] = fid;
        } else {
          overflow = true;
        }

        opp = sibhes[lid_prv + sibhes.size(1) * (fid - 1)];

        fid = static_cast<int>(static_cast<unsigned int>(opp) >> 4);
        if (static_cast<int>(static_cast<unsigned int>(opp) >> 4) == fid_in) {
          exitg1 = 1;
        } else {
          lid = static_cast<int>(opp & 15U);
        }
      } while (exitg1 == 0);

      //  Finished cycle
      if ((ring != 1.0) || ((*nverts < 1) && (*nfaces < 1024))) {
        double cur_ring;
        double ring_full;
        coder::SizeType nfaces_pre;
        coder::SizeType nverts_pre;

        // assert(nargin==9);
        vtags[vid - 1] = true;
        for (coder::SizeType i{0}; i < *nverts; i++) {
          vtags[ngbvs[i] - 1] = true;
        }

        for (coder::SizeType i{0}; i < *nfaces; i++) {
          ftags[ngbfs[i] - 1] = true;
        }

        //  Define buffers and prepare tags for further processing
        nverts_pre = 0;
        nfaces_pre = 0;

        //  Second, build full-size ring
        ring_full = std::trunc(ring);
        cur_ring = 1.0;
        boolean_T guard1{ false };

        do {
          coder::SizeType b_i;
          coder::SizeType nverts_last;
          coder::SizeType v;
          boolean_T guard2{ false };

          exitg1 = 0;
          guard1 = false;
          if ((cur_ring > ring_full) || ((cur_ring == ring_full) && (ring_full
                != ring))) {
            coder::SizeType nfaces_last;

            //  Collect halfring
            nfaces_last = *nfaces;
            nverts_last = *nverts;
            b_i = nfaces_pre + 1;
            for (coder::SizeType ii{b_i}; ii <= nfaces_last; ii++) {
              coder::SizeType jj;
              guard2 = false;
              if (elems.size(1) == 4) {
                c_tmp = ngbfs[ii - 1] - 1;
                if (elems[elems.size(1) * c_tmp + 3] != 0) {
                  boolean_T exitg4;

                  //  Insert missing vertices in quads
                  jj = 0;
                  exitg4 = false;
                  while ((!exitg4) && (jj < 4)) {
                    v = elems[jj + elems.size(1) * c_tmp];
                    if (!vtags[v - 1]) {
                      if (*nverts >= 1024) {
                        overflow = true;
                      } else {
                        (*nverts)++;
                        ngbvs[*nverts - 1] = v;
                        vtags[v - 1] = true;
                      }

                      exitg4 = true;
                    } else {
                      jj++;
                    }
                  }
                } else {
                  guard2 = true;
                }
              } else {
                guard2 = true;
              }

              if (guard2) {
                boolean_T exitg3;

                //  take opposite vertex in opposite face of triangle
                jj = 0;
                exitg3 = false;
                while ((!exitg3) && (jj < 3)) {
                  coder::SizeType oppe;
                  oppe = sibhes[jj + sibhes.size(1) * (ngbfs[ii - 1] - 1)];

                  c = static_cast<unsigned int>(oppe) >> 4;
                  if (oppe != 0) {
                    b = ftags[static_cast<coder::SizeType>(c) - 1];
                    if (!b) {
                      if ((elems.size(1) == 4) && (elems[elems.size(1) * (
                            static_cast<coder::SizeType>(c) - 1) + 3] != 0)) {
                        b = true;
                      } else {
                        b = false;
                      }

                      v = elems[(prv[static_cast<int>(oppe & 15U) + (b << 2)] +
                                 elems.size(1) * (static_cast<coder::SizeType>(c) - 1)) - 1]
                        - 1;
                      if (overflow || ((!vtags[v]) && (*nverts >= 1024)) ||
                          (*nfaces >= 1024)) {
                        overflow = true;
                      } else {
                        overflow = false;
                      }

                      if (!overflow) {
                        (*nfaces)++;
                        ngbfs[*nfaces - 1] = static_cast<int>(c);
                        ftags[static_cast<coder::SizeType>(c) - 1] = true;
                      }

                      if ((!vtags[v]) && (!overflow)) {
                        (*nverts)++;
                        ngbvs[*nverts - 1] = v + 1;
                        vtags[v] = true;
                      }

                      exitg3 = true;
                    } else {
                      jj++;
                    }
                  } else {
                    jj++;
                  }
                }
              }
            }

            if ((*nverts >= 1) || (*nfaces >= 1024) || (*nfaces == nfaces_last))
            {
              exitg1 = 1;
            } else {
              //  If needs to expand, then undo the last half ring
              *nverts = nverts_last;
              b_i = nfaces_last + 1;
              for (coder::SizeType i{b_i}; i <= *nfaces; i++) {
                ftags[ngbfs[i - 1] - 1] = false;
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
            b_i = nverts_pre + 1;
            for (coder::SizeType ii{b_i}; ii <= nverts_last; ii++) {
              coder::SizeType nedges;
              boolean_T allow_early_term;
              boolean_T isfirst;

              c = static_cast<unsigned int>(v2he[ngbvs[ii - 1] - 1]) >> 4;
              fid = static_cast<int>(static_cast<unsigned int>(v2he[ngbvs[ii - 1]
                - 1]) >> 4) - 1;

              lid = static_cast<int>(v2he[ngbvs[ii - 1] - 1] & 15U);
              if ((elems.size(1) == 4) && (elems[elems.size(1) * (static_cast<
                    coder::SizeType>(c) - 1) + 3] != 0)) {
                b = true;
              } else {
                b = false;
              }

              nedges = b;

              //  Allow early termination of the loop if an incident halfedge
              c_tmp = hebuf[ii - 1];
              if ((c_tmp != 0) && (sibhes[static_cast<int>(v2he[ngbvs[ii - 1] -
                    1] & 15U) + sibhes.size(1) * (static_cast<coder::SizeType>(c) - 1)] != 0))
              {
                allow_early_term = true;

                fid = static_cast<int>(static_cast<unsigned int>(hebuf[ii - 1]) >>
                  4) - 1;

                lid = static_cast<int>(c_tmp & 15U);
                if ((elems.size(1) == 4) && (elems[elems.size(1) * (static_cast<
                      int>(static_cast<coder::SizeType>(hebuf[ii - 1]) >> 4) - 1) +
                     3] != 0)) {
                  b = true;
                } else {
                  b = false;
                }

                nedges = b;
              } else {
                allow_early_term = false;
              }

              //  Starting point of counterclockwise rotation
              if (sibhes[lid + sibhes.size(1) * fid] != 0) {
                fid_in = fid;
              } else {
                fid_in = -1;
              }

              c_tmp = nxt[lid + (nedges << 2)];
              v = elems[(c_tmp + elems.size(1) * fid) - 1] - 1;
              if (overflow || ((!vtags[v]) && (*nverts >= 1024))) {
                overflow = true;
              } else {
                overflow = false;
              }

              if ((!overflow) && (!vtags[v])) {
                (*nverts)++;
                ngbvs[*nverts - 1] = v + 1;
                vtags[v] = true;

                //  Save starting position for next vertex
                hebuf[*nverts - 1] = (static_cast<int>(static_cast<unsigned int>
                  (fid + 1) << 4) + c_tmp) - 1;
              }

              //  Rotate counterclockwise around the vertex.
              isfirst = true;
              coder::SizeType exitg2;
              do {
                exitg2 = 0;

                //  Insert vertx into list
                if ((elems.size(1) == 4) && (elems[elems.size(1) * fid + 3] != 0))
                {
                  b = true;
                } else {
                  b = false;
                }

                lid_prv = prv[lid + (b << 2)] - 1;

                //  Insert face into list
                guard2 = false;
                if (ftags[fid]) {
                  if (allow_early_term && (!isfirst)) {
                    exitg2 = 1;
                  } else {
                    guard2 = true;
                  }
                } else {
                  //  If the face has already been inserted, then the vertex
                  v = elems[lid_prv + elems.size(1) * fid] - 1;
                  if (overflow || ((!vtags[v]) && (*nverts >= 1024)) ||
                      ((!ftags[fid]) && (*nfaces >= 1024))) {
                    overflow = true;
                  } else {
                    overflow = false;
                  }

                  if ((!vtags[v]) && (!overflow)) {
                    (*nverts)++;
                    ngbvs[*nverts - 1] = v + 1;
                    vtags[v] = true;

                    //  Save starting position for next ring
                    hebuf[*nverts - 1] = static_cast<int>(static_cast<unsigned
                      int>(fid + 1) << 4) + lid_prv;
                  }

                  if ((!ftags[fid]) && (!overflow)) {
                    (*nfaces)++;
                    ngbfs[*nfaces - 1] = fid + 1;
                    ftags[fid] = true;
                  }

                  isfirst = false;
                  guard2 = true;
                }

                if (guard2) {
                  opp = sibhes[lid_prv + sibhes.size(1) * fid];

                  fid = static_cast<int>(static_cast<unsigned int>(opp) >> 4) -
                    1;
                  if (static_cast<int>(static_cast<unsigned int>(opp) >> 4) ==
                      fid_in + 1) {
                    //  Finished cycle
                    exitg2 = 1;
                  } else {
                    lid = static_cast<int>(opp & 15U);
                  }
                }
              } while (exitg2 == 0);
            }

            cur_ring++;
            if (((*nverts >= 1) && (cur_ring >= ring)) || (*nfaces == nfaces_pre)
                || overflow) {
              exitg1 = 1;
            } else {
              nverts_pre = nverts_last;
            }
          }
        } while (exitg1 == 0);

        //  Reset flags
        vtags[vid - 1] = false;
        for (coder::SizeType i{0}; i < *nverts; i++) {
          vtags[ngbvs[i] - 1] = false;
        }

        for (coder::SizeType i{0}; i < *nfaces; i++) {
          ftags[ngbfs[i] - 1] = false;
        }
      }
    }
  }

  static void obtain_nring_vol(int vid, double ring, const ::coder::array<int,
    2U> &tets, const ::coder::array<int, 2U> &sibhfs, const ::coder::array<int,
    1U> &v2hf, int ngbvs[2048], int ngbes[2048], ::coder::array<boolean_T, 1U>
    &vtags, ::coder::array<boolean_T, 1U> &etags, ::coder::array<boolean_T, 1U>
    &etags_elem, int *nverts, int *nelems)
  {
    ::coder::array<int, 1U> temp_verts;
    int ngbes_lcl[128];
    boolean_T overflow;

    // OBTAIN_NRING_VOL Collect n-ring neighbor vertics and elements.
    overflow = false;

    //  For collecting the next layer of nodes before adding them to the
    vtags[vid - 1] = true;
    *nverts = 0;
    *nelems = 0;

    //  nverts is number of vertices in the neighborhood (not in entire mesh)
    if (v2hf[vid - 1] != 0) {
      coder::SizeType i;
      coder::SizeType nverts_idx_0;
      boolean_T guard1{ false };

      //  Obtain the 1-ring
      append_one_ring(vid, tets, sibhfs, v2hf, ngbvs, vtags, etags, ngbes,
                      nverts, nelems);

      // {
      guard1 = false;
      if (ring == 1.0) {
        if (*nverts >= 1) {
          vtags[vid - 1] = false;
          for (i = 0; i < *nverts; i++) {
            vtags[ngbvs[i] - 1] = false;
          }

          if (*nelems < 1) {
            nverts_idx_0 = 0;
          } else {
            nverts_idx_0 = *nelems;
          }

          for (i = 0; i < nverts_idx_0; i++) {
            etags[ngbes[i] - 1] = false;
          }
        } else {
          //  otherwise expand to a 1.3 ring
          ring = 1.3;
          guard1 = true;
        }
      } else {
        guard1 = true;
      }

      if (guard1) {
        double cur_ring;
        double ring_full;
        coder::SizeType nelems_pre;
        coder::SizeType nverts_pre;
        boolean_T onethird;
        boolean_T twothird;

        // }
        nverts_pre = 1;
        nelems_pre = 0;
        cur_ring = 1.0;
        ring = std::fmin(8.0, ring);
        ring_full = std::trunc(ring);
        onethird = false;
        twothird = false;
        if (ring != ring_full) {
          if (std::round(ring - ring_full) == 0.0) {
            onethird = true;
          } else {
            twothird = true;
          }
        }

        coder::SizeType exitg1;
        boolean_T b_guard1{ false };

        boolean_T guard2{ false };

        do {
          coder::SizeType MAX_temp_verts;
          coder::SizeType eid;
          coder::SizeType ii;
          coder::SizeType jj;
          coder::SizeType nelems_lcl;
          coder::SizeType ntemp;
          coder::SizeType nverts_last;
          boolean_T exitg2;
          exitg1 = 0;
          b_guard1 = false;
          guard2 = false;
          if ((cur_ring > ring_full) || ((cur_ring == ring_full) && onethird)) {
            //  Collect one-third-ring
            nverts_last = *nverts;
            MAX_temp_verts = 12 * *nelems;
            temp_verts.set_size(MAX_temp_verts);
            for (i = 0; i < MAX_temp_verts; i++) {
              temp_verts[i] = 0;
            }

            ntemp = -1;

            // Check currently listed neighboring elements
            i = nelems_pre + 1;
            for (ii = i; ii <= *nelems; ii++) {
              coder::SizeType i1;
              coder::SizeType i2;
              coder::SizeType i3;
              nverts_idx_0 = tets[tets.size(1) * (ii - 1)];
              i1 = tets[tets.size(1) * (ii - 1) + 1];
              i2 = tets[tets.size(1) * (ii - 1) + 2];
              i3 = tets[tets.size(1) * (ii - 1) + 3];
              if (((vtags[nverts_idx_0 - 1] + vtags[i1 - 1]) + vtags[i2 - 1]) +
                  vtags[i3 - 1] == 3) {
                temp_verts[ntemp + 1] = nverts_idx_0;
                temp_verts[ntemp + 2] = i1;
                temp_verts[ntemp + 3] = i2;
                temp_verts[ntemp + 4] = i3;
                ntemp += 4;
              }
            }

            // Check latest ring of neighboring elements
            ii = nverts_pre;
            exitg2 = false;
            while ((!exitg2) && (ii <= *nverts)) {
              //             ngbes_lcl = int32(zeros(128,1));
              nelems_lcl = elem_one_ring(ngbvs[ii - 1], tets, sibhfs, v2hf,
                etags_elem, ngbes_lcl);
              if (overflow) {
                exitg2 = true;
              } else {
                boolean_T exitg3;
                jj = 0;
                exitg3 = false;
                while ((!exitg3) && (jj <= nelems_lcl - 1)) {
                  eid = ngbes_lcl[jj] - 1;
                  etags_elem[ngbes_lcl[jj] - 1] = false;
                  if ((!etags[ngbes_lcl[jj] - 1]) && (((vtags[tets[tets.size(1) *
                         (ngbes_lcl[jj] - 1)] - 1] + vtags[tets[tets.size(1) *
                         (ngbes_lcl[jj] - 1) + 1] - 1]) + vtags[tets[tets.size(1)
                        * (ngbes_lcl[jj] - 1) + 2] - 1]) + vtags[tets[tets.size
                       (1) * (ngbes_lcl[jj] - 1) + 3] - 1] == 3)) {
                    if (ntemp + 5 < MAX_temp_verts) {
                      if (!vtags[tets[tets.size(1) * eid] - 1]) {
                        ntemp++;
                        temp_verts[ntemp] = tets[tets.size(1) * eid];
                      }

                      if (!vtags[tets[tets.size(1) * eid + 1] - 1]) {
                        ntemp++;
                        temp_verts[ntemp] = tets[tets.size(1) * eid + 1];
                      }

                      if (!vtags[tets[tets.size(1) * eid + 2] - 1]) {
                        ntemp++;
                        temp_verts[ntemp] = tets[tets.size(1) * eid + 2];
                      }

                      if (!vtags[tets[tets.size(1) * eid + 3] - 1]) {
                        ntemp++;
                        temp_verts[ntemp] = tets[tets.size(1) * eid + 3];
                      }

                      jj++;
                    } else {
                      overflow = true;
                      m2cPrintf("Overflow for variable temp_verts when obtaining 1/3 ring in obtain_n_third_ring_tet.m \n");
                      fflush(stdout);
                      exitg3 = true;
                    }
                  } else {
                    jj++;
                  }
                }

                ii++;
              }
            }

            for (ii = 0; ii <= ntemp; ii++) {
              if (*nverts >= 2048) {
                m2cPrintf("Overflow in 1/3 ring in obtain_n_third_ring_tet.m \n");
                fflush(stdout);
                overflow = true;
              } else if (!vtags[temp_verts[ii] - 1]) {
                (*nverts)++;
                ngbvs[*nverts - 1] = temp_verts[ii];
                vtags[temp_verts[ii] - 1] = true;
              }
            }

            if (*nverts >= 1) {
              exitg1 = 1;
            } else {
              //  If needs to expand, then undo the last third ring
              *nverts = nverts_last;
              twothird = true;

              // If the 1/3 ring is not big enough, we will try a 2/3 rings to
              guard2 = true;
            }
          } else {
            guard2 = true;
          }

          if (guard2) {
            if ((cur_ring > ring_full) || ((cur_ring == ring_full) && twothird))
            {
              coder::SizeType numv;

              //  Collect two-third-ring
              nverts_last = *nverts;
              MAX_temp_verts = 24 * *nelems;
              temp_verts.set_size(MAX_temp_verts);
              for (i = 0; i < MAX_temp_verts; i++) {
                temp_verts[i] = 0;
              }

              ntemp = -1;

              // Check currently listed neighboring elements
              ii = nelems_pre;
              exitg2 = false;
              while ((!exitg2) && (ii + 1 <= *nelems)) {
                numv = ((vtags[tets[tets.size(1) * ii] - 1] +
                         vtags[tets[tets.size(1) * ii + 1] - 1]) +
                        vtags[tets[tets.size(1) * ii + 2] - 1]) +
                  vtags[tets[tets.size(1) * ii + 3] - 1];
                if ((numv == 2) || (numv == 3)) {
                  if (ntemp + 5 <= MAX_temp_verts) {
                    temp_verts[ntemp + 1] = tets[tets.size(1) * ii];
                    temp_verts[ntemp + 2] = tets[tets.size(1) * ii + 1];
                    temp_verts[ntemp + 3] = tets[tets.size(1) * ii + 2];
                    temp_verts[ntemp + 4] = tets[tets.size(1) * ii + 3];
                    ntemp += 4;
                    ii++;
                  } else {
                    overflow = true;
                    m2cPrintf("Overflow for variable temp_verts when obtaining 2/3 ring in obtain_n_third_ring_tet.m \n");
                    fflush(stdout);
                    exitg2 = true;
                  }
                } else {
                  ii++;
                }
              }

              // Check latest ring of neighboring elements
              for (ii = nverts_pre; ii <= *nverts; ii++) {
                //             ngbes_lcl = int32(zeros(128,1));
                nelems_lcl = elem_one_ring(ngbvs[ii - 1], tets, sibhfs, v2hf,
                  etags_elem, ngbes_lcl);
                jj = 0;
                exitg2 = false;
                while ((!exitg2) && (jj <= nelems_lcl - 1)) {
                  eid = ngbes_lcl[jj] - 1;
                  etags_elem[ngbes_lcl[jj] - 1] = false;
                  if (!etags[ngbes_lcl[jj] - 1]) {
                    boolean_T b_numv_tmp;
                    boolean_T c_numv_tmp;
                    boolean_T d_numv_tmp;
                    boolean_T numv_tmp;
                    numv_tmp = vtags[tets[tets.size(1) * (ngbes_lcl[jj] - 1)] -
                      1];
                    b_numv_tmp = vtags[tets[tets.size(1) * (ngbes_lcl[jj] - 1) +
                      1] - 1];
                    c_numv_tmp = vtags[tets[tets.size(1) * (ngbes_lcl[jj] - 1) +
                      2] - 1];
                    d_numv_tmp = vtags[tets[tets.size(1) * (ngbes_lcl[jj] - 1) +
                      3] - 1];
                    numv = ((numv_tmp + b_numv_tmp) + c_numv_tmp) + d_numv_tmp;
                    if ((numv == 2) || (numv == 3)) {
                      if (ntemp + 5 <= MAX_temp_verts) {
                        if (!numv_tmp) {
                          ntemp++;
                          temp_verts[ntemp] = tets[tets.size(1) * eid];
                        }

                        if (!b_numv_tmp) {
                          ntemp++;
                          temp_verts[ntemp] = tets[tets.size(1) * eid + 1];
                        }

                        if (!c_numv_tmp) {
                          ntemp++;
                          temp_verts[ntemp] = tets[tets.size(1) * eid + 2];
                        }

                        if (!d_numv_tmp) {
                          ntemp++;
                          temp_verts[ntemp] = tets[tets.size(1) * eid + 3];
                        }

                        jj++;
                      } else {
                        overflow = true;
                        m2cPrintf("Overflow for variable temp_verts when obtaining 2/3 ring in obtain_nring_vol_edited.m \n");
                        fflush(stdout);
                        if (nelems_lcl < 1) {
                          nverts_idx_0 = 0;
                        } else {
                          nverts_idx_0 = nelems_lcl;
                        }

                        for (i = 0; i < nverts_idx_0; i++) {
                          etags_elem[ngbes_lcl[i] - 1] = false;
                        }

                        exitg2 = true;
                      }
                    } else {
                      jj++;
                    }
                  } else {
                    jj++;
                  }
                }
              }

              ii = 0;
              exitg2 = false;
              while ((!exitg2) && (ii <= ntemp)) {
                if (*nverts >= 2048) {
                  m2cPrintf("Overflow in 2/3 ring in obtain_n_third_ring_tet.m \n");
                  fflush(stdout);
                  overflow = true;
                  exitg2 = true;
                } else {
                  if (!vtags[temp_verts[ii] - 1]) {
                    (*nverts)++;
                    ngbvs[*nverts - 1] = temp_verts[ii];
                    vtags[temp_verts[ii] - 1] = true;
                  }

                  ii++;
                }
              }

              if (*nverts >= 1) {
                exitg1 = 1;
              } else {
                //  If needs to expand, then undo the last third ring
                *nverts = nverts_last;
                b_guard1 = true;
              }
            } else {
              b_guard1 = true;
            }
          }

          if (b_guard1) {
            coder::SizeType b_nverts_last;
 
            // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            nverts_last = *nverts + 1;
            nelems_pre = *nelems;
            b_nverts_last = *nverts;
            for (ii = nverts_pre; ii <= b_nverts_last; ii++) {
              //         ngbes_lcl = int32(zeros(128,1));
              nelems_lcl = b_append_one_ring(ngbvs[ii - 1], tets, sibhfs, v2hf,
                ngbvs, nverts, vtags, etags_elem, ngbes_lcl);
              jj = 0;
              exitg2 = false;
              while ((!exitg2) && (jj <= nelems_lcl - 1)) {
                etags_elem[ngbes_lcl[jj] - 1] = false;
                if (!etags[ngbes_lcl[jj] - 1]) {
                  if (*nelems < 2048) {
                    etags[ngbes_lcl[jj] - 1] = true;
                    (*nelems)++;
                    ngbes[*nelems - 1] = ngbes_lcl[jj];
                    jj++;
                  } else {
                    overflow = true;
                    m2cPrintf("Overflow in %g ring in obtain_nring_ring_tet.m \n",
                           cur_ring);
                    fflush(stdout);
                    exitg2 = true;
                  }
                } else {
                  jj++;
                }
              }
            }

            cur_ring++;
            if (((*nverts >= 1) && (cur_ring >= ring)) || (nelems_pre == *nelems)
                || overflow) {
              exitg1 = 1;
            } else {
              nverts_pre = nverts_last;
            }
          }
        } while (exitg1 == 0);

        //  Reset flags
        vtags[vid - 1] = false;
        if (*nverts < 1) {
          nverts_idx_0 = 0;
        } else {
          nverts_idx_0 = *nverts;
        }

        for (i = 0; i < nverts_idx_0; i++) {
          vtags[ngbvs[i] - 1] = false;
        }

        if (*nelems < 1) {
          nverts_idx_0 = 0;
        } else {
          nverts_idx_0 = *nelems;
        }

        for (i = 0; i < nverts_idx_0; i++) {
          etags[ngbes[i] - 1] = false;
        }
      }
    }
  }

  void rdi_compute_stencils(int n, const ::coder::array<int, 2U> &conn, const
    RdiParams *params, const ::coder::array<int, 1U> &nrange, ::coder::array<int,
    2U> &stcls)
  {
    double tEnd;
    double tStart;
    coder::SizeType i;
    coder::SizeType maxStcl;

    // Rrdi_compute_stencils - Compute stencils using AHF
    tStart = 0.0;

#ifdef _OPENMP

    tStart = omp_get_wtime();

#endif // _OPENMP

    //  topological dimension
    maxStcl = params->maxStclSize;

    //  maximum stencil per node
    if (params->maxStclSize == 0) {
      coder::SizeType c;
      c = static_cast<int>(std::ceil(params->ring));
      if (c < 5) {
        maxStcl = iv[(c + ((params->dim - 1) << 2)) - 1];
      } else {
        maxStcl = iv[((params->dim - 1) << 2) + 3];
      }
    }

    if (nrange.size(0) != 0) {
      maxStcl <<= 1;
      i = nrange.size(0);
    } else {
      i = n;
    }

    stcls.set_size(i, maxStcl + 1);
    if (params->verbose > 1) {
      if (nrange.size(0) == 0) {
        m2cPrintf(" Compute stencil with %g ring in %dD...\n", params->ring,
               params->dim);
        fflush(stdout);
      } else {
        m2cPrintf(" Update stencil for %d nodes with %g ring in %dD...\n",
               (int)nrange.size(0), params->ring, params->dim);
        fflush(stdout);
      }
    }

    //  build stencils
    if (params->dim == 1) {
      compute_stcl_kernel1(n, conn, params->ring, stcls, nrange);
    } else if (params->dim == 2) {
      compute_stcl_kernel2(n, conn, params->ring, stcls, nrange);
    } else {
      //  3D, only support tet mesh
      if (conn.size(1) != 4) {
        m2cErrMsgIdAndTxt("rdi_compute_stencils:unsupportMesh",
                          "Only linear tet meshes are supported in 3D");
      }

      compute_stcl_kernel_tet(n, conn, params->ring, stcls, nrange);
    }

    // Provides a portable wall clock timing routine.
    tEnd = 0.0;

#ifdef _OPENMP

    tEnd = omp_get_wtime();

#endif // _OPENMP

    if (params->verbose > 1) {
      m2cPrintf(" Stencil computation finished in %gs...\n", tEnd - tStart);
      fflush(stdout);
    }
  }

  void rdi_compute_stencils(int n, const ::coder::array<int, 2U> &conn, const
    RdiParams *params, ::coder::array<int, 2U> &stcls)
  {
    double tEnd;
    double tStart;
    coder::SizeType maxStcl;

    // Rrdi_compute_stencils - Compute stencils using AHF
    tStart = 0.0;

#ifdef _OPENMP

    tStart = omp_get_wtime();

#endif // _OPENMP

    //  topological dimension
    maxStcl = params->maxStclSize;

    //  maximum stencil per node
    if (params->maxStclSize == 0) {
      coder::SizeType c;
      c = static_cast<int>(std::ceil(params->ring));
      if (c < 5) {
        maxStcl = iv[(c + ((params->dim - 1) << 2)) - 1];
      } else {
        maxStcl = iv[((params->dim - 1) << 2) + 3];
      }
    }

    stcls.set_size(n, maxStcl + 1);
    if (params->verbose > 1) {
      m2cPrintf(" Compute stencil with %g ring in %dD...\n", params->ring,
             params->dim);
      fflush(stdout);
    }

    //  build stencils
    if (params->dim == 1) {
      compute_stcl_kernel1(n, conn, params->ring, stcls);
    } else if (params->dim == 2) {
      compute_stcl_kernel2(n, conn, params->ring, stcls);
    } else {
      //  3D, only support tet mesh
      if (conn.size(1) != 4) {
        m2cErrMsgIdAndTxt("rdi_compute_stencils:unsupportMesh",
                          "Only linear tet meshes are supported in 3D");
      }

      compute_stcl_kernel_tet(n, conn, params->ring, stcls);
    }

    // Provides a portable wall clock timing routine.
    tEnd = 0.0;

#ifdef _OPENMP

    tEnd = omp_get_wtime();

#endif // _OPENMP

    if (params->verbose > 1) {
      m2cPrintf(" Stencil computation finished in %gs...\n", tEnd - tStart);
      fflush(stdout);
    }
  }
}

// End of code generation (rdi_compute_stencils.cpp)
