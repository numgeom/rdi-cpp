// Copyright 2022 The NumGeom Group, Stony Brook University
//
// rdi_params.hpp
//

#ifndef RDI_PARAMS_H
#define RDI_PARAMS_H

#ifndef MATLAB_MEX_FILE
// definition of parameter structure
namespace rdi{
/*!
 * @struct RdiParams
 * @brief Parameter structure used in RDI library
 *
 */
#endif // MATLAB_MEX_FILE

struct RdiParams {
  int    dim;          ///< topological dimension
  double ring;         ///< ring size (2)
  int    maxStclSize;  ///< max stencil size per node (1-3D: 10,30,60)
  int    verbose;      ///< verbose (<1 silent, 1 default, >1 detailed)
  int    degree;       ///< monomial degree (2)
  int surfType;        ///< special surface type (0); surfType=1 means spherical (3D)
                       ///< or circle (2D) surface
  double epsBeta;      ///< threshold for computing @f$\beta@f$
  double hGlobal;      ///< global mesh size (0)
  double cGlobal;      ///< global threshold for OSUS detection (0.05)
  double cLocal;       ///< local threshold for OSUS detection (0.5)
  double kappa1;       ///< threshold for osc. detection for @f$C_1$@f disc. (0.3)
  double kappa0;       ///< threshold for osc. detection for @f$C_0$@f disc. (1.0)
  bool   markNearDis;  ///< mark near-by (1-ring) nodes as dis. (true)
  bool   wlsInterp0;   ///< use interpolation mode for WLS (false)
  bool   wlsUseDag;    ///< use dag for WLS (false)
  int    nThreads;     ///< number of parallel threads
  bool   parTask;      ///< use OpenMP task for building OSUS operators
};

#ifndef MATLAB_MEX_FILE

}   // namespace rdi

namespace rdi_stencils {
using rdi::RdiParams;
}

namespace rdi_kernel {
using rdi::RdiParams;
}
#endif

#endif // RDI_PARAMS_H
