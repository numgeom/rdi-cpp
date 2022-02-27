# Welcome to the _RDI_ Project #

_RDI_ is in alpha stage..

## Introduction ##

Robust Discontinuity Indicators (_RDI_) is a novel technique for detecting discontinuities of functions on a mesh (or grid). It is developed using a weighted-least-squares (WLS) reconstruction. Besides detection for discontinuities, RDI outputs useful numerical values, which can be used to compute weights in WENO-like methods (or WLS-ENO-based techniques in particular). So far, the RDI technique has been successfully applied to (1) seeking high-order numerical solutions from FVM-discretized hyperbolic conservation laws on unstructured meshes and (2) high-order and conservative data remapping across unmatched interface meshes in multiphysics coupling.

## Installation ##

## Third-Party Dependency: WLSLIB ##

RDI requires _wlslib-cpp_ for its headers during compilation. You can obtain _wlslib-cpp_ freely at <https://github.com/numgeom/wlslib-cpp>.

## Copyright and Licenses ##

The _RDI_ software is developed by the NumGeom Research Group at Stony Brook University, and it is released under BSD-3 license.

## How to Cite _RDI_ ##

If you use _RDI_ in your work for solving numerical PDEs, then please cite the following paper:

```bibtex
@Article{liu2016wls,
  author  = {Hongxu Liu and Xiangmin Jiao},
  journal = {J. Comput. Phys.},
  title   = {{WLS-ENO}: Weighted-least-squares based essentially non-oscillatory schemes for finite volume methods on unstructured meshes},
  year    = {2016},
  pages   = {749--773},
  volume  = {314},
  doi     = {10.1016/j.jcp.2016.03.039},
}
```

If you use _RDI_ for data remapping (or, in general, function reconstructions), then please cite the following paper:

```bibtex
@Article{li2020wls,
  author  = {Yipeng Li and Qiao Chen and Xuebin Wang and Xiangmin Jiao},
  journal = {J. Comput. Phys.},
  title   = {{WLS-ENO} Remap: Superconvergent and non-oscillatory weighted least squares data transfer on surfaces},
  year    = {2020},
  pages   = {109578},
  volume  = {417},
  doi     = {10.1016/j.jcp.2020.109578},
}
```

If you use _RDI_ in other fields, please consider citing both papers. In addition, the core of the RDI technique is based on the _WALF_ reconstruction, which was developed in:

```bibtex
@Article{jiao2012reconstructing,
  author  = {Xiangmin Jiao and Duo Wang},
  journal = {Engrg. Comput.},
  title   = {Reconstructing high-order surfaces for meshing},
  year    = {2012},
  number  = {4},
  pages   = {361--373},
  volume  = {28},
  doi     = {10.1007/s00366-011-0244-8},
}
```

