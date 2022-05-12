# Matriplex

Matriplex is a C++ header-only template library implementing data structures
and methods to enhance the efficiency of matrix operations involving sets of
identically-sized small matrices. Such matrices are typically too small to
benefit individually from SIMD processing, but when taken as a collection,
they are able to do so.

Matriplex's memory layout uses a matrix-major representation optimized for
loading vector registers for SIMD operations on a set of small matrices, using
the native vector-unit width on processors with vector units. Matriplex
includes a code generator for defining optimized matrix operations, with
support for symmetric matrices and on-the-fly matrix transposition. Patterns
of elements that are known by construction to be zero or one can be specified,
and the resulting code will be optimized to eliminate unnecessary register
loads and arithmetic operations. The generated code can be either standard C++
or macros that map to architecture-specific intrinsic functions.

Matriplex data structures and auto-generated code are at the core of the
[mkFit](https://iris-hep.org/projects/mkfit.html) package for the vectorized,
parallelized reconstruction of particle tracks in multi-layered detectors. The
mkFit algorithm is based on Kalman Filter (KF) techniques, and within mkFit,
Matriplex is used for all KF-related operations on tracks and hits, as well as
matrix operations in general. In 2022, the mkFit package was integrated into
the CMSSW analysis toolkit for the CMS detector in the LHC at CERN.
