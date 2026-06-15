# Updates from CMSSW Integration

This branch brings the standalone Matriplex library up to parity with the version in
`RecoTracker/MkFitCore/src/Matriplex` in CMSSW. The changes fall into four files.

---

## `MatriplexCommon.h`

- Added `MPLEX_ALIGN` macro that selects the alignment boundary based on the detected ISA
  at compile time: 64 bytes for AVX-512, 32 bytes for AVX/AVX2, 16 bytes for SSE3.
- Added `internal::sincos4()`: a fast 4th-order Taylor series approximation for sin and cos,
  used in track propagation.

---

## `Matriplex.h`

### Alignment
The alignment strategy has changed. The class itself now carries
`__attribute__((aligned(MPLEX_ALIGN)))`, and the per-member alignment attribute on `fArray`
has been removed. Alignment is enforced at the object level rather than the array member level.

### Element-wise access helpers
New inner types `QReduced` and `QAssigner` enable extracting or assigning a fixed `(i,j)`
element across all N matrices as a 1×1 Matriplex:

```cpp
m.aij(2, 3) = other.rij(0, 0);   // assign element (0,0) of other into (2,3) of m
```

### Arithmetic operators
Full in-place and free-function arithmetic operator suite:
- In-place: `+=`, `-=`, `*=`, `/=` — for both scalar `T` and `MPlex` operands
- Unary `-`
- Free-function: `+`, `-`, `*`, `/` — for MPlex×MPlex, MPlex×scalar, and scalar×MPlex

### Math methods
Element-wise math methods added as both member functions and free functions at namespace scope:

| Method | Notes |
|---|---|
| `abs()`, `sqr()` | magnitude helpers |
| `sqrt()`, `hypot()` | |
| `sin()`, `cos()`, `tan()`, `atan2()` | standard-library backed |
| `sincos()` | computes sin and cos together |
| `sincos4()` | fast Taylor series via `internal::sincos4()` |
| `negate()`, `negate_if_ltz(sign)` | sign-flip helpers |
| `min()`, `max()`, `min_max()` | element-wise min/max free functions |

### VDT fast-math (`#ifdef MPLEX_VDT`)
When CERN's [VDT library](https://github.com/dpiparo/vdt) is available, approximate but
faster transcendental variants are exposed:
`fast_isqrt()`, `fast_sin()`, `fast_cos()`, `fast_sincos()`, `fast_tan()`, `fast_atan2()`.
These dispatch to `vdt::fast_*f()` for `float` and `vdt::fast_*()` for `double`.

---

## `MatriplexSym.h`

- Same class-level alignment change as `Matriplex.h`.
- Added the D=5 symmetric offset table — previously only D=3 and D=6 were present.
- Added explicit `= default` copy constructor.
- Added `ReduceFixedIJ()` method, mirroring the `QReduced` accessor added to `Matriplex.h`.

---

## `GenMPlexOps.pl`

Added a new **"KalmanOps on Plane"** section that generates Jacobian and error-propagation
operations for a plane-based (local coordinate) Kalman filter variant. This is required to
handle CMS silicon strip and pixel geometry, where detector layers use local plane coordinates
rather than global cylindrical ones. The new generated files are:

| File | Operation |
|---|---|
| `JacCCS2Loc.ah` | Jacobian: Curvilinear → local coordinates |
| `PsErrLoc.ah` / `PsErrLocTransp.ah` | Propagate state error into local frame |
| `PsErrLocUpd.ah` | Apply Kalman gain update to local error |
| `JacLoc2CCS.ah` | Jacobian: local → Curvilinear coordinates |
| `OutErrCCS.ah` / `OutErrCCSTransp.ah` | Propagate updated error back to CCS frame |
| `MultHelixPlaneProp.ah` / `MultHelixPlanePropTransp.ah` | Helix propagation similarity transform |
| `JacErrPropCurv1.ah` / `JacErrPropCurv2.ah` | Error propagation Jacobian in curvature basis |

---

## Subsequent work (not in CMSSW)

`MatriplexSym` was extended to match the capabilities added to `Matriplex` in this update:
arithmetic operators, math methods, VDT fast-math support, `QReduced`/`QAssigner` element
access helpers, `negate()`, and namespace-scope free functions.
