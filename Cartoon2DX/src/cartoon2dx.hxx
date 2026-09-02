// cartoon2dx.hxx -- driver-independent math for the Cartoon2DX thorn.
//
// This header contains the pieces of the "cartoon" axisymmetry trick
// (see /home/sbrandt/Cactus/cartoon2dx.md and notes.md) that are pure
// math with no dependency on CarpetX's grid/loop machinery: the
// rho/rotation geometry, the tensor-rotation algebra (ported from the
// original Cartoon2D thorn's Cartoon2DBC.c), and 1D Lagrange
// interpolation (ported from its interpolate.c). Kept dependency-free so
// it can be exercised by a standalone host-compiled test
// (test/standalone_math_test.cxx) without building all of Cactus/AMReX.
//
// What this header deliberately does NOT do yet (see notes.md
// "Interpolation design" and the CarpetX-exploration notes for why):
//   - gather a stencil of samples from an actual CarpetX grid function
//     (that needs the tiled/box-based GF accessors and ghost-point loops,
//     and is the next piece of the port);
//   - handle the near-axis case where invert_rho_to_x() would need to
//     take the square root of a negative number (rho < |y_ref|). Per
//     project convention (CLAUDE.md: never silently ignore a problem --
///    make it fatal instead), this is a loud, fatal error for now rather
//     than a silent guess.
//   - ENO interpolation, or excision-aware interpolation. Both existed in
//     the original Cartoon2D thorn but are not needed for the vacuum
//     Brill-Lindquist head-on-collision test that motivates this port.
//     Revisit if a future use case needs them.

#ifndef CARTOON2DX_HXX
#define CARTOON2DX_HXX

#include <cmath>

// This header is used both inside Cactus (where CCTK_REAL/CCTK_ERROR/
// CCTK_DEVICE/CCTK_HOST are already defined via cctk.h/loop.hxx) and by
// the standalone test (which has no Cactus headers at all). Fall back to
// plain `double`, a thrown exception, and no-op device markers in the
// standalone case.
#ifndef CCTK_REAL
#define CARTOON2DX_STANDALONE
using CCTK_REAL = double;
#endif

#ifdef CARTOON2DX_STANDALONE
#include <stdexcept>
#define CARTOON2DX_FATAL(msg) (throw std::domain_error(msg))
#else
#define CARTOON2DX_FATAL(msg) CCTK_ERROR(msg)
#endif

// These functions are called from inside CarpetX's CCTK_DEVICE loop
// lambdas (see cartoon2dx_carpetx.hxx), so they need to be compilable for
// the device too. CCTK_DEVICE/CCTK_HOST are already defined by the time
// cctk.h has been included; the standalone test never defines them, so
// fall back to nothing.
#ifndef CCTK_DEVICE
#define CCTK_DEVICE
#endif
#ifndef CCTK_HOST
#define CCTK_HOST
#endif

namespace Cartoon2DX {

// The three tensor types the original Cartoon2D thorn supports (see its
// Cartoon2D_tensors.h): a plain scalar, a vector with one upper index
// ("u"), and a symmetric rank-2 tensor with two lower indices ("ddsym").
// ADM/BSSN-style variables (lapse, shift, metric, extrinsic curvature)
// are respectively scalar/u/ddsym.
enum class tensortype_t { scalar, u, ddsym };

// ---------------------------------------------------------------------
// Geometry: rho and the relative rotation between two axisymmetric points
// ---------------------------------------------------------------------

inline CCTK_DEVICE CCTK_HOST CCTK_REAL rho_of(CCTK_REAL x, CCTK_REAL y) { return std::sqrt(x * x + y * y); }

// The rotation R(delta) = [[c,-s],[s,c]] that carries tensor components
// sampled at a reference point (x_ref, y_ref) to their value at a target
// point (x_target, y_target), under the assumption that both points lie
// on the same circle about the z-axis (rho_of(x_ref,y_ref) ==
// rho_of(x_target,y_target)). Represented as (c,s) = (cos delta, sin
// delta) rather than the angle itself, computed from coordinate ratios
// rather than atan2/cos/sin -- this mirrors the original Cartoon2D
// thorn's approach (Cartoon2DBC.c), which favored exact ratios over
// transcendental calls, and generalizes it: the original thorn always
// had y_ref == 0 (its reference plane was exactly y=0), so its "delta"
// was simply the target's own angle. Cartoon2DX's vertex-centered
// reference layers sit at y_ref = +-dy/2, not 0, so the angle of the
// reference point itself must be subtracted off.
struct Rotation2D {
  CCTK_REAL c, s; // cos(delta), sin(delta)
};

inline CCTK_DEVICE CCTK_HOST Rotation2D relative_rotation(CCTK_REAL x_ref, CCTK_REAL y_ref,
                                     CCTK_REAL x_target, CCTK_REAL y_target) {
  const CCTK_REAL rho_ref = rho_of(x_ref, y_ref);
  const CCTK_REAL rho_target = rho_of(x_target, y_target);
  const CCTK_REAL cos_t = x_target / rho_target, sin_t = y_target / rho_target;
  const CCTK_REAL cos_r = x_ref / rho_ref, sin_r = y_ref / rho_ref;
  // cos(delta) = cos(phi_t - phi_r), sin(delta) = sin(phi_t - phi_r)
  return Rotation2D{cos_t * cos_r + sin_t * sin_r, sin_t * cos_r - cos_t * sin_r};
}

inline CCTK_DEVICE CCTK_HOST Rotation2D inverse(Rotation2D delta) { return Rotation2D{delta.c, -delta.s}; }

// Solve for the x-coordinate on a reference layer at fixed y = y_ref whose
// distance from the z-axis equals a target rho. Feeding this x back into a
// uniform-grid-in-x interpolator is how Cartoon2DX interpolates "in rho"
// without needing a non-uniform-grid interpolator (see notes.md
// "Interpolation design" for the reasoning and the alternative that was
// rejected). For y_ref == 0 (the cell-centered case) this is just rho
// itself, reproducing the original Cartoon2D thorn's behavior exactly.
//
// Fatal for rho < |y_ref|: the target point would be closer to the
// z-axis than the reference layer, which requires axis-reflection
// handling analogous to the original thorn's x<0 "zombie" case. Not yet
// implemented -- see the file-level comment above.
//
// This branch is called from CarpetX device loop code (see
// cartoon2dx_carpetx.hxx), so CARTOON2DX_FATAL must be device-safe if it
// can actually trigger there. It provably cannot: cartoon2dx_carpetx.hxx
// only ever calls this with y_ref = the canonical INTERIOR reference
// layer and (x,y) = an actual y-GHOST point (gated by is_y_ghost()). A
// genuine y-ghost point is by construction strictly farther from y=0 than
// the interior reference layer, so rho = rho_of(x,y) >= |y| > |y_ref|
// always holds there, and disc = rho^2 - y_ref^2 is always positive. This
// branch is only reachable from the standalone/general-purpose use of
// this function (e.g. the unit test), where CCTK_ERROR's device-safety is
// moot.
inline CCTK_DEVICE CCTK_HOST CCTK_REAL invert_rho_to_x(CCTK_REAL rho, CCTK_REAL y_ref) {
  const CCTK_REAL disc = rho * rho - y_ref * y_ref;
  if (disc < 0)
    CARTOON2DX_FATAL(
        "Cartoon2DX::invert_rho_to_x: rho < |y_ref| (near-axis case not yet implemented)");
  return std::sqrt(disc);
}

// Same inversion, but with the sign of `x_same_sign`. Ghost-fill must interpolate
// on the same side of the axis as the target point: the unsigned +sqrt result
// would send every x<0 query to x=+rho, which is typically a different AMReX
// box (and a different MPI rank). Axisymmetric scalars are even in x, so that
// jump is only correct if the +rho samples are actually readable; they are not
// in a local-box interpolator. Using the same-sign source keeps the stencil
// within ~|y_ghost| of the target, i.e. local.
inline CCTK_DEVICE CCTK_HOST CCTK_REAL invert_rho_to_x_signed(CCTK_REAL rho, CCTK_REAL y_ref,
                                                              CCTK_REAL x_same_sign) {
  return std::copysign(invert_rho_to_x(rho, y_ref), x_same_sign);
}

// ---------------------------------------------------------------------
// Tensor rotation algebra (ported from Cartoon2DBC.c's BndCartoon2DVI)
// ---------------------------------------------------------------------

// Rotates a one-upper-index vector's in-plane (x,y) components by delta.
// Applies to: the full tensortype_t::u case, and (see rotate_ddsym below)
// the (xz,yz) sub-block of a tensortype_t::ddsym tensor. The z-component
// of a "u" tensor is invariant and simply passed through by the caller.
inline CCTK_DEVICE CCTK_HOST void rotate_u(Rotation2D delta, CCTK_REAL fx, CCTK_REAL fy,
                     CCTK_REAL &tx, CCTK_REAL &ty) {
  tx = delta.c * fx - delta.s * fy;
  ty = delta.s * fx + delta.c * fy;
}

// Rotates the in-plane (xx,xy,yy) block of a symmetric two-lower-index
// tensor by delta. Two-lower-index (covariant) components transform with
// the *inverse* rotation relative to one-upper-index components -- this
// is exactly why the original thorn's DDSYM case applied its S matrix
// where its U case applied R, and vice versa for the opposite-sign ghost
// copy (compare the +dj/-dj blocks in Cartoon2DBC.c).
inline CCTK_DEVICE CCTK_HOST void rotate_ddsym_inplane(Rotation2D delta, CCTK_REAL fxx, CCTK_REAL fxy,
                                 CCTK_REAL fyy, CCTK_REAL &txx, CCTK_REAL &txy,
                                 CCTK_REAL &tyy) {
  const Rotation2D r = inverse(delta);
  const CCTK_REAL mxx = r.c, mxy = r.s, myx = -r.s, myy = r.c; // S(delta) = R(-delta)
  txx = fxx * mxx * mxx + 2 * fxy * mxx * myx + fyy * myx * myx;
  tyy = fxx * mxy * mxy + 2 * fxy * mxy * myy + fyy * myy * myy;
  txy = fxx * mxx * mxy + fxy * (mxy * myx + mxx * myy) + fyy * myx * myy;
}

// Convenience wrapper rotating all six components of a symmetric
// two-lower-index tensor (xx,xy,xz,yy,yz,zz) by delta. The (xz,yz)
// sub-block transforms like a one-upper-index vector (see rotate_u); zz
// is invariant.
inline CCTK_DEVICE CCTK_HOST void rotate_ddsym(Rotation2D delta, CCTK_REAL fxx, CCTK_REAL fxy, CCTK_REAL fxz,
                         CCTK_REAL fyy, CCTK_REAL fyz, CCTK_REAL fzz, CCTK_REAL &txx,
                         CCTK_REAL &txy, CCTK_REAL &txz, CCTK_REAL &tyy, CCTK_REAL &tyz,
                         CCTK_REAL &tzz) {
  rotate_ddsym_inplane(delta, fxx, fxy, fyy, txx, txy, tyy);
  rotate_u(delta, fxz, fyz, txz, tyz);
  tzz = fzz;
}

// ---------------------------------------------------------------------
// 1D Lagrange interpolation on a uniform grid
// (ported from interpolate.c's interpolate_local, generalized: the
// original had 5 separate hand-unrolled order1..order5 polynomials
// generated via Maple; a single loop-based evaluator covers any order
// without the duplication -- see notes.md.)
// ---------------------------------------------------------------------

// Interpolates the degree-`order` polynomial through the order+1
// uniformly spaced points (x0 + i*dx, y[i]), i = 0..order, evaluated at
// x. No excision support (see file-level comment).
inline CCTK_DEVICE CCTK_HOST CCTK_REAL interpolate_lagrange(int order, CCTK_REAL x0, CCTK_REAL dx,
                                      const CCTK_REAL *y, CCTK_REAL x) {
  CCTK_REAL result = 0;
  for (int i = 0; i <= order; ++i) {
    CCTK_REAL term = y[i];
    const CCTK_REAL xi = x0 + i * dx;
    for (int j = 0; j <= order; ++j) {
      if (j == i)
        continue;
      const CCTK_REAL xj = x0 + j * dx;
      term *= (x - xj) / (xi - xj);
    }
    result += term;
  }
  return result;
}

} // namespace Cartoon2DX

#endif // CARTOON2DX_HXX
