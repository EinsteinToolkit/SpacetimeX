// cartoon2dx_carpetx.hxx -- CarpetX-specific glue for the ghost-fill.
//
// Builds on cartoon2dx.hxx's driver-independent math (rho/rotation
// geometry, tensor rotation, 1D Lagrange interpolation) to actually fill
// the cartoon points of a CarpetX grid function.
//
// Mixed centering: CI/CJ/CK are 0 for vertex and 1 for cell on each
// axis (CarpetX CENTERING={vvv}/{ccc}/{vcc}/...). The loop and the
// GF3D5layout use those same flags. rho is always rho_of(p.x, p.y):
// PointDesc's .x/.y are already the centering-aware coordinates, so a
// cell-centered hydro point at the thin-slab interior has y=0 and a
// vertex-centered metric point has y=+-dy/2.
//
// Reference layer: j_ref = nghostzones[1] is the first interior y
// index for both 2-vertex (VVV) and 1-cell (CCC) slabs.

#ifndef CARTOON2DX_CARPETX_HXX
#define CARTOON2DX_CARPETX_HXX

#include "cartoon2dx.hxx"

#include <cctk.h>
#include <loop_device.hxx>

#include <cmath>

namespace Cartoon2DX {

using Loop::GF3D5;
using Loop::GF3D5index;
using Loop::GF3D5layout;
using Loop::PointDesc;

inline int reference_j(const Arith::vect<int, Loop::dim> &nghostzones) {
  return nghostzones[1];
}

// Local shape of a GF with this centering: one fewer point on each
// cell-centered axis than the vertex lsh stored on GridDescBase.
inline Arith::vect<int, Loop::dim>
centering_lsh(const Arith::vect<int, Loop::dim> &vertex_lsh, int CI, int CJ,
              int CK) {
  return {vertex_lsh[0] - CI, vertex_lsh[1] - CJ, vertex_lsh[2] - CK};
}

// True iff p is in the y-direction ghost zone of a GF whose local
// y-size is lsh[1] (already centering-adjusted).
inline bool is_y_ghost(const PointDesc &p,
                       const Arith::vect<int, Loop::dim> &nghostzones,
                       const Arith::vect<int, Loop::dim> &lsh) {
  return p.I[1] < nghostzones[1] || p.I[1] >= lsh[1] - nghostzones[1];
}

// y-ghosts, and (when fill_negative_x) the x<0 hang-over used as a
// derivative stencil at the axis. Upper-x and both z faces stay with
// CarpetX / NewRadX. p.x is the centering-aware coordinate.
inline bool needs_cartoon_fill(const PointDesc &p,
                               const Arith::vect<int, Loop::dim> &nghostzones,
                               const Arith::vect<int, Loop::dim> &lsh,
                               bool fill_negative_x) {
  if (is_y_ghost(p, nghostzones, lsh))
    return true;
  return fill_negative_x && p.x < 0;
}

// Gathers an (order+1)-point stencil along x at the fixed reference
// y-layer (same k as p) and interpolates it at x_src on that layer.
// Canonical data lives at x>=0: when fill_negative_x, invert unsigned
// (always +sqrt) and clamp the stencil out of the x<0 ghosts we are
// filling. Otherwise keep the same-sign invert so a full-x domain can
// interpolate x<0 locally.
//
// y_ref is reconstructed from this point's own centering-aware p.y,
// so cell-centered interiors (y=0) and vertex interiors (y=+-dy/2)
// both work.
inline CCTK_REAL interpolate_reference(const GF3D5<CCTK_REAL> &gf,
                                       const GF3D5layout &layout,
                                       const PointDesc &p, int j_ref, int order,
                                       int lsh_x, int nghost_x,
                                       bool fill_negative_x, CCTK_REAL rho) {
  const CCTK_REAL y_ref = p.y - (p.I[1] - j_ref) * p.DX[1];
  const CCTK_REAL x_src = fill_negative_x
                              ? invert_rho_to_x(rho, y_ref)
                              : invert_rho_to_x_signed(rho, y_ref, p.x);

  // Fractional index of x_src along the reference layer, referenced from
  // this point's own (index, coordinate) pair -- valid since the grid
  // spacing is uniform, and avoids needing the grid's origin separately.
  const CCTK_REAL i_real = p.I[0] + (x_src - p.X[0]) / p.DX[0];
  const int i_center = int(std::lround(i_real));
  const int i_lo = fill_negative_x ? nghost_x : 0;

  int offset = order / 2;
  if (i_center - offset < i_lo)
    offset = i_center - i_lo;
  if (i_center - offset + order >= lsh_x)
    offset = i_center + order - lsh_x + 1;
  if (i_center - offset < i_lo)
    offset = i_center - i_lo;
  const int i0 = i_center - offset;

  CCTK_REAL y[6]; // order <= 5 (see param.ccl), so at most 6 points
  for (int n = 0; n <= order; ++n) {
    const GF3D5index idx(layout,
                         Arith::vect<int, Loop::dim>{i0 + n, j_ref, p.I[2]});
    y[n] = gf(idx);
  }
  const CCTK_REAL x0 = p.X[0] + (i0 - p.I[0]) * p.DX[0];
  return interpolate_lagrange(order, x0, p.DX[0], y, x_src);
}

template <int CI, int CJ, int CK, typename Grid>
void fill_ghosts_scalar(const cGH *restrict cctkGH, const Grid &grid, int order,
                        bool fill_negative_x, const GF3D5<CCTK_REAL> &gf) {
  static_assert(CI == 0 || CI == 1);
  static_assert(CJ == 0 || CJ == 1);
  static_assert(CK == 0 || CK == 1);
  const GF3D5layout layout(cctkGH, {CI, CJ, CK});
  const int j_ref = reference_j(grid.nghostzones);
  const auto lsh = centering_lsh(grid.lsh, CI, CJ, CK);
  const int lsh_x = lsh[0];
  const int nghost_x = grid.nghostzones[0];
  const auto nghostzones = grid.nghostzones;
  grid.template loop_all_device<CI, CJ, CK>(
      grid.nghostzones,
      [=] CCTK_DEVICE(const PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
        if (!needs_cartoon_fill(p, nghostzones, lsh, fill_negative_x))
          return;
        const CCTK_REAL rho = rho_of(p.x, p.y);
        const GF3D5index idx(layout, p.I);
        gf.store(idx, interpolate_reference(gf, layout, p, j_ref, order, lsh_x,
                                            nghost_x, fill_negative_x, rho));
      });
}

template <int CI, int CJ, int CK, typename Grid>
void fill_ghosts_vector(const cGH *restrict cctkGH, const Grid &grid, int order,
                        bool fill_negative_x, const GF3D5<CCTK_REAL> &gfx,
                        const GF3D5<CCTK_REAL> &gfy,
                        const GF3D5<CCTK_REAL> &gfz) {
  static_assert(CI == 0 || CI == 1);
  static_assert(CJ == 0 || CJ == 1);
  static_assert(CK == 0 || CK == 1);
  const GF3D5layout layout(cctkGH, {CI, CJ, CK});
  const int j_ref = reference_j(grid.nghostzones);
  const auto lsh = centering_lsh(grid.lsh, CI, CJ, CK);
  const int lsh_x = lsh[0];
  const int nghost_x = grid.nghostzones[0];
  const auto nghostzones = grid.nghostzones;
  grid.template loop_all_device<CI, CJ, CK>(
      grid.nghostzones,
      [=] CCTK_DEVICE(const PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
        if (!needs_cartoon_fill(p, nghostzones, lsh, fill_negative_x))
          return;
        const CCTK_REAL rho = rho_of(p.x, p.y);
        const CCTK_REAL y_ref = p.y - (p.I[1] - j_ref) * p.DX[1];
        const CCTK_REAL x_src = fill_negative_x
                                    ? invert_rho_to_x(rho, y_ref)
                                    : invert_rho_to_x_signed(rho, y_ref, p.x);
        const Rotation2D delta = relative_rotation(x_src, y_ref, p.x, p.y);

        const CCTK_REAL fx = interpolate_reference(
            gfx, layout, p, j_ref, order, lsh_x, nghost_x, fill_negative_x, rho);
        const CCTK_REAL fy = interpolate_reference(
            gfy, layout, p, j_ref, order, lsh_x, nghost_x, fill_negative_x, rho);
        const CCTK_REAL fz = interpolate_reference(
            gfz, layout, p, j_ref, order, lsh_x, nghost_x, fill_negative_x, rho);
        CCTK_REAL tx, ty;
        rotate_u(delta, fx, fy, tx, ty);

        const GF3D5index idx(layout, p.I);
        gfx.store(idx, tx);
        gfy.store(idx, ty);
        gfz.store(idx, fz);
      });
}

template <int CI, int CJ, int CK, typename Grid>
void fill_ghosts_ddsym(const cGH *restrict cctkGH, const Grid &grid, int order,
                       bool fill_negative_x, const GF3D5<CCTK_REAL> &gfxx,
                       const GF3D5<CCTK_REAL> &gfxy,
                       const GF3D5<CCTK_REAL> &gfxz,
                       const GF3D5<CCTK_REAL> &gfyy,
                       const GF3D5<CCTK_REAL> &gfyz,
                       const GF3D5<CCTK_REAL> &gfzz) {
  static_assert(CI == 0 || CI == 1);
  static_assert(CJ == 0 || CJ == 1);
  static_assert(CK == 0 || CK == 1);
  const GF3D5layout layout(cctkGH, {CI, CJ, CK});
  const int j_ref = reference_j(grid.nghostzones);
  const auto lsh = centering_lsh(grid.lsh, CI, CJ, CK);
  const int lsh_x = lsh[0];
  const int nghost_x = grid.nghostzones[0];
  const auto nghostzones = grid.nghostzones;
  grid.template loop_all_device<CI, CJ, CK>(
      grid.nghostzones,
      [=] CCTK_DEVICE(const PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
        if (!needs_cartoon_fill(p, nghostzones, lsh, fill_negative_x))
          return;
        const CCTK_REAL rho = rho_of(p.x, p.y);
        const CCTK_REAL y_ref = p.y - (p.I[1] - j_ref) * p.DX[1];
        const CCTK_REAL x_src = fill_negative_x
                                    ? invert_rho_to_x(rho, y_ref)
                                    : invert_rho_to_x_signed(rho, y_ref, p.x);
        const Rotation2D delta = relative_rotation(x_src, y_ref, p.x, p.y);

        const CCTK_REAL fxx = interpolate_reference(
            gfxx, layout, p, j_ref, order, lsh_x, nghost_x, fill_negative_x, rho);
        const CCTK_REAL fxy = interpolate_reference(
            gfxy, layout, p, j_ref, order, lsh_x, nghost_x, fill_negative_x, rho);
        const CCTK_REAL fxz = interpolate_reference(
            gfxz, layout, p, j_ref, order, lsh_x, nghost_x, fill_negative_x, rho);
        const CCTK_REAL fyy = interpolate_reference(
            gfyy, layout, p, j_ref, order, lsh_x, nghost_x, fill_negative_x, rho);
        const CCTK_REAL fyz = interpolate_reference(
            gfyz, layout, p, j_ref, order, lsh_x, nghost_x, fill_negative_x, rho);
        const CCTK_REAL fzz = interpolate_reference(
            gfzz, layout, p, j_ref, order, lsh_x, nghost_x, fill_negative_x, rho);
        CCTK_REAL txx, txy, txz, tyy, tyz, tzz;
        rotate_ddsym(delta, fxx, fxy, fxz, fyy, fyz, fzz, txx, txy, txz, tyy,
                     tyz, tzz);

        const GF3D5index idx(layout, p.I);
        gfxx.store(idx, txx);
        gfxy.store(idx, txy);
        gfxz.store(idx, txz);
        gfyy.store(idx, tyy);
        gfyz.store(idx, tyz);
        gfzz.store(idx, tzz);
      });
}

} // namespace Cartoon2DX

#endif // CARTOON2DX_CARPETX_HXX
