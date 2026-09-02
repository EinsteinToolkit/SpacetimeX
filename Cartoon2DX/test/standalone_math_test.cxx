// standalone_math_test.cxx -- sanity tests for cartoon2dx.hxx.
//
// This is NOT part of the Cactus build (it is intentionally omitted from
// src/make.code.defn's SRCS) and is not a Cactus parameter-file test like
// the other files under test/. It exercises the driver-independent pure
// math in cartoon2dx.hxx directly, without needing to build all of
// Cactus/AMReX/CarpetX first.
//
// Build and run with (single line):
//   g++ -std=c++17 -Wall -Wextra -I../src -o standalone_math_test standalone_math_test.cxx && ./standalone_math_test

#include "cartoon2dx.hxx"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>

using namespace Cartoon2DX;

namespace {

int g_failures = 0;

void expect_near(double actual, double expected, double tol, const char *what) {
  if (std::fabs(actual - expected) > tol) {
    std::printf("FAIL: %s -- expected %.15g, got %.15g (diff %.3g)\n", what, expected,
                actual, actual - expected);
    ++g_failures;
  } else {
    std::printf("PASS: %s\n", what);
  }
}

void expect_true(bool cond, const char *what) {
  if (!cond) {
    std::printf("FAIL: %s\n", what);
    ++g_failures;
  } else {
    std::printf("PASS: %s\n", what);
  }
}

// --- rho / rotation geometry ---------------------------------------------

void test_rho_of() {
  expect_near(rho_of(3, 4), 5, 1e-14, "rho_of(3,4) == 5");
  expect_near(rho_of(-3, 4), 5, 1e-14, "rho_of(-3,4) == 5 (sign-independent)");
  expect_near(rho_of(0, 0.5), 0.5, 1e-14, "rho_of(0,0.5) == 0.5");
}

void test_relative_rotation_matches_old_cartoon2d_when_yref_zero() {
  // When the reference point sits at y_ref=0 (old Cartoon2D's assumption),
  // delta must reduce to exactly the target point's own angle: rotating
  // (rho,0) by delta should land exactly on (x_target,y_target).
  const double x_target = 3, y_target = 4, rho = rho_of(x_target, y_target);
  const auto delta = relative_rotation(rho, 0, x_target, y_target);
  double tx, ty;
  rotate_u(delta, rho, 0, tx, ty);
  expect_near(tx, x_target, 1e-13, "y_ref=0: rotate_u(delta, rho, 0) recovers x_target");
  expect_near(ty, y_target, 1e-13, "y_ref=0: rotate_u(delta, rho, 0) recovers y_target");
}

void test_relative_rotation_offaxis_reference() {
  // General case: reference point off the y=0 plane (the vertex-centered
  // case, y_ref = dy/2), same rho as the target. Rotating the reference
  // point by delta must still land exactly on the target.
  const double y_ref = 0.05; // dy/2
  const double x_target = 3, y_target = 4;
  const double rho = rho_of(x_target, y_target);
  const double x_ref = invert_rho_to_x(rho, y_ref);
  expect_near(rho_of(x_ref, y_ref), rho, 1e-13,
              "invert_rho_to_x round-trips through rho_of");

  const auto delta = relative_rotation(x_ref, y_ref, x_target, y_target);
  double tx, ty;
  rotate_u(delta, x_ref, y_ref, tx, ty);
  expect_near(tx, x_target, 1e-13, "off-axis ref: rotate_u recovers x_target");
  expect_near(ty, y_target, 1e-13, "off-axis ref: rotate_u recovers y_target");
}

void test_invert_rho_to_x_near_axis_is_fatal() {
  bool threw = false;
  try {
    invert_rho_to_x(0.01, 0.05); // rho < y_ref
  } catch (const std::domain_error &) {
    threw = true;
  }
  expect_true(threw, "invert_rho_to_x(rho < y_ref) throws (near-axis case not implemented)");
}

// --- tensor rotation algebra ----------------------------------------------

void test_rotate_u_is_orthogonal_round_trip() {
  const auto delta = relative_rotation(1, 0, 0.6, 0.8); // rho=1 in both cases
  const double fx = 2.0, fy = -3.5;
  double tx, ty, bx, by;
  rotate_u(delta, fx, fy, tx, ty);
  rotate_u(inverse(delta), tx, ty, bx, by);
  expect_near(bx, fx, 1e-13, "rotate_u forward+inverse round-trips fx");
  expect_near(by, fy, 1e-13, "rotate_u forward+inverse round-trips fy");
  // Rotation must preserve length.
  expect_near(tx * tx + ty * ty, fx * fx + fy * fy, 1e-12,
              "rotate_u preserves vector length");
}

void test_rotate_ddsym_inplane_round_trip_and_trace() {
  const auto delta = relative_rotation(1, 0, -0.6, 0.8);
  const double fxx = 1.5, fxy = 0.3, fyy = -2.0;
  double txx, txy, tyy, bxx, bxy, byy;
  rotate_ddsym_inplane(delta, fxx, fxy, fyy, txx, txy, tyy);
  rotate_ddsym_inplane(inverse(delta), txx, txy, tyy, bxx, bxy, byy);
  expect_near(bxx, fxx, 1e-12, "rotate_ddsym_inplane round-trips fxx");
  expect_near(bxy, fxy, 1e-12, "rotate_ddsym_inplane round-trips fxy");
  expect_near(byy, fyy, 1e-12, "rotate_ddsym_inplane round-trips fyy");
  // Trace of a symmetric 2-tensor is a rotation invariant.
  expect_near(txx + tyy, fxx + fyy, 1e-12, "rotate_ddsym_inplane preserves trace");
}

// Cross-check rotate_ddsym_inplane against a naive, independently-written
// congruence transform using explicit trig, to catch sign-convention bugs
// that a mere round-trip test could miss (a round trip is consistent with
// itself even if inverse() had the wrong sign).
void test_rotate_ddsym_inplane_matches_naive_congruence() {
  const double phi_ref = 0.2, phi_target = 1.1;
  const double x_ref = std::cos(phi_ref), y_ref = std::sin(phi_ref);
  const double x_target = std::cos(phi_target), y_target = std::sin(phi_target);
  const auto delta = relative_rotation(x_ref, y_ref, x_target, y_target);

  const double fxx = 0.7, fxy = -0.4, fyy = 1.3;
  double txx, txy, tyy;
  rotate_ddsym_inplane(delta, fxx, fxy, fyy, txx, txy, tyy);

  // Naive ground truth: T' = S^T F S where S = R(-delta_angle),
  // delta_angle = phi_target - phi_ref (matches the derivation in
  // cartoon2dx.hxx's comments), built directly from cos/sin.
  const double da = phi_target - phi_ref;
  const double c = std::cos(-da), s = std::sin(-da); // S = R(-da)
  const double naive_txx = fxx * c * c + 2 * fxy * c * (-s) + fyy * s * s;
  const double naive_tyy = fxx * s * s + 2 * fxy * c * s + fyy * c * c;
  const double naive_txy = fxx * c * s + fxy * (c * c - s * s) - fyy * s * c;

  expect_near(txx, naive_txx, 1e-12, "rotate_ddsym_inplane matches naive trig congruence (xx)");
  expect_near(tyy, naive_tyy, 1e-12, "rotate_ddsym_inplane matches naive trig congruence (yy)");
  expect_near(txy, naive_txy, 1e-12, "rotate_ddsym_inplane matches naive trig congruence (xy)");
}

void test_rotate_ddsym_full_passthrough_and_vector_block() {
  const auto delta = relative_rotation(1, 0, 0, 1); // quarter turn
  double txx, txy, txz, tyy, tyz, tzz;
  rotate_ddsym(delta, /*fxx*/ 1, /*fxy*/ 0, /*fxz*/ 5, /*fyy*/ 1, /*fyz*/ 7, /*fzz*/ 9,
               txx, txy, txz, tyy, tyz, tzz);
  expect_near(tzz, 9, 1e-14, "rotate_ddsym passes zz through unchanged");
  // (xz,yz) should behave exactly like rotate_u on the same delta.
  double vx, vy;
  rotate_u(delta, 5, 7, vx, vy);
  expect_near(txz, vx, 1e-13, "rotate_ddsym's xz block matches rotate_u");
  expect_near(tyz, vy, 1e-13, "rotate_ddsym's yz block matches rotate_u");
}

// --- Lagrange interpolation -----------------------------------------------

void test_interpolate_lagrange_reproduces_polynomial() {
  // A cubic sampled at 4 uniformly spaced points must be reproduced
  // exactly (up to roundoff) by order-3 Lagrange interpolation anywhere,
  // including well outside the sample points.
  auto f = [](double x) { return 2 * x * x * x - 3 * x * x + x - 7; };
  const double x0 = -1.0, dx = 0.5;
  double y[4];
  for (int i = 0; i < 4; ++i)
    y[i] = f(x0 + i * dx);
  for (double x : {-2.0, -1.0, -0.3, 0.0, 0.7, 1.5, 3.0}) {
    char what[128];
    std::snprintf(what, sizeof what, "interpolate_lagrange(order=3) reproduces cubic at x=%g", x);
    expect_near(interpolate_lagrange(3, x0, dx, y, x), f(x), 1e-11, what);
  }
}

void test_interpolate_lagrange_passes_through_samples() {
  const double x0 = 0.0, dx = 0.1;
  const double y[6] = {1.0, 2.0, -1.0, 0.5, 3.0, 2.5}; // order 5, not a polynomial by construction
  for (int i = 0; i <= 5; ++i) {
    char what[64];
    std::snprintf(what, sizeof what, "interpolate_lagrange(order=5) passes through sample %d", i);
    expect_near(interpolate_lagrange(5, x0, dx, y, x0 + i * dx), y[i], 1e-12, what);
  }
}

void test_signed_invert_stays_on_the_same_side_of_the_axis() {
  const double y_ref = 0.05;
  const double x_target = -3, y_target = 0.4;
  const double rho = rho_of(x_target, y_target);
  const double x_ref = invert_rho_to_x_signed(rho, y_ref, x_target);
  expect_true(x_ref < 0, "invert_rho_to_x_signed keeps x<0 on the negative-x side");
  expect_near(rho_of(x_ref, y_ref), rho, 1e-13,
              "signed invert still lands on the same rho-circle");

  const auto delta = relative_rotation(x_ref, y_ref, x_target, y_target);
  double tx, ty;
  rotate_u(delta, x_ref, y_ref, tx, ty);
  expect_near(tx, x_target, 1e-13, "negative-x ref: rotate_u recovers x_target");
  expect_near(ty, y_target, 1e-13, "negative-x ref: rotate_u recovers y_target");
}

// 1D Lagrange along a uniform x-row, with the same offset-clamping as
// cartoon2dx_carpetx.hxx::interpolate_reference. `row` has `nx` samples at
// x = x0 + i*dx, i = 0..nx-1.
double interpolate_row(int order, int nx, double x0, double dx, const double *row, double x) {
  const double i_real = (x - x0) / dx;
  const int i_center = int(std::lround(i_real));
  int offset = order / 2;
  if (i_center - offset < 0)
    offset = i_center;
  if (i_center - offset + order >= nx)
    offset = i_center + order - nx + 1;
  const int i0 = i_center - offset;
  double y[6];
  for (int n = 0; n <= order; ++n)
    y[n] = row[i0 + n];
  return interpolate_lagrange(order, x0 + i0 * dx, dx, y, x);
}

// Vertex-centered thin-y cartoon fill of an axisymmetric scalar Gaussian,
// including x<0, using only a LOCAL x-window around each target (the thing
// a real AMReX box can actually see). The unsigned +sqrt inversion would
// look up +rho, which is outside that window for every x<0 point.
void test_cartoon_fill_axisymmetric_gaussian_local_box() {
  const int ng = 3, nint_y = 2, nx = 81, order = 4;
  const int ny = nint_y + 2 * ng;
  const int j_ref = ng;
  const double dx = 0.05, dy = 0.05;
  const double xmin = -2.0;
  const int box_half = 8; // ~local box: 2*8+1 samples, like a small AMReX box

  auto x_of = [&](int i) { return xmin + i * dx; };
  auto y_of = [&](int j) { return -0.5 * dy + (j - ng) * dy; };
  auto f = [&](double x, double y) { return std::exp(-(x * x + y * y)); };

  double u[ny * nx];
  for (int j = 0; j < ny; ++j)
    for (int i = 0; i < nx; ++i)
      u[j * nx + i] = (j >= ng && j < ng + nint_y) ? f(x_of(i), y_of(j)) : 0.0;

  double max_abs_err = 0;
  int n_checked = 0;
  for (int j = 0; j < ny; ++j) {
    if (j >= ng && j < ng + nint_y)
      continue;
    const double y = y_of(j);
    const double y_ref = y_of(j_ref);
    for (int i = box_half; i < nx - box_half; ++i) {
      const double x = x_of(i);
      const double rho = rho_of(x, y);
      const double x_src = invert_rho_to_x_signed(rho, y_ref, x);

      // Local window around THIS point, not around +rho. For x<0 the
      // unsigned inversion would ask for a sample outside this window.
      const int i_lo = i - box_half;
      const int nwin = 2 * box_half + 1;
      double row[32];
      for (int k = 0; k < nwin; ++k)
        row[k] = u[j_ref * nx + (i_lo + k)];

      const double got = interpolate_row(order, nwin, x_of(i_lo), dx, row, x_src);
      u[j * nx + i] = got;
      const double err = std::fabs(got - f(x, y));
      if (err > max_abs_err)
        max_abs_err = err;
      ++n_checked;
    }
  }
  expect_true(n_checked > 0, "cartoon fill checked at least one ghost point");
  expect_true(max_abs_err < 1e-4, "local-box cartoon fill of exp(-rho^2) matches analytic in y-ghosts");
  if (max_abs_err >= 1e-4)
    std::printf("      max_abs_err = %.3g over %d points\n", max_abs_err, n_checked);

  // Same fill but with the OLD unsigned invert: must fail on this local-box
  // restriction, because x<0 queries jump to +rho, outside the window.
  double max_unsigned_err = 0;
  int n_neg = 0;
  for (int j = 0; j < ny; ++j) {
    if (j >= ng && j < ng + nint_y)
      continue;
    const double y = y_of(j);
    const double y_ref = y_of(j_ref);
    for (int i = box_half; i < nx - box_half; ++i) {
      const double x = x_of(i);
      if (x >= 0)
        continue;
      const double rho = rho_of(x, y);
      const double x_src_unsigned = invert_rho_to_x(rho, y_ref); // always >= 0
      const int i_lo = i - box_half;
      const int nwin = 2 * box_half + 1;
      double row[32];
      for (int k = 0; k < nwin; ++k)
        row[k] = u[j_ref * nx + (i_lo + k)];
      const double got = interpolate_row(order, nwin, x_of(i_lo), dx, row, x_src_unsigned);
      const double err = std::fabs(got - f(x, y));
      if (err > max_unsigned_err)
        max_unsigned_err = err;
      ++n_neg;
    }
  }
  expect_true(n_neg > 0, "unsigned-invert comparison saw x<0 points");
  expect_true(max_unsigned_err > 1e-2,
              "unsigned invert on a local x<0 box is WRONG (would hide behind a global array)");
  std::printf("      unsigned-invert max_abs_err on x<0 = %.3g (signed was %.3g)\n",
              max_unsigned_err, max_abs_err);
}

// Brill-Lindquist conformal-flat metric (the 2BH initial data):
//   psi = 1 + sum m_n / (2 r_n),  g_ij = psi^4 delta_ij.
// Filling y-ghosts by interpolating each component on the reference
// layer then rotate_ddsym must recover the analytic metric, including
// at x<0, on a local x-window. This is the tensor path Cartoon2DXZ4c4m
// actually uses, which the scalar Gaussian test does not exercise.
void test_cartoon_fill_brill_lindquist_metric_local_box() {
  const int ng = 3, nint_y = 2, nx = 81, order = 4;
  const int ny = nint_y + 2 * ng;
  const int j_ref = ng;
  const double dx = 0.05, dy = 0.05;
  const double xmin = -2.0;
  const int box_half = 8;
  const double z = 0.125; // the slice where the t=9 origin discrepancy showed up
  const double m = 0.5, z1 = 1.15, z2 = -1.15;

  auto x_of = [&](int i) { return xmin + i * dx; };
  auto y_of = [&](int j) { return -0.5 * dy + (j - ng) * dy; };
  auto psi = [&](double x, double y) {
    const double r1 = std::sqrt(x * x + y * y + (z - z1) * (z - z1));
    const double r2 = std::sqrt(x * x + y * y + (z - z2) * (z - z2));
    return 1.0 + m / (2.0 * r1) + m / (2.0 * r2);
  };
  auto gdiag = [&](double x, double y) {
    const double p = psi(x, y);
    return p * p * p * p;
  };

  double gxx[ny * nx], gxy[ny * nx], gxz[ny * nx], gyy[ny * nx], gyz[ny * nx], gzz[ny * nx];
  for (int j = 0; j < ny; ++j)
    for (int i = 0; i < nx; ++i) {
      const int idx = j * nx + i;
      if (j >= ng && j < ng + nint_y) {
        const double gd = gdiag(x_of(i), y_of(j));
        gxx[idx] = gd;
        gyy[idx] = gd;
        gzz[idx] = gd;
        gxy[idx] = gxz[idx] = gyz[idx] = 0.0;
      } else {
        gxx[idx] = gxy[idx] = gxz[idx] = gyy[idx] = gyz[idx] = gzz[idx] = 0.0;
      }
    }

  auto interp_comp = [&](const double *gf, int i, double x_src) {
    const int i_lo = i - box_half;
    const int nwin = 2 * box_half + 1;
    double row[32];
    for (int k = 0; k < nwin; ++k)
      row[k] = gf[j_ref * nx + (i_lo + k)];
    return interpolate_row(order, nwin, x_of(i_lo), dx, row, x_src);
  };

  double max_gxx = 0, max_gxy = 0, max_gyy = 0, max_gzz = 0, max_gxz = 0, max_gyz = 0;
  int n_checked = 0;
  for (int j = 0; j < ny; ++j) {
    if (j >= ng && j < ng + nint_y)
      continue;
    const double y = y_of(j);
    const double y_ref = y_of(j_ref);
    for (int i = box_half; i < nx - box_half; ++i) {
      const double x = x_of(i);
      const double rho = rho_of(x, y);
      const double x_src = invert_rho_to_x_signed(rho, y_ref, x);
      const auto delta = relative_rotation(x_src, y_ref, x, y);
      const double fxx = interp_comp(gxx, i, x_src);
      const double fxy = interp_comp(gxy, i, x_src);
      const double fxz = interp_comp(gxz, i, x_src);
      const double fyy = interp_comp(gyy, i, x_src);
      const double fyz = interp_comp(gyz, i, x_src);
      const double fzz = interp_comp(gzz, i, x_src);
      double txx, txy, txz, tyy, tyz, tzz;
      rotate_ddsym(delta, fxx, fxy, fxz, fyy, fyz, fzz, txx, txy, txz, tyy, tyz, tzz);
      const double gd = gdiag(x, y);
      max_gxx = std::max(max_gxx, std::fabs(txx - gd));
      max_gyy = std::max(max_gyy, std::fabs(tyy - gd));
      max_gzz = std::max(max_gzz, std::fabs(tzz - gd));
      max_gxy = std::max(max_gxy, std::fabs(txy));
      max_gxz = std::max(max_gxz, std::fabs(txz));
      max_gyz = std::max(max_gyz, std::fabs(tyz));
      ++n_checked;
    }
  }
  expect_true(n_checked > 0, "BL metric fill checked at least one ghost point");
  expect_true(max_gxx < 1e-4, "BL cartoon fill recovers gxx = psi^4");
  expect_true(max_gyy < 1e-4, "BL cartoon fill recovers gyy = psi^4");
  expect_true(max_gzz < 1e-4, "BL cartoon fill recovers gzz = psi^4");
  expect_true(max_gxy < 1e-4, "BL cartoon fill keeps gxy ~ 0");
  expect_true(max_gxz < 1e-4, "BL cartoon fill keeps gxz ~ 0");
  expect_true(max_gyz < 1e-4, "BL cartoon fill keeps gyz ~ 0");
  std::printf("      BL metric ghost maxabs err: gxx=%.3g gyy=%.3g gzz=%.3g gxy=%.3g gxz=%.3g gyz=%.3g\n",
              max_gxx, max_gyy, max_gzz, max_gxy, max_gxz, max_gyz);
}

// Radial in-plane vector V = (x, y, 0) * psi^4. Transforms as a U tensor
// under xy rotation; the z component stays 0.
void test_cartoon_fill_radial_vector_local_box() {
  const int ng = 3, nint_y = 2, nx = 81, order = 4;
  const int ny = nint_y + 2 * ng;
  const int j_ref = ng;
  const double dx = 0.05, dy = 0.05;
  const double xmin = -2.0;
  const int box_half = 8;
  const double z = 0.125, m = 0.5, z1 = 1.15, z2 = -1.15;

  auto x_of = [&](int i) { return xmin + i * dx; };
  auto y_of = [&](int j) { return -0.5 * dy + (j - ng) * dy; };
  auto amp = [&](double x, double y) {
    const double r1 = std::sqrt(x * x + y * y + (z - z1) * (z - z1));
    const double r2 = std::sqrt(x * x + y * y + (z - z2) * (z - z2));
    const double p = 1.0 + m / (2.0 * r1) + m / (2.0 * r2);
    return p * p * p * p;
  };

  double vx[ny * nx], vy[ny * nx], vz[ny * nx];
  for (int j = 0; j < ny; ++j)
    for (int i = 0; i < nx; ++i) {
      const int idx = j * nx + i;
      if (j >= ng && j < ng + nint_y) {
        const double a = amp(x_of(i), y_of(j));
        vx[idx] = a * x_of(i);
        vy[idx] = a * y_of(j);
        vz[idx] = 0.0;
      } else {
        vx[idx] = vy[idx] = vz[idx] = 0.0;
      }
    }

  auto interp_comp = [&](const double *gf, int i, double x_src) {
    const int i_lo = i - box_half;
    const int nwin = 2 * box_half + 1;
    double row[32];
    for (int k = 0; k < nwin; ++k)
      row[k] = gf[j_ref * nx + (i_lo + k)];
    return interpolate_row(order, nwin, x_of(i_lo), dx, row, x_src);
  };

  double max_vx = 0, max_vy = 0, max_vz = 0;
  int n_checked = 0;
  for (int j = 0; j < ny; ++j) {
    if (j >= ng && j < ng + nint_y)
      continue;
    const double y = y_of(j);
    const double y_ref = y_of(j_ref);
    for (int i = box_half; i < nx - box_half; ++i) {
      const double x = x_of(i);
      const double rho = rho_of(x, y);
      const double x_src = invert_rho_to_x_signed(rho, y_ref, x);
      const auto delta = relative_rotation(x_src, y_ref, x, y);
      const double fx = interp_comp(vx, i, x_src);
      const double fy = interp_comp(vy, i, x_src);
      const double fz = interp_comp(vz, i, x_src);
      double tx, ty;
      rotate_u(delta, fx, fy, tx, ty);
      const double a = amp(x, y);
      max_vx = std::max(max_vx, std::fabs(tx - a * x));
      max_vy = std::max(max_vy, std::fabs(ty - a * y));
      max_vz = std::max(max_vz, std::fabs(fz));
      ++n_checked;
    }
  }
  expect_true(n_checked > 0, "radial vector fill checked at least one ghost point");
  expect_true(max_vx < 1e-3, "radial vector cartoon fill recovers Vx");
  expect_true(max_vy < 1e-3, "radial vector cartoon fill recovers Vy");
  expect_true(max_vz < 1e-12, "radial vector cartoon fill keeps Vz = 0");
  std::printf("      radial vector ghost maxabs err: vx=%.3g vy=%.3g vz=%.3g\n", max_vx, max_vy, max_vz);
}

// Half-plane cartoon: interior lives at x>=0 only. x<0 (hang-over for
// derivatives) is filled by unsigned invert + rotation from the +x side,
// with the interpolation stencil clamped to x>=0.
void test_half_x_fill_from_positive_canonical() {
  const int ng = 3, nint_y = 2, nint_x = 40, order = 4;
  const int nx = nint_x + ng; // ghosts only on the lower-x (axis) side
  const int ny = nint_y + 2 * ng;
  const int j_ref = ng;
  const int i_axis = ng; // first interior index, x=0
  const double dx = 0.05, dy = 0.05;
  auto x_of = [&](int i) { return (i - i_axis) * dx; };
  auto y_of = [&](int j) { return -0.5 * dy + (j - ng) * dy; };
  auto f = [&](double x, double y) { return std::exp(-(x * x + y * y)); };

  double u[ny * nx];
  for (int j = 0; j < ny; ++j)
    for (int i = 0; i < nx; ++i)
      u[j * nx + i] =
          (j >= ng && j < ng + nint_y && i >= i_axis) ? f(x_of(i), y_of(j))
                                                      : 0.0;

  auto interp_pos = [&](int /*i_target*/, double x_src) {
    const double i_real = i_axis + x_src / dx;
    const int i_center = int(std::lround(i_real));
    int offset = order / 2;
    if (i_center - offset < i_axis)
      offset = i_center - i_axis;
    if (i_center - offset + order >= nx)
      offset = i_center + order - nx + 1;
    if (i_center - offset < i_axis)
      offset = i_center - i_axis;
    const int i0 = i_center - offset;
    double yv[6];
    for (int n = 0; n <= order; ++n)
      yv[n] = u[j_ref * nx + i0 + n];
    return interpolate_lagrange(order, x_of(i0), dx, yv, x_src);
  };

  double max_abs_err = 0;
  int n_checked = 0;
  for (int j = 0; j < ny; ++j) {
    const bool y_ghost = j < ng || j >= ng + nint_y;
    for (int i = 0; i < nx; ++i) {
      const bool x_neg = i < i_axis;
      if (!y_ghost && !x_neg)
        continue;
      const double x = x_of(i), y = y_of(j);
      const double y_ref = y_of(j_ref);
      const double rho = rho_of(x, y);
      const double x_src = invert_rho_to_x(rho, y_ref); // always +x
      const auto delta = relative_rotation(x_src, y_ref, x, y);
      const double fs = interp_pos(i, x_src);
      double tx, ty;
      rotate_u(delta, x_src, y_ref, tx, ty); // unused; scalar copies fs
      (void)tx;
      (void)ty;
      const double err = std::fabs(fs - f(x, y));
      if (err > max_abs_err)
        max_abs_err = err;
      ++n_checked;
    }
  }
  expect_true(n_checked > 0, "half-x fill checked y-ghosts and x<0 hang-over");
  expect_true(max_abs_err < 1e-4,
              "half-x fill of exp(-rho^2) from x>=0 matches analytic");
  if (max_abs_err >= 1e-4)
    std::printf("      max_abs_err = %.3g over %d points\n", max_abs_err,
                n_checked);
}

// Cell-centered thin-y cartoon (HydroBaseX/TOV primitives): one interior
// y-cell whose center sits at y=0, plus ghosts. PointDesc.x/.y for a CCC
// group are the cell centers, so y_ref=0 and invert_rho_to_x(rho,0)=rho.
void test_cell_centered_half_x_fill() {
  const int ng = 2, nint_y = 1, nint_x = 40, order = 4;
  const int nx = nint_x + ng; // ghosts only on the lower-x (axis) side
  const int ny = nint_y + 2 * ng;
  const int j_ref = ng;   // the single interior cell
  const int i_axis = ng;  // first interior cell, center at x=+dx/2
  const double dx = 0.05, dy = 0.05;
  auto x_of = [&](int i) { return (i - i_axis + 0.5) * dx; };
  auto y_of = [&](int j) { return (j - j_ref) * dy; }; // interior y=0
  auto f = [&](double x, double y) { return std::exp(-(x * x + y * y)); };

  expect_near(y_of(j_ref), 0, 1e-15, "cell-centered interior y_ref is 0");
  expect_true(x_of(i_axis) > 0, "first interior cell sits at x=+dx/2");
  expect_true(x_of(i_axis - 1) < 0, "first x-ghost cell sits at x=-dx/2");
  expect_near(invert_rho_to_x(0.3, y_of(j_ref)), 0.3, 1e-14,
              "cell-centered invert_rho_to_x(rho, 0) == rho");

  double u[ny * nx];
  for (int j = 0; j < ny; ++j)
    for (int i = 0; i < nx; ++i)
      u[j * nx + i] =
          (j == j_ref && i >= i_axis) ? f(x_of(i), y_of(j)) : 0.0;

  auto interp_pos = [&](double x_src) {
    const double i_real = i_axis - 0.5 + x_src / dx;
    const int i_center = int(std::lround(i_real));
    int offset = order / 2;
    if (i_center - offset < i_axis)
      offset = i_center - i_axis;
    if (i_center - offset + order >= nx)
      offset = i_center + order - nx + 1;
    if (i_center - offset < i_axis)
      offset = i_center - i_axis;
    const int i0 = i_center - offset;
    double yv[6];
    for (int n = 0; n <= order; ++n)
      yv[n] = u[j_ref * nx + i0 + n];
    return interpolate_lagrange(order, x_of(i0), dx, yv, x_src);
  };

  double max_abs_err = 0;
  int n_checked = 0;
  for (int j = 0; j < ny; ++j) {
    const bool y_ghost = j != j_ref;
    for (int i = 0; i < nx; ++i) {
      const double x = x_of(i), y = y_of(j);
      const bool x_neg = x < 0;
      if (!y_ghost && !x_neg)
        continue;
      const double y_ref = y_of(j_ref);
      const double rho = rho_of(x, y);
      const double x_src = invert_rho_to_x(rho, y_ref);
      const double fs = interp_pos(x_src);
      const double err = std::fabs(fs - f(x, y));
      if (err > max_abs_err)
        max_abs_err = err;
      ++n_checked;
    }
  }
  expect_true(n_checked > 0,
              "cell-centered fill checked y-ghosts and x<0 hang-over");
  expect_true(max_abs_err < 1e-4,
              "cell-centered fill of exp(-rho^2) from x>=0 matches analytic");
  if (max_abs_err >= 1e-4)
    std::printf("      cell-centered max_abs_err = %.3g over %d points\n",
                max_abs_err, n_checked);
}

} // namespace

int main() {
  test_rho_of();
  test_relative_rotation_matches_old_cartoon2d_when_yref_zero();
  test_relative_rotation_offaxis_reference();
  test_invert_rho_to_x_near_axis_is_fatal();
  test_rotate_u_is_orthogonal_round_trip();
  test_rotate_ddsym_inplane_round_trip_and_trace();
  test_rotate_ddsym_inplane_matches_naive_congruence();
  test_rotate_ddsym_full_passthrough_and_vector_block();
  test_interpolate_lagrange_reproduces_polynomial();
  test_interpolate_lagrange_passes_through_samples();
  test_signed_invert_stays_on_the_same_side_of_the_axis();
  test_cartoon_fill_axisymmetric_gaussian_local_box();
  test_cartoon_fill_brill_lindquist_metric_local_box();
  test_cartoon_fill_radial_vector_local_box();
  test_half_x_fill_from_positive_canonical();
  test_cell_centered_half_x_fill();

  if (g_failures > 0) {
    std::printf("\n%d CHECK(S) FAILED\n", g_failures);
    return EXIT_FAILURE;
  }
  std::printf("\nAll checks passed.\n");
  return EXIT_SUCCESS;
}
