# Cartoon2DX

Cactus Code Thorn Cartoon2DX
Author(s)    : Steven R. Brandt <sbrandt@cct.lsu.edu>
Licence      : LGPL (same as SpacetimeX)

--------------------------------------------------------------------------

Cactus/CarpetX thorn: like `cactusnumerical/Cartoon2D`, but for the
CarpetX driver instead of Carpet/PUGH. Implements the "cartoon" trick for
running an axisymmetric evolution as a thin-in-y 3D evolution, filling the
y ghost zones by interpolating a reference profile in cylindrical radius
(rho) and rotating the tensor components into place.

This thorn lives in the SpacetimeX arrangement (alongside NewRadX, another
CarpetX physical boundary condition). Activate it and set
`CarpetX::boundary_y` / `boundary_upper_y` (and `boundary_x` for half-x)
to `"cartoon"`. A GetComponents list that checks out this branch plus the
matching CarpetX/AsterX changes is `Docs/thornlist/cartoon2dx.th`.

Paper: Alcubierre, Brandt, Brügmann, Holz, Seidel, Takahashi, Thornburg
(2001), "Symmetry without symmetry", gr-qc/9908012. Original cartoon
idea due to Steve Brandt; Cactus4 Cartoon2D by Sai Iyer.

## This thorn is a consumer-agnostic library, not a CottonmouthZ4c4m thorn

`Cartoon2DX` itself has no `INHERITS:` of any kind. It provides:

- `src/cartoon2dx.hxx`: rho/rotation geometry, tensor-rotation algebra for
  scalar/vector("u")/symmetric-tensor("ddsym") types, and 1D Lagrange
  interpolation -- all ported from the original Cartoon2D thorn's
  `Cartoon2DBC.c`/`interpolate.c`, generalized where needed for CarpetX's
  mixed centering (see notes.md). Fully tested standalone (30/30 checks).
- `src/cartoon2dx_carpetx.hxx`: CarpetX-specific glue -- gathers a
  stencil from an actual `GF3D5<CCTK_REAL>` grid function and loops over
  y-ghost points with `grid.loop_all_device` (filtered to the y faces
  only), exposing `fill_ghosts_scalar`/`fill_ghosts_vector`/
  `fill_ghosts_ddsym`. Compiled and exercised against a real CarpetX/AMReX
  build.
- `src/check_parameters.cxx`: `Cartoon2DX_CheckParameters` at
  `CCTK_PARAMCHECK`.
- `src/apply.cxx`: `ApplyCartoonBoundary(cctkGH, gi, tl)`, called from
  CarpetX `apply_boundary_conditions` when a face is `boundary_t::cartoon`,
  for the timelevel that owns the MultiFab.
  Classifies each GF group from `tensortypealias` (Scalar / U or D /
  DD_sym -- the same Einstein Toolkit tag old Cartoon2D reads) or, if
  that is absent, from the CarpetX `parities` tag (1 / 3 / 6
  variables). Set `CarpetX::boundary_y` / `boundary_upper_y` (and
  `boundary_x` for half-x) to `"cartoon"`.
- `param.ccl`'s `order`/`verbose`/`fill_negative_x`.
- `test/standalone_math_test.cxx`: correctness tests for
  `cartoon2dx.hxx`, runnable directly with `g++` (no Cactus/AMReX build
  required -- see that file's header comment).

To use Cartoon2DX, activate the thorn and tag GF groups with
`tensortypealias` (or rely on existing `parities=` on generated
CarpetX thorns). Do not add a glue thorn.

### Why there is no `INHERITS:` (learned the hard way, not designed up front)

The very first version of this thorn had `INHERITS: CottonmouthZ4c4m`
directly in its own `interface.ccl`, with the 4 hook functions living
right here. That seemed reasonable until an actual smoke-test run of
`Cartoon2DXWaveToyX` failed activation outright:

```
Error: Implementation 'ADMBASEX' not activated.
       This implementation is required by activated thorn(s):
           Cartoon2DX (implementing Cartoon2DX)
Error: Implementation 'COTTONMOUTHZ4C4M' not activated.
Error: Implementation 'TMUNUBASEX' not activated.
```

`INHERITS:` in Cactus isn't "grant access if the other thorn happens to
be active" -- it's a hard requirement that the *implementation* be
active, transitively. So `Cartoon2DX` could never be activated for an
unrelated consumer without also activating ADMBaseX/CottonmouthZ4c4m/
TmunuBaseX. The generic after-sync fill needs no `INHERITS:`: it looks
up groups by index and writes through `CCTK_VarDataPtrI`.

The leftover `Cartoon2DXZ4c4m` / `Cartoon2DXWaveToyX` thorns are empty
stubs so old par files that still list them continue to activate.

## Status

The generic after-sync BC plus half-x (`xmin=0`, lower-x ghosts filled
by rotation) is the current integration. The wave-equation isolation
test tracks WaveToyX's exact Gaussian at 2nd-order FD accuracy. The
Brill-Lindquist head-on test evolves real punctures through t=10
(gxx.max stays O(100–600); the holes fall together along z). The old
origin-gxx=1 finding was NewRadX replacing the Einstein RHS on both
thin-y interior layers; cartoon pars set `NewRadX::apply_y=no` and
`NewRadX::apply_lower_x=no`. At t=1, cartoon vs full3d gxx agrees to a
few percent in the bulk (exact match at t=0). See notes.md.

Not yet implemented (see notes.md "Status" for reasoning):

- ENO interpolation and excision support (present in the original
  Cartoon2D thorn, not needed for the vacuum Brill-Lindquist
  head-on-collision correctness test that motivates this port).
- Automatic grid resizing (not needed: CarpetX grid extents are set
  directly via `CarpetX::ncells_x/y/z` etc.).

The near-axis special case (`rho < |y_ref|`) in `invert_rho_to_x()`
remains a loud, fatal error in general, but is now *proved* unreachable
via the actual ghost-fill call path (see cartoon2dx.hxx's comment on
`invert_rho_to_x`) -- not just deferred.
