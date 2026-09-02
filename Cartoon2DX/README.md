# Cartoon2DX

Cactus/CarpetX thorn: like `cactusnumerical/Cartoon2D`, but for CarpetX
instead of Carpet/PUGH. An axisymmetric evolution is a thin-in-y 3D
slab. Ghosts on the y-faces (and optional x<0 hang-over) are filled by
interpolating in cylindrical radius $\rho=\sqrt{x^2+y^2}$ and rotating
tensor components into place.

This thorn lives in the SpacetimeX arrangement, next to NewRadX (another
CarpetX physical BC). Licence: LGPL.

Paper: Alcubierre, Brandt, Brügmann, Holz, Seidel, Takahashi, Thornburg
(2001), [gr-qc/9908012](https://arxiv.org/abs/gr-qc/9908012). Original
cartoon idea: Steve Brandt. Cactus4 Cartoon2D: Sai Iyer.

Full ThornGuide write-up: `doc/documentation.tex`.

## Using it

Activate `Cartoon2DX` and set the cartoon faces:

```
CarpetX::boundary_y       = "cartoon"
CarpetX::boundary_upper_y = "cartoon"
CarpetX::boundary_x       = "cartoon"   # half-x; xmin is the z-axis
CarpetX::boundary_upper_x = "neumann"   # or linear extrapolation
CarpetX::ncells_y = 1                   # one cell; y=0 at the centre
Cartoon2DX::fill_negative_x = yes
Cartoon2DX::order = 4                   # Lagrange, order+1 points
```

`ghost_size` must be large enough for the interpolation stencil
(`order=4` needs 5 points, so `ghost_size >= 2`).

Tag each grid-function group with Einstein Toolkit `tensortypealias`
(`Scalar`, `U`/`D`, `DD_sym`). If that tag is missing, Cartoon2DX falls
back to the group's `parities` tag (1 / 3 / 6 variables). Packed groups
such as AsterX `cons_vector` are walked component-wise from parities.
Untagged 1-variable GFs are treated as scalars. Mixed CarpetX centering
is taken from each group's `CENTERING` table; $\rho$ uses
`PointDesc.x` / `.y`.

There is no `INHERITS:` and no per-consumer glue thorn.

GetComponents list (this branch plus the matching CarpetX/AsterX
changes): `Docs/thornlist/cartoon2dx.th` in SpacetimeX.

## Tests

| Test | What it is |
|---|---|
| `test/tov.par` | Cactus testsuite: low-res Cowling TOV, 2 RK4 steps. Meant to run on a laptop in well under 30 s. |
| `test/standalone_math_test.cxx` | Interpolation and tensor-rotation math, no Cactus/AMReX. |

Standalone math test:

```
cd arrangements/SpacetimeX/Cartoon2DX/test
g++ -std=c++17 -Wall -Wextra -I../src -o standalone_math_test standalone_math_test.cxx
./standalone_math_test
```

Cactus testsuite (from a configuration that includes Cartoon2DX, AsterX,
and TOVSolverX):

```
make sim-testsuite TESTPROCS=1 tests=tov
```

Longer comparison pars (`tov_cartoon2dx*.par`, `head_on_*.par`) live in
`test/` as well; they are not part of the testsuite.

## Parameters

| Parameter | Default | Meaning |
|---|---|---|
| `order` | 4 | Lagrange interpolation order (`order+1` points), 1–5 |
| `fill_negative_x` | yes | Fill x<0 from x>=0 (half-x domain, `xmin=0`) |
| `verbose` | no | Log each group the first time it is filled |
