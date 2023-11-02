#ifndef FIELD_HXX
#define FIELD_HXX

#include <loop.hxx>

#include <cassert>
#include <cstddef>
#include <vector>

namespace Z4c {
using namespace Loop;
using namespace std;

template <typename T> class field {
  const GridDescBase &restrict grid;
  vect<int, dim> imin, imax;
  vect<ptrdiff_t, dim> istr;
  vector<T> elts;

  constexpr ptrdiff_t idx(const vect<T, dim> &I) const {
#ifdef CCTK_DEBUG
    for (int d = 0; d < dim; ++d)
      assert(I[d] >= imin[d] && I[d] < imax[d]);
#endif
    ptrdiff_t n = 0;
    if (dim > 0) {
      n += I[0];
      for (int d = 1; d < dim; ++d)
        n += istr[d] * I[d];
    }
#ifdef CCTK_DEBUG
    assert(n >= 0 && n < ptrdiff_t(elts.size()));
#endif
    return n;
  }

public:
  field(const GridDescBase &grid) : grid(grid) {
    grid.box_all<1, 1, 1>(grid.nghostzones, imin, imax);
    ptrdiff_t str = 1;
    for (int d = 0; d < dim; ++d) {
      istr[d] = str;
      str *= imax[d] - imin[d];
    }
    elts.resize(str);
  }

  const T &operator()(const vect<T, dim> &I) const { return elts[idx(I)]; }
  T &operator()(const vect<T, dim> &I) { return elts[idx(I)]; }

  template <int dir>
  field deriv(const int order, const vect<T, dim> &dx) const {
    static_assert(dir >= 0 && dir < dim, "");
    constexpr auto DI = vect<int, dim>::unit(dir);
    const auto &f = *this;
    field r(grid);
    switch (order) {
    case 2:
      assert(grid.nghostzones[dir] >= 1);
      loop_int(grid.nghostzones, [&](const PointDesc &p) {
        r(p.I) = -1 / T(2) * (f(p.I - DI) - f(p.I + DI)) / dx[dir];
      });
      break;
    case 4:
      assert(grid.nghostzones[dir] >= 2);
      loop_int(grid.nghostzones, [&](const PointDesc &p) {
        r(p.I) = (1 / T(12) * (f(p.I - 2 * DI) - f(p.I + 2 * DI)) +
                  2 / T(3) * (f(p.I - DI) - f(p.I + DI))) /
                 dx[dir];
      });
      break;
    default:
      assert(0);
    }
    return r;
  }
};

} // namespace Z4c

#endif // #ifndef FIELD_HXX
