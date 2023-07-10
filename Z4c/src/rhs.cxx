#include <cctk.h>

#ifdef __CUDACC__
// Disable CCTK_DEBUG since the debug information takes too much
// parameter space to launch the kernels
#ifdef CCTK_DEBUG
#undef CCTK_DEBUG
#endif
#endif

#include "derivs.hxx"
#include "physics.hxx"
#include "z4c_vars.hxx"

#include <loop_device.hxx>
#include <mat.hxx>
#include <simd.hxx>
#include <vec.hxx>

#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>

#ifdef __CUDACC__
#include <nvToolsExt.h>
#endif

#include <cmath>

namespace Z4c {
using namespace Arith;
using namespace Loop;
using namespace std;

typedef simd<CCTK_REAL> vreal;
typedef simdl<CCTK_REAL> vbool;
constexpr size_t vsize = tuple_size_v<vreal>;

// helpers to access tile functions as a 4D array
template <typename T> constexpr auto next(const T&) { return T::end; };

template <int begin>
struct scalar_t {
  enum {idx = begin};
  enum {end = begin + 1};

  vreal operator()(const vbool mask, const GF3D5index& index,
    const GF3D5<CCTK_REAL>& tmp0, const int np) const CCTK_HOST CCTK_DEVICE {
    return tmp0(mask, index, idx*np);
  }
  vreal operator()(const vbool mask, const GF3D2index& index,
    const GF3D2<const CCTK_REAL>& tmp0, const int np) const CCTK_HOST CCTK_DEVICE {
    return tmp0(mask, index, idx*np);
  }
  void store(const vbool mask, const GF3D2index& index,
    const GF3D2<CCTK_REAL>& tmp0, const int np,
    const vreal& value) const CCTK_HOST CCTK_DEVICE {
    tmp0.store(mask, index, idx*np, value);
  }
};

template <int begin>
struct vec_t {
  enum {idx_x = begin, idx_y, idx_z};
  enum {end = begin + 3};

  vec<vreal, 3> operator()(const vbool mask, const GF3D5index& index,
    const GF3D5<CCTK_REAL>& tmp0, const int np) const CCTK_HOST CCTK_DEVICE {
    return vec<vreal, 3>({tmp0(mask, index, idx_x*np), tmp0(mask, index, idx_y*np), tmp0(mask, index, idx_z*np)});
  }
  vec<vreal, 3> operator()(const vbool mask, const GF3D2index& index,
    const GF3D2<const CCTK_REAL>& tmp0, const int np) const CCTK_HOST CCTK_DEVICE {
    return vec<vreal, 3>({tmp0(mask, index, idx_x*np), tmp0(mask, index, idx_y*np), tmp0(mask, index, idx_z*np)});
  }
  void store(const vbool mask, const GF3D2index& index,
    const GF3D2<CCTK_REAL>& tmp0, const int np,
    const vec<vreal, 3>& value) const CCTK_HOST CCTK_DEVICE {
    tmp0.store(mask, index, idx_x*np, value.elts[idx_x-idx_x]);
    tmp0.store(mask, index, idx_y*np, value.elts[idx_y-idx_x]);
    tmp0.store(mask, index, idx_z*np, value.elts[idx_z-idx_x]);
  }
};

template <int begin>
struct mat_t {
  enum {idx_xx = begin, idx_xy, idx_xz, idx_yy, idx_yz, idx_zz};
  enum {end = begin + 6};

  smat<vreal, 3> operator()(const vbool mask, const GF3D5index& index,
    const GF3D5<CCTK_REAL>& tmp0, const int np) const CCTK_HOST CCTK_DEVICE {
    return smat<vreal, 3>({tmp0(mask, index, idx_xx*np), tmp0(mask, index, idx_xy*np),
      tmp0(mask, index, idx_xz*np), tmp0(mask, index, idx_yy*np), tmp0(mask, index, idx_yz*np),
      tmp0(mask, index, idx_zz*np)});
  }
  smat<vreal, 3> operator()(const vbool mask, const GF3D2index& index,
    const GF3D2<const CCTK_REAL>& tmp0, const int np) const CCTK_HOST CCTK_DEVICE {
    return smat<vreal, 3>({tmp0(mask, index, idx_xx*np), tmp0(mask, index, idx_xy*np),
      tmp0(mask, index, idx_xz*np), tmp0(mask, index, idx_yy*np), tmp0(mask, index, idx_yz*np),
      tmp0(mask, index, idx_zz*np)});
  }
  void store(const vbool mask, const GF3D2index& index,
    const GF3D2<CCTK_REAL>& tmp0, const int np,
    const smat<vreal, 3>& value) const CCTK_HOST CCTK_DEVICE {
    tmp0.store(mask, index, idx_xx*np, value.elts[idx_xx-idx_xx]);
    tmp0.store(mask, index, idx_xy*np, value.elts[idx_xy-idx_xx]);
    tmp0.store(mask, index, idx_xz*np, value.elts[idx_xz-idx_xx]);
    tmp0.store(mask, index, idx_yy*np, value.elts[idx_yy-idx_xx]);
    tmp0.store(mask, index, idx_yz*np, value.elts[idx_yz-idx_xx]);
    tmp0.store(mask, index, idx_zz*np, value.elts[idx_zz-idx_xx]);
  }
};

template <int begin>
struct vec_vec_t {
  enum {idx_xx = begin, idx_xy, idx_xz, idx_yx, idx_yy, idx_yz,
        idx_zx, idx_zy, idx_zz};
  enum {end = begin + 9};

  vec<vec<vreal, 3>, 3> operator()(const vbool mask, const GF3D5index& index,
    const GF3D5<CCTK_REAL>& tmp0, const int np) const CCTK_HOST CCTK_DEVICE {
    return vec<vec<vreal, 3>, 3>(
      {vec<vreal, 3>({tmp0(mask, index, idx_xx*np), tmp0(mask, index, idx_xy*np), tmp0(mask, index, idx_xz*np)}),
       vec<vreal, 3>({tmp0(mask, index, idx_yx*np), tmp0(mask, index, idx_yy*np), tmp0(mask, index, idx_yz*np)}),
       vec<vreal, 3>({tmp0(mask, index, idx_zx*np), tmp0(mask, index, idx_zy*np), tmp0(mask, index, idx_zz*np)})});
  }
  vec<vec<vreal, 3>, 3> operator()(const vbool mask, const GF3D2index& index,
    const GF3D2<const CCTK_REAL>& tmp0, const int np) const CCTK_HOST CCTK_DEVICE {
    return vec<vec<vreal, 3>, 3>(
      {vec<vreal, 3>({tmp0(mask, index, idx_xx*np), tmp0(mask, index, idx_xy*np), tmp0(mask, index, idx_xz*np)}),
       vec<vreal, 3>({tmp0(mask, index, idx_yx*np), tmp0(mask, index, idx_yy*np), tmp0(mask, index, idx_yz*np)}),
       vec<vreal, 3>({tmp0(mask, index, idx_zx*np), tmp0(mask, index, idx_zy*np), tmp0(mask, index, idx_zz*np)})});
  }
  void store(const vbool mask, const GF3D2index& index,
    const GF3D2<CCTK_REAL>& tmp0, const int np,
    const vec<vec<vreal, 3>, 3>& value) const CCTK_HOST CCTK_DEVICE {
    tmp0.store(mask, index, idx_xx*np, value.elts[idx_xx-idx_xx]);
    tmp0.store(mask, index, idx_xy*np, value.elts[idx_xy-idx_xx]);
    tmp0.store(mask, index, idx_xz*np, value.elts[idx_xz-idx_xx]);
    tmp0.store(mask, index, idx_yx*np, value.elts[idx_yx-idx_xx]);
    tmp0.store(mask, index, idx_yy*np, value.elts[idx_yy-idx_xx]);
    tmp0.store(mask, index, idx_yz*np, value.elts[idx_yz-idx_xx]);
    tmp0.store(mask, index, idx_xz*np, value.elts[idx_zx-idx_xx]);
    tmp0.store(mask, index, idx_xz*np, value.elts[idx_zy-idx_xx]);
    tmp0.store(mask, index, idx_xz*np, value.elts[idx_zz-idx_xx]);
  }
};

template <int begin>
struct vec_mat_t {
  enum {idx_xxx = begin, idx_xxy, idx_xxz, idx_xyy, idx_xyz, idx_xzz,
        idx_yxx, idx_yxy, idx_yxz, idx_yyy, idx_yyz, idx_yzz,
        idx_zxx, idx_zxy, idx_zxz, idx_zyy, idx_zyz, idx_zzz};
  enum {end = begin + 18};

  vec<smat<vreal, 3>, 3> operator()(const vbool mask, const GF3D5index& index,
    const GF3D5<CCTK_REAL>& tmp0, const int np) const CCTK_HOST CCTK_DEVICE {
    return vec<smat<vreal, 3>, 3>(
      {smat<vreal, 3>({tmp0(mask, index, idx_xxx*np), tmp0(mask, index, idx_xxy*np), tmp0(mask, index, idx_xxz*np), tmp0(mask, index, idx_xyy*np), tmp0(mask, index, idx_xyz*np), tmp0(mask, index, idx_xzz*np)}),
       smat<vreal, 3>({tmp0(mask, index, idx_yxx*np), tmp0(mask, index, idx_yxy*np), tmp0(mask, index, idx_yxz*np), tmp0(mask, index, idx_yyy*np), tmp0(mask, index, idx_yyz*np), tmp0(mask, index, idx_yzz*np)}),
       smat<vreal, 3>({tmp0(mask, index, idx_zxx*np), tmp0(mask, index, idx_zxy*np), tmp0(mask, index, idx_zxz*np), tmp0(mask, index, idx_zyy*np), tmp0(mask, index, idx_zyz*np), tmp0(mask, index, idx_zzz*np)})});
  }
  vec<smat<vreal, 3>, 3> operator()(const vbool mask, const GF3D2index& index,
    const GF3D2<const CCTK_REAL>& tmp0, const int np) const CCTK_HOST CCTK_DEVICE {
    return vec<smat<vreal, 3>, 3>(
      {smat<vreal, 3>({tmp0(mask, index, idx_xxx*np), tmp0(mask, index, idx_xxy*np), tmp0(mask, index, idx_xxz*np), tmp0(mask, index, idx_xyy*np), tmp0(mask, index, idx_xyz*np), tmp0(mask, index, idx_xzz*np)}),
       smat<vreal, 3>({tmp0(mask, index, idx_yxx*np), tmp0(mask, index, idx_yxy*np), tmp0(mask, index, idx_yxz*np), tmp0(mask, index, idx_yyy*np), tmp0(mask, index, idx_yyz*np), tmp0(mask, index, idx_yzz*np)}),
       smat<vreal, 3>({tmp0(mask, index, idx_zxx*np), tmp0(mask, index, idx_zxy*np), tmp0(mask, index, idx_zxz*np), tmp0(mask, index, idx_zyy*np), tmp0(mask, index, idx_zyz*np), tmp0(mask, index, idx_zzz*np)})});
  }
  void store(const vbool mask, const GF3D2index& index,
    const GF3D2<CCTK_REAL>& tmp0, const int np,
    const vec<smat<vreal, 3>, 3>& value) const CCTK_HOST CCTK_DEVICE {
    tmp0.store(mask, index, idx_xxx*np, value.elts[idx_xxx-idx_xxx]);
    tmp0.store(mask, index, idx_xxy*np, value.elts[idx_xxy-idx_xxx]);
    tmp0.store(mask, index, idx_xxz*np, value.elts[idx_xxz-idx_xxx]);
    tmp0.store(mask, index, idx_xyy*np, value.elts[idx_xyy-idx_xxx]);
    tmp0.store(mask, index, idx_xyz*np, value.elts[idx_xyz-idx_xxx]);
    tmp0.store(mask, index, idx_xzz*np, value.elts[idx_xzz-idx_xxx]);
    tmp0.store(mask, index, idx_yxx*np, value.elts[idx_yxx-idx_xxx]);
    tmp0.store(mask, index, idx_yxy*np, value.elts[idx_yxy-idx_xxx]);
    tmp0.store(mask, index, idx_yxz*np, value.elts[idx_yxz-idx_xxx]);
    tmp0.store(mask, index, idx_yyy*np, value.elts[idx_yyy-idx_xxx]);
    tmp0.store(mask, index, idx_yyz*np, value.elts[idx_yyz-idx_xxx]);
    tmp0.store(mask, index, idx_yzz*np, value.elts[idx_yzz-idx_xxx]);
    tmp0.store(mask, index, idx_zxx*np, value.elts[idx_zxx-idx_xxx]);
    tmp0.store(mask, index, idx_zxy*np, value.elts[idx_zxy-idx_xxx]);
    tmp0.store(mask, index, idx_zxz*np, value.elts[idx_zxz-idx_xxx]);
    tmp0.store(mask, index, idx_zyy*np, value.elts[idx_zyy-idx_xxx]);
    tmp0.store(mask, index, idx_zyz*np, value.elts[idx_zyz-idx_xxx]);
    tmp0.store(mask, index, idx_zzz*np, value.elts[idx_zzz-idx_xxx]);
  }
};

template <int begin>
struct mat_vec_t {
  enum {idx_xxx = begin, idx_xxy, idx_xxz, idx_xyx, idx_xyy, idx_xyz,
        idx_xzx, idx_xzy, idx_xzz, idx_yyx, idx_yyy, idx_yyz,
        idx_yzx, idx_yzy, idx_yzz, idx_zzx, idx_zzy, idx_zzz};
  enum {end = begin + 18};

  smat<vec<vreal, 3>, 3> operator()(const vbool mask, const GF3D5index& index,
    const GF3D5<CCTK_REAL>& tmp0, const int np) const CCTK_HOST CCTK_DEVICE {
    return smat<vec<vreal, 3>, 3>(
      {vec<vreal, 3>({tmp0(mask, index, idx_xxx*np), tmp0(mask, index, idx_xxy*np), tmp0(mask, index, idx_xxz*np)}),
       vec<vreal, 3>({tmp0(mask, index, idx_xyx*np), tmp0(mask, index, idx_xyy*np), tmp0(mask, index, idx_xyz*np)}),
       vec<vreal, 3>({tmp0(mask, index, idx_xzx*np), tmp0(mask, index, idx_xzy*np), tmp0(mask, index, idx_xzz*np)}),
       vec<vreal, 3>({tmp0(mask, index, idx_yyx*np), tmp0(mask, index, idx_yyy*np), tmp0(mask, index, idx_yyz*np)}),
       vec<vreal, 3>({tmp0(mask, index, idx_yzx*np), tmp0(mask, index, idx_yzy*np), tmp0(mask, index, idx_yzz*np)}),
       vec<vreal, 3>({tmp0(mask, index, idx_zzx*np), tmp0(mask, index, idx_zzy*np), tmp0(mask, index, idx_zzz*np)})});
  }
  smat<vec<vreal, 3>, 3> operator()(const vbool mask, const GF3D2index& index,
    const GF3D2<const CCTK_REAL>& tmp0, const int np) const CCTK_HOST CCTK_DEVICE {
    return smat<vec<vreal, 3>, 3>(
      {vec<vreal, 3>({tmp0(mask, index, idx_xxx*np), tmp0(mask, index, idx_xxy*np), tmp0(mask, index, idx_xxz*np)}),
       vec<vreal, 3>({tmp0(mask, index, idx_xyx*np), tmp0(mask, index, idx_xyy*np), tmp0(mask, index, idx_xyz*np)}),
       vec<vreal, 3>({tmp0(mask, index, idx_xzx*np), tmp0(mask, index, idx_xzy*np), tmp0(mask, index, idx_xzz*np)}),
       vec<vreal, 3>({tmp0(mask, index, idx_yyx*np), tmp0(mask, index, idx_yyy*np), tmp0(mask, index, idx_yyz*np)}),
       vec<vreal, 3>({tmp0(mask, index, idx_yzx*np), tmp0(mask, index, idx_yzy*np), tmp0(mask, index, idx_yzz*np)}),
       vec<vreal, 3>({tmp0(mask, index, idx_zzx*np), tmp0(mask, index, idx_zzy*np), tmp0(mask, index, idx_zzz*np)})});
  }
  void store(const vbool mask, const GF3D2index& index,
    const GF3D2<CCTK_REAL>& tmp0, const int np,
    const smat<vec<vreal, 3>, 3>& value) const {
    tmp0.store(mask, index, idx_xxx*np, value.elts[idx_xxx-idx_xxx]);
    tmp0.store(mask, index, idx_xxy*np, value.elts[idx_xxy-idx_xxx]);
    tmp0.store(mask, index, idx_xxz*np, value.elts[idx_xxz-idx_xxx]);
    tmp0.store(mask, index, idx_xyx*np, value.elts[idx_xyx-idx_xxx]);
    tmp0.store(mask, index, idx_xyy*np, value.elts[idx_xyy-idx_xxx]);
    tmp0.store(mask, index, idx_xyz*np, value.elts[idx_xyz-idx_xxx]);
    tmp0.store(mask, index, idx_xzx*np, value.elts[idx_xzx-idx_xxx]);
    tmp0.store(mask, index, idx_xzy*np, value.elts[idx_xzy-idx_xxx]);
    tmp0.store(mask, index, idx_xzz*np, value.elts[idx_xzz-idx_xxx]);
    tmp0.store(mask, index, idx_yyx*np, value.elts[idx_yyx-idx_xxx]);
    tmp0.store(mask, index, idx_yyy*np, value.elts[idx_yyy-idx_xxx]);
    tmp0.store(mask, index, idx_yyz*np, value.elts[idx_yyz-idx_xxx]);
    tmp0.store(mask, index, idx_yzx*np, value.elts[idx_yzx-idx_xxx]);
    tmp0.store(mask, index, idx_yzy*np, value.elts[idx_yzy-idx_xxx]);
    tmp0.store(mask, index, idx_yzz*np, value.elts[idx_yzz-idx_xxx]);
    tmp0.store(mask, index, idx_zzx*np, value.elts[idx_zzx-idx_xxx]);
    tmp0.store(mask, index, idx_zzy*np, value.elts[idx_zzy-idx_xxx]);
    tmp0.store(mask, index, idx_zzz*np, value.elts[idx_zzz-idx_xxx]);
  }
};

template <int begin>
struct mat_mat_t {
  enum {idx_xxxx = begin, idx_xxxy, idx_xxxz, idx_xxyy, idx_xxyz, idx_xxzz,
        idx_xyxx, idx_xyxy, idx_xyxz, idx_xyyy, idx_xyyz, idx_xyzz,
        idx_xzxx, idx_xzxy, idx_xzxz, idx_xzyy, idx_xzyz, idx_xzzz,
        idx_yyxx, idx_yyxy, idx_yyxz, idx_yyyy, idx_yyyz, idx_yyzz,
        idx_yzxx, idx_yzxy, idx_yzxz, idx_yzyy, idx_yzyz, idx_yzzz,
        idx_zzxx, idx_zzxy, idx_zzxz, idx_zzyy, idx_zzyz, idx_zzzz};
  enum {end = begin + 36};

  smat<smat<vreal, 3>, 3> operator()(const vbool mask, const GF3D5index& index,
    const GF3D5<CCTK_REAL>& tmp0, const int np) const CCTK_HOST CCTK_DEVICE {
    return 
      smat<smat<vreal, 3>, 3>(
      {smat<vreal, 3>({tmp0(mask, index, idx_xxxx*np), tmp0(mask, index, idx_xxxy*np), tmp0(mask, index, idx_xxxz*np), tmp0(mask, index, idx_xxyy*np), tmp0(mask, index, idx_xxyz*np), tmp0(mask, index, idx_xxzz*np)}),
       smat<vreal, 3>({tmp0(mask, index, idx_xyxx*np), tmp0(mask, index, idx_xyxy*np), tmp0(mask, index, idx_xyxz*np), tmp0(mask, index, idx_xyyy*np), tmp0(mask, index, idx_xyyz*np), tmp0(mask, index, idx_xyzz*np)}),
       smat<vreal, 3>({tmp0(mask, index, idx_xzxx*np), tmp0(mask, index, idx_xzxy*np), tmp0(mask, index, idx_xzxz*np), tmp0(mask, index, idx_xzyy*np), tmp0(mask, index, idx_xzyz*np), tmp0(mask, index, idx_xzzz*np)}),
       smat<vreal, 3>({tmp0(mask, index, idx_yyxx*np), tmp0(mask, index, idx_yyxy*np), tmp0(mask, index, idx_yyxz*np), tmp0(mask, index, idx_yyyy*np), tmp0(mask, index, idx_yyyz*np), tmp0(mask, index, idx_yyzz*np)}),
       smat<vreal, 3>({tmp0(mask, index, idx_yzxx*np), tmp0(mask, index, idx_yzxy*np), tmp0(mask, index, idx_yzxz*np), tmp0(mask, index, idx_yzyy*np), tmp0(mask, index, idx_yzyz*np), tmp0(mask, index, idx_yzzz*np)}),
       smat<vreal, 3>({tmp0(mask, index, idx_zzxx*np), tmp0(mask, index, idx_zzxy*np), tmp0(mask, index, idx_zzxz*np), tmp0(mask, index, idx_zzyy*np), tmp0(mask, index, idx_zzyz*np), tmp0(mask, index, idx_zzzz*np)})});
  }
  smat<smat<vreal, 3>, 3> operator()(const vbool mask, const GF3D2index& index,
    const GF3D2<const CCTK_REAL>& tmp0, const int np) const CCTK_HOST CCTK_DEVICE {
    return 
      smat<smat<vreal, 3>, 3>(
      {smat<vreal, 3>({tmp0(mask, index, idx_xxxx*np), tmp0(mask, index, idx_xxxy*np), tmp0(mask, index, idx_xxxz*np), tmp0(mask, index, idx_xxyy*np), tmp0(mask, index, idx_xxyz*np), tmp0(mask, index, idx_xxzz*np)}),
       smat<vreal, 3>({tmp0(mask, index, idx_xyxx*np), tmp0(mask, index, idx_xyxy*np), tmp0(mask, index, idx_xyxz*np), tmp0(mask, index, idx_xyyy*np), tmp0(mask, index, idx_xyyz*np), tmp0(mask, index, idx_xyzz*np)}),
       smat<vreal, 3>({tmp0(mask, index, idx_xzxx*np), tmp0(mask, index, idx_xzxy*np), tmp0(mask, index, idx_xzxz*np), tmp0(mask, index, idx_xzyy*np), tmp0(mask, index, idx_xzyz*np), tmp0(mask, index, idx_xzzz*np)}),
       smat<vreal, 3>({tmp0(mask, index, idx_yyxx*np), tmp0(mask, index, idx_yyxy*np), tmp0(mask, index, idx_yyxz*np), tmp0(mask, index, idx_yyyy*np), tmp0(mask, index, idx_yyyz*np), tmp0(mask, index, idx_yyzz*np)}),
       smat<vreal, 3>({tmp0(mask, index, idx_yzxx*np), tmp0(mask, index, idx_yzxy*np), tmp0(mask, index, idx_yzxz*np), tmp0(mask, index, idx_yzyy*np), tmp0(mask, index, idx_yzyz*np), tmp0(mask, index, idx_yzzz*np)}),
       smat<vreal, 3>({tmp0(mask, index, idx_zzxx*np), tmp0(mask, index, idx_zzxy*np), tmp0(mask, index, idx_zzxz*np), tmp0(mask, index, idx_zzyy*np), tmp0(mask, index, idx_zzyz*np), tmp0(mask, index, idx_zzzz*np)})});
  }
  void store(const vbool mask, const GF3D2index& index,
    const GF3D2<CCTK_REAL>& tmp0, const int np,
    const smat<smat<vreal, 3>, 3>& value) const {
    tmp0.store(mask, index, idx_xxxx*np, value.elts[idx_xxxx-idx_xxxx]);
    tmp0.store(mask, index, idx_xxxy*np, value.elts[idx_xxxy-idx_xxxx]);
    tmp0.store(mask, index, idx_xxxz*np, value.elts[idx_xxxz-idx_xxxx]);
    tmp0.store(mask, index, idx_xxyy*np, value.elts[idx_xxyy-idx_xxxx]);
    tmp0.store(mask, index, idx_xxyz*np, value.elts[idx_xxyz-idx_xxxx]);
    tmp0.store(mask, index, idx_xxzz*np, value.elts[idx_xxzz-idx_xxxx]);

    tmp0.store(mask, index, idx_xyxx*np, value.elts[idx_xyxy-idx_xxxx]);
    tmp0.store(mask, index, idx_xyxy*np, value.elts[idx_xyxy-idx_xxxx]);
    tmp0.store(mask, index, idx_xyxz*np, value.elts[idx_xyxz-idx_xxxx]);
    tmp0.store(mask, index, idx_xyyy*np, value.elts[idx_xyyy-idx_xxxx]);
    tmp0.store(mask, index, idx_xyyz*np, value.elts[idx_xyyz-idx_xxxx]);
    tmp0.store(mask, index, idx_xyzz*np, value.elts[idx_xyzz-idx_xxxx]);

    tmp0.store(mask, index, idx_xzxx*np, value.elts[idx_xzxz-idx_xxxx]);
    tmp0.store(mask, index, idx_xzxy*np, value.elts[idx_xzxy-idx_xxxx]);
    tmp0.store(mask, index, idx_xzxz*np, value.elts[idx_xzxz-idx_xxxx]);
    tmp0.store(mask, index, idx_xzyy*np, value.elts[idx_xzyy-idx_xxxx]);
    tmp0.store(mask, index, idx_xzyz*np, value.elts[idx_xzyz-idx_xxxx]);
    tmp0.store(mask, index, idx_xzzz*np, value.elts[idx_xzzz-idx_xxxx]);

    tmp0.store(mask, index, idx_yyxx*np, value.elts[idx_yyyy-idx_xxxx]);
    tmp0.store(mask, index, idx_yyxy*np, value.elts[idx_yyxy-idx_xxxx]);
    tmp0.store(mask, index, idx_yyxz*np, value.elts[idx_yyxz-idx_xxxx]);
    tmp0.store(mask, index, idx_yyyy*np, value.elts[idx_yyyy-idx_xxxx]);
    tmp0.store(mask, index, idx_yyyz*np, value.elts[idx_yyyz-idx_xxxx]);
    tmp0.store(mask, index, idx_yyzz*np, value.elts[idx_yyzz-idx_xxxx]);

    tmp0.store(mask, index, idx_yzxx*np, value.elts[idx_yzyz-idx_xxxx]);
    tmp0.store(mask, index, idx_yzxy*np, value.elts[idx_yzxy-idx_xxxx]);
    tmp0.store(mask, index, idx_yzxz*np, value.elts[idx_yzxz-idx_xxxx]);
    tmp0.store(mask, index, idx_yzyy*np, value.elts[idx_yzyy-idx_xxxx]);
    tmp0.store(mask, index, idx_yzyz*np, value.elts[idx_yzyz-idx_xxxx]);
    tmp0.store(mask, index, idx_yzzz*np, value.elts[idx_yzzz-idx_xxxx]);

    tmp0.store(mask, index, idx_zzxx*np, value.elts[idx_zzzz-idx_xxxx]);
    tmp0.store(mask, index, idx_zzxy*np, value.elts[idx_zzxy-idx_xxxx]);
    tmp0.store(mask, index, idx_zzxz*np, value.elts[idx_zzxz-idx_xxxx]);
    tmp0.store(mask, index, idx_zzyy*np, value.elts[idx_zzyy-idx_xxxx]);
    tmp0.store(mask, index, idx_zzyz*np, value.elts[idx_zzyz-idx_xxxx]);
    tmp0.store(mask, index, idx_zzzz*np, value.elts[idx_zzzz-idx_xxxx]);
  }
};

extern "C" void Z4c_RHS(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_Z4c_RHS;
  DECLARE_CCTK_PARAMETERS;

  for (int d = 0; d < 3; ++d)
    if (cctk_nghostzones[d] < deriv_order / 2 + 1)
      CCTK_VERROR("Need at least %d ghost zones", deriv_order / 2 + 1);

  //

  const array<int, dim> indextype = {0, 0, 0};
  const array<int, dim> nghostzones = {cctk_nghostzones[0], cctk_nghostzones[1],
                                       cctk_nghostzones[2]};
  vect<int, dim> imin, imax;
  GridDescBase(cctkGH).box_int<0, 0, 0>(nghostzones, imin, imax);
  // Suffix 1: with ghost zones, suffix 0: without ghost zones
  const GF3D2layout layout1(cctkGH, indextype);
  const GF3D5layout layout0(imin, imax);

  const GF3D2<const CCTK_REAL> gf_chi1(layout1, chi);
  //const scalar_t<0> gf_chi1_v;

  const smat<GF3D2<const CCTK_REAL>, 3> gf_gammat1{
      GF3D2<const CCTK_REAL>(layout1, gammatxx),
      GF3D2<const CCTK_REAL>(layout1, gammatxy),
      GF3D2<const CCTK_REAL>(layout1, gammatxz),
      GF3D2<const CCTK_REAL>(layout1, gammatyy),
      GF3D2<const CCTK_REAL>(layout1, gammatyz),
      GF3D2<const CCTK_REAL>(layout1, gammatzz)};
  //const mat_t<0> gf_gammat1_v;

  const GF3D2<const CCTK_REAL> gf_Kh1(layout1, Kh);
  //const scalar_t<0> gf_Kh1_v;

  const smat<GF3D2<const CCTK_REAL>, 3> gf_At1{
      GF3D2<const CCTK_REAL>(layout1, Atxx),
      GF3D2<const CCTK_REAL>(layout1, Atxy),
      GF3D2<const CCTK_REAL>(layout1, Atxz),
      GF3D2<const CCTK_REAL>(layout1, Atyy),
      GF3D2<const CCTK_REAL>(layout1, Atyz),
      GF3D2<const CCTK_REAL>(layout1, Atzz)};
  //const mat_t<0> gf_At1_v;

  const vec<GF3D2<const CCTK_REAL>, 3> gf_Gamt1{
      GF3D2<const CCTK_REAL>(layout1, Gamtx),
      GF3D2<const CCTK_REAL>(layout1, Gamty),
      GF3D2<const CCTK_REAL>(layout1, Gamtz)};
  //const vec_t<0> gf_Gamt1_v;

  const GF3D2<const CCTK_REAL> gf_Theta1(layout1, Theta);
  //const scalar_t<0> gf_Theta1_v;

  const GF3D2<const CCTK_REAL> gf_alphaG1(layout1, alphaG);
  //const scalar_t<0> gf_alphaG1_v;

  const vec<GF3D2<const CCTK_REAL>, 3> gf_betaG1{
      GF3D2<const CCTK_REAL>(layout1, betaGx),
      GF3D2<const CCTK_REAL>(layout1, betaGy),
      GF3D2<const CCTK_REAL>(layout1, betaGz)};
  //const vec_t<0> gf_betaG1_v;

  //

  // Ideas:
  //
  // - Outline certain functions, e.g. `det` or `raise_index`. Ensure
  //   they are called with floating-point arguments, not tensor
  //   indices.

  const int ntmps = 154;
  GF3D5vector<CCTK_REAL> tmps(layout0, ntmps);
  int itmp = 0;

  const auto make_gf = [&]() { return GF3D5<CCTK_REAL>(tmps(itmp++)); };
  const auto make_vec = [&](const auto &f) {
    return vec<result_of_t<decltype(f)()>, 3>([&](int) { return f(); });
  };
  const auto make_mat = [&](const auto &f) {
    return smat<result_of_t<decltype(f)()>, 3>([&](int, int) { return f(); });
  };
  const auto make_vec_gf = [&]() { return make_vec(make_gf); };
  const auto make_mat_gf = [&]() { return make_mat(make_gf); };
  const auto make_vec_vec_gf = [&]() { return make_vec(make_vec_gf); };
  const auto make_vec_mat_gf = [&]() { return make_vec(make_mat_gf); };
  const auto make_mat_vec_gf = [&]() { return make_mat(make_vec_gf); };
  const auto make_mat_mat_gf = [&]() { return make_mat(make_mat_gf); };

  const GF3D5<CCTK_REAL> gf_chi0(make_gf());
  const scalar_t<0> gf_chi0_v;
  const vec<GF3D5<CCTK_REAL>, 3> gf_dchi0(make_vec_gf());
  const vec_t<next(gf_chi0_v)> gf_dchi0_v;
  const smat<GF3D5<CCTK_REAL>, 3> gf_ddchi0(make_mat_gf());
  const mat_t<next(gf_dchi0_v)> gf_ddchi0_v;
  calc_derivs2(cctkGH, gf_chi1, gf_chi0, gf_dchi0, gf_ddchi0, layout0);

  const smat<GF3D5<CCTK_REAL>, 3> gf_gammat0(make_mat_gf());
  const mat_t<next(gf_ddchi0_v)> gf_gammat0_v;
  const smat<vec<GF3D5<CCTK_REAL>, 3>, 3> gf_dgammat0(make_mat_vec_gf());
  const mat_vec_t<next(gf_gammat0_v)> gf_dgammat0_v;
  const smat<smat<GF3D5<CCTK_REAL>, 3>, 3> gf_ddgammat0(make_mat_mat_gf());
  const mat_mat_t<next(gf_dgammat0_v)> gf_ddgammat0_v;
  calc_derivs2(cctkGH, gf_gammat1, gf_gammat0, gf_dgammat0, gf_ddgammat0,
               layout0);

  const GF3D5<CCTK_REAL> gf_Kh0(make_gf());
  const scalar_t<next(gf_ddgammat0_v)> gf_Kh0_v;
  const vec<GF3D5<CCTK_REAL>, 3> gf_dKh0(make_vec_gf());
  const vec_t<next(gf_Kh0_v)> gf_dKh0_v;
  calc_derivs(cctkGH, gf_Kh1, gf_Kh0, gf_dKh0, layout0);

  const smat<GF3D5<CCTK_REAL>, 3> gf_At0(make_mat_gf());
  const mat_t<next(gf_dKh0_v)> gf_At0_v;
  const smat<vec<GF3D5<CCTK_REAL>, 3>, 3> gf_dAt0(make_mat_vec_gf());
  const mat_vec_t<next(gf_At0_v)> gf_dAt0_v;
  calc_derivs(cctkGH, gf_At1, gf_At0, gf_dAt0, layout0);

  const vec<GF3D5<CCTK_REAL>, 3> gf_Gamt0(make_vec_gf());
  const vec_t<next(gf_dAt0_v)> gf_Gamt0_v;
  const vec<vec<GF3D5<CCTK_REAL>, 3>, 3> gf_dGamt0(make_vec_vec_gf());
  const vec_vec_t<next(gf_Gamt0_v)> gf_dGamt0_v;
  calc_derivs(cctkGH, gf_Gamt1, gf_Gamt0, gf_dGamt0, layout0);

  const GF3D5<CCTK_REAL> gf_Theta0(make_gf());
  const scalar_t<next(gf_dGamt0_v)> gf_Theta0_v;
  const vec<GF3D5<CCTK_REAL>, 3> gf_dTheta0(make_vec_gf());
  const vec_t<next(gf_Theta0_v)> gf_dTheta0_v;
  calc_derivs(cctkGH, gf_Theta1, gf_Theta0, gf_dTheta0, layout0);

  const GF3D5<CCTK_REAL> gf_alphaG0(make_gf());
  const scalar_t<next(gf_dTheta0_v)> gf_alphaG0_v;
  const vec<GF3D5<CCTK_REAL>, 3> gf_dalphaG0(make_vec_gf());
  const vec_t<next(gf_alphaG0_v)> gf_dalphaG0_v;
  const smat<GF3D5<CCTK_REAL>, 3> gf_ddalphaG0(make_mat_gf());
  const mat_t<next(gf_dalphaG0_v)> gf_ddalphaG0_v;
  calc_derivs2(cctkGH, gf_alphaG1, gf_alphaG0, gf_dalphaG0, gf_ddalphaG0,
               layout0);

  const vec<GF3D5<CCTK_REAL>, 3> gf_betaG0(make_vec_gf());
  const vec_t<next(gf_ddalphaG0_v)> gf_betaG0_v;
  const vec<vec<GF3D5<CCTK_REAL>, 3>, 3> gf_dbetaG0(make_vec_vec_gf());
  const vec_vec_t<next(gf_betaG0_v)> gf_dbetaG0_v;
  const vec<smat<GF3D5<CCTK_REAL>, 3>, 3> gf_ddbetaG0(make_vec_mat_gf());
  const vec_mat_t<next(gf_dbetaG0_v)> gf_ddbetaG0_v;
  calc_derivs2(cctkGH, gf_betaG1, gf_betaG0, gf_dbetaG0, gf_ddbetaG0, layout0);

  if (itmp != ntmps)
    CCTK_VERROR("Wrong number of temporary variables: ntmps=%d itmp=%d", ntmps,
                itmp);
  itmp = -1;

  //

  const GF3D2<const CCTK_REAL> gf_eTtt1(layout1, eTtt);
  const scalar_t<0> gf_eTtt1_v;

  const vec<GF3D2<const CCTK_REAL>, 3> gf_eTti1{
      GF3D2<const CCTK_REAL>(layout1, eTtx),
      GF3D2<const CCTK_REAL>(layout1, eTty),
      GF3D2<const CCTK_REAL>(layout1, eTtz)};
  const vec_t<0> gf_eTti1_v;

  const smat<GF3D2<const CCTK_REAL>, 3> gf_eTij1{
      GF3D2<const CCTK_REAL>(layout1, eTxx),
      GF3D2<const CCTK_REAL>(layout1, eTxy),
      GF3D2<const CCTK_REAL>(layout1, eTxz),
      GF3D2<const CCTK_REAL>(layout1, eTyy),
      GF3D2<const CCTK_REAL>(layout1, eTyz),
      GF3D2<const CCTK_REAL>(layout1, eTzz)};
  const mat_t<0> gf_eTij1_v;

  //

  const GF3D2<CCTK_REAL> gf_chi_rhs1(layout1, chi_rhs);
  const scalar_t<0> gf_chi_rhs1_v;

  const smat<GF3D2<CCTK_REAL>, 3> gf_gammat_rhs1{
      GF3D2<CCTK_REAL>(layout1, gammatxx_rhs),
      GF3D2<CCTK_REAL>(layout1, gammatxy_rhs),
      GF3D2<CCTK_REAL>(layout1, gammatxz_rhs),
      GF3D2<CCTK_REAL>(layout1, gammatyy_rhs),
      GF3D2<CCTK_REAL>(layout1, gammatyz_rhs),
      GF3D2<CCTK_REAL>(layout1, gammatzz_rhs)};
  const mat_t<0> gf_gammat_rhs1_v;

  const GF3D2<CCTK_REAL> gf_Kh_rhs1(layout1, Kh_rhs);
  const scalar_t<0> gf_Kh_rhs1_v;

  const smat<GF3D2<CCTK_REAL>, 3> gf_At_rhs1{
      GF3D2<CCTK_REAL>(layout1, Atxx_rhs), GF3D2<CCTK_REAL>(layout1, Atxy_rhs),
      GF3D2<CCTK_REAL>(layout1, Atxz_rhs), GF3D2<CCTK_REAL>(layout1, Atyy_rhs),
      GF3D2<CCTK_REAL>(layout1, Atyz_rhs), GF3D2<CCTK_REAL>(layout1, Atzz_rhs)};
  const mat_t<0> gf_At_rhs1_v;

  const vec<GF3D2<CCTK_REAL>, 3> gf_Gamt_rhs1{
      GF3D2<CCTK_REAL>(layout1, Gamtx_rhs),
      GF3D2<CCTK_REAL>(layout1, Gamty_rhs),
      GF3D2<CCTK_REAL>(layout1, Gamtz_rhs)};
  const vec_t<0> gf_Gamt_rhs1_v;

  const GF3D2<CCTK_REAL> gf_Theta_rhs1(layout1, Theta_rhs);
  const scalar_t<0> gf_Theta_rhs1_v;

  const GF3D2<CCTK_REAL> gf_alphaG_rhs1(layout1, alphaG_rhs);
  const scalar_t<0> gf_alphaG_rhs1_v;

  const vec<GF3D2<CCTK_REAL>, 3> gf_betaG_rhs1{
      GF3D2<CCTK_REAL>(layout1, betaGx_rhs),
      GF3D2<CCTK_REAL>(layout1, betaGy_rhs),
      GF3D2<CCTK_REAL>(layout1, betaGz_rhs)};
  const vec_t<0> gf_betaG_rhs1_v;

  //

  const Loop::GridDescBaseDevice grid(cctkGH);

#if 1

#ifdef __CUDACC__
  const nvtxRangeId_t range = nvtxRangeStartA("Z4c_RHS::rhs");
#endif
  noinline([&]() __attribute__((__flatten__, __hot__)) {
    grid.loop_int_device<0, 0, 0, vsize>(
        grid.nghostzones, [=] ARITH_DEVICE(const PointDesc &p) ARITH_INLINE {
          const vbool mask = mask_for_loop_tail<vbool>(p.i, p.imax);
          const GF3D2index index1(layout1, p.I);
          const GF3D5index index0(layout0, p.I);
	  const int n0 = layout0.np;
	  const int n1 = layout1.np;

          // Load and calculate
          const z4c_vars<vreal> vars(
              set_Theta_zero, kappa1, kappa2, f_mu_L, f_mu_S, eta, //
              gf_chi0_v(mask, index0, gf_chi0, n0), gf_dchi0_v(mask, index0, gf_chi0, n0),
              gf_ddchi0_v(mask, index0, gf_chi0, n0), //
              gf_gammat0_v(mask, index0, gf_chi0, n0), gf_dgammat0_v(mask, index0, gf_chi0, n0),
              gf_ddgammat0_v(mask, index0, gf_chi0, n0),                        //
              gf_Kh0_v(mask, index0, gf_chi0, n0), gf_dKh0_v(mask, index0, gf_chi0, n0),       //
              gf_At0_v(mask, index0, gf_chi0, n0), gf_dAt0_v(mask, index0, gf_chi0, n0),       //
              gf_Gamt0_v(mask, index0, gf_chi0, n0), gf_dGamt0_v(mask, index0, gf_chi0, n0),   //
              gf_Theta0_v(mask, index0, gf_chi0, n0), gf_dTheta0_v(mask, index0, gf_chi0, n0), //
              gf_alphaG0_v(mask, index0, gf_chi0, n0), gf_dalphaG0_v(mask, index0, gf_chi0, n0),
              gf_ddalphaG0_v(mask, index0, gf_chi0, n0), //
              gf_betaG0_v(mask, index0, gf_chi0, n0), gf_dbetaG0_v(mask, index0, gf_chi0, n0),
              gf_ddbetaG0_v(mask, index0, gf_chi0, n0), //
              gf_eTtt1_v(mask, index1, gf_eTtt1, n1), gf_eTti1_v(mask, index1, gf_eTti1(0), n1),
              gf_eTij1_v(mask, index1, gf_eTij1(0,0), n1));

          gf_chi_rhs1_v.store(mask, index1, gf_chi_rhs1, n1, vars.chi_rhs);
          gf_gammat_rhs1_v.store(mask, index1, gf_gammat_rhs1(0,0), n1, vars.gammat_rhs);
          gf_Kh_rhs1_v.store(mask, index1, gf_Kh_rhs1, n1, vars.Kh_rhs);
          gf_At_rhs1_v.store(mask, index1, gf_At_rhs1(0,0), n1, vars.At_rhs);
          gf_Gamt_rhs1_v.store(mask, index1, gf_Gamt_rhs1(0), n1, vars.Gamt_rhs);
          gf_Theta_rhs1_v.store(mask, index1, gf_Theta_rhs1, n1, vars.Theta_rhs);
          gf_alphaG_rhs1_v.store(mask, index1, gf_alphaG_rhs1, n1, vars.alphaG_rhs);
          gf_betaG_rhs1_v.store(mask, index1, gf_betaG_rhs1(0), n1, vars.betaG_rhs);
        });
  });
#ifdef __CUDACC__
  nvtxRangeEnd(range);
#endif

#else

  noinline([&]() __attribute__((__flatten__, __hot__)) {
    grid.loop_int_device<0, 0, 0, vsize>(
        grid.nghostzones, [=] ARITH_DEVICE(const PointDesc &p) ARITH_INLINE {
          const vbool mask = mask_for_loop_tail<vbool>(p.i, p.imax);
          const GF3D2index index1(layout1, p.I);
          const GF3D5index index0(layout0, p.I);

          // Load and calculate
          const z4c_vars<vreal> vars(
              kappa1, kappa2, f_mu_L, f_mu_S, eta, //
              gf_chi0(mask, index0), gf_dchi0(mask, index0),
              gf_ddchi0(mask, index0), //
              gf_gammat0(mask, index0), gf_dgammat0(mask, index0),
              gf_ddgammat0(mask, index0),                        //
              gf_Kh0(mask, index0), gf_dKh0(mask, index0),       //
              gf_At0(mask, index0), gf_dAt0(mask, index0),       //
              gf_Gamt0(mask, index0), gf_dGamt0(mask, index0),   //
              gf_Theta0(mask, index0), gf_dTheta0(mask, index0), //
              gf_alphaG0(mask, index0), gf_dalphaG0(mask, index0),
              gf_ddalphaG0(mask, index0), //
              gf_betaG0(mask, index0), gf_dbetaG0(mask, index0),
              gf_ddbetaG0(mask, index0), //
              gf_eTtt1(mask, index1), gf_eTti1(mask, index1),
              gf_eTij1(mask, index1));

          // Store Kh_rhs, At_rhs, Gamt_rhs, Theta_rhs
          gf_Kh_rhs1.store(mask, index1, vars.Kh_rhs);
          gf_At_rhs1.store(mask, index1, vars.At_rhs);
          gf_Gamt_rhs1.store(mask, index1, vars.Gamt_rhs);
          gf_Theta_rhs1.store(mask, index1, vars.Theta_rhs);
        });
  });

  noinline([&]() __attribute__((__flatten__, __hot__)) {
    grid.loop_int_device<0, 0, 0, vsize>(
        grid.nghostzones, [=] ARITH_DEVICE(const PointDesc &p) ARITH_INLINE {
          const vbool mask = mask_for_loop_tail<vbool>(p.i, p.imax);
          const GF3D2index index1(layout1, p.I);
          const GF3D5index index0(layout0, p.I);

          // Load and calculate
          const z4c_vars<vreal> vars(
              kappa1, kappa2, f_mu_L, f_mu_S, eta, //
              gf_chi0(mask, index0), gf_dchi0(mask, index0),
              gf_ddchi0(mask, index0), //
              gf_gammat0(mask, index0), gf_dgammat0(mask, index0),
              gf_ddgammat0(mask, index0),                        //
              gf_Kh0(mask, index0), gf_dKh0(mask, index0),       //
              gf_At0(mask, index0), gf_dAt0(mask, index0),       //
              gf_Gamt0(mask, index0), gf_dGamt0(mask, index0),   //
              gf_Theta0(mask, index0), gf_dTheta0(mask, index0), //
              gf_alphaG0(mask, index0), gf_dalphaG0(mask, index0),
              gf_ddalphaG0(mask, index0), //
              gf_betaG0(mask, index0), gf_dbetaG0(mask, index0),
              gf_ddbetaG0(mask, index0), //
              gf_eTtt1(mask, index1), gf_eTti1(mask, index1),
              gf_eTij1(mask, index1));

          // Store chi_rhs, gammat_rhs, alphaG_rhs, betaG_rhs
          gf_chi_rhs1.store(mask, index1, vars.chi_rhs);
          gf_gammat_rhs1.store(mask, index1, vars.gammat_rhs);
          gf_alphaG_rhs1.store(mask, index1, vars.alphaG_rhs);
          gf_betaG_rhs1.store(mask, index1, vars.betaG_rhs);
        });
  });

#endif

  // Upwind and dissipation terms

  // TODO: Consider fusing the loops to reduce memory bandwidth

  apply_upwind_diss(cctkGH, gf_chi1, gf_betaG1, gf_chi_rhs1);

  for (int a = 0; a < 3; ++a)
    for (int b = a; b < 3; ++b)
      apply_upwind_diss(cctkGH, gf_gammat1(a, b), gf_betaG1,
                        gf_gammat_rhs1(a, b));

  apply_upwind_diss(cctkGH, gf_Kh1, gf_betaG1, gf_Kh_rhs1);

  for (int a = 0; a < 3; ++a)
    for (int b = a; b < 3; ++b)
      apply_upwind_diss(cctkGH, gf_At1(a, b), gf_betaG1, gf_At_rhs1(a, b));

  for (int a = 0; a < 3; ++a)
    apply_upwind_diss(cctkGH, gf_Gamt1(a), gf_betaG1, gf_Gamt_rhs1(a));

  if (!set_Theta_zero)
    apply_upwind_diss(cctkGH, gf_Theta1, gf_betaG1, gf_Theta_rhs1);

  apply_upwind_diss(cctkGH, gf_alphaG1, gf_betaG1, gf_alphaG_rhs1);

  for (int a = 0; a < 3; ++a)
    apply_upwind_diss(cctkGH, gf_betaG1(a), gf_betaG1, gf_betaG_rhs1(a));
}

} // namespace Z4c
