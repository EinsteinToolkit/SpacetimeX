#ifndef KRANC_DIFF_FUNCTIONS
#  define PDstandard2nd1(u) ((-KRANC_GFOFFSET3D(u,-1,0,0) + KRANC_GFOFFSET3D(u,1,0,0))*p1o2dx)
#else
#  define PDstandard2nd1(u) (PDstandard2nd1_impl(u,p1o2dx,cdj,cdk))
static CCTK_REAL PDstandard2nd1_impl(const CCTK_REAL* restrict const u, const CCTK_REAL p1o2dx, const ptrdiff_t cdj, const ptrdiff_t cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL PDstandard2nd1_impl(const CCTK_REAL* restrict const u, const CCTK_REAL p1o2dx, const ptrdiff_t cdj, const ptrdiff_t cdk)
{
  const ptrdiff_t cdi CCTK_ATTRIBUTE_UNUSED = sizeof(CCTK_REAL);
  return (-KRANC_GFOFFSET3D(u,-1,0,0) + KRANC_GFOFFSET3D(u,1,0,0))*p1o2dx;
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDstandard2nd2(u) ((-KRANC_GFOFFSET3D(u,0,-1,0) + KRANC_GFOFFSET3D(u,0,1,0))*p1o2dy)
#else
#  define PDstandard2nd2(u) (PDstandard2nd2_impl(u,p1o2dy,cdj,cdk))
static CCTK_REAL PDstandard2nd2_impl(const CCTK_REAL* restrict const u, const CCTK_REAL p1o2dy, const ptrdiff_t cdj, const ptrdiff_t cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL PDstandard2nd2_impl(const CCTK_REAL* restrict const u, const CCTK_REAL p1o2dy, const ptrdiff_t cdj, const ptrdiff_t cdk)
{
  const ptrdiff_t cdi CCTK_ATTRIBUTE_UNUSED = sizeof(CCTK_REAL);
  return (-KRANC_GFOFFSET3D(u,0,-1,0) + KRANC_GFOFFSET3D(u,0,1,0))*p1o2dy;
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDstandard2nd3(u) ((-KRANC_GFOFFSET3D(u,0,0,-1) + KRANC_GFOFFSET3D(u,0,0,1))*p1o2dz)
#else
#  define PDstandard2nd3(u) (PDstandard2nd3_impl(u,p1o2dz,cdj,cdk))
static CCTK_REAL PDstandard2nd3_impl(const CCTK_REAL* restrict const u, const CCTK_REAL p1o2dz, const ptrdiff_t cdj, const ptrdiff_t cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL PDstandard2nd3_impl(const CCTK_REAL* restrict const u, const CCTK_REAL p1o2dz, const ptrdiff_t cdj, const ptrdiff_t cdk)
{
  const ptrdiff_t cdi CCTK_ATTRIBUTE_UNUSED = sizeof(CCTK_REAL);
  return PDstandard2nd2_impl(u, p1o2dz, cdk, cdj);
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDstandard2nd11(u) ((-2*KRANC_GFOFFSET3D(u,0,0,0) + KRANC_GFOFFSET3D(u,-1,0,0) + KRANC_GFOFFSET3D(u,1,0,0))*p1odx2)
#else
#  define PDstandard2nd11(u) (PDstandard2nd11_impl(u,p1odx2,cdj,cdk))
static CCTK_REAL PDstandard2nd11_impl(const CCTK_REAL* restrict const u, const CCTK_REAL p1odx2, const ptrdiff_t cdj, const ptrdiff_t cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL PDstandard2nd11_impl(const CCTK_REAL* restrict const u, const CCTK_REAL p1odx2, const ptrdiff_t cdj, const ptrdiff_t cdk)
{
  const ptrdiff_t cdi CCTK_ATTRIBUTE_UNUSED = sizeof(CCTK_REAL);
  return (-2*KRANC_GFOFFSET3D(u,0,0,0) + KRANC_GFOFFSET3D(u,-1,0,0) + KRANC_GFOFFSET3D(u,1,0,0))*p1odx2;
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDstandard2nd22(u) ((-2*KRANC_GFOFFSET3D(u,0,0,0) + KRANC_GFOFFSET3D(u,0,-1,0) + KRANC_GFOFFSET3D(u,0,1,0))*p1ody2)
#else
#  define PDstandard2nd22(u) (PDstandard2nd22_impl(u,p1ody2,cdj,cdk))
static CCTK_REAL PDstandard2nd22_impl(const CCTK_REAL* restrict const u, const CCTK_REAL p1ody2, const ptrdiff_t cdj, const ptrdiff_t cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL PDstandard2nd22_impl(const CCTK_REAL* restrict const u, const CCTK_REAL p1ody2, const ptrdiff_t cdj, const ptrdiff_t cdk)
{
  const ptrdiff_t cdi CCTK_ATTRIBUTE_UNUSED = sizeof(CCTK_REAL);
  return (-2*KRANC_GFOFFSET3D(u,0,0,0) + KRANC_GFOFFSET3D(u,0,-1,0) + KRANC_GFOFFSET3D(u,0,1,0))*p1ody2;
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDstandard2nd33(u) ((-2*KRANC_GFOFFSET3D(u,0,0,0) + KRANC_GFOFFSET3D(u,0,0,-1) + KRANC_GFOFFSET3D(u,0,0,1))*p1odz2)
#else
#  define PDstandard2nd33(u) (PDstandard2nd33_impl(u,p1odz2,cdj,cdk))
static CCTK_REAL PDstandard2nd33_impl(const CCTK_REAL* restrict const u, const CCTK_REAL p1odz2, const ptrdiff_t cdj, const ptrdiff_t cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL PDstandard2nd33_impl(const CCTK_REAL* restrict const u, const CCTK_REAL p1odz2, const ptrdiff_t cdj, const ptrdiff_t cdk)
{
  const ptrdiff_t cdi CCTK_ATTRIBUTE_UNUSED = sizeof(CCTK_REAL);
  return PDstandard2nd22_impl(u, p1odz2, cdk, cdj);
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDstandard2nd12(u) ((KRANC_GFOFFSET3D(u,-1,-1,0) - KRANC_GFOFFSET3D(u,-1,1,0) - KRANC_GFOFFSET3D(u,1,-1,0) + KRANC_GFOFFSET3D(u,1,1,0))*p1o4dxdy)
#else
#  define PDstandard2nd12(u) (PDstandard2nd12_impl(u,p1o4dxdy,cdj,cdk))
static CCTK_REAL PDstandard2nd12_impl(const CCTK_REAL* restrict const u, const CCTK_REAL p1o4dxdy, const ptrdiff_t cdj, const ptrdiff_t cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL PDstandard2nd12_impl(const CCTK_REAL* restrict const u, const CCTK_REAL p1o4dxdy, const ptrdiff_t cdj, const ptrdiff_t cdk)
{
  const ptrdiff_t cdi CCTK_ATTRIBUTE_UNUSED = sizeof(CCTK_REAL);
  return (KRANC_GFOFFSET3D(u,-1,-1,0) - KRANC_GFOFFSET3D(u,-1,1,0) - KRANC_GFOFFSET3D(u,1,-1,0) + KRANC_GFOFFSET3D(u,1,1,0))*p1o4dxdy;
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDstandard2nd13(u) ((KRANC_GFOFFSET3D(u,-1,0,-1) - KRANC_GFOFFSET3D(u,-1,0,1) - KRANC_GFOFFSET3D(u,1,0,-1) + KRANC_GFOFFSET3D(u,1,0,1))*p1o4dxdz)
#else
#  define PDstandard2nd13(u) (PDstandard2nd13_impl(u,p1o4dxdz,cdj,cdk))
static CCTK_REAL PDstandard2nd13_impl(const CCTK_REAL* restrict const u, const CCTK_REAL p1o4dxdz, const ptrdiff_t cdj, const ptrdiff_t cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL PDstandard2nd13_impl(const CCTK_REAL* restrict const u, const CCTK_REAL p1o4dxdz, const ptrdiff_t cdj, const ptrdiff_t cdk)
{
  const ptrdiff_t cdi CCTK_ATTRIBUTE_UNUSED = sizeof(CCTK_REAL);
  return PDstandard2nd12_impl(u, p1o4dxdz, cdk, cdj);
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDstandard2nd21(u) ((KRANC_GFOFFSET3D(u,-1,-1,0) - KRANC_GFOFFSET3D(u,-1,1,0) - KRANC_GFOFFSET3D(u,1,-1,0) + KRANC_GFOFFSET3D(u,1,1,0))*p1o4dxdy)
#else
#  define PDstandard2nd21(u) (PDstandard2nd21_impl(u,p1o4dxdy,cdj,cdk))
static CCTK_REAL PDstandard2nd21_impl(const CCTK_REAL* restrict const u, const CCTK_REAL p1o4dxdy, const ptrdiff_t cdj, const ptrdiff_t cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL PDstandard2nd21_impl(const CCTK_REAL* restrict const u, const CCTK_REAL p1o4dxdy, const ptrdiff_t cdj, const ptrdiff_t cdk)
{
  const ptrdiff_t cdi CCTK_ATTRIBUTE_UNUSED = sizeof(CCTK_REAL);
  return PDstandard2nd12_impl(u, p1o4dxdy, cdj, cdk);
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDstandard2nd23(u) ((KRANC_GFOFFSET3D(u,0,-1,-1) - KRANC_GFOFFSET3D(u,0,-1,1) - KRANC_GFOFFSET3D(u,0,1,-1) + KRANC_GFOFFSET3D(u,0,1,1))*p1o4dydz)
#else
#  define PDstandard2nd23(u) (PDstandard2nd23_impl(u,p1o4dydz,cdj,cdk))
static CCTK_REAL PDstandard2nd23_impl(const CCTK_REAL* restrict const u, const CCTK_REAL p1o4dydz, const ptrdiff_t cdj, const ptrdiff_t cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL PDstandard2nd23_impl(const CCTK_REAL* restrict const u, const CCTK_REAL p1o4dydz, const ptrdiff_t cdj, const ptrdiff_t cdk)
{
  const ptrdiff_t cdi CCTK_ATTRIBUTE_UNUSED = sizeof(CCTK_REAL);
  return (KRANC_GFOFFSET3D(u,0,-1,-1) - KRANC_GFOFFSET3D(u,0,-1,1) - KRANC_GFOFFSET3D(u,0,1,-1) + KRANC_GFOFFSET3D(u,0,1,1))*p1o4dydz;
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDstandard2nd31(u) ((KRANC_GFOFFSET3D(u,-1,0,-1) - KRANC_GFOFFSET3D(u,-1,0,1) - KRANC_GFOFFSET3D(u,1,0,-1) + KRANC_GFOFFSET3D(u,1,0,1))*p1o4dxdz)
#else
#  define PDstandard2nd31(u) (PDstandard2nd31_impl(u,p1o4dxdz,cdj,cdk))
static CCTK_REAL PDstandard2nd31_impl(const CCTK_REAL* restrict const u, const CCTK_REAL p1o4dxdz, const ptrdiff_t cdj, const ptrdiff_t cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL PDstandard2nd31_impl(const CCTK_REAL* restrict const u, const CCTK_REAL p1o4dxdz, const ptrdiff_t cdj, const ptrdiff_t cdk)
{
  const ptrdiff_t cdi CCTK_ATTRIBUTE_UNUSED = sizeof(CCTK_REAL);
  return PDstandard2nd12_impl(u, p1o4dxdz, cdk, cdj);
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDstandard2nd32(u) ((KRANC_GFOFFSET3D(u,0,-1,-1) - KRANC_GFOFFSET3D(u,0,-1,1) - KRANC_GFOFFSET3D(u,0,1,-1) + KRANC_GFOFFSET3D(u,0,1,1))*p1o4dydz)
#else
#  define PDstandard2nd32(u) (PDstandard2nd32_impl(u,p1o4dydz,cdj,cdk))
static CCTK_REAL PDstandard2nd32_impl(const CCTK_REAL* restrict const u, const CCTK_REAL p1o4dydz, const ptrdiff_t cdj, const ptrdiff_t cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL PDstandard2nd32_impl(const CCTK_REAL* restrict const u, const CCTK_REAL p1o4dydz, const ptrdiff_t cdj, const ptrdiff_t cdk)
{
  const ptrdiff_t cdi CCTK_ATTRIBUTE_UNUSED = sizeof(CCTK_REAL);
  return PDstandard2nd23_impl(u, p1o4dydz, cdj, cdk);
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDstandard4th1(u) ((-8*KRANC_GFOFFSET3D(u,-1,0,0) + 8*KRANC_GFOFFSET3D(u,1,0,0) + KRANC_GFOFFSET3D(u,-2,0,0) - KRANC_GFOFFSET3D(u,2,0,0))*p1o12dx)
#else
#  define PDstandard4th1(u) (PDstandard4th1_impl(u,p1o12dx,cdj,cdk))
static CCTK_REAL PDstandard4th1_impl(const CCTK_REAL* restrict const u, const CCTK_REAL p1o12dx, const ptrdiff_t cdj, const ptrdiff_t cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL PDstandard4th1_impl(const CCTK_REAL* restrict const u, const CCTK_REAL p1o12dx, const ptrdiff_t cdj, const ptrdiff_t cdk)
{
  const ptrdiff_t cdi CCTK_ATTRIBUTE_UNUSED = sizeof(CCTK_REAL);
  return (-8*KRANC_GFOFFSET3D(u,-1,0,0) + 8*KRANC_GFOFFSET3D(u,1,0,0) + KRANC_GFOFFSET3D(u,-2,0,0) - KRANC_GFOFFSET3D(u,2,0,0))*p1o12dx;
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDstandard4th2(u) ((-8*KRANC_GFOFFSET3D(u,0,-1,0) + 8*KRANC_GFOFFSET3D(u,0,1,0) + KRANC_GFOFFSET3D(u,0,-2,0) - KRANC_GFOFFSET3D(u,0,2,0))*p1o12dy)
#else
#  define PDstandard4th2(u) (PDstandard4th2_impl(u,p1o12dy,cdj,cdk))
static CCTK_REAL PDstandard4th2_impl(const CCTK_REAL* restrict const u, const CCTK_REAL p1o12dy, const ptrdiff_t cdj, const ptrdiff_t cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL PDstandard4th2_impl(const CCTK_REAL* restrict const u, const CCTK_REAL p1o12dy, const ptrdiff_t cdj, const ptrdiff_t cdk)
{
  const ptrdiff_t cdi CCTK_ATTRIBUTE_UNUSED = sizeof(CCTK_REAL);
  return (-8*KRANC_GFOFFSET3D(u,0,-1,0) + 8*KRANC_GFOFFSET3D(u,0,1,0) + KRANC_GFOFFSET3D(u,0,-2,0) - KRANC_GFOFFSET3D(u,0,2,0))*p1o12dy;
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDstandard4th3(u) ((-8*KRANC_GFOFFSET3D(u,0,0,-1) + 8*KRANC_GFOFFSET3D(u,0,0,1) + KRANC_GFOFFSET3D(u,0,0,-2) - KRANC_GFOFFSET3D(u,0,0,2))*p1o12dz)
#else
#  define PDstandard4th3(u) (PDstandard4th3_impl(u,p1o12dz,cdj,cdk))
static CCTK_REAL PDstandard4th3_impl(const CCTK_REAL* restrict const u, const CCTK_REAL p1o12dz, const ptrdiff_t cdj, const ptrdiff_t cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL PDstandard4th3_impl(const CCTK_REAL* restrict const u, const CCTK_REAL p1o12dz, const ptrdiff_t cdj, const ptrdiff_t cdk)
{
  const ptrdiff_t cdi CCTK_ATTRIBUTE_UNUSED = sizeof(CCTK_REAL);
  return PDstandard4th2_impl(u, p1o12dz, cdk, cdj);
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDstandard4th11(u) ((30*KRANC_GFOFFSET3D(u,0,0,0) - 16*(KRANC_GFOFFSET3D(u,-1,0,0) + KRANC_GFOFFSET3D(u,1,0,0)) + KRANC_GFOFFSET3D(u,-2,0,0) + KRANC_GFOFFSET3D(u,2,0,0))*pm1o12dx2)
#else
#  define PDstandard4th11(u) (PDstandard4th11_impl(u,pm1o12dx2,cdj,cdk))
static CCTK_REAL PDstandard4th11_impl(const CCTK_REAL* restrict const u, const CCTK_REAL pm1o12dx2, const ptrdiff_t cdj, const ptrdiff_t cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL PDstandard4th11_impl(const CCTK_REAL* restrict const u, const CCTK_REAL pm1o12dx2, const ptrdiff_t cdj, const ptrdiff_t cdk)
{
  const ptrdiff_t cdi CCTK_ATTRIBUTE_UNUSED = sizeof(CCTK_REAL);
  return (30*KRANC_GFOFFSET3D(u,0,0,0) - 16*(KRANC_GFOFFSET3D(u,-1,0,0) + KRANC_GFOFFSET3D(u,1,0,0)) + KRANC_GFOFFSET3D(u,-2,0,0) + KRANC_GFOFFSET3D(u,2,0,0))*pm1o12dx2;
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDstandard4th22(u) ((30*KRANC_GFOFFSET3D(u,0,0,0) - 16*(KRANC_GFOFFSET3D(u,0,-1,0) + KRANC_GFOFFSET3D(u,0,1,0)) + KRANC_GFOFFSET3D(u,0,-2,0) + KRANC_GFOFFSET3D(u,0,2,0))*pm1o12dy2)
#else
#  define PDstandard4th22(u) (PDstandard4th22_impl(u,pm1o12dy2,cdj,cdk))
static CCTK_REAL PDstandard4th22_impl(const CCTK_REAL* restrict const u, const CCTK_REAL pm1o12dy2, const ptrdiff_t cdj, const ptrdiff_t cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL PDstandard4th22_impl(const CCTK_REAL* restrict const u, const CCTK_REAL pm1o12dy2, const ptrdiff_t cdj, const ptrdiff_t cdk)
{
  const ptrdiff_t cdi CCTK_ATTRIBUTE_UNUSED = sizeof(CCTK_REAL);
  return (30*KRANC_GFOFFSET3D(u,0,0,0) - 16*(KRANC_GFOFFSET3D(u,0,-1,0) + KRANC_GFOFFSET3D(u,0,1,0)) + KRANC_GFOFFSET3D(u,0,-2,0) + KRANC_GFOFFSET3D(u,0,2,0))*pm1o12dy2;
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDstandard4th33(u) ((30*KRANC_GFOFFSET3D(u,0,0,0) - 16*(KRANC_GFOFFSET3D(u,0,0,-1) + KRANC_GFOFFSET3D(u,0,0,1)) + KRANC_GFOFFSET3D(u,0,0,-2) + KRANC_GFOFFSET3D(u,0,0,2))*pm1o12dz2)
#else
#  define PDstandard4th33(u) (PDstandard4th33_impl(u,pm1o12dz2,cdj,cdk))
static CCTK_REAL PDstandard4th33_impl(const CCTK_REAL* restrict const u, const CCTK_REAL pm1o12dz2, const ptrdiff_t cdj, const ptrdiff_t cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL PDstandard4th33_impl(const CCTK_REAL* restrict const u, const CCTK_REAL pm1o12dz2, const ptrdiff_t cdj, const ptrdiff_t cdk)
{
  const ptrdiff_t cdi CCTK_ATTRIBUTE_UNUSED = sizeof(CCTK_REAL);
  return PDstandard4th22_impl(u, pm1o12dz2, cdk, cdj);
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDstandard4th12(u) ((-64*(KRANC_GFOFFSET3D(u,-1,1,0) + KRANC_GFOFFSET3D(u,1,-1,0)) + 64*(KRANC_GFOFFSET3D(u,-1,-1,0) + KRANC_GFOFFSET3D(u,1,1,0)) + 8*(KRANC_GFOFFSET3D(u,-1,2,0) + KRANC_GFOFFSET3D(u,1,-2,0) + KRANC_GFOFFSET3D(u,-2,1,0) + KRANC_GFOFFSET3D(u,2,-1,0)) - 8*(KRANC_GFOFFSET3D(u,-1,-2,0) + KRANC_GFOFFSET3D(u,1,2,0) + KRANC_GFOFFSET3D(u,-2,-1,0) + KRANC_GFOFFSET3D(u,2,1,0)) + KRANC_GFOFFSET3D(u,-2,-2,0) - KRANC_GFOFFSET3D(u,-2,2,0) - KRANC_GFOFFSET3D(u,2,-2,0) + KRANC_GFOFFSET3D(u,2,2,0))*p1o144dxdy)
#else
#  define PDstandard4th12(u) (PDstandard4th12_impl(u,p1o144dxdy,cdj,cdk))
static CCTK_REAL PDstandard4th12_impl(const CCTK_REAL* restrict const u, const CCTK_REAL p1o144dxdy, const ptrdiff_t cdj, const ptrdiff_t cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL PDstandard4th12_impl(const CCTK_REAL* restrict const u, const CCTK_REAL p1o144dxdy, const ptrdiff_t cdj, const ptrdiff_t cdk)
{
  const ptrdiff_t cdi CCTK_ATTRIBUTE_UNUSED = sizeof(CCTK_REAL);
  return (-64*(KRANC_GFOFFSET3D(u,-1,1,0) + KRANC_GFOFFSET3D(u,1,-1,0)) + 64*(KRANC_GFOFFSET3D(u,-1,-1,0) + KRANC_GFOFFSET3D(u,1,1,0)) + 8*(KRANC_GFOFFSET3D(u,-1,2,0) + KRANC_GFOFFSET3D(u,1,-2,0) + KRANC_GFOFFSET3D(u,-2,1,0) + KRANC_GFOFFSET3D(u,2,-1,0)) - 8*(KRANC_GFOFFSET3D(u,-1,-2,0) + KRANC_GFOFFSET3D(u,1,2,0) + KRANC_GFOFFSET3D(u,-2,-1,0) + KRANC_GFOFFSET3D(u,2,1,0)) + KRANC_GFOFFSET3D(u,-2,-2,0) - KRANC_GFOFFSET3D(u,-2,2,0) - KRANC_GFOFFSET3D(u,2,-2,0) + KRANC_GFOFFSET3D(u,2,2,0))*p1o144dxdy;
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDstandard4th13(u) ((-64*(KRANC_GFOFFSET3D(u,-1,0,1) + KRANC_GFOFFSET3D(u,1,0,-1)) + 64*(KRANC_GFOFFSET3D(u,-1,0,-1) + KRANC_GFOFFSET3D(u,1,0,1)) + 8*(KRANC_GFOFFSET3D(u,-1,0,2) + KRANC_GFOFFSET3D(u,1,0,-2) + KRANC_GFOFFSET3D(u,-2,0,1) + KRANC_GFOFFSET3D(u,2,0,-1)) - 8*(KRANC_GFOFFSET3D(u,-1,0,-2) + KRANC_GFOFFSET3D(u,1,0,2) + KRANC_GFOFFSET3D(u,-2,0,-1) + KRANC_GFOFFSET3D(u,2,0,1)) + KRANC_GFOFFSET3D(u,-2,0,-2) - KRANC_GFOFFSET3D(u,-2,0,2) - KRANC_GFOFFSET3D(u,2,0,-2) + KRANC_GFOFFSET3D(u,2,0,2))*p1o144dxdz)
#else
#  define PDstandard4th13(u) (PDstandard4th13_impl(u,p1o144dxdz,cdj,cdk))
static CCTK_REAL PDstandard4th13_impl(const CCTK_REAL* restrict const u, const CCTK_REAL p1o144dxdz, const ptrdiff_t cdj, const ptrdiff_t cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL PDstandard4th13_impl(const CCTK_REAL* restrict const u, const CCTK_REAL p1o144dxdz, const ptrdiff_t cdj, const ptrdiff_t cdk)
{
  const ptrdiff_t cdi CCTK_ATTRIBUTE_UNUSED = sizeof(CCTK_REAL);
  return PDstandard4th12_impl(u, p1o144dxdz, cdk, cdj);
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDstandard4th21(u) ((-64*(KRANC_GFOFFSET3D(u,-1,1,0) + KRANC_GFOFFSET3D(u,1,-1,0)) + 64*(KRANC_GFOFFSET3D(u,-1,-1,0) + KRANC_GFOFFSET3D(u,1,1,0)) + 8*(KRANC_GFOFFSET3D(u,-1,2,0) + KRANC_GFOFFSET3D(u,1,-2,0) + KRANC_GFOFFSET3D(u,-2,1,0) + KRANC_GFOFFSET3D(u,2,-1,0)) - 8*(KRANC_GFOFFSET3D(u,-1,-2,0) + KRANC_GFOFFSET3D(u,1,2,0) + KRANC_GFOFFSET3D(u,-2,-1,0) + KRANC_GFOFFSET3D(u,2,1,0)) + KRANC_GFOFFSET3D(u,-2,-2,0) - KRANC_GFOFFSET3D(u,-2,2,0) - KRANC_GFOFFSET3D(u,2,-2,0) + KRANC_GFOFFSET3D(u,2,2,0))*p1o144dxdy)
#else
#  define PDstandard4th21(u) (PDstandard4th21_impl(u,p1o144dxdy,cdj,cdk))
static CCTK_REAL PDstandard4th21_impl(const CCTK_REAL* restrict const u, const CCTK_REAL p1o144dxdy, const ptrdiff_t cdj, const ptrdiff_t cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL PDstandard4th21_impl(const CCTK_REAL* restrict const u, const CCTK_REAL p1o144dxdy, const ptrdiff_t cdj, const ptrdiff_t cdk)
{
  const ptrdiff_t cdi CCTK_ATTRIBUTE_UNUSED = sizeof(CCTK_REAL);
  return PDstandard4th12_impl(u, p1o144dxdy, cdj, cdk);
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDstandard4th23(u) ((-64*(KRANC_GFOFFSET3D(u,0,-1,1) + KRANC_GFOFFSET3D(u,0,1,-1)) + 64*(KRANC_GFOFFSET3D(u,0,-1,-1) + KRANC_GFOFFSET3D(u,0,1,1)) + 8*(KRANC_GFOFFSET3D(u,0,-1,2) + KRANC_GFOFFSET3D(u,0,1,-2) + KRANC_GFOFFSET3D(u,0,-2,1) + KRANC_GFOFFSET3D(u,0,2,-1)) - 8*(KRANC_GFOFFSET3D(u,0,-1,-2) + KRANC_GFOFFSET3D(u,0,1,2) + KRANC_GFOFFSET3D(u,0,-2,-1) + KRANC_GFOFFSET3D(u,0,2,1)) + KRANC_GFOFFSET3D(u,0,-2,-2) - KRANC_GFOFFSET3D(u,0,-2,2) - KRANC_GFOFFSET3D(u,0,2,-2) + KRANC_GFOFFSET3D(u,0,2,2))*p1o144dydz)
#else
#  define PDstandard4th23(u) (PDstandard4th23_impl(u,p1o144dydz,cdj,cdk))
static CCTK_REAL PDstandard4th23_impl(const CCTK_REAL* restrict const u, const CCTK_REAL p1o144dydz, const ptrdiff_t cdj, const ptrdiff_t cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL PDstandard4th23_impl(const CCTK_REAL* restrict const u, const CCTK_REAL p1o144dydz, const ptrdiff_t cdj, const ptrdiff_t cdk)
{
  const ptrdiff_t cdi CCTK_ATTRIBUTE_UNUSED = sizeof(CCTK_REAL);
  return (-64*(KRANC_GFOFFSET3D(u,0,-1,1) + KRANC_GFOFFSET3D(u,0,1,-1)) + 64*(KRANC_GFOFFSET3D(u,0,-1,-1) + KRANC_GFOFFSET3D(u,0,1,1)) + 8*(KRANC_GFOFFSET3D(u,0,-1,2) + KRANC_GFOFFSET3D(u,0,1,-2) + KRANC_GFOFFSET3D(u,0,-2,1) + KRANC_GFOFFSET3D(u,0,2,-1)) - 8*(KRANC_GFOFFSET3D(u,0,-1,-2) + KRANC_GFOFFSET3D(u,0,1,2) + KRANC_GFOFFSET3D(u,0,-2,-1) + KRANC_GFOFFSET3D(u,0,2,1)) + KRANC_GFOFFSET3D(u,0,-2,-2) - KRANC_GFOFFSET3D(u,0,-2,2) - KRANC_GFOFFSET3D(u,0,2,-2) + KRANC_GFOFFSET3D(u,0,2,2))*p1o144dydz;
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDstandard4th31(u) ((-64*(KRANC_GFOFFSET3D(u,-1,0,1) + KRANC_GFOFFSET3D(u,1,0,-1)) + 64*(KRANC_GFOFFSET3D(u,-1,0,-1) + KRANC_GFOFFSET3D(u,1,0,1)) + 8*(KRANC_GFOFFSET3D(u,-1,0,2) + KRANC_GFOFFSET3D(u,1,0,-2) + KRANC_GFOFFSET3D(u,-2,0,1) + KRANC_GFOFFSET3D(u,2,0,-1)) - 8*(KRANC_GFOFFSET3D(u,-1,0,-2) + KRANC_GFOFFSET3D(u,1,0,2) + KRANC_GFOFFSET3D(u,-2,0,-1) + KRANC_GFOFFSET3D(u,2,0,1)) + KRANC_GFOFFSET3D(u,-2,0,-2) - KRANC_GFOFFSET3D(u,-2,0,2) - KRANC_GFOFFSET3D(u,2,0,-2) + KRANC_GFOFFSET3D(u,2,0,2))*p1o144dxdz)
#else
#  define PDstandard4th31(u) (PDstandard4th31_impl(u,p1o144dxdz,cdj,cdk))
static CCTK_REAL PDstandard4th31_impl(const CCTK_REAL* restrict const u, const CCTK_REAL p1o144dxdz, const ptrdiff_t cdj, const ptrdiff_t cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL PDstandard4th31_impl(const CCTK_REAL* restrict const u, const CCTK_REAL p1o144dxdz, const ptrdiff_t cdj, const ptrdiff_t cdk)
{
  const ptrdiff_t cdi CCTK_ATTRIBUTE_UNUSED = sizeof(CCTK_REAL);
  return PDstandard4th12_impl(u, p1o144dxdz, cdk, cdj);
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDstandard4th32(u) ((-64*(KRANC_GFOFFSET3D(u,0,-1,1) + KRANC_GFOFFSET3D(u,0,1,-1)) + 64*(KRANC_GFOFFSET3D(u,0,-1,-1) + KRANC_GFOFFSET3D(u,0,1,1)) + 8*(KRANC_GFOFFSET3D(u,0,-1,2) + KRANC_GFOFFSET3D(u,0,1,-2) + KRANC_GFOFFSET3D(u,0,-2,1) + KRANC_GFOFFSET3D(u,0,2,-1)) - 8*(KRANC_GFOFFSET3D(u,0,-1,-2) + KRANC_GFOFFSET3D(u,0,1,2) + KRANC_GFOFFSET3D(u,0,-2,-1) + KRANC_GFOFFSET3D(u,0,2,1)) + KRANC_GFOFFSET3D(u,0,-2,-2) - KRANC_GFOFFSET3D(u,0,-2,2) - KRANC_GFOFFSET3D(u,0,2,-2) + KRANC_GFOFFSET3D(u,0,2,2))*p1o144dydz)
#else
#  define PDstandard4th32(u) (PDstandard4th32_impl(u,p1o144dydz,cdj,cdk))
static CCTK_REAL PDstandard4th32_impl(const CCTK_REAL* restrict const u, const CCTK_REAL p1o144dydz, const ptrdiff_t cdj, const ptrdiff_t cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL PDstandard4th32_impl(const CCTK_REAL* restrict const u, const CCTK_REAL p1o144dydz, const ptrdiff_t cdj, const ptrdiff_t cdk)
{
  const ptrdiff_t cdi CCTK_ATTRIBUTE_UNUSED = sizeof(CCTK_REAL);
  return PDstandard4th23_impl(u, p1o144dydz, cdj, cdk);
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDstandardfdOrder21(u) ((-KRANC_GFOFFSET3D(u,-1,0,0) + KRANC_GFOFFSET3D(u,1,0,0))*p1o2dx)
#else
#  define PDstandardfdOrder21(u) (PDstandardfdOrder21_impl(u,p1o2dx,cdj,cdk))
static CCTK_REAL PDstandardfdOrder21_impl(const CCTK_REAL* restrict const u, const CCTK_REAL p1o2dx, const ptrdiff_t cdj, const ptrdiff_t cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL PDstandardfdOrder21_impl(const CCTK_REAL* restrict const u, const CCTK_REAL p1o2dx, const ptrdiff_t cdj, const ptrdiff_t cdk)
{
  const ptrdiff_t cdi CCTK_ATTRIBUTE_UNUSED = sizeof(CCTK_REAL);
  return (-KRANC_GFOFFSET3D(u,-1,0,0) + KRANC_GFOFFSET3D(u,1,0,0))*p1o2dx;
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDstandardfdOrder22(u) ((-KRANC_GFOFFSET3D(u,0,-1,0) + KRANC_GFOFFSET3D(u,0,1,0))*p1o2dy)
#else
#  define PDstandardfdOrder22(u) (PDstandardfdOrder22_impl(u,p1o2dy,cdj,cdk))
static CCTK_REAL PDstandardfdOrder22_impl(const CCTK_REAL* restrict const u, const CCTK_REAL p1o2dy, const ptrdiff_t cdj, const ptrdiff_t cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL PDstandardfdOrder22_impl(const CCTK_REAL* restrict const u, const CCTK_REAL p1o2dy, const ptrdiff_t cdj, const ptrdiff_t cdk)
{
  const ptrdiff_t cdi CCTK_ATTRIBUTE_UNUSED = sizeof(CCTK_REAL);
  return (-KRANC_GFOFFSET3D(u,0,-1,0) + KRANC_GFOFFSET3D(u,0,1,0))*p1o2dy;
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDstandardfdOrder23(u) ((-KRANC_GFOFFSET3D(u,0,0,-1) + KRANC_GFOFFSET3D(u,0,0,1))*p1o2dz)
#else
#  define PDstandardfdOrder23(u) (PDstandardfdOrder23_impl(u,p1o2dz,cdj,cdk))
static CCTK_REAL PDstandardfdOrder23_impl(const CCTK_REAL* restrict const u, const CCTK_REAL p1o2dz, const ptrdiff_t cdj, const ptrdiff_t cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL PDstandardfdOrder23_impl(const CCTK_REAL* restrict const u, const CCTK_REAL p1o2dz, const ptrdiff_t cdj, const ptrdiff_t cdk)
{
  const ptrdiff_t cdi CCTK_ATTRIBUTE_UNUSED = sizeof(CCTK_REAL);
  return PDstandardfdOrder22_impl(u, p1o2dz, cdk, cdj);
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDstandardfdOrder41(u) ((-8*KRANC_GFOFFSET3D(u,-1,0,0) + 8*KRANC_GFOFFSET3D(u,1,0,0) + KRANC_GFOFFSET3D(u,-2,0,0) - KRANC_GFOFFSET3D(u,2,0,0))*p1o12dx)
#else
#  define PDstandardfdOrder41(u) (PDstandardfdOrder41_impl(u,p1o12dx,cdj,cdk))
static CCTK_REAL PDstandardfdOrder41_impl(const CCTK_REAL* restrict const u, const CCTK_REAL p1o12dx, const ptrdiff_t cdj, const ptrdiff_t cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL PDstandardfdOrder41_impl(const CCTK_REAL* restrict const u, const CCTK_REAL p1o12dx, const ptrdiff_t cdj, const ptrdiff_t cdk)
{
  const ptrdiff_t cdi CCTK_ATTRIBUTE_UNUSED = sizeof(CCTK_REAL);
  return (-8*KRANC_GFOFFSET3D(u,-1,0,0) + 8*KRANC_GFOFFSET3D(u,1,0,0) + KRANC_GFOFFSET3D(u,-2,0,0) - KRANC_GFOFFSET3D(u,2,0,0))*p1o12dx;
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDstandardfdOrder42(u) ((-8*KRANC_GFOFFSET3D(u,0,-1,0) + 8*KRANC_GFOFFSET3D(u,0,1,0) + KRANC_GFOFFSET3D(u,0,-2,0) - KRANC_GFOFFSET3D(u,0,2,0))*p1o12dy)
#else
#  define PDstandardfdOrder42(u) (PDstandardfdOrder42_impl(u,p1o12dy,cdj,cdk))
static CCTK_REAL PDstandardfdOrder42_impl(const CCTK_REAL* restrict const u, const CCTK_REAL p1o12dy, const ptrdiff_t cdj, const ptrdiff_t cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL PDstandardfdOrder42_impl(const CCTK_REAL* restrict const u, const CCTK_REAL p1o12dy, const ptrdiff_t cdj, const ptrdiff_t cdk)
{
  const ptrdiff_t cdi CCTK_ATTRIBUTE_UNUSED = sizeof(CCTK_REAL);
  return (-8*KRANC_GFOFFSET3D(u,0,-1,0) + 8*KRANC_GFOFFSET3D(u,0,1,0) + KRANC_GFOFFSET3D(u,0,-2,0) - KRANC_GFOFFSET3D(u,0,2,0))*p1o12dy;
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDstandardfdOrder43(u) ((-8*KRANC_GFOFFSET3D(u,0,0,-1) + 8*KRANC_GFOFFSET3D(u,0,0,1) + KRANC_GFOFFSET3D(u,0,0,-2) - KRANC_GFOFFSET3D(u,0,0,2))*p1o12dz)
#else
#  define PDstandardfdOrder43(u) (PDstandardfdOrder43_impl(u,p1o12dz,cdj,cdk))
static CCTK_REAL PDstandardfdOrder43_impl(const CCTK_REAL* restrict const u, const CCTK_REAL p1o12dz, const ptrdiff_t cdj, const ptrdiff_t cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL PDstandardfdOrder43_impl(const CCTK_REAL* restrict const u, const CCTK_REAL p1o12dz, const ptrdiff_t cdj, const ptrdiff_t cdk)
{
  const ptrdiff_t cdi CCTK_ATTRIBUTE_UNUSED = sizeof(CCTK_REAL);
  return PDstandardfdOrder42_impl(u, p1o12dz, cdk, cdj);
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDstandardfdOrder61(u) ((-45*KRANC_GFOFFSET3D(u,-1,0,0) + 45*KRANC_GFOFFSET3D(u,1,0,0) + 9*KRANC_GFOFFSET3D(u,-2,0,0) - 9*KRANC_GFOFFSET3D(u,2,0,0) - KRANC_GFOFFSET3D(u,-3,0,0) + KRANC_GFOFFSET3D(u,3,0,0))*p1o60dx)
#else
#  define PDstandardfdOrder61(u) (PDstandardfdOrder61_impl(u,p1o60dx,cdj,cdk))
static CCTK_REAL PDstandardfdOrder61_impl(const CCTK_REAL* restrict const u, const CCTK_REAL p1o60dx, const ptrdiff_t cdj, const ptrdiff_t cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL PDstandardfdOrder61_impl(const CCTK_REAL* restrict const u, const CCTK_REAL p1o60dx, const ptrdiff_t cdj, const ptrdiff_t cdk)
{
  const ptrdiff_t cdi CCTK_ATTRIBUTE_UNUSED = sizeof(CCTK_REAL);
  return (-45*KRANC_GFOFFSET3D(u,-1,0,0) + 45*KRANC_GFOFFSET3D(u,1,0,0) + 9*KRANC_GFOFFSET3D(u,-2,0,0) - 9*KRANC_GFOFFSET3D(u,2,0,0) - KRANC_GFOFFSET3D(u,-3,0,0) + KRANC_GFOFFSET3D(u,3,0,0))*p1o60dx;
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDstandardfdOrder62(u) ((-45*KRANC_GFOFFSET3D(u,0,-1,0) + 45*KRANC_GFOFFSET3D(u,0,1,0) + 9*KRANC_GFOFFSET3D(u,0,-2,0) - 9*KRANC_GFOFFSET3D(u,0,2,0) - KRANC_GFOFFSET3D(u,0,-3,0) + KRANC_GFOFFSET3D(u,0,3,0))*p1o60dy)
#else
#  define PDstandardfdOrder62(u) (PDstandardfdOrder62_impl(u,p1o60dy,cdj,cdk))
static CCTK_REAL PDstandardfdOrder62_impl(const CCTK_REAL* restrict const u, const CCTK_REAL p1o60dy, const ptrdiff_t cdj, const ptrdiff_t cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL PDstandardfdOrder62_impl(const CCTK_REAL* restrict const u, const CCTK_REAL p1o60dy, const ptrdiff_t cdj, const ptrdiff_t cdk)
{
  const ptrdiff_t cdi CCTK_ATTRIBUTE_UNUSED = sizeof(CCTK_REAL);
  return (-45*KRANC_GFOFFSET3D(u,0,-1,0) + 45*KRANC_GFOFFSET3D(u,0,1,0) + 9*KRANC_GFOFFSET3D(u,0,-2,0) - 9*KRANC_GFOFFSET3D(u,0,2,0) - KRANC_GFOFFSET3D(u,0,-3,0) + KRANC_GFOFFSET3D(u,0,3,0))*p1o60dy;
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDstandardfdOrder63(u) ((-45*KRANC_GFOFFSET3D(u,0,0,-1) + 45*KRANC_GFOFFSET3D(u,0,0,1) + 9*KRANC_GFOFFSET3D(u,0,0,-2) - 9*KRANC_GFOFFSET3D(u,0,0,2) - KRANC_GFOFFSET3D(u,0,0,-3) + KRANC_GFOFFSET3D(u,0,0,3))*p1o60dz)
#else
#  define PDstandardfdOrder63(u) (PDstandardfdOrder63_impl(u,p1o60dz,cdj,cdk))
static CCTK_REAL PDstandardfdOrder63_impl(const CCTK_REAL* restrict const u, const CCTK_REAL p1o60dz, const ptrdiff_t cdj, const ptrdiff_t cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL PDstandardfdOrder63_impl(const CCTK_REAL* restrict const u, const CCTK_REAL p1o60dz, const ptrdiff_t cdj, const ptrdiff_t cdk)
{
  const ptrdiff_t cdi CCTK_ATTRIBUTE_UNUSED = sizeof(CCTK_REAL);
  return PDstandardfdOrder62_impl(u, p1o60dz, cdk, cdj);
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDstandardfdOrder81(u) ((-672*KRANC_GFOFFSET3D(u,-1,0,0) + 672*KRANC_GFOFFSET3D(u,1,0,0) + 168*KRANC_GFOFFSET3D(u,-2,0,0) - 168*KRANC_GFOFFSET3D(u,2,0,0) - 32*KRANC_GFOFFSET3D(u,-3,0,0) + 32*KRANC_GFOFFSET3D(u,3,0,0) + 3*KRANC_GFOFFSET3D(u,-4,0,0) - 3*KRANC_GFOFFSET3D(u,4,0,0))*p1o840dx)
#else
#  define PDstandardfdOrder81(u) (PDstandardfdOrder81_impl(u,p1o840dx,cdj,cdk))
static CCTK_REAL PDstandardfdOrder81_impl(const CCTK_REAL* restrict const u, const CCTK_REAL p1o840dx, const ptrdiff_t cdj, const ptrdiff_t cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL PDstandardfdOrder81_impl(const CCTK_REAL* restrict const u, const CCTK_REAL p1o840dx, const ptrdiff_t cdj, const ptrdiff_t cdk)
{
  const ptrdiff_t cdi CCTK_ATTRIBUTE_UNUSED = sizeof(CCTK_REAL);
  return (-672*KRANC_GFOFFSET3D(u,-1,0,0) + 672*KRANC_GFOFFSET3D(u,1,0,0) + 168*KRANC_GFOFFSET3D(u,-2,0,0) - 168*KRANC_GFOFFSET3D(u,2,0,0) - 32*KRANC_GFOFFSET3D(u,-3,0,0) + 32*KRANC_GFOFFSET3D(u,3,0,0) + 3*KRANC_GFOFFSET3D(u,-4,0,0) - 3*KRANC_GFOFFSET3D(u,4,0,0))*p1o840dx;
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDstandardfdOrder82(u) ((-672*KRANC_GFOFFSET3D(u,0,-1,0) + 672*KRANC_GFOFFSET3D(u,0,1,0) + 168*KRANC_GFOFFSET3D(u,0,-2,0) - 168*KRANC_GFOFFSET3D(u,0,2,0) - 32*KRANC_GFOFFSET3D(u,0,-3,0) + 32*KRANC_GFOFFSET3D(u,0,3,0) + 3*KRANC_GFOFFSET3D(u,0,-4,0) - 3*KRANC_GFOFFSET3D(u,0,4,0))*p1o840dy)
#else
#  define PDstandardfdOrder82(u) (PDstandardfdOrder82_impl(u,p1o840dy,cdj,cdk))
static CCTK_REAL PDstandardfdOrder82_impl(const CCTK_REAL* restrict const u, const CCTK_REAL p1o840dy, const ptrdiff_t cdj, const ptrdiff_t cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL PDstandardfdOrder82_impl(const CCTK_REAL* restrict const u, const CCTK_REAL p1o840dy, const ptrdiff_t cdj, const ptrdiff_t cdk)
{
  const ptrdiff_t cdi CCTK_ATTRIBUTE_UNUSED = sizeof(CCTK_REAL);
  return (-672*KRANC_GFOFFSET3D(u,0,-1,0) + 672*KRANC_GFOFFSET3D(u,0,1,0) + 168*KRANC_GFOFFSET3D(u,0,-2,0) - 168*KRANC_GFOFFSET3D(u,0,2,0) - 32*KRANC_GFOFFSET3D(u,0,-3,0) + 32*KRANC_GFOFFSET3D(u,0,3,0) + 3*KRANC_GFOFFSET3D(u,0,-4,0) - 3*KRANC_GFOFFSET3D(u,0,4,0))*p1o840dy;
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDstandardfdOrder83(u) ((-672*KRANC_GFOFFSET3D(u,0,0,-1) + 672*KRANC_GFOFFSET3D(u,0,0,1) + 168*KRANC_GFOFFSET3D(u,0,0,-2) - 168*KRANC_GFOFFSET3D(u,0,0,2) - 32*KRANC_GFOFFSET3D(u,0,0,-3) + 32*KRANC_GFOFFSET3D(u,0,0,3) + 3*KRANC_GFOFFSET3D(u,0,0,-4) - 3*KRANC_GFOFFSET3D(u,0,0,4))*p1o840dz)
#else
#  define PDstandardfdOrder83(u) (PDstandardfdOrder83_impl(u,p1o840dz,cdj,cdk))
static CCTK_REAL PDstandardfdOrder83_impl(const CCTK_REAL* restrict const u, const CCTK_REAL p1o840dz, const ptrdiff_t cdj, const ptrdiff_t cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL PDstandardfdOrder83_impl(const CCTK_REAL* restrict const u, const CCTK_REAL p1o840dz, const ptrdiff_t cdj, const ptrdiff_t cdk)
{
  const ptrdiff_t cdi CCTK_ATTRIBUTE_UNUSED = sizeof(CCTK_REAL);
  return PDstandardfdOrder82_impl(u, p1o840dz, cdk, cdj);
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDstandardfdOrder211(u) ((-2*KRANC_GFOFFSET3D(u,0,0,0) + KRANC_GFOFFSET3D(u,-1,0,0) + KRANC_GFOFFSET3D(u,1,0,0))*p1odx2)
#else
#  define PDstandardfdOrder211(u) (PDstandardfdOrder211_impl(u,p1odx2,cdj,cdk))
static CCTK_REAL PDstandardfdOrder211_impl(const CCTK_REAL* restrict const u, const CCTK_REAL p1odx2, const ptrdiff_t cdj, const ptrdiff_t cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL PDstandardfdOrder211_impl(const CCTK_REAL* restrict const u, const CCTK_REAL p1odx2, const ptrdiff_t cdj, const ptrdiff_t cdk)
{
  const ptrdiff_t cdi CCTK_ATTRIBUTE_UNUSED = sizeof(CCTK_REAL);
  return (-2*KRANC_GFOFFSET3D(u,0,0,0) + KRANC_GFOFFSET3D(u,-1,0,0) + KRANC_GFOFFSET3D(u,1,0,0))*p1odx2;
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDstandardfdOrder222(u) ((-2*KRANC_GFOFFSET3D(u,0,0,0) + KRANC_GFOFFSET3D(u,0,-1,0) + KRANC_GFOFFSET3D(u,0,1,0))*p1ody2)
#else
#  define PDstandardfdOrder222(u) (PDstandardfdOrder222_impl(u,p1ody2,cdj,cdk))
static CCTK_REAL PDstandardfdOrder222_impl(const CCTK_REAL* restrict const u, const CCTK_REAL p1ody2, const ptrdiff_t cdj, const ptrdiff_t cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL PDstandardfdOrder222_impl(const CCTK_REAL* restrict const u, const CCTK_REAL p1ody2, const ptrdiff_t cdj, const ptrdiff_t cdk)
{
  const ptrdiff_t cdi CCTK_ATTRIBUTE_UNUSED = sizeof(CCTK_REAL);
  return (-2*KRANC_GFOFFSET3D(u,0,0,0) + KRANC_GFOFFSET3D(u,0,-1,0) + KRANC_GFOFFSET3D(u,0,1,0))*p1ody2;
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDstandardfdOrder233(u) ((-2*KRANC_GFOFFSET3D(u,0,0,0) + KRANC_GFOFFSET3D(u,0,0,-1) + KRANC_GFOFFSET3D(u,0,0,1))*p1odz2)
#else
#  define PDstandardfdOrder233(u) (PDstandardfdOrder233_impl(u,p1odz2,cdj,cdk))
static CCTK_REAL PDstandardfdOrder233_impl(const CCTK_REAL* restrict const u, const CCTK_REAL p1odz2, const ptrdiff_t cdj, const ptrdiff_t cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL PDstandardfdOrder233_impl(const CCTK_REAL* restrict const u, const CCTK_REAL p1odz2, const ptrdiff_t cdj, const ptrdiff_t cdk)
{
  const ptrdiff_t cdi CCTK_ATTRIBUTE_UNUSED = sizeof(CCTK_REAL);
  return PDstandardfdOrder222_impl(u, p1odz2, cdk, cdj);
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDstandardfdOrder411(u) ((30*KRANC_GFOFFSET3D(u,0,0,0) - 16*(KRANC_GFOFFSET3D(u,-1,0,0) + KRANC_GFOFFSET3D(u,1,0,0)) + KRANC_GFOFFSET3D(u,-2,0,0) + KRANC_GFOFFSET3D(u,2,0,0))*pm1o12dx2)
#else
#  define PDstandardfdOrder411(u) (PDstandardfdOrder411_impl(u,pm1o12dx2,cdj,cdk))
static CCTK_REAL PDstandardfdOrder411_impl(const CCTK_REAL* restrict const u, const CCTK_REAL pm1o12dx2, const ptrdiff_t cdj, const ptrdiff_t cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL PDstandardfdOrder411_impl(const CCTK_REAL* restrict const u, const CCTK_REAL pm1o12dx2, const ptrdiff_t cdj, const ptrdiff_t cdk)
{
  const ptrdiff_t cdi CCTK_ATTRIBUTE_UNUSED = sizeof(CCTK_REAL);
  return (30*KRANC_GFOFFSET3D(u,0,0,0) - 16*(KRANC_GFOFFSET3D(u,-1,0,0) + KRANC_GFOFFSET3D(u,1,0,0)) + KRANC_GFOFFSET3D(u,-2,0,0) + KRANC_GFOFFSET3D(u,2,0,0))*pm1o12dx2;
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDstandardfdOrder422(u) ((30*KRANC_GFOFFSET3D(u,0,0,0) - 16*(KRANC_GFOFFSET3D(u,0,-1,0) + KRANC_GFOFFSET3D(u,0,1,0)) + KRANC_GFOFFSET3D(u,0,-2,0) + KRANC_GFOFFSET3D(u,0,2,0))*pm1o12dy2)
#else
#  define PDstandardfdOrder422(u) (PDstandardfdOrder422_impl(u,pm1o12dy2,cdj,cdk))
static CCTK_REAL PDstandardfdOrder422_impl(const CCTK_REAL* restrict const u, const CCTK_REAL pm1o12dy2, const ptrdiff_t cdj, const ptrdiff_t cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL PDstandardfdOrder422_impl(const CCTK_REAL* restrict const u, const CCTK_REAL pm1o12dy2, const ptrdiff_t cdj, const ptrdiff_t cdk)
{
  const ptrdiff_t cdi CCTK_ATTRIBUTE_UNUSED = sizeof(CCTK_REAL);
  return (30*KRANC_GFOFFSET3D(u,0,0,0) - 16*(KRANC_GFOFFSET3D(u,0,-1,0) + KRANC_GFOFFSET3D(u,0,1,0)) + KRANC_GFOFFSET3D(u,0,-2,0) + KRANC_GFOFFSET3D(u,0,2,0))*pm1o12dy2;
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDstandardfdOrder433(u) ((30*KRANC_GFOFFSET3D(u,0,0,0) - 16*(KRANC_GFOFFSET3D(u,0,0,-1) + KRANC_GFOFFSET3D(u,0,0,1)) + KRANC_GFOFFSET3D(u,0,0,-2) + KRANC_GFOFFSET3D(u,0,0,2))*pm1o12dz2)
#else
#  define PDstandardfdOrder433(u) (PDstandardfdOrder433_impl(u,pm1o12dz2,cdj,cdk))
static CCTK_REAL PDstandardfdOrder433_impl(const CCTK_REAL* restrict const u, const CCTK_REAL pm1o12dz2, const ptrdiff_t cdj, const ptrdiff_t cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL PDstandardfdOrder433_impl(const CCTK_REAL* restrict const u, const CCTK_REAL pm1o12dz2, const ptrdiff_t cdj, const ptrdiff_t cdk)
{
  const ptrdiff_t cdi CCTK_ATTRIBUTE_UNUSED = sizeof(CCTK_REAL);
  return PDstandardfdOrder422_impl(u, pm1o12dz2, cdk, cdj);
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDstandardfdOrder611(u) ((-490*KRANC_GFOFFSET3D(u,0,0,0) + 270*(KRANC_GFOFFSET3D(u,-1,0,0) + KRANC_GFOFFSET3D(u,1,0,0)) - 27*(KRANC_GFOFFSET3D(u,-2,0,0) + KRANC_GFOFFSET3D(u,2,0,0)) + 2*(KRANC_GFOFFSET3D(u,-3,0,0) + KRANC_GFOFFSET3D(u,3,0,0)))*p1o180dx2)
#else
#  define PDstandardfdOrder611(u) (PDstandardfdOrder611_impl(u,p1o180dx2,cdj,cdk))
static CCTK_REAL PDstandardfdOrder611_impl(const CCTK_REAL* restrict const u, const CCTK_REAL p1o180dx2, const ptrdiff_t cdj, const ptrdiff_t cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL PDstandardfdOrder611_impl(const CCTK_REAL* restrict const u, const CCTK_REAL p1o180dx2, const ptrdiff_t cdj, const ptrdiff_t cdk)
{
  const ptrdiff_t cdi CCTK_ATTRIBUTE_UNUSED = sizeof(CCTK_REAL);
  return (-490*KRANC_GFOFFSET3D(u,0,0,0) + 270*(KRANC_GFOFFSET3D(u,-1,0,0) + KRANC_GFOFFSET3D(u,1,0,0)) - 27*(KRANC_GFOFFSET3D(u,-2,0,0) + KRANC_GFOFFSET3D(u,2,0,0)) + 2*(KRANC_GFOFFSET3D(u,-3,0,0) + KRANC_GFOFFSET3D(u,3,0,0)))*p1o180dx2;
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDstandardfdOrder622(u) ((-490*KRANC_GFOFFSET3D(u,0,0,0) + 270*(KRANC_GFOFFSET3D(u,0,-1,0) + KRANC_GFOFFSET3D(u,0,1,0)) - 27*(KRANC_GFOFFSET3D(u,0,-2,0) + KRANC_GFOFFSET3D(u,0,2,0)) + 2*(KRANC_GFOFFSET3D(u,0,-3,0) + KRANC_GFOFFSET3D(u,0,3,0)))*p1o180dy2)
#else
#  define PDstandardfdOrder622(u) (PDstandardfdOrder622_impl(u,p1o180dy2,cdj,cdk))
static CCTK_REAL PDstandardfdOrder622_impl(const CCTK_REAL* restrict const u, const CCTK_REAL p1o180dy2, const ptrdiff_t cdj, const ptrdiff_t cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL PDstandardfdOrder622_impl(const CCTK_REAL* restrict const u, const CCTK_REAL p1o180dy2, const ptrdiff_t cdj, const ptrdiff_t cdk)
{
  const ptrdiff_t cdi CCTK_ATTRIBUTE_UNUSED = sizeof(CCTK_REAL);
  return (-490*KRANC_GFOFFSET3D(u,0,0,0) + 270*(KRANC_GFOFFSET3D(u,0,-1,0) + KRANC_GFOFFSET3D(u,0,1,0)) - 27*(KRANC_GFOFFSET3D(u,0,-2,0) + KRANC_GFOFFSET3D(u,0,2,0)) + 2*(KRANC_GFOFFSET3D(u,0,-3,0) + KRANC_GFOFFSET3D(u,0,3,0)))*p1o180dy2;
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDstandardfdOrder633(u) ((-490*KRANC_GFOFFSET3D(u,0,0,0) + 270*(KRANC_GFOFFSET3D(u,0,0,-1) + KRANC_GFOFFSET3D(u,0,0,1)) - 27*(KRANC_GFOFFSET3D(u,0,0,-2) + KRANC_GFOFFSET3D(u,0,0,2)) + 2*(KRANC_GFOFFSET3D(u,0,0,-3) + KRANC_GFOFFSET3D(u,0,0,3)))*p1o180dz2)
#else
#  define PDstandardfdOrder633(u) (PDstandardfdOrder633_impl(u,p1o180dz2,cdj,cdk))
static CCTK_REAL PDstandardfdOrder633_impl(const CCTK_REAL* restrict const u, const CCTK_REAL p1o180dz2, const ptrdiff_t cdj, const ptrdiff_t cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL PDstandardfdOrder633_impl(const CCTK_REAL* restrict const u, const CCTK_REAL p1o180dz2, const ptrdiff_t cdj, const ptrdiff_t cdk)
{
  const ptrdiff_t cdi CCTK_ATTRIBUTE_UNUSED = sizeof(CCTK_REAL);
  return PDstandardfdOrder622_impl(u, p1o180dz2, cdk, cdj);
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDstandardfdOrder811(u) ((-14350*KRANC_GFOFFSET3D(u,0,0,0) + 8064*(KRANC_GFOFFSET3D(u,-1,0,0) + KRANC_GFOFFSET3D(u,1,0,0)) - 1008*(KRANC_GFOFFSET3D(u,-2,0,0) + KRANC_GFOFFSET3D(u,2,0,0)) + 128*(KRANC_GFOFFSET3D(u,-3,0,0) + KRANC_GFOFFSET3D(u,3,0,0)) - 9*(KRANC_GFOFFSET3D(u,-4,0,0) + KRANC_GFOFFSET3D(u,4,0,0)))*p1o5040dx2)
#else
#  define PDstandardfdOrder811(u) (PDstandardfdOrder811_impl(u,p1o5040dx2,cdj,cdk))
static CCTK_REAL PDstandardfdOrder811_impl(const CCTK_REAL* restrict const u, const CCTK_REAL p1o5040dx2, const ptrdiff_t cdj, const ptrdiff_t cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL PDstandardfdOrder811_impl(const CCTK_REAL* restrict const u, const CCTK_REAL p1o5040dx2, const ptrdiff_t cdj, const ptrdiff_t cdk)
{
  const ptrdiff_t cdi CCTK_ATTRIBUTE_UNUSED = sizeof(CCTK_REAL);
  return (-14350*KRANC_GFOFFSET3D(u,0,0,0) + 8064*(KRANC_GFOFFSET3D(u,-1,0,0) + KRANC_GFOFFSET3D(u,1,0,0)) - 1008*(KRANC_GFOFFSET3D(u,-2,0,0) + KRANC_GFOFFSET3D(u,2,0,0)) + 128*(KRANC_GFOFFSET3D(u,-3,0,0) + KRANC_GFOFFSET3D(u,3,0,0)) - 9*(KRANC_GFOFFSET3D(u,-4,0,0) + KRANC_GFOFFSET3D(u,4,0,0)))*p1o5040dx2;
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDstandardfdOrder822(u) ((-14350*KRANC_GFOFFSET3D(u,0,0,0) + 8064*(KRANC_GFOFFSET3D(u,0,-1,0) + KRANC_GFOFFSET3D(u,0,1,0)) - 1008*(KRANC_GFOFFSET3D(u,0,-2,0) + KRANC_GFOFFSET3D(u,0,2,0)) + 128*(KRANC_GFOFFSET3D(u,0,-3,0) + KRANC_GFOFFSET3D(u,0,3,0)) - 9*(KRANC_GFOFFSET3D(u,0,-4,0) + KRANC_GFOFFSET3D(u,0,4,0)))*p1o5040dy2)
#else
#  define PDstandardfdOrder822(u) (PDstandardfdOrder822_impl(u,p1o5040dy2,cdj,cdk))
static CCTK_REAL PDstandardfdOrder822_impl(const CCTK_REAL* restrict const u, const CCTK_REAL p1o5040dy2, const ptrdiff_t cdj, const ptrdiff_t cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL PDstandardfdOrder822_impl(const CCTK_REAL* restrict const u, const CCTK_REAL p1o5040dy2, const ptrdiff_t cdj, const ptrdiff_t cdk)
{
  const ptrdiff_t cdi CCTK_ATTRIBUTE_UNUSED = sizeof(CCTK_REAL);
  return (-14350*KRANC_GFOFFSET3D(u,0,0,0) + 8064*(KRANC_GFOFFSET3D(u,0,-1,0) + KRANC_GFOFFSET3D(u,0,1,0)) - 1008*(KRANC_GFOFFSET3D(u,0,-2,0) + KRANC_GFOFFSET3D(u,0,2,0)) + 128*(KRANC_GFOFFSET3D(u,0,-3,0) + KRANC_GFOFFSET3D(u,0,3,0)) - 9*(KRANC_GFOFFSET3D(u,0,-4,0) + KRANC_GFOFFSET3D(u,0,4,0)))*p1o5040dy2;
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDstandardfdOrder833(u) ((-14350*KRANC_GFOFFSET3D(u,0,0,0) + 8064*(KRANC_GFOFFSET3D(u,0,0,-1) + KRANC_GFOFFSET3D(u,0,0,1)) - 1008*(KRANC_GFOFFSET3D(u,0,0,-2) + KRANC_GFOFFSET3D(u,0,0,2)) + 128*(KRANC_GFOFFSET3D(u,0,0,-3) + KRANC_GFOFFSET3D(u,0,0,3)) - 9*(KRANC_GFOFFSET3D(u,0,0,-4) + KRANC_GFOFFSET3D(u,0,0,4)))*p1o5040dz2)
#else
#  define PDstandardfdOrder833(u) (PDstandardfdOrder833_impl(u,p1o5040dz2,cdj,cdk))
static CCTK_REAL PDstandardfdOrder833_impl(const CCTK_REAL* restrict const u, const CCTK_REAL p1o5040dz2, const ptrdiff_t cdj, const ptrdiff_t cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL PDstandardfdOrder833_impl(const CCTK_REAL* restrict const u, const CCTK_REAL p1o5040dz2, const ptrdiff_t cdj, const ptrdiff_t cdk)
{
  const ptrdiff_t cdi CCTK_ATTRIBUTE_UNUSED = sizeof(CCTK_REAL);
  return PDstandardfdOrder822_impl(u, p1o5040dz2, cdk, cdj);
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDstandardfdOrder212(u) ((KRANC_GFOFFSET3D(u,-1,-1,0) - KRANC_GFOFFSET3D(u,-1,1,0) - KRANC_GFOFFSET3D(u,1,-1,0) + KRANC_GFOFFSET3D(u,1,1,0))*p1o4dxdy)
#else
#  define PDstandardfdOrder212(u) (PDstandardfdOrder212_impl(u,p1o4dxdy,cdj,cdk))
static CCTK_REAL PDstandardfdOrder212_impl(const CCTK_REAL* restrict const u, const CCTK_REAL p1o4dxdy, const ptrdiff_t cdj, const ptrdiff_t cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL PDstandardfdOrder212_impl(const CCTK_REAL* restrict const u, const CCTK_REAL p1o4dxdy, const ptrdiff_t cdj, const ptrdiff_t cdk)
{
  const ptrdiff_t cdi CCTK_ATTRIBUTE_UNUSED = sizeof(CCTK_REAL);
  return (KRANC_GFOFFSET3D(u,-1,-1,0) - KRANC_GFOFFSET3D(u,-1,1,0) - KRANC_GFOFFSET3D(u,1,-1,0) + KRANC_GFOFFSET3D(u,1,1,0))*p1o4dxdy;
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDstandardfdOrder213(u) ((KRANC_GFOFFSET3D(u,-1,0,-1) - KRANC_GFOFFSET3D(u,-1,0,1) - KRANC_GFOFFSET3D(u,1,0,-1) + KRANC_GFOFFSET3D(u,1,0,1))*p1o4dxdz)
#else
#  define PDstandardfdOrder213(u) (PDstandardfdOrder213_impl(u,p1o4dxdz,cdj,cdk))
static CCTK_REAL PDstandardfdOrder213_impl(const CCTK_REAL* restrict const u, const CCTK_REAL p1o4dxdz, const ptrdiff_t cdj, const ptrdiff_t cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL PDstandardfdOrder213_impl(const CCTK_REAL* restrict const u, const CCTK_REAL p1o4dxdz, const ptrdiff_t cdj, const ptrdiff_t cdk)
{
  const ptrdiff_t cdi CCTK_ATTRIBUTE_UNUSED = sizeof(CCTK_REAL);
  return PDstandardfdOrder212_impl(u, p1o4dxdz, cdk, cdj);
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDstandardfdOrder221(u) ((KRANC_GFOFFSET3D(u,-1,-1,0) - KRANC_GFOFFSET3D(u,-1,1,0) - KRANC_GFOFFSET3D(u,1,-1,0) + KRANC_GFOFFSET3D(u,1,1,0))*p1o4dxdy)
#else
#  define PDstandardfdOrder221(u) (PDstandardfdOrder221_impl(u,p1o4dxdy,cdj,cdk))
static CCTK_REAL PDstandardfdOrder221_impl(const CCTK_REAL* restrict const u, const CCTK_REAL p1o4dxdy, const ptrdiff_t cdj, const ptrdiff_t cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL PDstandardfdOrder221_impl(const CCTK_REAL* restrict const u, const CCTK_REAL p1o4dxdy, const ptrdiff_t cdj, const ptrdiff_t cdk)
{
  const ptrdiff_t cdi CCTK_ATTRIBUTE_UNUSED = sizeof(CCTK_REAL);
  return PDstandardfdOrder212_impl(u, p1o4dxdy, cdj, cdk);
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDstandardfdOrder223(u) ((KRANC_GFOFFSET3D(u,0,-1,-1) - KRANC_GFOFFSET3D(u,0,-1,1) - KRANC_GFOFFSET3D(u,0,1,-1) + KRANC_GFOFFSET3D(u,0,1,1))*p1o4dydz)
#else
#  define PDstandardfdOrder223(u) (PDstandardfdOrder223_impl(u,p1o4dydz,cdj,cdk))
static CCTK_REAL PDstandardfdOrder223_impl(const CCTK_REAL* restrict const u, const CCTK_REAL p1o4dydz, const ptrdiff_t cdj, const ptrdiff_t cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL PDstandardfdOrder223_impl(const CCTK_REAL* restrict const u, const CCTK_REAL p1o4dydz, const ptrdiff_t cdj, const ptrdiff_t cdk)
{
  const ptrdiff_t cdi CCTK_ATTRIBUTE_UNUSED = sizeof(CCTK_REAL);
  return (KRANC_GFOFFSET3D(u,0,-1,-1) - KRANC_GFOFFSET3D(u,0,-1,1) - KRANC_GFOFFSET3D(u,0,1,-1) + KRANC_GFOFFSET3D(u,0,1,1))*p1o4dydz;
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDstandardfdOrder231(u) ((KRANC_GFOFFSET3D(u,-1,0,-1) - KRANC_GFOFFSET3D(u,-1,0,1) - KRANC_GFOFFSET3D(u,1,0,-1) + KRANC_GFOFFSET3D(u,1,0,1))*p1o4dxdz)
#else
#  define PDstandardfdOrder231(u) (PDstandardfdOrder231_impl(u,p1o4dxdz,cdj,cdk))
static CCTK_REAL PDstandardfdOrder231_impl(const CCTK_REAL* restrict const u, const CCTK_REAL p1o4dxdz, const ptrdiff_t cdj, const ptrdiff_t cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL PDstandardfdOrder231_impl(const CCTK_REAL* restrict const u, const CCTK_REAL p1o4dxdz, const ptrdiff_t cdj, const ptrdiff_t cdk)
{
  const ptrdiff_t cdi CCTK_ATTRIBUTE_UNUSED = sizeof(CCTK_REAL);
  return PDstandardfdOrder212_impl(u, p1o4dxdz, cdk, cdj);
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDstandardfdOrder232(u) ((KRANC_GFOFFSET3D(u,0,-1,-1) - KRANC_GFOFFSET3D(u,0,-1,1) - KRANC_GFOFFSET3D(u,0,1,-1) + KRANC_GFOFFSET3D(u,0,1,1))*p1o4dydz)
#else
#  define PDstandardfdOrder232(u) (PDstandardfdOrder232_impl(u,p1o4dydz,cdj,cdk))
static CCTK_REAL PDstandardfdOrder232_impl(const CCTK_REAL* restrict const u, const CCTK_REAL p1o4dydz, const ptrdiff_t cdj, const ptrdiff_t cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL PDstandardfdOrder232_impl(const CCTK_REAL* restrict const u, const CCTK_REAL p1o4dydz, const ptrdiff_t cdj, const ptrdiff_t cdk)
{
  const ptrdiff_t cdi CCTK_ATTRIBUTE_UNUSED = sizeof(CCTK_REAL);
  return PDstandardfdOrder223_impl(u, p1o4dydz, cdj, cdk);
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDstandardfdOrder412(u) ((-64*(KRANC_GFOFFSET3D(u,-1,1,0) + KRANC_GFOFFSET3D(u,1,-1,0)) + 64*(KRANC_GFOFFSET3D(u,-1,-1,0) + KRANC_GFOFFSET3D(u,1,1,0)) + 8*(KRANC_GFOFFSET3D(u,-1,2,0) + KRANC_GFOFFSET3D(u,1,-2,0) + KRANC_GFOFFSET3D(u,-2,1,0) + KRANC_GFOFFSET3D(u,2,-1,0)) - 8*(KRANC_GFOFFSET3D(u,-1,-2,0) + KRANC_GFOFFSET3D(u,1,2,0) + KRANC_GFOFFSET3D(u,-2,-1,0) + KRANC_GFOFFSET3D(u,2,1,0)) + KRANC_GFOFFSET3D(u,-2,-2,0) - KRANC_GFOFFSET3D(u,-2,2,0) - KRANC_GFOFFSET3D(u,2,-2,0) + KRANC_GFOFFSET3D(u,2,2,0))*p1o144dxdy)
#else
#  define PDstandardfdOrder412(u) (PDstandardfdOrder412_impl(u,p1o144dxdy,cdj,cdk))
static CCTK_REAL PDstandardfdOrder412_impl(const CCTK_REAL* restrict const u, const CCTK_REAL p1o144dxdy, const ptrdiff_t cdj, const ptrdiff_t cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL PDstandardfdOrder412_impl(const CCTK_REAL* restrict const u, const CCTK_REAL p1o144dxdy, const ptrdiff_t cdj, const ptrdiff_t cdk)
{
  const ptrdiff_t cdi CCTK_ATTRIBUTE_UNUSED = sizeof(CCTK_REAL);
  return (-64*(KRANC_GFOFFSET3D(u,-1,1,0) + KRANC_GFOFFSET3D(u,1,-1,0)) + 64*(KRANC_GFOFFSET3D(u,-1,-1,0) + KRANC_GFOFFSET3D(u,1,1,0)) + 8*(KRANC_GFOFFSET3D(u,-1,2,0) + KRANC_GFOFFSET3D(u,1,-2,0) + KRANC_GFOFFSET3D(u,-2,1,0) + KRANC_GFOFFSET3D(u,2,-1,0)) - 8*(KRANC_GFOFFSET3D(u,-1,-2,0) + KRANC_GFOFFSET3D(u,1,2,0) + KRANC_GFOFFSET3D(u,-2,-1,0) + KRANC_GFOFFSET3D(u,2,1,0)) + KRANC_GFOFFSET3D(u,-2,-2,0) - KRANC_GFOFFSET3D(u,-2,2,0) - KRANC_GFOFFSET3D(u,2,-2,0) + KRANC_GFOFFSET3D(u,2,2,0))*p1o144dxdy;
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDstandardfdOrder413(u) ((-64*(KRANC_GFOFFSET3D(u,-1,0,1) + KRANC_GFOFFSET3D(u,1,0,-1)) + 64*(KRANC_GFOFFSET3D(u,-1,0,-1) + KRANC_GFOFFSET3D(u,1,0,1)) + 8*(KRANC_GFOFFSET3D(u,-1,0,2) + KRANC_GFOFFSET3D(u,1,0,-2) + KRANC_GFOFFSET3D(u,-2,0,1) + KRANC_GFOFFSET3D(u,2,0,-1)) - 8*(KRANC_GFOFFSET3D(u,-1,0,-2) + KRANC_GFOFFSET3D(u,1,0,2) + KRANC_GFOFFSET3D(u,-2,0,-1) + KRANC_GFOFFSET3D(u,2,0,1)) + KRANC_GFOFFSET3D(u,-2,0,-2) - KRANC_GFOFFSET3D(u,-2,0,2) - KRANC_GFOFFSET3D(u,2,0,-2) + KRANC_GFOFFSET3D(u,2,0,2))*p1o144dxdz)
#else
#  define PDstandardfdOrder413(u) (PDstandardfdOrder413_impl(u,p1o144dxdz,cdj,cdk))
static CCTK_REAL PDstandardfdOrder413_impl(const CCTK_REAL* restrict const u, const CCTK_REAL p1o144dxdz, const ptrdiff_t cdj, const ptrdiff_t cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL PDstandardfdOrder413_impl(const CCTK_REAL* restrict const u, const CCTK_REAL p1o144dxdz, const ptrdiff_t cdj, const ptrdiff_t cdk)
{
  const ptrdiff_t cdi CCTK_ATTRIBUTE_UNUSED = sizeof(CCTK_REAL);
  return PDstandardfdOrder412_impl(u, p1o144dxdz, cdk, cdj);
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDstandardfdOrder421(u) ((-64*(KRANC_GFOFFSET3D(u,-1,1,0) + KRANC_GFOFFSET3D(u,1,-1,0)) + 64*(KRANC_GFOFFSET3D(u,-1,-1,0) + KRANC_GFOFFSET3D(u,1,1,0)) + 8*(KRANC_GFOFFSET3D(u,-1,2,0) + KRANC_GFOFFSET3D(u,1,-2,0) + KRANC_GFOFFSET3D(u,-2,1,0) + KRANC_GFOFFSET3D(u,2,-1,0)) - 8*(KRANC_GFOFFSET3D(u,-1,-2,0) + KRANC_GFOFFSET3D(u,1,2,0) + KRANC_GFOFFSET3D(u,-2,-1,0) + KRANC_GFOFFSET3D(u,2,1,0)) + KRANC_GFOFFSET3D(u,-2,-2,0) - KRANC_GFOFFSET3D(u,-2,2,0) - KRANC_GFOFFSET3D(u,2,-2,0) + KRANC_GFOFFSET3D(u,2,2,0))*p1o144dxdy)
#else
#  define PDstandardfdOrder421(u) (PDstandardfdOrder421_impl(u,p1o144dxdy,cdj,cdk))
static CCTK_REAL PDstandardfdOrder421_impl(const CCTK_REAL* restrict const u, const CCTK_REAL p1o144dxdy, const ptrdiff_t cdj, const ptrdiff_t cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL PDstandardfdOrder421_impl(const CCTK_REAL* restrict const u, const CCTK_REAL p1o144dxdy, const ptrdiff_t cdj, const ptrdiff_t cdk)
{
  const ptrdiff_t cdi CCTK_ATTRIBUTE_UNUSED = sizeof(CCTK_REAL);
  return PDstandardfdOrder412_impl(u, p1o144dxdy, cdj, cdk);
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDstandardfdOrder423(u) ((-64*(KRANC_GFOFFSET3D(u,0,-1,1) + KRANC_GFOFFSET3D(u,0,1,-1)) + 64*(KRANC_GFOFFSET3D(u,0,-1,-1) + KRANC_GFOFFSET3D(u,0,1,1)) + 8*(KRANC_GFOFFSET3D(u,0,-1,2) + KRANC_GFOFFSET3D(u,0,1,-2) + KRANC_GFOFFSET3D(u,0,-2,1) + KRANC_GFOFFSET3D(u,0,2,-1)) - 8*(KRANC_GFOFFSET3D(u,0,-1,-2) + KRANC_GFOFFSET3D(u,0,1,2) + KRANC_GFOFFSET3D(u,0,-2,-1) + KRANC_GFOFFSET3D(u,0,2,1)) + KRANC_GFOFFSET3D(u,0,-2,-2) - KRANC_GFOFFSET3D(u,0,-2,2) - KRANC_GFOFFSET3D(u,0,2,-2) + KRANC_GFOFFSET3D(u,0,2,2))*p1o144dydz)
#else
#  define PDstandardfdOrder423(u) (PDstandardfdOrder423_impl(u,p1o144dydz,cdj,cdk))
static CCTK_REAL PDstandardfdOrder423_impl(const CCTK_REAL* restrict const u, const CCTK_REAL p1o144dydz, const ptrdiff_t cdj, const ptrdiff_t cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL PDstandardfdOrder423_impl(const CCTK_REAL* restrict const u, const CCTK_REAL p1o144dydz, const ptrdiff_t cdj, const ptrdiff_t cdk)
{
  const ptrdiff_t cdi CCTK_ATTRIBUTE_UNUSED = sizeof(CCTK_REAL);
  return (-64*(KRANC_GFOFFSET3D(u,0,-1,1) + KRANC_GFOFFSET3D(u,0,1,-1)) + 64*(KRANC_GFOFFSET3D(u,0,-1,-1) + KRANC_GFOFFSET3D(u,0,1,1)) + 8*(KRANC_GFOFFSET3D(u,0,-1,2) + KRANC_GFOFFSET3D(u,0,1,-2) + KRANC_GFOFFSET3D(u,0,-2,1) + KRANC_GFOFFSET3D(u,0,2,-1)) - 8*(KRANC_GFOFFSET3D(u,0,-1,-2) + KRANC_GFOFFSET3D(u,0,1,2) + KRANC_GFOFFSET3D(u,0,-2,-1) + KRANC_GFOFFSET3D(u,0,2,1)) + KRANC_GFOFFSET3D(u,0,-2,-2) - KRANC_GFOFFSET3D(u,0,-2,2) - KRANC_GFOFFSET3D(u,0,2,-2) + KRANC_GFOFFSET3D(u,0,2,2))*p1o144dydz;
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDstandardfdOrder431(u) ((-64*(KRANC_GFOFFSET3D(u,-1,0,1) + KRANC_GFOFFSET3D(u,1,0,-1)) + 64*(KRANC_GFOFFSET3D(u,-1,0,-1) + KRANC_GFOFFSET3D(u,1,0,1)) + 8*(KRANC_GFOFFSET3D(u,-1,0,2) + KRANC_GFOFFSET3D(u,1,0,-2) + KRANC_GFOFFSET3D(u,-2,0,1) + KRANC_GFOFFSET3D(u,2,0,-1)) - 8*(KRANC_GFOFFSET3D(u,-1,0,-2) + KRANC_GFOFFSET3D(u,1,0,2) + KRANC_GFOFFSET3D(u,-2,0,-1) + KRANC_GFOFFSET3D(u,2,0,1)) + KRANC_GFOFFSET3D(u,-2,0,-2) - KRANC_GFOFFSET3D(u,-2,0,2) - KRANC_GFOFFSET3D(u,2,0,-2) + KRANC_GFOFFSET3D(u,2,0,2))*p1o144dxdz)
#else
#  define PDstandardfdOrder431(u) (PDstandardfdOrder431_impl(u,p1o144dxdz,cdj,cdk))
static CCTK_REAL PDstandardfdOrder431_impl(const CCTK_REAL* restrict const u, const CCTK_REAL p1o144dxdz, const ptrdiff_t cdj, const ptrdiff_t cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL PDstandardfdOrder431_impl(const CCTK_REAL* restrict const u, const CCTK_REAL p1o144dxdz, const ptrdiff_t cdj, const ptrdiff_t cdk)
{
  const ptrdiff_t cdi CCTK_ATTRIBUTE_UNUSED = sizeof(CCTK_REAL);
  return PDstandardfdOrder412_impl(u, p1o144dxdz, cdk, cdj);
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDstandardfdOrder432(u) ((-64*(KRANC_GFOFFSET3D(u,0,-1,1) + KRANC_GFOFFSET3D(u,0,1,-1)) + 64*(KRANC_GFOFFSET3D(u,0,-1,-1) + KRANC_GFOFFSET3D(u,0,1,1)) + 8*(KRANC_GFOFFSET3D(u,0,-1,2) + KRANC_GFOFFSET3D(u,0,1,-2) + KRANC_GFOFFSET3D(u,0,-2,1) + KRANC_GFOFFSET3D(u,0,2,-1)) - 8*(KRANC_GFOFFSET3D(u,0,-1,-2) + KRANC_GFOFFSET3D(u,0,1,2) + KRANC_GFOFFSET3D(u,0,-2,-1) + KRANC_GFOFFSET3D(u,0,2,1)) + KRANC_GFOFFSET3D(u,0,-2,-2) - KRANC_GFOFFSET3D(u,0,-2,2) - KRANC_GFOFFSET3D(u,0,2,-2) + KRANC_GFOFFSET3D(u,0,2,2))*p1o144dydz)
#else
#  define PDstandardfdOrder432(u) (PDstandardfdOrder432_impl(u,p1o144dydz,cdj,cdk))
static CCTK_REAL PDstandardfdOrder432_impl(const CCTK_REAL* restrict const u, const CCTK_REAL p1o144dydz, const ptrdiff_t cdj, const ptrdiff_t cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL PDstandardfdOrder432_impl(const CCTK_REAL* restrict const u, const CCTK_REAL p1o144dydz, const ptrdiff_t cdj, const ptrdiff_t cdk)
{
  const ptrdiff_t cdi CCTK_ATTRIBUTE_UNUSED = sizeof(CCTK_REAL);
  return PDstandardfdOrder423_impl(u, p1o144dydz, cdj, cdk);
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDstandardfdOrder612(u) ((-2025*(KRANC_GFOFFSET3D(u,-1,1,0) + KRANC_GFOFFSET3D(u,1,-1,0)) + 2025*(KRANC_GFOFFSET3D(u,-1,-1,0) + KRANC_GFOFFSET3D(u,1,1,0)) + 405*(KRANC_GFOFFSET3D(u,-1,2,0) + KRANC_GFOFFSET3D(u,1,-2,0) + KRANC_GFOFFSET3D(u,-2,1,0) + KRANC_GFOFFSET3D(u,2,-1,0)) - 405*(KRANC_GFOFFSET3D(u,-1,-2,0) + KRANC_GFOFFSET3D(u,1,2,0) + KRANC_GFOFFSET3D(u,-2,-1,0) + KRANC_GFOFFSET3D(u,2,1,0)) - 81*(KRANC_GFOFFSET3D(u,-2,2,0) + KRANC_GFOFFSET3D(u,2,-2,0)) + 81*(KRANC_GFOFFSET3D(u,-2,-2,0) + KRANC_GFOFFSET3D(u,2,2,0)) - 45*(KRANC_GFOFFSET3D(u,-1,3,0) + KRANC_GFOFFSET3D(u,1,-3,0) + KRANC_GFOFFSET3D(u,-3,1,0) + KRANC_GFOFFSET3D(u,3,-1,0)) + 45*(KRANC_GFOFFSET3D(u,-1,-3,0) + KRANC_GFOFFSET3D(u,1,3,0) + KRANC_GFOFFSET3D(u,-3,-1,0) + KRANC_GFOFFSET3D(u,3,1,0)) + 9*(KRANC_GFOFFSET3D(u,-2,3,0) + KRANC_GFOFFSET3D(u,2,-3,0) + KRANC_GFOFFSET3D(u,-3,2,0) + KRANC_GFOFFSET3D(u,3,-2,0)) - 9*(KRANC_GFOFFSET3D(u,-2,-3,0) + KRANC_GFOFFSET3D(u,2,3,0) + KRANC_GFOFFSET3D(u,-3,-2,0) + KRANC_GFOFFSET3D(u,3,2,0)) + KRANC_GFOFFSET3D(u,-3,-3,0) - KRANC_GFOFFSET3D(u,-3,3,0) - KRANC_GFOFFSET3D(u,3,-3,0) + KRANC_GFOFFSET3D(u,3,3,0))*p1o3600dxdy)
#else
#  define PDstandardfdOrder612(u) (PDstandardfdOrder612_impl(u,p1o3600dxdy,cdj,cdk))
static CCTK_REAL PDstandardfdOrder612_impl(const CCTK_REAL* restrict const u, const CCTK_REAL p1o3600dxdy, const ptrdiff_t cdj, const ptrdiff_t cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL PDstandardfdOrder612_impl(const CCTK_REAL* restrict const u, const CCTK_REAL p1o3600dxdy, const ptrdiff_t cdj, const ptrdiff_t cdk)
{
  const ptrdiff_t cdi CCTK_ATTRIBUTE_UNUSED = sizeof(CCTK_REAL);
  return (-2025*(KRANC_GFOFFSET3D(u,-1,1,0) + KRANC_GFOFFSET3D(u,1,-1,0)) + 2025*(KRANC_GFOFFSET3D(u,-1,-1,0) + KRANC_GFOFFSET3D(u,1,1,0)) + 405*(KRANC_GFOFFSET3D(u,-1,2,0) + KRANC_GFOFFSET3D(u,1,-2,0) + KRANC_GFOFFSET3D(u,-2,1,0) + KRANC_GFOFFSET3D(u,2,-1,0)) - 405*(KRANC_GFOFFSET3D(u,-1,-2,0) + KRANC_GFOFFSET3D(u,1,2,0) + KRANC_GFOFFSET3D(u,-2,-1,0) + KRANC_GFOFFSET3D(u,2,1,0)) - 81*(KRANC_GFOFFSET3D(u,-2,2,0) + KRANC_GFOFFSET3D(u,2,-2,0)) + 81*(KRANC_GFOFFSET3D(u,-2,-2,0) + KRANC_GFOFFSET3D(u,2,2,0)) - 45*(KRANC_GFOFFSET3D(u,-1,3,0) + KRANC_GFOFFSET3D(u,1,-3,0) + KRANC_GFOFFSET3D(u,-3,1,0) + KRANC_GFOFFSET3D(u,3,-1,0)) + 45*(KRANC_GFOFFSET3D(u,-1,-3,0) + KRANC_GFOFFSET3D(u,1,3,0) + KRANC_GFOFFSET3D(u,-3,-1,0) + KRANC_GFOFFSET3D(u,3,1,0)) + 9*(KRANC_GFOFFSET3D(u,-2,3,0) + KRANC_GFOFFSET3D(u,2,-3,0) + KRANC_GFOFFSET3D(u,-3,2,0) + KRANC_GFOFFSET3D(u,3,-2,0)) - 9*(KRANC_GFOFFSET3D(u,-2,-3,0) + KRANC_GFOFFSET3D(u,2,3,0) + KRANC_GFOFFSET3D(u,-3,-2,0) + KRANC_GFOFFSET3D(u,3,2,0)) + KRANC_GFOFFSET3D(u,-3,-3,0) - KRANC_GFOFFSET3D(u,-3,3,0) - KRANC_GFOFFSET3D(u,3,-3,0) + KRANC_GFOFFSET3D(u,3,3,0))*p1o3600dxdy;
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDstandardfdOrder613(u) ((-2025*(KRANC_GFOFFSET3D(u,-1,0,1) + KRANC_GFOFFSET3D(u,1,0,-1)) + 2025*(KRANC_GFOFFSET3D(u,-1,0,-1) + KRANC_GFOFFSET3D(u,1,0,1)) + 405*(KRANC_GFOFFSET3D(u,-1,0,2) + KRANC_GFOFFSET3D(u,1,0,-2) + KRANC_GFOFFSET3D(u,-2,0,1) + KRANC_GFOFFSET3D(u,2,0,-1)) - 405*(KRANC_GFOFFSET3D(u,-1,0,-2) + KRANC_GFOFFSET3D(u,1,0,2) + KRANC_GFOFFSET3D(u,-2,0,-1) + KRANC_GFOFFSET3D(u,2,0,1)) - 81*(KRANC_GFOFFSET3D(u,-2,0,2) + KRANC_GFOFFSET3D(u,2,0,-2)) + 81*(KRANC_GFOFFSET3D(u,-2,0,-2) + KRANC_GFOFFSET3D(u,2,0,2)) - 45*(KRANC_GFOFFSET3D(u,-1,0,3) + KRANC_GFOFFSET3D(u,1,0,-3) + KRANC_GFOFFSET3D(u,-3,0,1) + KRANC_GFOFFSET3D(u,3,0,-1)) + 45*(KRANC_GFOFFSET3D(u,-1,0,-3) + KRANC_GFOFFSET3D(u,1,0,3) + KRANC_GFOFFSET3D(u,-3,0,-1) + KRANC_GFOFFSET3D(u,3,0,1)) + 9*(KRANC_GFOFFSET3D(u,-2,0,3) + KRANC_GFOFFSET3D(u,2,0,-3) + KRANC_GFOFFSET3D(u,-3,0,2) + KRANC_GFOFFSET3D(u,3,0,-2)) - 9*(KRANC_GFOFFSET3D(u,-2,0,-3) + KRANC_GFOFFSET3D(u,2,0,3) + KRANC_GFOFFSET3D(u,-3,0,-2) + KRANC_GFOFFSET3D(u,3,0,2)) + KRANC_GFOFFSET3D(u,-3,0,-3) - KRANC_GFOFFSET3D(u,-3,0,3) - KRANC_GFOFFSET3D(u,3,0,-3) + KRANC_GFOFFSET3D(u,3,0,3))*p1o3600dxdz)
#else
#  define PDstandardfdOrder613(u) (PDstandardfdOrder613_impl(u,p1o3600dxdz,cdj,cdk))
static CCTK_REAL PDstandardfdOrder613_impl(const CCTK_REAL* restrict const u, const CCTK_REAL p1o3600dxdz, const ptrdiff_t cdj, const ptrdiff_t cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL PDstandardfdOrder613_impl(const CCTK_REAL* restrict const u, const CCTK_REAL p1o3600dxdz, const ptrdiff_t cdj, const ptrdiff_t cdk)
{
  const ptrdiff_t cdi CCTK_ATTRIBUTE_UNUSED = sizeof(CCTK_REAL);
  return PDstandardfdOrder612_impl(u, p1o3600dxdz, cdk, cdj);
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDstandardfdOrder621(u) ((-2025*(KRANC_GFOFFSET3D(u,-1,1,0) + KRANC_GFOFFSET3D(u,1,-1,0)) + 2025*(KRANC_GFOFFSET3D(u,-1,-1,0) + KRANC_GFOFFSET3D(u,1,1,0)) + 405*(KRANC_GFOFFSET3D(u,-1,2,0) + KRANC_GFOFFSET3D(u,1,-2,0) + KRANC_GFOFFSET3D(u,-2,1,0) + KRANC_GFOFFSET3D(u,2,-1,0)) - 405*(KRANC_GFOFFSET3D(u,-1,-2,0) + KRANC_GFOFFSET3D(u,1,2,0) + KRANC_GFOFFSET3D(u,-2,-1,0) + KRANC_GFOFFSET3D(u,2,1,0)) - 81*(KRANC_GFOFFSET3D(u,-2,2,0) + KRANC_GFOFFSET3D(u,2,-2,0)) + 81*(KRANC_GFOFFSET3D(u,-2,-2,0) + KRANC_GFOFFSET3D(u,2,2,0)) - 45*(KRANC_GFOFFSET3D(u,-1,3,0) + KRANC_GFOFFSET3D(u,1,-3,0) + KRANC_GFOFFSET3D(u,-3,1,0) + KRANC_GFOFFSET3D(u,3,-1,0)) + 45*(KRANC_GFOFFSET3D(u,-1,-3,0) + KRANC_GFOFFSET3D(u,1,3,0) + KRANC_GFOFFSET3D(u,-3,-1,0) + KRANC_GFOFFSET3D(u,3,1,0)) + 9*(KRANC_GFOFFSET3D(u,-2,3,0) + KRANC_GFOFFSET3D(u,2,-3,0) + KRANC_GFOFFSET3D(u,-3,2,0) + KRANC_GFOFFSET3D(u,3,-2,0)) - 9*(KRANC_GFOFFSET3D(u,-2,-3,0) + KRANC_GFOFFSET3D(u,2,3,0) + KRANC_GFOFFSET3D(u,-3,-2,0) + KRANC_GFOFFSET3D(u,3,2,0)) + KRANC_GFOFFSET3D(u,-3,-3,0) - KRANC_GFOFFSET3D(u,-3,3,0) - KRANC_GFOFFSET3D(u,3,-3,0) + KRANC_GFOFFSET3D(u,3,3,0))*p1o3600dxdy)
#else
#  define PDstandardfdOrder621(u) (PDstandardfdOrder621_impl(u,p1o3600dxdy,cdj,cdk))
static CCTK_REAL PDstandardfdOrder621_impl(const CCTK_REAL* restrict const u, const CCTK_REAL p1o3600dxdy, const ptrdiff_t cdj, const ptrdiff_t cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL PDstandardfdOrder621_impl(const CCTK_REAL* restrict const u, const CCTK_REAL p1o3600dxdy, const ptrdiff_t cdj, const ptrdiff_t cdk)
{
  const ptrdiff_t cdi CCTK_ATTRIBUTE_UNUSED = sizeof(CCTK_REAL);
  return PDstandardfdOrder612_impl(u, p1o3600dxdy, cdj, cdk);
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDstandardfdOrder623(u) ((-2025*(KRANC_GFOFFSET3D(u,0,-1,1) + KRANC_GFOFFSET3D(u,0,1,-1)) + 2025*(KRANC_GFOFFSET3D(u,0,-1,-1) + KRANC_GFOFFSET3D(u,0,1,1)) + 405*(KRANC_GFOFFSET3D(u,0,-1,2) + KRANC_GFOFFSET3D(u,0,1,-2) + KRANC_GFOFFSET3D(u,0,-2,1) + KRANC_GFOFFSET3D(u,0,2,-1)) - 405*(KRANC_GFOFFSET3D(u,0,-1,-2) + KRANC_GFOFFSET3D(u,0,1,2) + KRANC_GFOFFSET3D(u,0,-2,-1) + KRANC_GFOFFSET3D(u,0,2,1)) - 81*(KRANC_GFOFFSET3D(u,0,-2,2) + KRANC_GFOFFSET3D(u,0,2,-2)) + 81*(KRANC_GFOFFSET3D(u,0,-2,-2) + KRANC_GFOFFSET3D(u,0,2,2)) - 45*(KRANC_GFOFFSET3D(u,0,-1,3) + KRANC_GFOFFSET3D(u,0,1,-3) + KRANC_GFOFFSET3D(u,0,-3,1) + KRANC_GFOFFSET3D(u,0,3,-1)) + 45*(KRANC_GFOFFSET3D(u,0,-1,-3) + KRANC_GFOFFSET3D(u,0,1,3) + KRANC_GFOFFSET3D(u,0,-3,-1) + KRANC_GFOFFSET3D(u,0,3,1)) + 9*(KRANC_GFOFFSET3D(u,0,-2,3) + KRANC_GFOFFSET3D(u,0,2,-3) + KRANC_GFOFFSET3D(u,0,-3,2) + KRANC_GFOFFSET3D(u,0,3,-2)) - 9*(KRANC_GFOFFSET3D(u,0,-2,-3) + KRANC_GFOFFSET3D(u,0,2,3) + KRANC_GFOFFSET3D(u,0,-3,-2) + KRANC_GFOFFSET3D(u,0,3,2)) + KRANC_GFOFFSET3D(u,0,-3,-3) - KRANC_GFOFFSET3D(u,0,-3,3) - KRANC_GFOFFSET3D(u,0,3,-3) + KRANC_GFOFFSET3D(u,0,3,3))*p1o3600dydz)
#else
#  define PDstandardfdOrder623(u) (PDstandardfdOrder623_impl(u,p1o3600dydz,cdj,cdk))
static CCTK_REAL PDstandardfdOrder623_impl(const CCTK_REAL* restrict const u, const CCTK_REAL p1o3600dydz, const ptrdiff_t cdj, const ptrdiff_t cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL PDstandardfdOrder623_impl(const CCTK_REAL* restrict const u, const CCTK_REAL p1o3600dydz, const ptrdiff_t cdj, const ptrdiff_t cdk)
{
  const ptrdiff_t cdi CCTK_ATTRIBUTE_UNUSED = sizeof(CCTK_REAL);
  return (-2025*(KRANC_GFOFFSET3D(u,0,-1,1) + KRANC_GFOFFSET3D(u,0,1,-1)) + 2025*(KRANC_GFOFFSET3D(u,0,-1,-1) + KRANC_GFOFFSET3D(u,0,1,1)) + 405*(KRANC_GFOFFSET3D(u,0,-1,2) + KRANC_GFOFFSET3D(u,0,1,-2) + KRANC_GFOFFSET3D(u,0,-2,1) + KRANC_GFOFFSET3D(u,0,2,-1)) - 405*(KRANC_GFOFFSET3D(u,0,-1,-2) + KRANC_GFOFFSET3D(u,0,1,2) + KRANC_GFOFFSET3D(u,0,-2,-1) + KRANC_GFOFFSET3D(u,0,2,1)) - 81*(KRANC_GFOFFSET3D(u,0,-2,2) + KRANC_GFOFFSET3D(u,0,2,-2)) + 81*(KRANC_GFOFFSET3D(u,0,-2,-2) + KRANC_GFOFFSET3D(u,0,2,2)) - 45*(KRANC_GFOFFSET3D(u,0,-1,3) + KRANC_GFOFFSET3D(u,0,1,-3) + KRANC_GFOFFSET3D(u,0,-3,1) + KRANC_GFOFFSET3D(u,0,3,-1)) + 45*(KRANC_GFOFFSET3D(u,0,-1,-3) + KRANC_GFOFFSET3D(u,0,1,3) + KRANC_GFOFFSET3D(u,0,-3,-1) + KRANC_GFOFFSET3D(u,0,3,1)) + 9*(KRANC_GFOFFSET3D(u,0,-2,3) + KRANC_GFOFFSET3D(u,0,2,-3) + KRANC_GFOFFSET3D(u,0,-3,2) + KRANC_GFOFFSET3D(u,0,3,-2)) - 9*(KRANC_GFOFFSET3D(u,0,-2,-3) + KRANC_GFOFFSET3D(u,0,2,3) + KRANC_GFOFFSET3D(u,0,-3,-2) + KRANC_GFOFFSET3D(u,0,3,2)) + KRANC_GFOFFSET3D(u,0,-3,-3) - KRANC_GFOFFSET3D(u,0,-3,3) - KRANC_GFOFFSET3D(u,0,3,-3) + KRANC_GFOFFSET3D(u,0,3,3))*p1o3600dydz;
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDstandardfdOrder631(u) ((-2025*(KRANC_GFOFFSET3D(u,-1,0,1) + KRANC_GFOFFSET3D(u,1,0,-1)) + 2025*(KRANC_GFOFFSET3D(u,-1,0,-1) + KRANC_GFOFFSET3D(u,1,0,1)) + 405*(KRANC_GFOFFSET3D(u,-1,0,2) + KRANC_GFOFFSET3D(u,1,0,-2) + KRANC_GFOFFSET3D(u,-2,0,1) + KRANC_GFOFFSET3D(u,2,0,-1)) - 405*(KRANC_GFOFFSET3D(u,-1,0,-2) + KRANC_GFOFFSET3D(u,1,0,2) + KRANC_GFOFFSET3D(u,-2,0,-1) + KRANC_GFOFFSET3D(u,2,0,1)) - 81*(KRANC_GFOFFSET3D(u,-2,0,2) + KRANC_GFOFFSET3D(u,2,0,-2)) + 81*(KRANC_GFOFFSET3D(u,-2,0,-2) + KRANC_GFOFFSET3D(u,2,0,2)) - 45*(KRANC_GFOFFSET3D(u,-1,0,3) + KRANC_GFOFFSET3D(u,1,0,-3) + KRANC_GFOFFSET3D(u,-3,0,1) + KRANC_GFOFFSET3D(u,3,0,-1)) + 45*(KRANC_GFOFFSET3D(u,-1,0,-3) + KRANC_GFOFFSET3D(u,1,0,3) + KRANC_GFOFFSET3D(u,-3,0,-1) + KRANC_GFOFFSET3D(u,3,0,1)) + 9*(KRANC_GFOFFSET3D(u,-2,0,3) + KRANC_GFOFFSET3D(u,2,0,-3) + KRANC_GFOFFSET3D(u,-3,0,2) + KRANC_GFOFFSET3D(u,3,0,-2)) - 9*(KRANC_GFOFFSET3D(u,-2,0,-3) + KRANC_GFOFFSET3D(u,2,0,3) + KRANC_GFOFFSET3D(u,-3,0,-2) + KRANC_GFOFFSET3D(u,3,0,2)) + KRANC_GFOFFSET3D(u,-3,0,-3) - KRANC_GFOFFSET3D(u,-3,0,3) - KRANC_GFOFFSET3D(u,3,0,-3) + KRANC_GFOFFSET3D(u,3,0,3))*p1o3600dxdz)
#else
#  define PDstandardfdOrder631(u) (PDstandardfdOrder631_impl(u,p1o3600dxdz,cdj,cdk))
static CCTK_REAL PDstandardfdOrder631_impl(const CCTK_REAL* restrict const u, const CCTK_REAL p1o3600dxdz, const ptrdiff_t cdj, const ptrdiff_t cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL PDstandardfdOrder631_impl(const CCTK_REAL* restrict const u, const CCTK_REAL p1o3600dxdz, const ptrdiff_t cdj, const ptrdiff_t cdk)
{
  const ptrdiff_t cdi CCTK_ATTRIBUTE_UNUSED = sizeof(CCTK_REAL);
  return PDstandardfdOrder612_impl(u, p1o3600dxdz, cdk, cdj);
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDstandardfdOrder632(u) ((-2025*(KRANC_GFOFFSET3D(u,0,-1,1) + KRANC_GFOFFSET3D(u,0,1,-1)) + 2025*(KRANC_GFOFFSET3D(u,0,-1,-1) + KRANC_GFOFFSET3D(u,0,1,1)) + 405*(KRANC_GFOFFSET3D(u,0,-1,2) + KRANC_GFOFFSET3D(u,0,1,-2) + KRANC_GFOFFSET3D(u,0,-2,1) + KRANC_GFOFFSET3D(u,0,2,-1)) - 405*(KRANC_GFOFFSET3D(u,0,-1,-2) + KRANC_GFOFFSET3D(u,0,1,2) + KRANC_GFOFFSET3D(u,0,-2,-1) + KRANC_GFOFFSET3D(u,0,2,1)) - 81*(KRANC_GFOFFSET3D(u,0,-2,2) + KRANC_GFOFFSET3D(u,0,2,-2)) + 81*(KRANC_GFOFFSET3D(u,0,-2,-2) + KRANC_GFOFFSET3D(u,0,2,2)) - 45*(KRANC_GFOFFSET3D(u,0,-1,3) + KRANC_GFOFFSET3D(u,0,1,-3) + KRANC_GFOFFSET3D(u,0,-3,1) + KRANC_GFOFFSET3D(u,0,3,-1)) + 45*(KRANC_GFOFFSET3D(u,0,-1,-3) + KRANC_GFOFFSET3D(u,0,1,3) + KRANC_GFOFFSET3D(u,0,-3,-1) + KRANC_GFOFFSET3D(u,0,3,1)) + 9*(KRANC_GFOFFSET3D(u,0,-2,3) + KRANC_GFOFFSET3D(u,0,2,-3) + KRANC_GFOFFSET3D(u,0,-3,2) + KRANC_GFOFFSET3D(u,0,3,-2)) - 9*(KRANC_GFOFFSET3D(u,0,-2,-3) + KRANC_GFOFFSET3D(u,0,2,3) + KRANC_GFOFFSET3D(u,0,-3,-2) + KRANC_GFOFFSET3D(u,0,3,2)) + KRANC_GFOFFSET3D(u,0,-3,-3) - KRANC_GFOFFSET3D(u,0,-3,3) - KRANC_GFOFFSET3D(u,0,3,-3) + KRANC_GFOFFSET3D(u,0,3,3))*p1o3600dydz)
#else
#  define PDstandardfdOrder632(u) (PDstandardfdOrder632_impl(u,p1o3600dydz,cdj,cdk))
static CCTK_REAL PDstandardfdOrder632_impl(const CCTK_REAL* restrict const u, const CCTK_REAL p1o3600dydz, const ptrdiff_t cdj, const ptrdiff_t cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL PDstandardfdOrder632_impl(const CCTK_REAL* restrict const u, const CCTK_REAL p1o3600dydz, const ptrdiff_t cdj, const ptrdiff_t cdk)
{
  const ptrdiff_t cdi CCTK_ATTRIBUTE_UNUSED = sizeof(CCTK_REAL);
  return PDstandardfdOrder623_impl(u, p1o3600dydz, cdj, cdk);
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDstandardfdOrder812(u) ((-451584*(KRANC_GFOFFSET3D(u,-1,1,0) + KRANC_GFOFFSET3D(u,1,-1,0)) + 451584*(KRANC_GFOFFSET3D(u,-1,-1,0) + KRANC_GFOFFSET3D(u,1,1,0)) + 112896*(KRANC_GFOFFSET3D(u,-1,2,0) + KRANC_GFOFFSET3D(u,1,-2,0) + KRANC_GFOFFSET3D(u,-2,1,0) + KRANC_GFOFFSET3D(u,2,-1,0)) - 112896*(KRANC_GFOFFSET3D(u,-1,-2,0) + KRANC_GFOFFSET3D(u,1,2,0) + KRANC_GFOFFSET3D(u,-2,-1,0) + KRANC_GFOFFSET3D(u,2,1,0)) - 28224*(KRANC_GFOFFSET3D(u,-2,2,0) + KRANC_GFOFFSET3D(u,2,-2,0)) + 28224*(KRANC_GFOFFSET3D(u,-2,-2,0) + KRANC_GFOFFSET3D(u,2,2,0)) - 21504*(KRANC_GFOFFSET3D(u,-1,3,0) + KRANC_GFOFFSET3D(u,1,-3,0) + KRANC_GFOFFSET3D(u,-3,1,0) + KRANC_GFOFFSET3D(u,3,-1,0)) + 21504*(KRANC_GFOFFSET3D(u,-1,-3,0) + KRANC_GFOFFSET3D(u,1,3,0) + KRANC_GFOFFSET3D(u,-3,-1,0) + KRANC_GFOFFSET3D(u,3,1,0)) + 5376*(KRANC_GFOFFSET3D(u,-2,3,0) + KRANC_GFOFFSET3D(u,2,-3,0) + KRANC_GFOFFSET3D(u,-3,2,0) + KRANC_GFOFFSET3D(u,3,-2,0)) - 5376*(KRANC_GFOFFSET3D(u,-2,-3,0) + KRANC_GFOFFSET3D(u,2,3,0) + KRANC_GFOFFSET3D(u,-3,-2,0) + KRANC_GFOFFSET3D(u,3,2,0)) - 1024*(KRANC_GFOFFSET3D(u,-3,3,0) + KRANC_GFOFFSET3D(u,3,-3,0)) + 1024*(KRANC_GFOFFSET3D(u,-3,-3,0) + KRANC_GFOFFSET3D(u,3,3,0)) + 2016*(KRANC_GFOFFSET3D(u,-1,4,0) + KRANC_GFOFFSET3D(u,1,-4,0) + KRANC_GFOFFSET3D(u,-4,1,0) + KRANC_GFOFFSET3D(u,4,-1,0)) - 2016*(KRANC_GFOFFSET3D(u,-1,-4,0) + KRANC_GFOFFSET3D(u,1,4,0) + KRANC_GFOFFSET3D(u,-4,-1,0) + KRANC_GFOFFSET3D(u,4,1,0)) - 504*(KRANC_GFOFFSET3D(u,-2,4,0) + KRANC_GFOFFSET3D(u,2,-4,0) + KRANC_GFOFFSET3D(u,-4,2,0) + KRANC_GFOFFSET3D(u,4,-2,0)) + 504*(KRANC_GFOFFSET3D(u,-2,-4,0) + KRANC_GFOFFSET3D(u,2,4,0) + KRANC_GFOFFSET3D(u,-4,-2,0) + KRANC_GFOFFSET3D(u,4,2,0)) + 96*(KRANC_GFOFFSET3D(u,-3,4,0) + KRANC_GFOFFSET3D(u,3,-4,0) + KRANC_GFOFFSET3D(u,-4,3,0) + KRANC_GFOFFSET3D(u,4,-3,0)) - 96*(KRANC_GFOFFSET3D(u,-3,-4,0) + KRANC_GFOFFSET3D(u,3,4,0) + KRANC_GFOFFSET3D(u,-4,-3,0) + KRANC_GFOFFSET3D(u,4,3,0)) - 9*(KRANC_GFOFFSET3D(u,-4,4,0) + KRANC_GFOFFSET3D(u,4,-4,0)) + 9*(KRANC_GFOFFSET3D(u,-4,-4,0) + KRANC_GFOFFSET3D(u,4,4,0)))*p1o705600dxdy)
#else
#  define PDstandardfdOrder812(u) (PDstandardfdOrder812_impl(u,p1o705600dxdy,cdj,cdk))
static CCTK_REAL PDstandardfdOrder812_impl(const CCTK_REAL* restrict const u, const CCTK_REAL p1o705600dxdy, const ptrdiff_t cdj, const ptrdiff_t cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL PDstandardfdOrder812_impl(const CCTK_REAL* restrict const u, const CCTK_REAL p1o705600dxdy, const ptrdiff_t cdj, const ptrdiff_t cdk)
{
  const ptrdiff_t cdi CCTK_ATTRIBUTE_UNUSED = sizeof(CCTK_REAL);
  return (-451584*(KRANC_GFOFFSET3D(u,-1,1,0) + KRANC_GFOFFSET3D(u,1,-1,0)) + 451584*(KRANC_GFOFFSET3D(u,-1,-1,0) + KRANC_GFOFFSET3D(u,1,1,0)) + 112896*(KRANC_GFOFFSET3D(u,-1,2,0) + KRANC_GFOFFSET3D(u,1,-2,0) + KRANC_GFOFFSET3D(u,-2,1,0) + KRANC_GFOFFSET3D(u,2,-1,0)) - 112896*(KRANC_GFOFFSET3D(u,-1,-2,0) + KRANC_GFOFFSET3D(u,1,2,0) + KRANC_GFOFFSET3D(u,-2,-1,0) + KRANC_GFOFFSET3D(u,2,1,0)) - 28224*(KRANC_GFOFFSET3D(u,-2,2,0) + KRANC_GFOFFSET3D(u,2,-2,0)) + 28224*(KRANC_GFOFFSET3D(u,-2,-2,0) + KRANC_GFOFFSET3D(u,2,2,0)) - 21504*(KRANC_GFOFFSET3D(u,-1,3,0) + KRANC_GFOFFSET3D(u,1,-3,0) + KRANC_GFOFFSET3D(u,-3,1,0) + KRANC_GFOFFSET3D(u,3,-1,0)) + 21504*(KRANC_GFOFFSET3D(u,-1,-3,0) + KRANC_GFOFFSET3D(u,1,3,0) + KRANC_GFOFFSET3D(u,-3,-1,0) + KRANC_GFOFFSET3D(u,3,1,0)) + 5376*(KRANC_GFOFFSET3D(u,-2,3,0) + KRANC_GFOFFSET3D(u,2,-3,0) + KRANC_GFOFFSET3D(u,-3,2,0) + KRANC_GFOFFSET3D(u,3,-2,0)) - 5376*(KRANC_GFOFFSET3D(u,-2,-3,0) + KRANC_GFOFFSET3D(u,2,3,0) + KRANC_GFOFFSET3D(u,-3,-2,0) + KRANC_GFOFFSET3D(u,3,2,0)) - 1024*(KRANC_GFOFFSET3D(u,-3,3,0) + KRANC_GFOFFSET3D(u,3,-3,0)) + 1024*(KRANC_GFOFFSET3D(u,-3,-3,0) + KRANC_GFOFFSET3D(u,3,3,0)) + 2016*(KRANC_GFOFFSET3D(u,-1,4,0) + KRANC_GFOFFSET3D(u,1,-4,0) + KRANC_GFOFFSET3D(u,-4,1,0) + KRANC_GFOFFSET3D(u,4,-1,0)) - 2016*(KRANC_GFOFFSET3D(u,-1,-4,0) + KRANC_GFOFFSET3D(u,1,4,0) + KRANC_GFOFFSET3D(u,-4,-1,0) + KRANC_GFOFFSET3D(u,4,1,0)) - 504*(KRANC_GFOFFSET3D(u,-2,4,0) + KRANC_GFOFFSET3D(u,2,-4,0) + KRANC_GFOFFSET3D(u,-4,2,0) + KRANC_GFOFFSET3D(u,4,-2,0)) + 504*(KRANC_GFOFFSET3D(u,-2,-4,0) + KRANC_GFOFFSET3D(u,2,4,0) + KRANC_GFOFFSET3D(u,-4,-2,0) + KRANC_GFOFFSET3D(u,4,2,0)) + 96*(KRANC_GFOFFSET3D(u,-3,4,0) + KRANC_GFOFFSET3D(u,3,-4,0) + KRANC_GFOFFSET3D(u,-4,3,0) + KRANC_GFOFFSET3D(u,4,-3,0)) - 96*(KRANC_GFOFFSET3D(u,-3,-4,0) + KRANC_GFOFFSET3D(u,3,4,0) + KRANC_GFOFFSET3D(u,-4,-3,0) + KRANC_GFOFFSET3D(u,4,3,0)) - 9*(KRANC_GFOFFSET3D(u,-4,4,0) + KRANC_GFOFFSET3D(u,4,-4,0)) + 9*(KRANC_GFOFFSET3D(u,-4,-4,0) + KRANC_GFOFFSET3D(u,4,4,0)))*p1o705600dxdy;
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDstandardfdOrder813(u) ((-451584*(KRANC_GFOFFSET3D(u,-1,0,1) + KRANC_GFOFFSET3D(u,1,0,-1)) + 451584*(KRANC_GFOFFSET3D(u,-1,0,-1) + KRANC_GFOFFSET3D(u,1,0,1)) + 112896*(KRANC_GFOFFSET3D(u,-1,0,2) + KRANC_GFOFFSET3D(u,1,0,-2) + KRANC_GFOFFSET3D(u,-2,0,1) + KRANC_GFOFFSET3D(u,2,0,-1)) - 112896*(KRANC_GFOFFSET3D(u,-1,0,-2) + KRANC_GFOFFSET3D(u,1,0,2) + KRANC_GFOFFSET3D(u,-2,0,-1) + KRANC_GFOFFSET3D(u,2,0,1)) - 28224*(KRANC_GFOFFSET3D(u,-2,0,2) + KRANC_GFOFFSET3D(u,2,0,-2)) + 28224*(KRANC_GFOFFSET3D(u,-2,0,-2) + KRANC_GFOFFSET3D(u,2,0,2)) - 21504*(KRANC_GFOFFSET3D(u,-1,0,3) + KRANC_GFOFFSET3D(u,1,0,-3) + KRANC_GFOFFSET3D(u,-3,0,1) + KRANC_GFOFFSET3D(u,3,0,-1)) + 21504*(KRANC_GFOFFSET3D(u,-1,0,-3) + KRANC_GFOFFSET3D(u,1,0,3) + KRANC_GFOFFSET3D(u,-3,0,-1) + KRANC_GFOFFSET3D(u,3,0,1)) + 5376*(KRANC_GFOFFSET3D(u,-2,0,3) + KRANC_GFOFFSET3D(u,2,0,-3) + KRANC_GFOFFSET3D(u,-3,0,2) + KRANC_GFOFFSET3D(u,3,0,-2)) - 5376*(KRANC_GFOFFSET3D(u,-2,0,-3) + KRANC_GFOFFSET3D(u,2,0,3) + KRANC_GFOFFSET3D(u,-3,0,-2) + KRANC_GFOFFSET3D(u,3,0,2)) - 1024*(KRANC_GFOFFSET3D(u,-3,0,3) + KRANC_GFOFFSET3D(u,3,0,-3)) + 1024*(KRANC_GFOFFSET3D(u,-3,0,-3) + KRANC_GFOFFSET3D(u,3,0,3)) + 2016*(KRANC_GFOFFSET3D(u,-1,0,4) + KRANC_GFOFFSET3D(u,1,0,-4) + KRANC_GFOFFSET3D(u,-4,0,1) + KRANC_GFOFFSET3D(u,4,0,-1)) - 2016*(KRANC_GFOFFSET3D(u,-1,0,-4) + KRANC_GFOFFSET3D(u,1,0,4) + KRANC_GFOFFSET3D(u,-4,0,-1) + KRANC_GFOFFSET3D(u,4,0,1)) - 504*(KRANC_GFOFFSET3D(u,-2,0,4) + KRANC_GFOFFSET3D(u,2,0,-4) + KRANC_GFOFFSET3D(u,-4,0,2) + KRANC_GFOFFSET3D(u,4,0,-2)) + 504*(KRANC_GFOFFSET3D(u,-2,0,-4) + KRANC_GFOFFSET3D(u,2,0,4) + KRANC_GFOFFSET3D(u,-4,0,-2) + KRANC_GFOFFSET3D(u,4,0,2)) + 96*(KRANC_GFOFFSET3D(u,-3,0,4) + KRANC_GFOFFSET3D(u,3,0,-4) + KRANC_GFOFFSET3D(u,-4,0,3) + KRANC_GFOFFSET3D(u,4,0,-3)) - 96*(KRANC_GFOFFSET3D(u,-3,0,-4) + KRANC_GFOFFSET3D(u,3,0,4) + KRANC_GFOFFSET3D(u,-4,0,-3) + KRANC_GFOFFSET3D(u,4,0,3)) - 9*(KRANC_GFOFFSET3D(u,-4,0,4) + KRANC_GFOFFSET3D(u,4,0,-4)) + 9*(KRANC_GFOFFSET3D(u,-4,0,-4) + KRANC_GFOFFSET3D(u,4,0,4)))*p1o705600dxdz)
#else
#  define PDstandardfdOrder813(u) (PDstandardfdOrder813_impl(u,p1o705600dxdz,cdj,cdk))
static CCTK_REAL PDstandardfdOrder813_impl(const CCTK_REAL* restrict const u, const CCTK_REAL p1o705600dxdz, const ptrdiff_t cdj, const ptrdiff_t cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL PDstandardfdOrder813_impl(const CCTK_REAL* restrict const u, const CCTK_REAL p1o705600dxdz, const ptrdiff_t cdj, const ptrdiff_t cdk)
{
  const ptrdiff_t cdi CCTK_ATTRIBUTE_UNUSED = sizeof(CCTK_REAL);
  return PDstandardfdOrder812_impl(u, p1o705600dxdz, cdk, cdj);
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDstandardfdOrder821(u) ((-451584*(KRANC_GFOFFSET3D(u,-1,1,0) + KRANC_GFOFFSET3D(u,1,-1,0)) + 451584*(KRANC_GFOFFSET3D(u,-1,-1,0) + KRANC_GFOFFSET3D(u,1,1,0)) + 112896*(KRANC_GFOFFSET3D(u,-1,2,0) + KRANC_GFOFFSET3D(u,1,-2,0) + KRANC_GFOFFSET3D(u,-2,1,0) + KRANC_GFOFFSET3D(u,2,-1,0)) - 112896*(KRANC_GFOFFSET3D(u,-1,-2,0) + KRANC_GFOFFSET3D(u,1,2,0) + KRANC_GFOFFSET3D(u,-2,-1,0) + KRANC_GFOFFSET3D(u,2,1,0)) - 28224*(KRANC_GFOFFSET3D(u,-2,2,0) + KRANC_GFOFFSET3D(u,2,-2,0)) + 28224*(KRANC_GFOFFSET3D(u,-2,-2,0) + KRANC_GFOFFSET3D(u,2,2,0)) - 21504*(KRANC_GFOFFSET3D(u,-1,3,0) + KRANC_GFOFFSET3D(u,1,-3,0) + KRANC_GFOFFSET3D(u,-3,1,0) + KRANC_GFOFFSET3D(u,3,-1,0)) + 21504*(KRANC_GFOFFSET3D(u,-1,-3,0) + KRANC_GFOFFSET3D(u,1,3,0) + KRANC_GFOFFSET3D(u,-3,-1,0) + KRANC_GFOFFSET3D(u,3,1,0)) + 5376*(KRANC_GFOFFSET3D(u,-2,3,0) + KRANC_GFOFFSET3D(u,2,-3,0) + KRANC_GFOFFSET3D(u,-3,2,0) + KRANC_GFOFFSET3D(u,3,-2,0)) - 5376*(KRANC_GFOFFSET3D(u,-2,-3,0) + KRANC_GFOFFSET3D(u,2,3,0) + KRANC_GFOFFSET3D(u,-3,-2,0) + KRANC_GFOFFSET3D(u,3,2,0)) - 1024*(KRANC_GFOFFSET3D(u,-3,3,0) + KRANC_GFOFFSET3D(u,3,-3,0)) + 1024*(KRANC_GFOFFSET3D(u,-3,-3,0) + KRANC_GFOFFSET3D(u,3,3,0)) + 2016*(KRANC_GFOFFSET3D(u,-1,4,0) + KRANC_GFOFFSET3D(u,1,-4,0) + KRANC_GFOFFSET3D(u,-4,1,0) + KRANC_GFOFFSET3D(u,4,-1,0)) - 2016*(KRANC_GFOFFSET3D(u,-1,-4,0) + KRANC_GFOFFSET3D(u,1,4,0) + KRANC_GFOFFSET3D(u,-4,-1,0) + KRANC_GFOFFSET3D(u,4,1,0)) - 504*(KRANC_GFOFFSET3D(u,-2,4,0) + KRANC_GFOFFSET3D(u,2,-4,0) + KRANC_GFOFFSET3D(u,-4,2,0) + KRANC_GFOFFSET3D(u,4,-2,0)) + 504*(KRANC_GFOFFSET3D(u,-2,-4,0) + KRANC_GFOFFSET3D(u,2,4,0) + KRANC_GFOFFSET3D(u,-4,-2,0) + KRANC_GFOFFSET3D(u,4,2,0)) + 96*(KRANC_GFOFFSET3D(u,-3,4,0) + KRANC_GFOFFSET3D(u,3,-4,0) + KRANC_GFOFFSET3D(u,-4,3,0) + KRANC_GFOFFSET3D(u,4,-3,0)) - 96*(KRANC_GFOFFSET3D(u,-3,-4,0) + KRANC_GFOFFSET3D(u,3,4,0) + KRANC_GFOFFSET3D(u,-4,-3,0) + KRANC_GFOFFSET3D(u,4,3,0)) - 9*(KRANC_GFOFFSET3D(u,-4,4,0) + KRANC_GFOFFSET3D(u,4,-4,0)) + 9*(KRANC_GFOFFSET3D(u,-4,-4,0) + KRANC_GFOFFSET3D(u,4,4,0)))*p1o705600dxdy)
#else
#  define PDstandardfdOrder821(u) (PDstandardfdOrder821_impl(u,p1o705600dxdy,cdj,cdk))
static CCTK_REAL PDstandardfdOrder821_impl(const CCTK_REAL* restrict const u, const CCTK_REAL p1o705600dxdy, const ptrdiff_t cdj, const ptrdiff_t cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL PDstandardfdOrder821_impl(const CCTK_REAL* restrict const u, const CCTK_REAL p1o705600dxdy, const ptrdiff_t cdj, const ptrdiff_t cdk)
{
  const ptrdiff_t cdi CCTK_ATTRIBUTE_UNUSED = sizeof(CCTK_REAL);
  return PDstandardfdOrder812_impl(u, p1o705600dxdy, cdj, cdk);
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDstandardfdOrder823(u) ((-451584*(KRANC_GFOFFSET3D(u,0,-1,1) + KRANC_GFOFFSET3D(u,0,1,-1)) + 451584*(KRANC_GFOFFSET3D(u,0,-1,-1) + KRANC_GFOFFSET3D(u,0,1,1)) + 112896*(KRANC_GFOFFSET3D(u,0,-1,2) + KRANC_GFOFFSET3D(u,0,1,-2) + KRANC_GFOFFSET3D(u,0,-2,1) + KRANC_GFOFFSET3D(u,0,2,-1)) - 112896*(KRANC_GFOFFSET3D(u,0,-1,-2) + KRANC_GFOFFSET3D(u,0,1,2) + KRANC_GFOFFSET3D(u,0,-2,-1) + KRANC_GFOFFSET3D(u,0,2,1)) - 28224*(KRANC_GFOFFSET3D(u,0,-2,2) + KRANC_GFOFFSET3D(u,0,2,-2)) + 28224*(KRANC_GFOFFSET3D(u,0,-2,-2) + KRANC_GFOFFSET3D(u,0,2,2)) - 21504*(KRANC_GFOFFSET3D(u,0,-1,3) + KRANC_GFOFFSET3D(u,0,1,-3) + KRANC_GFOFFSET3D(u,0,-3,1) + KRANC_GFOFFSET3D(u,0,3,-1)) + 21504*(KRANC_GFOFFSET3D(u,0,-1,-3) + KRANC_GFOFFSET3D(u,0,1,3) + KRANC_GFOFFSET3D(u,0,-3,-1) + KRANC_GFOFFSET3D(u,0,3,1)) + 5376*(KRANC_GFOFFSET3D(u,0,-2,3) + KRANC_GFOFFSET3D(u,0,2,-3) + KRANC_GFOFFSET3D(u,0,-3,2) + KRANC_GFOFFSET3D(u,0,3,-2)) - 5376*(KRANC_GFOFFSET3D(u,0,-2,-3) + KRANC_GFOFFSET3D(u,0,2,3) + KRANC_GFOFFSET3D(u,0,-3,-2) + KRANC_GFOFFSET3D(u,0,3,2)) - 1024*(KRANC_GFOFFSET3D(u,0,-3,3) + KRANC_GFOFFSET3D(u,0,3,-3)) + 1024*(KRANC_GFOFFSET3D(u,0,-3,-3) + KRANC_GFOFFSET3D(u,0,3,3)) + 2016*(KRANC_GFOFFSET3D(u,0,-1,4) + KRANC_GFOFFSET3D(u,0,1,-4) + KRANC_GFOFFSET3D(u,0,-4,1) + KRANC_GFOFFSET3D(u,0,4,-1)) - 2016*(KRANC_GFOFFSET3D(u,0,-1,-4) + KRANC_GFOFFSET3D(u,0,1,4) + KRANC_GFOFFSET3D(u,0,-4,-1) + KRANC_GFOFFSET3D(u,0,4,1)) - 504*(KRANC_GFOFFSET3D(u,0,-2,4) + KRANC_GFOFFSET3D(u,0,2,-4) + KRANC_GFOFFSET3D(u,0,-4,2) + KRANC_GFOFFSET3D(u,0,4,-2)) + 504*(KRANC_GFOFFSET3D(u,0,-2,-4) + KRANC_GFOFFSET3D(u,0,2,4) + KRANC_GFOFFSET3D(u,0,-4,-2) + KRANC_GFOFFSET3D(u,0,4,2)) + 96*(KRANC_GFOFFSET3D(u,0,-3,4) + KRANC_GFOFFSET3D(u,0,3,-4) + KRANC_GFOFFSET3D(u,0,-4,3) + KRANC_GFOFFSET3D(u,0,4,-3)) - 96*(KRANC_GFOFFSET3D(u,0,-3,-4) + KRANC_GFOFFSET3D(u,0,3,4) + KRANC_GFOFFSET3D(u,0,-4,-3) + KRANC_GFOFFSET3D(u,0,4,3)) - 9*(KRANC_GFOFFSET3D(u,0,-4,4) + KRANC_GFOFFSET3D(u,0,4,-4)) + 9*(KRANC_GFOFFSET3D(u,0,-4,-4) + KRANC_GFOFFSET3D(u,0,4,4)))*p1o705600dydz)
#else
#  define PDstandardfdOrder823(u) (PDstandardfdOrder823_impl(u,p1o705600dydz,cdj,cdk))
static CCTK_REAL PDstandardfdOrder823_impl(const CCTK_REAL* restrict const u, const CCTK_REAL p1o705600dydz, const ptrdiff_t cdj, const ptrdiff_t cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL PDstandardfdOrder823_impl(const CCTK_REAL* restrict const u, const CCTK_REAL p1o705600dydz, const ptrdiff_t cdj, const ptrdiff_t cdk)
{
  const ptrdiff_t cdi CCTK_ATTRIBUTE_UNUSED = sizeof(CCTK_REAL);
  return (-451584*(KRANC_GFOFFSET3D(u,0,-1,1) + KRANC_GFOFFSET3D(u,0,1,-1)) + 451584*(KRANC_GFOFFSET3D(u,0,-1,-1) + KRANC_GFOFFSET3D(u,0,1,1)) + 112896*(KRANC_GFOFFSET3D(u,0,-1,2) + KRANC_GFOFFSET3D(u,0,1,-2) + KRANC_GFOFFSET3D(u,0,-2,1) + KRANC_GFOFFSET3D(u,0,2,-1)) - 112896*(KRANC_GFOFFSET3D(u,0,-1,-2) + KRANC_GFOFFSET3D(u,0,1,2) + KRANC_GFOFFSET3D(u,0,-2,-1) + KRANC_GFOFFSET3D(u,0,2,1)) - 28224*(KRANC_GFOFFSET3D(u,0,-2,2) + KRANC_GFOFFSET3D(u,0,2,-2)) + 28224*(KRANC_GFOFFSET3D(u,0,-2,-2) + KRANC_GFOFFSET3D(u,0,2,2)) - 21504*(KRANC_GFOFFSET3D(u,0,-1,3) + KRANC_GFOFFSET3D(u,0,1,-3) + KRANC_GFOFFSET3D(u,0,-3,1) + KRANC_GFOFFSET3D(u,0,3,-1)) + 21504*(KRANC_GFOFFSET3D(u,0,-1,-3) + KRANC_GFOFFSET3D(u,0,1,3) + KRANC_GFOFFSET3D(u,0,-3,-1) + KRANC_GFOFFSET3D(u,0,3,1)) + 5376*(KRANC_GFOFFSET3D(u,0,-2,3) + KRANC_GFOFFSET3D(u,0,2,-3) + KRANC_GFOFFSET3D(u,0,-3,2) + KRANC_GFOFFSET3D(u,0,3,-2)) - 5376*(KRANC_GFOFFSET3D(u,0,-2,-3) + KRANC_GFOFFSET3D(u,0,2,3) + KRANC_GFOFFSET3D(u,0,-3,-2) + KRANC_GFOFFSET3D(u,0,3,2)) - 1024*(KRANC_GFOFFSET3D(u,0,-3,3) + KRANC_GFOFFSET3D(u,0,3,-3)) + 1024*(KRANC_GFOFFSET3D(u,0,-3,-3) + KRANC_GFOFFSET3D(u,0,3,3)) + 2016*(KRANC_GFOFFSET3D(u,0,-1,4) + KRANC_GFOFFSET3D(u,0,1,-4) + KRANC_GFOFFSET3D(u,0,-4,1) + KRANC_GFOFFSET3D(u,0,4,-1)) - 2016*(KRANC_GFOFFSET3D(u,0,-1,-4) + KRANC_GFOFFSET3D(u,0,1,4) + KRANC_GFOFFSET3D(u,0,-4,-1) + KRANC_GFOFFSET3D(u,0,4,1)) - 504*(KRANC_GFOFFSET3D(u,0,-2,4) + KRANC_GFOFFSET3D(u,0,2,-4) + KRANC_GFOFFSET3D(u,0,-4,2) + KRANC_GFOFFSET3D(u,0,4,-2)) + 504*(KRANC_GFOFFSET3D(u,0,-2,-4) + KRANC_GFOFFSET3D(u,0,2,4) + KRANC_GFOFFSET3D(u,0,-4,-2) + KRANC_GFOFFSET3D(u,0,4,2)) + 96*(KRANC_GFOFFSET3D(u,0,-3,4) + KRANC_GFOFFSET3D(u,0,3,-4) + KRANC_GFOFFSET3D(u,0,-4,3) + KRANC_GFOFFSET3D(u,0,4,-3)) - 96*(KRANC_GFOFFSET3D(u,0,-3,-4) + KRANC_GFOFFSET3D(u,0,3,4) + KRANC_GFOFFSET3D(u,0,-4,-3) + KRANC_GFOFFSET3D(u,0,4,3)) - 9*(KRANC_GFOFFSET3D(u,0,-4,4) + KRANC_GFOFFSET3D(u,0,4,-4)) + 9*(KRANC_GFOFFSET3D(u,0,-4,-4) + KRANC_GFOFFSET3D(u,0,4,4)))*p1o705600dydz;
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDstandardfdOrder831(u) ((-451584*(KRANC_GFOFFSET3D(u,-1,0,1) + KRANC_GFOFFSET3D(u,1,0,-1)) + 451584*(KRANC_GFOFFSET3D(u,-1,0,-1) + KRANC_GFOFFSET3D(u,1,0,1)) + 112896*(KRANC_GFOFFSET3D(u,-1,0,2) + KRANC_GFOFFSET3D(u,1,0,-2) + KRANC_GFOFFSET3D(u,-2,0,1) + KRANC_GFOFFSET3D(u,2,0,-1)) - 112896*(KRANC_GFOFFSET3D(u,-1,0,-2) + KRANC_GFOFFSET3D(u,1,0,2) + KRANC_GFOFFSET3D(u,-2,0,-1) + KRANC_GFOFFSET3D(u,2,0,1)) - 28224*(KRANC_GFOFFSET3D(u,-2,0,2) + KRANC_GFOFFSET3D(u,2,0,-2)) + 28224*(KRANC_GFOFFSET3D(u,-2,0,-2) + KRANC_GFOFFSET3D(u,2,0,2)) - 21504*(KRANC_GFOFFSET3D(u,-1,0,3) + KRANC_GFOFFSET3D(u,1,0,-3) + KRANC_GFOFFSET3D(u,-3,0,1) + KRANC_GFOFFSET3D(u,3,0,-1)) + 21504*(KRANC_GFOFFSET3D(u,-1,0,-3) + KRANC_GFOFFSET3D(u,1,0,3) + KRANC_GFOFFSET3D(u,-3,0,-1) + KRANC_GFOFFSET3D(u,3,0,1)) + 5376*(KRANC_GFOFFSET3D(u,-2,0,3) + KRANC_GFOFFSET3D(u,2,0,-3) + KRANC_GFOFFSET3D(u,-3,0,2) + KRANC_GFOFFSET3D(u,3,0,-2)) - 5376*(KRANC_GFOFFSET3D(u,-2,0,-3) + KRANC_GFOFFSET3D(u,2,0,3) + KRANC_GFOFFSET3D(u,-3,0,-2) + KRANC_GFOFFSET3D(u,3,0,2)) - 1024*(KRANC_GFOFFSET3D(u,-3,0,3) + KRANC_GFOFFSET3D(u,3,0,-3)) + 1024*(KRANC_GFOFFSET3D(u,-3,0,-3) + KRANC_GFOFFSET3D(u,3,0,3)) + 2016*(KRANC_GFOFFSET3D(u,-1,0,4) + KRANC_GFOFFSET3D(u,1,0,-4) + KRANC_GFOFFSET3D(u,-4,0,1) + KRANC_GFOFFSET3D(u,4,0,-1)) - 2016*(KRANC_GFOFFSET3D(u,-1,0,-4) + KRANC_GFOFFSET3D(u,1,0,4) + KRANC_GFOFFSET3D(u,-4,0,-1) + KRANC_GFOFFSET3D(u,4,0,1)) - 504*(KRANC_GFOFFSET3D(u,-2,0,4) + KRANC_GFOFFSET3D(u,2,0,-4) + KRANC_GFOFFSET3D(u,-4,0,2) + KRANC_GFOFFSET3D(u,4,0,-2)) + 504*(KRANC_GFOFFSET3D(u,-2,0,-4) + KRANC_GFOFFSET3D(u,2,0,4) + KRANC_GFOFFSET3D(u,-4,0,-2) + KRANC_GFOFFSET3D(u,4,0,2)) + 96*(KRANC_GFOFFSET3D(u,-3,0,4) + KRANC_GFOFFSET3D(u,3,0,-4) + KRANC_GFOFFSET3D(u,-4,0,3) + KRANC_GFOFFSET3D(u,4,0,-3)) - 96*(KRANC_GFOFFSET3D(u,-3,0,-4) + KRANC_GFOFFSET3D(u,3,0,4) + KRANC_GFOFFSET3D(u,-4,0,-3) + KRANC_GFOFFSET3D(u,4,0,3)) - 9*(KRANC_GFOFFSET3D(u,-4,0,4) + KRANC_GFOFFSET3D(u,4,0,-4)) + 9*(KRANC_GFOFFSET3D(u,-4,0,-4) + KRANC_GFOFFSET3D(u,4,0,4)))*p1o705600dxdz)
#else
#  define PDstandardfdOrder831(u) (PDstandardfdOrder831_impl(u,p1o705600dxdz,cdj,cdk))
static CCTK_REAL PDstandardfdOrder831_impl(const CCTK_REAL* restrict const u, const CCTK_REAL p1o705600dxdz, const ptrdiff_t cdj, const ptrdiff_t cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL PDstandardfdOrder831_impl(const CCTK_REAL* restrict const u, const CCTK_REAL p1o705600dxdz, const ptrdiff_t cdj, const ptrdiff_t cdk)
{
  const ptrdiff_t cdi CCTK_ATTRIBUTE_UNUSED = sizeof(CCTK_REAL);
  return PDstandardfdOrder812_impl(u, p1o705600dxdz, cdk, cdj);
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDstandardfdOrder832(u) ((-451584*(KRANC_GFOFFSET3D(u,0,-1,1) + KRANC_GFOFFSET3D(u,0,1,-1)) + 451584*(KRANC_GFOFFSET3D(u,0,-1,-1) + KRANC_GFOFFSET3D(u,0,1,1)) + 112896*(KRANC_GFOFFSET3D(u,0,-1,2) + KRANC_GFOFFSET3D(u,0,1,-2) + KRANC_GFOFFSET3D(u,0,-2,1) + KRANC_GFOFFSET3D(u,0,2,-1)) - 112896*(KRANC_GFOFFSET3D(u,0,-1,-2) + KRANC_GFOFFSET3D(u,0,1,2) + KRANC_GFOFFSET3D(u,0,-2,-1) + KRANC_GFOFFSET3D(u,0,2,1)) - 28224*(KRANC_GFOFFSET3D(u,0,-2,2) + KRANC_GFOFFSET3D(u,0,2,-2)) + 28224*(KRANC_GFOFFSET3D(u,0,-2,-2) + KRANC_GFOFFSET3D(u,0,2,2)) - 21504*(KRANC_GFOFFSET3D(u,0,-1,3) + KRANC_GFOFFSET3D(u,0,1,-3) + KRANC_GFOFFSET3D(u,0,-3,1) + KRANC_GFOFFSET3D(u,0,3,-1)) + 21504*(KRANC_GFOFFSET3D(u,0,-1,-3) + KRANC_GFOFFSET3D(u,0,1,3) + KRANC_GFOFFSET3D(u,0,-3,-1) + KRANC_GFOFFSET3D(u,0,3,1)) + 5376*(KRANC_GFOFFSET3D(u,0,-2,3) + KRANC_GFOFFSET3D(u,0,2,-3) + KRANC_GFOFFSET3D(u,0,-3,2) + KRANC_GFOFFSET3D(u,0,3,-2)) - 5376*(KRANC_GFOFFSET3D(u,0,-2,-3) + KRANC_GFOFFSET3D(u,0,2,3) + KRANC_GFOFFSET3D(u,0,-3,-2) + KRANC_GFOFFSET3D(u,0,3,2)) - 1024*(KRANC_GFOFFSET3D(u,0,-3,3) + KRANC_GFOFFSET3D(u,0,3,-3)) + 1024*(KRANC_GFOFFSET3D(u,0,-3,-3) + KRANC_GFOFFSET3D(u,0,3,3)) + 2016*(KRANC_GFOFFSET3D(u,0,-1,4) + KRANC_GFOFFSET3D(u,0,1,-4) + KRANC_GFOFFSET3D(u,0,-4,1) + KRANC_GFOFFSET3D(u,0,4,-1)) - 2016*(KRANC_GFOFFSET3D(u,0,-1,-4) + KRANC_GFOFFSET3D(u,0,1,4) + KRANC_GFOFFSET3D(u,0,-4,-1) + KRANC_GFOFFSET3D(u,0,4,1)) - 504*(KRANC_GFOFFSET3D(u,0,-2,4) + KRANC_GFOFFSET3D(u,0,2,-4) + KRANC_GFOFFSET3D(u,0,-4,2) + KRANC_GFOFFSET3D(u,0,4,-2)) + 504*(KRANC_GFOFFSET3D(u,0,-2,-4) + KRANC_GFOFFSET3D(u,0,2,4) + KRANC_GFOFFSET3D(u,0,-4,-2) + KRANC_GFOFFSET3D(u,0,4,2)) + 96*(KRANC_GFOFFSET3D(u,0,-3,4) + KRANC_GFOFFSET3D(u,0,3,-4) + KRANC_GFOFFSET3D(u,0,-4,3) + KRANC_GFOFFSET3D(u,0,4,-3)) - 96*(KRANC_GFOFFSET3D(u,0,-3,-4) + KRANC_GFOFFSET3D(u,0,3,4) + KRANC_GFOFFSET3D(u,0,-4,-3) + KRANC_GFOFFSET3D(u,0,4,3)) - 9*(KRANC_GFOFFSET3D(u,0,-4,4) + KRANC_GFOFFSET3D(u,0,4,-4)) + 9*(KRANC_GFOFFSET3D(u,0,-4,-4) + KRANC_GFOFFSET3D(u,0,4,4)))*p1o705600dydz)
#else
#  define PDstandardfdOrder832(u) (PDstandardfdOrder832_impl(u,p1o705600dydz,cdj,cdk))
static CCTK_REAL PDstandardfdOrder832_impl(const CCTK_REAL* restrict const u, const CCTK_REAL p1o705600dydz, const ptrdiff_t cdj, const ptrdiff_t cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL PDstandardfdOrder832_impl(const CCTK_REAL* restrict const u, const CCTK_REAL p1o705600dydz, const ptrdiff_t cdj, const ptrdiff_t cdk)
{
  const ptrdiff_t cdi CCTK_ATTRIBUTE_UNUSED = sizeof(CCTK_REAL);
  return PDstandardfdOrder823_impl(u, p1o705600dydz, cdj, cdk);
}
#endif

