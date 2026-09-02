// Generic cartoon boundary condition: CarpetX calls
// Cartoon2DX_ApplyBoundary from apply_boundary_conditions for each
// group whose faces are boundary_t::cartoon. Classify the group as
// scalar / U / DD_sym from the Einstein Toolkit `tensortypealias` tag
// (same tag old Cartoon2D reads) or from the group's `parities` tag,
// and fill y-ghosts (and x<0 when fill_negative_x) by
// rotation+interpolation.
//
// Mixed centering is dispatched from the group's CENTERING table onto
// fill_ghosts_*<CI,CJ,CK>. Packed groups (AsterX cons_vector: scalars
// plus a momentum triple in one GF group) are walked component-wise
// from the parities array.

#include <cartoon2dx_carpetx.hxx>

#include <loop_device.hxx>

#include <cctk.h>
#include <cctk_Parameters.h>
#include <util_Table.h>
#include <cctk_Groups.h>

#include <cstdlib>
#include <cstring>
#include <vector>

namespace Cartoon2DX {
namespace {

const char *tensortype_name(tensortype_t tt) {
  switch (tt) {
  case tensortype_t::scalar:
    return "scalar";
  case tensortype_t::u:
    return "U";
  case tensortype_t::ddsym:
    return "DD_sym";
  }
  return "?";
}

bool read_centering(int gi, int &CI, int &CJ, int &CK) {
  CI = CJ = CK = 0;
  const int centering = CCTK_GroupCenteringTableI(gi);
  if (centering < 0)
    return true;
  CCTK_INT index[3] = {0, 0, 0};
  const int iret = Util_TableGetIntArray(centering, 3, index, "centering");
  if (iret == UTIL_ERROR_TABLE_NO_SUCH_KEY)
    return true;
  if (iret != 3)
    return false;
  CI = int(index[0]);
  CJ = int(index[1]);
  CK = int(index[2]);
  return (CI == 0 || CI == 1) && (CJ == 0 || CJ == 1) && (CK == 0 || CK == 1);
}

// Same strings old Cartoon2D/ApplyCartoon.c accepts, plus the mixed-case
// ADMBase spellings. CCTK_EQUALS is case-insensitive; the extra branches
// are documentation.
bool classify_alias(const char *alias, int numvars, tensortype_t &tt) {
  if (CCTK_EQUALS(alias, "scalar")) {
    tt = tensortype_t::scalar;
    return true;
  }
  if (CCTK_EQUALS(alias, "u") || CCTK_EQUALS(alias, "d")) {
    tt = tensortype_t::u;
    return numvars == 3;
  }
  if (CCTK_EQUALS(alias, "dd_sym") || CCTK_EQUALS(alias, "ddsym")) {
    tt = tensortype_t::ddsym;
    return numvars == 6;
  }
  return false;
}

int group_parities(int gi, std::vector<CCTK_INT> &par) {
  par.clear();
  const int tags = CCTK_GroupTagsTableI(gi);
  if (tags < 0)
    return 0;
  const int npar = Util_TableGetIntArray(tags, 0, nullptr, "parities");
  if (npar == UTIL_ERROR_TABLE_NO_SUCH_KEY || npar <= 0)
    return 0;
  par.resize(npar);
  const int iret = Util_TableGetIntArray(tags, npar, par.data(), "parities");
  if (iret != npar)
    CCTK_VERROR("Cartoon2DX: parities array for group %s has length %d, "
                "read %d",
                CCTK_FullGroupName(gi), npar, iret);
  return npar;
}

// tensortypealias is the Einstein Toolkit tag (ADMBase, old Cartoon2D).
// parities is what CarpetX/Cottonmouth already put on every tensor group.
bool classify_group(int gi, int numvars, tensortype_t &tt) {
  const int tags = CCTK_GroupTagsTableI(gi);
  if (tags < 0)
    return false;

  char alias[64] = "";
  const int ierr =
      Util_TableGetString(tags, sizeof alias, alias, "tensortypealias");
  if (ierr >= 0)
    return classify_alias(alias, numvars, tt);

  // Query length the same way CarpetX::get_group_parities does.
  const int npar = Util_TableGetIntArray(tags, 0, nullptr, "parities");
  if (npar == UTIL_ERROR_TABLE_NO_SUCH_KEY || npar <= 0)
    return false;
  if (numvars == 1 && npar == 3) {
    tt = tensortype_t::scalar;
    return true;
  }
  if (numvars == 3 && npar == 9) {
    tt = tensortype_t::u;
    return true;
  }
  if (numvars == 6 && npar == 18) {
    tt = tensortype_t::ddsym;
    return true;
  }
  return false;
}

// Polar (shift/vel/mom) or axial (B) Cartesian triple in a packed group.
bool is_vector_triple(const std::vector<CCTK_INT> &par, int c) {
  if (3 * (c + 3) > int(par.size()))
    return false;
  const CCTK_INT *p0 = par.data() + 3 * c;
  const CCTK_INT *p1 = par.data() + 3 * (c + 1);
  const CCTK_INT *p2 = par.data() + 3 * (c + 2);
  const bool polar = p0[0] == -1 && p0[1] == +1 && p0[2] == +1 &&
                     p1[0] == +1 && p1[1] == -1 && p1[2] == +1 &&
                     p2[0] == +1 && p2[1] == +1 && p2[2] == -1;
  const bool axial = p0[0] == +1 && p0[1] == -1 && p0[2] == -1 &&
                     p1[0] == -1 && p1[1] == +1 && p1[2] == -1 &&
                     p2[0] == -1 && p2[1] == -1 && p2[2] == +1;
  return polar || axial;
}

void log_group_once(int gi, const char *gname, const char *how, int numvars,
                    bool verbose) {
  if (!verbose)
    return;
  static std::vector<char> logged;
  const int ng = CCTK_NumGroups();
  if (gi < 0 || gi >= ng)
    return;
  if (int(logged.size()) < ng)
    logged.resize(ng, 0);
  if (logged[gi])
    return;
  logged[gi] = 1;
  CCTK_VINFO("Cartoon2DX: filling group %s as %s (%d vars)", gname, how,
             numvars);
}

template <int CI, int CJ, int CK>
void fill_as(const cGH *cctkGH, const Loop::GridDescBaseDevice &grid, int gi,
             int tl, int order, bool fill_negative_x, tensortype_t tt, int v0) {
  const Loop::GF3D5layout layout(cctkGH, {CI, CJ, CK});
  auto gf = [&](int c) {
    CCTK_REAL *ptr =
        static_cast<CCTK_REAL *>(CCTK_VarDataPtrI(cctkGH, tl, v0 + c));
    if (!ptr)
      CCTK_VERROR("Cartoon2DX: missing data pointer for %s component %d tl=%d",
                  CCTK_FullGroupName(gi), c, tl);
    return Loop::GF3D5<CCTK_REAL>(layout, ptr);
  };

  switch (tt) {
  case tensortype_t::scalar:
    fill_ghosts_scalar<CI, CJ, CK>(cctkGH, grid, order, fill_negative_x, gf(0));
    break;
  case tensortype_t::u:
    fill_ghosts_vector<CI, CJ, CK>(cctkGH, grid, order, fill_negative_x, gf(0),
                                   gf(1), gf(2));
    break;
  case tensortype_t::ddsym:
    fill_ghosts_ddsym<CI, CJ, CK>(cctkGH, grid, order, fill_negative_x, gf(0),
                                  gf(1), gf(2), gf(3), gf(4), gf(5));
    break;
  }
}

template <int CI, int CJ, int CK>
void fill_packed(const cGH *cctkGH, const Loop::GridDescBaseDevice &grid,
                 int gi, int tl, int order, bool fill_negative_x, int v0,
                 int numvars, const std::vector<CCTK_INT> &par) {
  const Loop::GF3D5layout layout(cctkGH, {CI, CJ, CK});
  auto gf = [&](int c) {
    CCTK_REAL *ptr =
        static_cast<CCTK_REAL *>(CCTK_VarDataPtrI(cctkGH, tl, v0 + c));
    if (!ptr)
      CCTK_VERROR("Cartoon2DX: missing data pointer for %s component %d tl=%d",
                  CCTK_FullGroupName(gi), c, tl);
    return Loop::GF3D5<CCTK_REAL>(layout, ptr);
  };

  int c = 0;
  while (c < numvars) {
    if (c + 2 < numvars && is_vector_triple(par, c)) {
      fill_ghosts_vector<CI, CJ, CK>(cctkGH, grid, order, fill_negative_x,
                                     gf(c), gf(c + 1), gf(c + 2));
      c += 3;
    } else {
      fill_ghosts_scalar<CI, CJ, CK>(cctkGH, grid, order, fill_negative_x,
                                     gf(c));
      c += 1;
    }
  }
}

template <int CI, int CJ, int CK>
void apply_group_c(const cGH *cctkGH, const Loop::GridDescBaseDevice &grid,
                   int gi, int tl, int order, bool fill_negative_x, bool verbose,
                   const char *gname, int numvars) {
  tensortype_t tt;
  const int v0 = CCTK_FirstVarIndexI(gi);
  if (v0 < 0)
    return;
  if (!CCTK_VarDataPtrI(cctkGH, tl, v0))
    return;

  if (classify_group(gi, numvars, tt)) {
    log_group_once(gi, gname, tensortype_name(tt), numvars, verbose);
    if (tt == tensortype_t::scalar) {
      for (int c = 0; c < numvars; ++c)
        fill_as<CI, CJ, CK>(cctkGH, grid, gi, tl, order, fill_negative_x, tt,
                            v0 + c);
    } else {
      fill_as<CI, CJ, CK>(cctkGH, grid, gi, tl, order, fill_negative_x, tt, v0);
    }
    return;
  }

  // Untagged 1-variable GF (HydroBaseX::rho, Z4c::chi, ADM lapse, ...).
  if (numvars == 1) {
    log_group_once(gi, gname, "scalar (untagged)", numvars, verbose);
    fill_as<CI, CJ, CK>(cctkGH, grid, gi, tl, order, fill_negative_x,
                        tensortype_t::scalar, v0);
    return;
  }

  std::vector<CCTK_INT> par;
  const int npar = group_parities(gi, par);
  if (npar == 3 * numvars) {
    log_group_once(gi, gname, "packed (parities walk)", numvars, verbose);
    fill_packed<CI, CJ, CK>(cctkGH, grid, gi, tl, order, fill_negative_x, v0,
                            numvars, par);
    return;
  }

  log_group_once(gi, gname, "SKIP (no tensortypealias/parities)", numvars,
                 verbose);
}

void apply_group(const cGH *cctkGH, const Loop::GridDescBaseDevice &grid,
                 int gi, int tl, int order, bool fill_negative_x, bool verbose) {
  cGroup gdata;
  if (CCTK_GroupData(gi, &gdata) < 0)
    return;
  if (gdata.grouptype != CCTK_GF)
    return;
  // CarpetX does not overload QueryGroupStorageB (the flesh dummy
  // always returns 0). Storage is implied by a non-null data pointer.

  char *gname = CCTK_GroupName(gi);
  if (!gname)
    return;

  int CI = 0, CJ = 0, CK = 0;
  if (!read_centering(gi, CI, CJ, CK)) {
    CCTK_VERROR("Cartoon2DX: group %s has an invalid CENTERING table", gname);
  }

  const int key = CI + 2 * CJ + 4 * CK;
  switch (key) {
  case 0:
    apply_group_c<0, 0, 0>(cctkGH, grid, gi, tl, order, fill_negative_x, verbose,
                           gname, gdata.numvars);
    break;
  case 1:
    apply_group_c<1, 0, 0>(cctkGH, grid, gi, tl, order, fill_negative_x, verbose,
                           gname, gdata.numvars);
    break;
  case 2:
    apply_group_c<0, 1, 0>(cctkGH, grid, gi, tl, order, fill_negative_x, verbose,
                           gname, gdata.numvars);
    break;
  case 3:
    apply_group_c<1, 1, 0>(cctkGH, grid, gi, tl, order, fill_negative_x, verbose,
                           gname, gdata.numvars);
    break;
  case 4:
    apply_group_c<0, 0, 1>(cctkGH, grid, gi, tl, order, fill_negative_x, verbose,
                           gname, gdata.numvars);
    break;
  case 5:
    apply_group_c<1, 0, 1>(cctkGH, grid, gi, tl, order, fill_negative_x, verbose,
                           gname, gdata.numvars);
    break;
  case 6:
    apply_group_c<0, 1, 1>(cctkGH, grid, gi, tl, order, fill_negative_x, verbose,
                           gname, gdata.numvars);
    break;
  case 7:
    apply_group_c<1, 1, 1>(cctkGH, grid, gi, tl, order, fill_negative_x, verbose,
                           gname, gdata.numvars);
    break;
  default:
    CCTK_VERROR("Cartoon2DX: group %s centering {%d,%d,%d} is not 0/1", gname,
                CI, CJ, CK);
  }
  std::free(gname);
}

} // namespace

extern "C" void Cartoon2DX_ApplyBoundary(CCTK_POINTER_TO_CONST cctkGH_,
                                         CCTK_INT gi, CCTK_INT tl) {
  DECLARE_CCTK_PARAMETERS;
  const cGH *const cctkGH = static_cast<const cGH *>(cctkGH_);
  if (tl < 0)
    CCTK_VERROR("Cartoon2DX: ApplyCartoonBoundary got tl=%d for group %s", tl,
                CCTK_FullGroupName(gi));
  const Loop::GridDescBaseDevice grid(cctkGH);
  apply_group(cctkGH, grid, gi, tl, order, fill_negative_x, verbose);
}

} // namespace Cartoon2DX
