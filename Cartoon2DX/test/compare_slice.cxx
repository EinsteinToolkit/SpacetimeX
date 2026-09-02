// compare_slice.cxx -- pointwise comparison tool for the Cartoon2DX
// correctness test (requirements.md). Extracts a y=const 2D (x,z) slice
// of a named variable from two CarpetX Silo output files (one full-3D,
// one thin-y Cartoon2DX) and reports the pointwise difference at
// matching (x,z) grid points. Not part of any thorn -- a one-off
// analysis tool, kept here rather than in a thorn's src/ since it has
// nothing to do with the simulation itself.
//
// Build (single line):
//   g++ -std=c++17 -O2 -I/home/sbrandt/Cactus/configs/sim/scratch/external/Silo/include -I/usr/include/hdf5/openmpi compare_slice.cxx -o compare_slice -L/home/sbrandt/Cactus/configs/sim/scratch/external/Silo/lib -lsiloh5 -L/usr/lib/x86_64-linux-gnu/hdf5/openmpi -lhdf5 -lz -lm -Wl,-rpath,/usr/lib/x86_64-linux-gnu/hdf5/openmpi
//
// Usage:
//   compare_slice <file1.silo> <file2.silo> <target_y> <exclude_radius> <varname> [varname...]
//
// exclude_radius: points within this distance (in the x-z plane) of
// either puncture (hardcoded below at x=0,z=+-1.15, matching
// BrillLindquistMulti's defaults) are excluded from the reported
// statistics -- pass 0 to include everything.
//
// For each variable, opens the multimesh/multivar in each file, scans
// all quadmesh blocks, and for any block whose y-coordinate array
// contains target_y (within tolerance), extracts the (x,z) slice at
// that y-index -- restricted to each block's own interior (min_index/
// max_index) to avoid double-counting the ghost-zone overlap between
// neighboring blocks. Points are matched between the two files by exact
// (x,z) coordinate (both runs share dx=dz and domain extent, so this
// should be an exact match, not an interpolation).

#include <silo.h>

#include <unistd.h>

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <map>
#include <string>
#include <vector>

namespace {

struct Key {
  long xi, zi; // coordinates rounded to a fixed grid to use as exact map keys
  bool operator<(const Key &o) const {
    return xi < o.xi || (xi == o.xi && zi < o.zi);
  }
};

Key make_key(double x, double z) {
  // Round to the nearest 1e-6 -- both runs share dx=dz to full double
  // precision (same formula, same xmin/xmax/ncells), so this is just
  // insurance against last-bit roundoff, not real quantization.
  return Key{std::lround(x * 1e6), std::lround(z * 1e6)};
}

// Multimesh/multivar block references look like "subdir/file.silo:objname"
// -- DBGetQuadmesh/DBGetQuadvar on the TOP-LEVEL DBfile handle do not
// auto-follow the "file:" prefix (that's only done by higher-level tools
// like VisIt), so each referenced sub-file must be opened directly and
// queried by its bare object name. Caches open sub-files by name since
// blocks commonly share one (e.g. all blocks in a single-process run).
struct SubfileCache {
  std::map<std::string, DBfile *> open_files;

  // Splits "file:obj" on the LAST colon, opens `file` (cached), and
  // returns (that DBfile*, "obj"). Exits the process on any failure.
  std::pair<DBfile *, std::string> resolve(const char *ref) {
    const std::string s(ref);
    const size_t colon = s.find_last_of(':');
    if (colon == std::string::npos) {
      std::fprintf(stderr, "ERROR: reference has no ':' separator: %s\n", ref);
      std::exit(1);
    }
    const std::string file = s.substr(0, colon);
    const std::string obj = s.substr(colon + 1);
    auto it = open_files.find(file);
    if (it != open_files.end()) return {it->second, obj};
    DBfile *db = DBOpen(file.c_str(), DB_UNKNOWN, DB_READ);
    if (!db) {
      std::fprintf(stderr, "ERROR: could not open referenced sub-file %s\n", file.c_str());
      std::exit(1);
    }
    open_files[file] = db;
    return {db, obj};
  }

  ~SubfileCache() {
    for (auto &kv : open_files) DBClose(kv.second);
  }
};

// Extracts a y=target_y slice of `varname` from `filename` into `out`
// (keyed by rounded (x,z)). Returns the number of points extracted.
int extract_slice(const char *filename, double target_y, const char *varname,
                  std::map<Key, double> &out) {
  // The multimesh's block references are relative paths ("subdir/file:obj"),
  // which Silo resolves relative to the CURRENT WORKING DIRECTORY, not the
  // top-level file's own location -- so chdir into the file's directory
  // and open it by basename, restoring cwd when done.
  char saved_cwd[4096];
  if (!getcwd(saved_cwd, sizeof saved_cwd)) { std::perror("getcwd"); std::exit(1); }
  std::string path(filename);
  const size_t slash = path.find_last_of('/');
  std::string dir = (slash == std::string::npos) ? "." : path.substr(0, slash);
  std::string base = (slash == std::string::npos) ? path : path.substr(slash + 1);
  if (chdir(dir.c_str()) != 0) { std::perror("chdir"); std::exit(1); }

  DBfile *db = DBOpen(base.c_str(), DB_UNKNOWN, DB_READ);
  if (!db) {
    std::fprintf(stderr, "ERROR: could not open %s\n", filename);
    std::exit(1);
  }

  char **toc_names = nullptr;
  DBtoc *toc = DBGetToc(db);
  if (!toc || toc->nmultimesh < 1) {
    std::fprintf(stderr, "ERROR: no multimesh found in %s\n", filename);
    std::exit(1);
  }
  const char *meshname = toc->multimesh_names[0];

  DBmultimesh *mm = DBGetMultimesh(db, meshname);
  if (!mm) {
    std::fprintf(stderr, "ERROR: DBGetMultimesh(%s) failed in %s\n", meshname, filename);
    std::exit(1);
  }
  DBmultivar *mv = DBGetMultivar(db, varname);
  if (!mv) {
    std::fprintf(stderr, "ERROR: DBGetMultivar(%s) failed in %s\n", varname, filename);
    std::exit(1);
  }
  if (mv->nvars != mm->nblocks) {
    std::fprintf(stderr, "ERROR: block count mismatch (mesh=%d, var=%d) in %s\n",
                mm->nblocks, mv->nvars, filename);
    std::exit(1);
  }

  SubfileCache cache;
  int n_extracted = 0;
  for (int b = 0; b < mm->nblocks; ++b) {
    auto [mesh_db, mesh_obj] = cache.resolve(mm->meshnames[b]);
    DBquadmesh *qm = DBGetQuadmesh(mesh_db, mesh_obj.c_str());
    if (!qm) continue;
    if (qm->coordtype != DB_COLLINEAR || qm->datatype != DB_DOUBLE) {
      std::fprintf(stderr, "ERROR: unexpected mesh format in block %d of %s\n", b, filename);
      std::exit(1);
    }
    const double *xs = static_cast<const double *>(qm->coords[0]);
    const double *ys = static_cast<const double *>(qm->coords[1]);
    const double *zs = static_cast<const double *>(qm->coords[2]);

    // Find the y-index matching target_y, if this block's y-range covers it.
    int jy = -1;
    for (int j = 0; j < qm->dims[1]; ++j) {
      if (std::fabs(ys[j] - target_y) < 1e-9) { jy = j; break; }
    }
    if (jy < 0) { DBFreeQuadmesh(qm); continue; }

    auto [var_db, var_obj] = cache.resolve(mv->varnames[b]);
    DBquadvar *qv = DBGetQuadvar(var_db, var_obj.c_str());
    if (!qv) {
      std::fprintf(stderr, "ERROR: DBGetQuadvar(%s) failed in %s\n", mv->varnames[b], filename);
      std::exit(1);
    }
    const double *vals = static_cast<const double *>(qv->vals[0]);

    for (int i = qm->min_index[0]; i <= qm->max_index[0]; ++i) {
      for (int k = qm->min_index[2]; k <= qm->max_index[2]; ++k) {
        const long idx = i * (long)qv->stride[0] + jy * (long)qv->stride[1] +
                         k * (long)qv->stride[2];
        out[make_key(xs[i], zs[k])] = vals[idx];
        ++n_extracted;
      }
    }
    DBFreeQuadvar(qv);
    DBFreeQuadmesh(qm);
  }

  DBFreeMultivar(mv);
  DBFreeMultimesh(mm);
  DBClose(db);
  if (chdir(saved_cwd) != 0) { std::perror("chdir back"); std::exit(1); }
  return n_extracted;
}

} // namespace

int main(int argc, char **argv) {
  if (argc < 5) {
    std::fprintf(stderr,
                "Usage: %s <file1.silo> <file2.silo> <target_y> <varname> [varname...]\n",
                argv[0]);
    return 1;
  }
  const char *file1 = argv[1];
  const char *file2 = argv[2];
  const double target_y = std::atof(argv[3]);
  const double exclude_radius = std::atof(argv[4]);
  // Puncture positions, matching BrillLindquistMulti's defaults
  // (z0_1/z0_2; see repos/SpacetimeX/BrillLindquistMulti/param.ccl).
  const double punctures_xz[2][2] = {{0.0, 1.15}, {0.0, -1.15}};

  for (int a = 5; a < argc; ++a) {
    const char *varname = argv[a];
    std::map<Key, double> s1, s2;
    const int n1 = extract_slice(file1, target_y, varname, s1);
    const int n2 = extract_slice(file2, target_y, varname, s2);

    int n_matched = 0, n_excluded = 0;
    double max_abs_diff = 0, rms_diff_sq_sum = 0;
    Key max_key{0, 0};
    double max_v1 = 0, max_v2 = 0;
    for (const auto &kv : s1) {
      auto it = s2.find(kv.first);
      if (it == s2.end()) continue;
      const double x = kv.first.xi / 1e6, z = kv.first.zi / 1e6;
      if (exclude_radius > 0) {
        bool near_puncture = false;
        for (const auto &p : punctures_xz) {
          const double dx = x - p[0], dz = z - p[1];
          if (std::sqrt(dx * dx + dz * dz) < exclude_radius) { near_puncture = true; break; }
        }
        if (near_puncture) { ++n_excluded; continue; }
      }
      const double diff = std::fabs(kv.second - it->second);
      rms_diff_sq_sum += diff * diff;
      if (diff > max_abs_diff) {
        max_abs_diff = diff;
        max_key = kv.first;
        max_v1 = kv.second;
        max_v2 = it->second;
      }
      ++n_matched;
    }
    const double rms = n_matched > 0 ? std::sqrt(rms_diff_sq_sum / n_matched) : 0.0;

    std::printf("%-28s file1_pts=%-6d file2_pts=%-6d matched=%-6d excluded=%-6d "
               "max_abs_diff=%.6e (at x=%.6g,z=%.6g: v1=%.6e v2=%.6e) rms_diff=%.6e\n",
               varname, n1, n2, n_matched, n_excluded, max_abs_diff, max_key.xi / 1e6,
               max_key.zi / 1e6, max_v1, max_v2, rms);

    // Explicit symmetric boundary-point check: is the discrepancy at one
    // z-edge mirrored at the other, or one-sided? Prints both faces along
    // x=0 (the z-axis) plus, for context, one grid point inward from each.
    for (double z : {-8.0, -7.875, 0.0, 0.125, 1.15, 7.875, 8.0}) {
      const Key k = make_key(0.0, z);
      auto i1 = s1.find(k), i2 = s2.find(k);
      if (i1 == s1.end() || i2 == s2.end()) {
        std::printf("  %-26s x=0,z=%-8g  (point not found in one or both files)\n", "", z);
        continue;
      }
      std::printf("  %-26s x=0,z=%-8g  v1=%.6e v2=%.6e diff=%.6e\n", "", z, i1->second,
                 i2->second, std::fabs(i1->second - i2->second));
    }
    // Near-origin x-scan at z=0.125 (the slice that previously went flat
    // to gxx~1 in cartoon while full3d developed a bump).
    for (double x : {-0.5, -0.25, -0.125, 0.0, 0.125, 0.25, 0.5}) {
      const Key k = make_key(x, 0.125);
      auto i1 = s1.find(k), i2 = s2.find(k);
      if (i1 == s1.end() || i2 == s2.end()) {
        std::printf("  %-26s x=%-8g,z=0.125  (point not found in one or both files)\n", "", x);
        continue;
      }
      std::printf("  %-26s x=%-8g,z=0.125  v1=%.6e v2=%.6e diff=%.6e\n", "", x, i1->second,
                 i2->second, std::fabs(i1->second - i2->second));
    }
  }
  return 0;
}
