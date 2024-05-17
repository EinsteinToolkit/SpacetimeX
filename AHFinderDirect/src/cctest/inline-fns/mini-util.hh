// mini-util.hh -- cut-down version of src/jtutil/util.hh for testing inline fns
// $Header$

namespace jtutil
          {
// how many integers are in the closed interval [low,high]
inline int how_many_in_range(int low, int high) { return high - low + 1; }
          }
