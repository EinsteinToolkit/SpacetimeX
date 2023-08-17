// stl_vector.hh -- define STL vector class in global namespace.
// $Header$

//
// prerequisites:
//	"cctk.h"
//

//
// This header file defines the STL vector class in the global namespace.
// It makes use of Cactus's configuration-time probing of what header
// files are available.
//

#if   defined(HAVE_VECTOR)
  #include <vector>
  using std::vector;
#elif defined(HAVE_VECTOR_H)
  #include <vector.h>
#else
  #error "Cactus couldn't find the STL vector class!"
#endif
