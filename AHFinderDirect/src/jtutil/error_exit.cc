// error_exit.cc -- print an error message and exit (= CCTK_VWarn() wrapper)
// $Header$
//
// error_exit -- print an error message and exit (= CCTK_VWarn() wrapper)
//

#include <stdio.h>
#include <stdarg.h>
#include <stdlib.h>
#include <string.h>

#ifdef STANDALONE_TEST
  #include "fake_cctk.h"
#else
  #include "cctk.h"
#endif

#include "config.h"
#include "stdc.h"

//******************************************************************************

//
// This function prints an error message (formatted via vfprintf(3S)),
// then exits.  It does *not* return to its caller.
#ifdef STANDALONE_TEST
// The exit is done via either  exit()  (for ERROR_EXIT)  or  abort()
// (for PANIC_EXIT).
#else
// The exit is done via  CCTK_VWarn()  with the caller-supplied
// "message level" argument; if this is one of
//	ERROR_EXIT
//	PANIC_EXIT
// defined in "stdc.h", then the  CCTK_VWarn()  call itself will terminate
// the Cactus run; if not we do an  abort() .
#endif
//
// Despite not actually returning, this function is declared as returning
// an  int , so it may be easily used in conditional expressions like
//      foo = bar_ok ? baz(bar) : error_exit(...);
//
// Usage:
//	error_exit(exit_code, message, args...)
//
// Arguments:
// msg_level = (in) The "message level" for  CCTK_VWarn() .
// format = (in) vprintf(3S) format string for error message to be printed.
// args... = (in) Any additional arguments are (presumably) used in formatting
//		  the error message string.
//
namespace jtutil
	  {
extern "C"
/*VARARGS*/
int error_exit(int msg_level, const char *format, ...)
{
const int N_buffer = 2000;
char buffer[N_buffer];

va_list ap;
va_start(ap, format);
// FIXME: We should really do something if msg was truncated due to
//        overflowing the buffer.  But what to do?
vsnprintf(buffer, N_buffer, format, ap);
va_end(ap);

// delete trailing '\n' if present, since CCTK_VWarn() doesn't like this
const int len = strlen(buffer);
if ((len > 0) && (buffer[len-1] == '\n'))
   then buffer[len-1] = '\0';

#ifdef STANDALONE_TEST
  fprintf(stderr, "%s\n", buffer);
  if (msg_level == PANIC_EXIT)
     then abort();
     else exit(msg_level);
#else
  CCTK_VWarn(0, __LINE__, __FILE__, CCTK_THORNSTRING, "%s", buffer);
#endif

// if we got here, evidently  msg_level  wasn't drastic enough
abort();							/*NOTREACHED*/
}
	  }	// namespace jtutil::
