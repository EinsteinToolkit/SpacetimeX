/* stdc.h -- JT standard C/C++ header file, cut-down for modern code only */
/* $Header$ */

#ifndef AHFINDERDIRECT__STDC_H
#define AHFINDERDIRECT__STDC_H

/*
 * ***** THIS VERSION IS FOR ANSI/ISO C AND C++ *****
 */

/******************************************************************************/

/*
 * I use this macro in all if statements -- I like the symmetry between
 * the two branches of an if statement.
 *
 *	if (foo)
 *	   then bar;
 *
 *	if (foo)
 *	   then bar;
 *	   else baz;
 */
#define then	/* empty */

/******************************************************************************/

/*
 * misc stuff
 */

#ifdef M_PI		/* usually defined in <math.h> */
  #define PI	M_PI
#endif

#define STRING_EQUAL(s_,t_)	(strcmp(s_,t_) == 0)

/* C library uses abs() for integer absolute value; I prefer iabs() */
#define iabs(x_)	abs(x_)

/******************************************************************************/

/*
 * misc sentinal values
 */

/* n.b. this value is usually 10...0 binary for 2's-complement arithmetic */
#define INT_SENTINAL_VALUE	(~ INT_MAX)	/* from <limits.h> */

#define FLOAT_SENTINAL_VALUE	(- FLT_MAX)	/* from <float.h> */
#define DOUBLE_SENTINAL_VALUE	(- DBL_MAX)	/* from <float.h> */

/******************************************************************************/

#ifdef __cplusplus
namespace jtutil
	  {
#endif

/*
 * Low-level code in this thorn is done with  error_exit() .  In this
 * this thorn it's a wrapper around  CCTK_VWarn(msg_level, ...) .  It's
 * declared to return  int  so it may be easily used in conditional
 * expressions like
 *	foo = bar_ok ? baz(bar) : error_exit(...);
 */
#ifdef __cplusplus
  extern "C"
#endif
int error_exit(int msg_level, const char *format, ...)
  #ifdef __GNUC__
  __attribute__ ((noreturn))
  __attribute__ ((format(printf,2,3)))
  #endif
;

/*
 * error_exit() uses the following  msg_level  codes:
 *	ERROR_EXIT	==> something bad has happened
 *			    (eg inconsistent user input)
 *	PANIC_EXIT	==> diasaster (eg internal consistency check failed)
 *			==> we'd like to force (eg) a core dump or stack trace
 *			    (but at present this isn't implemented; we just
 *			     do a  CCTK_VWarn()  just like for  ERROR_EXIT)
 */
#define ERROR_EXIT	(-1)
#define PANIC_EXIT	(-2)

#ifdef __cplusplus
	  }	/* namespace jtutil */
#endif

/******************************************************************************/

#endif	/* AHFINDERDIRECT__STDC_H */
