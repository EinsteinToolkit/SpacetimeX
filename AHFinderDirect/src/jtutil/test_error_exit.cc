// test_error_exit.cc -- test driver for  error_exit()  function
// $Header$

#include <stdio.h>

#include "stdc.h"
using jtutil::error_exit;

int main()
{
error_exit(ERROR_EXIT, "two+two=%.3f", 4.0);			/*NOTREACHED*/
return 0;
}
