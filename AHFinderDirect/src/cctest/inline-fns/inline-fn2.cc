// inline-fn2.cc -- test non-class inline function in header file
// $Header$
#include <stdio.h>
#include "mini-util.hh"

int main(void)
{
printf("jtutil::how_many_in_range(69,105) = %d\n",
       jtutil::how_many_in_range(69,105));
return 0;
}
