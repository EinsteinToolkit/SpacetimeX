// inline-fn.cc -- test non-class inline function
// $Header$
#include <stdio.h>

inline int add3(int x) { return x+3; }

int main(void)
{
printf("add3(42) = %d\n", add3(42));
return 0;
}
