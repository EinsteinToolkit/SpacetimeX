// $Header$

#include <cstdio>

int main()
{
printf("testing <cstdio> functions in global namespace (THIS SHOULD FAIL):\n");
printf("==> #include <cstdio>; printf() is ok (THIS SHOULD FAIL)\n");
return 0;
}
