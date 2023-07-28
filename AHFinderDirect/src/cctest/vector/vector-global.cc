// $Header$

#include <stdio.h>
#include <vector>

int main()
{
printf("testing <vector> functions in global namespace (THIS SHOULD FAIL):\n");
vector<int> v(3);
v[0] = 42;
v[1] = 69;
v[2] = 105;
printf("%d %d %d should be 42 69 105... ", v[0], v[1], v[2]);
printf(((v[0] == 42) && (v[1] == 69) && (v[2] == 105)) ? "ok\n" : "FAIL\n");
printf("==> #include <vector>; vector is ok (THIS SHOULD FAIL)\n");
return 0;
}
