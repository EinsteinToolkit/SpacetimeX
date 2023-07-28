// inline-class.cc -- test simple class with inline function
// $Header$
#include <stdio.h>

class	addx
	{
public:
	int operator()(int i) { return i+x_; }
	addx(int x);
private:
	int x_;
	};

addx::addx(int x)
	: x_(x)
{ }

int main(void)
{
addx add3(3);
printf("add3(42) = %d\n", add3(42));
return 0;
}
