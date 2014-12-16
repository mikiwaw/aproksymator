#include "aproksymator_legendre.c"
#include <stdio.h>

int main()
{
	int i, j;
	double x;
	i = 2;
	x = 1;
	
	printf("%d\t%f\t%f\n", i, x, legendre(i, x));
	
	for (i = 0; i <= 10; i++)
		for (j = 0; j <= 10; j++)
			printf("%d\t%d\t%f\n", i, j, legendre(i, j));

	system("pause");
	return 0;
}
