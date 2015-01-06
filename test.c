//#include "aproksymator_legendre.c"
#include "gaus/piv_ge_solver.h"
#include "gaus/piv_ge_solver.c"
#include "gaus/pivot.c"
#include "gaus/matrix.h"
#include "gaus/matrix.c"
#include <stdio.h>
#include "gsl/matrix/gsl_matrix.h"

int main()
{
	matrix_t       *eqs = NULL;
	int x, y, v;

	printf("podaj wiersze:\n");
	scanf("%d", &x);
	printf("podaj kolumny:\n");
	scanf("%d", &y);

	eqs = make_matrix(x, y);
	for (x = 0; x < eqs->rn; x++)
		for (y = 0; y < eqs->cn; y++)
		{
			printf("%dx%d :", x, y);
			scanf("%d", &v);
			put_entry_matrix(eqs, x, y, v);

		}


	//put_entry_matrix(eqs, 0, 0, 1.0);
	//put_entry_matrix(eqs, 0, 1, 1.0);
	//put_entry_matrix(eqs, 0, 2, 1.0);
	////put_entry_matrix(eqs, 0, 3, 1.0);
	//put_entry_matrix(eqs, 1, 0, 2.0);
	//put_entry_matrix(eqs, 1, 1, 1.0);
	//put_entry_matrix(eqs, 1, 2, 1.0);
	////put_entry_matrix(eqs, 1, 3, 5.0);
	//put_entry_matrix(eqs, 2, 0, 1.0);
	//put_entry_matrix(eqs, 2, 1, 2.0);
	//put_entry_matrix(eqs, 2, 2, 5.0);
	////put_entry_matrix(eqs, 2, 3, 11.0);

	write_matrix(eqs, stdout);
	piv_ge_solver(eqs);
	write_matrix(eqs, stdout);

	system("pause");
	return 0;
}
