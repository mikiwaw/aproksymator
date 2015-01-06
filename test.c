//#include "aproksymator_legendre.c"
#include "gaus/piv_ge_solver.h"
#include "gaus/piv_ge_solver.c"
#include "gaus/pivot.c"
#include "gaus/matrix.h"
#include "gaus/matrix.c"
#include <stdio.h>
#include "gsl/matrix/gsl_matrix.h"
#include "gsl/linalg/gsl_linalg.h"
#include "gsl/cblas/gsl_cblas.h"



int main()
{

	int n = 2;
	int s;
	int i, j, v;
	matrix_t *X;

	gsl_matrix * m = gsl_matrix_alloc(n, n);
	gsl_matrix * inverse = gsl_matrix_alloc(n, n);
	gsl_permutation * perm = gsl_permutation_alloc(n);

	X = make_martix(n, n);
	for (i = 0; i < X->rn; i++)
		for (j = 0; j < X->cn; j++)
		{
			printf("%dx%d: ", i, j);
			scnaf("%d", &v);
			put_entry_matrix(X, i, j, v);
			gsl_matrix_set(m, i, j, v);
		}



	gsl_matrix_fprintf(stdout, m, "%g");
	gsl_matrix_fprintf(stdout, m, "%f");

	gsl_linalg_LU_decomp(m, perm, &s);


	gsl_linalg_LU_invert(m, perm, inverse);

	for (i = 0; i < X->rn; i++)
		for (j = 0; j < X->cn; j++)
		{
			put_entry_matrix(X, i, j, gsl_matrix_get(m, i, j));
		}


	//matrix_t       *eqs = NULL;
	//int x, y, v;

	//printf("podaj wiersze:\n");
	//scanf("%d", &x);
	//printf("podaj kolumny:\n");
	//scanf("%d", &y);

	//eqs = make_matrix(x, y);
	//for (x = 0; x < eqs->rn; x++)
	//	for (y = 0; y < eqs->cn; y++)
	//	{
	//		printf("%dx%d :", x, y);
	//		scanf("%d", &v);
	//		put_entry_matrix(eqs, x, y, v);

	//	}


	////put_entry_matrix(eqs, 0, 0, 1.0);
	////put_entry_matrix(eqs, 0, 1, 1.0);
	////put_entry_matrix(eqs, 0, 2, 1.0);
	//////put_entry_matrix(eqs, 0, 3, 1.0);
	////put_entry_matrix(eqs, 1, 0, 2.0);
	////put_entry_matrix(eqs, 1, 1, 1.0);
	////put_entry_matrix(eqs, 1, 2, 1.0);
	//////put_entry_matrix(eqs, 1, 3, 5.0);
	////put_entry_matrix(eqs, 2, 0, 1.0);
	////put_entry_matrix(eqs, 2, 1, 2.0);
	////put_entry_matrix(eqs, 2, 2, 5.0);
	//////put_entry_matrix(eqs, 2, 3, 11.0);

	//write_matrix(eqs, stdout);
	//piv_ge_solver(eqs);
	//write_matrix(eqs, stdout);

	//system("pause");
	return 0;
}
