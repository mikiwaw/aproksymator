#include "makespl.h"
#include "math.h"
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_cblas.h>

double legendre(int degree, double x) {
	double tmp = (double)degree;
	switch (degree) {
	case 0:
		return 1.0;
	case 1:
		return x;
	case 2:
		return (3 * pow(x, 2) - 1.0) / 2.0;
	case 3:
		return (5 * pow(x, 3) - 3 * x) / 2.0;
	case 4:
		return (35 * pow(x, 4) - 30 * pow(x, 2) + 3) / 8.0;
	case 5:
		return (63 * pow(x, 5) - 70 * pow(x, 3) + 15 * x) / 8.0;
	case 6:
		return (231 * pow(x, 6) - 315 * pow(x, 4) + 105 * pow(x, 2) - 5) / 16.0;
	case 7:
		return (429 * pow(x, 7) - 693 * pow(x, 5) + 315 * pow(x, 3) - 35 * x) / 16.0;
	case 8:
		return (6435 * pow(x, 8) - 12012 * pow(x, 6) + 6930 * pow(x, 4) - 1260 * pow(x, 2) + 35) / 128.0;
	case 9:
		return (6435 * pow(x, 9) - 25740 * pow(x, 7) + 18018 * pow(x, 5) - 4620 * pow(x, 3) + 315 * x) / 128.0;
	case 10:
		return (46189 * pow(x, 10) - 109395 * pow(x, 8) + 90090 * pow(x, 6) - 30030 * pow(x, 4) + 3465 * pow(x, 2) - 63) / 256.0;
	case 11:
		return (88179 * pow(x, 11) - 230945 * pow(x, 9) + 218790 * pow(x, 7) - 90090 * pow(x, 5) + 15015 * pow(x, 3) - 693 * x) / 256.0;
	default:
		return (2.0 * tmp + 1.0) / (tmp + 1.0) * x * legendre(degree - 1, x) - tmp / (tmp + 1) * legendre(degree - 2, x);
	}
}

//void rows_swap(matrix_t *X, int a, int b)
//{
//	int i;
//	double tmp;
//	for (i = 0; i < X->cn; i++)
//	{
//
//	}
//}

//gsl_matrix * inverse = gsl_matrix_alloc(n, n);
//gsl_permutation * perm = gsl_permutation_alloc(n);


void  make_spl(points_t *pts, spline_t *spl) {
//nazwy
	gsl_matrix *A, *A_inverse, *BF;
	gsl_permutation * perm;

	double		*x = pts->x;
	double		*y = pts->y;
	double		a = x[0];
	double		b = x[pts->n - 1];
	int			m = pts->n - 2 > 10 ? 10 : pts->n - 2;
	int			i, j, k;
	char		*mEnv = getenv("APPROX_BASE_SIZE");

	if (mEnv != NULL && atoi(mEnv) > 0)
		m = atoi(mEnv);

	m++;

	A = gsl_matrix_alloc(m, m);

	for (i = 0; i < m; i++)
	{
		for (j = 0; j < m; j++)
		{
			for (k = 0; k < pts->n; k++)
			{
				gsl_matrix_set(A, i, j, legendre(i, x[k]) * legendre(j, x[k]));
			}
		}
	}

	BF = gsl_matrix_alloc(m, 1);

	for (i = 0; i < m; i++)
	{
		for (k = 0; k < pts->n; k++)
			gsl_matrix_set(BF, i, 0, legendre(i, x[k]) * y[k]);
	}

	gsl_linalg_LU_decomp(A, perm, &i);
	gsl_linalg_LU_invert(A, perm, A_inverse);


    gsl_matrix_fprintf(stdout, A, "%g");
    gsl_matrix_fprintf(stdout, BF, "%g");
    gsl_matrix_fprintf(stdout, A_inverse, "%g");
}
