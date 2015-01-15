#include "makecof.h"

#include <math.h>
//#include <float.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>
//#include <gsl/gsl_cblas.h>
//#include <gsl/gsl_blas.h>

void  make_cof(points_t *pts, coefficients_t *cof)
{

    gsl_matrix *A;
    gsl_permutation *perm;
    gsl_vector *wsp, *BF;

    double		*x = pts->x;
    double		*y = pts->y;
    double		a = x[0];
    double		b = x[pts->n - 1];
    int			m = pts->n - 2 > 10 ? 10 : pts->n - 2;
    int			i, j, k;
    double      tmp;
    char		*mEnv = getenv("APPROX_BASE_SIZE");

    if (mEnv != NULL && atoi(mEnv) > 0)
        m = atoi(mEnv);
    m++;

    A = gsl_matrix_alloc(m, m);
    wsp = gsl_vector_alloc(m);
    BF = gsl_vector_alloc(m);
    perm = gsl_permutation_alloc(m);

    gsl_matrix_set_zero (A);
    for (i = 0; i < m; i++)
    {
        for (j = 0; j < m; j++)
        {
            for (k = 0; k < pts->n; k++)
            {
                A->data[i * A->tda + j] += legendre(i, x[k]) * legendre(j, x[k]);
            }
        }
    }

    gsl_vector_set_zero (BF);
    for (i = 0; i < m; i++)
    {
        tmp = 0.0;
        for (k = 0; k < pts->n; k++)
            tmp += legendre(i, x[k]) * y[k];
        gsl_vector_set(BF, i, tmp);
    }

    gsl_linalg_LU_decomp(A, perm, &i);
    gsl_linalg_LU_solve(A, perm, BF, wsp);

    gsl_matrix_free(A);
    gsl_permutation_free(perm);
    gsl_vector_free(BF);

    if (alloc_cof(cof, m) == 0)
    {
        for (i = 0; i < cof->base; i++)
            cof->coefficients[i] = gsl_vector_get(wsp, i);
    }

    gsl_vector_free(wsp);

    return;
}
