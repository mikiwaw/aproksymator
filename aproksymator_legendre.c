#include "makespl.h"
#include <math.h>
#include <float.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_cblas.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_deriv.h>


#define h 1.0




double legendre(int degree, double x)
{
    double tmp = (double)degree;
    switch (degree)
    {
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
        return (2.0 * (tmp - 1) + 1.0) / tmp * x * legendre(degree - 1, x) - (tmp - 1.0) / tmp * legendre(degree - 2, x);
    }
}

struct legendre_wrap_params
{
    int degree;
};

double legendre_wrap(double x, void *p)
{
    struct legendre_wrap_params * params = (struct legendre_wrap_params *)p;
    return legendre(params->degree, x);
}

double d1legendre(double x, void *p)
{
    struct legendre_wrap_params * params = (struct legendre_wrap_params *)p;
    double result = 0, tmp;

    gsl_function F;

    F.function = &legendre_wrap;
    F.params = params;

    gsl_deriv_forward(&F, x, h, &result, &tmp);

    return result;
}

double d2legendre(double x, void *p)
{
    struct legendre_wrap_params * params = (struct legendre_wrap_params *)p;
    double result = 0, tmp;

    gsl_function F;

    F.function = &d1legendre;
    F.params = params;

    gsl_deriv_forward(&F, x, h, &result, &tmp);

    return result;

}

double d3legendre(double x, void *p)
{
    struct legendre_wrap_params * params = (struct legendre_wrap_params *)p;
    double result = 0, tmp;

    gsl_function F;

    F.function = &d2legendre;
    F.params = params;

    gsl_deriv_forward(&F, x, h, &result, &tmp);

    return result;

}

double fi(const gsl_matrix *a, const int m, const double x)
{
    int i;
    double result = 0.0;

    for (i = 0; i < m; i++)
    {
        result += gsl_matrix_get(a, i, 0) * legendre(i, x);
    }

    return result;
}

double d1fi(const gsl_matrix *a, const int m, const double x)
{
    struct legendre_wrap_params params;
    int i;
    double result = 0.0;

    for (i = 0; i < m; i++)
    {
        params.degree = i;
        result += gsl_matrix_get(a, i, 0) * d1legendre(x, &params);
    }

    return result;
}

double d2fi(const gsl_matrix *a, const int m, const double x)
{
    struct legendre_wrap_params params;
    int i;
    double result = 0.0;

    for (i = 0; i < m; i++)
    {
        params.degree = i;
        result += gsl_matrix_get(a, i, 0) * d2legendre(x, &params);
    }

    return result;
}

double d3fi(const gsl_matrix *a, const int m, const double x)
{
    struct legendre_wrap_params params;
    int i;
    double result = 0.0;

    for (i = 0; i < m; i++)
    {
        params.degree = i;
        result += gsl_matrix_get(a, i, 0) * d3legendre(x, &params);
    }

    return result;
}


void  make_spl(points_t *pts, spline_t *spl)
{

    gsl_matrix *A, *A_inverse, *BF, *wsp;
    gsl_permutation * perm;

    double		*x = pts->x;
    double		*y = pts->y;
    double		a = x[0];
    double		b = x[pts->n - 1];
    int			m = pts->n - 1 > 10 ? 10 : pts->n - 1;


    int			i, j, k;
    char		*mEnv = getenv("APPROX_BASE_SIZE");

    if (mEnv != NULL && atoi(mEnv) > 0)
        m = atoi(mEnv);


    m++;

    A = gsl_matrix_alloc(m, m);
    A_inverse = gsl_matrix_alloc(m, m);
    perm = gsl_permutation_alloc(m);
    BF = gsl_matrix_alloc(m, 1);
    wsp = gsl_matrix_alloc(m, 1);

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
//macierz odwrotna
    gsl_linalg_LU_decomp(A, perm, &i);
    gsl_linalg_LU_invert(A, perm, A_inverse);
    gsl_matrix_free(A);
    gsl_permutation_free(perm);

    gsl_matrix_set_zero (BF);
    for (i = 0; i < m; i++)
    {
        for (k = 0; k < pts->n; k++)
            BF->data[i * BF->tda + 0] += legendre(i, x[k]) * y[k];
    }
//mnozenie macierzy
    gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, A_inverse, BF, 0.0, wsp);
    gsl_matrix_free(A_inverse);
    gsl_matrix_free(BF);

    if (alloc_spl(spl, m) == 0)
    {
        for (i = 0; i < spl->n; i++)
        {
            double xx = spl->x[i] = a + i*(b - a) / (spl->n - 1);
            xx += 10.0*DBL_EPSILON;
            spl->f[i] = fi(wsp, m, xx);
            spl->f1[i] = d1fi(wsp, m, xx);
            spl->f2[i] = d2fi(wsp, m, xx);
            spl->f3[i] = d3fi(wsp, m, xx);
        }
    }

    return;
}
