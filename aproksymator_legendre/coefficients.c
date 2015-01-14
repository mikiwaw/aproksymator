#include "coefficients.h"

#include <math.h>
#include <stdlib.h>

#define MALLOC_FAILED( P, SIZE ) (((P)=malloc( (SIZE)*sizeof( *(P))))==NULL)

double legendre(const int degree, const double x)
{
    long tmp = (double)degree;
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
        return (12155 * pow(x, 9) - 25740 * pow(x, 7) + 18018 * pow(x, 5) - 4620 * pow(x, 3) + 315 * x) / 128.0;
    case 10:
        return (46189 * pow(x, 10) - 109395 * pow(x, 8) + 90090 * pow(x, 6) - 30030 * pow(x, 4) + 3465 * pow(x, 2) - 63) / 256.0;
    case 11:
        return (88179 * pow(x, 11) - 230945 * pow(x, 9) + 218790 * pow(x, 7) - 90090 * pow(x, 5) + 15015 * pow(x, 3) - 693 * x) / 256.0;
    default:
        return (2.0 * (tmp - 1) + 1.0) / tmp * x * legendre(degree - 1, x) - (tmp - 1.0) / tmp * legendre(degree - 2, x);
    }
}

int
alloc_cof (coefficients_t * cof, int base)
{
    cof->base = base;
    return MALLOC_FAILED (cof->coefficients, cof->base);
}

int
read_cof (FILE * inf, coefficients_t * cof)
{
    int i;
    if (fscanf (inf, "%d", &(cof->base)) != 1 || cof->base < 0)
        return 1;

    if (alloc_cof (cof, cof->base))
        return 1;

    for (i = 0; i < cof->base; i++)
        if (fscanf (inf, "%lf", cof->coefficients + i) != 1)
            return 1;

    return 0;
}

void
write_cof (coefficients_t * cof, FILE * ouf)
{
    int i;
    fprintf (ouf, "%d\n", cof->base);
    for (i = 0; i < cof->base; i++)
        fprintf (ouf, "%g\n", cof->coefficients[i]);
}

double
value_cof (coefficients_t * cof, double x)
{
    int i;
    double result = 0.0;



    for (i = 0; i < cof->base; i++)
    {
        result += cof->coefficients[i] * legendre(i, x);
    }

    return result;
}
