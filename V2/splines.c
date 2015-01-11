#include "splines.h"

#include <stdlib.h>

#define MALLOC_FAILED( P, SIZE ) (((P)=malloc( (SIZE)*sizeof( *(P))))==NULL)

int
alloc_spl (spline_t * spl, int n)
{
    spl->n = n;
    return MALLOC_FAILED (spl->x, spl->n);
}

int
read_spl (FILE * inf, spline_t * spl)
{
    int i;
    if (fscanf (inf, "%d", &(spl->n)) != 1 || spl->n < 0)
        return 1;

    if (alloc_spl (spl, spl->n))
        return 1;

    for (i = 0; i < spl->n; i++)
        if (fscanf(inf, "%lf", spl->x + i) != 1)
            return 1;

    return 0;
}

void
write_spl (spline_t * spl, FILE * ouf)
{
    int i;
    fprintf (ouf, "%d\n", spl->n);
    for (i = 0; i < spl->n; i++)
        fprintf (ouf, "%g\n", spl->x[i]);
}

double
value_spl (spline_t * spl, double x)
{
    return fi_wrap(spl, x);
}
