#ifndef COEFFICIENTS_H
#define COEFFICIENTS_H

#include <stdio.h>

typedef struct
{
    int base;
    double *coefficients;
} coefficients_t;

double legendre(const int degree, const double x);

int alloc_cof (coefficients_t *cof, int base );

int read_cof(FILE *inf, coefficients_t *cof );

void write_cof(coefficients_t *cof, FILE * ouf );

double value_cof(coefficients_t *cof, double x);

#endif
