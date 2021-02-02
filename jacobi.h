#ifndef JACOBI_H
#define JACOBI_H

#include <openssl/bn.h>

struct Point{
    BIGNUM * x;
    BIGNUM * y;
    BIGNUM * t;
    BIGNUM * z;
};

struct Curve{
    BIGNUM * p;
    BIGNUM * q;
    BIGNUM * b;
};

void init_point(struct Point * p, BIGNUM * X, BIGNUM * Y, BIGNUM * T, BIGNUM * Z);

void init_curve(struct Curve * C);

struct Point add_points(struct Point p_1, struct Point  p_2);

struct Point double_point(struct Point p);

struct Point multiple_point(struct Point p, BIGNUM * k);

int if_contains(struct Point P);

#endif
