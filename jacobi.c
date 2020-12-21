#include "jacobi.h"
#include <stdlib.h>

void init_point(struct Point* P, BIGNUM* X, BIGNUM* Y, BIGNUM* T, BIGNUM* Z) {
    BN_bn2dec(P->x, X);
    BN_bn2dec(P->y, Y);
    BN_bn2dec(P->t, T);
    BN_bn2dec(P->z, Z);
}

void init_curve(struct Curve* C) {
    BN_set_word(C->p, p_s);
    BN_set_word(C->q, q_s);
    BN_set_word(C->b, b_s);
}

struct Point add_points(struct Point p_1, struct Point p_2) {
    struct Point new_p;
    struct Curve C;
    init_curve(&C);
    BIGNUM* help;
    help = BN_new();
    init_point(&new_p, help, help, help, help);

    //x3 = x1*z1*y2*t2 + y1*t1*x2*z2
    BN_mul(new_p.x, p_1.x, p_1.z);
    BN_mul(new_p.x, new_p.x, p_2.y);
    BN_mul(new_p.x, new_p.x, p_2.t);
    BN_mul(help, p_1.y, p_1.t);
    BN_mul(help, help, p_2.x);
    BN_mul(help, help, p_2.z);
    BN_add(new_p.x, new_p.x, help);
    BN_mod(new_p.x, new_p.x, C.p);

    //y3 = y1*z1*y2*z2 - x1*t1*x2*t2
    BN_set_word(help, "0", 10);
    BN_mul(new_p.y, p_1.y, p_1.z);
    BN_mul(new_p.y, new_p.y, p_2.y);
    BN_mul(new_p.y, new_p.y, p_2.z);
    BN_mul(help, p_1.x, p_1.t);
    BN_mul(help, help, p_2.x);
    BN_mul(help, help, p_2.t);
    BN_sub(new_p.y, new_p.y, help);
    BN_mod(new_p.y, new_p.y, C.p);

    //t3 = t1*z1*t2*z2 - b*x1*y1*x2*y2
    BN_set_word(help, "0", 10);
    BN_mul(new_p.t, p_1.t, p_1.z);
    BN_mul(new_p.t, new_p.t, p_2.t);
    BN_mul(new_p.t, new_p.t, p_2.z);
    BN_mul(help, C.b, p_1.x);
    BN_mul(help, help, p_1.y);
    BN_mul(help, help, p_2.x);
    BN_mul(help, help, p_2.y);
    BN_sub(new_p.t, new_p.t, help);
    BN_mod(new_p.t, new_p.t, C.p);

    //z3 = z1^2*y2^2 + x2^2*t1^2
    BN_set_word(help, "0", 10);
    BN_mul(new_p.z, p_1.z, p_1.z);
    BN_mul(new_p.z, new_p.z, p_2.y);
    BN_mul(new_p.z, new_p.z, p_2.y);
    BN_mul(help, p_2.x, p_2.x);
    BN_mul(help, help, p_1.t);
    BN_mul(help, help, p_1.t);
    BN_add(new_p.z, new_p.z, help);
    BN_mod(new_p.z, new_p.z, C.p);

    BN_clear(help);

    return new_p;
}

struct Point double_point(struct Point P) {
    BIGNUM* yz, tz, yt;
    yz = BN_new();
    tz = BN_new();
    yt = BN_new();
    struct Point new_p;
    struct Curve C;
    init_curve(&C);
    init_point(&new_p, yz, yz, yz, yz);

    //x3 = 2*x*y*t*z
    BN_set_word(yz, "2", 10);
    BN_mul(new_p.x, P.x, P.y);
    BN_mul(new_p.x, new_p.x, P.t);
    BN_mul(new_p.x, new_p.x, P.z);
    BN_mul(new_p.x, new_p.x, yz);
    BN_mod(new_p.x, new_p.x, C.p);

    //yz = (y*z)^2
    BN_mul(yz, P.y, P.y);
    BN_mul(yz, yz, P.z);
    BN_mul(yz, yz, P.z);

    //tz = (t*z)^2
    BN_mul(tz, P.t, P.t);
    BN_mul(tz, tz, P.z);
    BN_mul(tz, tz, P.z);

    //yt = (y*t)^2
    BN_mul(yt, P.y, P.y);
    BN_mul(yt, yt, P.t);
    BN_mul(yt, yt, P.t);

    //y3 = (y*z)^2 - (t*z)^2 + (y*t)^2
    BN_sub(new_p.y, yz, tz);
    BN_add(new_p.y, new_p.y, yt);
    BN_mod(new_p.y, new_p.y, C.p);

    //t3 = (t*z)^2 - (y*z)^2 + (y*t)^2
    BN_sub(new_p.t, tz, yz);
    BN_add(new_p.t, new_p.t, yt);
    BN_mod(new_p.t, new_p.t, C.p);

    //z3 = (t*z)^2 + (y*z)^2 - (y*t)^2
    BN_add(new_p.z, tz, yz);
    BN_sub(new_p.z, new_p.z, yt);
    BN_mod(new_p.z, new_p.z, C.p);

    BN_clear(yz);
    BN_clear(tz);
    BN_clear(yt);

    return new_p;
}

struct Point multiple_point(struct Point P, BIGNUM* k) {
    struct Point Q, R;
    struct Curve C;
    BIGNUM* help;
    init_curve(&C);
    BN_set_word(help, "1", 10);
    init_point(&Q, 0, help, help, help);
    init_point(&R, P.x, P.y, P.t, P.z);
    int n = BN_mask_bits(k, 2);
    for (int i = n - 1; i >= 0; i--) {
        if (BN_is_bit_set(k, i) == 0) {     
            R = add_points(R, Q);
            Q = double_point(Q);
        }
        else {
            Q = add_points(Q, R);
            R = double_point(R);
        }
    }

    BN_clear(help);

    return Q;
}

int if_contains(struct Point P) {
    BIGNUM* left1, right, left2, help;
    BN_init(left1);
    BN_init(right);
    BN_init(left2);
    BN_init(help);
    struct Curve C;
    init_curve(&C);

    //left1
    BN_mul(left1, P.x, P.x);
    BN_set_bit(left2, left1);
    BN_mul(help, P.y, P.y);
    BN_add(left1, left1, help);
    BN_mod(left1, left1, C.p);

    //left2
    BN_mul(left2, left2, C.b);
    BN_mul(help, P.t, P.t);
    BN_add(left2, left2, help);
    BN_mod(left2, left2, C.p);

    //right
    BN_mul(right, P.z, P.z);
    BN_mod(right, right, C.p);

    if (BN_cmp(left1, right) == 0 && BN_cmp(left2, right) == 0)
        return 1;
    else
        return 0;

    BN_clear(left1);
    BN_clear(right);
    BN_clear(left2);
    BN_clear(help);
}