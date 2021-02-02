#include "jacobi.h"
#include <stdlib.h>

const char * p_s = "115792089237316195423570985008687907853269984665640564039457584007913111864739";
const char * q_s = "4824670384888174809315457708695329493868526062531865828358260708486391273721";
const char * b_s = "115792089237316195423570985008687907853269984665640564039457584007913111864737";

void init_point(struct Point * P, BIGNUM * X, BIGNUM * Y, BIGNUM * T, BIGNUM * Z){
    P->x = BN_new();
    P->y = BN_new();
    P->t = BN_new();
    P->z = BN_new();
    BN_copy(P->x, X);
    BN_copy(P->y, Y);
    BN_copy(P->t, T);
    BN_copy(P->z, Z);
}

void init_curve(struct Curve * C){
    C->q = BN_new();
    C->p = BN_new();
    C->b = BN_new();
    BN_dec2bn(&C->p, p_s);
    BN_dec2bn(&C->q, q_s);
    BN_dec2bn(&C->b, b_s);
}

struct Point add_points(struct Point p_1, struct Point p_2){
    struct Point new_p;
    struct Curve C;

    C.q = BN_new();
    C.p = BN_new();
    C.b = BN_new();

    new_p.t = BN_new();
    new_p.x = BN_new();
    new_p.y = BN_new();
    new_p.z = BN_new();

    init_curve(&C);
    BIGNUM * help = BN_new();
    init_point(&new_p, help, help, help, help);

    BN_CTX *ctx = BN_CTX_new();
    //x3 = x1*z1*y2*t2 + y1*t1*x2*z2
    BN_mod_mul(new_p.x, p_1.x, p_1.z, C.p, ctx);
    BN_mod_mul(new_p.x, new_p.x, p_2.y, C.p, ctx);
    BN_mod_mul(new_p.x, new_p.x, p_2.t, C.p, ctx);
    BN_mod_mul(help, p_1.y, p_1.t, C.p, ctx);
    BN_mod_mul(help, help, p_2.x, C.p, ctx);
    BN_mod_mul(help, help, p_2.z, C.p, ctx);
    BN_mod_add(new_p.x, new_p.x, help, C.p, ctx);
    //BN_mod(new_p.x, new_p.x, C.p, ctx);

    //y3 = y1*z1*y2*z2 - x1*t1*x2*t2
    BN_dec2bn(&help, "0");
    BN_mod_mul(new_p.y, p_1.y, p_1.z, C.p, ctx);
    BN_mod_mul(new_p.y, new_p.y, p_2.y, C.p, ctx);
    BN_mod_mul(new_p.y, new_p.y, p_2.z, C.p, ctx);
    BN_mod_mul(help, p_1.x, p_1.t, C.p, ctx);
    BN_mod_mul(help, help, p_2.x, C.p, ctx);
    BN_mod_mul(help, help, p_2.t, C.p, ctx);
    BN_mod_sub(new_p.y, new_p.y, help, C.p, ctx);
    //BN_mod(new_p.y, new_p.y, C.p, ctx);

    //t3 = t1*z1*t2*z2 - b*x1*y1*x2*y2
    BN_dec2bn(&help, "0");
    BN_mod_mul(new_p.t, p_1.t, p_1.z, C.p, ctx);
    BN_mod_mul(new_p.t, new_p.t, p_2.t, C.p, ctx);
    BN_mod_mul(new_p.t, new_p.t, p_2.z, C.p, ctx);
    BN_mod_mul(help, C.b, p_1.x, C.p, ctx);
    BN_mod_mul(help, help, p_1.y, C.p, ctx);
    BN_mod_mul(help, help, p_2.x, C.p, ctx);
    BN_mod_mul(help, help, p_2.y, C.p, ctx);
    BN_mod_sub(new_p.t, new_p.t, help, C.p, ctx);
    //BN_mod(new_p.t, new_p.t, C.p, ctx);

    //z3 = z1^2*y2^2 + x2^2*t1^2
    BN_dec2bn(&help, "0");
    BN_mod_mul(new_p.z, p_1.z, p_1.z, C.p, ctx);
    BN_mod_mul(new_p.z, new_p.z, p_2.y, C.p, ctx);
    BN_mod_mul(new_p.z, new_p.z, p_2.y, C.p, ctx);
    BN_mod_mul(help, p_2.x, p_2.x, C.p, ctx);
    BN_mod_mul(help, help, p_1.t, C.p, ctx);
    BN_mod_mul(help, help, p_1.t, C.p, ctx);
    BN_mod_add(new_p.z, new_p.z, help, C.p, ctx);
    //BN_mod(new_p.z, new_p.z, C.p, ctx);
    BN_CTX_free(ctx);

    BN_free(help);

    return new_p;
}

struct Point double_point(struct Point P){
    BIGNUM * yz;
    BIGNUM * tz;
    BIGNUM * yt;
    yz = BN_new();
    tz = BN_new();
    yt = BN_new();
    struct Point new_p;
    struct Curve C;
    init_curve(&C);
    init_point(&new_p, yz, yz, yz, yz);

    //x3 = 2*x*y*t*z
    BN_CTX *ctx = BN_CTX_new();
    BN_dec2bn(&yz, "2");
    BN_mod_mul(new_p.x, P.x, P.y, C.p, ctx);
    BN_mod_mul(new_p.x, new_p.x, P.t, C.p, ctx);
    BN_mod_mul(new_p.x, new_p.x, P.z, C.p, ctx);
    BN_mod_mul(new_p.x, new_p.x, yz, C.p, ctx);

    //yz = (y*z)^2
    BN_mod_mul(yz, P.y, P.y, C.p, ctx);
    BN_mod_mul(yz, yz, P.z, C.p, ctx);
    BN_mod_mul(yz, yz, P.z, C.p, ctx);

    //tz = (t*z)^2
    BN_mod_mul(tz, P.t, P.t, C.p, ctx);
    BN_mod_mul(tz, tz, P.z, C.p, ctx);
    BN_mod_mul(tz, tz, P.z, C.p, ctx);

    //yt = (y*t)^2
    BN_mod_mul(yt, P.y, P.y, C.p, ctx);
    BN_mod_mul(yt, yt, P.t, C.p, ctx);
    BN_mod_mul(yt, yt, P.t, C.p, ctx);

    //y3 = (y*z)^2 - (t*z)^2 + (y*t)^2
    BN_mod_sub(new_p.y, yz, tz, C.p, ctx);
    BN_mod_add(new_p.y, new_p.y, yt, C.p, ctx);

    //t3 = (t*z)^2 - (y*z)^2 + (y*t)^2
    BN_mod_sub(new_p.t, tz, yz, C.p, ctx);
    BN_mod_add(new_p.t, new_p.t, yt, C.p, ctx);

    //z3 = (t*z)^2 + (y*z)^2 - (y*t)^2
    BN_mod_add(new_p.z, tz, yz, C.p, ctx);
    BN_mod_sub(new_p.z, new_p.z, yt, C.p, ctx);

    BN_free(yz);
    BN_free(tz);
    BN_free(yt);

    BN_CTX_free(ctx);

    return new_p;
}

struct Point multiple_point(struct Point P, BIGNUM * k){
    struct Point Q, R;
    struct Curve C;
    BIGNUM * help = BN_new();
    BIGNUM * zero = BN_new();

    //printf("\n>> MulP : k = %s\n", BN_bn2dec(k));

    C.b = BN_new();
    C.p = BN_new();
    C.q = BN_new();

    init_curve(&C);
    BN_dec2bn(&help, "1");
    BN_dec2bn(&zero, "0");
    init_point(&Q, zero, help, help, help);
    init_point(&R, P.x, P.y, P.t, P.z);
    int n = BN_num_bytes(k);
    for (int i = n - 1; i >= 0; i--){
        if (BN_is_bit_set(k, i) == 0){
            R = add_points(R, Q);
            Q = double_point(Q);
        }
        else {
            Q = add_points(Q, R);
            R = double_point(R);
        }
    }

    BN_free(help);
    return Q;
}

int if_contains(struct Point P){
    BIGNUM * left1;
    BIGNUM * right;
    BIGNUM * left2;
    BIGNUM * help;

    left1 = BN_new();
    left2 = BN_new();
    right = BN_new();
    help = BN_new();

    struct Curve C;

    C.b = BN_new();
    C.p = BN_new();
    C.q = BN_new();

    init_curve(&C);

    BN_CTX *ctx = BN_CTX_new();
    //printf("Before calc :\n> left1 = %s\nleft2 = %s\nright = %s\n", BN_bn2dec(left1), BN_bn2dec(left2), BN_bn2dec(right));
    //left1
    BN_mul(left1, P.x, P.x, ctx);
    BN_copy(left2, left1);
    BN_mul(help, P.y, P.y, ctx);
    BN_add(left1, left1, help);
    BN_mod(left1, left1, C.p, ctx);
    //printf("LEFT1 CHANGED :\n> left1 = %s\nleft2 = %s\nright = %s\n", BN_bn2dec(left1), BN_bn2dec(left2), BN_bn2dec(right));

    //left2
    BN_mul(left2, left2, C.b, ctx);
    BN_mul(help, P.t, P.t, ctx);
    BN_add(left2, left2, help);
    BN_mod(left2, left2, C.p, ctx);
    //printf("LEFT2 CHANGED :\n> left1 = %s\nleft2 = %s\nright = %s\n", BN_bn2dec(left1), BN_bn2dec(left2), BN_bn2dec(right));

    //right
    BN_mul(right, P.z, P.z,ctx);
    BN_mod(right, right, C.p, ctx);
    BN_CTX_free(ctx);
    //printf("RIGHT CHANGED :\n> left1 = %s\nleft2 = %s\nright = %s\n", BN_bn2dec(left1), BN_bn2dec(left2), BN_bn2dec(right));

    //printf(">> if_con : \n> left1 = %s\n> left2 = %s\n> right = %s\n", BN_bn2dec(left1), BN_bn2dec(left2), BN_bn2dec(right));
    //printf(">> if_con : \n> CMP left1 & right = %d\n", BN_cmp(left1, right));
    //printf(">> if_con : \n> CMP left2 & right = %d\n", BN_cmp(left2, right));

    if(BN_cmp(left1, right) == 0 && BN_cmp(left2, right) == 0) {
        BN_free(left1);
        BN_free(right);
        BN_free(left2);
        BN_free(help);
        return 1;
    }
    else {
        BN_free(left1);
        BN_free(right);
        BN_free(left2);
        BN_free(help);
        return 0;
    }
}
