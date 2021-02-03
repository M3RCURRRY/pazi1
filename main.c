#include <stdio.h>
#include "jacobi.h"

int main()
{
    const char * x_s = "86773867670483746082709700497332781598034016226947998080571074369332058445893";
    const char * y_s = "112365783073842736181650671369471258524588261898469739668735228897200995462257";
    const char * t_s = "34969497036250894576602400300903320615522974325083826128152497107656115683806";
    const char * z_s = "53269948314907700562817814682691121416875131851387647397949804668930162464909";

    struct Curve C;
    struct Point P;

    P.x = BN_new();
    P.y = BN_new();
    P.t = BN_new();
    P.z = BN_new();

    C.q = BN_new();
    C.p = BN_new();
    C.b = BN_new();

    init_curve(&C);

    /*
    struct Point QQ;
    BIGNUM * zro = BN_new();
    BIGNUM * one = BN_new();
    BN_dec2bn(&zro, "0");
    BN_dec2bn(&one, "1");
    init_point(&QQ, zro, one, one, one);
    */

    BIGNUM * x1 = BN_new();
    BIGNUM * y1 = BN_new();
    BIGNUM * t1 = BN_new();
    BIGNUM * z1 = BN_new();
    BN_dec2bn(&x1, x_s);
    BN_dec2bn(&y1, y_s);
    BN_dec2bn(&z1, t_s);
    BN_dec2bn(&t1, z_s);
    init_point(&P, x1, y1, t1, z1);

    //Tест 1
    BIGNUM * k = BN_new();
    BIGNUM * n = BN_new();
    BN_dec2bn(&n, "100000000000000000000000000");
    BN_rand_range(k,n);
    printf("1 тест. Является ли точка [k]P частью кривой: ");
    if (if_contains(multiple_point(P, k)) == 1) {
        printf("да\n");
    }
    else
        printf("нет\n");

    //Tест 2
    struct Point N;

    N.x = BN_new();
    N.y = BN_new();
    N.t = BN_new();
    N.z = BN_new();

    BN_dec2bn(&x1, "1");
    BN_dec2bn(&y1, "0");
    init_point(&N, y1, x1, x1, x1);

    struct Point Test_2;

    Test_2.x = BN_new();
    Test_2.y = BN_new();
    Test_2.t = BN_new();
    Test_2.z = BN_new();
    Test_2 = multiple_point(P, C.q);

    printf("2 тест. Равны ли нейтральная точка и [q]P: ");
    if (BN_cmp(N.x, Test_2.x) && BN_cmp(N.y, Test_2.y) && BN_cmp(N.t, Test_2.t) && BN_cmp(N.z, Test_2.z))
        //if (N.x == Test_2.x && N.y == Test_2.y && N.t == Test_2.t && N.z == Test_2.z)
        printf("да\n");
    else
        printf("нет\n");

    //Tест 3
    printf("3.1 тест. Верно ли, что [q+1]P = P: ");
    BIGNUM * qq = BN_new();
    BN_dec2bn(&qq, "1");
    BN_add(qq, qq, C.q);
    struct Point Test_3;

    Test_3.x = BN_new();
    Test_3.y = BN_new();
    Test_3.t = BN_new();
    Test_3.z = BN_new();

    Test_3 = multiple_point(P, qq);
    if (BN_cmp(P.x,Test_3.x) && BN_cmp(P.y, Test_3.y) && BN_cmp(P.t, Test_3.t) && BN_cmp(P.z, Test_3.z))
        printf("да\n");
    else
        printf("нет\n");

    printf("3.2 тест. Верно ли, что [q-1]P = -P: ");
    BN_dec2bn(&qq, "1");
    BN_sub(qq, C.q, qq);
    Test_3 = multiple_point(P, qq);
    BN_dec2bn(&qq, "-1");
    BN_CTX *m_ctx = BN_CTX_new();
    BN_mul(qq, qq, P.x, m_ctx);
    if (BN_cmp(qq, Test_3.x) && BN_cmp(P.y, Test_3.y) && BN_cmp(P.t, Test_3.t) && BN_cmp(P.z, Test_3.z))
        printf("да\n");
    else
        printf("нет\n");

    //Tест 4
    BIGNUM * k1 = BN_new();
    BIGNUM * k2 = BN_new();
    BN_rand_range(k1, n);
    BN_rand_range(k2, n);
    struct Point Right, Left;

    Right.x = BN_new();
    Right.y = BN_new();
    Right.t = BN_new();
    Right.z = BN_new();

    Left.x = BN_new();
    Left.y = BN_new();
    Left.t = BN_new();
    Left.z = BN_new();

    Right = multiple_point(P, k2);
    Left = multiple_point(P, k1);
    Left = add_points(Left, Right);
    BN_add(k1, k1, k2);
    Right = multiple_point(P, k1);
    printf("4 тест. Верно ли, что [k1]P + [k2]P = [k1+k2]P: ");
    if (BN_cmp(Left.x, Right.x) && BN_cmp(Left.y, Right.y) && BN_cmp(Left.t, Right.t) && BN_cmp(Left.z, Right.z))
        printf("да\n");
    else
        printf("нет\n");

    BN_free(x1);
    BN_free(y1);
    BN_free(t1);
    BN_free(z1);
    BN_free(k);
    BN_free(n);
    BN_free(qq);
    BN_free(k1);
    BN_free(k2);

    BN_CTX_free(m_ctx);
}
