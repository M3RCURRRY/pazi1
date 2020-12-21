#include <stdio.h>
#include "jacobi.h"

int main()
{
    struct Curve C;
    struct Point P;
    init_curve(&C);

    BIGNUM* x1, y1, t1, z1;
    x1 = BN_new();
    y1 = BN_new();
    t1 = BN_new();
    z1 = BN_new();
    BN_set_word(x1, x_s);
    BN_set_word(y1, y_s);
    BN_set_word(z1, t_s);
    BN_set_word(t1, z_s);
    init_point(&P, x1, y1, t1, z1);

    //���� 1
    BIGNUM* k, n;
    BN_init(k);
    BN_init(n);
    BN_set_word(n, "100000000000000000000000000", 10);
    BN_CTX* state;
    BN_rand(state, n);

    printf("1 ����. �������� �� ����� [k]P ������ ������: ");
    if (if_contains(multiple_point(P, k)) == 1)
        printf("��\n");
    else
        printf("���\n");

    //���� 2
    struct Point N;
    BN_set_word(x1, "1", 10);
    BN_set_word(y1, "0", 10);
    init_point(&N, y1, x1, x1, x1);

    struct Point Test_2 = multiple_point(P, C.q);

    printf("2 ����. ����� �� ����������� ����� � [q]P: ");
    if (N.x == Test_2.x && N.y == Test_2.y && N.t == Test_2.t && N.z == Test_2.z)
        printf("��\n");
    else
        printf("���\n");

    //���� 3
    printf("3.1 ����. ����� ��, ��� [q+1]P = P: ");
    BIGNUM* qq;
    BN_init(qq);
    BN_set_word(qq, "1", 10);
    BN_add(qq, qq, C.q);
    struct Point Test_3 = multiple_point(P, qq);
    if (P.x == Test_3.x && P.y == Test_3.y && P.t == Test_3.t && P.z == Test_3.z)
        printf("��\n");
    else
        printf("���\n");

    printf("3.2 ����. ����� ��, ��� [q-1]P = -P: ");
    BN_set_word(qq, "1", 10);
    BN_sub(qq, C.q, qq);
    Test_3 = multiple_point(P, qq);
    BN_set_word(qq, "-1", 10);
    BN_mul(qq, qq, P.x);
    if (qq == Test_3.x && P.y == Test_3.y && P.t == Test_3.t && P.z == Test_3.z)
        printf("��\n");
    else
        printf("���\n");

    //���� 4
    BIGNUM* k1, k2;
    BN_init(k1);
    BN_init(k2);
    BN_pseudo_rand(k1, state, n);
    BN_pseudo_rand(k2, state, n);
    struct Point Right, Left;
    Right = multiple_point(P, k2);
    Left = multiple_point(P, k1);
    Left = add_points(Left, Right);
    BN_add(k1, k1, k2);
    Right = multiple_point(P, k1);
    printf("4 ����. ����� ��, ��� [k1]P + [k2]P = [k1+k2]P: ");
    if (Left.x == Right.x && Left.y == Right.y && Left.t == Right.t && Left.z == Right.z)
        printf("��\n");
    else
        printf("���\n");

    BN_clear(x1);
    BN_clear(y1);
    BN_clear(t1);
    BN_clear(z1);
    BN_clear(k);
    BN_clear(n);
    BN_clear(qq);
    BN_clear(k1);
    BN_clear(k2);
}