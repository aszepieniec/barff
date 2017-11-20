#include "gfp.h"
#include "gfpm.h"
#include "hqs.h"
#include "csprng.h"
#include "gf2x.h"
#include "gf256x.h"
#include "gf65536x.h"
#include "gf16777216x.h"
#include "gf4096x.h"

#include <stdio.h>
#include <stdlib.h>

int test_csprng( )
{
    unsigned short int i;
    unsigned int r;
    int success, allequal;
    unsigned char o1[50];
    unsigned char o2[50];
    csprng rng;

    printf("testing csprng ...");
    success = 1;
    csprng_init(&rng);
    r = rand();
    csprng_seed(&rng, sizeof(unsigned int), (unsigned char *)&r);
    csprng_generate(&rng, 50, o1);
    csprng_generate(&rng, 50, o2);

    allequal = 1;
    for( i = 0 ; i < 50 ; ++i )
    {
        if( o1[i] != o2[i] )
        {
            allequal = 0;
        }
    }
    if( allequal != 0 )
    {
        success = 0;
    }

    if( success == 0 )
    {
        printf( "fail! Did not produce 2 different output streams.\n");
        return 0;
    }

    csprng_init(&rng);
    csprng_seed(&rng, sizeof(unsigned int), (unsigned char *)&r);
    csprng_generate(&rng, 50, o2);

    allequal = 1;
    for( i = 0 ; i < 50 ; ++i )
    {
        if( o1[i] != o2[i] )
        {
            allequal = 0;
        }
    }
    if( allequal == 0 )
    {
        success = 0;
    }

    if( success == 1 )
    {
        printf("success!\n");
        return 1;
    }
    else
    {
        printf("fail!\n");
        return 0;
    }
}

int test_matrix_inverse( )
{
    FILE * fh;
    csprng rng;

    unsigned char * randomness;
    char * prand;
    int invertible;
    int num_limbs;

    gfpmatrix A, B, C;

    unsigned int random;
    random = rand();
    csprng_init(&rng);
    csprng_seed(&rng, sizeof(unsigned int), (unsigned char*)&random);

    A = gfpm_init(10,10);
    B = gfpm_init(10,10);
    C = gfpm_init(10,10);


    printf("testing matrix inverse ... ");

    num_limbs = (GFP_NUMBITS + sizeof(unsigned long int) * 8 - 1) / (sizeof(unsigned long int) * 8);
    randomness = malloc(sizeof(unsigned long int)*num_limbs*10*10);
    csprng_generate(&rng, sizeof(unsigned long int)*num_limbs*10*10, randomness);
    gfpm_random(A, randomness);

    invertible = gfpm_inverse(B, A);

    gfpm_multiply(&C, A, B);

    free(randomness);

    if( invertible == 0 || gfpm_is_eye(C) )
    {
        printf("success!\n");
        gfpm_destroy(A);
        gfpm_destroy(B);
        gfpm_destroy(C);
        return 1;
    }
    else
    {
        printf("fail!\n");
        printf("A:\n");
        gfpm_print(A);
        printf("B:\n");
        gfpm_print(B);
        printf("C:\n");
        gfpm_print(C);
        printf("random seed: %u\n", random);
        gfpm_destroy(A);
        gfpm_destroy(B);
        gfpm_destroy(C);
        return 0;
    }
}

int test_multiply_transpose( )
{
    gfpmatrix A, B, C1, C2;
    unsigned short int n, m;
    csprng rng;
    unsigned int random;
    unsigned char * randomness;
    int num_limbs;

    random = rand();
    csprng_init(&rng);
    csprng_seed(&rng, sizeof(unsigned int), (unsigned char *)&random);

    printf("testing transpose-multiplication ... ");

    m = 10;
    n = 5;

    A = gfpm_init(n, m);
    B = gfpm_init(n, m);
    C1 = gfpm_init(n, n);
    C2 = gfpm_init(n, n);

    num_limbs = (n*m*(GFP_NUMBITS + sizeof(unsigned long int)*8 - 1)) / (sizeof(unsigned long int) * 8);
    randomness = malloc(num_limbs*sizeof(unsigned long int));
    csprng_generate(&rng, (num_limbs*sizeof(unsigned long int)), randomness);
    gfpm_random(A, randomness);
    csprng_generate(&rng, (num_limbs*sizeof(unsigned long int)), randomness);
    gfpm_random(B, randomness);
    free(randomness);

    gfpm_multiply_transpose(&C1, A, B);
    gfpm_transpose(&B);
    gfpm_multiply(&C2, A, B);
    gfpm_transpose(&B);

    if( !gfpm_equals(C1, C2) )
    {
        printf("fail!\n");
        printf("C1 =/= C2 after gfpm_multiply_transpose\n");
        printf("rand seed: %i\n", random);

        gfpm_destroy(A);
        gfpm_destroy(B);
         gfpm_destroy(C1);
         gfpm_destroy(C2);
        return 0;
    }
    
    gfpm_destroy(C1);
    gfpm_destroy(C2);
    C1 = gfpm_init(m, m);
    C2 = gfpm_init(m, m);

    gfpm_transpose_multiply(&C1, A, B);
    gfpm_transpose(&A);
    gfpm_multiply(&C2, A, B);
    gfpm_transpose(&A);

    if( !gfpm_equals(C1, C2) )
    {
        printf("fail!\n");
        printf("C1 =/= C2 after gfpm_transpose_multiply\n");
        printf("rand seed: %i\n", random);

        gfpm_destroy(A);
        gfpm_destroy(B);
         gfpm_destroy(C1);
         gfpm_destroy(C2);
        return 0;
    }

    printf("success!\n");

    gfpm_destroy(A);
    gfpm_destroy(B);
     gfpm_destroy(C1);
     gfpm_destroy(C2);
    return 1;
}

int test_solve( )
{
    unsigned short int m, n;
    unsigned int random;
    gfpmatrix A, b, b2, x, s, K, k, v;
    csprng rng;
    int success;
    unsigned char * randomness;
    int num_limbs;

    success = 1;

    random = rand();
    csprng_init(&rng);
    csprng_seed(&rng, sizeof(unsigned int), (unsigned char *)&random);

    m = 5;
    n = 6;

    printf("testing matrix equation solver ... ");

    A = gfpm_init(m, n);
    num_limbs = (GFP_NUMBITS + sizeof(unsigned long int) * 8 - 1) / (sizeof(unsigned long int) * 8);
    randomness = malloc(m*n*sizeof(unsigned long int)*num_limbs);
    csprng_generate(&rng, m*n*sizeof(unsigned long int)*num_limbs, randomness);
    gfpm_random(A, randomness);
    free(randomness);
    b = gfpm_init(m, 1);
    x = gfpm_init(n, 1);
    
    randomness = malloc(n*1*sizeof(unsigned long int)*num_limbs);
    csprng_generate(&rng, n*1*sizeof(unsigned long int)*num_limbs, randomness);
    gfpm_random(x, randomness);
    free(randomness);
    gfpm_multiply(&b, A, x);
    gfpm_zeros(x);


    s = gfpm_init(n, 1);

    gfpm_solve(A, b, s, &K);

    /* get a kernel vector */
    k = gfpm_init(n, 1);
    if( K.width > 0 )
    {
        v = gfpm_init(K.width, 1);
        randomness = malloc(K.width*1*sizeof(unsigned long int)*num_limbs);
        csprng_generate(&rng, K.width*1*sizeof(unsigned long int)*num_limbs, randomness);
        gfpm_random(v, randomness);
        free(randomness);

        gfpm_multiply(&k, K, v);

        gfpm_destroy(v);
    }
    else
    {
        gfpm_zeros(k);
    }


    gfpm_add(x, s, k);
    b2 = gfpm_init(m,1);
    gfpm_multiply(&b2, A, x);

    if( !gfpm_equals(b, b2) )
    {
        printf("fail!\n"); getchar();
        printf("A:\n");
        gfpm_print(A);
        printf("s: ");
        gfpm_print_transpose(s);
        printf("K:\n");
        gfpm_print(K);

        printf("k^T: ");
        gfpm_print_transpose(k);
        printf("b^T:  ");
        gfpm_print_transpose(b);
        printf("b2^T: ");
        gfpm_print_transpose(b2);
        printf("x: ");
        gfpm_print_transpose(x);
        printf("random seed: %i\n", random);
        success = 0;
    }

    randomness = malloc(b.height*b.width*sizeof(unsigned long int)*num_limbs);
    csprng_generate(&rng, b.height*b.width*sizeof(unsigned long int)*num_limbs, randomness);
    gfpm_random(b, randomness);
    free(randomness);
    gfpm_destroy(K);
    if( gfpm_solve(A, b, x, &K) && success == 1 )
    {
        gfpm_multiply(&b2, A, x);
        if( !gfpm_equals(b, b2) )
        {
            printf("fail!\n");
            printf("Ax = b is supposed to have solution but --\n");
            printf("x =  ");
            gfpm_print_transpose(x);
            printf("b =  ");
            gfpm_print_transpose(b);
            printf("Ax = ");
            gfpm_print_transpose(b2);
            printf("random seed: %i\n", random);
            success = 0;
        }
    }

    if( success == 1 )
    {
        printf("success!\n");
    }

    gfpm_destroy(A);
    gfpm_destroy(b);
    gfpm_destroy(b2);
    gfpm_destroy(K);
    gfpm_destroy(x);
    gfpm_destroy(s);
    gfpm_destroy(k);

    return success;
}

int test_composition( )
{
    hqsystem F, P;
    gfpmatrix T, S;
    unsigned short int m, n;
    gfpmatrix x, Sx, FoSx, ToFoSx, y;
    unsigned int i;
    int equal;
    csprng rng;
    unsigned int random;
    unsigned char * randomness;
    int num_limbs;

    random = rand();
    csprng_init(&rng);
    csprng_seed(&rng, sizeof(unsigned int), (unsigned char *)&random);


    printf("testing composition of linear transforms with homogeneous quadratic systems ... ");


    m = 3;
    n = 4;

    F = hqs_init(n, m);
    num_limbs = (GFP_NUMBITS + sizeof(unsigned long int) * 8 - 1) / (sizeof(unsigned long int) * 8);
    randomness = malloc(num_limbs * sizeof(unsigned long int) * m * n * n);
    csprng_generate(&rng, num_limbs * sizeof(unsigned long int)* m * n * n, randomness);
    hqs_random(F, randomness);
    free(randomness);
    P = hqs_clone(F);

    T = gfpm_init(m, m);
    S = gfpm_init(n, n);


    randomness = malloc(num_limbs * sizeof(unsigned long int) * m * m);
    csprng_generate(&rng, num_limbs * sizeof(unsigned long int) * m * m, randomness);
    gfpm_random_invertible(T, randomness);
    free(randomness);
    randomness = malloc(num_limbs * sizeof(unsigned long int) * n * n);
    csprng_generate(&rng, num_limbs * sizeof(unsigned long int)* n * n, randomness);
    gfpm_random_invertible(S, randomness);
    free(randomness);

    hqs_compose_output(T, P);

    hqs_compose_input(P, S);

    x = gfpm_init(n, 1);
    y = gfpm_init(m, 1);
    Sx = gfpm_init(n, 1);
    FoSx = gfpm_init(m, 1);
    ToFoSx = gfpm_init(m, 1);

    equal = 1;
    for( i = 0 ; i < 1 && equal == 1 ; ++i )
    {
        randomness = malloc(num_limbs * sizeof(unsigned long int) * x.height * x.width);
        csprng_generate(&rng, num_limbs * sizeof(unsigned long int) * x.height * x.width, randomness);
        gfpm_random(x, randomness);
        free(randomness);

        gfpm_multiply(&Sx, S, x);

        hqs_eval(FoSx, F, Sx);

        gfpm_multiply(&ToFoSx, T, FoSx);

        hqs_eval(y, P, x);

        equal = gfpm_equals(y, ToFoSx);
    }

    gfpm_destroy(x);
    gfpm_destroy(y);
    gfpm_destroy(Sx);
    gfpm_destroy(FoSx);
    gfpm_destroy(ToFoSx);
    hqs_destroy(F);
    hqs_destroy(P);
    gfpm_destroy(T);
    gfpm_destroy(S);

    if( equal )
    {
        printf("success!\n");
        return 1;
    }
    else
    {
        printf("fail!\n");
        printf("random seed: %u\n", random);
        return 0;
    }
}

int test_gf256x_add( )
{
    unsigned short int m, n, o;
    unsigned int i;
    int equal;
    gf256x a, b, c, ab, abc1, bc, abc2;
    csprng rng;
    unsigned int random;
    unsigned char * randomness;
    int equals;

    random = rand();
    csprng_init(&rng);
    csprng_seed(&rng, sizeof(unsigned int), (unsigned char *)&random);


    m = csprng_generate_ulong(&rng) % 100;
    n = csprng_generate_ulong(&rng) % 100;
    o = csprng_generate_ulong(&rng) % 100;

    printf("testing addition of GF(256)[x] elements of degrees %i, %i, and %i... ", m, n, o);

    a = gf256x_init(m);
    csprng_generate(&rng, a.degree+1, a.data);

    b = gf256x_init(n);
    csprng_generate(&rng, b.degree+1, b.data);

    c = gf256x_init(o);
    csprng_generate(&rng, c.degree+1, c.data);

    ab = gf256x_init(0);
    gf256x_add(&ab, a, b);
    abc1 = gf256x_init(0);
    gf256x_add(&abc1, ab, c);

    bc = gf256x_init(0);
    gf256x_add(&bc, b, c);
    abc2 = gf256x_init(0);
    gf256x_add(&abc2, ab, c);

    equals = gf256x_equals(abc1, abc2);

    gf256x_destroy(a);
    gf256x_destroy(b);
    gf256x_destroy(c);
    gf256x_destroy(ab);
    gf256x_destroy(abc1);
    gf256x_destroy(bc);
    gf256x_destroy(abc2);

    if( equals == 1 )
    {
        printf("success!\n");
        return 1;
    }
    else
    {
        printf("fail!\n");
        return 0;
    }
}

int test_gf256x_multiply( )
{
    unsigned short int m, n, o;
    unsigned int i;
    int equal;
    gf256x a, b, c, ab, abc1, bc, abc2;
    csprng rng;
    unsigned int random;
    unsigned char * randomness;
    int equals;

    random = rand();
    csprng_init(&rng);
    csprng_seed(&rng, sizeof(unsigned int), (unsigned char *)&random);


    m = csprng_generate_ulong(&rng) % 100;
    n = csprng_generate_ulong(&rng) % 100;
    o = csprng_generate_ulong(&rng) % 100;

    printf("testing multiplication of GF(256)[x] elements of degrees %i, %i, and %i... ", m, n, o);

    a = gf256x_init(m);
    csprng_generate(&rng, a.degree+1, a.data);

    b = gf256x_init(n);
    csprng_generate(&rng, b.degree+1, b.data);

    c = gf256x_init(o);
    csprng_generate(&rng, c.degree+1, c.data);

    ab = gf256x_init(0);
    gf256x_multiply(&ab, a, b);
    abc1 = gf256x_init(0);
    gf256x_multiply(&abc1, ab, c);

    bc = gf256x_init(0);
    gf256x_multiply(&bc, b, c);
    abc2 = gf256x_init(0);
    gf256x_multiply(&abc2, ab, c);

    equals = gf256x_equals(abc1, abc2);

    gf256x_destroy(a);
    gf256x_destroy(b);
    gf256x_destroy(c);
    gf256x_destroy(ab);
    gf256x_destroy(abc1);
    gf256x_destroy(bc);
    gf256x_destroy(abc2);

    if( equals == 1 )
    {
        printf("success!\n");
        return 1;
    }
    else
    {
        printf("fail!\n");
        return 0;
    }
}

int test_gf256x_divide( )
{
    unsigned short int m, n, o;
    unsigned int i;
    int equal;
    gf256x numerator, divisor, quotient, remainder, product, sum;
    csprng rng;
    unsigned int random;
    unsigned char * randomness;
    int equals;

    random = rand();
    csprng_init(&rng);
    csprng_seed(&rng, sizeof(unsigned int), (unsigned char *)&random);


    m = csprng_generate_ulong(&rng) % 200;
    n = csprng_generate_ulong(&rng) % m;

    printf("testing division of GF(256)[x] elements of degrees %i, and %i... ", m, n);

    numerator = gf256x_init(m);
    csprng_generate(&rng, numerator.degree+1, numerator.data);

    divisor = gf256x_init(n);
    csprng_generate(&rng, divisor.degree+1, divisor.data);

    quotient = gf256x_init(0);
    remainder = gf256x_init(0);
    gf256x_divide(&quotient, &remainder, numerator, divisor);

    product = gf256x_init(0);
    gf256x_multiply(&product, divisor, quotient);

    sum = gf256x_init(0);
    gf256x_add(&sum, product, remainder);

    equals = gf256x_equals(numerator, sum);
    equals &= remainder.degree <= divisor.degree; /* includes division by constant */

    if( equals == 1 )
    {
        printf("success!\n");
    }
    else
    {
        printf("\n");
        printf("numerator: "); gf256x_print(numerator); printf("\n");
        printf("divisor: "); gf256x_print(divisor); printf("\n");
        printf("quotient: "); gf256x_print(quotient); printf("\n");
        printf("remainder: "); gf256x_print(remainder); printf("\n");
        printf("product: "); gf256x_print(product); printf("\n");
        printf("sum: "); gf256x_print(sum); printf("\n");
        printf("fail!\n");
    }

    gf256x_destroy(numerator);
    gf256x_destroy(divisor);
    gf256x_destroy(remainder);
    gf256x_destroy(quotient);
    gf256x_destroy(product);
    gf256x_destroy(sum);

    return equals;
}

int test_gf256x_xgcd( )
{
    unsigned short int m, n;
    unsigned int i;
    int equal;
    gf256x x, y, a, b, g, ax, by, sum, quotient, remainder;
    csprng rng;
    unsigned int random;
    unsigned char * randomness;
    int equals;

    random = rand();
    csprng_init(&rng);
    csprng_seed(&rng, sizeof(unsigned int), (unsigned char *)&random);


    m = csprng_generate_ulong(&rng) % 200;
    n = csprng_generate_ulong(&rng) % 200;

    printf("testing xgcd GF(256)[x] elements of degrees %i, and %i... ", m, n);

    x = gf256x_init(m);
    csprng_generate(&rng, x.degree+1, x.data);

    y = gf256x_init(n);
    csprng_generate(&rng, y.degree+1, y.data);

    quotient = gf256x_init(0);
    remainder = gf256x_init(0);
    a = gf256x_init(0);
    b = gf256x_init(0);
    g = gf256x_init(0);
    gf256x_xgcd(&a, &b, &g, x, y);

    ax = gf256x_init(0);
    by = gf256x_init(0);
    sum = gf256x_init(0);
    gf256x_multiply(&ax, a, x);
    gf256x_multiply(&by, b, y);
    gf256x_add(&sum, ax, by);
    equals = gf256x_equals(g, sum);

    gf256x_divide(&quotient, &remainder, x, g);
    equals &= gf256x_is_zero(remainder);

    gf256x_divide(&quotient, &remainder, y, g);
    equals &= gf256x_is_zero(remainder);

    if( equals == 1 )
    {
        printf("success!\n");
    }
    else
    {
        printf("fail!\n");
    }

    gf256x_destroy(x);
    gf256x_destroy(y);
    gf256x_destroy(a);
    gf256x_destroy(b);
    gf256x_destroy(g);
    gf256x_destroy(ax);
    gf256x_destroy(by);
    gf256x_destroy(sum);
    gf256x_destroy(remainder);
    gf256x_destroy(quotient);

    return equals;
}

int test_gf65536_inverse( )
{
    csprng rng;
    unsigned int random;
    int equals;
    unsigned int a, b, c;

    random = rand();
    csprng_init(&rng);
    csprng_seed(&rng, sizeof(unsigned int), (unsigned char *)&random);


    printf("testing inverse computation of GF(65536) elements ... ");

    a = csprng_generate_ulong(&rng) % 0xffff;

    b = gf65536_inverse(a);

    c = gf65536_multiply(a, b);

    if( c == 1 || a == 0 )
    {
        printf("success!\n");
        return 1;
    }
    else
    {
        printf("fail!\n");
        return 0;
    }
}

int test_gf65536x_add( )
{
    unsigned short int m, n, o;
    unsigned int i;
    int equal;
    gf65536x a, b, c, ab, abc1, bc, abc2;
    csprng rng;
    unsigned int random;
    unsigned char * randomness;
    int equals;

    random = rand();
    csprng_init(&rng);
    csprng_seed(&rng, sizeof(unsigned int), (unsigned char *)&random);


    m = csprng_generate_ulong(&rng) % 100;
    n = csprng_generate_ulong(&rng) % 100;
    o = csprng_generate_ulong(&rng) % 100;

    printf("testing addition of GF(65536)[x] elements of degrees %i, %i, and %i... ", m, n, o);

    a = gf65536x_init(m);
    csprng_generate(&rng, 2*a.degree+2, a.data);

    b = gf65536x_init(n);
    csprng_generate(&rng, 2*b.degree+2, b.data);

    c = gf65536x_init(o);
    csprng_generate(&rng, 2*c.degree+2, c.data);

    ab = gf65536x_init(0);
    gf65536x_add(&ab, a, b);
    abc1 = gf65536x_init(0);
    gf65536x_add(&abc1, ab, c);

    bc = gf65536x_init(0);
    gf65536x_add(&bc, b, c);
    abc2 = gf65536x_init(0);
    gf65536x_add(&abc2, ab, c);

    equals = gf65536x_equals(abc1, abc2);

    gf65536x_destroy(a);
    gf65536x_destroy(b);
    gf65536x_destroy(c);
    gf65536x_destroy(ab);
    gf65536x_destroy(abc1);
    gf65536x_destroy(bc);
    gf65536x_destroy(abc2);

    if( equals == 1 )
    {
        printf("success!\n");
        return 1;
    }
    else
    {
        printf("fail!\n");
        return 0;
    }
}

int test_gf65536x_multiply( )
{
    unsigned short int m, n, o;
    unsigned int i;
    int equal;
    gf65536x a, b, c, ab, abc1, bc, abc2;
    csprng rng;
    unsigned int random;
    unsigned char * randomness;
    int equals;

    random = rand();
    csprng_init(&rng);
    csprng_seed(&rng, sizeof(unsigned int), (unsigned char *)&random);


    m = csprng_generate_ulong(&rng) % 100;
    n = csprng_generate_ulong(&rng) % 100;
    o = csprng_generate_ulong(&rng) % 100;

    printf("testing multiplication of GF(65536)[x] elements of degrees %i, %i, and %i... ", m, n, o);

    a = gf65536x_init(m);
    csprng_generate(&rng, 2*a.degree+2, a.data);

    b = gf65536x_init(n);
    csprng_generate(&rng, 2*b.degree+2, b.data);

    c = gf65536x_init(o);
    csprng_generate(&rng, 2*c.degree+2, c.data);

    ab = gf65536x_init(0);
    gf65536x_multiply(&ab, a, b);
    abc1 = gf65536x_init(0);
    gf65536x_multiply(&abc1, ab, c);

    bc = gf65536x_init(0);
    gf65536x_multiply(&bc, b, c);
    abc2 = gf65536x_init(0);
    gf65536x_multiply(&abc2, ab, c);

    equals = gf65536x_equals(abc1, abc2);
    equals &= (a.degree + b.degree + c.degree) == abc1.degree;

    gf65536x_destroy(a);
    gf65536x_destroy(b);
    gf65536x_destroy(c);
    gf65536x_destroy(ab);
    gf65536x_destroy(abc1);
    gf65536x_destroy(bc);
    gf65536x_destroy(abc2);

    if( equals == 1 )
    {
        printf("success!\n");
        return 1;
    }
    else
    {
        printf("fail!\n");
        return 0;
    }
}

int test_gf65536x_divide( )
{
    unsigned short int m, n, o;
    unsigned int i;
    int equal;
    gf65536x numerator, divisor, quotient, remainder, product, sum;
    csprng rng;
    unsigned int random;
    unsigned char * randomness;
    int equals;

    random = rand();
    csprng_init(&rng);
    csprng_seed(&rng, sizeof(unsigned int), (unsigned char *)&random);


    m = csprng_generate_ulong(&rng) % 200;
    n = csprng_generate_ulong(&rng) % m;

    printf("testing division of GF(65536)[x] elements of degrees %i, and %i... ", m, n);

    numerator = gf65536x_init(m);
    csprng_generate(&rng, 2*numerator.degree+2, numerator.data);

    divisor = gf65536x_init(n);
    csprng_generate(&rng, 2*divisor.degree+2, divisor.data);

    quotient = gf65536x_init(0);
    remainder = gf65536x_init(0);
    gf65536x_divide(&quotient, &remainder, numerator, divisor);
    printf("\n");

    product = gf65536x_init(0);
    gf65536x_multiply(&product, divisor, quotient);

    sum = gf65536x_init(0);
    gf65536x_add(&sum, product, remainder);

    equals = gf65536x_equals(numerator, sum);
    equals &= remainder.degree <= divisor.degree; /* includes division by constant */

    if( equals == 1 )
    {
        printf("success!\n");
    }
    else
    {
        printf("\n");
        printf("numerator: "); gf65536x_print(numerator); printf("\n");
        printf("divisor: "); gf65536x_print(divisor); printf("\n");
        printf("quotient: "); gf65536x_print(quotient); printf("\n");
        printf("remainder: "); gf65536x_print(remainder); printf("\n");
        printf("product: "); gf65536x_print(product); printf("\n");
        printf("sum: "); gf65536x_print(sum); printf("\n");
        gf65536x_add(&sum, sum, numerator);
        printf("difference: "); gf65536x_print(sum); printf("\n");
        printf("fail!\n");
    }

    gf65536x_destroy(numerator);
    gf65536x_destroy(divisor);
    gf65536x_destroy(remainder);
    gf65536x_destroy(quotient);
    gf65536x_destroy(product);
    gf65536x_destroy(sum);

    return equals;
}

int test_gf65536x_xgcd( )
{
    unsigned short int m, n;
    unsigned int i;
    int equal;
    gf65536x x, y, a, b, g, ax, by, sum, quotient, remainder;
    csprng rng;
    unsigned int random;
    unsigned char * randomness;
    int equals;

    random = rand();
    csprng_init(&rng);
    csprng_seed(&rng, sizeof(unsigned int), (unsigned char *)&random);


    m = csprng_generate_ulong(&rng) % 200;
    n = csprng_generate_ulong(&rng) % 200;

    printf("randomness: %lu ...", random);

    printf("testing xgcd GF(65536)[x] elements of degrees %i, and %i... ", m, n);

    x = gf65536x_init(m);
    csprng_generate(&rng, 2*x.degree+2, x.data);

    y = gf65536x_init(n);
    csprng_generate(&rng, 2*y.degree+2, y.data);

    quotient = gf65536x_init(0);
    remainder = gf65536x_init(0);
    a = gf65536x_init(0);
    b = gf65536x_init(0);
    g = gf65536x_init(0);
    gf65536x_xgcd(&a, &b, &g, x, y);

    ax = gf65536x_init(0);
    by = gf65536x_init(0);
    sum = gf65536x_init(0);
    gf65536x_multiply(&ax, a, x);
    gf65536x_multiply(&by, b, y);
    gf65536x_add(&sum, ax, by);
    equals = gf65536x_equals(g, sum);

    gf65536x_divide(&quotient, &remainder, x, g);
    equals &= gf65536x_is_zero(remainder);

    gf65536x_divide(&quotient, &remainder, y, g);
    equals &= gf65536x_is_zero(remainder);

    if( equals == 1 )
    {
        printf("success!\n");
    }
    else
    {
        printf("fail!\n");
        printf("x: "); gf65536x_print(x); printf("\n");
        printf("y: "); gf65536x_print(y); printf("\n");
        printf("a: "); gf65536x_print(a); printf("\n");
        printf("b: "); gf65536x_print(b); printf("\n");
        printf("g: "); gf65536x_print(g); printf("\n");

    gf65536x_multiply(&ax, a, x);
    gf65536x_multiply(&by, b, y);
    gf65536x_add(&sum, ax, by);

        if( gf65536x_equals(g, sum) == 1 )
        {
            printf("sum equals\n");
        }
        else
        {
            printf("sum is different\n");
            printf("ax: "); gf65536x_print(ax); printf("\n");
            printf("by: "); gf65536x_print(by); printf("\n");
            printf("sum: "); gf65536x_print(sum); printf("\n");
            printf("g = "); gf65536x_print(g); printf("\n");
        }

    gf65536x_divide(&quotient, &remainder, x, g);
        if( gf65536x_is_zero(remainder) == 1 )
        {
            printf("x is zero mod g\n");
        }
        else
        {
            printf("x is nonzero mod g\n");
        }

    gf65536x_divide(&quotient, &remainder, y, g);
        if( gf65536x_is_zero(remainder) == 1 )
        {
            printf("y is zero mod g\n");
        }
        else
        {
            printf("y is nonzero mod g\n");
        }
    }

    gf65536x_destroy(x);
    gf65536x_destroy(y);
    gf65536x_destroy(a);
    gf65536x_destroy(b);
    gf65536x_destroy(g);
    gf65536x_destroy(ax);
    gf65536x_destroy(by);
    gf65536x_destroy(sum);
    gf65536x_destroy(remainder);
    gf65536x_destroy(quotient);

    return equals;
}


int test_gf16777216_inverse( )
{
    csprng rng;
    unsigned int random;
    int equals;
    unsigned int a, b, c;

    random = rand();
    csprng_init(&rng);
    csprng_seed(&rng, sizeof(unsigned int), (unsigned char *)&random);


    printf("testing inverse computation of GF(16777216) elements ... ");

    a = csprng_generate_ulong(&rng) % 0xffffff;

    b = gf16777216_inverse(a);

    c = gf16777216_multiply(a, b);

    if( c == 1 || a == 0 )
    {
        printf("success!\n");
        return 1;
    }
    else
    {
        printf("a: %02x%02x%02x\n", a&0xff, (a>>8)&0xff, (a>>16)&0xff);
        printf("b: %02x%02x%02x\n", b&0xff, (b>>8)&0xff, (b>>16)&0xff);
        printf("c: %02x%02x%02x\n", c&0xff, (c>>8)&0xff, (c>>16)&0xff);
        printf("fail!\n");
        return 0;
    }
}

int test_gf16777216x_add( )
{
    unsigned short int m, n, o;
    unsigned int i;
    int equal;
    gf16777216x a, b, c, ab, abc1, bc, abc2;
    csprng rng;
    unsigned int random;
    unsigned char * randomness;
    int equals;

    random = rand();
    csprng_init(&rng);
    csprng_seed(&rng, sizeof(unsigned int), (unsigned char *)&random);


    m = csprng_generate_ulong(&rng) % 100;
    n = csprng_generate_ulong(&rng) % 100;
    o = csprng_generate_ulong(&rng) % 100;

    printf("testing addition of GF(16777216)[x] elements of degrees %i, %i, and %i... ", m, n, o);

    a = gf16777216x_init(m);
    csprng_generate(&rng, 3*a.degree+3, a.data);

    b = gf16777216x_init(n);
    csprng_generate(&rng, 3*b.degree+3, b.data);

    c = gf16777216x_init(o);
    csprng_generate(&rng, 3*c.degree+3, c.data);

    ab = gf16777216x_init(0);
    gf16777216x_add(&ab, a, b);
    abc1 = gf16777216x_init(0);
    gf16777216x_add(&abc1, ab, c);

    bc = gf16777216x_init(0);
    gf16777216x_add(&bc, b, c);
    abc2 = gf16777216x_init(0);
    gf16777216x_add(&abc2, ab, c);

    equals = gf16777216x_equals(abc1, abc2);

    gf16777216x_destroy(a);
    gf16777216x_destroy(b);
    gf16777216x_destroy(c);
    gf16777216x_destroy(ab);
    gf16777216x_destroy(abc1);
    gf16777216x_destroy(bc);
    gf16777216x_destroy(abc2);

    if( equals == 1 )
    {
        printf("success!\n");
        return 1;
    }
    else
    {
        printf("fail!\n");
        return 0;
    }
}

int test_gf16777216x_multiply( )
{
    unsigned short int m, n, o;
    unsigned int i;
    int equal;
    gf16777216x a, b, c, ab, abc1, bc, abc2;
    csprng rng;
    unsigned int random;
    unsigned char * randomness;
    int equals;

    random = rand();
    csprng_init(&rng);
    csprng_seed(&rng, sizeof(unsigned int), (unsigned char *)&random);


    m = csprng_generate_ulong(&rng) % 100;
    n = csprng_generate_ulong(&rng) % 100;
    o = csprng_generate_ulong(&rng) % 100;

    printf("testing multiplication of GF(16777216)[x] elements of degrees %i, %i, and %i... ", m, n, o);

    a = gf16777216x_init(m);
    csprng_generate(&rng, 3*a.degree+3, a.data);

    b = gf16777216x_init(n);
    csprng_generate(&rng, 3*b.degree+3, b.data);

    c = gf16777216x_init(o);
    csprng_generate(&rng, 3*c.degree+3, c.data);

    ab = gf16777216x_init(0);
    gf16777216x_multiply(&ab, a, b);
    abc1 = gf16777216x_init(0);
    gf16777216x_multiply(&abc1, ab, c);

    bc = gf16777216x_init(0);
    gf16777216x_multiply(&bc, b, c);
    abc2 = gf16777216x_init(0);
    gf16777216x_multiply(&abc2, ab, c);

    equals = gf16777216x_equals(abc1, abc2);
    equals &= (a.degree + b.degree + c.degree) == abc1.degree;

    gf16777216x_destroy(a);
    gf16777216x_destroy(b);
    gf16777216x_destroy(c);
    gf16777216x_destroy(ab);
    gf16777216x_destroy(abc1);
    gf16777216x_destroy(bc);
    gf16777216x_destroy(abc2);

    if( equals == 1 )
    {
        printf("success!\n");
        return 1;
    }
    else
    {
        printf("fail!\n");
        return 0;
    }
}

int test_gf16777216x_divide( )
{
    unsigned short int m, n, o;
    unsigned int i;
    int equal;
    gf16777216x numerator, divisor, quotient, remainder, product, sum;
    csprng rng;
    unsigned int random;
    unsigned char * randomness;
    int equals;

    random = rand();
    csprng_init(&rng);
    csprng_seed(&rng, sizeof(unsigned int), (unsigned char *)&random);


    m = csprng_generate_ulong(&rng) % 200;
    n = csprng_generate_ulong(&rng) % m;

    printf("testing division of GF(16777216)[x] elements of degrees %i, and %i... ", m, n);

    numerator = gf16777216x_init(m);
    csprng_generate(&rng, 3*numerator.degree+3, numerator.data);

    divisor = gf16777216x_init(n);
    csprng_generate(&rng, 3*divisor.degree+3, divisor.data);

    quotient = gf16777216x_init(0);
    remainder = gf16777216x_init(0);
    gf16777216x_divide(&quotient, &remainder, numerator, divisor);

    product = gf16777216x_init(0);
    gf16777216x_multiply(&product, divisor, quotient);

    sum = gf16777216x_init(0);
    gf16777216x_add(&sum, product, remainder);

    equals = gf16777216x_equals(numerator, sum);
    equals &= remainder.degree <= divisor.degree; /* includes division by constant */

    if( equals == 1 )
    {
        printf("success!\n");
    }
    else
    {
        printf("\n");
        printf("numerator: "); gf16777216x_print(numerator); printf("\n");
        printf("divisor: "); gf16777216x_print(divisor); printf("\n");
        printf("quotient: "); gf16777216x_print(quotient); printf("\n");
        printf("remainder: "); gf16777216x_print(remainder); printf("\n");
        printf("product: "); gf16777216x_print(product); printf("\n");
        printf("sum: "); gf16777216x_print(sum); printf("\n");
        gf16777216x_add(&sum, sum, numerator);
        printf("difference: "); gf16777216x_print(sum); printf("\n");
        printf("fail!\n");
    }

    gf16777216x_destroy(numerator);
    gf16777216x_destroy(divisor);
    gf16777216x_destroy(remainder);
    gf16777216x_destroy(quotient);
    gf16777216x_destroy(product);
    gf16777216x_destroy(sum);

    return equals;
}

int test_gf16777216x_xgcd( )
{
    unsigned short int m, n;
    unsigned int i;
    int equal;
    gf16777216x x, y, a, b, g, ax, by, sum, quotient, remainder;
    csprng rng;
    unsigned int random;
    unsigned char * randomness;
    int equals;

    random = rand();
    csprng_init(&rng);
    csprng_seed(&rng, sizeof(unsigned int), (unsigned char *)&random);


    m = csprng_generate_ulong(&rng) % 200;
    n = csprng_generate_ulong(&rng) % 200;

    printf("randomness: %lu ...", random);

    printf("testing xgcd GF(16777216)[x] elements of degrees %i, and %i... ", m, n);

    x = gf16777216x_init(m);
    csprng_generate(&rng, 3*x.degree+3, x.data);

    y = gf16777216x_init(n);
    csprng_generate(&rng, 3*y.degree+3, y.data);

    quotient = gf16777216x_init(0);
    remainder = gf16777216x_init(0);
    a = gf16777216x_init(0);
    b = gf16777216x_init(0);
    g = gf16777216x_init(0);
    gf16777216x_xgcd(&a, &b, &g, x, y);

    ax = gf16777216x_init(0);
    by = gf16777216x_init(0);
    sum = gf16777216x_init(0);
    gf16777216x_multiply(&ax, a, x);
    gf16777216x_multiply(&by, b, y);
    gf16777216x_add(&sum, ax, by);
    equals = gf16777216x_equals(g, sum);

    gf16777216x_divide(&quotient, &remainder, x, g);
    equals &= gf16777216x_is_zero(remainder);

    gf16777216x_divide(&quotient, &remainder, y, g);
    equals &= gf16777216x_is_zero(remainder);

    if( equals == 1 )
    {
        printf("success!\n");
    }
    else
    {
        printf("fail!\n");
        printf("x: "); gf16777216x_print(x); printf("\n");
        printf("y: "); gf16777216x_print(y); printf("\n");
        printf("a: "); gf16777216x_print(a); printf("\n");
        printf("b: "); gf16777216x_print(b); printf("\n");
        printf("g: "); gf16777216x_print(g); printf("\n");

    gf16777216x_multiply(&ax, a, x);
    gf16777216x_multiply(&by, b, y);
    gf16777216x_add(&sum, ax, by);

        if( gf16777216x_equals(g, sum) == 1 )
        {
            printf("sum equals\n");
        }
        else
        {
            printf("sum is different\n");
            printf("ax: "); gf16777216x_print(ax); printf("\n");
            printf("by: "); gf16777216x_print(by); printf("\n");
            printf("sum: "); gf16777216x_print(sum); printf("\n");
            printf("g = "); gf16777216x_print(g); printf("\n");
        }

    gf16777216x_divide(&quotient, &remainder, x, g);
        if( gf16777216x_is_zero(remainder) == 1 )
        {
            printf("x is zero mod g\n");
        }
        else
        {
            printf("x is nonzero mod g\n");
        }

    gf16777216x_divide(&quotient, &remainder, y, g);
        if( gf16777216x_is_zero(remainder) == 1 )
        {
            printf("y is zero mod g\n");
        }
        else
        {
            printf("y is nonzero mod g\n");
        }
    }

    gf16777216x_destroy(x);
    gf16777216x_destroy(y);
    gf16777216x_destroy(a);
    gf16777216x_destroy(b);
    gf16777216x_destroy(g);
    gf16777216x_destroy(ax);
    gf16777216x_destroy(by);
    gf16777216x_destroy(sum);
    gf16777216x_destroy(remainder);
    gf16777216x_destroy(quotient);

    return equals;
}

int test_gf16777216x_modexp( )
{
    unsigned short int m, n, o;
    unsigned int i;
    int equal;
    gf16777216x g, p, ga, gb;
    csprng rng;
    unsigned long int a, b;
    unsigned int random;
    unsigned char * randomness;
    int equals;

    random = rand();
    csprng_init(&rng);
    csprng_seed(&rng, sizeof(unsigned int), (unsigned char *)&random);


    m = 10 + (csprng_generate_ulong(&rng) % 100);
    n = csprng_generate_ulong(&rng) % m;

    printf("testing modexp of GF(16777216)[x] elements for modulus degree %i ... ", m);

    p = gf16777216x_init(m);
    csprng_generate(&rng, 3*p.degree+3, p.data);

    g = gf16777216x_init(n);
    csprng_generate(&rng, 3*g.degree+3, g.data);

    ga = gf16777216x_init(0);
    gb = gf16777216x_init(0);

    a = csprng_generate_ulong(&rng) % (1 << (sizeof(unsigned long int)*3-1));
    b = csprng_generate_ulong(&rng) % (1 << (sizeof(unsigned long int)*3-1));

    gf16777216x_modexp(&ga, g, a, p);
    gf16777216x_modexp(&ga, ga, b, p);

    gf16777216x_modexp(&gb, g, a*b, p);

    equals = gf16777216x_equals(ga, gb);

    if( equals == 1 )
    {
        printf("success!\n");
    }
    else
    {
        printf("fail!\n");
        printf("ga: "); gf16777216x_print(ga); printf("\n");
        printf("gb: "); gf16777216x_print(gb); printf("\n");
    }

    gf16777216x_destroy(g);
    gf16777216x_destroy(p);
    gf16777216x_destroy(ga);
    gf16777216x_destroy(gb);

    return equals;
}

int test_gf2x_add( )
{
    unsigned short int m, n, o;
    unsigned int i;
    int equal;
    gf2x a, b, c, ab, abc1, bc, abc2;
    csprng rng;
    unsigned int random;
    unsigned char * randomness;
    int equals;

    random = rand();
    csprng_init(&rng);
    csprng_seed(&rng, sizeof(unsigned int), (unsigned char *)&random);


    m = csprng_generate_ulong(&rng) % 100;
    n = csprng_generate_ulong(&rng) % 100;
    o = csprng_generate_ulong(&rng) % 100;

    printf("testing addition of GF(2)[x] elements of degrees %i, %i, and %i... (randomness: %u) ", m, n, o, random);

    a = gf2x_init(m);
    csprng_generate(&rng, (a.degree+1+7)/8, a.data);

    b = gf2x_init(n);
    csprng_generate(&rng, (b.degree+1+7)/8, b.data);

    c = gf2x_init(o);
    csprng_generate(&rng, (c.degree+1+7)/8, c.data);

    ab = gf2x_init(0);
    gf2x_add(&ab, a, b);
    abc1 = gf2x_init(0);
    gf2x_add(&abc1, ab, c);

    bc = gf2x_init(0);
    gf2x_add(&bc, b, c);
    abc2 = gf2x_init(0);
    gf2x_add(&abc2, ab, c);

    equals = gf2x_equals(abc1, abc2);

    gf2x_destroy(a);
    gf2x_destroy(b);
    gf2x_destroy(c);
    gf2x_destroy(ab);
    gf2x_destroy(abc1);
    gf2x_destroy(bc);
    gf2x_destroy(abc2);

    if( equals == 1 )
    {
        printf("success!\n");
        return 1;
    }
    else
    {
        printf("fail!\n");
        return 0;
    }
}

int test_gf2x_multiply( )
{
    unsigned short int m, n, o;
    unsigned int i;
    int equal;
    gf2x a, b, c, ab, abc1, bc, abc2;
    csprng rng;
    unsigned int random;
    unsigned char * randomness;
    int equals;

    random = rand();
    csprng_init(&rng);
    csprng_seed(&rng, sizeof(unsigned int), (unsigned char *)&random);


    m = 10 + (csprng_generate_ulong(&rng) % 100);
    n = 10 + (csprng_generate_ulong(&rng) % 100);
    o = 10 + (csprng_generate_ulong(&rng) % 100);


    a = gf2x_init(m);
    csprng_generate(&rng, (a.degree+1+7)/8, a.data);
    gf2x_trim(&a);

    b = gf2x_init(n);
    csprng_generate(&rng, (b.degree+1+7)/8, b.data);
    gf2x_trim(&b);

    c = gf2x_init(o);
    csprng_generate(&rng, (c.degree+1+7)/8, c.data);
    gf2x_trim(&c);

    printf("testing multiplication of GF(2)[x] elements of degrees %i, %i, and %i... (randomness: %u) ", a.degree, b.degree, c.degree, random);

    ab = gf2x_init(0);
    gf2x_multiply(&ab, a, b);
    abc1 = gf2x_init(0);
    gf2x_multiply(&abc1, ab, c);

    bc = gf2x_init(0);
    gf2x_multiply(&bc, b, c);
    abc2 = gf2x_init(0);
    gf2x_multiply(&abc2, ab, c);

    equals = gf2x_equals(abc1, abc2);

    gf2x_destroy(a);
    gf2x_destroy(b);
    gf2x_destroy(c);
    gf2x_destroy(ab);
    gf2x_destroy(abc1);
    gf2x_destroy(bc);
    gf2x_destroy(abc2);

    if( equals == 1 )
    {
        printf("success!\n");
        return 1;
    }
    else
    {
        printf("fail!\n");
        return 0;
    }
}

int test_gf2x_divide( )
{
    unsigned short int m, n, o;
    unsigned int i;
    int equal;
    gf2x numerator, divisor, quotient, remainder, product, sum;
    csprng rng;
    unsigned int random;
    unsigned char * randomness;
    int equals;

    random = rand();
    csprng_init(&rng);
    csprng_seed(&rng, sizeof(unsigned int), (unsigned char *)&random);


    m = 10 + (csprng_generate_ulong(&rng) % 200);
    n = 10 + (csprng_generate_ulong(&rng) % (m-9));


    numerator = gf2x_init(m);
    csprng_generate(&rng, (numerator.degree+1+7)/8, numerator.data);
    gf2x_trim(&numerator);

    divisor = gf2x_init(n);
    csprng_generate(&rng, (divisor.degree+1+7)/8, divisor.data);
    gf2x_trim(&divisor);

    printf("testing division of GF(2)[x] elements of degrees %i, and %i... (randomness: %u) ", numerator.degree, divisor.degree, random);

    quotient = gf2x_init(0);
    remainder = gf2x_init(0);

    gf2x_divide(&quotient, &remainder, numerator, divisor);

    product = gf2x_init(0);
    gf2x_multiply(&product, divisor, quotient);

    sum = gf2x_init(0);
    gf2x_add(&sum, product, remainder);

    equals = gf2x_equals(numerator, sum);
    equals &= remainder.degree <= divisor.degree; /* includes division by constant */

    if( equals == 1 )
    {
        printf("success!\n");
    }
    else
    {
        printf("\n");
        printf("numerator: "); gf2x_print(numerator); printf("\n");
        printf("divisor: "); gf2x_print(divisor); printf("\n");
        printf("quotient: "); gf2x_print(quotient); printf("\n");
        printf("remainder: "); gf2x_print(remainder); printf("\n");
        printf("product: "); gf2x_print(product); printf("\n");
        printf("sum:     "); gf2x_print(sum); printf("\n");
        printf("fail!\n");
    }

    gf2x_destroy(numerator);
    gf2x_destroy(divisor);
    gf2x_destroy(remainder);
    gf2x_destroy(quotient);
    gf2x_destroy(product);
    gf2x_destroy(sum);

    return equals;
}

int test_gf2x_xgcd( )
{
    unsigned short int m, n;
    unsigned int i;
    int equal;
    gf2x x, y, a, b, g, ax, by, sum, quotient, remainder;
    csprng rng;
    unsigned int random;
    unsigned char * randomness;
    int equals;

    random = rand();
    csprng_init(&rng);
    csprng_seed(&rng, sizeof(unsigned int), (unsigned char *)&random);


    m = 10 + (csprng_generate_ulong(&rng) % 200);
    n = 10 + (csprng_generate_ulong(&rng) % 200);


    x = gf2x_init(m);
    csprng_generate(&rng, (x.degree+1+7)/8, x.data);
    gf2x_trim(&x);

    y = gf2x_init(n);
    csprng_generate(&rng, (y.degree+1+7)/8, y.data);
    gf2x_trim(&y);

    printf("testing xgcd of GF(2)[x] elements of degrees %i, and %i... (randomness: %u) ", x.degree, y.degree, random);

    quotient = gf2x_init(0);
    remainder = gf2x_init(0);
    a = gf2x_init(0);
    b = gf2x_init(0);
    g = gf2x_init(0);
    gf2x_xgcd(&a, &b, &g, x, y);

    ax = gf2x_init(0);
    by = gf2x_init(0);
    sum = gf2x_init(0);
    gf2x_multiply(&ax, a, x);
    gf2x_multiply(&by, b, y);
    gf2x_add(&sum, ax, by);
    equals = gf2x_equals(g, sum);

    gf2x_divide(&quotient, &remainder, x, g);
    equals &= gf2x_is_zero(remainder);

    gf2x_divide(&quotient, &remainder, y, g);
    equals &= gf2x_is_zero(remainder);

    if( equals == 1 )
    {
        printf("success!\n");
    }
    else
    {
        printf("fail!\n");
        printf("x: "); gf2x_print(x); printf("\n");
        printf("y: "); gf2x_print(y); printf("\n");
        printf("a: "); gf2x_print(a); printf("\n");
        printf("ax: "); gf2x_print(ax); printf("\n");
        printf("b: "); gf2x_print(b); printf("\n");
        printf("by: "); gf2x_print(by); printf("\n");
        printf("sum: "); gf2x_print(sum); printf("\n");
        printf("g =  "); gf2x_print(g); printf("\n");

        gf2x_divide(&quotient, &remainder, x, g);
        printf("x // g = "); gf2x_print(quotient); printf("\n");
        printf("x %% g = "); gf2x_print(remainder); printf("\n");
    }

    gf2x_destroy(x);
    gf2x_destroy(y);
    gf2x_destroy(a);
    gf2x_destroy(b);
    gf2x_destroy(g);
    gf2x_destroy(ax);
    gf2x_destroy(by);
    gf2x_destroy(sum);
    gf2x_destroy(remainder);
    gf2x_destroy(quotient);

    return equals;
}

int test_gf2x_lcm()
{
    unsigned short int m, n;
    unsigned int i;
    int equal;
    gf2x x, y, l, g, temp;
    csprng rng;
    unsigned int random;
    unsigned char * randomness;
    int success;

    random = rand();
    csprng_init(&rng);
    csprng_seed(&rng, sizeof(unsigned int), (unsigned char *)&random);


    m = 10 + (csprng_generate_ulong(&rng) % 200);
    n = 10 + (csprng_generate_ulong(&rng) % 200);


    x = gf2x_init(m);
    csprng_generate(&rng, (x.degree+1+7)/8, x.data);
    gf2x_trim(&x);

    y = gf2x_init(n);
    csprng_generate(&rng, (y.degree+1+7)/8, y.data);
    gf2x_trim(&y);

    printf("testing lcm of GF(2)[x] elements of degrees %i, and %i... (randomness: %u) ", x.degree, y.degree, random);


    l = gf2x_init(0);
    g = gf2x_init(0);

    gf2x_gcd(&g, x, y);
    gf2x_lcm(&l, x, y);

    temp = gf2x_init(0);
    gf2x_multiply(&temp, x, y);

    success = gf2x_is_one(g) || (!gf2x_equals(temp, l));

    if( success == 1 )
    {
        printf("success!\n");
    }
    else
    {
        printf("failure!\n");
        printf("x: "); gf2x_print(x); printf("\n");
        printf("y: "); gf2x_print(y); printf("\n");
        printf("l: "); gf2x_print(l); printf("\n");
        printf("g: "); gf2x_print(g); printf("\n");
        printf("p: "); gf2x_print(temp); printf("\n");
    }

    gf2x_destroy(x);
    gf2x_destroy(y);
    gf2x_destroy(g);
    gf2x_destroy(l);
    gf2x_destroy(temp);

    return success;
}

int test_gf2x_modinv()
{
    int m, n;
    unsigned int random;
    csprng rng;
    gf2x inv, elm, mod;
    int success;

    random = rand();
    csprng_init(&rng);
    csprng_seed(&rng, sizeof(unsigned int), (unsigned char *)&random);


    m = 10 + (csprng_generate_ulong(&rng) % 200);
    n = 9 + (csprng_generate_ulong(&rng) % (m-9));

    printf("testing gf2x_modinv on polynomials of degree at most %i and %i ... (randomness: %u) ", n, m, random);

    elm = gf2x_init(n);
    csprng_generate(&rng, n/8+1, elm.data);
    gf2x_trim(&elm);
    mod = gf2x_init(m);
    csprng_generate(&rng, m/8+1, mod.data);
    gf2x_trim(&mod);

    inv = gf2x_init(0);
    gf2x_zero(&inv);
    gf2x_gcd(&inv, elm, mod);

    while( elm.degree >= mod.degree || gf2x_is_one(inv) == 0 )
    {
        //printf("inside loop because elm.degree = %i >= mod.degree = %i or inv = "); gf2x_print(inv); printf(" =/= 1\n");
        elm.degree = n;
        csprng_generate(&rng, n/8+1, elm.data);
        gf2x_trim(&elm);
        mod.degree = m;
        csprng_generate(&rng, m/8+1, mod.data);
        gf2x_trim(&mod);
        gf2x_gcd(&inv, elm, mod);
    }

    gf2x_modinv(&inv, elm, mod);
    gf2x_multiply(&inv, inv, elm);
    gf2x_mod(&inv, inv, mod);

    success = gf2x_is_one(inv);

    if( success == 1 )
    {
        printf("success!\n");
    }
    else
    {
        printf("failure!\n");
        printf("elm: "); gf2x_print(elm); printf("\n");
        printf("mod: "); gf2x_print(mod); printf("\n");
        gf2x_modinv(&inv, elm, mod);
        printf("inv: "); gf2x_print(inv); printf("\n");
        gf2x_multiply(&inv, inv, elm);
        printf("pro: "); gf2x_print(inv); printf("\n");
        gf2x_mod(&inv, inv, mod);
        printf("rem: "); gf2x_print(inv); printf("\n");
    }

    gf2x_destroy(mod);
    gf2x_destroy(elm);
    gf2x_destroy(inv);

    return 1;
}

int test_gf2x_modexp()
{
    int m, n;
    long int a, b;
    unsigned int random;
    csprng rng;
    gf2x A, B, inv, elm, mod;
    int success;

    random = rand();
    csprng_init(&rng);
    csprng_seed(&rng, sizeof(unsigned int), (unsigned char *)&random);


    m = 10 + (csprng_generate_ulong(&rng) % 200);
    n = 9 + (csprng_generate_ulong(&rng) % (m-9));

    printf("testing gf2x_modexp on polynomials of degree at most %i and %i ... (randomness: %u) ", n, m, random);

    elm = gf2x_init(n);
    csprng_generate(&rng, n/8+1, elm.data);
    gf2x_trim(&elm);
    mod = gf2x_init(m);
    csprng_generate(&rng, m/8+1, mod.data);
    gf2x_trim(&mod);

    inv = gf2x_init(0);
    gf2x_zero(&inv);
    gf2x_gcd(&inv, elm, mod);

    while( elm.degree >= mod.degree || gf2x_is_one(inv) == 0 )
    {
        //printf("inside loop because elm.degree = %i >= mod.degree = %i or inv = "); gf2x_print(inv); printf(" =/= 1\n");
        elm.degree = n;
        csprng_generate(&rng, n/8+1, elm.data);
        gf2x_trim(&elm);
        mod.degree = m;
        csprng_generate(&rng, m/8+1, mod.data);
        gf2x_trim(&mod);
        gf2x_gcd(&inv, elm, mod);
    }

    a = csprng_generate_ulong(&rng) >> 2;
    b = csprng_generate_ulong(&rng) >> 2;
    
    if( csprng_generate_ulong(&rng) % 2 == 1 )
    {
        a = -a;
    }
    if( csprng_generate_ulong(&rng) % 2 == 1 )
    {
        b = -b;
    }

    A = gf2x_init(0);
    B = gf2x_init(0);

    gf2x_modexp(&A, elm, a, mod);
    gf2x_modexp(&A, A, b, mod);

    gf2x_modexp(&B, elm, b, mod);
    gf2x_modexp(&B, B, a, mod);

    success = gf2x_equals(A, B);

    if( success == 1 )
    {
        printf("success!\n");
    }
    else
    {
        printf("failure!\n");
    }

    gf2x_destroy(inv);
    gf2x_destroy(elm);
    gf2x_destroy(mod);
    gf2x_destroy(A);
    gf2x_destroy(B);

    return success;
}

int test_gf2x_karatsuba( )
{
    unsigned short int m, n, o;
    unsigned int i;
    int equal;
    gf2x a, b, c, ab, abc1, bc, abc2;
    csprng rng;
    unsigned int random;
    unsigned char * randomness;
    int equals;

    random = rand();
    csprng_init(&rng);
    csprng_seed(&rng, sizeof(unsigned int), (unsigned char *)&random);


    m = 10 + (csprng_generate_ulong(&rng) % 10000);
    n = 10 + (csprng_generate_ulong(&rng) % 10000);
    o = 10 + (csprng_generate_ulong(&rng) % 10000);


    a = gf2x_init(m);
    csprng_generate(&rng, (a.degree+1+7)/8, a.data);
    gf2x_trim(&a);

    b = gf2x_init(n);
    csprng_generate(&rng, (b.degree+1+7)/8, b.data);
    gf2x_trim(&b);

    c = gf2x_init(o);
    csprng_generate(&rng, (c.degree+1+7)/8, c.data);
    gf2x_trim(&c);

    printf("testing karatsuba of GF(2)[x] elements of degrees %i, %i, and %i... (randomness: %u) ", a.degree, b.degree, c.degree, random);

    ab = gf2x_init(0);
    gf2x_karatsuba(&ab, a, b);
    abc1 = gf2x_init(0);
    gf2x_karatsuba(&abc1, ab, c);

    bc = gf2x_init(0);
    gf2x_karatsuba(&bc, b, c);
    abc2 = gf2x_init(0);
    gf2x_karatsuba(&abc2, ab, c);

    equals = gf2x_equals(abc1, abc2);
    equals &= gf2x_divides(a, abc1);
    equals &= gf2x_divides(b, abc1);
    equals &= gf2x_divides(c, abc1);
    equals &= gf2x_divides(a, abc2);
    equals &= gf2x_divides(b, abc2);
    equals &= gf2x_divides(c, abc2);
    if( equals == 1 )
    {
        printf("success!\n");
    }
    else
    {
        printf("fail!\n");
        printf("a: "); gf2x_print(a); printf(" (%i) \n", a.degree);
        printf("b: "); gf2x_print(b); printf(" (%i) \n", b.degree);
        printf("c: "); gf2x_print(c); printf(" (%i) \n", c.degree);
        printf("a*b: "); gf2x_print(ab); printf(" (%i) \n", ab.degree);
        printf("b*c: "); gf2x_print(bc); printf(" (%i) \n", bc.degree);
        printf("a*b*c: "); gf2x_print(abc1); printf(" (%i) \n", abc1.degree);
        printf("a*b*c: "); gf2x_print(abc2); printf(" (%i) \n", abc2.degree);
    }

    gf2x_destroy(a);
    gf2x_destroy(b);
    gf2x_destroy(c);
    gf2x_destroy(ab);
    gf2x_destroy(abc1);
    gf2x_destroy(bc);
    gf2x_destroy(abc2);

    return equals;
}

int test_gf2x_minpoly()
{
    int m, n;
    long int a, b;
    unsigned int random;
    csprng rng;
    gf2x min, elm, mod, ev, raised;
    int success;
    int i;

    random = rand();
    csprng_init(&rng);
    csprng_seed(&rng, sizeof(unsigned int), (unsigned char *)&random);


    m = 10 + (csprng_generate_ulong(&rng) % 200);
    n = 9 + (csprng_generate_ulong(&rng) % (m-9));

    printf("testing gf2x_minpoly on polynomials of degree at most %i and %i ... (randomness: %u) ", n, m, random);

    elm = gf2x_init(n);
    csprng_generate(&rng, n/8+1, elm.data);
    gf2x_trim(&elm);
    mod = gf2x_init(m);
    csprng_generate(&rng, m/8+1, mod.data);
    gf2x_trim(&mod);

    min = gf2x_init(0);
    gf2x_zero(&min);
    gf2x_gcd(&min, elm, mod);

    while( elm.degree >= mod.degree || gf2x_is_one(min) == 0 )
    {
        //printf("inside loop because elm.degree = %i >= mod.degree = %i or inv = "); gf2x_print(inv); printf(" =/= 1\n");
        elm.degree = n;
        csprng_generate(&rng, n/8+1, elm.data);
        gf2x_trim(&elm);
        mod.degree = m;
        csprng_generate(&rng, m/8+1, mod.data);
        gf2x_trim(&mod);
        gf2x_gcd(&min, elm, mod);
    }

    gf2x_minpoly(&min, elm, mod);

    /* test if minpoly evaluates to zero in elm */
    ev = gf2x_init(0);
    raised = gf2x_init(0);
    gf2x_one(&raised);
    gf2x_zero(&ev);
    for( i = 0 ; i <= min.degree ; ++i )
    {
        if( (min.data[i/8] & (1 << (i % 8))) != 0 )
        {
            gf2x_add(&ev, ev, raised);
        }
        gf2x_multiply(&raised, raised, elm);
        gf2x_mod(&raised, raised, mod);
    }

    success = gf2x_is_zero(ev);

    if( success == 1 )
    {
        printf("success!\n");
    }
    else
    {
        printf("failure!\n");
        gf2x_print(min); printf("("); gf2x_print(elm); printf(") = "); gf2x_print(ev); printf(" mod "); gf2x_print(mod); printf("\n");
    }

    gf2x_destroy(raised);
    gf2x_destroy(ev);
    gf2x_destroy(min);
    gf2x_destroy(elm);
    gf2x_destroy(mod);

    return success;
}


int test_gf4096_inverse( )
{
    csprng rng;
    unsigned int random;
    int equals;
    unsigned int a, b, c;

    random = rand();
    csprng_init(&rng);
    csprng_seed(&rng, sizeof(unsigned int), (unsigned char *)&random);


    printf("testing inverse computation of GF(4096) elements ... ");

    a = csprng_generate_ulong(&rng) % 0xffff;

    b = gf4096_inverse(a);

    c = gf4096_multiply(a, b);

    if( c == 1 || a == 0 )
    {
        printf("success!\n");
        return 1;
    }
    else
    {
        printf("fail!\n");
        return 0;
    }
}

int test_gf4096x_add( )
{
    unsigned short int m, n, o;
    unsigned int i;
    int equal;
    gf4096x a, b, c, ab, abc1, bc, abc2;
    csprng rng;
    unsigned int random;
    unsigned char * randomness;
    int equals;

    random = rand();
    csprng_init(&rng);
    csprng_seed(&rng, sizeof(unsigned int), (unsigned char *)&random);


    m = csprng_generate_ulong(&rng) % 100;
    n = csprng_generate_ulong(&rng) % 100;
    o = csprng_generate_ulong(&rng) % 100;

    printf("testing addition of GF(4096)[x] elements of degrees %i, %i, and %i... ", m, n, o);

    a = gf4096x_init(m);
    csprng_generate(&rng, 2*a.degree+2, a.data);

    b = gf4096x_init(n);
    csprng_generate(&rng, 2*b.degree+2, b.data);

    c = gf4096x_init(o);
    csprng_generate(&rng, 2*c.degree+2, c.data);

    ab = gf4096x_init(0);
    gf4096x_add(&ab, a, b);
    abc1 = gf4096x_init(0);
    gf4096x_add(&abc1, ab, c);

    bc = gf4096x_init(0);
    gf4096x_add(&bc, b, c);
    abc2 = gf4096x_init(0);
    gf4096x_add(&abc2, ab, c);

    equals = gf4096x_equals(abc1, abc2);

    gf4096x_destroy(a);
    gf4096x_destroy(b);
    gf4096x_destroy(c);
    gf4096x_destroy(ab);
    gf4096x_destroy(abc1);
    gf4096x_destroy(bc);
    gf4096x_destroy(abc2);

    if( equals == 1 )
    {
        printf("success!\n");
        return 1;
    }
    else
    {
        printf("fail!\n");
        return 0;
    }
}

int test_gf4096x_multiply( )
{
    unsigned short int m, n, o;
    unsigned int i;
    int equal;
    gf4096x a, b, c, ab, abc1, bc, abc2;
    csprng rng;
    unsigned int random;
    unsigned char * randomness;
    int equals;

    random = rand();
    csprng_init(&rng);
    csprng_seed(&rng, sizeof(unsigned int), (unsigned char *)&random);


    m = csprng_generate_ulong(&rng) % 100;
    n = csprng_generate_ulong(&rng) % 100;
    o = csprng_generate_ulong(&rng) % 100;

    printf("testing multiplication of GF(4096)[x] elements of degrees %i, %i, and %i... ", m, n, o);

    a = gf4096x_init(m);
    csprng_generate(&rng, 2*a.degree+2, a.data);

    b = gf4096x_init(n);
    csprng_generate(&rng, 2*b.degree+2, b.data);

    c = gf4096x_init(o);
    csprng_generate(&rng, 2*c.degree+2, c.data);

    ab = gf4096x_init(0);
    gf4096x_multiply(&ab, a, b);
    abc1 = gf4096x_init(0);
    gf4096x_multiply(&abc1, ab, c);

    bc = gf4096x_init(0);
    gf4096x_multiply(&bc, b, c);
    abc2 = gf4096x_init(0);
    gf4096x_multiply(&abc2, ab, c);

    equals = gf4096x_equals(abc1, abc2);
    equals &= (a.degree + b.degree + c.degree) == abc1.degree;

    gf4096x_destroy(a);
    gf4096x_destroy(b);
    gf4096x_destroy(c);
    gf4096x_destroy(ab);
    gf4096x_destroy(abc1);
    gf4096x_destroy(bc);
    gf4096x_destroy(abc2);

    if( equals == 1 )
    {
        printf("success!\n");
        return 1;
    }
    else
    {
        printf("fail!\n");
        return 0;
    }
}

int test_gf4096x_divide( )
{
    unsigned short int m, n, o;
    unsigned int i;
    int equal;
    gf4096x numerator, divisor, quotient, remainder, product, sum;
    csprng rng;
    unsigned int random;
    unsigned char * randomness;
    int equals;

    random = rand();
    csprng_init(&rng);
    csprng_seed(&rng, sizeof(unsigned int), (unsigned char *)&random);


    m = csprng_generate_ulong(&rng) % 200;
    n = csprng_generate_ulong(&rng) % m;

    printf("testing division of GF(4096)[x] elements of degrees %i, and %i... (randomness: %u) ", m, n, random);

    numerator = gf4096x_init(m);
    csprng_generate(&rng, 2*numerator.degree+2, numerator.data);

    divisor = gf4096x_init(n);
    csprng_generate(&rng, 2*divisor.degree+2, divisor.data);

    quotient = gf4096x_init(0);
    remainder = gf4096x_init(0);
    gf4096x_divide(&quotient, &remainder, numerator, divisor);

    product = gf4096x_init(0);
    gf4096x_multiply(&product, divisor, quotient);

    sum = gf4096x_init(0);
    gf4096x_add(&sum, product, remainder);

    equals = gf4096x_equals(numerator, sum);
    equals &= remainder.degree <= divisor.degree; /* includes division by constant */

    if( equals == 1 )
    {
        printf("success!\n");
    }
    else
    {
        printf("\n");
        printf("numerator: "); gf4096x_print(numerator); printf("\n");
        printf("(numerator degree: %i)\n", numerator.degree);
        printf("divisor: "); gf4096x_print(divisor); printf("\n");
        printf("(divisor degree: %i)\n", divisor.degree);
        printf("quotient: "); gf4096x_print(quotient); printf("\n");
        printf("remainder: "); gf4096x_print(remainder); printf("\n");
        printf("(remainder degree: %i)\n", remainder.degree);
        printf("product: "); gf4096x_print(product); printf("\n");
        printf("sum: "); gf4096x_print(sum); printf("\n");
        printf("(sum degree: %i)\n", sum.degree);
        gf4096x_add(&sum, sum, numerator);
        printf("difference: "); gf4096x_print(sum); printf("\n");
        printf("(difference degree: %i)\n", sum.degree);
        printf("fail!\n");
    }

    gf4096x_destroy(numerator);
    gf4096x_destroy(divisor);
    gf4096x_destroy(remainder);
    gf4096x_destroy(quotient);
    gf4096x_destroy(product);
    gf4096x_destroy(sum);

    return equals;
}

int test_gf4096x_xgcd( )
{
    unsigned short int m, n;
    unsigned int i;
    int equal;
    gf4096x x, y, a, b, g, ax, by, sum, quotient, remainder;
    csprng rng;
    unsigned int random;
    unsigned char * randomness;
    int equals;

    random = rand();
    csprng_init(&rng);
    csprng_seed(&rng, sizeof(unsigned int), (unsigned char *)&random);


    m = csprng_generate_ulong(&rng) % 200;
    n = csprng_generate_ulong(&rng) % 200;

    printf("testing xgcd GF(4096)[x] elements of degrees %i, and %i... ", m, n);

    x = gf4096x_init(m);
    csprng_generate(&rng, 2*x.degree+2, x.data);

    y = gf4096x_init(n);
    csprng_generate(&rng, 2*y.degree+2, y.data);

    quotient = gf4096x_init(0);
    remainder = gf4096x_init(0);
    a = gf4096x_init(0);
    b = gf4096x_init(0);
    g = gf4096x_init(0);
    gf4096x_xgcd(&a, &b, &g, x, y);

    ax = gf4096x_init(0);
    by = gf4096x_init(0);
    sum = gf4096x_init(0);
    gf4096x_multiply(&ax, a, x);
    gf4096x_multiply(&by, b, y);
    gf4096x_add(&sum, ax, by);
    equals = gf4096x_equals(g, sum);

    gf4096x_divide(&quotient, &remainder, x, g);
    equals &= gf4096x_is_zero(remainder);

    gf4096x_divide(&quotient, &remainder, y, g);
    equals &= gf4096x_is_zero(remainder);

    if( equals == 1 )
    {
        printf("success!\n");
    }
    else
    {
        printf("fail!\n");
        printf("x: "); gf4096x_print(x); printf("\n");
        printf("y: "); gf4096x_print(y); printf("\n");
        printf("a: "); gf4096x_print(a); printf("\n");
        printf("b: "); gf4096x_print(b); printf("\n");
        printf("g: "); gf4096x_print(g); printf("\n");

    gf4096x_multiply(&ax, a, x);
    gf4096x_multiply(&by, b, y);
    gf4096x_add(&sum, ax, by);

        if( gf4096x_equals(g, sum) == 1 )
        {
            printf("sum equals\n");
        }
        else
        {
            printf("sum is different\n");
            printf("ax: "); gf4096x_print(ax); printf("\n");
            printf("by: "); gf4096x_print(by); printf("\n");
            printf("sum: "); gf4096x_print(sum); printf("\n");
            printf("g = "); gf4096x_print(g); printf("\n");
        }

    gf4096x_divide(&quotient, &remainder, x, g);
        if( gf4096x_is_zero(remainder) == 1 )
        {
            printf("x is zero mod g\n");
        }
        else
        {
            printf("x is nonzero mod g\n");
        }

    gf4096x_divide(&quotient, &remainder, y, g);
        if( gf4096x_is_zero(remainder) == 1 )
        {
            printf("y is zero mod g\n");
        }
        else
        {
            printf("y is nonzero mod g\n");
        }
    }

    gf4096x_destroy(x);
    gf4096x_destroy(y);
    gf4096x_destroy(a);
    gf4096x_destroy(b);
    gf4096x_destroy(g);
    gf4096x_destroy(ax);
    gf4096x_destroy(by);
    gf4096x_destroy(sum);
    gf4096x_destroy(remainder);
    gf4096x_destroy(quotient);

    return equals;
}




int main( int argc, char ** argv )
{
    unsigned int i;
    unsigned int b;
#ifdef BIG
    bi ninetythree;
    ninetythree = bi_cast(93);
    prime_modulus = bi_cast(1);
    bi_shift_left(&prime_modulus, prime_modulus, 257);
    bi_subtract(&prime_modulus, prime_modulus, ninetythree);
#endif
    b = 1;

    printf("testing basic algebra routines for finite fields ...\n");

    for( i = 0 ; i < 0 && b == 1 ; ++i ) b = b & test_csprng();
    for( i = 0 ; i < 0 && b == 1 ; ++i ) b = b & test_matrix_inverse();
    for( i = 0 ; i < 0 && b == 1 ; ++i ) b = b & test_multiply_transpose();
    for( i = 0 ; i < 0 && b == 1 ; ++i ) b = b & test_solve();
    for( i = 0 ; i < 0 && b == 1 ; ++i ) b = b & test_composition();
    for( i = 0 ; i < 0 && b == 1 ; ++i ) b = b & test_gf256x_add();
    for( i = 0 ; i < 0 && b == 1 ; ++i ) b = b & test_gf256x_multiply();
    for( i = 0 ; i < 0 && b == 1 ; ++i ) b = b & test_gf256x_divide();
    for( i = 0 ; i < 0 && b == 1 ; ++i ) b = b & test_gf256x_xgcd();
    for( i = 0 ; i < 0 && b == 1 ; ++i ) b = b & test_gf65536_inverse();
    for( i = 0 ; i < 0 && b == 1 ; ++i ) b = b & test_gf65536x_add();
    for( i = 0 ; i < 0 && b == 1 ; ++i ) b = b & test_gf65536x_multiply();
    for( i = 0 ; i < 0 && b == 1 ; ++i ) b = b & test_gf65536x_divide();
    for( i = 0 ; i < 0 && b == 1 ; ++i ) b = b & test_gf65536x_xgcd();
    for( i = 0 ; i < 0 && b == 1 ; ++i ) b = b & test_gf16777216_inverse();
    for( i = 0 ; i < 0 && b == 1 ; ++i ) b = b & test_gf16777216x_add();
    for( i = 0 ; i < 0 && b == 1 ; ++i ) b = b & test_gf16777216x_multiply();
    for( i = 0 ; i < 0 && b == 1 ; ++i ) b = b & test_gf16777216x_divide();
    for( i = 0 ; i < 0 && b == 1 ; ++i ) b = b & test_gf16777216x_xgcd();
    for( i = 0 ; i < 0 && b == 1 ; ++i ) b = b & test_gf16777216x_modexp();
    for( i = 0 ; i < 0 && b == 1 ; ++i ) b = b & test_gf2x_add();
    for( i = 0 ; i < 0 && b == 1 ; ++i ) b = b & test_gf2x_multiply();
    for( i = 0 ; i < 0 && b == 1 ; ++i ) b = b & test_gf2x_divide();
    for( i = 0 ; i < 0 && b == 1 ; ++i ) b = b & test_gf2x_xgcd();
    for( i = 0 ; i < 0 && b == 1 ; ++i ) b = b & test_gf2x_lcm();
    for( i = 0 ; i < 0 && b == 1 ; ++i ) b = b & test_gf2x_modinv();
    for( i = 0 ; i < 0 && b == 1 ; ++i ) b = b & test_gf2x_modexp();
    for( i = 0 ; i < 0 && b == 1 ; ++i ) b = b & test_gf2x_karatsuba();
    for( i = 0 ; i < 0 && b == 1 ; ++i ) b = b & test_gf2x_minpoly();
    for( i = 0 ; i < 10 && b == 1 ; ++i ) b = b & test_gf4096_inverse();
    for( i = 0 ; i < 10 && b == 1 ; ++i ) b = b & test_gf4096x_add();
    for( i = 0 ; i < 10 && b == 1 ; ++i ) b = b & test_gf4096x_multiply();
    for( i = 0 ; i < 10 && b == 1 ; ++i ) b = b & test_gf4096x_divide();
    for( i = 0 ; i < 10 && b == 1 ; ++i ) b = b & test_gf4096x_xgcd();

#ifdef BIG
    bi_destroy(ninetythree);
    bi_destroy(prime_modulus);
#endif

    if( b == 1 )
    {
        printf("success.\n");
        return 1;
    }
    else
    {
        printf("failure.\n");
        return 0;
    }
}

