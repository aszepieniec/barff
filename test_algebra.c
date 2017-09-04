#include "gfp.h"
#include "gfpm.h"
#include "hqs.h"
#include "csprng.h"
#include "gf256x.h"

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
    equals &= (a.degree + b.degree + c.degree) == abc1.degree;

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

    for( i = 0 ; i < 10 && b == 1 ; ++i ) b = b & test_csprng();
    for( i = 0 ; i < 10 && b == 1 ; ++i ) b = b & test_matrix_inverse();
    for( i = 0 ; i < 10 && b == 1 ; ++i ) b = b & test_multiply_transpose();
    for( i = 0 ; i < 10 && b == 1 ; ++i ) b = b & test_solve();
    for( i = 0 ; i < 10 && b == 1 ; ++i ) b = b & test_composition();
    for( i = 0 ; i < 10 && b == 1 ; ++i ) b = b & test_gf256x_add();
    for( i = 0 ; i < 10 && b == 1 ; ++i ) b = b & test_gf256x_multiply();
    for( i = 0 ; i < 10 && b == 1 ; ++i ) b = b & test_gf256x_divide();
    for( i = 0 ; i < 10 && b == 1 ; ++i ) b = b & test_gf256x_xgcd();

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

