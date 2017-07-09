#include "gfp.h"
#include "gfpm.h"
#include "hqs.h"
#include "csprng.h"

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

    m = 15;
    n = 10;

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

    m = 20;
    n = 21;

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


    m = 4;
    n = 6;

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

