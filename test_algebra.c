#include "algebra.h"
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

int test_matrix_inverse()
{
    FILE * fh;
    csprng rng;
    unsigned char buf1[10*10];
    unsigned char buf2[10*10];
    unsigned char buf3[10*10];

    unsigned char randomness[10*10+1];
    char * prand;
    int invertible;

    gfmatrix A, B, C;

    unsigned int random;
    random = rand();
    csprng_init(&rng);
    csprng_seed(&rng, sizeof(unsigned int), (unsigned char*)&random);

    A = gfm(10, 10, buf1);
    B = gfm(10,10,buf2);
    C = gfm(10,10,buf3);


    printf("testing matrix inverse ... ");

    gfm_random(A, &rng);

    invertible = gfm_inverse(B, A);

    gfm_multiply(C, A, B);

    if( invertible == 0 || gfm_is_eye(C) )
    {
        printf("success!\n");
        return 1;
    }
    else
    {
        printf("fail!\n");
        printf("A:\n");
        gfm_print(A);
        printf("B:\n");
        gfm_print(B);
        printf("C:\n");
        gfm_print(C);
        printf("random seed: %u\n", random);
        return 0;
    }
}

int test_composition()
{
    hqsystem F, P;
    gfmatrix T, S;
    unsigned short int m, n;
    gfmatrix x, Sx, FoSx, ToFoSx, y;
    unsigned int i;
    int equal;
    csprng rng;
    unsigned int random;

    random = rand();
    csprng_init(&rng);
    csprng_seed(&rng, sizeof(unsigned int), (unsigned char *)&random);


    printf("testing composition of linear transforms with homogeneous quadratic systems ... ");


    m = 10;
    n = 15;

    F = hqs_init(n, m);
    hqs_random(F, &rng);
    P = hqs_copy_new(F);



    T = gfm_init(m, m);
    S = gfm_init(n, n);


    gfm_random_invertible(T, &rng);
    gfm_random_invertible(S, &rng);


    hqs_compose_output(T, P);
    hqs_compose_input(P, S);


    x = gfm_init(n, 1);
    y = gfm_init(m, 1);
    Sx = gfm_init(n, 1);
    FoSx = gfm_init(m, 1);
    ToFoSx = gfm_init(m, 1);


    equal = 1;
    for( i = 0 ; i < 1 && equal == 1 ; ++i )
    {
        gfm_random(x, &rng);

        gfm_multiply(Sx, S, x);

        hqs_eval(FoSx, F, Sx);

        gfm_multiply(ToFoSx, T, FoSx);

        hqs_eval(y, P, x);

        equal = gfm_equals(y, ToFoSx);
        getchar();
    }

    gfm_destroy(x);
    gfm_destroy(y);
    gfm_destroy(Sx);
    gfm_destroy(FoSx);
    gfm_destroy(ToFoSx);
    hqs_destroy(F);
    hqs_destroy(P);
    gfm_destroy(T);
    gfm_destroy(S);

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
    b = 1;

    printf("testing basic algebra routines for finite fields ...\n");

    for( i = 0 ; i < 10 && b == 1 ; ++i ) b = b & test_csprng();
    for( i = 0 ; i < 10 && b == 1 ; ++i ) b = b & test_matrix_inverse();
    for( i = 0 ; i < 10 && b == 1 ; ++i ) b = b & test_composition();

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

