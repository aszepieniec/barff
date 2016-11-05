#include "algebra.h"

#include <stdio.h>
#include <stdlib.h>

int test_matrix_inverse()
{
    FILE * fh;
    unsigned char buf1[10*10];
    unsigned char buf2[10*10];
    unsigned char buf3[10*10];

    unsigned char randomness[10*10+1];
    char * prand;
    int invertible;

    gfmatrix A, B, C;

    A = gfm(10, 10, buf1);
    B = gfm(10,10,buf2);
    C = gfm(10,10,buf3);


    printf("testing matrix inverse ... ");
    fh = fopen("/dev/urandom", "rb");

    fread(randomness, 1, sizeof(randomness), fh);
    gfm_random(A, randomness);

    fclose(fh);

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
    unsigned char * buf1;
    unsigned char * buf2;
    unsigned char * buf3;
    unsigned char * buf4;
    FILE * fh;

    printf("testing composition of linear transforms with homogeneous quadratic systems ... ");

    m = 10;
    n = 15;

    fh = fopen("/dev/urandom", "rb");

    buf1 = malloc(n+1);
    buf2 = malloc(n*n*m+1);
    buf3 = malloc(n*n);
    buf4 = malloc(m*m);

    fread(buf2, 1, m*n*n, fh);

    F = hqs_init(n, m);
    hqs_random(F, buf2);
    P = hqs_copy_new(F);


    fread(buf3, 1, n*n, fh);
    fread(buf4, 1, m*m, fh);
    T = gfm_init(m, m);
    S = gfm_init(n, n);


    gfm_random_invertible(T, buf3);
    gfm_random_invertible(S, buf4);

    hqs_compose_output(T, P);
    hqs_compose_input(P, S);


    x = gfm_init(n, 1);
    y = gfm_init(m, 1);
    Sx = gfm_init(n, 1);
    FoSx = gfm_init(m, 1);
    ToFoSx = gfm_init(m, 1);

    equal = 1;
    for( i = 0 ; i < 100 && equal == 1 ; ++i )
    {
        fread(buf1, 1, n, fh);
        gfm_random(x, buf1);

        gfm_multiply(Sx, S, x);

        hqs_eval(FoSx, F, Sx);

        gfm_multiply(ToFoSx, T, FoSx);

        hqs_eval(y, P, x);

        if( gfm_equals(y, ToFoSx) == 0 )
        {
            equal = 0;
        }
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
    free(buf1);
    free(buf2);

    if( equal )
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

int main( int argc, char ** argv )
{
    unsigned int i;
    unsigned int b;
    b = 1;

    printf("testing basic algebra routines for finite fields ...\n");

    for( i = 0 ; i < 100 && b == 1 ; ++i ) b = b & test_matrix_inverse();
    for( i = 0 ; i < 100 && b == 1 ; ++i ) b = b & test_composition();

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

