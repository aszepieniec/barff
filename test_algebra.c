#include "algebra.h"

#include <stdio.h>

int main( int argc, char ** argv )
{
    FILE * fh;
    unsigned char buf1[10*10];
    unsigned char buf2[10*10];
    unsigned char buf3[10*10];
    gfmatrix A = gfm(10, 10, buf1), B=gfm(10,10,buf2), C=gfm(10,10,buf3);

    unsigned char randomness[10*10+1];
    fh = fopen("/dev/urandom", "rb");

    fgets(randomness, sizeof(randomness), fh);
    gfm_random(A, randomness);

    fclose(fh);

    gfm_inverse(B, A);

    gfm_multiply(C, A, B);

    printf("A:\n");
    gfm_print(A);
    printf("B:\n");
    gfm_print(B);
    printf("C:\n");
    gfm_print(C);

    return 1;
};

