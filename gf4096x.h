#ifndef GF4096X_H
#define GF4096X_H

unsigned int gf4096_multiply( unsigned int lhs, unsigned int rhs );
unsigned int gf4096_inverse( unsigned int elm );
unsigned int gf4096_exp( unsigned int elm, int exponent );

typedef struct
{
    unsigned char * data;
    /* a single GF(2^12) element takes up 2 bytes; the 4 msb's
     * are ignored but for good practice should be set to zero
     * always. */
    int degree;
} gf4096x;

gf4096x gf4096x_init( int deg );
int gf4096x_copy( gf4096x* dest, gf4096x source );
int gf4096x_destroy( gf4096x p );

unsigned int gf4096x_eval( gf4096x polynomial, unsigned int point );
int gf4096x_add( gf4096x* dest, gf4096x lhs, gf4096x rhs );
int gf4096x_multiply( gf4096x* dest, gf4096x lhs, gf4096x rhs );
int gf4096x_multiply_constant_shift( gf4096x* dest, gf4096x poly, unsigned int constant, int shift );
int gf4096x_divide( gf4096x* quo, gf4096x* rem, gf4096x num, gf4096x den );
int gf4096x_xgcd( gf4096x* a, gf4096x* b, gf4096x* g, gf4096x x, gf4096x y );

int gf4096x_print( gf4096x p );

#endif

