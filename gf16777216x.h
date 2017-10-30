#ifndef GF16777216X_H
#define GF16777216X_H

unsigned int gf16777216_multiply( unsigned int lhs, unsigned int rhs );
unsigned int gf16777216_inverse( unsigned int elm );
unsigned int gf16777216_exp( unsigned int elm, int exponent );

typedef struct
{
    unsigned char * data;
    int degree;
} gf16777216x;

gf16777216x gf16777216x_init( int deg );
int gf16777216x_copy( gf16777216x* dest, gf16777216x source );
int gf16777216x_destroy( gf16777216x p );

unsigned int gf16777216x_eval( gf16777216x polynomial, unsigned int point );
int gf16777216x_add( gf16777216x* dest, gf16777216x lhs, gf16777216x rhs );
int gf16777216x_multiply( gf16777216x* dest, gf16777216x lhs, gf16777216x rhs );
int gf16777216x_multiply_constant_shift( gf16777216x* dest, gf16777216x poly, unsigned int constant, int shift );
int gf16777216x_divide( gf16777216x* quo, gf16777216x* rem, gf16777216x num, gf16777216x den );
int gf16777216x_xgcd( gf16777216x* a, gf16777216x* b, gf16777216x* g, gf16777216x x, gf16777216x y );
int gf16777216x_modexp( gf16777216x* dest, gf16777216x base, unsigned long int exponent, gf16777216x modulus );

int gf16777216x_equals( gf16777216x lhs, gf16777216x rhs );
int gf16777216x_is_zero( gf16777216x lhs );

int gf16777216x_print( gf16777216x p );

#endif

