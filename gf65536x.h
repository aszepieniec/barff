#ifndef GF65536X_H
#define GF65536X_H

unsigned int gf65536_multiply( unsigned int lhs, unsigned int rhs );
unsigned int gf65536_inverse( unsigned int elm );
unsigned int gf65536_exp( unsigned int elm, int exponent );

typedef struct
{
    unsigned char * data;
    int degree;
} gf65536x;

gf65536x gf65536x_init( int deg );
int gf65536x_copy( gf65536x* dest, gf65536x source );
int gf65536x_destroy( gf65536x p );

unsigned int gf65536x_eval( gf65536x polynomial, unsigned int point );
int gf65536x_add( gf65536x* dest, gf65536x lhs, gf65536x rhs );
int gf65536x_multiply( gf65536x* dest, gf65536x lhs, gf65536x rhs );
int gf65536x_multiply_constant_shift( gf65536x* dest, gf65536x poly, unsigned int constant, int shift );
int gf65536x_divide( gf65536x* quo, gf65536x* rem, gf65536x num, gf65536x den );
int gf65536x_xgcd( gf65536x* a, gf65536x* b, gf65536x* g, gf65536x x, gf65536x y );

int gf65536x_print( gf65536x p );

#endif

