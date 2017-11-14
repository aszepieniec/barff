#ifndef GF2X_H
#define GF2X_H

typedef struct
{
    unsigned char * data;
    int degree;
} gf2x;

gf2x gf2x_init( int deg );
int gf2x_copy( gf2x* dest, gf2x source );
int gf2x_destroy( gf2x p );

int gf2x_is_zero( gf2x poly );
int gf2x_is_one( gf2x poly );

int gf2x_trim( gf2x* poly );
int gf2x_add( gf2x* dest, gf2x lhs, gf2x rhs );
int gf2x_shift_left( gf2x* dest, gf2x poly, unsigned int shift );
int gf2x_multiply( gf2x* dest, gf2x lhs, gf2x rhs );
#define GF2X_KARATSUBA_CUTOFF 320
int gf2x_karatsuba( gf2x* dest, gf2x lhs, gf2x rhs );
int gf2x_divide( gf2x* quo, gf2x* rem, gf2x num, gf2x den );
int gf2x_xgcd( gf2x* a, gf2x* b, gf2x* g, gf2x x, gf2x y );
int gf2x_gcd( gf2x* g, gf2x x, gf2x y);
int gf2x_lcm( gf2x* l, gf2x x, gf2x y);
int gf2x_mod( gf2x* res, gf2x poly, gf2x modulus );
int gf2x_modinv( gf2x * inv, gf2x base, gf2x modulus );
int gf2x_modexp( gf2x * res, gf2x base, long int exp, gf2x modulus );

int gf2x_minpoly( gf2x * res, gf2x elm, gf2x modulus );


int gf2x_print( gf2x p );

#endif

