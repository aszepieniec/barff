#include "gfp.h"

#include <stdio.h>

gfp_element gfp( int castee )
{
    gfp_element e;
    e = bi_cast(castee);
    return e;
}

gfp_element gfp_init( unsigned int size )
{
    return bi_init(size);
}

gfp_element gfp_clone( gfp_element elm )
{
    return bi_clone(elm);
}

int gfp_destroy( gfp_element elm )
{
    return bi_destroy(elm);
}

int gfp_copy( gfp_element * dest, gfp_element source )
{
    return bi_copy(dest, source);
}

int gfp_zero( gfp_element* elm )
{
    return bi_zero(elm);
}

int gfp_one( gfp_element* elm )
{
    return bi_one(elm);
}

int gfp_random( gfp_element* elm, unsigned char * randomness )
{
    bi_random(elm, GFP_NUMBITS, randomness);
    bi_modulo(elm, *elm, prime_modulus);
    return 1;
}

int gfp_random_invertible( gfp_element* elm, unsigned char * randomness )
{
    bi phi, one;
    phi = bi_init(0);
    one = bi_cast(1);
    bi_subtract(&phi, prime_modulus, one);
    bi_random(elm, GFP_NUMBITS, randomness);
    bi_modulo(elm, *elm, phi);
    bi_add(elm, *elm, one);
    bi_destroy(one);
    bi_destroy(phi);
}

int gfp_compare( gfp_element lhs, gfp_element rhs )
{
    return !bi_compare(lhs, rhs);
}

int gfp_is_zero( gfp_element elm )
{
    return bi_is_zero(elm);
}

int gfp_is_one( gfp_element elm )
{
    return bi_is_one(elm);
}

int gfp_add( gfp_element * res, gfp_element lhs, gfp_element rhs )
{
    bi_add(res, lhs, rhs);
    bi_modulo(res, *res, prime_modulus);
    return 1;
}

int gfp_subtract( gfp_element * res, gfp_element lhs, gfp_element rhs )
{
    bi_subtract(res, lhs, rhs);
    bi_modulo(res, *res, prime_modulus);
    return 1;
}

int gfp_negate( gfp_element * res, gfp_element elm )
{
    return bi_subtract(res, prime_modulus, elm);
}

int gfp_multiply( gfp_element * res, gfp_element lhs, gfp_element rhs )
{
    bi_multiply(res, lhs, rhs);
    bi_modulo(res, *res, prime_modulus);
    return 1;
}

int gfp_divide( gfp_element * quo, gfp_element numerator, gfp_element divisor )
{
    gfp_element e;
    e = gfp_clone(divisor);
    gfp_inverse(&e, divisor);
    gfp_multiply(quo, numerator, e);
    return 1;
}

int gfp_inverse( gfp_element * res, gfp_element elm )
{
    bi x, y, g;
    x = bi_init(0);
    y = bi_init(0);
    g = bi_init(0);
    bi_xgcd(&x, &y, &g, elm, prime_modulus);
    bi_modulo(res, x, prime_modulus);
    bi_destroy(x);
    bi_destroy(y);
    bi_destroy(g);
    return 1;
}

int gfp_print( gfp_element elm )
{
    bi_print(elm);
    return 1;
}


