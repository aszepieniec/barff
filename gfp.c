#include "gfp.h"

#include <stdio.h>

#if GFP_NUMBYTES <= 4

int xgcd( int a, int b, int * x, int * y, int * gcd )
{
    int q, r;
    int u;
    int v;
    int m, n;

    *x = 0;
    *y = 1;
    u = 1;
    v = 0;
    while( a != 0 )
    {
        q = b / a;
        r = b % a;
        m = *x - u*q;
        n = *y - v*q;
        b = a;
        a = r;
        *x = u;
        *y = v;
        u = m;
        v = n;
    }
    *gcd = b;
    
    return *gcd;
}

gfp_element gfp( int castee )
{
    gfp_element e;
    e = gfp_init(sizeof(castee));
    e = ((castee % GF_PRIME_MODULUS) + GF_PRIME_MODULUS) % GF_PRIME_MODULUS;
    return e;
}

gfp_element gfp_init( unsigned int size )
{
    return 0;
}

gfp_element gfp_clone( gfp_element elm )
{
    return elm;
}

int gfp_destroy( gfp_element elm )
{
    return 1;
}

int gfp_copy( gfp_element * dest, gfp_element source )
{
    *dest = source;
    return 1;
}

int gfp_zero( gfp_element* elm )
{
    *elm = 0;
    return 1;
}

int gfp_one( gfp_element* elm )
{
    *elm = 1;
    return 1;
}

int gfp_random( gfp_element* elm, unsigned char * randomness )
{
    int i;
    unsigned int r = 0;
    for( i = 0 ; i < GFP_NUMBYTES + 1 ; ++i )
    {
        r = r * 256 + randomness[i];
    }
    *elm = r % GF_PRIME_MODULUS;
    return 1;
}

int gfp_random_invertible( gfp_element* elm, unsigned char * randomness )
{
    int i;
    unsigned int r = 0;
    for( i = 0 ; i < GFP_NUMBYTES + 1 ; ++i )
    {
        r = r * 256 + randomness[i];
    }
    *elm = 1 + (r % (GF_PRIME_MODULUS-1));
    return 1;
}

int gfp_compare( gfp_element lhs, gfp_element rhs )
{
    if( (int)(lhs) == (int)(rhs) )
    {
        return 1;
    }
    return 0;
}

int gfp_is_one( gfp_element elm )
{
    if( elm == 1 )
    {
        return 1;
    }
    return 0;
}

int gfp_is_zero( gfp_element elm )
{
    if( elm == 0 )
    {
        return 1;
    }
    return 0;
}

int gfp_add( gfp_element * res, gfp_element lhs, gfp_element rhs )
{
    *res = ((int)(lhs) + (int)(rhs)) % GF_PRIME_MODULUS;
    return 1;
}

int gfp_subtract( gfp_element * res, gfp_element lhs, gfp_element rhs )
{
    *res = (GF_PRIME_MODULUS + (int)(lhs) - (int)(lhs)) % GF_PRIME_MODULUS;
    return 1;
}

int gfp_negate( gfp_element * res, gfp_element elm )
{
    *res = GF_PRIME_MODULUS - (int)(elm);
    return 1;
}

int gfp_multiply( gfp_element * res, gfp_element lhs, gfp_element rhs )
{
    *res = ((int)(lhs) * (int)(rhs)) % GF_PRIME_MODULUS;
    return 1;
}

int gfp_divide( gfp_element * quo, gfp_element numerator, gfp_element divisor )
{
    gfp_element e;
    e = gfp_clone(divisor);
    gfp_inverse(&e, divisor);
    *quo = ((int)(numerator) * e) % GF_PRIME_MODULUS;
    return 1;
}

int gfp_inverse( gfp_element * res, gfp_element elm )
{
    int a, b;
    int x, y, g;
    a = elm;
    b = GF_PRIME_MODULUS;
    xgcd(a, b, &x, &y, &g);
    x = (GF_PRIME_MODULUS + (x % GF_PRIME_MODULUS)) % GF_PRIME_MODULUS;
    *res = x;
    return 1;
}

int gfp_print( gfp_element elm )
{
    if( (int)elm > 10 )
    {
        printf(" ");
    }
    printf("%i", (int)elm);
    return 1;
}

#endif


