#include "gf2x.h"
#include <stdlib.h>
#include <stdio.h>

/**
 * gf2x_init
 * Initialize a GF(256)[x] object of given degree. Allocate memory
 * and set to zero.
 */
gf2x gf2x_init( int deg )
{
    gf2x elm;
    int num_bytes;

    num_bytes = 1;
    if( (deg+1+7)/8 > num_bytes )
    {
        num_bytes = (deg+1+7)/8;
    }
    elm.data = malloc(num_bytes);
    elm.degree = deg;
    return elm;
}

/**
 * gf2x_zero
 * Set the given polynomial to zero.
 */
int gf2x_zero( gf2x* p )
{
    free(p->data);
    p->degree = -1;
    p->data = malloc(1);
    p->data[0] = 0;
    return 1;
}

/**
 * gf2x_one
 * Set the given polynomial to one.
 */
int gf2x_one( gf2x* p )
{
    free(p->data);
    p->degree = 0;
    p->data = malloc(1);
    p->data[0] = 1;
    return 1;
}

/**
 * gf2x_copy
 * Copy a GF(256)[x] element from one container to another, and
 * reinitialize as necessary.
 */
int gf2x_copy( gf2x* dest, gf2x source )
{
    int i;
    gf2x temp;

    temp = gf2x_init(source.degree);

    for( i = 0 ; i < (source.degree+1+7)/8 ; ++i )
    {
        temp.data[i] = source.data[i];
    }

    gf2x_destroy(*dest);
    *dest = temp;
    return 1;
}

/**
 * gf2x_destroy
 * Destroy a GF(256)[x] object. Free memory.
 */
int gf2x_destroy( gf2x p )
{
    free(p.data);
}

/**
 * gf2x_add
 * Add two GF(2)[x] elements together.
 */
int gf2x_add( gf2x* dest, gf2x lhs, gf2x rhs )
{
    int i;
    unsigned char * data;

    if( rhs.degree < 0 )
    {
        gf2x_copy(dest, lhs);
        return 1;
    }

    if( lhs.degree < 0 )
    {
        gf2x_copy(dest, rhs);
        return 1;
    }

    if( rhs.degree > lhs.degree )
    {
        return gf2x_add(dest, rhs, lhs);
    }

    data = malloc((lhs.degree+1+7)/8);

    for( i = 0 ; i <= lhs.degree/8 ; ++ i )
    {
        data[i] = 0;
    }
    for( i = 0 ; i < (rhs.degree)/8 ; ++i )
    {
        data[i] ^= lhs.data[i] ^ rhs.data[i];
    }
    data[i] ^= lhs.data[i] ^ (rhs.data[i] & (0xff >> (7 - (rhs.degree % 8))));
    i += 1;
    for( ; i < (lhs.degree+1+7)/8 ; ++i )
    {
        data[i] ^= lhs.data[i];
    }

    free(dest->data);
    dest->degree = lhs.degree;
    dest->data = data;

    //printf("lhs: "); gf2x_print(lhs); printf("\n");
    //printf("rhs: "); gf2x_print(rhs); printf("\n");
    //printf("res: "); gf2x_print(*dest); printf("\n");

    gf2x_trim(dest);

    return 1;
}

/**
 * gf2x_multiply
 * Multiple two GF(256)[x] elements together.
 */
int gf2x_multiply( gf2x* dest, gf2x lhs, gf2x rhs )
{
    int i, j;
    int degree;
    unsigned char * data;
    gf2x temp, shifted;

    if( gf2x_is_zero(rhs) == 1 || gf2x_is_zero(lhs) == 1 )
    {
        gf2x_zero(dest);
        return 1;
    }

    degree = lhs.degree + rhs.degree;
    temp = gf2x_init(0);
    shifted = gf2x_init(0);
    gf2x_zero(&temp);

    for( i = 0 ; i <= lhs.degree ; ++i )
    {
        if( (lhs.data[i/8] & (1 << (i%8))) != 0 )
        {
            gf2x_shift_left(&shifted, rhs, i);
            gf2x_add(&temp, temp, shifted);
        }
    }

    free(dest->data);
    dest->data = temp.data;
    dest->degree = temp.degree;
    gf2x_destroy(shifted);

    return 1;
}

/**
 * gf2x_equals
 * Decide if two elements of GF(256)[x] are equal, and return 1 if so.
 * (Return 0 otherwise.)
 */
int gf2x_equals( gf2x lhs, gf2x rhs )
{
    int i;
    int equal;
    if( lhs.degree != rhs.degree )
    {
        //printf("degrees are different.\n");
        return 0;
    }
    equal = 1;
    for( i = 0 ; i < lhs.degree/8 ; ++i )
    {
        equal = equal & ((lhs.data[i] ^ rhs.data[i]) == 0);
        //if( (lhs.data[i] ^ rhs.data[i]) != 0 )
        //    printf("found difference.\n");
    }
    equal = equal & ((unsigned char)(lhs.data[lhs.degree/8] << (8 - (lhs.degree % 8))) ^ (unsigned char)(rhs.data[rhs.degree/8] << (8 - (lhs.degree % 8)))) == 0;
    //if( (unsigned char)(((lhs.data[lhs.degree/8] ^ rhs.data[rhs.degree/8]) << (8 - (lhs.degree % 8)))) != 0 )
    //{
    //    printf("%02x vs %02x for degree mod 8 %i\n", (unsigned char)(lhs.data[lhs.degree/8] << (8 - (lhs.degree % 8))), (unsigned char)(rhs.data[rhs.degree/8] << (8 - (rhs.degree % 8))), lhs.degree);
    //    printf("found difference at end.\n");
    //}
    return equal;
}

/**
 * gf2x_is_zero
 * Determine if the given polynomial is equal to zero. Return one if
 * so, zero otherwise.
 */
int gf2x_is_zero( gf2x p )
{
    int zero;
    int i;
    zero = 1;
    for( i = 0 ; i < (p.degree+1+7)/8 ; ++i )
    {
        zero &= (p.data[i] == 0);
    }
    return zero;
}

/**
 * gf2x_is_one
 * Determine if the given polynomial is equal to one. Return one if
 * so, zero otherwise.
 */
int gf2x_is_one( gf2x p )
{
    int one;
    int i;
    one = 1;
    for( i = 1 ; i < (p.degree+1+7)/8 ; ++i )
    {
        one &= (p.data[i] == 0);
    }
    one &= (p.data[0] == 1);
    return one;
}

/**
 * gf2x_shift_left
 * Apply a shift it to the polynomial (to the left,
 * i.e., towards higher degree). Satisfies:
 *  dest == x^shift * poly
 */
int gf2x_shift_left( gf2x* dest, gf2x poly, unsigned int shift )
{
    unsigned char * data;
    int i;
    int degree;

    degree = shift + poly.degree;

    data = malloc((degree+1+7)/8+1);
    if( shift % 8 == 0 )
    {
        for( i = 0 ; i < shift/8 ; ++i )
        {
            data[i] = 0;
        }
        for( i = 0 ; i <= poly.degree/8 ; ++ i)
        {
            data[i + (shift/8)] = poly.data[i];
        }
    }
    else
    {
        for( i = 0 ; i <= degree/8 ; ++i )
        {
            data[i] = 0;
        }
        for( i = 0 ; i < (poly.degree+1+7)/8 ; ++i )
        {
            data[i+(shift/8)] ^= poly.data[i] << ((shift) % 8);
            data[i+(shift/8)+1] ^= poly.data[i] >> ((8-shift) % 8);
        }
    }

    free(dest->data);
    dest->data = data;
    dest->degree = degree;

    //printf("pre-shift: ");
    //for( i = 0 ; i < shift ; ++i ) printf("-");
    //gf2x_print(poly); printf("\n");


    //printf("postshift: ");
    //gf2x_print(*dest); printf("\n");

    return 1;
}

/**
 * gf2x_trim
 * Do some housekeeping: reduce the degree so that the coefficient
 * of the resulting degree is not zero. (But reduce no further than
 * necessary, and don't change any coefficients.)
 */
int gf2x_trim( gf2x* poly )
{
    int i;

    /* first, find out if we really have to do something at all */
    if( poly->degree < 0 || (poly->data[poly->degree/8] & (1 << (poly->degree % 8))) != 0 )
    {
        return 1;
    }
    else
    {
        //printf("trimming "); gf2x_print(*poly); printf("\n");
    }

    for( i = poly->degree % 8 ; i >= 0 ; --i )
    {
        if( poly->degree == 0 )
        {
            break;
        }
        if( (poly->data[poly->degree/8] & (1 << (poly->degree % 8))) != 0 )
        {
            return 1;
        }
        poly->degree -= 1;
    }

    for( i = poly->degree / 8 ; i > 0 ; --i )
    {
        if( poly->degree == 0 )
        {
            break;
        }
        if( poly->data[i] != 0 )
        {
            break;
        }
        poly->degree -= 8;
    }

    for( i = 7 ; i >= 0 ; --i )
    {
        if( poly->degree == 0 ||  (poly->data[poly->degree/8] & (1 << i)) != 0 )
        {
            break;
        }
        poly->degree -= 1;
    }

    if( poly->degree == 0 && (poly->data[0] & 1) == 0 )
    {
        poly->degree = -1;
    }

    return 1;
}

/**
 * gf2x_divide
 * Divide one GF(256)[x] element by another and record the quotient
 * and remainder.
 */
int gf2x_divide( gf2x* quo, gf2x* rem, gf2x num, gf2x divisor )
{
    int i, j;
    unsigned char inv, compl;
    gf2x remainder, poly;
    gf2x quotient;
    gf2x temp;

    /* make sure divisor leading coefficient is not zero */
    if( divisor.data[(divisor.degree/8)] & (1 << (divisor.degree%8)) == 0 )
    {
        //printf("lc of divisor is not zero\n");
        poly = gf2x_init(0);
        gf2x_copy(&poly, divisor);
        gf2x_trim(&poly);
        gf2x_divide(quo, rem, num, poly);
        gf2x_destroy(poly);
        return 1;
    }

    /* make sure numerator leading coefficient is not zero */
    if( num.data[num.degree/8] & (1 << (num.degree % 8)) == 0 )
    {
        //printf("lc of num is not zero\n");
        poly = gf2x_init(0);
        gf2x_copy(&poly, num);
        gf2x_trim(&poly);
        gf2x_divide(quo, rem, poly, divisor);
        gf2x_destroy(poly);
        return 1;
    }

    /* make sure deg(divisor) > deg(numerator) */
    if( divisor.degree > num.degree )
    {
        gf2x_zero(quo);
        gf2x_copy(rem, num);
        return 1;
    }

    /* filtered out edge cases, proceed with division already */
    remainder = gf2x_init(0);
    poly = gf2x_init(0);
    gf2x_copy(&remainder, num);
    quotient = gf2x_init(num.degree - divisor.degree);
    for( i = 0 ; i < (quotient.degree + 1 + 7)/8 ; ++i )
    {
        quotient.data[i] = 0;
    }

    temp = gf2x_init(0);

    for( i = remainder.degree - divisor.degree ; i >= 0 ; --i )
    {
        if( remainder.degree < divisor.degree + i )
        {
            continue;
        }

        gf2x_shift_left(&poly, divisor, i);

        quotient.data[i/8] ^= (1 << (i%8));
        gf2x_add(&remainder, remainder, poly);
    }

    free(quo->data);
    quo->data = quotient.data;
    quo->degree = num.degree - divisor.degree;

    gf2x_copy(rem, remainder);
    gf2x_destroy(remainder);
    gf2x_destroy(poly);
    gf2x_destroy(temp);

    return 1;
}

/**
 * gf2x_xgcd
 * Compute the greatest common divisor g and Bezout coefficients a
 * and b for x and y using the extended Euclidean algorithm.
 */
int gf2x_xgcd( gf2x* a, gf2x* b, gf2x* g, gf2x x, gf2x y )
{
    gf2x s, old_s;
    gf2x t, old_t;
    gf2x r, old_r;
    gf2x quotient, remainder;
    gf2x temp;

    s = gf2x_init(0);
    old_s = gf2x_init(0);
    t = gf2x_init(0);
    old_t = gf2x_init(0);
    r = gf2x_init(0);
    old_r = gf2x_init(0);
    quotient = gf2x_init(0);
    remainder = gf2x_init(0);
    temp = gf2x_init(0);

    gf2x_zero(&temp);
    gf2x_zero(&s);
    gf2x_one(&old_s);
    gf2x_one(&t);
    gf2x_zero(&old_t);
    gf2x_copy(&r, y);
    gf2x_copy(&old_r, x);

    while( gf2x_is_zero(r) == 0 )
    {
        gf2x_divide(&quotient, &remainder, old_r, r);

        gf2x_copy(&old_r, r);
        gf2x_copy(&r, remainder);

        gf2x_multiply(&temp, quotient, s);
        gf2x_add(&temp, temp, old_s);
        gf2x_copy(&old_s, s);
        gf2x_copy(&s, temp);

        gf2x_multiply(&temp, quotient, t);
        gf2x_add(&temp, temp, old_t);
        gf2x_copy(&old_t, t);
        gf2x_copy(&t, temp);
    }

    gf2x_copy(a, old_s);
    gf2x_copy(b, old_t);
    gf2x_copy(g, old_r);

    gf2x_destroy(s);
    gf2x_destroy(old_s);
    gf2x_destroy(t);
    gf2x_destroy(old_t);
    gf2x_destroy(r);
    gf2x_destroy(old_r);
    gf2x_destroy(quotient);
    gf2x_destroy(remainder);
    gf2x_destroy(temp);

    return 1;
}

/**
 * gf2x_gcd
 * Get the greatest common divisor but we don't care about the Bezout
 * coefficients.
 */
int gf2x_gcd( gf2x* g, gf2x x, gf2x y)
{
    gf2x a, b;
    a = gf2x_init(0);
    b = gf2x_init(0);
    gf2x_xgcd(&a, &b, g, x, y);
    gf2x_destroy(a);
    gf2x_destroy(b);
    return 1;
}


/**
 * gf2x_lcm
 * Get the least common multiple of the two given polynomials.
 */
int gf2x_lcm( gf2x* l, gf2x x, gf2x y)
{
    gf2x g;
    g = gf2x_init(0);
    gf2x_gcd(&g, x, y);
    gf2x_multiply(l, x, y);
    gf2x_divide(l, &g, *l, g);
    gf2x_destroy(g);
    return 1;
}

/**
 * gf2x_mod
 * Compute the remainder of the given polynomial after division by
 * modulus.
 */
int gf2x_mod( gf2x* res, gf2x poly, gf2x modulus )
{
    gf2x quotient;
    quotient = gf2x_init(0);
    gf2x_divide(&quotient, res, poly, modulus);
    gf2x_destroy(quotient);
    return 1;
}

/**
 * gf2x_modinv
 * Compute the modular inverse of the given polynomial (base) modulo
 * the modulus.
 */
int gf2x_modinv( gf2x * inv, gf2x base, gf2x modulus )
{
    int success;
    gf2x g, a, b;
    g = gf2x_init(0);
    a = gf2x_init(0);
    b = gf2x_init(0);

    gf2x_xgcd(&a, &b, &g, base, modulus);

    success = gf2x_is_one(g);

    gf2x_mod(inv, a, modulus);

    gf2x_destroy(g);
    gf2x_destroy(a);
    gf2x_destroy(b);

    return success;
}

/**
 * gf2x_modexp
 * Compute base to the power exp modulo modulus.
 */
int gf2x_modexp( gf2x * res, gf2x base, long int exp, gf2x modulus )
{
    gf2x temp;
    int i;

    if( exp == 0 )
    {
        gf2x_one(res);
        return 1;
    }

    if( exp < 0 )
    {
        temp = gf2x_init(0);
        gf2x_modinv(&temp, base, modulus);
        gf2x_modexp(res, base, -exp, modulus);
        gf2x_destroy(temp);
        return 1;
    }

    temp = gf2x_init(0);
    gf2x_one(&temp);

    for( i = 8*sizeof(long int)-1 ; i >= 0 ; --i )
    {
        gf2x_multiply(&temp, temp, temp);
        gf2x_mod(&temp, temp, modulus);

        if( (exp & (1 << i)) != 0 )
        {
            //gf2x_multiply(&temp, temp, base);
            gf2x_mod(&temp, temp, modulus);
        }
    }

    gf2x_copy(res, temp);
    gf2x_destroy(temp);
    return 1;
}

/**
 * gf2x_print
 */
int gf2x_print( gf2x p )
{
    int i;

    for( i = 0 ; i <= p.degree ; ++i )
    {
        printf("%i", (p.data[i/8] & (1 << (i%8))) != 0);
    }

    return 1;
}

