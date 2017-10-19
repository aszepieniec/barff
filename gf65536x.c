#include "gf65536x.h"
#include "gf256x.h"
#include <stdlib.h>
#include <stdio.h>

/**
 * gf65536_multiply
 * Multiply two GF(65536) elements as GF(256)[z] polynomials modulo
 * m(z) = z^2 + 0x01*z + 0x2d.
 * @params:
 *  * lhs, rhs : GF(65536) elements to multiply
 * @return:
 *  * product of lhs and rhs
 */
unsigned int gf65536_multiply( unsigned int lhs, unsigned int rhs )
{
    gf256x modulus, quotient, remainder, left, right, product;
    unsigned int result;

    if( (lhs & 0xffff) == 0 || (rhs & 0xffff) == 0 )
    {
        return 0;
    }

    modulus = gf256x_init(2);
    modulus.data[2] = 0x01;
    modulus.data[1] = 0x01;
    modulus.data[0] = 0x2d;

    quotient = gf256x_init(0);
    remainder = gf256x_init(0);

    left = gf256x_init(1);
    left.data[0] = lhs & 0xff;
    left.data[1] = (lhs >> 8) & 0xff;

    right = gf256x_init(1);
    right.data[0] = rhs & 0xff;
    right.data[1] = (rhs >> 8) & 0xff;

    product = gf256x_init(0);
    gf256x_multiply(&product, left, right);

    if( product.degree >= modulus.degree )
    {
        //printf("product: "); gf256x_print(product); printf("\n");
        //printf("modulus: "); gf256x_print(modulus); printf("\n");
        gf256x_divide(&quotient, &remainder, product, modulus);
    }
    else
    {
        gf256x_copy(&remainder, product);
    }
    //printf("remainder: "); gf256x_print(remainder); printf("\n");

    if( gf256x_is_zero(remainder) == 1 )
    {
        result = 0;
    }
    else if( remainder.degree == 0 )
    {
        result = remainder.data[0];
    }
    else
    {
        result = remainder.data[0] | (remainder.data[1] << 8);
    }

    //printf("((hex2Fx('"); gf256x_print(left); printf("') * hex2Fx('"); gf256x_print(right); printf("')) %% hex2Fx('"); gf256x_print(modulus); printf("')) == hex2Fx('"); gf256x_print(remainder); printf("')\n");

    gf256x_destroy(modulus);
    gf256x_destroy(quotient);
    gf256x_destroy(remainder);
    gf256x_destroy(product);
    gf256x_destroy(left);
    gf256x_destroy(right);

    return result;
}

/**
 * gf65536_inverse
 * Find the inverse of a GF(65536) element as a GF(256)[z] polynomial
 * modulo m(z) = z^2 + 0x01*z + 0x2d.
 * @params:
 *  * elm : GF(65536) element whose inverse is to be found
 * @return:
 *  * inverse of elm
 */
unsigned int gf65536_inverse( unsigned int elm )
{
    unsigned int inv;
    gf256x a, b, g, x, y;

    if( elm & 65535 == 0 )
    {
        return 0;
    }

    a = gf256x_init(0);
    b = gf256x_init(0);
    g = gf256x_init(0);
    x = gf256x_init(1);
    y = gf256x_init(2);

    x.data[0] = elm & 0xff;
    x.data[1] = (elm >> 8) & 0xff;

    y.data[0] = 0x2d;
    y.data[1] = 0x01;
    y.data[2] = 0x01;

    gf256x_xgcd(&a, &b, &g, x, y);

    if( a.degree == 0 )
    {
        inv = a.data[0];
    }
    else
    {
        inv = a.data[0] | (a.data[1] << 8);
    }

    gf256x_destroy(a);
    gf256x_destroy(b);
    gf256x_destroy(g);
    gf256x_destroy(x);
    gf256x_destroy(y);

    return inv;
}

/**
 * gf65536_exp
 * Raise a given GF(65536) element to the given power.
 */
unsigned int gf65536_exp( unsigned int element, int exponent )
{
    unsigned int acc;
    int i;

    if( exponent == 0 )
    {
        return 1;
    }

    if( exponent < 0 )
    {
        return gf65536_exp(gf65536_inverse(element), -exponent);
    }

    acc = 0;
    for( i = sizeof(int) * 8 - 1 ; i >= 0 ; ++i )
    {
        acc = gf65536_multiply(acc, acc);
        if( exponent & (1 << i) != 0 )
        {
            acc = gf65536_multiply(acc, element);
        }
    }

    return acc;
}

/**
 * gf65536x_init
 * Initialize a GF(65536)[x] object of given degree. Allocate memory
 * and set to zero.
 */
gf65536x gf65536x_init( int deg )
{
    gf65536x elm;
    elm.data = malloc(2*deg+2);
    elm.degree = deg;
    return elm;
}

/**
 * gf65536x_zero
 * Set the given polynomial to zero.
 */
int gf65536x_zero( gf65536x* p )
{
    free(p->data);
    p->degree = 0;
    p->data = malloc(2);
    p->data[0] = 0;
    p->data[1] = 0;
    return 1;
}

/**
 * gf65536x_one
 * Set the given polynomial to one.
 */
int gf65536x_one( gf65536x* p )
{
    free(p->data);
    p->degree = 0;
    p->data = malloc(2);
    p->data[0] = 1;
    p->data[1] = 0;
    return 1;
}

/**
 * gf65536x_copy
 * Copy a GF(65536)[x] element from one container to another, and
 * reinitialize as necessary.
 */
int gf65536x_copy( gf65536x* dest, gf65536x source )
{
    int i;
    if( dest->degree != source.degree )
    {
        free(dest->data);
        dest->data = malloc(2*source.degree+2);
        dest->degree = source.degree;
    }
    for( i = 0 ; i < 2 + 2*source.degree ; ++i )
    {
        dest->data[i] = source.data[i];
    }
    return 1;
}

/**
 * gf65536x_destroy
 * Destroy a GF(65536)[x] object. Free memory.
 */
int gf65536x_destroy( gf65536x p )
{
    free(p.data);
}

/**
 * gf65536x_add
 * Add two GF(65536)[x] elements together.
 */
int gf65536x_add( gf65536x* dest, gf65536x lhs, gf65536x rhs )
{
    int i;
    unsigned char * data;
    if( rhs.degree > lhs.degree )
    {
        return gf65536x_add(dest, rhs, lhs);
    }

    data = malloc(2*lhs.degree+2);

    for( i = 0 ; i < 2 + 2*rhs.degree ; ++i )
    {
        data[i] = lhs.data[i] ^ rhs.data[i];
    }
    for( ; i < 2 + 2*lhs.degree ; ++i )
    {
        data[i] = lhs.data[i];
    }

    free(dest->data);
    dest->degree = lhs.degree;
    dest->data = data;

    while( data[2*dest->degree] == 0 && data[2*dest->degree+1] == 0 && dest->degree > 0 )
    {
        dest->degree -= 1;
    }

    return 1;
}

/**
 * gf65536x_multiply
 * Multiply two GF(65536)[x] elements together.
 */
int gf65536x_multiply( gf65536x* dest, gf65536x lhs, gf65536x rhs )
{
    int i, j;
    int degree;
    unsigned int product, left, right;
    unsigned char * data;

    degree = lhs.degree + rhs.degree;
    data = malloc(2*degree + 2);

    for( i = 0 ; i < 2 + 2*degree ; ++i )
    {
        data[i] = 0;
    }

    product = 0;
    for( i = 0 ; i < 1 + lhs.degree ; ++i )
    {
        for( j = 0 ; j < 1 + rhs.degree ; ++j )
        {
            left = lhs.data[2*i] | (lhs.data[2*i + 1] << 8);
            right = rhs.data[2*j] | (rhs.data[2*j + 1] << 8);
            product = gf65536_multiply(left, right);
            data[2*(i+j)] ^= product & 0xff;
            data[2*(i+j) + 1] ^= (product >> 8) & 0xff;
        }
    }

    free(dest->data);
    dest->data = data;
    dest->degree = degree;

    return 1;
}

/**
 * gf65536x_equals
 * Decide if two elements of GF(65536)[x] are equal, and return 1 if so.
 * (Return 0 otherwise.)
 */
int gf65536x_equals( gf65536x lhs, gf65536x rhs )
{
    int i;
    int equal;
    if( lhs.degree != rhs.degree )
    {
        return 0;
    }
    equal = 1;
    for( i = 0 ; i < 2 + 2*lhs.degree ; ++i )
    {
        equal = equal & (lhs.data[i] == rhs.data[i]);
    }
    return equal;
}

/**
 * gf65536x_is_zero
 * Determine if the given polynomial is equal to zero. Return one if
 * so, zero otherwise.
 */
int gf65536x_is_zero( gf65536x p )
{
    int zero;
    int i;
    zero = 1;
    for( i = 0 ; i < 2 + 2*p.degree ; ++i )
    {
        zero &= (p.data[i] == 0);
    }
    return zero;
}

/**
 * gf65536x_multiply_constant_shift
 * Multiply the polynomial with a constant and shift it (to the left,
 * i.e., towards higher degree). Satisfies:
 *  dest == constant * x^shift * poly
 */
int gf65536x_multiply_constant_shift( gf65536x* dest, gf65536x poly, unsigned int constant, int shift )
{
    unsigned char * data;
    int i;
    int degree;
    unsigned int lhs, product;

    degree = shift + poly.degree;

    data = malloc(2*degree+2);
    for( i = 0 ; i < 2*shift ; ++i )
    {
        data[i] = 0;
    }

    for( i = shift ; i < 1 + degree ; ++i )
    {
        lhs = poly.data[2*(i-shift)] | (poly.data[2*(i-shift) + 1] << 8);
        product = gf65536_multiply(lhs, constant);
        data[2*i] = product & 0xff;
        data[2*i + 1] = (product >> 8) & 0xff;
        //data[i] = gf65536_multiply(poly.data[i-shift], constant);
    }

    free(dest->data);
    dest->data = data;
    dest->degree = degree;

    return 1;
}

/**
 * gf65536x_divide
 * Divide one GF(65536)[x] element by another and record the quotient
 * and remainder.
 */
int gf65536x_divide( gf65536x* quo, gf65536x* rem, gf65536x num, gf65536x divisor )
{
    int i, j;
    unsigned int lc, inv, compl, r;
    gf65536x remainder, poly, quotient;
    unsigned char * quotient_data;

    /* make sure divisor leading coefficient is not zero */
    if( divisor.data[2*divisor.degree] == 0 && divisor.data[2*divisor.degree+1] == 0 )
    {
        poly.data = malloc(2*divisor.degree + 2);
        for( i = 0 ; i < 1 + divisor.degree ; ++i )
        {
            poly.data[2*i] = divisor.data[2*i];
            poly.data[2*i+1] = divisor.data[2*i+1];
        }
        for( poly.degree = divisor.degree-1 ; poly.degree > 0 ; --poly.degree )
        {
            if( poly.data[2*poly.degree] != 0 || poly.data[2*poly.degree + 1] != 0 )
            {
                break;
            }
        }
        gf65536x_divide(quo, rem, num, poly);
        free(poly.data);
        return 1;
    }

    /* make sure numerator leading coefficient is not zero */
    if( num.data[2*num.degree] == 0 && num.data[2*num.degree+1] == 0 )
    {
        poly.data = malloc(2*num.degree + 2);
        for( i = 0 ; i < 2 + 2*num.degree ; ++i )
        {
            poly.data[i] = num.data[i];
        }
        for( poly.degree = num.degree-1 ; poly.degree > 0 ; --poly.degree )
        {
            if( poly.data[2*poly.degree] != 0 || poly.data[2*poly.degree+1] != 0 )
            {
                break;
            }
        }
        gf65536x_divide(quo, rem, poly, divisor);
        free(poly.data);
        return 1;
    }

    /* make sure deg(divisor) > deg(numerator) */
    if( divisor.degree > num.degree )
    {
        gf65536x_zero(quo);
        gf65536x_copy(rem, num);
        return 1;
    }

    /* filtered out edge cases, proceed with division already */
    remainder = gf65536x_init(0);
    poly = gf65536x_init(0);
    gf65536x_copy(&remainder, num);
    quotient = gf65536x_init(num.degree - divisor.degree);
    quotient_data = malloc(2*(num.degree - divisor.degree + 1));
    for( i = 0 ; i < 2 + 2*(num.degree - divisor.degree) ; ++i )
    {
        quotient_data[i] = 0;
    }

    lc = divisor.data[2*divisor.degree] | (divisor.data[2*divisor.degree + 1] << 8);
    inv = gf65536_inverse(lc);
    
    for( i = remainder.degree - divisor.degree ; i >= 0 ; --i )
    {
        if( remainder.degree < divisor.degree + i )
        {
            continue;
        }

        r = remainder.data[2*remainder.degree] | (remainder.data[2*remainder.degree + 1] << 8);
        compl = gf65536_multiply(r, inv);

        gf65536x_multiply_constant_shift(&poly, divisor, compl, i);

        quotient_data[2*i] = compl & 0xff;
        quotient_data[2*i+1] = (compl >> 8) & 0xff;
        //quotient.data[2*i] = compl & 0xff;
        //quotient.data[2*i+1] = (compl >> 8) & 0xff;
    
        gf65536x_add(&remainder, remainder, poly);

    }

    free(quo->data);
    quo->data = quotient_data;
    quo->degree = num.degree - divisor.degree;

    gf65536x_copy(rem, remainder);
    gf65536x_destroy(remainder);
    gf65536x_destroy(poly);
    gf65536x_destroy(quotient);

    return 1;
}

/**
 * gf65536x_xgcd
 * Compute the greatest common divisor g and Bezout coefficients a
 * and b for x and y using the extended Euclidean algorithm.
 */
int gf65536x_xgcd( gf65536x* a, gf65536x* b, gf65536x* g, gf65536x x, gf65536x y )
{
    gf65536x s, old_s;
    gf65536x t, old_t;
    gf65536x r, old_r;
    gf65536x quotient, remainder;
    gf65536x temp;
    gf65536x temp2;
    unsigned int lc;

    s = gf65536x_init(0);
    old_s = gf65536x_init(0);
    t = gf65536x_init(0);
    old_t = gf65536x_init(0);
    r = gf65536x_init(0);
    old_r = gf65536x_init(0);
    quotient = gf65536x_init(0);
    remainder = gf65536x_init(0);
    temp = gf65536x_init(0);
    temp2 = gf65536x_init(0);

    gf65536x_zero(&s);
    gf65536x_one(&old_s);
    gf65536x_one(&t);
    gf65536x_zero(&old_t);
    gf65536x_copy(&r, y);
    gf65536x_copy(&old_r, x);

    while( gf65536x_is_zero(r) == 0 ) /* while r =/= 0 */
    {
        gf65536x_divide(&quotient, &remainder, old_r, r);

        gf65536x_copy(&old_r, r);
        gf65536x_copy(&r, remainder);

        gf65536x_multiply(&temp, quotient, s);
        gf65536x_add(&temp, temp, old_s);
        gf65536x_copy(&old_s, s);
        gf65536x_copy(&s, temp);

        gf65536x_multiply(&temp, quotient, t);
        gf65536x_add(&temp, temp, old_t);
        gf65536x_copy(&old_t, t);
        gf65536x_copy(&t, temp);


        gf65536x_multiply(&temp, old_s, x);
        gf65536x_multiply(&temp2, old_t, y);
        gf65536x_add(&temp, temp, temp2);
        if( gf65536x_equals(temp, old_r) != 1 )
        {
            printf("ab + xy != g \n");
            getchar();
        }
    }

    gf65536x_copy(a, old_s);
    gf65536x_copy(b, old_t);
    gf65536x_copy(g, old_r);

    lc = g->data[2*g->degree] | (g->data[2*g->degree+1] << 8);
    lc = gf65536_inverse(lc);
    gf65536x_multiply_constant_shift(g, *g, lc, 0);
    gf65536x_multiply_constant_shift(a, *a, lc, 0);
    gf65536x_multiply_constant_shift(b, *b, lc, 0);

    gf65536x_destroy(s);
    gf65536x_destroy(old_s);
    gf65536x_destroy(t);
    gf65536x_destroy(old_t);
    gf65536x_destroy(r);
    gf65536x_destroy(old_r);
    gf65536x_destroy(quotient);
    gf65536x_destroy(remainder);
    gf65536x_destroy(temp);
    gf65536x_destroy(temp2);

    return 1;
}

/**
 * gf65536x_eval
 * Evaluate the given polynomial in a given point.
 */
unsigned int gf65536x_eval( gf65536x polynomial, unsigned int point )
{
    int i;
    unsigned int acc;
    unsigned int xi;
    acc = 0;
    xi = 1;
    for( i = 0 ; i < 1 + polynomial.degree ; ++i )
    {
        acc = acc ^ gf65536_multiply(polynomial.data[i], xi);
        xi = gf65536_multiply(xi, point);
    }

//    printf("evaluating polynomial "); gf65536x_print(polynomial); printf(" int point %02x; result: %02x\n", point, acc);

    return acc;
}

/**
 * gf65536x_print
 * Cast the polynomial's coefficients to hex number and throw them to
 * stdout.
 */
int gf65536x_print( gf65536x p )
{
    int i;

    for( i = 0 ; i < 1 + p.degree ; ++i )
    {
        printf("%02x%02x", p.data[2*i+1], p.data[2*i]);
    }

    return 1;
}

