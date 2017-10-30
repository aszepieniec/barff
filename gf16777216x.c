#include "gf16777216x.h"
#include "gf256x.h"
#include <stdlib.h>
#include <stdio.h>

/**
 * gf16777216_multiply
 * Multiply two GF(16777216) elements as GF(256)[z] polynomials modulo
 * m(z) = z^2 + 0x01*z + 0x2d.
 * @params:
 *  * lhs, rhs : GF(16777216) elements to multiply
 * @return:
 *  * product of lhs and rhs
 */
unsigned int gf16777216_multiply( unsigned int lhs, unsigned int rhs )
{
    gf256x modulus, quotient, remainder, left, right, product;
    unsigned int result;

    if( (lhs & 0xffffff) == 0 || (rhs & 0xffffff) == 0 )
    {
        return 0;
    }

    modulus = gf256x_init(3);
    modulus.data[0] = 0x01;
    modulus.data[1] = 0xe0;
    modulus.data[2] = 0x00;
    modulus.data[3] = 0x01;

    quotient = gf256x_init(0);
    remainder = gf256x_init(0);

    left = gf256x_init(2);
    left.data[0] = lhs & 0xff;
    left.data[1] = (lhs >> 8) & 0xff;
    left.data[2] = (lhs >> 16) & 0xff;

    right = gf256x_init(2);
    right.data[0] = rhs & 0xff;
    right.data[1] = (rhs >> 8) & 0xff;
    right.data[2] = (rhs >> 16) & 0xff;

    product = gf256x_init(0);
    gf256x_multiply(&product, left, right);

    if( product.degree >= modulus.degree )
    {
        gf256x_divide(&quotient, &remainder, product, modulus);
    }
    else
    {
        gf256x_copy(&remainder, product);
    }

    if( gf256x_is_zero(remainder) == 1 )
    {
        result = 0;
    }
    else if( remainder.degree == 0 )
    {
        result = remainder.data[0];
    }
    else if( remainder.degree == 1 )
    {
        result = remainder.data[0] | (remainder.data[1] << 8);
    }
    else
    {
        result = remainder.data[0] | (remainder.data[1] << 8) | (remainder.data[2] << 16);
    }

    gf256x_destroy(modulus);
    gf256x_destroy(quotient);
    gf256x_destroy(remainder);
    gf256x_destroy(product);
    gf256x_destroy(left);
    gf256x_destroy(right);

    return result;
}

/**
 * gf16777216_inverse
 * Find the inverse of a GF(16777216) element as a GF(256)[z] polynomial
 * modulo m(z) = z^3 + 0xe0*z + 0x01.
 * @params:
 *  * elm : GF(16777216) element whose inverse is to be found
 * @return:
 *  * inverse of elm
 */
unsigned int gf16777216_inverse( unsigned int elm )
{
    unsigned int inv;
    gf256x a, b, g, x, y;

    if( elm & 0xffffff == 0 )
    {
        return 0;
    }

    a = gf256x_init(0);
    b = gf256x_init(0);
    g = gf256x_init(0);
    x = gf256x_init(2);
    y = gf256x_init(3);

    x.data[0] = elm & 0xff;
    x.data[1] = (elm >> 8) & 0xff;
    x.data[2] = (elm >> 16) & 0xff;

    y.data[0] = 0x01;
    y.data[1] = 0xe0;
    y.data[2] = 0x00;
    y.data[3] = 0x01;

    gf256x_xgcd(&a, &b, &g, x, y);

    if( a.degree == 0 )
    {
        inv = a.data[0];
    }
    else if( a.degree == 1 )
    {
        inv = a.data[0] | (a.data[1] << 8);
    }
    else
    {
        inv = a.data[0] | (a.data[1] << 8) | (a.data[2] << 16);
    }

    gf256x_destroy(a);
    gf256x_destroy(b);
    gf256x_destroy(g);
    gf256x_destroy(x);
    gf256x_destroy(y);

    return inv;
}

/**
 * gf16777216_exp
 * Raise a given GF(16777216) element to the given power.
 */
unsigned int gf16777216_exp( unsigned int element, int exponent )
{
    unsigned int acc;
    int i;

    if( exponent == 0 )
    {
        return 1;
    }

    if( exponent < 0 )
    {
        return gf16777216_exp(gf16777216_inverse(element), -exponent);
    }

    acc = 0;
    for( i = sizeof(int) * 8 - 1 ; i >= 0 ; ++i )
    {
        acc = gf16777216_multiply(acc, acc);
        if( exponent & (1 << i) != 0 )
        {
            acc = gf16777216_multiply(acc, element);
        }
    }

    return acc;
}

/**
 * gf16777216x_init
 * Initialize a GF(16777216)[x] object of given degree. Allocate memory
 * and set to zero.
 */
gf16777216x gf16777216x_init( int deg )
{
    gf16777216x elm;
    elm.data = malloc(3*deg+3);
    elm.degree = deg;
    return elm;
}

/**
 * gf16777216x_zero
 * Set the given polynomial to zero.
 */
int gf16777216x_zero( gf16777216x* p )
{
    free(p->data);
    p->degree = 0;
    p->data = malloc(3);
    p->data[0] = 0;
    p->data[1] = 0;
    p->data[2] = 0;
    return 1;
}

/**
 * gf16777216x_one
 * Set the given polynomial to one.
 */
int gf16777216x_one( gf16777216x* p )
{
    free(p->data);
    p->degree = 0;
    p->data = malloc(3);
    p->data[0] = 1;
    p->data[1] = 0;
    p->data[2] = 0;
    return 1;
}

/**
 * gf16777216x_copy
 * Copy a GF(16777216)[x] element from one container to another, and
 * reinitialize as necessary.
 */
int gf16777216x_copy( gf16777216x* dest, gf16777216x source )
{
    int i;
    if( dest->degree != source.degree )
    {
        free(dest->data);
        dest->data = malloc(3*source.degree+3);
        dest->degree = source.degree;
    }
    for( i = 0 ; i < 3 + 3*source.degree ; ++i )
    {
        dest->data[i] = source.data[i];
    }
    return 1;
}

/**
 * gf16777216x_destroy
 * Destroy a GF(16777216)[x] object. Free memory.
 */
int gf16777216x_destroy( gf16777216x p )
{
    free(p.data);
}

/**
 * gf16777216x_add
 * Add two GF(16777216)[x] elements together.
 */
int gf16777216x_add( gf16777216x* dest, gf16777216x lhs, gf16777216x rhs )
{
    int i;
    unsigned char * data;
    if( rhs.degree > lhs.degree )
    {
        return gf16777216x_add(dest, rhs, lhs);
    }

    data = malloc(3*lhs.degree+3);

    for( i = 0 ; i < 3 + 3*rhs.degree ; ++i )
    {
        data[i] = lhs.data[i] ^ rhs.data[i];
    }
    for( ; i < 3 + 3*lhs.degree ; ++i )
    {
        data[i] = lhs.data[i];
    }

    free(dest->data);
    dest->degree = lhs.degree;
    dest->data = data;

    while( data[3*dest->degree] == 0 && data[3*dest->degree+1] == 0 && data[3*dest->degree+2] == 0 && dest->degree > 0 )
    {
        dest->degree -= 1;
    }

    return 1;
}

/**
 * gf16777216x_multiply
 * Multiply two GF(16777216)[x] elements together.
 */
int gf16777216x_multiply( gf16777216x* dest, gf16777216x lhs, gf16777216x rhs )
{
    int i, j;
    int degree;
    unsigned int product, left, right;
    unsigned char * data;

    degree = lhs.degree + rhs.degree;
    data = malloc(3*degree + 3);

    for( i = 0 ; i < 3 + 3*degree ; ++i )
    {
        data[i] = 0;
    }

    product = 0;
    for( i = 0 ; i < 1 + lhs.degree ; ++i )
    {
        for( j = 0 ; j < 1 + rhs.degree ; ++j )
        {
            left = lhs.data[3*i] | (lhs.data[3*i + 1] << 8) | (lhs.data[3*i + 2] << 16);
            right = rhs.data[3*j] | (rhs.data[3*j + 1] << 8) | (rhs.data[3*j + 2] << 16);
            product = gf16777216_multiply(left, right);
            data[3*(i+j)] ^= product & 0xff;
            data[3*(i+j) + 1] ^= (product >> 8) & 0xff;
            data[3*(i+j) + 2] ^= (product >> 16) & 0xff;
        }
    }

    free(dest->data);
    dest->data = data;
    dest->degree = degree;

    return 1;
}

/**
 * gf16777216x_equals
 * Decide if two elements of GF(16777216)[x] are equal, and return 1 if so.
 * (Return 0 otherwise.)
 */
int gf16777216x_equals( gf16777216x lhs, gf16777216x rhs )
{
    int i;
    int equal;
    if( lhs.degree != rhs.degree )
    {
        return 0;
    }
    equal = 1;
    for( i = 0 ; i < 3 + 3*lhs.degree ; ++i )
    {
        equal = equal & (lhs.data[i] == rhs.data[i]);
    }
    return equal;
}

/**
 * gf16777216x_is_zero
 * Determine if the given polynomial is equal to zero. Return one if
 * so, zero otherwise.
 */
int gf16777216x_is_zero( gf16777216x p )
{
    int zero;
    int i;
    zero = 1;
    for( i = 0 ; i < 3 + 3*p.degree ; ++i )
    {
        zero &= (p.data[i] == 0);
    }
    return zero;
}

/**
 * gf16777216x_multiply_constant_shift
 * Multiply the polynomial with a constant and shift it (to the left,
 * i.e., towards higher degree). Satisfies:
 *  dest == constant * x^shift * poly
 */
int gf16777216x_multiply_constant_shift( gf16777216x* dest, gf16777216x poly, unsigned int constant, int shift )
{
    unsigned char * data;
    int i;
    int degree;
    unsigned int lhs, product;

    degree = shift + poly.degree;

    data = malloc(3*degree+3);
    for( i = 0 ; i < 3*shift ; ++i )
    {
        data[i] = 0;
    }

    for( i = shift ; i < 1 + degree ; ++i )
    {
        lhs = poly.data[3*(i-shift)] | (poly.data[3*(i-shift) + 1] << 8) | (poly.data[3*(i-shift) + 2] << 16);
        product = gf16777216_multiply(lhs, constant);
        data[3*i] = product & 0xff;
        data[3*i + 1] = (product >> 8) & 0xff;
        data[3*i + 2] = (product >> 16) & 0xff;
    }

    free(dest->data);
    dest->data = data;
    dest->degree = degree;

    return 1;
}

/**
 * gf16777216x_divide
 * Divide one GF(16777216)[x] element by another and record the quotient
 * and remainder.
 */
int gf16777216x_divide( gf16777216x* quo, gf16777216x* rem, gf16777216x num, gf16777216x divisor )
{
    int i, j;
    unsigned int lc, inv, compl, r;
    gf16777216x remainder, poly, quotient;
    unsigned char * quotient_data;

    /* make sure divisor leading coefficient is not zero */
    if( divisor.data[3*divisor.degree] == 0 && divisor.data[3*divisor.degree+1] && divisor.data[3*divisor.degree+2] == 0 )
    {
        poly.data = malloc(3*divisor.degree + 3);
        for( i = 0 ; i < 1 + divisor.degree ; ++i )
        {
            poly.data[3*i] = divisor.data[3*i];
            poly.data[3*i+1] = divisor.data[3*i+1];
            poly.data[3*i+2] = divisor.data[3*i+2];
        }
        for( poly.degree = divisor.degree-1 ; poly.degree > 0 ; --poly.degree )
        {
            if( poly.data[3*poly.degree] != 0 || poly.data[3*poly.degree + 1] != 0 || poly.data[3*poly.degree + 2] != 0 )
            {
                break;
            }
        }
        gf16777216x_divide(quo, rem, num, poly);
        free(poly.data);
        return 1;
    }

    /* make sure numerator leading coefficient is not zero */
    if( num.data[3*num.degree] == 0 && num.data[3*num.degree+1] == 0 && num.data[3*num.degree+1] )
    {
        poly.data = malloc(3*num.degree + 3);
        for( i = 0 ; i < 3 + 3*num.degree ; ++i )
        {
            poly.data[i] = num.data[i];
        }
        for( poly.degree = num.degree-1 ; poly.degree > 0 ; --poly.degree )
        {
            if( poly.data[3*poly.degree] != 0 || poly.data[3*poly.degree+1] != 0 || poly.data[3*poly.degree+2] != 0 )
            {
                break;
            }
        }
        gf16777216x_divide(quo, rem, poly, divisor);
        free(poly.data);
        return 1;
    }

    /* make sure deg(divisor) > deg(numerator) */
    if( divisor.degree > num.degree )
    {
        gf16777216x_zero(quo);
        gf16777216x_copy(rem, num);
        return 1;
    }

    /* filtered out edge cases, proceed with division already */
    remainder = gf16777216x_init(0);
    poly = gf16777216x_init(0);
    gf16777216x_copy(&remainder, num);
    quotient = gf16777216x_init(num.degree - divisor.degree);
    quotient_data = malloc(3*(num.degree - divisor.degree + 1));
    for( i = 0 ; i < 3 + 3*(num.degree - divisor.degree) ; ++i )
    {
        quotient_data[i] = 0;
    }

    lc = divisor.data[3*divisor.degree] | (divisor.data[3*divisor.degree + 1] << 8) | (divisor.data[3*divisor.degree + 2] << 16);
    inv = gf16777216_inverse(lc);
    
    for( i = remainder.degree - divisor.degree ; i >= 0 ; --i )
    {
        if( remainder.degree < divisor.degree + i )
        {
            continue;
        }

        r = remainder.data[3*remainder.degree] | (remainder.data[3*remainder.degree + 1] << 8) | (remainder.data[3*remainder.degree + 2] << 16);
        compl = gf16777216_multiply(r, inv);

        gf16777216x_multiply_constant_shift(&poly, divisor, compl, i);

        quotient_data[3*i] = compl & 0xff;
        quotient_data[3*i+1] = (compl >> 8) & 0xff;
        quotient_data[3*i+2] = (compl >> 16) & 0xff;

        gf16777216x_add(&remainder, remainder, poly);

    }

    free(quo->data);
    quo->data = quotient_data;
    quo->degree = num.degree - divisor.degree;

    gf16777216x_copy(rem, remainder);
    gf16777216x_destroy(remainder);
    gf16777216x_destroy(poly);
    gf16777216x_destroy(quotient);

    return 1;
}

/**
 * gf16777216x_xgcd
 * Compute the greatest common divisor g and Bezout coefficients a
 * and b for x and y using the extended Euclidean algorithm.
 */
int gf16777216x_xgcd( gf16777216x* a, gf16777216x* b, gf16777216x* g, gf16777216x x, gf16777216x y )
{
    gf16777216x s, old_s;
    gf16777216x t, old_t;
    gf16777216x r, old_r;
    gf16777216x quotient, remainder;
    gf16777216x temp;
    gf16777216x temp2;
    unsigned int lc;

    s = gf16777216x_init(0);
    old_s = gf16777216x_init(0);
    t = gf16777216x_init(0);
    old_t = gf16777216x_init(0);
    r = gf16777216x_init(0);
    old_r = gf16777216x_init(0);
    quotient = gf16777216x_init(0);
    remainder = gf16777216x_init(0);
    temp = gf16777216x_init(0);
    temp2 = gf16777216x_init(0);

    gf16777216x_zero(&s);
    gf16777216x_one(&old_s);
    gf16777216x_one(&t);
    gf16777216x_zero(&old_t);
    gf16777216x_copy(&r, y);
    gf16777216x_copy(&old_r, x);

    while( gf16777216x_is_zero(r) == 0 ) /* while r =/= 0 */
    {
        gf16777216x_divide(&quotient, &remainder, old_r, r);

        gf16777216x_copy(&old_r, r);
        gf16777216x_copy(&r, remainder);

        gf16777216x_multiply(&temp, quotient, s);
        gf16777216x_add(&temp, temp, old_s);
        gf16777216x_copy(&old_s, s);
        gf16777216x_copy(&s, temp);

        gf16777216x_multiply(&temp, quotient, t);
        gf16777216x_add(&temp, temp, old_t);
        gf16777216x_copy(&old_t, t);
        gf16777216x_copy(&t, temp);


        gf16777216x_multiply(&temp, old_s, x);
        gf16777216x_multiply(&temp2, old_t, y);
        gf16777216x_add(&temp, temp, temp2);
        if( gf16777216x_equals(temp, old_r) != 1 )
        {
            printf("ab + xy != g \n");
            getchar();
        }
    }

    gf16777216x_copy(a, old_s);
    gf16777216x_copy(b, old_t);
    gf16777216x_copy(g, old_r);

    lc = g->data[2*g->degree] | (g->data[2*g->degree+1] << 8);
    lc = gf16777216_inverse(lc);
    gf16777216x_multiply_constant_shift(g, *g, lc, 0);
    gf16777216x_multiply_constant_shift(a, *a, lc, 0);
    gf16777216x_multiply_constant_shift(b, *b, lc, 0);

    gf16777216x_destroy(s);
    gf16777216x_destroy(old_s);
    gf16777216x_destroy(t);
    gf16777216x_destroy(old_t);
    gf16777216x_destroy(r);
    gf16777216x_destroy(old_r);
    gf16777216x_destroy(quotient);
    gf16777216x_destroy(remainder);
    gf16777216x_destroy(temp);
    gf16777216x_destroy(temp2);

    return 1;
}

/**
 * gf16777216x_modexp
 * Raise a given polynomial to a given power modulo another
 * polynomial.
 */
int gf16777216x_modexp( gf16777216x* dest, gf16777216x base, unsigned long int exponent, gf16777216x modulus )
{
    gf16777216x raised, quo;
    int i;

    raised = gf16777216x_init(0);
    quo = gf16777216x_init(0);
    raised.data[0] = 1;
    raised.data[1] = 0;
    raised.data[2] = 0;

    for( i = sizeof(unsigned long int) * 8 - 1 ; i >= 0 ; --i )
    {
        gf16777216x_multiply(&raised, raised, raised);
        if( (exponent >> i) & 1 != 0 )
        {
            gf16777216x_multiply(&raised, raised, base);
        }
        gf16777216x_divide(&quo, &raised, raised, modulus);
    }

    gf16777216x_destroy(quo);
    gf16777216x_destroy(*dest);
    dest->degree = raised.degree;
    dest->data = raised.data;

    return 1;
}

/**
 * gf16777216x_eval
 * Evaluate the given polynomial in a given point.
 */
unsigned int gf16777216x_eval( gf16777216x polynomial, unsigned int point )
{
    int i;
    unsigned int acc;
    unsigned int xi;
    unsigned int coeff;
    acc = 0;
    xi = 1;
    for( i = 0 ; i < 1 + polynomial.degree ; ++i )
    {
        coeff = polynomial.data[3*i] | (polynomial.data[3*i+1] << 8) | (polynomial.data[3*i+2] << 16);
        acc = acc ^ gf16777216_multiply(coeff, xi);
        xi = gf16777216_multiply(xi, point);
    }

//    printf("evaluating polynomial "); gf16777216x_print(polynomial); printf(" int point %02x; result: %02x\n", point, acc);

    return acc;
}

/**
 * gf16777216x_print
 * Cast the polynomial's coefficients to hex number and throw them to
 * stdout.
 */
int gf16777216x_print( gf16777216x p )
{
    int i;

    for( i = 0 ; i < 1 + p.degree ; ++i )
    {
        printf("%02x%02x%02x", p.data[2*i], p.data[2*i+1], p.data[2*i+2]);
    }

    return 1;
}

