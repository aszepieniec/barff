#include "gf4096x.h"
#include "gf4096_tables.h"
#include <stdlib.h>
#include <stdio.h>

/**
 * gf4096_multiply
 * Multiply two GF(4096) elements as GF(4096)[z] polynomials modulo
 * m(z) = x^12 + x^7 + x^6 + x^5 + x^3 + x + 1.
 * @params:
 *  * lhs, rhs : GF(4096) elements to multiply
 * @return:
 *  * product of lhs and rhs
 */
unsigned int gf4096_multiply( unsigned int lhs, unsigned int rhs )
{
    int a, b;
    if( lhs == 0 || rhs == 0 )
    {
        return 0;
    }
    a = gf4096_dlogs[(lhs & 0xfff)];
    b = gf4096_dlogs[(rhs & 0xfff)];
    return gf4096_antilogs[(a + b) % 4095];
}

/**
 * gf4096_inverse
 * Find the inverse of a GF(4096) element.
 * @params:
 *  * elm : GF(4096) element whose inverse is to be found
 * @return:
 *  * inverse of elm
 */
unsigned int gf4096_inverse( unsigned int elm )
{
    int a;
    unsigned int inv;
    if( elm == 0 )
    {
        return 0;
    }
    a = gf4096_dlogs[(elm & 0xfff)];
    inv = gf4096_antilogs[4095 - a];
    return inv;
}

/**
 * gf4096_exp
 * Raise a given GF(4096) element to the given power.
 */
unsigned int gf4096_exp( unsigned int element, int exponent )
{
    int index;
    index = (4095 + ((gf4096_dlogs[(element & 0xfff)] * exponent) % 4095)) % 4095;
    return gf4096_antilogs[index];
}

/**
 * gf4096x_init
 * Initialize a GF(4096)[x] object of given degree. Allocate memory
 * and set to zero.
 */
gf4096x gf4096x_init( int deg )
{
    gf4096x elm;
    elm.data = malloc(2*deg+2);
    elm.degree = deg;
    return elm;
}

/**
 * gf4096x_zero
 * Set the given polynomial to zero.
 */
int gf4096x_zero( gf4096x* p )
{
    free(p->data);
    p->degree = 0;
    p->data = malloc(2);
    p->data[0] = 0;
    p->data[1] = 0;
    return 1;
}

/**
 * gf4096x_one
 * Set the given polynomial to one.
 */
int gf4096x_one( gf4096x* p )
{
    free(p->data);
    p->degree = 0;
    p->data = malloc(2);
    p->data[0] = 1;
    p->data[1] = 0;
    return 1;
}

/**
 * gf4096x_copy
 * Copy a GF(4096)[x] element from one container to another, and
 * reinitialize as necessary.
 */
int gf4096x_copy( gf4096x* dest, gf4096x source )
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
 * gf4096x_destroy
 * Destroy a GF(4096)[x] object. Free memory.
 */
int gf4096x_destroy( gf4096x p )
{
    free(p.data);
}

/**
 * gf4096x_add
 * Add two GF(4096)[x] elements together.
 */
int gf4096x_add( gf4096x* dest, gf4096x lhs, gf4096x rhs )
{
    int i;
    gf4096x res;

    if( rhs.degree > lhs.degree )
    {
        return gf4096x_add(dest, rhs, lhs);
    }

    res = gf4096x_init(lhs.degree);

    for( i = 0 ; i < 2 + 2*rhs.degree ; ++i )
    {
        res.data[i] = lhs.data[i] ^ rhs.data[i];
    }
    
    for( ; i < 2 + 2*lhs.degree ; ++i )
    {
        res.data[i] = lhs.data[i];
    }

    free(dest->data);
    dest->degree = lhs.degree;
    dest->data = res.data;

    while( res.data[2*dest->degree] == 0 && (res.data[2*dest->degree+1] & 0xf) == 0 && dest->degree > 0 )
    {
        dest->degree -= 1;
    }

    return 1;
}

/**
 * gf4096x_multiply
 * Multiply two GF(4096)[x] elements together.
 */
int gf4096x_multiply( gf4096x* dest, gf4096x lhs, gf4096x rhs )
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
        left = lhs.data[2*i] | ((lhs.data[2*i + 1] & 0xf) << 8);
        for( j = 0 ; j < 1 + rhs.degree ; ++j )
        {
            right = rhs.data[2*j] | ((rhs.data[2*j + 1] & 0xf) << 8);
            product = gf4096_multiply(left, right);
            data[2*(i+j)] ^= product & 0xff;
            data[2*(i+j) + 1] ^= (product >> 8) & 0xf;
        }
    }

    free(dest->data);
    dest->data = data;
    dest->degree = degree;

    return 1;
}

/**
 * gf4096x_equals
 * Decide if two elements of GF(4096)[x] are equal, and return 1 if so.
 * (Return 0 otherwise.)
 */
int gf4096x_equals( gf4096x lhs, gf4096x rhs )
{
    int i;
    int equal;
    if( lhs.degree != rhs.degree )
    {
        return 0;
    }
    equal = 1;
    for( i = 0 ; i < 1 + lhs.degree ; ++i )
    {
        equal = equal & (lhs.data[2*i] == rhs.data[2*i]);
        equal = equal & (((unsigned char)(lhs.data[2*i+1] ^ rhs.data[2*i+1]) & 0xf) == 0);
    }
    return equal;
}

/**
 * gf4096x_is_zero
 * Determine if the given polynomial is equal to zero. Return one if
 * so, zero otherwise.
 */
int gf4096x_is_zero( gf4096x p )
{
    int zero;
    int i;
    zero = 1;
    for( i = 0 ; i < 1 + p.degree ; ++i )
    {
        zero &= (p.data[2*i] == 0);
        zero &= ((p.data[2*i+1] & 0xf) == 0);
    }
    return zero;
}

/**
 * gf4096x_multiply_constant_shift
 * Multiply the polynomial with a constant and shift it (to the left,
 * i.e., towards higher degree). Satisfies:
 *  dest == constant * x^shift * poly
 */
int gf4096x_multiply_constant_shift( gf4096x* dest, gf4096x poly, unsigned int constant, int shift )
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
        lhs = poly.data[2*(i-shift)] | ((unsigned char)(poly.data[2*(i-shift) + 1] & 0xf) << 8);
        product = gf4096_multiply(lhs, constant);
        data[2*i] = product & 0xff;
        data[2*i + 1] = (product >> 8) & 0xf;
    }

    free(dest->data);
    dest->data = data;
    dest->degree = degree;

    return 1;
}

/**
 * gf4096x_divide
 * Divide one GF(4096)[x] element by another and record the quotient
 * and remainder.
 */
int gf4096x_divide( gf4096x* quo, gf4096x* rem, gf4096x num, gf4096x divisor )
{
    int i, j;
    unsigned int lc, inv, compl, r;
    gf4096x remainder, poly, quotient;
    gf4096x temp;

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
        gf4096x_divide(quo, rem, num, poly);
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
        gf4096x_divide(quo, rem, poly, divisor);
        free(poly.data);
        return 1;
    }

    /* make sure deg(divisor) > deg(numerator) */
    if( divisor.degree > num.degree )
    {
        gf4096x_zero(quo);
        gf4096x_copy(rem, num);
        return 1;
    }

    /* filtered out edge cases, proceed with division already */
    remainder = gf4096x_init(0);
    poly = gf4096x_init(0);
    gf4096x_copy(&remainder, num);
    quotient = gf4096x_init(num.degree - divisor.degree);
    for( i = 0 ; i < 2 + 2*(num.degree - divisor.degree) ; ++i )
    {
        quotient.data[i] = 0;
    }

    lc = divisor.data[2*divisor.degree] | ((divisor.data[2*divisor.degree + 1] & 0xf) << 8);
    inv = gf4096_inverse(lc);
    
    for( i = remainder.degree - divisor.degree ; i >= 0 ; --i )
    {
        if( remainder.degree < divisor.degree + i )
        {
            continue;
        }

        r = remainder.data[2*remainder.degree] | ((remainder.data[2*remainder.degree + 1] & 0xf) << 8);
        compl = gf4096_multiply(r, inv);

        gf4096x_multiply_constant_shift(&poly, divisor, compl, i);

        quotient.data[2*i] = compl & 0xff;
        quotient.data[2*i+1] = (compl >> 8) & 0xf;

        gf4096x_add(&remainder, remainder, poly);

        temp = gf4096x_init(0);
        gf4096x_multiply(&temp, quotient, divisor);
        gf4096x_add(&temp, temp, remainder);
        gf4096x_destroy(temp);

    }

    free(quo->data);
    quo->data = quotient.data;
    quo->degree = num.degree - divisor.degree;

    gf4096x_copy(rem, remainder);
    gf4096x_destroy(remainder);
    gf4096x_destroy(poly);

    return 1;
}

/**
 * gf4096x_xgcd
 * Compute the greatest common divisor g and Bezout coefficients a
 * and b for x and y using the extended Euclidean algorithm.
 */
int gf4096x_xgcd( gf4096x* a, gf4096x* b, gf4096x* g, gf4096x x, gf4096x y )
{
    gf4096x s, old_s;
    gf4096x t, old_t;
    gf4096x r, old_r;
    gf4096x quotient, remainder;
    gf4096x temp;
    gf4096x temp2;
    unsigned int lc;

    s = gf4096x_init(0);
    old_s = gf4096x_init(0);
    t = gf4096x_init(0);
    old_t = gf4096x_init(0);
    r = gf4096x_init(0);
    old_r = gf4096x_init(0);
    quotient = gf4096x_init(0);
    remainder = gf4096x_init(0);
    temp = gf4096x_init(0);
    temp2 = gf4096x_init(0);

    gf4096x_zero(&s);
    gf4096x_one(&old_s);
    gf4096x_one(&t);
    gf4096x_zero(&old_t);
    gf4096x_copy(&r, y);
    gf4096x_copy(&old_r, x);

    while( gf4096x_is_zero(r) == 0 ) /* while r =/= 0 */
    {
        gf4096x_divide(&quotient, &remainder, old_r, r);
    
        gf4096x_copy(&old_r, r);
        gf4096x_copy(&r, remainder);

        gf4096x_multiply(&temp, quotient, s);
        gf4096x_add(&temp, temp, old_s);
        gf4096x_copy(&old_s, s);
        gf4096x_copy(&s, temp);

        gf4096x_multiply(&temp, quotient, t);
        gf4096x_add(&temp, temp, old_t);
        gf4096x_copy(&old_t, t);
        gf4096x_copy(&t, temp);


        gf4096x_multiply(&temp, old_s, x);
        gf4096x_multiply(&temp2, old_t, y);
        gf4096x_add(&temp, temp, temp2);
        if( gf4096x_equals(temp, old_r) != 1 )
        {
            printf("ab + xy != g \n");
            getchar();
        }
    }

    gf4096x_copy(a, old_s);
    gf4096x_copy(b, old_t);
    gf4096x_copy(g, old_r);

    lc = g->data[2*g->degree] | (g->data[2*g->degree+1] << 8);
    lc = gf4096_inverse(lc);
    gf4096x_multiply_constant_shift(g, *g, lc, 0);
    gf4096x_multiply_constant_shift(a, *a, lc, 0);
    gf4096x_multiply_constant_shift(b, *b, lc, 0);

    gf4096x_destroy(s);
    gf4096x_destroy(old_s);
    gf4096x_destroy(t);
    gf4096x_destroy(old_t);
    gf4096x_destroy(r);
    gf4096x_destroy(old_r);
    gf4096x_destroy(quotient);
    gf4096x_destroy(remainder);
    gf4096x_destroy(temp);
    gf4096x_destroy(temp2);

    return 1;
}

/**
 * gf4096x_eval
 * Evaluate the given polynomial in a given point.
 */
unsigned int gf4096x_eval( gf4096x polynomial, unsigned int point )
{
    int i;
    unsigned int acc;
    unsigned int xi;
    unsigned int coeff;
    acc = 0;
    xi = 1;
    for( i = 0 ; i < 1 + polynomial.degree ; ++i )
    {
        coeff = polynomial.data[2*i] | (polynomial.data[2*i+1] << 8);
        acc = acc ^ gf4096_multiply(coeff, xi);
        xi = gf4096_multiply(xi, point);
    }

//    printf("evaluating polynomial "); gf4096x_print(polynomial); printf(" int point %02x; result: %02x\n", point, acc);

    return acc;
}

/**
 * gf4096x_print
 * Cast the polynomial's coefficients to hex number and throw them to
 * stdout.
 */
int gf4096x_print( gf4096x p )
{
    int i;

    for( i = 0 ; i < 1 + p.degree ; ++i )
    {
        printf("%x%02x", (p.data[2*i+1] & 0xf), p.data[2*i]);
    }

    return 1;
}

