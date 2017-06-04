#include "bi.h"
#include <stdlib.h>
#include <string.h>
#include <stdio.h>

/**
 * bi_init
 * Initialize (allocate) space for a new integer.
 * @params:
 *  * numberofbits : int -- number of bits to allocate for. One
 *    limb contains sizeof(unsigned long int)*8 bits.
 * @return:
 *  * big integer : binteger -- with freshly mallocated memory all
 *    set to zero.
 */
bi bi_init( int numberofbits )
{
    int i;
    bi binteger;
    unsigned int numberoflimbs;
    numberoflimbs = (numberofbits + sizeof(unsigned long int) * 8 - 1) / (sizeof(unsigned long int)*8);
    if( numberoflimbs == 0 )
    {
        numberoflimbs = numberoflimbs + 1;
    }
    binteger.sign = 1;
    binteger.num_limbs = numberoflimbs;
    binteger.data = malloc(sizeof(unsigned long int) * numberoflimbs);
    for( i = 0 ; i < binteger.num_limbs ; ++i )
    {
        binteger.data[i] = 0;
    }
    return binteger;
}

/**
 * bi_copy
 * Copy a big integer into another big integer object.
 */
int bi_copy( bi * dest, bi source )
{
    int i;
    if( dest->num_limbs < source.num_limbs )
    {
        dest->num_limbs = source.num_limbs;
        free(dest->data);
        dest->data = malloc(sizeof(unsigned long int)*dest->num_limbs);
    }
    dest->sign = source.sign;
    for( i = 0 ; i < source.num_limbs ; ++i )
    {
        dest->data[i] = source.data[i];
    }
    for( i = source.num_limbs ; i < dest->num_limbs ; ++i )
    {
        dest->data[i] = 0;
    }
    dest->num_limbs = source.num_limbs;
    return 1;
}

/**
 * bi_cast
 * Initialize a big integer and set its value to the given (signed)
 * integer.
 * @params
 *  * integer : long int (signed) -- value to set the resulting big
 *    integer to
 * @return
 *  * binteger : bi -- with freshly mallocated memory set to the
 *    value of integer
 */
bi bi_cast( long int integer )
{
    bi binteger;
    
    binteger.num_limbs = 1;
    binteger.data = malloc(sizeof(unsigned long int));
    if( integer >= 0 )
    {
        binteger.sign = 1;
        binteger.data[0] = integer;
    }
    else
    {
        binteger.sign = -1;
        binteger.data[0] = -integer;
    }

    return binteger;
}

/**
 * bi_local
 * Group together all information into a big integer object with,
 * importantly, local (stack-based) memory. No mallocation!
 * @params:
 *  * sign : char -- +1 if positive, -1 if negative
 *  * num_limbs : unsigned int -- number of limbs
 *  * data : (pointer to) array of unsigned long ints -- represents
 *    the mantissa of the integer
 * @return:
 *  * binteger : bi -- big integer object built from the given ar-
 *    guments.
 */
bi bi_local( char sign, int num_limbs, unsigned long int * limbs )
{
    bi binteger;
    binteger.sign = sign;
    binteger.num_limbs = num_limbs;
    binteger.data = limbs;
    return binteger;
}

/**
 * bi_clone
 * Clone a big integer; create a new object with allocated memory.
 * @params:
 *  * source : bit integer -- object to clone
 * @return:
 *  * binteger : big integer -- exact clone of source
 */
bi bi_clone( bi source )
{
    bi binteger;
    binteger.sign = source.sign;
    binteger.num_limbs = source.num_limbs;
    binteger.data = malloc(binteger.num_limbs*sizeof(unsigned long int));
    strncpy((char*)binteger.data, (char*)source.data, binteger.num_limbs*sizeof(unsigned long int));
    return binteger;
}

/**
 * bi_destroy
 * Deallocate memory that was allocated to the given object.
 */
int bi_destroy( bi integer )
{
    free(integer.data);
    return 1;
}

/**
 * bi_one
 * Set the given big integer to hold the value 1.
 */
int bi_one( bi * integer )
{
    int i;
    if( integer->num_limbs < 1 )
    {
        free(integer->data);
        integer->data = malloc(sizeof(unsigned long int));
        integer->num_limbs = 1;
    }
    for( i = 1 ; i < integer->num_limbs ; ++i )
    {
        integer->data[i] = 0;
    }
    integer->data[0] = 1;
    integer->sign = 1;

    return 1;
}

/**
 * bi_zero
 * Set the given big integer to hold the value 1.
 */
int bi_zero( bi * integer )
{
    int i;
    if( integer->num_limbs < 1 )
    {
        free(integer->data);
        integer->data = malloc(sizeof(unsigned long int));
        integer->num_limbs = 1;
    }
    for( i = 0 ; i < integer->num_limbs ; ++i )
    {
        integer->data[i] = 0;
    }
    integer->sign = 1;

    return 1;
}

/**
 * bi_add_limbs
 * Add two limbs together, and record if overflow happens. Speci-
 * fically, the sum of the left and right hand sides is stored in
 * sum; whereas carry will be increased by one if overflow happens.
 * @params:
 *  * sum : pointer to unsigned long int -- designated container for
 *    the sum
 *  * carry : pointer to unsigned long int -- to be incremented by
 *    one if overload happens
 *  * lhs, rhs : unsigned long ints -- operands to sum together
 */
int bi_add_limbs( unsigned long int * sum, unsigned long int * carry, unsigned long int lhs, unsigned long int rhs )
{
    *sum = lhs + rhs;
    if( *sum < lhs || *sum < rhs )
    {
        *carry = *carry + 1;
    }

    return 1;
}

/**
 * bi_subtract_limbs
 * Subtract two limbs, and record if borrowing happens. Specifically,
 * the difference of the left and right hand sides is stored in dif-
 * borrow; whereas borrow will be increased by one if borrow happens.
 * @params:
 *  * difference : pointer to unsigned long int -- designated
 *    container for the difference
 *  * borrow : pointer to unsigned long int -- to be incremented by
 *    one if overload happens
 *  * lhs, rhs : unsigned long ints -- operands to subtract from one
 *    another
 */
int bi_subtract_limbs( unsigned long int * difference, unsigned long int * borrow, unsigned long int lhs, unsigned long int rhs )
{
    *difference = lhs - rhs;
    if( *difference > lhs )
    {
        *borrow = *borrow + 1;
    }

    return 1;
}

/**
 * bi_multiply_limbs
 * Multiply two limbs; record the overflow by increasing the carry
 * container by that amount.
 * @params:
 *  * produce : pointer to unsigned long int -- designated container
 *    for the product
 *  * carry : pointer to unsigend long int -- to be set to
 *    whatever overload occurs
 *  * lhs, rhs : unsigned long ints -- operands to multiply together
 */
int bi_multiply_limbs( unsigned long int * product, unsigned long int * carry, unsigned long int lhs, unsigned long int rhs )
{
    unsigned long int lchar, rchar;
    int i, j, k;
    unsigned long int lhs_uchars[sizeof(unsigned long int)], rhs_uchars[sizeof(unsigned long int)];
    unsigned long int prod_ulongs[2*sizeof(unsigned long int)+1];

    for( i = 0 ; i < sizeof(unsigned long int) ; ++i )
    {
        lhs_uchars[i] = (lhs >> (i*8)) & 255;
        rhs_uchars[i] = (rhs >> (i*8)) & 255;
    }
    prod_ulongs[0] = 0;
    for( i = 0 ; i < 2*sizeof(unsigned long int) ; ++i )
    {
        for( j = 0 ; j < sizeof(unsigned long int) ; ++j )
        {
            k = i-j;
            if( k < 0 || k >= sizeof(unsigned long int) )
            {
                continue;
            }
            prod_ulongs[i] = prod_ulongs[i] + lhs_uchars[j]* rhs_uchars[k];
        }
        prod_ulongs[i+1] = prod_ulongs[i] >> 8;
    }

    *product = (unsigned long int)(lhs * rhs);
    
    *carry = 0;
    for( i = 0 ; i < sizeof(unsigned long int) ; ++i )
    {
        *carry = *carry | ((prod_ulongs[sizeof(unsigned long int)+i] & 255) << (8*i));
    }

    return 1;
}

/**
 * bi_add
 * Add two bigintegers. This function is oblivious of the arguments'
 * signs and furthermore assumes that |lhs| >= |rhs|.
 */
int bi_add_ignoresign( bi * res, bi lhs, bi rhs )
{
    unsigned long int * data;
    unsigned long int carry, temp1, temp2;
    int max, min, i;

    max = lhs.num_limbs;
    min = rhs.num_limbs;

    if( min > max )
    {
        return bi_add_ignoresign(res, rhs, lhs);
    }


    data = malloc(sizeof(unsigned long int) * (max + 1));
    for( i = 0 ; i < max+1 ; ++i )
    {
        data[i] = 0;
    }

    carry = 0;
    for( i = 0 ; i < min ; ++i )
    {
        temp2 = carry;
        carry = 0;
        bi_add_limbs(&temp1, &carry, temp2, lhs.data[i]);
        bi_add_limbs(&data[i], &carry, rhs.data[i], temp1);
    }
    for( i = min ; i < max ; ++i )
    {
        temp1 = 0;
        bi_add_limbs(&data[i], &temp1, lhs.data[i], carry);
        carry = temp1;
    }
    if( carry != 0 )
    {
        data[max] = carry;
    }
    else
    {
        data[max] = 0;
    }

    free(res->data);
    res->data = data;
    if( res->data[max] == 0 )
    {
        res->num_limbs = max;
    }
    else
    {
        res->num_limbs = max + 1;
    }
    res->sign = 1;


    return 1;
}

/**
 * bi_add
 * Adds two big integers.
 * Internally, this function merely decides which subprocedure to
 * call, depending on the sign and relative absolute magnitude (ie,
 * whether |lhs| > |rhs| or |lhs| < |rhs|).
 */
int bi_add( bi * res, bi lhs, bi rhs )
{
    int cmp;

    if( bi_is_zero(lhs) )
    {
        bi_copy(res, rhs);
        return 1;
    }
    else if( bi_is_zero(rhs) )
    {
        bi_copy(res, lhs);
        return 1;
    }

    cmp =  bi_compare_absolute(lhs, rhs);

    if( lhs.sign == 1 && rhs.sign == 1 )
    {
        if( cmp < 0 )
        {
            return bi_add_ignoresign(res, rhs, lhs);
        }
        else
        {
            return bi_add_ignoresign(res, lhs, rhs);
        }
    }

    else if( lhs.sign == -1 && rhs.sign == -1 )
    {
        if( cmp < 0 )
        {
            bi_add_ignoresign(res, rhs, lhs);
        }
        else
        {
            bi_add_ignoresign(res, lhs, rhs);
        }
        bi_negate(res);
        return 1;
    }

    else if( lhs.sign == 1 && rhs.sign == -1 )
    {
        if( cmp < 0 )
        {
            bi_subtract_ignoresign(res, rhs, lhs);
            bi_negate(res);
            return 1;
        }
        else if( cmp > 0 )
        {
            bi_subtract_ignoresign(res, lhs, rhs);
            return 1;
        }
        else
        {
            bi_zero(res);
            return 1;
        }
    }

    else if( lhs.sign == -1 && rhs.sign == 1 )
    {
        cmp =  bi_compare_absolute(lhs, rhs);
        if( cmp < 0 )
        {
            bi_subtract_ignoresign(res, rhs, lhs);
            return 1;
        }
        else if( cmp > 0 )
        {
            bi_subtract_ignoresign(res, lhs, rhs);
            bi_negate(res);
            return 1;
        }
        else
        {
            bi_zero(res);
            return 1;
        }
    }
    else
    {
        return 1;
    }
}

/**
 * bi_subtract_ignoresign
 * Subtract two integers. This function ignores the operands' signs
 * and furthermore assumes |lhs| > |rhs|.
 */
int bi_subtract_ignoresign( bi * res, bi lhs, bi rhs )
{
    int min, max;
    unsigned long int temp, borrow;
    unsigned long int * borrows;
    int i, j;
    unsigned long int * data;

    max = lhs.num_limbs;
    min = rhs.num_limbs;

    if( max < min )
    {
        for( min = min ; min > 1 ; --min )
        {
            if( rhs.data[min-1] != 0 )
            {
                break;
            }
        }
    }

    res->num_limbs = max;
    data = malloc(sizeof(unsigned long int)*max);
    for( i = 0 ; i < res->num_limbs ; ++i )
    {
        data[i] = 0;
    }
    res->sign = 1;

    borrows = malloc(sizeof(unsigned long int)*(max+1));
    for( i = 0 ; i < max+1 ; ++i )
    {
        borrows[i] = 0;
    }

    for( i = 0 ; i < min ; ++i )
    {
        bi_subtract_limbs(&data[i], &borrows[i+1], lhs.data[i], rhs.data[i]);
    }
    for( i = min ; i < max ; ++i )
    {
        data[i] = lhs.data[i];
    }

    for( i = 1 ; i < max ; ++i )
    {
        temp = data[i];
        bi_subtract_limbs(&data[i], &borrows[i+1], temp, borrows[i]);
    }

    for( res->num_limbs = max ; res->num_limbs > 1 ; --res->num_limbs )
    {
        if( data[res->num_limbs-1] != 0 )
        {
            break;
        }
    }

    free(borrows);
    free(res->data);
    res->data = data;

    return 1;
}

/**
 * bi_subtract
 * Subtract one big integer from another.
 * Internally, this function just decides which subprocedure to call.
 */
int bi_subtract( bi * res, bi lhs, bi rhs )
{
    int cmp;
    cmp = bi_compare_absolute(lhs, rhs);

    if( lhs.sign == 1 && rhs.sign == 1 )
    {
        if( cmp > 0 )
        {
            return bi_subtract_ignoresign(res, lhs, rhs);
        }
        else if( cmp < 0 )
        {
            bi_subtract_ignoresign(res, rhs, lhs);
            bi_negate(res);
            return 1;
        }
        else
        {
            bi_zero(res);
            return 1;
        }
    }
    else if( lhs.sign == -1 && rhs.sign == -1 )
    {
        if( cmp > 0 )
        {
            bi_subtract_ignoresign(res, lhs, rhs);
            bi_negate(res);
            return 1;
        }
        else if( cmp < 0 )
        {
            return bi_subtract_ignoresign(res, rhs, lhs);
        }
        else
        {
            bi_zero(res);
            return 1;
        }
    }
    else if( lhs.sign == 1 && rhs.sign == -1 )
    {
        if( cmp >= 0 )
        {
            return bi_add_ignoresign(res, lhs, rhs);
        }
        else
        {
            return bi_add_ignoresign(res, rhs, lhs);
        }
    }
    else if( lhs.sign == -1 && rhs.sign == 1 )
    {
        if( cmp <= 0 )
        {
            bi_add_ignoresign(res, rhs, lhs);
            bi_negate(res);
            return 1;
        }
        else
        {
            bi_add_ignoresign(res, lhs, rhs);
            bi_negate(res);
            return 1;
        }
    }
    else
    {
        return 1;
    }
}


/**
 * bi_negate
 * Flip the sign of the big integer.
 */
int bi_negate( bi * op )
{
    op->sign = -op->sign;
    return 1;
}

/**
 * bi_multiply_shift_limb
 * Multiply a big integer by a limb and shift it left by a given
 * number of limbs.
 * @params:
 *  * res : big integer -- stores the product
 *  * op : big integer -- the operand to be multiplied by the limb
 *  * limb : unsigned long int -- the operand to multiply with
 *  * int shamt : number of limbs to shift by
 */
int bi_multiply_shift_limb( bi * res, bi op, unsigned long int limb, int shamt )
{
    unsigned long int carry, prod;
    int i;
    bi hi, lo;

    hi = bi_init((op.num_limbs+1) * sizeof(unsigned long int) * 8);
    lo = bi_init((op.num_limbs) * sizeof(unsigned long int) * 8);

    for( i = 0 ; i < op.num_limbs ; ++i )
    {
        hi.data[i] = 0;
        lo.data[i] = 0;
    }
    hi.data[op.num_limbs] = 0;

    for( i = 0 ; i < op.num_limbs ; ++i )
    {
        bi_multiply_limbs(&prod, &carry, limb, op.data[i]);
        lo.data[i] = prod;
        hi.data[i+1] = carry;
    }

    bi_add(&hi, hi, lo);

    if( res->num_limbs < hi.num_limbs + shamt )
    {
        res->num_limbs = hi.num_limbs + shamt;
        free(res->data);
        res->data = malloc(res->num_limbs * sizeof(unsigned long int));
    }

    for( i = 0 ; i < shamt ; ++i )
    {
        res->data[i] = 0;
    }

    for( i = 0 ; i < hi.num_limbs ; ++i )
    {
        res->data[i+shamt] = hi.data[i];
    }

    bi_destroy(hi);
    bi_destroy(lo);

    res->sign = op.sign;

    return 1;
}

/**
 * bi_multiply
 * Multiply two big integers.
 */
int bi_multiply( bi * res, bi lhs, bi rhs )
{
    bi curr;
    unsigned int i;

    if( rhs.num_limbs > lhs.num_limbs )
    {
        return bi_multiply(res, rhs, lhs);
    }

    curr = bi_init((lhs.num_limbs + rhs.num_limbs + 1)*sizeof(unsigned long int)*8);
    bi_zero(res);

    for( i = 0 ; i < rhs.num_limbs ; ++i )
    {
        bi_multiply_shift_limb(&curr, lhs, rhs.data[i], i);
        bi_add(res, *res, curr);
    }

    res->sign = lhs.sign * rhs.sign;

    bi_destroy(curr);

    return 1;
}

/**
 * bi_divide
 * Divide the numerator by the denominator and store the quotient in
 * quo and the remainder in rem. The remainder is positive while the
 * sign of the quotient adapts to the signs of the numerator and the
 * denominator.
 * @params:
 *  * quo, rem : pointers to big integrs -- containers for quotient
 *    and remainder, respectively
 *  * numerator, denominator : big integers -- the integers to divide
 */
int bi_divide( bi * quo, bi * rem, bi numerator, bi denominator )
{
    int bitsize_difference;
    int bitsize_numerator;
    int bitsize_denominator;
    int i;
    int cmp;
    bi temp;
    bi shifted;
    unsigned char * randomness;
    bi positive_numerator, positive_denominator, negative_quotient, negative_remainder;
    bi one;


    if( numerator.sign == -1 && denominator.sign == -1 )
    {
        positive_numerator = bi_init(0);
        positive_denominator = bi_init(0);
        negative_quotient = bi_init(0);
        negative_remainder = bi_init(0);

        bi_copy(&positive_numerator, numerator);
        bi_copy(&positive_denominator, denominator);
        bi_negate(&positive_numerator);
        bi_negate(&positive_denominator);
        bi_divide(quo, &negative_remainder, positive_numerator, positive_denominator);

        bi_negate(&negative_remainder);
        bi_subtract(rem, negative_remainder, denominator);
        
        one = bi_init(0);
        bi_one(&one);
        temp = bi_init(0);
        bi_copy(&temp, *quo);
        bi_add(quo, temp, one);

        bi_destroy(temp);
        bi_destroy(one);
        
        bi_destroy(positive_numerator);
        bi_destroy(positive_denominator);
        bi_destroy(negative_quotient);
        bi_destroy(negative_remainder);
        return 1;
    }
    else if( numerator.sign == -1 )
    {
        positive_numerator = bi_init(0);
        negative_quotient = bi_init(0);
        negative_remainder = bi_init(0);
        bi_copy(&positive_numerator, numerator);
        bi_negate(&positive_numerator);

        bi_divide(&negative_quotient, &negative_remainder, positive_numerator, denominator);
        bi_negate(&negative_quotient);
        bi_copy(quo, negative_quotient);

        bi_negate(&negative_remainder);
        bi_add(rem, negative_remainder, denominator);

        one = bi_init(0);
        bi_one(&one);
        temp = bi_init(0);
        bi_copy(&temp, *quo);
        bi_subtract(quo, temp, one);
        bi_destroy(temp);
        bi_destroy(one);

        bi_destroy(positive_numerator);
        bi_destroy(negative_remainder);
        bi_destroy(negative_quotient);
        return 1;
    }
    else if( denominator.sign == -1 )
    {
        positive_denominator = bi_init(0);
        bi_copy(&positive_denominator, denominator);
        bi_negate(&positive_denominator);

        bi_divide(quo, rem, numerator, positive_denominator);
        bi_negate(quo);

        bi_destroy(positive_denominator);
        return 1;
    }
    else if( bi_compare(numerator, denominator) < 0 )
    {
        bi_copy(rem, numerator);
        bi_zero(quo);
        return 1;
    }

    if( bi_is_one(denominator) )
    {
        bi_copy(quo, denominator);
        bi_zero(rem);
        return 1;
    }

    if( bi_is_zero(denominator) )
    {
        bi_copy(rem, denominator);
        return 1;
    }

    bitsize_numerator = bi_bitsize(numerator);
    bitsize_denominator = bi_bitsize(denominator);
    bitsize_difference = bitsize_numerator - bitsize_denominator;
    bi_copy(rem, numerator);
    shifted = bi_init(0);

    bi_shift_left(&shifted, denominator, bitsize_difference);

    if( quo->num_limbs < (bitsize_difference + sizeof(unsigned long int)*8 - 1)/(sizeof(unsigned long int)*8) )
    {
        quo->num_limbs = (bitsize_difference + sizeof(unsigned long int)*8 - 1)/(sizeof(unsigned long int)*8);
        free(quo->data);
        quo->data = malloc(quo->num_limbs * sizeof(unsigned long int));
    }
    for( i = 0 ; i < quo->num_limbs ; ++i )
    {
        quo->data[i] = 0;
    }
    quo->sign = 1;

    for( i = bitsize_difference ; i >= 0 ; --i )
    {
        cmp = bi_compare(*rem, shifted);
        if( cmp < 0 )
        {
            bi_setbit(quo, i, 0 );
        }
        else if( cmp > 0 )
        {
            bi_setbit(quo, i, 1);
            bi_subtract(rem, *rem, shifted);
        }
        else
        {
            bi_setbit(quo, i, 1);
            bi_subtract(rem, *rem, shifted);
            bi_destroy(shifted);
            return 1;
        }
        bi_shift_right(&shifted, 1);
    }

    bi_destroy(shifted);

    return 1;
}

/**
 * bi_setbit
 * Set the ith bit to one or zero.
 */
int bi_setbit( bi * integer, int bit_index, char bit )
{
    unsigned int limb_index;
    unsigned int rem;
    unsigned long int mask;
    int num_bits;
    unsigned long int * data;
    int i;

    num_bits = integer->num_limbs * sizeof(unsigned long int) * 8;
    limb_index = bit_index / (sizeof(unsigned long int) * 8);
    rem = bit_index % (sizeof(unsigned long int) * 8);

    if( num_bits <= bit_index )
    {
        data = malloc(sizeof(unsigned long int) * (limb_index + 1));
        for( i = 0 ; i < integer->num_limbs ; ++i )
        {
            data[i] = integer->data[i];
        }
        for( i = integer->num_limbs ; i < limb_index + 1 ; ++i )
        {
            data[i] = 0;
        }
        free(integer->data);
        integer->data = data;
        integer->num_limbs = limb_index + 1;
    }

    mask = (unsigned long int)1 << bit_index;

    if( bit == 1 )
    {
        integer->data[limb_index] = integer->data[limb_index] | mask;
    }
    else if( bit == 0 )
    {
        integer->data[limb_index] = integer->data[limb_index] & (~mask);
    }

    return 1;
}

/**
 * bi_compare_absolute
 * Compare two integers, ignoring their signs.
 * @return
 * * cmp : integer -- same sign as |lhs|-|rhs|
 */
int bi_compare_absolute( bi lhs, bi rhs )
{
    int i, j;

    if( rhs.num_limbs > lhs.num_limbs )
    {
        return -bi_compare_absolute(rhs, lhs);
    }

    for( i = lhs.num_limbs - 1 ; i >= rhs.num_limbs ; --i )
    {
        if( lhs.data[i] != 0 )
        {
            return 1;
        }
    }
    for( i = rhs.num_limbs - 1 ; i >= 0 ; --i )
    {
        if( lhs.data[i] > rhs.data[i] )
        {
            return 1;
        }
        else if( lhs.data[i] < rhs.data[i] )
        {
            return -1;
        }
    }
    return 0;
}

/**
 * bi_compare
 * Find out which one of two big integers is larger.
 * @returns:
 *  * 1 if lhs > rhs; -1 if lhs < rhs; 0 if lhs = rhs
 */
int bi_compare( bi lhs, bi rhs )
{
    int abscmp;

    abscmp = bi_compare_absolute(lhs,rhs);

    if( abscmp == 0 )
    {
        if( lhs.sign == rhs.sign )
        {
            return 0;
        }
        if( bi_is_zero(lhs) == 1 && bi_is_zero(rhs) == 1 )
        {
            return 0;
        }
        if( lhs.sign > 0 )
        {
            return 1;
        }
        return -1;
    }
    else if( abscmp > 0 )
    {
        if( lhs.sign > 0 )
        {
            return 1;
        }
        return -1;
    }
    if( rhs.sign > 0 )
    {
        return -1;
    }
    return 1;
}

/**
 * bi_is_one
 * Decide if this big integer is equal to one or not.
 * @returns:
 *  * 1 if equal to one; 0 otherwise
 */
int bi_is_one( bi test )
{
    int i;
    if( test.sign != 1 )
    {
        return 0;
    }
    for( i = test.num_limbs-1 ; i >= 1 ; --i )
    {
        if( test.data[i] != 0 )
        {
            return 0;
        }
    }
    if( test.data[0] != 1 )
    {
        return 0;
    }
    return 1;
}

/**
 * bi_is_zero
 * Decide if this big integer is equal to zero or not.
 * returns:
 *  * 1 if equal to zero; 0 otherwise
 */
int bi_is_zero( bi test )
{
    int i;
    int ret;
    ret = 1;
    for( i = test.num_limbs-1 ; i >= 0 ; --i )
    {
        if( test.data[i] != 0 )
        {
            ret = 0;
        }
    }
    return ret;
}

/**
 * bi_bitsize
 * Determine the size of this integer in terms of number of bit.
 */
int bi_bitsize( bi integer )
{
    int bitsize;
    int i;
    unsigned long int msd;

    for( i = integer.num_limbs - 1 ; i >= 0 ; --i )
    {
        if( integer.data[i] != 0 )
        {
            break;
        }
    }
    if( i == -1 )
    {
        return 0;
    }

    bitsize = i * sizeof(unsigned long int) * 8;

    msd = integer.data[i];
    for( i = 0 ; i < sizeof(unsigned long int) * 8 ; ++i )
    {
        bitsize = bitsize + 1;
        msd = msd >> 1;
        if( msd == 0 )
        {
            break;
        }
    }

    return bitsize;
}

/**
 * bi_shift_left
 * Shift the given integer to the left as if multiplied by the given
 * power of two.
 * @params:
 *  * dest : pointer to big integer -- container for the shifted
 *    version
 *  * source : big integer -- the integer to copy and shift
 *  * shamt : int -- amount to shift by
 */
int bi_shift_left( bi * dest, bi source, int shamt )
{
    int num_bits;
    int num_limbs;
    int i;
    int shamtmod, shamtdiv;
    unsigned long int  * data;

    if( shamt < 0 )
    {
        return 1;
    }

    if( shamt == 0 )
    {
        return bi_copy(dest, source);
    }

    num_bits = bi_bitsize(source);
    shamtmod = shamt % (sizeof(unsigned long int)*8);
    shamtdiv = shamt / (sizeof(unsigned long int)*8);
    num_limbs = (shamt + num_bits + sizeof(unsigned long int)*8 - 1) / (sizeof(unsigned long int)*8);

        dest->num_limbs = num_limbs;
        data = malloc(sizeof(unsigned long int) * num_limbs);
    dest->sign = source.sign;
    for( i = 0 ; i < dest->num_limbs ; ++i )
    {
        data[i] = 0;
    }
    for( i = 0 ; i < source.num_limbs && i < dest->num_limbs ; ++i )
    {
        data[shamtdiv+i] = source.data[i] << shamtmod;
    }
    for( i = 0 ; i < source.num_limbs && i+shamtdiv+1 < dest->num_limbs ; ++i )
    {
        data[shamtdiv+i+1] = data[shamtdiv+1+i] | (source.data[i] >> (sizeof(unsigned long int)*8 - shamtmod));
    }

    for( dest->num_limbs = dest->num_limbs ; dest->num_limbs > 0 ; --dest->num_limbs )
    {
        if( data[dest->num_limbs-1] != 0 )
        {
            break;
        }
    }

    free(dest->data);
    dest->data = data;

    return 1;
}

/**
 * bi_shift_right
 * Shift the given to the right, and drop the rightmost bits.
 */
int bi_shift_right( bi * op, int shamt )
{
    int shamtdiv, shamtmod;
    int i;

    shamtdiv = shamt / (sizeof(unsigned long int)*8);
    shamtmod = shamt % (sizeof(unsigned long int)*8);

    for( i = 0 ; i < op->num_limbs - shamtdiv ; ++i )
    {
        op->data[i] = op->data[i+shamtdiv];
    }
    for( i = op->num_limbs - shamtdiv ; i < op->num_limbs ; ++i )
    {
        op->data[i] = 0;
    }
    op->num_limbs = op->num_limbs - shamtdiv;

    for( i = 0 ; i < op->num_limbs-1 ; ++i )
    {
        op->data[i] = (op->data[i] >> shamtmod) | (op->data[i+1] << (sizeof(unsigned long int)*8 - shamtmod));
    }
    op->data[op->num_limbs-1] = op->data[op->num_limbs-1] >> shamtmod;

    return 1;
}

/**
 * bi_random
 * Generate a random integer from the given string of random bytes.
 */
int bi_random( bi * dest, unsigned int num_bits, unsigned char * randomness )
{
    unsigned int num_limbs;
    int rem;
    int i;
    num_limbs = (num_bits + sizeof(unsigned long int)*8 - 1) / (sizeof(unsigned long int)*8);
   
    if( dest->num_limbs != num_limbs )
    { 
        free(dest->data);
        dest->data = malloc(num_limbs * sizeof(unsigned long int));
        dest->num_limbs = num_limbs;
    }
    dest->sign = 1;

    for( i = 0 ; i < dest->num_limbs ; ++i )
    {
        dest->data[i] = ((unsigned long int*)randomness)[i];
    }

    rem = num_bits % (sizeof(unsigned long int) * 8);
    if( rem == 0 )
    {
        rem = rem + sizeof(unsigned long int) * 8;
    }
    dest->data[num_limbs-1] = dest->data[num_limbs-1] >> (sizeof(unsigned long int) * 8 - rem);

    return 1;
}

/**
 * bi_modulo
 * Compute the remainder of the given integer after division by the
 * modulus.
 */
int bi_modulo( bi * res, bi integer, bi modulus )
{
    bi quotient;
    quotient = bi_init(0);
    bi_divide(&quotient, res, integer, modulus);
    bi_destroy(quotient);
    return 1;
}

/**
 * bi_xgcd
 * Compute the greatest common divisor and Bezout coefficients using
 * the extended Euclidean algorithm. In other words, find (x,y,g)
 * such that  x*a + y*b = g = gcd(a,b).
 * @returns:
 *  x : big integer -- Bezout coefficient for a
 *  y : big integer -- Bezout coefficient for b
 *  g : bit integer -- greatest common divisor
 */
int bi_xgcd( bi * x, bi * y, bi * g, bi a, bi b )
{
    bi r, s, t, old_r, old_s, old_t, quotient, temp, temp2;

    r = bi_init(0);
    s = bi_init(0);
    t = bi_init(0);
    quotient = bi_init(0);
    old_r = bi_init(0);
    old_s = bi_init(0);
    old_t = bi_init(0);
    temp = bi_init(0);
    temp2 = bi_init(0);

    bi_zero(&s);
    bi_one(&old_s);
    bi_one(&t);
    bi_zero(&old_t);
    bi_copy(&r, b);
    bi_copy(&old_r, a);

    while( !bi_is_zero(r) )
    {
        bi_divide(&quotient, &temp, old_r, r);

        bi_copy(&old_r, r);
        bi_copy(&r, temp);

        bi_copy(&temp2, old_s);
        bi_copy(&old_s, s);
        bi_multiply(&temp, quotient, s);
        bi_subtract(&s, temp2, temp);

        bi_copy(&temp2, old_t);
        bi_copy(&old_t, t);
        bi_multiply(&temp, quotient, t);
        bi_subtract(&t, temp2, temp);
    }

    bi_copy(x, old_s);
    bi_copy(y, old_t);
    bi_copy(g, old_r);

    bi_destroy(s);
    bi_destroy(t);
    bi_destroy(r);
    bi_destroy(old_s);
    bi_destroy(old_t);
    bi_destroy(old_r);
    bi_destroy(quotient);
    bi_destroy(temp);
    bi_destroy(temp2);

    return 1;
}

/**
 * bi_modcent
 * Compute the representative of the given integer modulo the modulus
 * inside the range [-floor(modulus/2):floor(modulus/2)].
 */
int bi_modcent( bi * res, bi integer, bi modulus )
{
    bi modover2;
    modover2 = bi_init(0);
    bi_copy(&modover2, modulus);
    bi_shift_right(&modover2, 1);

    bi_modulo(res, integer, modulus);

    if( bi_compare(integer, modover2) > 0 )
    {
        bi_copy(&modover2, *res);
        bi_subtract(res, modover2, modulus);
    }

    bi_destroy(modover2);

    return 1;
}

/**
 * bi_naf
 * Convert the given integer to non-adjacent form (NAF), namely to a
 * pair of integers (pos, neg) such that source = pos - neg and both
 * pos and neg have low hamming weight.
 */
int bi_naf( bi * respos, bi * resneg, bi source )
{
    int i, j;

    bi oneshift;

    oneshift = bi_init(1);
    bi_one(&oneshift);

    bi_copy(respos, source);
    bi_zero(resneg);

    for( i = 0 ; i < bi_bitsize(source) - 4 ; ++i )
    {
        if( bi_getbit(*respos, i) == 1 && bi_getbit(*respos, i+1) == 1  && bi_getbit(*respos, i+2) == 1 && bi_getbit(*respos, i+3) == 1 && bi_getbit(*respos, i+4) == 1 )
        {
            bi_add(respos, *respos, oneshift);
            bi_add(resneg, *resneg, oneshift);
            i = i + 4;
        }

        bi_shift_left(&oneshift, oneshift, 1);
    }

    bi_destroy(oneshift);
    return 1;
}

/**
 * bi_modexp
 * Raise the base to the given power but compute only the residue of
 * that function modulo the modulus. Uses the square-and-multiply
 * routine.
 */
int bi_modexp( bi * res, bi base, bi power, bi modulus )
{
    int bitsize;
    int i;
    bi temp, temp2;
    bi posones, negones;
    bi inverse;

    if( power.sign == -1 )
    {
        inverse = bi_init(0);
        bi_modinverse(&inverse, base, modulus);
        temp2 = bi_init(0);
        bi_copy(&temp2, power);
        bi_negate(&temp2);
        bi_modexp(res, inverse, temp2, modulus);
        bi_destroy(inverse);
        bi_destroy(temp2);
        return 1;
    }

    temp = bi_init(0);
    posones = bi_init(0);
    negones = bi_init(0);
    bi_naf(&posones, &negones, power);
    inverse = bi_init(0);
    bi_modinverse(&inverse, base, modulus);
    bi_one(res);

    bitsize = bi_bitsize(power);
    for( i = bitsize ; i >= 0 ; --i )
    {
        bi_multiply(&temp, *res, *res);
        bi_modulo(res, temp, modulus);

        /*
        if( bi_getbit(power, i) == 1 )
        {
            bi_multiply(&temp, *res, base);
            bi_modulo(res, temp, modulus);
        }
        */

        
        if( bi_getbit(posones, i) == 1 )
        {
            bi_multiply(&temp, *res, base);
            bi_modulo(res, temp, modulus);
        }
        

        
        if( bi_getbit(negones, i) == 1 )
        {
            bi_multiply(&temp, *res, inverse);
            bi_modulo(res, temp, modulus);
        }
        
    }

    bi_destroy(temp);
    bi_destroy(inverse);
    bi_destroy(posones);
    bi_destroy(negones);

    return 1;
}

/**
 * bi_modinverse
 * Compute the multiplicative inverse of some element modulo the
 * modulus. If no multiplicative inverse exists, this function will
 * return zero and offers no guarantee on the contents of res.
 */
int bi_modinverse( bi * res, bi element, bi modulus )
{
    int gcd_return;
    bi g, y;
    g = bi_init(0);
    y = bi_init(1);

    bi_xgcd(res, &y, &g, element, modulus);

    gcd_return = (g.data[0] == 1);

    bi_destroy(g);
    bi_destroy(y);

    return gcd_return;
}

/**
 * bi_miller_rabin_trial
 * Perform a single Miller-Rabin trial of the integer with base as
 * test.
 */
int bi_miller_rabin_trial( bi integer, bi base )
{
    bi nm1, one, two;
    bi d, x, temp;
    int i, r;

    d = bi_init(0);
    nm1 = bi_init(0);
    one = bi_cast(1);
    two = bi_cast(2);
    bi_subtract(&nm1, integer, one);

    for( r = 0 ; r < bi_bitsize(nm1) ; ++r )
    {
        if( bi_getbit(nm1, r) == 1 )
        {
            break;
        }
    }
    bi_copy(&d, nm1);
    bi_shift_right(&d, r-1);

    x = bi_init(0);
    bi_modexp(&x, base, d, integer);

    if( bi_is_one(x) == 1 )
    {
        bi_destroy(x);
        bi_destroy(nm1);
        bi_destroy(d);
        bi_destroy(one);
        bi_destroy(two);
        return 1;
    }

    if( bi_compare(x, nm1) == 0 )
    {
        bi_destroy(x);
        bi_destroy(nm1);
        bi_destroy(d);
        bi_destroy(one);
        bi_destroy(two);
        return 1;
    }

    temp = bi_init(0);
    for( i = 0 ; i < r-1 ; ++i )
    {
        bi_copy(&temp, x);
        bi_modexp(&x, temp, two, integer);
        if( bi_is_one(x) )
        {
            bi_destroy(x);
            bi_destroy(nm1);
            bi_destroy(d);
            bi_destroy(one);
            bi_destroy(two);
            return 0;
        }
    }

    bi_destroy(x);
    bi_destroy(nm1);
    bi_destroy(d);
    bi_destroy(one);
    bi_destroy(two);
    bi_destroy(temp);
    return 1;
}

/**
 * is_prime
 * Test probabilistically if the given integer is a prime or not.
 * Use the Miller-Rabin primality test, which passes composites with
 * probability at most 1/4 per trial for the number of trials equal
 * to certainty.
 * @params:
 *  * integer : big integer -- determine if this is a prime or not
 *  * random_ints : pointer to list of unsigned long int -- list of
 *    random challenge numbers; this list contains certainty
 *    elements.
 *  * certainty : int -- number of trials to run
 * @returns:
 *  * 1 if integer is prime; 0 otherwise. With probability
 *    no more than (1/4)^certainty this algorithm return 1 for a
 *    composite number.
 */
int bi_is_prime( bi integer, unsigned long int * random_ints, int certainty )
{
    bi base, mod6;
    int i;
    int composite;

    /* negative numbers can't be prime */
    if( integer.sign == -1 )
    {
        return 0;
    }

    /* filter out small integers */
    if( integer.num_limbs == 1 )
    {
        if( integer.data[0] == 1 )
        {
            return 0;
        }
        if( integer.data[0] == 2 || integer.data[0] == 3 )
        {
            return 1;
        }
    }

    /* make sure the integer is 1 or -1 modulo 6 */
    base = bi_cast(6);
    mod6 = bi_init(3);
    bi_modulo(&mod6, integer, base);
    bi_destroy(base);
    if( mod6.data[0] != 1 && mod6.data[0] != 5 )
    {
        bi_destroy(mod6);
        return 0;
    }
    bi_destroy(mod6);

    /* run Miller-Rabin trials */
    composite = 0;
    for( i = 0 ; i < certainty ; ++i )
    {
        base = bi_cast(random_ints[i]);
        composite = bi_miller_rabin_trial(integer, base);
        bi_destroy(base);
        if( composite == 0 )
        {
            return 0;
        }
    }

    return 1;
}

/**
 * bi_getbit
 * Get the ith bit of the given integer.
 */
int bi_getbit( bi integer, int bit_index )
{
    int limb_index;
    int shamt;
    unsigned long int mask;

    shamt = bit_index % (sizeof(unsigned long int) * 8);
    limb_index = bit_index / (sizeof(unsigned long int) * 8);

    mask = (unsigned long int)(1) << shamt;

    return (integer.data[limb_index] & mask) >> shamt;
}

/**
 * bi_print
 * Send the decimal representation of this integer to stdout.
 */
int bi_print( bi integer )
{
    bi quo, rem1, rem2, ten;
    int i;
    char * expansion;
    int numdigits;
    int bitsize;

    bitsize = bi_bitsize(integer);
    expansion = malloc(bi_bitsize(integer));
    numdigits = 0;
    quo = bi_init(0);
    rem1 = bi_init(0);
    rem2 = bi_init(0);
    ten = bi_cast(10);
    bi_copy(&rem2, integer);

    if( integer.sign == -1 )
    {
        printf("-");
    }


    for( i = 0 ; i < bitsize && bi_is_zero(rem2) == 0 ; ++i )
    {
        bi_divide(&quo, &rem1, rem2, ten);
        bi_copy(&rem2, rem1);
        expansion[numdigits++] = rem1.data[0];
    }

    for( numdigits = numdigits - 1 ; numdigits >= 0 ; --numdigits )
    {
        printf("%i", expansion[numdigits]);
    }

    bi_destroy(rem1);
    bi_destroy(rem2);
    bi_destroy(quo);
    bi_destroy(ten);
    free(expansion);
}

/**
 * bi_print_bitstring
 * Print the integer as a bitstring and send it to stdout.
 */
int bi_print_bitstring( bi integer )
{
    int num_bits;

    num_bits = bi_bitsize(integer);
    if( integer.sign == -1 )
    {
        printf("-");
    }
    if( num_bits == 0 )
    {
        printf("0");
    }
    for( num_bits = num_bits - 1 ; num_bits >= 0 ; --num_bits )
    {
        printf("%i", bi_getbit(integer, num_bits));
    }
    return 1;
}

