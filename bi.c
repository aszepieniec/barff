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
 *  * carry : pointer to unsigend long int -- to be incremented by
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
            if( k < 0 || j < 0 || k >= sizeof(unsigned long int) || j >= sizeof(unsigned long int) )
            {
                continue;
            }
            prod_ulongs[i] = prod_ulongs[i] + lhs_uchars[j]* rhs_uchars[k];
        }
        //printf("prod_ulongs[%i] = %u\n", i, prod_ulongs[i]);
        prod_ulongs[i+1] = prod_ulongs[i] >> 8;
    }
    //getchar();

    *product = (unsigned long int)(lhs * rhs);
    
    *carry = 0;
    for( i = 0 ; i < sizeof(unsigned long int) ; ++i )
    {
        *carry = *carry | ((prod_ulongs[sizeof(unsigned long int)+i] & 255) << (8*i));
    }

    //printf(" bin(%lu * %lu - %lu * 2^%i - %lu)\n", lhs, rhs, *carry, sizeof(unsigned long int)*8, *product);

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
    int min, max, i;
    unsigned long int borrow, temp;
    unsigned long int * data;
    unsigned long int * borrows;

    min = rhs.num_limbs;
    max = lhs.num_limbs;

    borrows = malloc((max+1)*sizeof(unsigned long int));
    data = malloc(max * sizeof(unsigned long int));
    for( i = 0 ; i < max ; ++i )
    {
        data[i] = 0;
    }

    borrows[0] = 0;
    for( i = 0 ; i < min ; ++i )
    {
        borrow = 0;
        bi_subtract_limbs(&temp, &borrow, lhs.data[i], borrows[i]);
        bi_subtract_limbs(&data[i], &borrow, temp, rhs.data[i]);
        borrows[i+1] = borrow;
    }
    for( i = min ; i < max ; ++i )
    {
        borrow = 0;
        bi_subtract_limbs(&data[i], &borrow, lhs.data[i], borrows[i]);
        borrows[i+1] = borrow;
    }

    free(res->data);
    res->data = data;
    res->sign = 1;
    for( res->num_limbs = max ; res->num_limbs > 1 ; --res->num_limbs )
    {
        if( res->data[res->num_limbs-1] != 0 )
        {
            break;
        }
    }

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
 * Multiply a big integer by a limb and shift it left bu a given
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
    bi hi, lo, temp;

    hi = bi_init((op.num_limbs+1) * sizeof(unsigned long int) * 8);
    lo = bi_init((op.num_limbs) * sizeof(unsigned long int) * 8);
    temp = bi_init(0);

    for( i = 0 ; i < op.num_limbs ; ++i )
    {
        hi.data[i] = 0;
        lo.data[i] = 0;
    }
    hi.data[op.num_limbs] = 0;

    for( i = 0 ; i < op.num_limbs ; ++i )
    {
        carry = 0;
        bi_multiply_limbs(&prod, &carry, limb, op.data[i]);
        lo.data[i] = prod;
        hi.data[i+1] = carry;
    }

    bi_add(&temp, hi, lo);

    if( res->num_limbs < temp.num_limbs + shamt )
    {
        res->num_limbs = temp.num_limbs + shamt;
        free(res->data);
        res->data = malloc(res->num_limbs * sizeof(unsigned long int));
    }

    for( i = 0 ; i < res->num_limbs ; ++i )
    {
        res->data[i] = 0;
    }

    for( i = 0 ; i < temp.num_limbs ; ++i )
    {
        res->data[i+shamt] = temp.data[i];
    }

    bi_destroy(hi);
    bi_destroy(lo);
    bi_destroy(temp);

    res->sign = op.sign;

    if( limb == 2529350679851345645 )
    {
        printf("int('");
        bi_print_bitstring(op);
        printf("', 2) * %lu * 2^(%i) - int('", limb, sizeof(unsigned long int) * shamt * 8);
        bi_print_bitstring(*res);
        printf("', 2)\n");
    }

    return 1;
}

/**
 * bi_multiply
 * Multiply two big integers.
 */
int bi_multiply( bi * res, bi lhs, bi rhs )
{
    bi curr, temp;
    unsigned int i;

    temp = bi_init(res->num_limbs);
    curr = bi_init((lhs.num_limbs + rhs.num_limbs + 1)*sizeof(unsigned long int)*8);

    for( i = 0 ; i < rhs.num_limbs ; ++i )
    {
        bi_multiply_shift_limb(&curr, lhs, rhs.data[i], i);
        bi_copy(&temp, *res);
        bi_add(res, temp, curr);
    }

    res->sign = lhs.sign * rhs.sign;

    bi_destroy(curr);
    bi_destroy(temp);

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
    bi shifted;
    unsigned char * randomness;


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
    shifted = bi_init(numerator.num_limbs * sizeof(unsigned long int) * 8);

    for( i = bitsize_difference ; i >= 0 ; --i )
    {
        bi_shift_left(&shifted, denominator, i);
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
            return 1;
        }
    }

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

    mask = 1 << bit_index;

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

int bi_xgcd( bi * g, bi * a, bi * b, bi x, bi y );

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
    int num_limbs;
    int i;
    int shamtmod, shamtdiv;

    if( shamt < 0 )
    {
        return 1;
    }

    shamtmod = shamt % (sizeof(unsigned long int)*8);
    shamtdiv = shamt / (sizeof(unsigned long int)*8);
    num_limbs = 1+source.num_limbs+shamtdiv;

    free(dest->data);
    dest->num_limbs = num_limbs;
    dest->data = malloc(sizeof(unsigned long int) * num_limbs);

    for( i = 0 ; i < num_limbs ; ++i )
    {
        dest->data[i] = 0;
    }
    for( i = 0 ; i < source.num_limbs ; ++i )
    {
        dest->data[shamtdiv+i] = source.data[i] << shamtmod;
    }
    for( i = 0 ; i < source.num_limbs ; ++i )
    {
        dest->data[shamtdiv+i+1] = dest->data[shamtdiv+1+i] | (source.data[i] >> (sizeof(unsigned long int)*8 - shamtmod));
    }
    dest->sign = source.sign;

    for( dest->num_limbs = dest->num_limbs ; dest->num_limbs > 0 ; --dest->num_limbs )
    {
        if( dest->data[dest->num_limbs-1] != 0 )
        {
            break;
        }
    }

    return 1;
}

/**
 * bi_random
 * Generate a random integer from the given string of random bytes.
 */
int bi_random( bi * dest, unsigned int num_bits, char * randomness )
{
    unsigned int num_limbs;
    int rem;
    int i;
    num_limbs = (num_bits + sizeof(unsigned long int)*8 - 1) / (sizeof(unsigned long int)*8);
    
    free(dest->data);
    dest->data = malloc(num_limbs * sizeof(unsigned long int));
    dest->num_limbs = num_limbs;
    dest->sign = 1;

    //memcpy((unsigned char*)dest->data, randomness, (num_bits + 7) / 8);
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

int bi_modulo( bi * res, bi integer, bi modulus );
int bi_modcent( bi * res, bi integer, bi modulus );
int bi_naf( bi * respos, bi * resneg, bi source );
int bi_modexp( bi * res, bi base, bi power, bi modulus );
int bi_modinverse( bi * res, bi element, bi modulus );
int bi_isprime( bi integer, char randomseed[SEED_LENGTH], int certainty );

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

