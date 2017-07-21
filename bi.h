#ifndef BI_H
#define BI_H

typedef struct
{
    char sign;
    unsigned int num_limbs;
    unsigned long int * data;
} bi;

bi bi_init( int number_of_limbs );
bi bi_cast( long int integer );
bi bi_cast_unsigned( unsigned long int integer );
bi bi_cast_bytestring( int num_bytes, unsigned char * bytes );
bi bi_local( char sign, int num_limbs, unsigned long int * limbs );
bi bi_clone( bi source );
int bi_destroy( bi integer );

int bi_one( bi * integer );
int bi_zero( bi * integer );
int bi_add_limbs( unsigned long int * sum, unsigned long int * carry, unsigned long int lhs, unsigned long int rhs );
int bi_subtract_limbs( unsigned long int * difference, unsigned long int * borrow, unsigned long int lhs, unsigned long int rhs );
int bi_multiply_limbs( unsigned long int * product, unsigned long int * carry, unsigned long int lhs, unsigned long int rhs );
int bi_add_ignoresign( bi * res, bi lhs, bi rhs );
int bi_add( bi * res, bi lhs, bi rhs );
int bi_increase( bi * lhs, bi rhs );
int bi_subtract( bi * res, bi lhs, bi rhs );
int bi_subtract_ignoresign( bi * res, bi lhs, bi rhs );
int bi_multiply( bi * res, bi lhs, bi rhs );
int bi_divide( bi * quo, bi * rem, bi numerator, bi denominator );
int bi_xgcd( bi * g, bi * a, bi * b, bi x, bi y );
int bi_negate( bi * integer );
int bi_compare( bi lhs, bi rhs );
int bi_compare_absolute( bi lhs, bi rhs );
int bi_is_zero( bi test );
int bi_is_one( bi test );
int bi_modulo( bi * res, bi integer, bi modulus );
int bi_modcent( bi * res, bi integer, bi modulus );
int bi_modinverse( bi * res, bi integer, bi modulus );
int bi_naf( bi * respos, bi * resneg, bi source );
int bi_modexp( bi * res, bi base, bi power, bi modulus );
int bi_bitsize( bi integer );
int bi_setbit( bi * integer, int bit_index, char bit );
int bi_getbit( bi integer, int bit_index );
int bi_random( bi * dest, unsigned int num_bits, unsigned char * randomness );
int bi_is_prime( bi integer, unsigned long int * random_ints, int certainty );
int bi_print( bi integer );
int bi_shift_left( bi * dest, bi src, int shamt );
int bi_shift_right( bi * op, int shamt );
unsigned char bi_getbyte( bi integer, int byte_index );
int bi_serialize( unsigned char * bytes, bi integer );
int bi_deserialize( bi * integer, unsigned char * bytes, int num_bytes );

#endif

