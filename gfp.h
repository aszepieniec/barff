#ifndef GFP
#define GFP

#ifndef BIG
    
    #ifndef GF_PRIME_MODULUS
    #define GF_PRIME_MODULUS 31
    #endif
    
    #ifndef GFP_NUMBITS
    #define GFP_NUMBITS 5
    #endif
    
    #ifndef GFP_NUMBYTES
    #define GFP_NUMBYTES ((GFP_NUMBITS+7)/8)
    #endif
    
    #ifndef GFP_CONTAINER
    #define GFP_CONTAINER
    #if GFP_NUMBYTES == 1
    typedef unsigned char gfp_element;
    #elif GFP_NUMBYTES <= 4
    typedef unsigned int gfp_element;
    #else
    typedef unsigned char[NUMBYTES] gfp_element;
    #endif
    #endif
    
#else
    
    #include "bi.h"
    bi prime_modulus;
    #define GF_PRIME_MODULUS prime_modulus
    #define GFP_NUMBITS bi_bitsize(prime_modulus)
    #define GFP_NUMBYTES ((GFP_NUMBITS+7)/8)
    typedef bi gfp_element;

#endif

gfp_element gfp( int castee );
gfp_element gfp_init( unsigned int size );
gfp_element gfp_clone( gfp_element elm );
int gfp_destroy( gfp_element elm );
int gfp_copy( gfp_element* dest, gfp_element source );
int gfp_zero( gfp_element* elm );
int gfp_one( gfp_element* elm );
int gfp_random( gfp_element* elm, unsigned char * randomness );
int gfp_random_invertible( gfp_element* elm, unsigned char * randomness );
int gfp_compare( gfp_element lhs, gfp_element rhs );
int gfp_copy( gfp_element * dest, gfp_element source );
int gfp_add( gfp_element * res, gfp_element lhs, gfp_element rhs );
int gfp_subtract( gfp_element * res, gfp_element lhs, gfp_element rhs );
int gfp_negate( gfp_element * res, gfp_element elm );
int gfp_multiply( gfp_element * res, gfp_element lhs, gfp_element rhs );
int gfp_divide( gfp_element * quo, gfp_element numerator, gfp_element divisor );
int gfp_inverse( gfp_element * res, gfp_element elm );
int gfp_print( gfp_element elm );
int gfp_is_one( gfp_element elm );
int gfp_is_zero( gfp_element elm );

#endif

