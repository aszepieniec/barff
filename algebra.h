#ifndef ALGEBRA_H
#define ALGEBRA_H

#include "csprng.h"

#ifndef MOD
#define MOD 11
#endif


/**
 * finite field arithmetic
 */
int xgcd( int x, int y, int* a, int* b, int* c );
unsigned char gf_inverse( unsigned char element );

/**
 * matrix 
 */

typedef struct 
{
    unsigned short int width;
    unsigned short int height;
    unsigned char * data;
} gfmatrix;

gfmatrix gfm( unsigned short int height, unsigned short int width, unsigned char * pdata );

gfmatrix gfm_init( unsigned short int height, unsigned short int width );
int gfm_destroy( gfmatrix fm );
int gfm_copy( gfmatrix dest, gfmatrix source );
gfmatrix gfm_clone( gfmatrix source );

int gfm_eye( gfmatrix mat );
int gfm_is_eye( gfmatrix mat );
int gfm_equals( gfmatrix lhs, gfmatrix rhs );
int gfm_zeros( gfmatrix mat );
int gfm_random( gfmatrix mat, csprng * rng );
int gfm_random_upper_triangular( gfmatrix mat, csprng * rng  );
int gfm_random_invertible( gfmatrix mat, csprng * rng  );
int gfm_transpose( gfmatrix * mat );
int gfm_multiply( gfmatrix dest, gfmatrix left, gfmatrix right );
int gfm_multiply_transpose( gfmatrix dest, gfmatrix left, gfmatrix rightT );
int gfm_transpose_multiply( gfmatrix dest, gfmatrix leftT, gfmatrix right );
int gfm_multiply_constant( gfmatrix dest, gfmatrix source, unsigned char constant );
int gfm_sum( gfmatrix dest, gfmatrix left_matrix, gfmatrix right_matrix );
int gfm_weighted_sum( gfmatrix dest, unsigned char left_constant, gfmatrix left_matrix, unsigned char right_constant, gfmatrix right_matrix );
int gfm_rowop( gfmatrix mat, unsigned short int destrow, unsigned short int sourcerow, unsigned char constant, unsigned short int offset );
int gfm_fliprows( gfmatrix mat, unsigned short int destrow, unsigned short int sourcerow );
int gfm_scalerow( gfmatrix mat, unsigned short int rowidx, unsigned char constant );
int gfm_redech( gfmatrix mat );
int gfm_solve( gfmatrix coeffs, gfmatrix target, gfmatrix solution, gfmatrix * kernel );

int gfm_stack( gfmatrix dest, gfmatrix top, gfmatrix bottom );
int gfm_cat( gfmatrix dest, gfmatrix left, gfmatrix right );
int gfm_slice( gfmatrix dest, gfmatrix mat, unsigned short int row_start, unsigned short int col_start );

int gfm_inverse( gfmatrix dest, gfmatrix mat );

int gfm_print( gfmatrix mat );

/**
 * univariate polynomials
 */

typedef struct
{
    unsigned short int degree;
    unsigned char * data;
} gfpolynomial;

/**
 * (homogeneous) quadratic systems
 */
typedef struct
{
    gfmatrix * quadratic_forms;
    unsigned short int n;
    unsigned short int m;
} hqsystem;

hqsystem hqs( gfmatrix* qfs, unsigned short int n, unsigned short int m );
hqsystem hqs_init( unsigned short int n, unsigned short int m );
int hqs_destroy( hqsystem sys );
int hqs_random( hqsystem sys, csprng * rng );
int hqs_copy( hqsystem dest, hqsystem source );
hqsystem hqs_copy_new( hqsystem source );
int hqs_compose_output( gfmatrix T, hqsystem P );
int hqs_compose_input( hqsystem P, gfmatrix S );
int hqs_eval( gfmatrix y, hqsystem sys, gfmatrix x );

#endif

