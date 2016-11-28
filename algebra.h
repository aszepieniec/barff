#ifndef ALGEBRA_H
#define ALGEBRA_H

#include "csprng.h"

#ifndef MOD
#define MOD 199
#endif

typedef unsigned short int field_element;

/**
 * finite field arithmetic
 */
int xgcd( int x, int y, int* a, int* b, int* c );
field_element gf_inverse( field_element element );

/**
 * matrix 
 */

typedef struct 
{
    unsigned  int width;
    unsigned  int height;
    field_element * data;
} gfmatrix;

gfmatrix gfm( unsigned  int height, unsigned  int width, field_element * pdata );

gfmatrix gfm_init( unsigned  int height, unsigned  int width );
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
int gfm_multiply_constant( gfmatrix dest, gfmatrix source, field_element constant );
int gfm_sum( gfmatrix dest, gfmatrix left_matrix, gfmatrix right_matrix );
int gfm_weighted_sum( gfmatrix dest, field_element left_constant, gfmatrix left_matrix, field_element right_constant, gfmatrix right_matrix );
int gfm_rowop( gfmatrix mat, unsigned  int destrow, unsigned  int sourcerow, field_element constant, unsigned  int offset );
int gfm_fliprows( gfmatrix mat, unsigned  int destrow, unsigned  int sourcerow );
int gfm_scalerow( gfmatrix mat, unsigned  int rowidx, field_element constant );
int gfm_redech( gfmatrix mat );
int gfm_solve( gfmatrix coeffs, gfmatrix target, gfmatrix solution, gfmatrix * kernel );

int gfm_stack( gfmatrix dest, gfmatrix top, gfmatrix bottom );
int gfm_cat( gfmatrix dest, gfmatrix left, gfmatrix right );
int gfm_slice( gfmatrix dest, gfmatrix mat, unsigned  int row_start, unsigned  int col_start );

int gfm_inverse( gfmatrix dest, gfmatrix mat );

int gfm_print( gfmatrix mat );

/**
 * univariate polynomials
 */

typedef struct
{
    unsigned  int degree;
    field_element * data;
} gfpolynomial;

/**
 * (homogeneous) quadratic systems
 */
typedef struct
{
    gfmatrix * quadratic_forms;
    unsigned  int n;
    unsigned  int m;
} hqsystem;

hqsystem hqs( gfmatrix* qfs, unsigned  int n, unsigned  int m );
hqsystem hqs_init( unsigned  int n, unsigned  int m );
int hqs_destroy( hqsystem sys );
int hqs_random( hqsystem sys, csprng * rng );
int hqs_copy( hqsystem dest, hqsystem source );
hqsystem hqs_copy_new( hqsystem source );
int hqs_compose_output( gfmatrix T, hqsystem P );
int hqs_compose_input( hqsystem P, gfmatrix S );
int hqs_eval( gfmatrix y, hqsystem sys, gfmatrix x );

#endif

