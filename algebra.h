#ifndef ALGEBRA_H
#define ALGEBRA_H

#include "gfp.h"

/**
 * finite field arithmetic
 */
int xgcd( int x, int y, int* a, int* b, int* c );
gfp_element gf_inverse( gfp_element element );

/**
 * matrix 
 */

typedef struct 
{
    unsigned  int width;
    unsigned  int height;
    gfp_element * data;
} gfpmatrix;

gfpmatrix gfpm( unsigned  int height, unsigned  int width, gfp_element * pdata );

gfpmatrix gfpm_init( unsigned  int height, unsigned  int width );
int gfpm_destroy( gfpmatrix fm );
int gfpm_copy( gfpmatrix dest, gfpmatrix source );
gfpmatrix gfpm_clone( gfpmatrix source );

int gfpm_eye( gfpmatrix mat );
int gfpm_is_eye( gfpmatrix mat );
int gfpm_equals( gfpmatrix lhs, gfpmatrix rhs );
int gfpm_zeros( gfpmatrix mat );
int gfpm_random( gfpmatrix mat, unsigned char * randomness );
int gfpm_random_upper_triangular( gfpmatrix mat, unsigned char * randomness  );
int gfpm_random_invertible( gfpmatrix mat, unsigned char * randomness  );
int gfpm_transpose( gfpmatrix * mat );
int gfpm_multiply( gfpmatrix dest, gfpmatrix left, gfpmatrix right );
int gfpm_multiply_transpose( gfpmatrix dest, gfpmatrix left, gfpmatrix rightT );
int gfpm_transpose_multiply( gfpmatrix dest, gfpmatrix leftT, gfpmatrix right );
int gfpm_multiply_constant( gfpmatrix dest, gfpmatrix source, gfp_element constant );
int gfpm_sum( gfpmatrix dest, gfpmatrix left_matrix, gfpmatrix right_matrix );
int gfpm_weighted_sum( gfpmatrix dest, gfp_element left_constant, gfpmatrix left_matrix, gfp_element right_constant, gfpmatrix right_matrix );
int gfpm_rowop( gfpmatrix mat, unsigned  int destrow, unsigned  int sourcerow, gfp_element constant, unsigned  int offset );
int gfpm_fliprows( gfpmatrix mat, unsigned  int destrow, unsigned  int sourcerow );
int gfpm_scalerow( gfpmatrix mat, unsigned  int rowidx, gfp_element constant );
int gfpm_redech( gfpmatrix mat );
int gfpm_solve( gfpmatrix coeffs, gfpmatrix target, gfpmatrix solution, gfpmatrix * kernel );
int gfpm_inspan( gfpmatrix vec, gfpmatrix mat );

int gfpm_stack( gfpmatrix dest, gfpmatrix top, gfpmatrix bottom );
int gfpm_cat( gfpmatrix dest, gfpmatrix left, gfpmatrix right );
int gfpm_slice( gfpmatrix dest, gfpmatrix mat, unsigned  int row_start, unsigned  int col_start );

int gfpm_inverse( gfpmatrix dest, gfpmatrix mat );

int gfpm_print( gfpmatrix mat );

/**
 * univariate polynomials
 */

typedef struct
{
    unsigned  int degree;
    gfp_element * data;
} gfpolynomial;

/**
 * (homogeneous) quadratic systems
 */
typedef struct
{
    gfpmatrix * quadratic_forms;
    unsigned  int n;
    unsigned  int m;
} hqsystem;

hqsystem hqs( gfpmatrix* qfs, unsigned  int n, unsigned  int m );
hqsystem hqs_init( unsigned  int n, unsigned  int m );
int hqs_destroy( hqsystem sys );
int hqs_random( hqsystem sys, unsigned char * randomness );
int hqs_copy( hqsystem dest, hqsystem source );
hqsystem hqs_clone( hqsystem source );
int hqs_compose_output( gfpmatrix T, hqsystem P );
int hqs_compose_input( hqsystem P, gfpmatrix S );
int hqs_eval( gfpmatrix y, hqsystem sys, gfpmatrix x );

#endif

