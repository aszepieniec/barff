#ifndef ALGEBRA_H
#define ALGEBRA_H

#ifndef MOD
#define MOD 61
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

int gfm_eye( gfmatrix mat );
int gfm_zeros( gfmatrix mat );
int gfm_random( gfmatrix mat, unsigned char * randomness );
int gfm_random_upper_triangular( gfmatrix mat, unsigned char * randomness );
int gfm_transpose( gfmatrix * mat );
int gfm_multiply( gfmatrix dest, gfmatrix left, gfmatrix right );
int gfm_multiply_constant( gfmatrix dest, unsigned char constant );
int gfm_sum( gfmatrix dest, gfmatrix left_matrix, gfmatrix right_matrix );
int gfm_weighted_sum( gfmatrix dest, unsigned char left_constant, gfmatrix left_matrix, unsigned char right_constant, gfmatrix right_matrix );
int gfm_rowop( gfmatrix mat, unsigned short int destrow, unsigned short int sourcerow, unsigned char constant, unsigned short int offset );
int gfm_fliprows( gfmatrix mat, unsigned short int destrow, unsigned short int sourcerow );
int gfm_scalerow( gfmatrix mat, unsigned short int rowidx, unsigned char constant );
int gfm_redech( gfmatrix mat );

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
int hqs_copy( hqsystem dest, hqsystem source );
int hqs_compose_output( gfmatrix T, hqsystem P );
int hqs_compose_input( hqsystem P, gfmatrix S );

#endif

