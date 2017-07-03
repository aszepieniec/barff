#ifndef HQS_H
#define HQS_H

/**
 * HQS.H
 * Routines and structures for homogeneous quadratic systems.
 */

#include "gfp.h"
#include "gfpm.h"

/**
 *  * univariate polynomials
 *   */

typedef struct
{
    unsigned  int degree;
    gfp_element * data;
} gfpolynomial;

/**
 *  * (homogeneous) quadratic systems
 *   */
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

