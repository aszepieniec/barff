#include "hqs.h"

#include <stdlib.h>
#include <stdio.h>

/**
 * hqs
 * Create homogeneous quadratic system object from a list of
 * of quadratic forms and the system's dimensions.
 */
hqsystem hqs( gfpmatrix* qfs, unsigned  int n, unsigned  int m )
{
    hqsystem hqs;
    hqs.quadratic_forms = qfs;
    hqs.n = n;
    hqs.m = m;
    return hqs;
}

/**
 * hqs_init
 * Creates a homogeneous quadratic system object from the given
 * dimensions and allocates memory for it. Don't forget to call
 * hqs_destroy afterwards.
 */
hqsystem hqs_init( unsigned  int n, unsigned  int m )
{
    unsigned int i;
    hqsystem hqs;
    hqs.n = n;
    hqs.m = m;
    hqs.quadratic_forms = malloc(sizeof(gfpmatrix) * m);
    for( i = 0 ; i < m ; ++i )
    {
        hqs.quadratic_forms[i] = gfpm_init(n,n);
    }
    return hqs;
}

/**
 * hqs_destroy
 * Deallocates space allocated for a homogeneous quadratic system
 * object.
 */
int hqs_destroy( hqsystem sys )
{
    unsigned int i;
    for( i = 0 ; i < sys.m ; ++i )
    {
        gfpm_destroy(sys.quadratic_forms[i]);
    }
    free(sys.quadratic_forms);
    return 1;
}

/**
 * hqs_copy
 * Copies the data associated to a homogeneous quadratic system. Does
 * not allocate the space necessary -- you have to do that yourself.
 * (And consequently you have the option of keeping everything in the
 * stack.)
 */
int hqs_copy( hqsystem dest, hqsystem source )
{
    unsigned int i;
#ifdef DEBUG
    if( dest.n != source.n || dest.m != source.m )
    {
        printf("in hqs_copy: cannot copy hqs object because dimensions do not match! dest: %i -> %i vs. source: %i -> %i\n", dest.n, dest.m, source.n, source.m);
        return 0;
    }
#endif

    for( i = 0 ; i < dest.m ; ++i )
    {
        gfpm_copy(dest.quadratic_forms[i], source.quadratic_forms[i]);
    }
    return 1;
}

/**
 * hqs_clone
 * Copy one homogeneous quadratic system to a new one. Remember to
 * destroy it when scope ends!
 * @params
 *  * source : the homogeneous quadratic system to copy
 * @return
 *  * dest : a new homogeneous quadratic system identical to source
 */
hqsystem hqs_clone( hqsystem source )
{
    hqsystem dest;
    dest = hqs_init(source.n, source.m);
    hqs_copy(dest, source);
    return dest;
}

/**
 * hqs_random
 * Given an empty homogeneous quadratic system, assign random values
 * to its coefficients.
 * The amount of randomness required is m*n*n bytes.
 * @params
 *  * sys : a homogeneous quadratic system; the dimensions of this
 *    system will be retained; the coefficients (entries of the
 *    matrices) will be forgotten
 *  * randomness : a pointer to a large enough string of random bytes
 *    "large enough" means n*n*m*(1+GF_NUMBYTES)
 * @result
 *  * sys will contain random coefficients
 * @return
 *  * 1 if success, 0 otherwise
 */
int hqs_random( hqsystem sys, unsigned char * randomness )
{
    unsigned int i, j, k, l;

    l = 0;
    for( k = 0 ; k < sys.m ; ++k )
    {
        for( i = 0 ; i < sys.n ; ++i )
        {
            for( j = 0 ; j < sys.n ; ++j )
            {
                gfp_random(&sys.quadratic_forms[k].data[i*sys.n + j], &randomness[l]);
                l = l + (GFP_NUMBITS + sizeof(unsigned long int) * 8 - 1) / (sizeof(unsigned long int) * 8);
            }
        }
    }

    return 1;
}

/**
 * hqs_compose_output
 * Compose a homogeneous quadratic system on the left with a
 * linear transform.
 * @params
 *  * T : a gfpmatrix object of dimensions mxm
 *  * F : an hqsystem object of dimensions n -> m
 * @promise
 *  * T is square
 *  * has the same number of columns as the number of quadratic
 *    forms in F
 * @result
 *  new.F = T * old.F
 * @return
 * 1 if success, 0 otherwise
 */
int hqs_compose_output( gfpmatrix T, hqsystem F )
{
    unsigned int i, j;
    /* declare helper variables */
    hqsystem P;
    gfpmatrix temp;

#ifdef DEBUG
    if( T.width != F.m )
    {
        printf("in hqs_compose_output: cannot compose homogeneous quadratic system with linear transform on output side because dimension mismatch! T: %ix%i vs F: %i -> %i\n", T.height, T.width, F.n, F.m);
        return 0;
    }
#endif

    /* init helper variables */
    P = hqs_init(F.n, T.height);
    temp = gfpm_init(F.n, F.n);

    /* perform multiplication */
    for( i = 0 ; i < T.height ; ++i )
    {
        gfpm_zeros(P.quadratic_forms[i]);
        for( j = 0 ; j < T.width ; ++j )
        {
            gfpm_multiply_constant(temp, F.quadratic_forms[j], T.data[i*T.width + j]);
            gfpm_add(P.quadratic_forms[i], P.quadratic_forms[i], temp);
        }
    }

    /* copy to argument */
    hqs_copy(F, P);

    /* destroy helper variables */
    hqs_destroy(P);
    gfpm_destroy(temp);

    return 1;
}

/**
 * hqs_compose_input
 * Compose a homogeneous quadratic system with a linear transform on
 * the input variables.
 * @params
 *  * F : homogeneous quadratic system object to which the input
 *    transform is to be applied
 *  * S : matrix object that represents the linear transform
 * @promise
 *  * F.n = S.width
 *  * S.width = S.height
 * @result
 *  * new.F = old.F o S
 * @return
 *  * 1 if success, 0 otherwise
 */
int hqs_compose_input( hqsystem F, gfpmatrix S )
{
    unsigned int i;

    /* declare helper variables */
    gfpmatrix temp;

    /* debug stuff */
#ifdef DEBUG
    if( F.n != S.height || S.height != S.width )
    {
        printf("hqs_compose_input: cannot compose homogeneous quadratic system with linear transform on input side because of dimension mismatch! F : %i -> %i vs. S : %ix%i\n", F.n, F.m, S.height, S.width);
        return 0;
    }
#endif

    /* init helper variables */

    temp = gfpm_init(S.height, S.height);

    /* perform multiplication */
    for( i = 0 ; i < F.m ; ++i )
    {
        gfpm_transpose_multiply(temp, S, F.quadratic_forms[i]);
        gfpm_multiply(F.quadratic_forms[i], temp, S);
    }

    /* destroy helper variables */
    gfpm_destroy(temp);

    return 1;
}

/**
 * hqs_eval
 * Evaluate a homogeneous quadratic system in a vector or, by
 * treating the columns as a list of vectors, as a matrix.
 */
int hqs_eval( gfpmatrix y, hqsystem sys, gfpmatrix x )
{
    unsigned int i, j, k;
    gfpmatrix vector, transposed_vector, temp, e;
    gfp_element edata;

    edata = gfp_init(1);

#ifdef DEBUG
    if( y.height != sys.m || sys.n != x.height || y.width != x.width )
    {
        printf("hqs_eval: cannot evaluate quadratic system because of dimension mismatch! F: %i -> %i vs. in: %ix%i, out: %ix%i\n", sys.n, sys.m, x.height, x.width, y.height, y.width);
        return 0;
    }
#endif

    vector = gfpm_init(sys.n, 1);
    transposed_vector = gfpm(1, sys.n, vector.data);
    temp = gfpm_init(sys.n, 1);
    e = gfpm(1, 1, &edata);

    for( j = 0 ; j < x.width ; ++j )
    {
        gfpm_slice(vector, x, 0, j);
        for( i = 0 ; i < x.height ; ++i )
        {
            for( k = 0 ; k < sys.m ; ++k )
            {
                gfpm_multiply(temp, sys.quadratic_forms[k], vector);
                gfpm_multiply(e, transposed_vector, temp);
                gfp_copy(&y.data[k*y.width + j], e.data[0]);
            }
        }
    }

    gfpm_destroy(vector);
    gfpm_destroy(temp);
    gfp_destroy(edata);

    return 1;
}

