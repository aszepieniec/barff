#include "algebra.h"

#include <stdlib.h>
#include <stdio.h>

/**
 * xgcd
 * Computes the Bezout relation
 *   a * x  +  b * y  =  gcd
 * from a and b using standard ints.
 */
int xgcd( int a, int b, int* x, int* y, int* gcd )
{
    int q, r;
    int u;
    int v;
    int m, n;

    *x = 0;
    *y = 1;
    u = 1;
    v = 0;
    while( a != 0 )
    {
        q = b / a;
        r = b % a;
        m = *x - u*q;
        n = *y - v*q;
        b = a;
        a = r;
        *x = u;
        *y = v;
        u = m;
        v = n;
    }
    *gcd = b;
    
    return *gcd;
}

/**
 * gf_inverse
 * Computes the multiplicative inverse of a field element.
 */
gfp_element gf_inverse( gfp_element element )
{
    int a;
    int b;
    int x, y, g;

    a = element;
    b = GF_PRIME_MODULUS;
    xgcd(a, b, &x, &y, &g);
    x = (GF_PRIME_MODULUS + (x % GF_PRIME_MODULUS)) % GF_PRIME_MODULUS;
    return x;
}

/**
 * gfpm
 * Create gfpmatrix object with given buffer.
 */
gfpmatrix gfpm( unsigned  int height, unsigned  int width, gfp_element* pdata )
{
    gfpmatrix mat;
    mat.height = height;
    mat.width = width;
    mat.data = pdata;
    return mat;
}

/**
 * gfpm_init
 * Create a field matrix object and allocate space for it.
 */
gfpmatrix gfpm_init( unsigned  int height, unsigned  int width )
{
    gfpmatrix mat;
    mat.width = width;
    mat.height = height;
    mat.data = malloc(width*height*sizeof(gfp_element));
    /*
    printf("created gfpm object with data member set to memory address %#010x\n", mat.data);
    */
    return mat;
}

/**
 * gfpm_destroy
 * Deallocates space allocated to a field matrix object. Call this
 * function before closing the scope where the field matrix object
 * was initialized.
 * @return
 *  * 1 if success
 */
int gfpm_destroy( gfpmatrix fm )
{
    /*
    printf("destroying gfpm object with data member set to memory address %#010x\n", fm.data);
    */
    free(fm.data);
    fm.width = 0;
    fm.height = 0;
    return 1;
}

/**
 * gfpm_copy
 * Copy the contents of one matrix to another. Does not allocate
 * memory for the new object; you must do that yourself! (Or use
 * gfpm_clone instead.)
 * @promise
 *  * dest.width >= source.width
 *  * dest.height >= source.height
 * @return
 *  * 1 if success
 */
int gfpm_copy( gfpmatrix dest, gfpmatrix source )
{
    unsigned int i, j;
    for( i = 0 ; i < source.height ; ++i )
    {
        for( j = 0 ; j < source.width ; ++j )
        {
            dest.data[i*dest.width + j] = source.data[i*source.width + j];
        }
    }

    return 1;
}

/**
 * gfpm_clone
 * Copy one matrix into a new object. Don't forget to call
 * gfpm_destroy at the end of scope!
 */
gfpmatrix gfpm_clone( gfpmatrix source )
{
    gfpmatrix mat;
    mat = gfpm_init(source.height, source.width);
    gfpm_copy(mat, source);
    return mat;
}

/**
 * gfpm_eye
 * Set a matrix to the identity matrix.
 * @param
 *  * eye : matrix object to set to identity
 * @returns
 *  * 1 if success, 0 otherwise
 */
int gfpm_eye( gfpmatrix eye )
{
    unsigned  int i, j;
    for( i = 0 ; i < eye.height ; ++i )
    {
        for( j = 0 ; j < eye.width ; ++j )
        {
            eye.data[i*eye.width + j] = 0;
        }
        eye.data[i*eye.width + i] = 1;
    }
    return 1;
}

/**
 * Decide whether the given matrix is an identity matrix or, for
 * rectangular matrices, whether the main diagonal has ones and all
 * the rest is zero.
 * @return
 *  * 1 if identity, 0 otherwise
 */
int gfpm_is_eye( gfpmatrix eye )
{
    unsigned int i, j;
    int b = 1;
    for( i = 0 ; i < eye.height ; ++i )
    {
        for( j = 0 ; j < eye.width ; ++j )
        {
            if( i == j )
            {
                if( eye.data[i*eye.width + j] != 1 )
                {
                    b = 0;
                }
            }
            else
            {
                if( eye.data[i*eye.width + j] != 0 )
                {
                    b = 0;
                }
            }
        }
    }

    return b;
}

/**
 * gfpm_equals
 * Test two matrices for equality.
 * @return
 *  * 1 if equal, 0 otherwise
 */
int gfpm_equals( gfpmatrix lhs, gfpmatrix rhs )
{
    unsigned int i, j;
    int b;

    if( lhs.width != rhs.width || lhs.height != rhs.height )
    {
        return 0;
    }

    b = 1;
    for( i = 0 ; i < lhs.height ; ++i )
    {
        for( j = 0 ; j < lhs.width ; ++j )
        {
            b = b & (lhs.data[lhs.width*i + j] == rhs.data[rhs.width*i + j]);
        }
    }

    return b;
}

/**
 * gfpm_zero
 * Sets a matrix to all zeros.
 * @param
 *  * zero : matrix object to set to zero
 * @returns
 *  * 1 if success
 */
int gfpm_zeros( gfpmatrix zero )
{
    unsigned  int i, j;
    for( i = 0 ; i < zero.height ; ++i )
    {
        for( j = 0 ; j < zero.width ; ++j )
        {
            zero.data[i*zero.width + j] = 0;
        }
    }
    return 1;
}

/**
 * gfpm_random
 * Put random values into the matrix.
 * @params
 *  * random : matrix objects with data already allocated and whose
 *    elements are to be assigned random values
 *  * randomness : pointer to large-enough string of random bytes
 *    "large enough" means n*n*sizeof(unsigned int)
 * @result
 *  * random <-$- matrix_space(random.height, random.width)
 */
int gfpm_random( gfpmatrix random, unsigned char * randomness )
{
    unsigned  int i, j;
    unsigned int l;
    unsigned int * rand = (unsigned int *)randomness;
   
    l = 0;
    for( i = 0 ; i < random.height ; ++i )
    {
        for( j = 0 ; j < random.width ; ++j )
        {
            random.data[i*random.width + j] = rand[l++] % GF_PRIME_MODULUS;
        }
    }

    return 1;
}

/**
 * gfpm_random_upper_triangular
 * Set the matrix to a random upper triangular with ones on the
 * diagonal.
 * @params
 *  * random : matrix objects with data already allocated and whose
 *    elements above the diagonal are to be assigned random values;
 *    whereas the elements above the diagonal are to be 0 and the
 *    elements on the diagonal are to be 1.
 *  * randomness : pointer to large-enough string of random bytes
 *    "large enough" means n*(n-1)/2*sizeof(unsigned int)
 * @result
 *  * random <-$- matrix_space(random.height, random.width)
 *    subject to
 *    forall i . random[i,i] = 1
 *    forall i, j>i . random[i,j] = 0
 */
int gfpm_random_upper_triangular( gfpmatrix random, unsigned char * randomness )
{
    unsigned  int i, j;
    unsigned int l;
    unsigned  int * rand = (unsigned int *)randomness;

    l = 0;
    for( i = 0 ; i < random.height ; ++i )
    {
        for( j = 0 ; j < i ; ++j )
        {
            random.data[i*random.width + j] = rand[l++] % GF_PRIME_MODULUS;
        }
        random.data[i*random.width + i] = 1;
        for( j = i+1 ; j < random.width ; ++j )
        {
            random.data[i*random.width + j] = 0;
        }
    }

    return 1;
}

/**
 * gfpm_random_invertible
 * Generate a random invertible matrix. This is accomplished by first
 * generating two random triangular matrices, where one has ones on
 * the diagonal is upper-triangular and the other has random nonzero
 * elements on the diagonal in addition to being lower-triangular.
 * The generated matrix is the product of the two triangular ones.
 * @params
 *  * mat : the matrix object in which to store the generated random
 *    invertible matrix
 *  * randomness : pointer to large enough string of random bytes
 *    "large enough" means n*n*sizeof(unsigned int)
 * @pre
 *  * mat.height = mat.width
 * @result
 *  * mat <-$- GL(mat.height)
 * @return
 *  * 1 if success, 0 otherwise
 */
int gfpm_random_invertible( gfpmatrix mat, unsigned char * randomness )
{
    gfpmatrix utm, ltm;
    unsigned int offset;
    unsigned int i;

    unsigned  int * rand = (unsigned int *) randomness;

#ifdef DEBUG
    if( mat.height != mat.width )
    {
        printf("gfpm_random_invertible: cannot generate random invertible matrix because matrix not square! %ix%i\n", mat.height, mat.width);
        return 0;
    }
#endif

    utm = gfpm_init(mat.height, mat.width);
    ltm = gfpm_init(mat.height, mat.width);

    /* generate triangular matrices with ones on the diagonal */
    gfpm_random_upper_triangular(utm, randomness);
    gfpm_random_upper_triangular(ltm, randomness + sizeof(unsigned int)*(mat.height * (mat.width-1))/2);
    gfpm_transpose(&ltm);

    /* set the diagonal elements of one matrix to random nonzero
     * elements */
    for( i = 0 ; i < utm.height ; ++i )
    {
        utm.data[i*utm.width + i] = 1 + (rand[i + mat.height * (mat.width-1)] % (GF_PRIME_MODULUS - 1));
    }

    /* multiply L * U to get the random invertible matrix */
    gfpm_multiply(mat, ltm, utm);

    gfpm_destroy(ltm);
    gfpm_destroy(utm);

    return 1;
}

/**
 * gfpm_transpose_square
 * Perform a matrix transposition in situ.
 */
int gfpm_transpose( gfpmatrix * trans )
{
    gfp_element a;
    unsigned  int i, j;

    gfpmatrix T;

    T = gfpm_init(trans->height, trans->width);
    gfpm_copy(T, *trans);

    a = trans->width;
    trans->width = trans->height;
    trans->height = a;

    for( i = 0 ; i < trans->height ; ++i )
    {
        for( j = 0 ; j < trans->width ; ++j )
        {
            trans->data[i*trans->width + j] = T.data[j*T.width + i];
        }
    }

    gfpm_destroy(T);

    return 1;
}

/**
 * gfpm_multiply
 * Multiplies two matrices and stores the result in the third, which
 * should be pre-allocated.
 * @params
 *  * dest : field matrix object to store the matrix product
 *  * left, right : field matrix object representing left- and right-
 *    hand-sides respectively.
 * @return
 *  * 1 if success, 0 otherwise
 * NOTE. Modular reduction is applied only once, after computing the
 * inner product between left row and right column and storing the
 * intermediate result in an int. If the int is four bytes and the
 * field element takes up a full byte, then in the worst case every
 * product requires two bytes and every 2^k additions requires k
 * extra bits. Since we have 16 bits to spare we have at most 2^16
 * additions which is anyway the max. height and width of matrices
 * that can be stored in a  unsigned int.
 */
int gfpm_multiply( gfpmatrix dest, gfpmatrix left, gfpmatrix right )
{
    unsigned  int i, j, k;
    unsigned int acc;

    #ifdef DEBUG
        if( dest.height != left.height || dest.width != right.width || left.width != right.height )
        {
            printf("in gfpm_multiply: trying to multiply matrices with unmatched dimensions: %ix%i * %ix%i = %ix%i\n", left.height, left.width, right.height, right.width, dest.height, dest.width);
            return 0;
        }
    #endif

    for( i = 0 ; i < left.height ; ++i )
    {
        for( j = 0 ; j < right.width ; ++j )
        {
            acc = 0;
            for( k = 0 ; k < left.width ; ++k )
            {
                acc = acc + left.data[i*left.width + k] * right.data[k*right.width + j];
            }
            dest.data[i*dest.width + j] = acc % GF_PRIME_MODULUS;
        }
    }

    return 1;
}

/**
 * gfpm_multiply_transpose
 * Multiplies the left hand side matrix with the transpose of the
 * hand side matrix and stores the result in the third, which
 * should be pre-allocated.
 * @params
 *  * dest : field matrix object to store the matrix product
 *  * left, rightT : field matrix object representing left- and right-
 *    hand-sides respectively.
 * @return
 *  * 1 if success, 0 otherwise
 * NOTE. Modular reduction is applied only once, after computing the
 * inner product between left row and right column and storing the
 * intermediate result in an int. If the int is four bytes and the
 * field element takes up a full byte, then in the worst case every
 * product requires two bytes and every 2^k additions requires k
 * extra bits. Since we have 16 bits to spare we have at most 2^16
 * additions which is anyway the max. height and width of matrices
 * that can be stored in a  unsigned int.
 */
int gfpm_multiply_transpose( gfpmatrix dest, gfpmatrix left, gfpmatrix rightT )
{
    unsigned  int i, j, k;
    unsigned int acc;

    #ifdef DEBUG
        if( dest.height != left.height || dest.width != rightT.height || left.width != rightT.width )
        {
            printf("in gfpm_multiply_transpose: trying to multiply matrices with unmatched dimensions: %ix%i * (%ix%i)^T = %ix%i\n", left.height, left.width, rightT.height, rightT.width, dest.height, dest.width);
            return 0;
        }
    #endif

    for( i = 0 ; i < left.height ; ++i )
    {
        for( j = 0 ; j < rightT.height ; ++j )
        {
            acc = 0;
            for( k = 0 ; k < left.width ; ++k )
            {
                acc = acc + left.data[i*left.width + k] * rightT.data[j*rightT.width + k];
            }
            dest.data[i*dest.width + j] = acc % GF_PRIME_MODULUS;
        }
    }

    return 1;
}

/**
 * gfpm_transpose_multiply
 * Multiplies the transpose of the left hand side matrix with the
 * hand side matrix and stores the result in the third, which
 * should be pre-allocated.
 * @params
 *  * dest : field matrix object to store the matrix product
 *  * leftT, right : field matrix object representing left- and right-
 *    hand-sides respectively.
 * @return
 *  * 1 if success, 0 otherwise
 * NOTE. Modular reduction is applied only once, after computing the
 * inner product between left row and right column and storing the
 * intermediate result in an int. If the int is four bytes and the
 * field element takes up a full byte, then in the worst case every
 * product requires two bytes and every 2^k additions requires k
 * extra bits. Since we have 16 bits to spare we have at most 2^16
 * additions which is anyway the max. height and width of matrices
 * that can be stored in a  unsigned int.
 */
int gfpm_transpose_multiply( gfpmatrix dest, gfpmatrix leftT, gfpmatrix right )
{
    unsigned  int i, j, k;
    unsigned int acc;

    #ifdef DEBUG
        if( dest.height != leftT.width || dest.width != right.width || leftT.height != right.height )
        {
            printf("in gfpm_transpose_multiply: trying to multiply matrices with unmatched dimensions: (%ix%i)^T * %ix%i = %ix%i\n", leftT.height, leftT.width, right.height, right.width, dest.height, dest.width);
            return 0;
        }
    #endif

    for( i = 0 ; i < leftT.width ; ++i )
    {
        for( j = 0 ; j < right.width ; ++j )
        {
            acc = 0;
            for( k = 0 ; k < leftT.height ; ++k )
            {
                acc = acc + leftT.data[k*leftT.width + i] * right.data[k*right.width + j];
            }
            dest.data[i*dest.width + j] = acc % GF_PRIME_MODULUS;
        }
    }

    return 1;
}

/**
 * gfpm_multiply_constant
 * Multiply the matrix with a constant.
 * @return
 *  * 1 if success
 */
int gfpm_multiply_constant( gfpmatrix dest, gfpmatrix source, gfp_element constant )
{
    unsigned  int i, j;

#ifdef DEBUG
    if( dest.width != source.width || dest.height != source.height )
    {
        printf("gfpm_multiply_constant: cannot multiply matrix with constant because dimensions of destination do not match those of source! %ix%i <- %ix%i\n", dest.height, dest.width, source.height, source.width);
        return 0;
    }
#endif

    for( i = 0 ; i < dest.height ; ++i )
    {
        for( j = 0 ; j < dest.width ; ++j )
        {
            dest.data[i*dest.width + j] = (source.data[i*dest.width + j] * constant) % GF_PRIME_MODULUS;
        }
    }
    return 1;
}

/**
 * gfpm_add
 * Add one matrix to another and store the result in a third. The
 * third matrix must be preallocated.
 * @params
 *  * dest : the matrix object to store the result in
 *  * left, right : the matrix objects to add together; they must 
 *    have the same dimensions!
 * @return
 *  * 1 if success, 0 otherwise
 */
int gfpm_add( gfpmatrix dest, gfpmatrix left, gfpmatrix right )
{
    unsigned  int i, j;

    #ifdef DEBUG
        if( dest.width != left.width || left.width != right.width || dest.height != left.height || left.height != right.height )
        {
            printf("in gfpm_add: trying to add matrices of incompatible dimensions! %ix%i + %ix%i = %ix%i\n", left.height, left.width, right.height, right.width, dest.height, dest.width);
            return 0;
        }
    #endif
    
    for( i = 0 ; i < dest.height ; ++i )
    {
        for( j = 0 ; j < dest.width ; ++j )
        {
            dest.data[i*dest.width + j] = (left.data[i*left.width + j] + right.data[i*right.width + j]) % GF_PRIME_MODULUS;
        }
    }

    return 1;
}

/**
 * gfpm_weighted_sum
 * Compute the weighted sum of two matrices, and store the result in
 * a third one. This third matrix must be pre-allocated.
 * @params:
 *  * dest : the matrix object to store the result into
 *  * left_constant, right_constant : gfp_elements that represent
 *    the field elements to weight the left and right matrices with
 *  * left_matrix, right_matrix : the two matrix objects to add
 *    together.
 * @return
 * 1 if success, 0 otherwise
 */
int gfpm_weighted_sum( gfpmatrix dest, gfp_element left_constant, gfpmatrix left_matrix, gfp_element right_constant, gfpmatrix right_matrix )
{
    unsigned  int i, j;
    unsigned int a;

    #ifdef DEBUG
        if( dest.width != left_matrix.width || left_matrix.width != right_matrix.width || dest.height != left_matrix.height || left_matrix.height != right_matrix.height )
        {
            printf("in gfpm_weighted_sum: trying to add matrices of incompatible dimensions! %ix%i + %ix%i = %ix%i\n", left_matrix.height, left_matrix.width, right_matrix.height, right_matrix.width, dest.height, dest.width);
            return 0;
        }
    #endif
    
    for( i = 0 ; i < dest.height ; ++i )
    {
        for( j = 0 ; j < dest.width ; ++j )
        {
            a = left_matrix.data[i*left_matrix.width + j] * left_constant + right_matrix.data[i*right_matrix.width + j] * right_constant;
            dest.data[i*dest.width + j] = a % GF_PRIME_MODULUS;
        }
    }

    return 1;
}

/**
 * gfpm_rowop
 * Perform a row operation on the given matrix, i.e., add one row,
 * weighted by a constant, to another.
 * @params
 *  * mat : the matrix object to operate on
 *  * destrow, sourcerow : unsigned  ints representing indices of the rows to operate on
 *  * constant : gfp_element representing the right constant
 *  * offset : unsigned  int, represents the number of zeros to skip before applying the row operation
 * @returns
 *  * 1 if success
 */
int gfpm_rowop( gfpmatrix mat, unsigned  int destrow, unsigned  int sourcerow, gfp_element constant, unsigned  int offset )
{
    unsigned  int j;
    for( j = offset ; j < mat.width ; ++j )
    {
        mat.data[destrow*mat.width + j] = (mat.data[destrow*mat.width + j] + mat.data[sourcerow*mat.width + j] * constant) % GF_PRIME_MODULUS;
    }

    return 1;
}

/**
 * gfpm_scalerow
 * Scales a single row in the matrix with a given constant.
 * @params
 *  * mat : the matrix object to operate on
 *  * rowidx : index of the row to scale
 *  * constant : gfp_element -- the field element to multiply the
 *    row with
 * @returns
 *  * 1 if success, 0 otherwise
 */
int gfpm_scalerow( gfpmatrix mat, unsigned  int rowidx, gfp_element constant )
{
    unsigned  int j;
    for( j = 0 ; j < mat.width ; ++j )
    {
        mat.data[rowidx*mat.width + j] = (mat.data[rowidx*mat.width + j] * constant) % GF_PRIME_MODULUS;
    }
    return 1;
}

/**
 * gfpm_fliprows
 * Flip two rows in the given matrix.
 * @params
 *  * mat : the matrix object to operate on
 *  * destrow, sourcerow : the indices of the rows to flip
 * @return
 *  * 1 if success
 */
int gfpm_fliprows( gfpmatrix mat, unsigned  int destrow, unsigned  int sourcerow )
{
    unsigned  int j;
    gfp_element a;
    for( j = 0 ; j < mat.width ; ++j )
    {
        a = mat.data[destrow*mat.width + j];
        mat.data[destrow*mat.width + j] = mat.data[sourcerow*mat.width + j];
        mat.data[sourcerow*mat.width + j] = a;
    }

    return 1;
}

/**
 * gfpm_redech
 * Reduce the given matrix to reduced row echelon form using row
 * operations.
 * @params
 *  * mat : the matrix object to work on
 * @return
 *  1 if success
 */
int gfpm_redech( gfpmatrix mat )
{
    unsigned  int col, row, i;
    gfp_element inv;
    row = 0;
    for( col = 0 ; col < mat.width ; ++col )
    {
        for( i = row ; i < mat.height ; ++i )
        {
            if( mat.data[i*mat.width + col] != 0 )
            {
                if( i != row )
                {
                    gfpm_fliprows(mat, i, row);
                }
                break;
            }
        }

        if( i == mat.height )
        {
            continue;
        }

        inv = gf_inverse(mat.data[row*mat.width + col]);
        
        if( inv != 1 )
        {
            gfpm_scalerow(mat, row, inv);
        }

        for( i = 0 ; i < mat.height ; ++i )
        {
            if( i == row )
            {
                continue;
            }
            gfpm_rowop(mat, i, row, (GF_PRIME_MODULUS - mat.data[i*mat.width + col]) % GF_PRIME_MODULUS, col);
        }

        row = row + 1;

        if( row == mat.height )
        {
            break;
        }
    }

    return 1;
}

/**
 * gfpm_solve
 * Solve a matrix equation of the form Ax = b for x up to a term in
 * the kernel of A. This routine also initializes a kernel matrix,
 * whose rows form a basis for the kernel of A.
 * @params
 *  * coeffs : a mxn gfpmatrix object representing the coefficient
 *    matrix A
 *  * target : a mx1 gfpmatrix object representing the b vector
 *  * solution : a nx1 gfpmatrix object to store one solution into
 *  * kernel : an uninitialized gfpmatrix object whose columns will
 *    span the kernel of A
 * @post
 *  * for all kernel.width x 1 vectors "random" holds:
 *          coeffs * (solution + kernel * random) = target
 * @return
 *  * 1 if a solution exists, 0 otherwise
 */
int gfpm_solve( gfpmatrix coeffs, gfpmatrix target, gfpmatrix solution, gfpmatrix * kernel )
{
    /* declare variables for echelon reduction */
    unsigned  int col, row, i, j;
    gfp_element inv;
    gfpmatrix mat;

    /* declare variables for pivot tracking */
    unsigned  int *pivots;
    unsigned  int *npivots;
    unsigned  int num_pivots;
    unsigned  int num_npivots;
    int have_solution;

    /* initialize variables for pivot tracking */
    num_pivots = 0;
    num_npivots = 0;
    pivots = malloc(sizeof(unsigned  int) * (coeffs.width+1));
    npivots = malloc(sizeof(unsigned  int) * (coeffs.width+1));

    /* initialize mat and copy coeffs and target to it */
    mat = gfpm_init(coeffs.height, coeffs.width+1);
    /*gfpm_copy(mat, coeffs);*/
    for( i = 0 ; i < mat.height ; ++i )
    {
        for( j = 0 ; j < coeffs.width ; ++j )
        {
            mat.data[i*mat.width + j] = coeffs.data[i*coeffs.width + j];
        }
    }
    for( i = 0 ; i < mat.height ; ++i )
    {
        mat.data[i*mat.width + mat.width - 1] = target.data[i*target.width + 0];
    }

    /* perform row echelon reduction */
    row = 0;
    for( col = 0 ; col < mat.width ; ++col )
    {
        for( i = row ; i < mat.height ; ++i )
        {
            if( mat.data[i*mat.width + col] != 0 )
            {
                if( i != row )
                {
                    gfpm_fliprows(mat, i, row);
                }
                break;
            }
        }

        if( i == mat.height )
        {
            /* non pivot */
            npivots[num_npivots++] = col;
            continue;
        }

        /* pivot */
        pivots[num_pivots++] = col;

        inv = gf_inverse(mat.data[row*mat.width + col]);
        
        if( inv != 1 )
        {
            gfpm_scalerow(mat, row, inv);
        }

        for( i = 0 ; i < mat.height ; ++i )
        {
            if( i == row )
            {
                continue;
            }
            gfpm_rowop(mat, i, row, (GF_PRIME_MODULUS - mat.data[i*mat.width + col]) % GF_PRIME_MODULUS, col);
        }

        row = row + 1;

        if( row == mat.height )
        {
            for( i = col+1 ; i < mat.width ; ++i )
            {
                npivots[num_npivots++] = i;
            }
            break;
        }

    }

    /* read out solution if the system is consistent */
    have_solution = (pivots[num_pivots-1] != mat.width-1);
    for( i = 0 ; i < mat.width-1 ; ++i )
    {
        solution.data[i*solution.width] = 0;
    }
    if( have_solution == 1 )
    {
        for( i = 0 ; i < num_pivots ; ++i )
        {
            solution.data[pivots[i]*solution.width] = mat.data[i*mat.width + mat.width - 1];
        }
    }

    /* read out kernel, if it exists */
    if( num_npivots > 1 )
    {
        *kernel = gfpm_init(mat.width-1, num_npivots-1);
        gfpm_zeros(*kernel);
        for( j = 0 ; j < num_npivots-1 ; ++j )
        {
            kernel->data[npivots[j]*kernel->width + j] = GF_PRIME_MODULUS - 1;
            for( i = 0 ; i < num_pivots && pivots[i] < npivots[j] ; ++i )
            {
                kernel->data[pivots[i]*kernel->width + j] = mat.data[i*mat.width + npivots[j]];
            }
        }
    }
    else
    {
        kernel->width = 0;
    }

    /* free allocated memory */
    gfpm_destroy(mat);
    free(pivots);
    free(npivots);

    return have_solution;
}

/**
 * gfpm_inspan
 * Decide if a vector is in the column span of a matrix.
 * @params:
 *  * vec : nx1 matrix
 *  * mat : nxm matrix
 * @returns:
 *  * 1 if vec is in colspan(mat); 0 otherwise
 */
int gfpm_inspan( gfpmatrix vec, gfpmatrix mat )
{
    gfpmatrix solution, kernel;
    int i, j;
    int success;

    #ifdef DEBUG
        if( mat.height != vec.height || vec.width != 1 )
        {
            printf("cannot decide if %ix%i vector is in colspan of %ix%i matrix because of dimension mismatch.\n", vec.height, vec.width, mat.height, mat.width);
            return 0;
        }
    #endif

    /* try and solve the system; if the system is consistent, the vector
     * lies in the span of the matrix */
    solution = gfpm_init(mat.width, 1);
    success = gfpm_solve(mat, vec, solution, &kernel);
    gfpm_destroy(solution);
    if( kernel.width > 0 )
        gfpm_destroy(kernel);

    return success;
}

/**
 * gfpm_stack
 * Stacks one matrix on top of another, and stores the result in the third
 * @params
 *  * mat : matrix object to store the result into
 *  * top, bottom : matrix objects to stack
 * @return
 *  * 1 if success, 0 otherwise
 */
int gfpm_stack( gfpmatrix mat, gfpmatrix top, gfpmatrix bottom )
{
    unsigned  int i, j;

    #ifdef DEBUG
        if( mat.width != top.width || top.width != bottom.width || mat.height != top.height + bottom.height )
        {
            printf("in gfpm_stack: cannot stack matrices of conflicting dimensions! %ix%i stack %ix%i = %ix%i\n", top.height, top.width, bottom.height, bottom.width, mat.height, mat.width);
            return 0;
        }
    #endif
    for( i = 0 ; i < top.height ; ++i )
    {
        for( j = 0 ; j < top.width ; ++j )
        {
            mat.data[i*top.width + j] = top.data[i*top.width + j];
        }
    }
    for( i = 0 ; i < bottom.height ; ++i )
    {
        for( j = 0 ; j < bottom.width ; ++j )
        {
            mat.data[(i+top.height)*mat.width + j] = bottom.data[i*bottom.width + j];
        }
    }

    return 1;
}

/**
 * gfpm_cat
 * Concatenates one matrix to another, and stores the result in a
 * third matrix.
 * @params
 *  * res : matrix object that contains the result
 *  * left, right : matrix objects to be stacked on the left, and
 *    right, respectively
 * @return
 *  * 1 if success, 0 otherwise
 */
int gfpm_cat( gfpmatrix res, gfpmatrix left, gfpmatrix right )
{
    unsigned  int i, j;

    #ifdef DEBUG
        if( res.height != left.height || left.height != right.height || res.width != left.width + right.width )
        {
            printf("in gfpm_cat: cannot concatenate two matrices of conflicting dimensions! %ix%i cat %ix%i = %ix%i\n", left.height, left.width, right.height, right.width, res.height, res.width);
            return 0;
        }
    #endif

    for( i = 0 ; i < res.height ; ++i )
    {
        for( j = 0 ; j < left.width ; ++j )
        {
            res.data[i*res.width + j] = left.data[i*left.width + j];
        }
        for( j = 0 ; j < right.width ; ++j )
        {
            res.data[i*res.width + left.width + j] = right.data[i*right.width + j];
        }
    }
    return 1;
}

/**
 * gfpm_slice
 * Slice a submatrix out of another matrix.
 * @params:
 *  * dest : the matrix to store the result in; the height and width
 *    of this matrix determine the area of elements to be copied
 *  * source : the matrix where the slice comes from
 *  * row_start, col_start : the indices of the row and column at
 *    which to start the slice
 * @return
 *  * 1 if success, 0 otherwise
 */
int gfpm_slice( gfpmatrix dest, gfpmatrix source, unsigned  int row_start, unsigned  int col_start )
{
    unsigned  int i, j;
    #ifdef DEBUG
        if( source.width < col_start + dest.width || source.height < row_start + dest.height )
        {
            printf("in gfpm_slice: cannot grab slice because slice size exceeds bounds! slicing %ix%i submatrix starting at (%i,%i) from %ix%i matrix\n", dest.height, dest.width, row_start, col_start, source.height, source.width);
            return 0;
        }
    #endif
    for( i = 0 ; i < dest.height ; ++i )
    {
        for( j = 0 ; j < dest.width ; ++j )
        {
            dest.data[i*dest.width + j] = source.data[(i+row_start)*source.width + col_start + j];
        }
    }
    return 0;
}

/**
 * gfpm_inverse
 * Compute the matrix inverse of mat, store the result in inv.
 * @return
 *  * 1 if success
 */
int gfpm_inverse( gfpmatrix inv, gfpmatrix mat )
{
    unsigned int i, j;
    unsigned  int catwidth;
    int invertible;
    gfpmatrix concat;
   
    catwidth = inv.width + mat.width;

    /* Set inv to the identity matrix. */
    for( i = 0 ; i < inv.height ; ++i )
    {
        for( j = 0 ; j < inv.width ; ++j )
        {
            inv.data[i*inv.width + j] = 0;
        }
        inv.data[i*inv.width + i] = 1;
    }

    /* Concatenate mat with identity */
    concat = gfpm_init(mat.height, catwidth);
    gfpm_cat(concat, mat, inv);

    /* row-reduce concat to echelon form */
    gfpm_redech(concat);

    /* test if main diagonal has only ones, because otherwise the
     * matrix is not invertible */
    invertible = 1;
    for( i = 0 ; i < inv.height ; ++i )
    {
        invertible = invertible & (int)(concat.data[i*concat.width + i]);
    }

    if( 0 == invertible )
    {
        gfpm_destroy(concat);
        return 0;
    }

    /* select rightmost square from concat */
    for( i = 0 ; i < inv.height ; ++i )
    {
        for( j = 0 ; j < inv.width ; ++j )
        {
            inv.data[i*inv.width + j] = concat.data[i*concat.width + mat.width + j];
        }
    }

    /* free concat */
    gfpm_destroy(concat);

    return 1;
}

/**
 * gfpm_print
 * Use printf to print the matrix to stdout.
 */
int gfpm_print( gfpmatrix mat )
{
    unsigned int i, j;
    for( i = 0 ; i < mat.height ; ++i )
    {
        for( j = 0 ; j < mat.width ; ++j )
        {
            if( mat.data[i*mat.width + j] < 10 )
                printf(" ");
            printf(" %i", mat.data[i*mat.width+j]);
        }
        printf("\n");
    }
}

/**
 * gfpm_print_transpose
 * Use printf to print the transpose of the matrix to stdout.
 */
int gfpm_print_transpose( gfpmatrix mat )
{
    gfpmatrix temp;
    temp = gfpm_clone(mat);
    gfpm_transpose(&temp);
    gfpm_print(temp);
    gfpm_destroy(temp);
    return 1;
}

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
        hqs.quadratic_forms[i].data = malloc(n*n*sizeof(gfp_element));
        hqs.quadratic_forms[i].width = n;
        hqs.quadratic_forms[i].height = n;
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
        free(sys.quadratic_forms[i].data);
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
 *    "large enough" means n*n*m*sizeof(unsigned int)
 * @result
 *  * sys will contain random coefficients
 * @return
 *  * 1 if success, 0 otherwise
 */
int hqs_random( hqsystem sys, unsigned char * randomness )
{
    unsigned int i, j, k, l;
    unsigned  int * rand = (unsigned int *)randomness;

    l = 0;
    for( k = 0 ; k < sys.m ; ++k )
    {
        for( i = 0 ; i < sys.n ; ++i )
        {
            for( j = 0 ; j < sys.n ; ++j )
            {
                sys.quadratic_forms[k].data[i*sys.n + j] = randomness[l] % GF_PRIME_MODULUS;
                l++;
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
                y.data[k*y.width + j] = e.data[0];
            }
        }
    }

    gfpm_destroy(vector);
    gfpm_destroy(temp);

    return 1;
}

