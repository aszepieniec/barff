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
unsigned char gf_inverse( unsigned char element )
{
    int a;
    int b;
    int x, y, g;

    a = element;
    b = MOD;
    xgcd(a, b, &x, &y, &g);
    x = (MOD + (x % MOD)) % MOD;
    return x;
}

/**
 * gfm
 * Create gfmatrix object with given buffer.
 */
gfmatrix gfm( unsigned short int height, unsigned short int width, unsigned char* pdata )
{
    gfmatrix mat;
    mat.height = height;
    mat.width = width;
    mat.data = pdata;
    return mat;
}

/**
 * gfm_init
 * Create a field matrix object and allocate space for it.
 */
gfmatrix gfm_init( unsigned short int height, unsigned short int width )
{
    gfmatrix mat;
    mat.width = width;
    mat.height = height;
    mat.data = malloc(width*height);
    return mat;
}

/**
 * gfm_destroy
 * Deallocates space allocated to a field matrix object. Call this
 * function before closing the scope where the field matrix object
 * was initialized.
 * @return
 *  * 1 if success
 */
int gfm_destroy( gfmatrix fm )
{
    free(fm.data);
    fm.width = 0;
    fm.height = 0;
    return 1;
}

/**
 * gfm_copy
 * Copy the contents of one matrix to another. Does not allocate
 * memory for the new object; you must do that yourself! (Or use
 * gfm_copy_new instead.)
 * @promise
 *  * dest and source have the same dimensions
 * @return
 *  * 1 if success
 */
int gfm_copy( gfmatrix dest, gfmatrix source )
{
    unsigned int i, j;
    for( i = 0 ; i < dest.height ; ++i )
    {
        for( j = 0 ; j < dest.width ; ++j )
        {
            dest.data[i*dest.width + j] = source.data[i*source.width + j];
        }
    }

    return 1;
}

/**
 * gfm_copy_new
 * Copy one matrix into a new object. Don't forget to call
 * gfm_destroy at the end of scope!
 */
gfmatrix gfm_copy_new( gfmatrix source )
{
    gfmatrix mat;
    mat = gfm_init(source.height, source.width);
    gfm_copy(mat, source);
    return mat;
}

/**
 * gfm_eye
 * Set a matrix to the identity matrix.
 * @param
 *  * eye : matrix object to set to identity
 * @returns
 *  * 1 if success, 0 otherwise
 */
int gfm_eye( gfmatrix eye )
{
    unsigned short int i, j;
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
int gfm_is_eye( gfmatrix eye )
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
 * gfm_equals
 * Test two matrices for equality.
 * @return
 *  * 1 if equal, 0 otherwise
 */
int gfm_equals( gfmatrix lhs, gfmatrix rhs )
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
 * gfm_zero
 * Sets a matrix to all zeros.
 * @param
 *  * zero : matrix object to set to zero
 * @returns
 *  * 1 if success
 */
int gfm_zeros( gfmatrix zero )
{
    unsigned short int i, j;
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
 * gfm_random
 * Put random values into the matrix.
 * @params
 *  * random : matrix objects with data already allocated and whose
 *    elements are to be assigned random values
 *  * rng : pointer to the csprng object from which to draw the
 *    random numbers
 * @result
 *  * random <-$- matrix_space(random.height, random.width)
 */
int gfm_random( gfmatrix random, csprng * rng )
{
    unsigned short int i, j;
    unsigned int l;
    unsigned short int * randomness;
   
    randomness = malloc(sizeof(unsigned short int) * random.width * random.height);
    csprng_generate(rng, sizeof(unsigned short int) * random.height * random.width, (unsigned char *)randomness);

    l = 0;
    for( i = 0 ; i < random.height ; ++i )
    {
        for( j = 0 ; j < random.width ; ++j )
        {
            random.data[i*random.width + j] = randomness[l++] % MOD;
        }
    }

    free(randomness);

    return 1;
}

/**
 * gfm_random_upper_triangular
 * Set the matrix to a random upper triangular with ones on the
 * diagonal.
 * @params
 *  * random : matrix objects with data already allocated and whose
 *    elements above the diagonal are to be assigned random values;
 *    whereas the elements above the diagonal are to be 0 and the
 *    elements on the diagonal are to be 1.
 *  * rng : pointer to the csprng object from which to draw the
 *    random numbers
 * @result
 *  * random <-$- matrix_space(random.height, random.width)
 *    subject to
 *    forall i . random[i,i] = 1
 *    forall i, j>i . random[i,j] = 0
 */
int gfm_random_upper_triangular( gfmatrix random, csprng * rng )
{
    unsigned short int i, j;
    unsigned int l;
    unsigned short int * randomness;

    randomness = malloc(sizeof(unsigned short int) * random.height * random.width);
    csprng_generate(rng, sizeof(unsigned short int) * random.height * random.width, (unsigned char *)randomness);

    l = 0;
    for( i = 0 ; i < random.height ; ++i )
    {
        for( j = 0 ; j < i ; ++j )
        {
            random.data[i*random.width + j] = randomness[l++] % MOD;
        }
        random.data[i*random.width + i] = 1;
        for( j = i+1 ; j < random.width ; ++j )
        {
            random.data[i*random.width + j] = randomness[l++] % MOD;
        }
    }

    free(randomness);

    return 1;
}

/**
 * gfm_random_invertible
 * Generate a random invertible matrix. This is accomplished by first
 * generating two random triangular matrices, where one has ones on
 * the diagonal is upper-triangular and the other has random nonzero
 * elements on the diagonal in addition to being lower-triangular.
 * The generated matrix is the product of the two triangular ones.
 * @params
 *  * mat : the matrix object in which to store the generated random
 *    invertible matrix
 *  * rng : a pointer to the csprng object to generate the random
 *    numbers from
 * @pre
 *  * mat.height = mat.width
 * @result
 *  * mat <-$- GL(mat.height)
 * @return
 *  * 1 if success, 0 otherwise
 */
int gfm_random_invertible( gfmatrix mat, csprng * rng )
{
    gfmatrix utm, ltm;
    unsigned int offset;
    unsigned int i;

    unsigned short int * randomness;

    randomness = malloc(sizeof(unsigned short int) * mat.height);
    csprng_generate(rng, sizeof(unsigned short int) * mat.height, (unsigned char *)randomness);

#ifdef DEBUG
    if( mat.height != mat.width )
    {
        printf("gfm_random_invertible: cannot generate random invertible matrix because matrix not square! %ix%i\n", mat.height, mat.width);
        return 0;
    }
#endif

    utm = gfm_init(mat.height, mat.width);
    ltm = gfm_init(mat.height, mat.width);

    /* generate triangular matrices with ones on the diagonal */
    gfm_random_upper_triangular(utm, rng);
    gfm_random_upper_triangular(ltm, rng);
    gfm_transpose(&ltm);

    /* set the diagonal elements of one matrix to random nonzero
     * elements */
    for( i = 0 ; i < utm.height ; ++i )
    {
        utm.data[i*utm.width + i] = 1 + (randomness[i] % (MOD - 1));
    }

    /* multiply L * U to get the random invertible matrix */
    gfm_multiply(mat, ltm, utm);

    gfm_destroy(ltm);
    gfm_destroy(utm);

    free(randomness);

    return 1;
}

/**
 * gfm_transpose_square
 * Perform a matrix transposition in situ.
 */
int gfm_transpose( gfmatrix * trans )
{
    unsigned char a;
    unsigned short int i, j;

    if( trans->width > trans->height )
    {
        for( i = 0 ; i < trans->height ; ++i )
        {
            for( j = i+1 ; j < trans->width ; ++j )
            {
                a = trans->data[i*trans->width + j];
                trans->data[i*trans->width + j] = trans->data[j*trans->width + i];
                trans->data[j*trans->width + i] = a;
            }
        }
    }
    else
    {
        for( i = 0 ; i < trans->height ; ++i )
        {
            for( j = 0 ; j < i && j < trans->width ; ++j )
            {
                a = trans->data[i*trans->width + j];
                trans->data[i*trans->width + j] = trans->data[j*trans->width + i];
                trans->data[j*trans->width + i] = a;
            }
        }
    }

    a = trans->width;
    trans->width = trans->height;
    trans->height = a;

    return 1;
}

/**
 * gfm_multiply
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
 * that can be stored in a short unsigned int.
 */
int gfm_multiply( gfmatrix dest, gfmatrix left, gfmatrix right )
{
    unsigned short int i, j, k;
    unsigned int acc;

    #ifdef DEBUG
        if( dest.height != left.height || dest.width != right.width || left.width != right.height )
        {
            printf("in gfm_multiply: trying to multiply matrices with unmatched dimensions: %ix%i * %ix%i = %ix%i\n", left.height, left.width, right.height, right.width, dest.height, dest.width);
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
            dest.data[i*dest.width + j] = acc % MOD;
        }
    }

    return 1;
}

/**
 * gfm_multiply_constant
 * Multiply the matrix with a constant.
 * @return
 *  * 1 if success
 */
int gfm_multiply_constant( gfmatrix dest, gfmatrix source, unsigned char constant )
{
    unsigned short int i, j;

#ifdef DEBUG
    if( dest.width != source.width || dest.height != source.height )
    {
        printf("gfm_multiply_constant: cannot multiply matrix with constant because dimensions of destination do not match those of source! %ix%i <- %ix%i\n", dest.height, dest.width, source.height, source.width);
        return 0;
    }
#endif

    for( i = 0 ; i < dest.height ; ++i )
    {
        for( j = 0 ; j < dest.width ; ++j )
        {
            dest.data[i*dest.width + j] = (source.data[i*dest.width + j] * constant) % MOD;
        }
    }
    return 1;
}

/**
 * gfm_add
 * Add one matrix to another and store the result in a third. The
 * third matrix must be preallocated.
 * @params
 *  * dest : the matrix object to store the result in
 *  * left, right : the matrix objects to add together; they must 
 *    have the same dimensions!
 * @return
 *  * 1 if success, 0 otherwise
 */
int gfm_add( gfmatrix dest, gfmatrix left, gfmatrix right )
{
    unsigned short int i, j;

    #ifdef DEBUG
        if( dest.width != left.width || left.width != right.width || dest.height != left.height || left.height != right.height )
        {
            printf("in gfm_add: trying to add matrices of incompatible dimensions! %ix%i + %ix%i = %ix%i\n", left.height, left.width, right.height, right.width, dest.height, dest.width);
            return 0;
        }
    #endif
    
    for( i = 0 ; i < dest.height ; ++i )
    {
        for( j = 0 ; j < dest.width ; ++j )
        {
            dest.data[i*dest.width + j] = (left.data[i*left.width + j] + right.data[i*right.width + j]) % MOD;
        }
    }

    return 1;
}

/**
 * gfm_weighted_sum
 * Compute the weighted sum of two matrices, and store the result in
 * a third one. This third matrix must be pre-allocated.
 * @params:
 *  * dest : the matrix object to store the result into
 *  * left_constant, right_constant : unsigned chars that represent
 *    the field elements to weight the left and right matrices with
 *  * left_matrix, right_matrix : the two matrix objects to add
 *    together.
 * @return
 * 1 if success, 0 otherwise
 */
int gfm_weighted_sum( gfmatrix dest, unsigned char left_constant, gfmatrix left_matrix, unsigned char right_constant, gfmatrix right_matrix )
{
    unsigned short int i, j;
    unsigned int a;

    #ifdef DEBUG
        if( dest.width != left_matrix.width || left_matrix.width != right_matrix.width || dest.height != left_matrix.height || left_matrix.height != right_matrix.height )
        {
            printf("in gfm_weighted_sum: trying to add matrices of incompatible dimensions! %ix%i + %ix%i = %ix%i\n", left_matrix.height, left_matrix.width, right_matrix.height, right_matrix.width, dest.height, dest.width);
            return 0;
        }
    #endif
    
    for( i = 0 ; i < dest.height ; ++i )
    {
        for( j = 0 ; j < dest.width ; ++j )
        {
            a = left_matrix.data[i*left_matrix.width + j] * left_constant + right_matrix.data[i*right_matrix.width + j] * right_constant;
            dest.data[i*dest.width + j] = a % MOD;
        }
    }

    return 1;
}

/**
 * gfm_rowop
 * Perform a row operation on the given matrix, i.e., add one row,
 * weighted by a constant, to another.
 * @params
 *  * mat : the matrix object to operate on
 *  * destrow, sourcerow : unsigned short ints representing indices of the rows to operate on
 *  * constant : unsigned char representing the right constant
 *  * offset : unsigned short int, represents the number of zeros to skip before applying the row operation
 * @returns
 *  * 1 if success
 */
int gfm_rowop( gfmatrix mat, unsigned short int destrow, unsigned short int sourcerow, unsigned char constant, unsigned short int offset )
{
    unsigned short int j;
    for( j = offset ; j < mat.width ; ++j )
    {
        mat.data[destrow*mat.width + j] = (mat.data[destrow*mat.width + j] + mat.data[sourcerow*mat.width + j] * constant) % MOD;
    }

    return 1;
}

/**
 * gfm_scalerow
 * Scales a single row in the matrix with a given constant.
 * @params
 *  * mat : the matrix object to operate on
 *  * rowidx : index of the row to scale
 *  * constant : unsigned char -- the field element to multiply the
 *    row with
 * @returns
 *  * 1 if success, 0 otherwise
 */
int gfm_scalerow( gfmatrix mat, unsigned short int rowidx, unsigned char constant )
{
    unsigned short int j;
    for( j = 0 ; j < mat.width ; ++j )
    {
        mat.data[rowidx*mat.width + j] = (mat.data[rowidx*mat.width + j] * constant) % MOD;
    }
    return 1;
}

/**
 * gfm_fliprows
 * Flip two rows in the given matrix.
 * @params
 *  * mat : the matrix object to operate on
 *  * destrow, sourcerow : the indices of the rows to flip
 * @return
 *  * 1 if success
 */
int gfm_fliprows( gfmatrix mat, unsigned short int destrow, unsigned short int sourcerow )
{
    unsigned short int j;
    unsigned char a;
    for( j = 0 ; j < mat.width ; ++j )
    {
        a = mat.data[destrow*mat.width + j];
        mat.data[destrow*mat.width + j] = mat.data[sourcerow*mat.width + j];
        mat.data[sourcerow*mat.width + j] = a;
    }

    return 1;
}

/**
 * gfm_redech
 * Reduce the given matrix to reduced row echelon form using row
 * operations.
 * @params
 *  * mat : the matrix object to work on
 * @return
 *  1 if success
 */
int gfm_redech( gfmatrix mat )
{
    unsigned short int col, row, i;
    unsigned char inv;
    row = 0;
    for( col = 0 ; col < mat.width ; ++col )
    {
        for( i = row ; i < mat.height ; ++i )
        {
            if( mat.data[i*mat.width + col] != 0 )
            {
                if( i != row )
                {
                    gfm_fliprows(mat, i, row);
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
            gfm_scalerow(mat, row, inv);
        }

        for( i = 0 ; i < mat.height ; ++i )
        {
            if( i == row )
            {
                continue;
            }
            gfm_rowop(mat, i, row, (MOD - mat.data[i*mat.width + col]) % MOD, col);
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
 * gfm_stack
 * Stacks one matrix on top of another, and stores the result in the third
 * @params
 *  * mat : matrix object to store the result into
 *  * top, bottom : matrix objects to stack
 * @return
 *  * 1 if success, 0 otherwise
 */
int gfm_stack( gfmatrix mat, gfmatrix top, gfmatrix bottom )
{
    unsigned short int i, j;

    #ifdef DEBUG
        if( mat.width != top.width || top.width != bottom.width || mat.height != top.height + bottom.height )
        {
            printf("in gfm_stack: cannot stack matrices of conflicting dimensions! %ix%i stack %ix%i = %ix%i\n", top.height, top.width, bottom.height, bottom.width, mat.height, mat.width);
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
 * gfm_cat
 * Concatenates one matrix to another, and stores the result in a
 * third matrix.
 * @params
 *  * res : matrix object that contains the result
 *  * left, right : matrix objects to be stacked on the left, and
 *    right, respectively
 * @return
 *  * 1 if success, 0 otherwise
 */
int gfm_cat( gfmatrix res, gfmatrix left, gfmatrix right )
{
    unsigned short int i, j;

    #ifdef DEBUG
        if( res.height != left.height || left.height != right.height || res.width != left.width + right.width )
        {
            printf("in gfm_cat: cannot concatenate two matrices of conflicting dimensions! %ix%i cat %ix%i = %ix%i\n", left.height, left.width, right.height, right.width, res.height, res.width);
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
 * gfm_slice
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
int gfm_slice( gfmatrix dest, gfmatrix source, unsigned short int row_start, unsigned short int col_start )
{
    unsigned short int i, j;
    #ifdef DEBUG
        if( source.width < col_start + dest.width || source.height < row_start + dest.height )
        {
            printf("in gfm_slice: cannot grab slice because slice size exceeds bounds! slicing %ix%i submatrix starting at (%i,%i) from %ix%i matrix\n", dest.height, dest.width, row_start, col_start, source.height, source.width);
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
 * gfm_inverse
 * Compute the matrix inverse of mat, store the result in inv.
 * @return
 *  * 1 if success
 */
int gfm_inverse( gfmatrix inv, gfmatrix mat )
{
    unsigned int i, j;
    unsigned short int catwidth;
    int invertible;
    gfmatrix concat;
   
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
    concat = gfm_init(mat.height, catwidth);
    gfm_cat(concat, mat, inv);

    /* row-reduce concat to echelon form */
    gfm_redech(concat);

    /* test if main diagonal has only ones, because otherwise the
     * matrix is not invertible */
    invertible = 1;
    for( i = 0 ; i < inv.height ; ++i )
    {
        invertible = invertible & (int)(concat.data[i*concat.width + i]);
    }

    if( !invertible )
    {
        gfm_destroy(concat);
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
    gfm_destroy(concat);

    return 1;
}

/**
 * gfm_print
 * Use printf to print the matrix to stdout.
 */
int gfm_print( gfmatrix mat )
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
 * gfm_print_transpose
 * Use printf to print the transpose of the matrix to stdout.
 */
int gfm_print_transpose( gfmatrix mat )
{
    gfmatrix temp;
    temp = gfm_copy_new(mat);
    gfm_transpose(&temp);
    gfm_print(temp);
    gfm_destroy(temp);
    return 1;
}

/**
 * hqs
 * Create homogeneous quadratic system object from a list of
 * of quadratic forms and the system's dimensions.
 */
hqsystem hqs( gfmatrix* qfs, unsigned short int n, unsigned short int m )
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
hqsystem hqs_init( unsigned short int n, unsigned short int m )
{
    unsigned int i;
    hqsystem hqs;
    hqs.n = n;
    hqs.m = m;
    hqs.quadratic_forms = malloc(sizeof(gfmatrix) * m);
    for( i = 0 ; i < m ; ++i )
    {
        hqs.quadratic_forms[i].data = malloc(n*n);
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
        gfm_copy(dest.quadratic_forms[i], source.quadratic_forms[i]);
    }
    return 1;
}

/**
 * hqs_copy_new
 * Copy one homogeneous quadratic system to a new one. Remember to
 * destroy it when scope ends!
 * @params
 *  * source : the homogeneous quadratic system to copy
 * @return
 *  * dest : a new homogeneous quadratic system identical to source
 */
hqsystem hqs_copy_new( hqsystem source )
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
 *  * rng : pointer to the csprng object to draw the random numbers
 *    from
 * @result
 *  * sys will contain random coefficients
 * @return
 *  * 1 if success, 0 otherwise
 */
int hqs_random( hqsystem sys, csprng * rng )
{
    unsigned int i, j, k, l;
    unsigned short int * randomness;

    randomness = malloc(sizeof(unsigned short int) * sys.m * sys.n * sys.n);
    csprng_generate(rng, sizeof(unsigned short int) * sys.m * sys.n * sys.n, (unsigned char *)randomness);

    l = 0;
    for( k = 0 ; k < sys.m ; ++k )
    {
        for( i = 0 ; i < sys.n ; ++i )
        {
            for( j = 0 ; j < sys.n ; ++j )
            {
                sys.quadratic_forms[k].data[i*sys.n + j] = randomness[l++] % MOD;
            }
        }
    }

    free(randomness);

    return 1;
}

/**
 * hqs_compose_output
 * Compose a homogeneous quadratic system on the left with a
 * linear transform.
 * @params
 *  * T : a gfmatrix object of dimensions mxm
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
int hqs_compose_output( gfmatrix T, hqsystem F )
{
    unsigned int i, j;
    /* declare helper variables */
    hqsystem P;
    gfmatrix temp;

#ifdef DEBUG
    if( T.width != F.m )
    {
        printf("in hqs_compose_output: cannot compose homogeneous quadratic system with linear transform on output side because dimension mismatch! T: %ix%i vs F: %i -> %i\n", T.height, T.width, F.n, F.m);
        return 0;
    }
#endif

    /* init helper variables */
    P = hqs_init(F.n, T.height);
    temp = gfm_init(F.n, F.n);

    /* perform multiplication */
    for( i = 0 ; i < T.height ; ++i )
    {
        gfm_zeros(P.quadratic_forms[i]);
        for( j = 0 ; j < T.width ; ++j )
        {
            gfm_multiply_constant(temp, F.quadratic_forms[j], T.data[i*T.width + j]);
            gfm_add(P.quadratic_forms[i], P.quadratic_forms[i], temp);
        }
    }

    /* copy to argument */
    hqs_copy(F, P);

    /* destroy helper variables */
    hqs_destroy(P);
    gfm_destroy(temp);

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
int hqs_compose_input( hqsystem F, gfmatrix S )
{
    unsigned int i;

    /* declare helper variables */
    gfmatrix temp;
    gfmatrix ST;

    /* debug stuff */
#ifdef DEBUG
    if( F.n != S.height || S.height != S.width )
    {
        printf("hqs_compose_input: cannot compose homogeneous quadratic system with linear transform on input side because of dimension mismatch! F : %i -> %i vs. S : %ix%i\n", F.n, F.m, S.height, S.width);
        return 0;
    }
#endif

    /* init helper variables */

    temp = gfm_init(S.height, S.height);
    ST = gfm_init(S.height, S.width);
    gfm_copy(ST, S);
    gfm_transpose(&ST);

    /* perform multiplication */
    for( i = 0 ; i < F.m ; ++i )
    {
        gfm_multiply(temp, ST, F.quadratic_forms[i]);
        gfm_multiply(F.quadratic_forms[i], temp, S);
    }

    /* destroy helper variables */
    gfm_destroy(temp);
    gfm_destroy(ST);

    return 1;
}

/**
 * hqs_eval
 * Evaluate a homogeneous quadratic system in a vector or, by
 * treating the columns as a list of vectors, as a matrix.
 */
int hqs_eval( gfmatrix y, hqsystem sys, gfmatrix x )
{
    unsigned int i, j, k;
    gfmatrix vector, transposed_vector, temp, e;
    unsigned char edata;

#ifdef DEBUG
    if( y.height != sys.m || sys.n != x.height || y.width != x.width )
    {
        printf("hqs_eval: cannot evaluate quadratic system because of dimension mismatch! F: %i -> %i vs. in: %ix%i, out: %ix%i\n", sys.n, sys.m, x.height, x.width, y.height, y.width);
        return 0;
    }
#endif

    vector = gfm_init(sys.n, 1);
    transposed_vector = gfm(1, sys.n, vector.data);
    temp = gfm_init(sys.n, 1);
    e = gfm(1, 1, &edata);

    for( j = 0 ; j < x.width ; ++j )
    {
        gfm_slice(vector, x, 0, j);
        for( i = 0 ; i < x.height ; ++i )
        {
            for( k = 0 ; k < sys.m ; ++k )
            {
                gfm_multiply(temp, sys.quadratic_forms[k], vector);
                gfm_multiply(e, transposed_vector, temp);
                y.data[k*y.width + j] = e.data[0];
            }
        }
    }

    gfm_destroy(vector);
    gfm_destroy(temp);

    return 1;
}

