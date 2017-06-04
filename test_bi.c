#include <stdio.h>
#include <stdlib.h>
#include "csprng.h"
#include "bi.h"

int test_addition( unsigned int * random )
{
    unsigned int k, l, m;
    bi a, b, c, d, e, f, g;
    csprng rng;
    unsigned char * randomness;
    unsigned int randint;
    int cmp;

    printf("testing addition (%u) ...", *random);

    csprng_init(&rng);
    csprng_seed(&rng, sizeof(unsigned int), (unsigned char*)random);

    csprng_generate(&rng, sizeof(unsigned int), (unsigned char*)&randint);
    k = 100 + (randint%1000);

    csprng_generate(&rng, sizeof(unsigned int), (unsigned char*)&randint);
    l = 100 + (randint%1000);

    csprng_generate(&rng, sizeof(unsigned int), (unsigned char*)&randint);
    m = 100 + (randint%1000);

    a = bi_init(k);
    randomness = malloc(k);
    csprng_generate(&rng, k, randomness);
    bi_random(&a, k, randomness);
    if(randomness[0]%2 == 1)
    {
        a.sign = -1;
    }
    free(randomness);

    b = bi_init(l);
    randomness = malloc(l);
    csprng_generate(&rng, l, randomness);
    bi_random(&b, l, randomness);
    if(randomness[0]%2 == 1)
    {
        b.sign = -1;
    }
    free(randomness);

    c = bi_init(m);
    randomness = malloc(m);
    csprng_generate(&rng, m, randomness);
    bi_random(&c, m, randomness);
    if(randomness[0]%2 == 1)
    {
        c.sign = -1;
    }
    free(randomness);

    d = bi_init(0);
    e = bi_init(0);
    f = bi_init(0);
    g = bi_init(0);

    bi_add(&d, a, b);
    bi_add(&f, d, c);


    bi_add(&e, b, c);
    bi_add(&g, a, e);

    cmp = bi_compare(f,g);

    bi_destroy(a);
    bi_destroy(b);
    bi_destroy(c);
    bi_destroy(d);
    bi_destroy(e);
    bi_destroy(f);
    bi_destroy(g);

    *random = csprng_generate_ulong(&rng);

    if( cmp == 0 )
    {
        printf("success.\n");
        return 1;
    }
    printf("cmp = %i\n", cmp);
    printf("failure.\n");
    return 0;
}

int test_multiplication( unsigned int * random )
{
    csprng rng;
    unsigned long int k, l, m;
    unsigned char * randomness;
    bi a, b, c, d, e, f;
    int success;

    printf("testing multiplication (%u) ...  ", *random);
    success = 0;

    csprng_init(&rng);
    csprng_seed(&rng, sizeof(unsigned int), (unsigned char*)random);
    csprng_generate(&rng, sizeof(unsigned int), (unsigned char*) random);

    k = 100 + csprng_generate_ulong(&rng) % 1000;
    l = 100 + csprng_generate_ulong(&rng) % 1000;
    m = 100 + csprng_generate_ulong(&rng) % 1000;

    printf("k: %i, l: %i, m: %i  ", k, l, m);

    a = bi_init(k);
    b = bi_init(l);
    c = bi_init(m);
    d = bi_init(0);
    e = bi_init(0);
    f = bi_init(0);

    randomness = malloc(k/8+sizeof(unsigned long int));
    csprng_generate(&rng, k/8+sizeof(unsigned long int), randomness);
    bi_random(&a, k, randomness);
    free(randomness);
    if( csprng_generate_ulong(&rng) % 2 == 1 )
    {
        bi_negate(&a);
    }

    randomness = malloc(l/8+sizeof(unsigned long int));
    csprng_generate(&rng, l/8+sizeof(unsigned long int), randomness);
    bi_random(&b, l, randomness);
    free(randomness);
    if( csprng_generate_ulong(&rng) % 2 == 1 )
    {
        bi_negate(&b);
    }

    randomness = malloc(m/8+sizeof(unsigned long int));
    csprng_generate(&rng, m/8+sizeof(unsigned long int), randomness);
    bi_random(&c, m, randomness);
    free(randomness);
    if( csprng_generate_ulong(&rng) % 2 == 1 )
    {
        bi_negate(&c);
    }

    bi_multiply(&d, a, b);
    bi_multiply(&e, a, c);
    bi_add(&f, d, e);

    bi_add(&d, b, c);
    bi_multiply(&e, a, d);

    if( bi_compare(e, f) == 0 )
    {
        success = 1;
        printf("success.\n");
    }
    else
    {
        printf("failure.\n");
    }

    bi_destroy(a);
    bi_destroy(b);
    bi_destroy(c);
    bi_destroy(d);
    bi_destroy(e);
    bi_destroy(f);

    return success;
}

int test_divide( unsigned int * random )
{
    csprng rng;
    bi divisor, numerator;
    bi quotient, remainder;
    bi product, sum;
    int k, l;
    int cmp, success;
    unsigned char * randomness;

    success = 0;
    printf("testing divide (%u) ... ", *random);

    csprng_init(&rng);
    csprng_seed(&rng, sizeof(unsigned int), (unsigned char*)random);
    csprng_generate(&rng, sizeof(unsigned int), (unsigned char*) random);

    k = 100 + (csprng_generate_ulong(&rng) % 1000);
    l = 100 + (csprng_generate_ulong(&rng) % (k-100));

    printf("k = %i  l = %i ... ", k, l);

    randomness = malloc(k/8 + sizeof(unsigned long int));
    csprng_generate(&rng, k/8 + sizeof(unsigned long int), randomness);
    numerator = bi_init(0);
    bi_random(&numerator, k, randomness);
    free(randomness);
    if( csprng_generate_ulong(&rng) % 2 == 1 )
    {
        bi_negate(&numerator);
    }

    randomness = malloc(l/8 + sizeof(unsigned long int));
    csprng_generate(&rng, l/8 + sizeof(unsigned long int), randomness);
    divisor = bi_init(0);
    bi_random(&divisor, l, randomness);
    free(randomness);
    if( csprng_generate_ulong(&rng) % 2 == 1 )
    {
        bi_negate(&divisor);
    }

    quotient = bi_init(0);
    remainder = bi_init(0);
    product = bi_init(0);
    sum = bi_init(0);

    bi_divide(&quotient, &remainder, numerator, divisor);

    bi_multiply(&product, quotient, divisor);
    bi_add(&sum, product, remainder);

    cmp = bi_compare(sum, numerator);

    if( cmp != 0 || remainder.sign == -1 )
    {
        printf("failure.\n");
        printf("cmp: %i\n", cmp);
        printf("remainder sign: %i\n", remainder.sign);
        printf("sum: "); bi_print_bitstring(sum); printf("\n");
        printf("num: "); bi_print_bitstring(numerator); printf("\n");
    }
    else
    {
        success = 1;
        printf("success.\n");
    }

    bi_destroy(numerator);
    bi_destroy(divisor);
    bi_destroy(quotient);
    bi_destroy(remainder);
    bi_destroy(sum);
    bi_destroy(product);

    return success;
}

int test_xgcd( unsigned int * random )
{
    csprng rng;
    bi a, b, x, y, g;
    bi xa, yb, xayb;
    int k, l;
    int cmp, success;
    unsigned char * randomness;

    success = 0;
    printf("testing xgcd (%u) ... ", *random);

    csprng_init(&rng);
    csprng_seed(&rng, sizeof(unsigned int), (unsigned char*)random);
    csprng_generate(&rng, sizeof(unsigned int), (unsigned char*) random);

    k = 100 + (csprng_generate_ulong(&rng) % 1000);
    l = 100 + (csprng_generate_ulong(&rng) % (k-100));

    printf("k = %i  l = %i ... ", k, l);

    randomness = malloc(k/8 + sizeof(unsigned long int));
    csprng_generate(&rng, k/8 + sizeof(unsigned long int), randomness);
    a = bi_init(0);
    bi_random(&a, k, randomness);
    free(randomness);
    if( csprng_generate_ulong(&rng) % 2 == 1 )
    {
        bi_negate(&a);
    }

    randomness = malloc(l/8 + sizeof(unsigned long int));
    csprng_generate(&rng, l/8 + sizeof(unsigned long int), randomness);
    b = bi_init(0);
    bi_random(&b, l, randomness);
    free(randomness);
    if( csprng_generate_ulong(&rng) % 2 == 1 )
    {
        bi_negate(&b);
    }

    x = bi_init(0);
    y = bi_init(0);
    g = bi_init(0);
    xa = bi_init(0);
    yb = bi_init(0);
    xayb = bi_init(0);

    bi_xgcd(&x, &y, &g, a, b);

    bi_multiply(&xa, x, a);
    bi_multiply(&yb, y, b);
    bi_add(&xayb, xa, yb);

    cmp = bi_compare(xayb, g);

    if( cmp != 0 )
    {
        printf("failure.\n");
    }
    else
    {
        success = 1;
        printf("success.\n");
    }

    bi_destroy(x);
    bi_destroy(y);
    bi_destroy(a);
    bi_destroy(b);
    bi_destroy(g);
    bi_destroy(xa);
    bi_destroy(yb);
    bi_destroy(xayb);

    return success;
}

int test_modexp( unsigned int * random )
{
    csprng rng;
    bi a, b, g, ga, gb, gab, gba, p;
    int k, l;
    int cmp, success;
    unsigned char * randomness;

    success = 0;
    printf("testing modpow (%u) ... ", *random);

    csprng_init(&rng);
    csprng_seed(&rng, sizeof(unsigned int), (unsigned char*)random);
    csprng_generate(&rng, sizeof(unsigned int), (unsigned char*) random);

    p = bi_init(1);
    g = bi_init(1);
    bi_one(&p);
    bi_shift_left(&g, p, 256);
    a = bi_cast(1539);
    bi_subtract(&p, g, a); /* sets p to 2^256 - 1539 */

    k = 100 + (csprng_generate_ulong(&rng) % 100);
    l = 100 + (csprng_generate_ulong(&rng) % (k-100));

    printf("k = %i  l = %i ... ", k, l);

    randomness = malloc(k/8 + sizeof(unsigned long int));
    csprng_generate(&rng, k/8 + sizeof(unsigned long int), randomness);
    bi_random(&a, k, randomness);
    free(randomness);
    if( csprng_generate_ulong(&rng) % 2 == 1 )
    {
        bi_negate(&a);
    }

    randomness = malloc(k/8 + sizeof(unsigned long int));
    csprng_generate(&rng, k/8 + sizeof(unsigned long int), randomness);
    b = bi_init(0);
    bi_random(&b, k, randomness);
    free(randomness);
    if( csprng_generate_ulong(&rng) % 2 == 1 )
    {
        bi_negate(&b);
    }

    randomness = malloc(k/8 + sizeof(unsigned long int));
    csprng_generate(&rng, k/8 + sizeof(unsigned long int), randomness);
    bi_random(&g, k, randomness);
    free(randomness);
    if( csprng_generate_ulong(&rng) % 2 == 1 )
    {
        bi_negate(&g);
    }

    ga = bi_init(0);
    bi_modexp(&ga, g, a, p);

    gb = bi_init(0);
    bi_modexp(&gb, g, b, p);

    gab = bi_init(0);
    bi_modexp(&gab, ga, b, p);

    gba = bi_init(0);
    bi_modexp(&gba, gb, a, p);

    cmp = bi_compare(gab, gba);

    if( cmp == 0 )
    {
        printf("success.\n");
        success = 1;
    }
    else
    {
        printf("p: "); bi_print_bitstring(p); printf("\n");
        printf("a: "); bi_print_bitstring(a); printf("\n");
        printf("b: "); bi_print_bitstring(b); printf("\n");
        printf("g: "); bi_print_bitstring(g); printf("\n");
        printf("g^a: "); bi_print_bitstring(ga); printf("\n");
        printf("g^b: "); bi_print_bitstring(gb); printf("\n");
        printf("g^ab: "); bi_print_bitstring(gab); printf("\n");
        printf("g^ba: "); bi_print_bitstring(gba); printf("\n");
        printf("failure.\n");
    }

    bi_destroy(g);
    bi_destroy(ga);
    bi_destroy(gb);
    bi_destroy(gab);
    bi_destroy(gba);
    bi_destroy(p);
    bi_destroy(a);
    bi_destroy(b);

    return success;
}

int test_naf( unsigned int * random )
{
    csprng rng;
    int k, l;
    int cmp, success;
    bi a, b, c, temp;
    unsigned char * randomness;

    success = 0;
    printf("testing naf (%u) ... ", *random);

    csprng_init(&rng);
    csprng_seed(&rng, sizeof(unsigned int), (unsigned char*)random);
    csprng_generate(&rng, sizeof(unsigned int), (unsigned char*) random);

    k = 100 + (csprng_generate_ulong(&rng) % 1000);

    randomness = malloc(k/8 + sizeof(unsigned long int));
    csprng_generate(&rng, k/8 + sizeof(unsigned long int), randomness);
    a = bi_init(k);
    bi_random(&a, k, randomness);
    free(randomness);

    b = bi_init(0);
    c = bi_init(0);

    bi_naf(&b, &c, a);
    temp = bi_init(0);
    bi_subtract(&temp, b, c);

    cmp = bi_compare(temp, a);

    if( cmp == 0 )
    {
        printf("success.\n");
        success = 1;
    }
    else
    {
        printf("failure.\n");
    }

    bi_destroy(a);
    bi_destroy(b);
    bi_destroy(c);
    bi_destroy(temp);

    return 1;
}

int test_primality( unsigned int * random )
{
    csprng rng;
    bi a, b, p;
    int k, l;
    int cmp, success;
    unsigned char * randomness;
    unsigned long int * random_ints;

    success = 1;
    printf("testing primality test (%u) ... ", *random);

    csprng_init(&rng);
    csprng_seed(&rng, sizeof(unsigned int), (unsigned char*)random);
    csprng_generate(&rng, sizeof(unsigned int), (unsigned char*) random);

    p = bi_init(1);
    b = bi_init(1);
    bi_one(&p);
    bi_shift_left(&b, p, 256);
    a = bi_cast(1539);
    bi_subtract(&p, b, a); /* sets p to 2^256 - 1539 */

    k = 100 + (csprng_generate_ulong(&rng) % 100);
    l = 100 + (csprng_generate_ulong(&rng) % (k-100));

    printf("k = %i  l = %i ... ", k, l);

    randomness = malloc(k/8 + sizeof(unsigned long int));
    csprng_generate(&rng, k/8 + sizeof(unsigned long int), randomness);
    bi_random(&a, k, randomness);
    free(randomness);
    if( csprng_generate_ulong(&rng) % 2 == 1 )
    {
        bi_negate(&a);
    }

    randomness = malloc(l/8 + sizeof(unsigned long int));
    csprng_generate(&rng, l/8 + sizeof(unsigned long int), randomness);
    bi_random(&b, l, randomness);
    free(randomness);
    if( csprng_generate_ulong(&rng) % 2 == 1 )
    {
        bi_negate(&b);
    }

    random_ints = malloc(sizeof(unsigned long int) * 50);
    csprng_generate(&rng, 50 * sizeof(unsigned long int), (unsigned char *)random_ints);
    if( bi_is_prime(p, random_ints, 50) == 0 )
    {
        printf("failure. Did not recognize legit prime.\n");
        success = 0;
    }

    bi_multiply(&p, a, b);
    csprng_generate(&rng, 50 * sizeof(unsigned long int), (unsigned char *)random_ints);
    if( bi_is_prime(p, random_ints, 50) == 1 )
    {
        printf("Failure. Failed to recognize legit composite.\n");
        success = 0;
    }

    bi_destroy(a);
    bi_destroy(b);
    bi_destroy(p);
    free(random_ints);

    if( success == 1 )
    {
        printf("success.\n");
    }

    return success;
}

int main( int argc, char ** argv )
{

    int i;
    int b;

    unsigned int random;

    random = 1804289383;

    printf("randomness: %u\n", random);

    b = 1;
    for( i = 0 ; i < 0 && b == 1 ; ++i ) b = b & test_addition(&random);
    for( i = 0 ; i < 100 && b == 1 ; ++i ) b = b & test_multiplication(&random);
    for( i = 0 ; i < 0 && b == 1 ; ++i ) b = b & test_divide(&random);
    for( i = 0 ; i < 0 && b == 1 ; ++i ) b = b & test_xgcd(&random);
    for( i = 0 ; i < 0 && b == 1 ; ++i ) b = b & test_modexp(&random);
    for( i = 0 ; i < 0 && b == 1 ; ++i ) b = b & test_naf(&random);
    for( i = 0 ; i < 0 && b == 1 ; ++i ) b = b & test_primality(&random);

    if( b == 1 )
    {
        printf("success!\n");
    }
    else
    {
        printf("failure!\n");
    }
    return 1;
}

