#include <stdio.h>
#include <stdlib.h>
#include "csprng.h"
#include "bi.h"

int test_shift_left( unsigned int * random )
{
    unsigned int k;
    bi a, b;
    csprng rng;
    unsigned char * randomness;

    printf("testing shift left (%u) ...\n", *random);

    csprng_init(&rng);
    csprng_seed(&rng, sizeof(unsigned int), (unsigned char*)random);
    csprng_generate(&rng, sizeof(unsigned int), (unsigned char*)random);
    k = 20 + (*random%500);
    randomness = malloc(k);
    csprng_generate(&rng, k, randomness);

    a = bi_init(k);
    bi_random(&a, k, randomness);
    bi_print_bitstring(a);

    b = bi_init(1);
    bi_shift_left(&b, a, 121);
    bi_print_bitstring(b);

    bi_destroy(a);
    bi_destroy(b);
    free(randomness);
    return 1;
}

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

    *random = 1361598189;

    printf("testing multiplication (%u) ...", *random);
    success = 0;

    csprng_init(&rng);
    csprng_seed(&rng, sizeof(unsigned int), (unsigned char*)random);
    csprng_generate(&rng, sizeof(unsigned int), (unsigned char*) random);

    k = 100 + csprng_generate_ulong(&rng) % 1000;
    l = 100 + csprng_generate_ulong(&rng) % 1000;
    m = 100 + csprng_generate_ulong(&rng) % 1000;

    printf("k: %i, l: %i, m: %i\n", k, l, m);

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

    printf("a: "); bi_print_bitstring(a); printf("\n");
    printf("b: "); bi_print_bitstring(b); printf("\n");
    printf("c: "); bi_print_bitstring(c); printf("\n");

    bi_multiply(&d, a, b);
    printf("a*b: "); bi_print_bitstring(d); printf("\n");

    bi_multiply(&e, a, c);
    printf("a*c: "); bi_print_bitstring(e); printf("\n");

    bi_add(&f, d, e);

    printf("a*b + a*c: "); bi_print_bitstring(f); printf("\n");

    bi_add(&d, b, c);
    printf("b + c: "); bi_print_bitstring(d); printf("\n");

    bi_multiply(&e, a, d);

    printf("a*(b + c): ");bi_print_bitstring(e); printf("\n");

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

int main( int argc, char ** argv )
{

    int i;
    int b;

    unsigned int random;

    random = 1804289383;

    printf("randomness: %u\n", random);

    b = 1;
    for( i = 0 ; i < 0 && b == 1 ; ++i ) b = b & test_shift_left(&random);
    for( i = 0 ; i < 0 && b == 1 ; ++i ) b = b & test_addition(&random);
    for( i = 0 ; i < 100 && b == 1 ; ++i ) b = b & test_multiplication(&random);

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

