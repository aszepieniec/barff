#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>

#include "bi.h"
#include "csprng.h"

int main( int argc, char ** argv )
{

    int num_trials;
    int bitsize, num_limbs;
    int i;

    unsigned int random;

    unsigned long int * randomness;
    unsigned long int * random_ints;
    int certainty;

    csprng rng;

    bi p;

    if( argc < 3 )
    {
        printf("usage: ./benchmark_primality 10 10\n");
        return 0;
    }

    num_trials = atoi(argv[1]);
    bitsize = atoi(argv[2]);

    random = rand();

    printf("benchmarking: sampling %i prime numbers of %i bits with randomness %u \n", num_trials, bitsize, random);

    csprng_init(&rng);
    csprng_seed(&rng, sizeof(unsigned int), (unsigned char*)&random);

    num_limbs = (bitsize + sizeof(unsigned long int)*8 - 1) / (sizeof(unsigned long int)*8);
    randomness = malloc(sizeof(unsigned long int) * num_limbs);
    certainty = 50;
    random_ints = malloc(sizeof(unsigned long int) * certainty);

    p = bi_init(0);

    for( i = 0 ; i < num_trials ; ++i )
    {
        csprng_generate(&rng, sizeof(unsigned long int)*certainty, (unsigned char *)random_ints);
        do
        {
            csprng_generate(&rng, sizeof(unsigned long int) * num_limbs, (unsigned char*)randomness);
            bi_random(&p, bitsize, (unsigned char *)randomness);
            bi_setbit(&p, 0, 1);
        }
        while( bi_is_prime(p, random_ints, certainty) == 0 );
    }

    printf("done.\n");

    if( num_trials == 1 )
    {
        printf("p = int('"); bi_print_bitstring(p); printf("', 2)\n");
    }

    bi_destroy(p);
    free(randomness);

    return 1;
}
