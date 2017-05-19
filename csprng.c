#include "csprng.h"

void KeccakF1600_StatePermute(void *state);

/**
 * csprng_init
 * Set the csprng state to all zeros.
 */
int csprng_init( csprng* rng )
{
    unsigned short int i;
    for( i = 0 ; i < csprng_state_size ; ++i )
    {
        rng->state[i] = 0;
    }

    return 1;
}

/**
 * csprng_seed
 * Randomize the csprng's state by feeding it a seed.
 * @params
 *  * rng : the csprng object, which contains a state
 *  * seed_length : unsigned short int describing the length (in
 *    bytes) of the seed
 *  * seed : pointer to the byte array of the seed
 */
int csprng_seed( csprng* rng, unsigned short int seed_length, unsigned char * seed )
{
    unsigned short int i, j;

    /* Absorb all the full input blocks */
    for( i = 0 ; i < seed_length / csprng_seed_rate ; ++i )
    {
        for( j = 0 ; j < csprng_seed_rate ; ++j )
        {
            rng->state[j] ^= seed[i*csprng_seed_rate + j];
        }
        KeccakF1600_StatePermute(rng->state);
    }

    /* Absorb the remainder of the last input block */
    for( j = 0 ; j < seed_length % csprng_seed_rate ; ++j )
    {
        rng->state[j] ^= seed[i*csprng_seed_rate + j];
    }
    KeccakF1600_StatePermute(rng->state);

    return 1;
}

/**
 * csprng_generate
 * Generate a string of random bytes of the given length from the
 * csprng.
 * @params
 *  * rng : the csprng object which contains a randomized state
 *  * buffer_length : the number of random output bytes to produce
 *  * buffer : the memory address to store the desired random bytes
 *    to
 */
int csprng_generate( csprng* rng, unsigned int buffer_length, unsigned char * buffer )
{
    unsigned short int i, j;

    /* squeeze out all the full output blocks */
    for( i = 0 ; i < buffer_length / csprng_gen_rate ; ++i )
    {
        for( j = 0 ; j < csprng_gen_rate ; ++j )
        {
            buffer[i*csprng_gen_rate + j] = rng->state[j];
        }
        KeccakF1600_StatePermute(rng->state);
    }
    /* squeeze out the remaining bytes of the last output block */
    for( j = 0 ; j < buffer_length % csprng_gen_rate ; ++j )
    {
        buffer[i*csprng_gen_rate + j] = rng->state[j];
    }
    KeccakF1600_StatePermute(rng->state);

    return 1;
}

/**
 * csprng_generate_ulong
 * Use the csprng to generate only an unsigned long int.
 */
unsigned long int csprng_generate_ulong( csprng * rng )
{
    unsigned long int ulong;
    csprng_generate(rng, sizeof(unsigned long int), (unsigned char*)&ulong);
    return ulong;
}

