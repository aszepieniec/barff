# barff
Ansi C library for basic algebra routines for prime finite fields, for instance for MQ cryptography; or for error correcting codes

# Features
 * matrix routines including initialization to random values, multiplication, echelon reduction, and system solving
 * homogeneous quadratic systems, evaluation and composition with linear transforms
 * csprng based on an ANSI C89 compliant version of Keccak
 * big integer routines

# TODO
 * polynomial and extension field arithmetic
 * make primality test faster

# How To
 * compile:
   * standard unit tests: `gcc -o test_algebra test_algebra.c gfp.c csprng.c Keccak-readable-and-compact-c89.c gf*x.c gfpm.c hqs.c -ansi -Wpedantic`
   * big field unit tests: `gcc -o test_algebra test_algebra.c gfbi.c bi.c csprng.c Keccak-readable-and-compact-c89.c gf*x.c gfpm.c hqs.c -DBIG -ansi -Wpedantic`
   * big number arithmetic: `gcc -o test_bi test_bi.c bi.c csprng.c Keccak-readable-and-compact-c89.c -ansi -Wpedantic`
   * prime sampling benchmark: `gcc -o benchmark_primality benchmark_primality.c bi.c csprng.c Keccak-readable-and-compact-c89.c -ansi -Wpedantic -O3`

