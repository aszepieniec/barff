# barff
Ansi C library for basic algebra routines for prime finite fields, chiefly targeting MQ cryptography.

# Features
 * matrix routines including initialization to random values, multiplication, echelon reduction, and system solving
 * homogeneous quadratic systems, evaluation and composition with linear transforms
 * csprng based on an ANSI C89 compliant version of Keccak
 * big integer routines

# TODO
 * polynomial and extension field arithmetic
 * link big integer routines with finite field arithmetic

# How To
 * compile: `gcc -o test_algebra test_algebra.c bi.c csprng.c Keccak-readable-and-compact-c89.c -ansi -Wpedantic`

