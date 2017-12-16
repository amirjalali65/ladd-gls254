#include <x86intrin.h>
#include <stdint.h>

typedef __m128i gf_t[2];
typedef __m128i gd_t[4];

void gf_out(__m128i a);
void gf254_rnd(gf_t a);
void gf254_out(gf_t a);
void gf254_cpy(gf_t a, gf_t b);
void gf254_inp(gf_t a, uint64_t t[]);
int  gf254_cmp(gf_t a, gf_t b);

#define BETA_E1 0x8CCD57A68C1BF773
#define BETA_E2 0X2826AEC5683DD7BF