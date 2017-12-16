#include <stdio.h>

#include "gf254.h"
#include "ffa.h"

void gf_out(__m128i a) {
    uint64_t t[2];
    _mm_storeu_si128((__m128i *)(t+0), a);
    printf("%.16lX %.16lX\n", t[1], t[0]);
}

void gf254_out(gf_t a) {
    uint64_t t[4];
    low_rdc(&a[0], &a[1]);
    _mm_storeu_si128((__m128i *)(t+0), a[0]);
    _mm_storeu_si128((__m128i *)(t+2), a[1]);

    printf("0x%.16lX%.16lX, ", t[2], t[0]);
    printf("0x%.16lX%.16lX\n", t[3], t[1]);
}

void gf254_rnd(gf_t a) {
    unsigned long long int t[4];

    __builtin_ia32_rdrand64_step(t + 0);
    __builtin_ia32_rdrand64_step(t + 1);
    __builtin_ia32_rdrand64_step(t + 2);
    __builtin_ia32_rdrand64_step(t + 3);

    t[1] &= 0x7FFFFFFFFFFFFFFF;
    t[3] &= 0x7FFFFFFFFFFFFFFF;
    a[0] = _mm_loadu_si128((__m128i *)(t+0));
    a[1] = _mm_loadu_si128((__m128i *)(t+2));
}

int gf254_cmp(gf_t a, gf_t b) {
    uint64_t t[4], u[4];
    int result = 0;

    _mm_storeu_si128((__m128i *)(t+0), a[0]);
    _mm_storeu_si128((__m128i *)(t+2), a[1]);
    _mm_storeu_si128((__m128i *)(u+0), b[0]);
    _mm_storeu_si128((__m128i *)(u+2), b[1]);

    for (int i = 0; i < 4; i++) {
        result |= (t[i] != u[i]);
    }

    return result != 0;
}

void gf254_cpy(gf_t a, gf_t b) {
    a[0] = b[0];
    a[1] = b[1];
}

void gf254_inp(gf_t a, uint64_t t[]) {
    a[0] = _mm_set_epi64x(t[2], t[0]);
    a[1] = _mm_set_epi64x(t[3], t[1]);
}
