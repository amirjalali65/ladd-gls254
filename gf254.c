#include "gf254.h"
#include "ffa.h"

INLINE void gf254_add(gf_t c, gf_t a, gf_t b) {
    c[0] = _mm_xor_si128(a[0], b[0]);
    c[1] = _mm_xor_si128(a[1], b[1]);
}

INLINE void gf254_add2(gd_t c, gd_t a, gd_t b) {
    c[0] = _mm_xor_si128(a[0], b[0]);
    c[1] = _mm_xor_si128(a[1], b[1]);
    c[2] = _mm_xor_si128(a[2], b[2]);
    c[3] = _mm_xor_si128(a[3], b[3]);
}

INLINE void gf254_mul(gf_t c, gf_t a, gf_t b) {
    low_mul(&c[0], &c[1], a[0], a[1], b[0], b[1]);
}

INLINE void gf254_mul_0u(gf_t c, gf_t a) {
    low_mul_00u(&c[0], &c[1], a[0], a[1]);
}

INLINE void gf254_mul_1u(gf_t c, gf_t a) {
    low_mul_01u(&c[0], &c[1], a[0], a[1]);
}

INLINE void gf254_muln(gd_t c, gf_t a, gf_t b) {
    low_muln(&c[0], &c[1], &c[2], &c[3], a[0], a[1], b[0], b[1]);
}

INLINE void gf254_rdcn(gf_t a, gd_t c) {
    low_rdcn(&a[0], &a[1], c[0], c[1], c[2], c[3]);
}

INLINE void gf254_rdc(gf_t c) {
    low_rdc(&c[0], &c[1]);
}

INLINE void gf254_mul1(gf_t c, gf_t a, gf_t b) {
    low_mul_fq1(&c[0], &c[1], a[0], a[1], b[0]);
}

INLINE void gf254_mul_gamma(gf_t c, gf_t a) {
    low_mul_gamma(&c[0], &c[1], a[0], a[1]);
}

INLINE void gf254_lsh(gf_t c, gf_t a) {
    low_lsh(&c[0], &c[1], a[0], a[1]);
}

INLINE void gf254_sqr(gf_t c, gf_t a) {
    low_sqr(&c[0], &c[1], a[0], a[1]);
}

INLINE void gf254_inv(gf_t c, gf_t a) {
    low_inv(&c[0], &c[1], a[0], a[1]);
}

INLINE void gf254_cswap(gf_t a, gf_t b, uint64_t cond) {
    __m128i t, mask = _mm_set_epi64x(-cond, -cond);
    t = a[0];
    a[0] = _mm_blendv_epi8(a[0], b[0], mask);
    b[0] = _mm_blendv_epi8(b[0], t, mask);
    t = a[1];
    a[1] = _mm_blendv_epi8(a[1], b[1], mask);
    b[1] = _mm_blendv_epi8(b[1], t, mask);
    return;
}

INLINE void gf254_sel(gf_t c, gf_t a, gf_t b, uint64_t cond) {
    __m128i t, mask = _mm_set_epi64x(-cond, -cond);
    c[0] = _mm_blendv_epi8(a[0], b[0], mask);
    c[1] = _mm_blendv_epi8(a[1], b[1], mask);
    return;
}

INLINE void gf254_ccopy(gf_t c, gf_t a, uint64_t cond) {
    __m128i t, mask = _mm_set_epi64x(-cond, -cond);
    c[0] = _mm_blendv_epi8(c[0], a[0], mask);
    c[1] = _mm_blendv_epi8(c[1], a[1], mask);
    return;
}
