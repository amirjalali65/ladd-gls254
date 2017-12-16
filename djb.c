#include <stdio.h>

#include "gf254.h"
#define INLINE static inline
#include "gf254.c"
#include "ladd.c"

void djb_recoding(int *S0, int *S1, int *S2, int *S3, int *fa, int *fi,
		uint64_t k1[], uint64_t k2[], int s1, int s2) {
	uint64_t mm, nn, mn = (k1[0] & 0x1) ^ (k2[0] & 0x1);

	int i, first_add = 0, final_index = 0;

	if (mn == 1) {
		final_index = 2;
	} else if (k2[0] & 0x1) {
		final_index = 0;
	} else {
		final_index = 1;
	}

	first_add = (k1[0] & 0x1);
	for (i = 0; i < 63; i++) {
		mm = (k1[0] >> i) & 0x1;
		mm ^= (k1[0] >> (i + 1)) & 0x1;

		nn = (k2[0] >> i) & 0x1;
		nn ^= (k2[0] >> (i + 1)) & 0x1;

		mn = mm ^ nn;
		S0[i] = mn;
		S1[i] = mm;
		S2[i] = (k1[0] >> (i + 1)) & 0x1;
		S2[i] ^= (k2[0] >> (i + 1)) & 0x1;
		S3[i] = first_add;
		first_add = mm ^ (mn ^ 1) * first_add;
	}
	mm = (k1[0] >> 63) & 0x1;
	mm ^= k1[1] & 0x1;

	nn = (k2[0] >> 63) & 0x1;
	nn ^= k2[1] & 0x1;

	mn = mm ^ nn;
	S0[i] = mn;
	S1[i] = mm;
	S2[i] = k1[1] & 0x1;
	S2[i] ^= k2[1] & 0x1;
	S3[i] = first_add;
	first_add = mm ^ (mn ^ 1) * first_add;
	for (i = 64; i < 127; i++) {
		mm = (k1[1] >> i) & 0x1;
		mm ^= (k1[1] >> (i + 1)) & 0x1;

		nn = (k2[1] >> i) & 0x1;
		nn ^= (k2[1] >> (i + 1)) & 0x1;

		mn = mm ^ nn;
		S0[i] = mn;
		S1[i] = mm;
		S2[i] = (k1[1] >> (i + 1)) & 0x1;
		S2[i] ^= (k2[1] >> (i + 1)) & 0x1;
		S3[i] = first_add;
		first_add = mm ^ (mn ^ 1) * first_add;
	}

	*fa = first_add;
	*fi = final_index;
}

INLINE void ec_doub_weis(gf_t x3, gf_t z3, gf_t x1, gf_t z1, gf_t b) {
	gf_t t1;
	gf254_sqr(z3, z1);
	gf254_sqr(t1, x1);
	gf254_mul1(x3, z3, b);
	gf254_add(x3, x3, t1);
	gf254_mul(z3, t1, z3);
	gf254_sqr(x3, x3);
}

INLINE void ec_dadd_weis(gf_t x3, gf_t z3, gf_t x1, gf_t z1, gf_t x) {
	gf_t t1, t2;
	gd_t u0, u1;

	gf254_mul(t1, x1, z3);
	gf254_mul(t2, x3, z1);
	gf254_add(z3, t1, t2);
	gf254_sqr(z3, z3);
	gf254_muln(u0, t1, t2);
	gf254_muln(u1, z3, x);
	gf254_add2(u0, u0, u1);
	gf254_rdcn(x3, u0);
}

INLINE void ec_dda_weis(gf_t x0, gf_t z0, gf_t x1, gf_t z1, gf_t x2, gf_t z2, gf_t x, gf_t b) {
    gf_t x3, z3, t1, t2;
    gd_t u0, u1;

    gf254_mul(t1, x0, z1);
    gf254_mul(t2, x1, z0);
    gf254_add(z1, t1, t2);
    gf254_sqr(z1, z1);
    gf254_muln(u0, t1, t2);
    gf254_muln(u1, z1, x);
    gf254_add2(u0, u0, u1);
    gf254_rdcn(x1, u0);

    gf254_sqr(z0, z2);
    gf254_sqr(t1, x2);
    gf254_mul1(x0, z0, b);
    gf254_add(x0, x0, t1);
    gf254_mul(z0, t1, z0);
    gf254_sqr(x0, x0);
}

void djb_weis(gf_t v, gf_t xp, gf_t yp, gf_t b, uint64_t k[]) {
	gf_t t, one, xq, xpq, xpmq, yq, zp, _w[2], tabx[2], taby[2], tabw[5], tabz[5];
	int S0[128], S1[128], S2[128], S3[128];
	int i, j, s1, s2, s, fa, fi, _t;
	uint64_t k1[2], k2[2];
	uint64_t _a[] = { 0, 0, 1, 0 };
	uint64_t _one[] = { 1, 0, 0, 0 };

	gls_recoding(k, k1, k2, &s1, &s2, BETA_E1);
	djb_recoding(S0, S1, S2, S3, &fa, &fi, k1, k2, s1, s2);

	/* PSI */
	tabx[1][0] = _mm_unpacklo_epi64(xp[0], xp[1]);
	tabx[1][1] = _mm_unpackhi_epi64(xp[0], xp[1]);
	taby[1][0] = _mm_unpacklo_epi64(yp[0], yp[1]);
	taby[1][1] = _mm_unpackhi_epi64(yp[0], yp[1]);

	taby[1][0] = _mm_xor_si128(taby[1][0], taby[1][1]);
	taby[1][0] = _mm_xor_si128(tabx[1][0], taby[1][0]);
	tabx[1][0] = _mm_xor_si128(tabx[1][0], tabx[1][1]);
	taby[1][1] = _mm_xor_si128(tabx[1][0], taby[1][1]);

	xq[0] = _mm_unpacklo_epi64(tabx[1][0], tabx[1][1]);
	xq[1] = _mm_unpackhi_epi64(tabx[1][0], tabx[1][1]);
	yq[0] = _mm_unpacklo_epi64(taby[1][0], taby[1][1]);
	yq[1] = _mm_unpackhi_epi64(taby[1][0], taby[1][1]);

	/* Compute P+Q. */
	gf_t t0, t1, l;
	gf254_inp(one, _a);
	gf254_add(t0, xp, xq);
	gf254_add(l, yp, yq);
	gf254_inv(xpmq, t0);
	gf254_mul(l, l, xpmq);
	gf254_sqr(t1, l);
	gf254_add(t1, t1, one);
	gf254_add(t0, t0, l);
	gf254_add(xpq, t0, t1);
	/* Compute P-Q. */
	gf254_add(t0, xp, xq);
	gf254_add(l, yp, yq);
	gf254_add(l, l, xq);
	gf254_mul(l, l, xpmq);
	gf254_sqr(t1, l);
	gf254_add(t1, t1, one);
	gf254_add(t0, t0, l);
	gf254_add(xpmq, t0, t1);

	gf254_inp(one, _one);
	gf254_cpy(tabw[0], xpq);
	gf254_cpy(tabz[0], one);
	gf254_cpy(tabw[1], xpq);
	gf254_cpy(tabz[1], one);
	gf254_cpy(tabw[2], xpq);
	gf254_cpy(tabz[2], one);
	gf254_cpy(zp, one);

	gf254_cswap(xp, xq, fa);
	ec_dadd_weis(tabw[2], tabz[2], xp, zp, xq);
	ec_doub_weis(tabw[1], tabz[1], tabw[1], tabz[1], b);
	gf254_cswap(xp, xq, fa);

	for (j = 126; j >= 0; j--) {
		_t = S1[j] ^ (S3[j] & S0[j]);
		gf254_sel(_w[0], xq, xp, S3[j]);
		gf254_sel(_w[1], xpq, xpmq, S2[j]);
		gf254_sel(tabw[4], tabw[1], tabw[0], _t);
		gf254_sel(tabz[4], tabz[1], tabz[0], _t);
		gf254_sel(tabw[3], tabw[4], tabw[2], S0[j]);
		gf254_sel(tabz[3], tabz[4], tabz[2], S0[j]);
		ec_dadd_weis(tabw[2], tabz[2], tabw[4], tabz[4], _w[0]);
		ec_dda_weis(tabw[1], tabz[1], tabw[0], tabz[0], tabw[3], tabz[3], _w[1],
				b);
	}

	gf254_cpy(v, tabw[fi]);
	gf254_cpy(t, tabz[fi]);
	gf254_inv(t, t);
	gf254_mul(v, v, t);
	gf254_rdc(v);
}

INLINE void ec_add_aff(gf_t x3, gf_t y3, gf_t x1, gf_t y1, gf_t x2,
		gf_t y2, gf_t one) {
	gf_t lam, t0, t1;

	gf254_add(t0, x1, x2);
	gf254_add(lam, y1, y2);

	gf254_inv(t1, t0);
	gf254_mul(lam, lam, t1);

	gf254_sqr(t1, lam);
	gf254_add(t1, t1, one);
	gf254_add(t0, t0, lam);
	gf254_add(x3, t0, t1);
}

INLINE void ec_dadd_huff(gf_t w0, gf_t z0, gf_t w1, gf_t z1, gf_t wp) {
	gf_t t1, t2;
	gd_t u0, u1;
	//B := W0*Z1;
	gf254_mul(t1, w0, z1);
	//C := W1*Z0;
	gf254_mul(t2, w1, z0);
	//W4 := (B + C)^2;
	gf254_add(w0, t1, t2);
	gf254_sqr(w0, w0);
	//Z4 := B*C + wP*W4;
	gf254_muln(u0, t1, t2);
	gf254_muln(u1, wp, w0);
	gf254_add2(u0, u0, u1);
	gf254_rdcn(z0, u0);
}

INLINE void ec_doub_huff(gf_t w0, gf_t z0, gf_t w1, gf_t z1, gf_t gamma) {
	gf_t t0, t1, t2, t3;
	//A := W0*Z0;
	gf254_mul(t0, w1, z1);
	//Z3 := (W0/gamma + Z0)^4;
	gf254_mul1(t3, w1, gamma);
	gf254_add(z0, z1, t3);
	gf254_sqr(z0, z0);
	gf254_sqr(z0, z0);
	//W3 := A^2;
	gf254_sqr(w0, t0);
}

INLINE void ec_dda_huff(gf_t w0, gf_t z0, gf_t w1, gf_t z1, gf_t w2,
		gf_t z2, gf_t wp, gf_t gamma) {
	gf_t t0, t1, t2, t3;
	gd_t u0, u1;

	/* t0 = W0*W1 + Z0*Z1 */
	gf254_muln(u0, w0, w1);
	gf254_muln(u1, z0, z1);
	gf254_add2(u0, u0, u1);
	gf254_rdcn(t0, u0);

	/* t1 = (W0 + Z0)*(W1 + Z1) */
	gf254_add(t2, w0, z0);
	gf254_add(t1, w1, z1);
	gf254_mul(t1, t1, t2);

	/* W3 := gamma * (W2*Z2)^2 */
	gf254_mul(t3, w2, z2);
	gf254_mul1(t3, t3, gamma);

	/* W4 := ((W0 + Z0)*(W1 + Z1) + T)^2 */
	gf254_add(t1, t1, t0);
	gf254_sqr(w1, t1);

	/* Z3 := (W2 + Z2)^4 */
	gf254_add(t2, w2, z2);
	gf254_sqr(t2, t2);
	gf254_sqr(z0, t2);

	/* Z4 := wP * T^2 */
	gf254_sqr(t0, t0);
	gf254_mul(z1, wp, t0);

	gf254_sqr(w0, t3);
}

void djb_huff(gf_t v, gf_t xp, gf_t yp, gf_t gamma, gf_t b, uint64_t k[]) {
	gf_t t, one, xq, yq, tabx[2], taby[2], tabw[5], tabz[5], _w[2];
	gf_t wp, wq, wpq, wpmq, zp, zq, zpq;
	int S0[128], S1[128], S2[128], S3[128];
	int i, j, l, s1, s2, s, fa, fi, _t;
	uint64_t k1[2], k2[2];
	uint64_t _a[] = { 0, 0, 1, 0 };
	uint64_t _one[] = { 1, 0, 0, 0 };

	gls_recoding(k, k1, k2, &s1, &s2, BETA_E2);
	djb_recoding(S0, S1, S2, S3, &fa, &fi, k1, k2, s1, s2);

	/* PSI */
	tabx[1][0] = _mm_unpacklo_epi64(xp[0], xp[1]);
	tabx[1][1] = _mm_unpackhi_epi64(xp[0], xp[1]);
	taby[1][0] = _mm_unpacklo_epi64(yp[0], yp[1]);
	taby[1][1] = _mm_unpackhi_epi64(yp[0], yp[1]);

	taby[1][0] = _mm_xor_si128(taby[1][0], taby[1][1]);
	taby[1][0] = _mm_xor_si128(tabx[1][0], taby[1][0]);
	tabx[1][0] = _mm_xor_si128(tabx[1][0], tabx[1][1]);
	taby[1][1] = _mm_xor_si128(tabx[1][0], taby[1][1]);

	xq[0] = _mm_unpacklo_epi64(tabx[1][0], tabx[1][1]);
	xq[1] = _mm_unpackhi_epi64(tabx[1][0], tabx[1][1]);
	yq[0] = _mm_unpacklo_epi64(taby[1][0], taby[1][1]);
	yq[1] = _mm_unpackhi_epi64(taby[1][0], taby[1][1]);

	gf254_inp(one, _a);
	ec_add_aff(tabx[0], taby[0], xp, yp, xq, yq, one);
	gf254_inp(one, _one);

	gf254_cpy(wp, xp);
	gf254_cpy(wq, xq);
	gf254_cpy(wpq, tabx[0]);

	t[0] = _mm_unpacklo_epi64(wpq[0], wpq[1]);
	t[1] = _mm_unpackhi_epi64(wpq[0], wpq[1]);
	t[0] = _mm_xor_si128(t[0], t[1]);
	wpmq[0] = _mm_unpacklo_epi64(t[0], t[1]);
	wpmq[1] = _mm_unpackhi_epi64(t[0], t[1]);

	gf254_cpy(tabw[0], one);
	gf254_cpy(tabz[0], wpq);
	gf254_cpy(tabw[1], one);
	gf254_cpy(tabz[1], wpq);
	gf254_cpy(tabw[2], one);
	gf254_cpy(tabz[2], wpq);
	gf254_cpy(zp, one);

	gf254_cswap(wp, wq, fa);
	ec_dadd_huff(tabw[2], tabz[2], one, wp, wq);
	ec_doub_huff(tabw[1], tabz[1], tabw[1], tabz[1], gamma);
	gf254_cswap(wp, wq, fa);

	for (j = 126; j >= 0; j--) {
		_t = S1[j] ^ (S3[j] & S0[j]);
		gf254_sel(_w[0], wq, wp, S3[j]);
		gf254_sel(_w[1], wpq, wpmq, S2[j]);
		gf254_sel(tabw[4], tabw[1], tabw[0], _t);
		gf254_sel(tabz[4], tabz[1], tabz[0], _t);
		gf254_sel(tabw[3], tabw[2], tabw[4], S0[j] ^ 1);
		gf254_sel(tabz[3], tabz[2], tabz[4], S0[j] ^ 1);
		ec_dadd_huff(tabw[2], tabz[2], tabw[4], tabz[4], _w[0]);
		ec_dda_huff(tabw[1], tabz[1], tabw[0], tabz[0], tabw[3], tabz[3], _w[1], gamma);
	}
	gf254_cpy(v, tabw[fi]);
	gf254_cpy(t, tabz[fi]);
	gf254_inv(t, t);
	gf254_mul(v, v, t);
	gf254_rdc(v);
}
