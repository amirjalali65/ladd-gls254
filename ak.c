#include <stdio.h>

#include "gf254.h"
#define INLINE static inline
#include "gf254.c"
#include "ladd.c"

#define CMP_128(cmp, c0, c1, a0, a1)\
        asm ("subq %3, %0 \n\t"\
             "sbbq %4, %1 \n\t"\
             "sbbq $0, %2 \n\t"\
        : "+r" (c0), "+r" (c1), "+r" (cmp)\
        : "r" (a0), "r" (a1) \
        );

#define R1 0
#define R2 1
#define _R1 2
#define _R2 3

static inline void hswap(uint64_t * c, uint64_t * a, uint64_t cond) {
	uint64_t mask, t;

	mask = -cond;
	t = (a[0] ^ c[0]) & mask;
	a[0] ^= t;
	c[0] ^= t;
	t = (a[1] ^ c[1]) & mask;
	a[1] ^= t;
	c[1] ^= t;
}

#define SEL(a, b, c) ((-c & ((a) ^ (b))) ^ (a))

void ak_recoding(int *S, int *l, uint64_t _k1[], uint64_t _k2[], int s1, int s2) {
	uint64_t k1[2], k2[2], tmp[2], c, d, i = 0;;
	int t, _t;
    k1[0] = _k1[0];
    k1[1] = _k1[1];
    k2[0] = _k2[0];
    k2[1] = _k2[1];

	while ((k1[1] ^ k2[1]) | (k1[0] ^ k2[0])) {
		/* k1 = d, k2 = e. */
		t = (k1[0] ^ k2[0]) & 0x1;
		_t = k1[0] & 0x1;

		c = 1;
		tmp[0] = k1[0];
		tmp[1] = k1[1];
		CMP_128(c, tmp[0], tmp[1], k2[0], k2[1]);
		tmp[0] = SEL(-tmp[0], tmp[0], c);
		tmp[1] = SEL(~tmp[1], tmp[1], c);

		// c = (d > e), t = !(d \equiv e) mod 2, _t  = d mod 2
		d = (t & _t) | ((c ^ 1) & (t ^ 1));
		hswap(k1, k2, d);
		//types_copy(tmp, k1);
		//SUB_128(tmp[0], tmp[1], k2[0], k2[1]);
		//tmp[0] = (c^1)*-tmp[0] + c*tmp[0];
		//tmp[1] = (c^1)*~tmp[1] + c*tmp[1];
		//ccopy(k1, tmp, t);
		k1[0] = SEL(tmp[0], k1[0], t);
		k1[1] = SEL(tmp[1], k1[1], t);
		k1[0] = (k1[0] >> 1) ^ (k1[1] << 63);
		k1[1] >>= 1;
		hswap(k1, k2, d);
		//S[i++] = t * 2 * (c^1) + (t ^ 1)*(2 * _t + R2); /* 2 ou 1,3 */
		S[i++] = SEL((c ^ 1) << (uint64_t) 1, (_t << (uint64_t) 1) + R2, t);
		//printf("%d %X %X %X %X\n", i, k1[0], k1[1], k2[0], k2[1]);
	}
	*l = i;
	//return k1[0];
}

void ec_add_mix_LD(gf_t x3, gf_t y3, gf_t z3, gf_t x2, gf_t y2, gf_t z2,
		gf_t x1, gf_t y1) {
	gf_t t0, t1, t2;
	gf_t x2odd, y2odd;

	gf254_cpy(x2odd, x2);
	gf254_cpy(y2odd, y2);

	gf254_mul(t0, z2, x1);
	gf254_sqr(t1, z2);
	gf254_add(x3, x2odd, t0);
	gf254_mul(t0, z2, x3);
	gf254_mul(t2, t1, y1);
	gf254_add(y3, y2odd, t2);

	gf254_sqr(z3, t0);
	gf254_mul(t2, t0, y3);

	gf254_mul_0u(t1, t1);
	gf254_add(t0, t0, t1);

	gf254_sqr(t1, x3);
	gf254_mul(x3, t1, t0);
	gf254_sqr(t1, y3);
	gf254_add(x3, x3, t1);
	gf254_add(x3, x3, t2);
}

INLINE void ec_dadd_weis_proj(gf_t x3, gf_t z3, gf_t x1, gf_t z1, gf_t x2,
        gf_t z2, gf_t x0, gf_t z0, gf_t a1) {
    gf_t t1, t2, t3, t4;
    gf254_add(t1, x1, z1);
    gf254_add(t2, x2, z2);

    gd_t u0, u1;
    gf254_muln(u0, x1, x2);
    gf254_muln(u1, z1, z2);
    gf254_add2(u0, u0, u1);
    gf254_rdcn(t3, u0);

    gf254_mul(t2, t2, t1);
    gf254_add(t2, t2, t3);

    gf254_sqr(t2, t2);
    gf254_sqr(t3, t3);

    gf254_mul(t4, x1, z1);
    gf254_mul(z3, t2, x0);
    gf254_mul(x3, t3, z0);
}

void ak_weis(gf_t x2, gf_t x1, gf_t y1, gf_t a1, uint64_t k[]) {
	int S[256], i, j, l, s1, s2, s;
	uint64_t k1[2], k2[2];
	gf_t tabx[3], taby[3], tabz[3], one;
	uint64_t _a[] = { 0, 0, 1, 0 };
	uint64_t _one[] = { 1, 0, 0, 0 };
	__m128i t0, t1, t2, t3;

    gls_recoding(k, k1, k2, &s1, &s2, BETA_E1);
	ak_recoding(S, &i, k1, k2, s1, s2);
	//     for (j = 0; j < i; j++) {
	//         printf("%d ", S[j]);
	//     }
	//     printf("%d\n", i);

	gf254_cpy(tabx[0], x1);
	gf254_cpy(taby[0], y1);

	//gf254_out(x1);
	//gf254_out(y1);

	/* PSI */
	t0 = _mm_unpacklo_epi64(x1[0], x1[1]);
	t1 = _mm_unpackhi_epi64(x1[0], x1[1]);
	t2 = _mm_unpacklo_epi64(y1[0], y1[1]);
	t3 = _mm_unpackhi_epi64(y1[0], y1[1]);

	tabx[1][1] = t1;
	tabx[1][0] = _mm_xor_si128(t0, t1);
	taby[1][1] = _mm_xor_si128(tabx[1][0], t3);
	taby[1][0] = _mm_xor_si128(t0, t3);
	taby[1][0] = _mm_xor_si128(taby[1][0], t2);

	t0 = tabx[1][0];
	t1 = tabx[1][1];
	t2 = taby[1][0];
	t3 = taby[1][1];
	tabx[1][0] = _mm_unpacklo_epi64(t0, t1);
	tabx[1][1] = _mm_unpackhi_epi64(t0, t1);
	taby[1][0] = _mm_unpacklo_epi64(t2, t3);
	taby[1][1] = _mm_unpackhi_epi64(t2, t3);

	if (s1 == 1) {
		gf254_add(taby[0], taby[0], tabx[0]);
	}

	if (s2 == 1) {
		gf254_add(taby[1], taby[1], tabx[1]);
	}

	gf254_inp(tabz[0], _one);
	gf254_inp(tabz[1], _one);

	gf254_inp(one, _a);
	gf254_add(taby[1], taby[1], tabx[1]);
	//ec_add_aff(tabx[2], taby[2], tabx[0], taby[0], tabx[1], taby[1], one);
	ec_add_mix_LD(tabx[2], taby[2], tabz[2], tabx[0], taby[0], tabz[0], tabx[1],
			taby[1]);
	gf254_add(taby[1], taby[1], tabx[1]);

	gf254_sqr(taby[0], a1);
	gf254_mul(tabx[0], tabx[0], taby[0]);
	gf254_mul(tabx[1], tabx[1], taby[0]);
	gf254_mul(tabx[2], tabx[2], taby[0]);

	int b0 = 0, b1 = 0, b2 = 0, b3 = 0, _b1 = 0, _b2 = 0, _b3 = 0, c1 = 1, c2 =
			0, c3 = 0, _c2 = 0;

	for (j = 0; j < i; j++) {
		b1 = (S[j] == R2);
		b2 = (S[j] == _R1);
		b3 = (S[j] == _R2);

		c1 = _b1 ^ b1;
		c3 = _b3 ^ b3;
		_c2 = _b2 | _b3;
		c2 = (_b2 | _b3) ^ (b2 | b3);
		gf254_cswap(tabx[1], tabx[2], (c1 & ~_c2) | (c3 & _c2));
		gf254_cswap(tabz[1], tabz[2], (c1 & ~_c2) | (c3 & _c2));
		gf254_cswap(tabx[0], tabx[2], (c1 & _c2) | (c3 & ~_c2));
		gf254_cswap(tabz[0], tabz[2], (c1 & _c2) | (c3 & ~_c2));
		gf254_cswap(tabx[0], tabx[1], c2);
		gf254_cswap(tabz[0], tabz[1], c2);
		ladd_weis_proj(tabx[0], tabz[0], tabx[1], tabz[1], tabx[2], tabz[2],
				a1);
		/*if (S[j] == R1)
		 * //X0, Z0, X1, Z1 := diffAdd(X0, Z0, X1, Z1, X2, Z2, a1);
		 * if (S[j] == R2)
		 * //X0, Z0, X2, Z2 := diffAdd(X0, Z0, X2, Z2, X1, Z1, a1);
		 * if (S[j] == _R1)
		 * //X1, Z1, X0, Z0 := diffAdd(X1, Z1, X0, Z0, X2, Z2, a1);
		 * if (S[j] == _R2);
		 * //X1, Z1, X2, Z2 := diffAdd(X1, Z1, X2, Z2, X0, Z0, a1); */
		_b1 = b1;
		_b2 = b2;
		_b3 = b3;
	}

	gf254_cswap(tabx[0], tabx[1], _b2 | _b3);
	gf254_cswap(tabz[0], tabz[1], _b2 | _b3);
	gf254_cswap(tabx[0], tabx[2], _b3);
	gf254_cswap(tabz[0], tabz[2], _b3);
	gf254_cswap(tabx[1], tabx[2], _b1);
	gf254_cswap(tabz[1], tabz[2], _b1);

	ec_dadd_weis_proj(tabx[0], tabz[0], tabx[0], tabz[0], tabx[1], tabz[1], tabx[2],
			tabz[2], a1);

	//ec_add_aff(tabx[0], taby[0], tabx[0], taby[0], tabx[1], taby[1], c);
	gf254_inv(taby[2], tabz[0]);
	gf254_mul(x2, tabx[0], taby[2]);
	gf254_rdc(x2);
	//gf254_mul(x2, x2, a1);
	return;
}
