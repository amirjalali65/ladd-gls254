#include <stdio.h>

#include "gf254.h"
#define INLINE static inline
#include "gf254.c"
#include "ladd.c"

void monty_weis(gf_t v, gf_t w, gf_t x, gf_t y, gf_t b, uint64_t k[], int recover) {
	uint64_t p, c, bit, one[] = { 1, 0, 0, 0 };
	gf_t _x1, _z1, _x0, _z0, t0, t1;
	int i;

	gf254_cpy(_x0, x);
	gf254_inp(_z0, one);

	gf254_sqr(_z1, x);
	gf254_add(_x1, _z1, b);
	gf254_sqr(_x1, _x1);

	p = 0;
	for (i = 251; i >= 0; i--) {
		bit = (k[i / 64] >> (i % 64)) & 1;
		c = bit ^ p;
		p = bit;
		gf254_cswap(_x0, _x1, c);
		gf254_cswap(_z0, _z1, c);
		ladd_weis(_x0, _z0, _x1, _z1, x, b);
	}
	gf254_cswap(_x0, _x1, p);
	gf254_cswap(_z0, _z1, p);

	if (recover) {
		gf254_mul(t0, _z0, _z1);
		gf254_mul(_z0, _z0, x);
		gf254_add(_z0, _z0, _x0);
		gf254_mul(_z1, _z1, x);
		gf254_mul(_x0, _x0, _z1);
		gf254_add(_z1, _z1, _x1);
		gf254_mul(_z1, _z1, _z0);

		gf254_sqr(t1, x);
		gf254_add(t1, t1, y);
		gf254_mul(t1, t1, t0);
		gf254_add(t1, t1, _z1);

		gf254_mul(t0, t0, x);
		gf254_inv(t0, t0);
		gf254_mul(t1, t1, t0);
		gf254_mul(_x1, _x0, t0);
		gf254_add(_z1, _x1, x);
		gf254_mul(_z1, _z1, t1);
		gf254_add(_z1, _z1, y);

		/* init with zero and sort the special cases */
		one[0] ^= recover;
		gf254_inp(_z0, one);
		gf254_sel(v, _x1, x, gf254_cmp(_z1, _z0));
		gf254_add(w, x, y);
		gf254_sel(w, _z1, w, gf254_cmp(_z1, _z0));
		gf254_rdc(v);
		gf254_rdc(w);
	} else {
		gf254_inv(_z0, _z0);
		gf254_mul(v, _x0, _z0);
		gf254_rdc(v);
	}
}

void monty_huff(gf_t v, gf_t w, gf_t gamma, uint64_t k[]) {
	uint64_t p, c, bit, one[] = { 1, 0, 0, 0 };
	int i;
	gf_t _w1, _z1, _w0, _z0, _w;

	gf254_inp(_z0, one);
	gf254_cpy(_w0, w);

	gf254_sqr(_w1, w);
	gf254_mul1(_z1, w, gamma);
	gf254_add(_z1, _z0, _z1);
	gf254_sqr(_z1, _z1);
	gf254_sqr(_z1, _z1);

	gf254_inv(_w, w);

	p = 0;
	for (i = 251; i >= 0; i--) {
		bit = (k[i / 64] >> (i % 64)) & 1;
		c = bit ^ p;
		p = bit;
		gf254_cswap(_w0, _w1, c);
		gf254_cswap(_z0, _z1, c);
		ladd_huff(_w0, _z0, _w1, _z1, _w, gamma);
	}
	gf254_cswap(_w0, _w1, p);
	gf254_cswap(_z0, _z1, p);

	gf254_inv(_z0, _z0);
	gf254_mul(v, _w0, _z0);
	gf254_rdc(v);
}

void monty_edws(gf_t v, gf_t w, gf_t c1, uint64_t k[]) {
	uint64_t p, c, bit, one[] = { 1, 0, 0, 0 };
	int i;
	gf_t _w1, _z1, _w0, _z0, t0;

	gf254_cpy(_w0, w);
	gf254_sqr(_w1, _w0);
	gf254_lsh(_z1, w);
	gf254_mul1(_z1, _z1, c1);
	gf254_sqr(_z1, _z1);
	gf254_sqr(_z1, _z1);
	gf254_inp(_z0, one);
	gf254_add(_z1, _z1, _z0);

	gf254_inv(t0, w);

	p = 0;
	for (i = 251; i >= 0; i--) {
		bit = (k[i / 64] >> (i % 64)) & 1;
		c = bit ^ p;
		p = bit;
		gf254_cswap(_w0, _w1, c);
		gf254_cswap(_z0, _z1, c);
		ladd_edws(_w0, _z0, _w1, _z1, t0, c1);
	}
	gf254_cswap(_w0, _w1, p);
	gf254_cswap(_z0, _z1, p);

	gf254_inv(_z0, _z0);
	gf254_mul(v, _w0, _z0);
	gf254_rdc(v);
}
