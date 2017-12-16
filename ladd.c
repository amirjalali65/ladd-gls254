#include <stdio.h>

#include "gf254.h"

INLINE void ladd_weis(gf_t x0, gf_t z0, gf_t x1, gf_t z1, gf_t x, gf_t b) {
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

	gf254_sqr(z0, z0);
	gf254_sqr(t1, x0);
	gf254_mul1(x0, z0, b);
	gf254_add(x0, x0, t1);
	gf254_sqr(x0, x0);
	gf254_mul(z0, t1, z0);
}

INLINE void ladd_weis_proj(gf_t x1, gf_t z1, gf_t x2, gf_t z2, gf_t x0, gf_t z0, gf_t a1) {
	gf_t t1, t2, t3, t4;
	gf254_add(t1, x1, z1);
	gf254_add(t2, x2, z2);

	gd_t u0, u1, u2;
	gf254_muln(u0, x1, x2);
	gf254_muln(u1, z1, z2);
	gf254_add2(u0, u0, u1);
	gf254_rdcn(t3, u0);

	gf254_mul(t2, t2, t1);
	gf254_add(t2, t2, t3);

	gf254_sqr(t1, t1);
	gf254_sqr(t2, t2);
	gf254_sqr(t3, t3);

	gf254_mul(t4, x1, z1);
	gf254_mul(z2, t2, x0);
	gf254_mul(x2, t3, z0);

	gf254_mul1(z1, t4, a1);

	gf254_sqr(x1, t1);
	gf254_sqr(z1, z1);
}

/* Older formula. */
INLINE void ladd_huff_old(gf_t w0, gf_t z0, gf_t w1, gf_t z1, gf_t wp, gf_t gamma) {
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

	/* W3 := gamma * (W0*Z0)^2 */
	gf254_mul(t3, w0, z0);
	gf254_mul1(t3, t3, gamma);
	//gf254_mul_gamma(t3, t3);

	/* W4 := ((W0 + Z0)*(W1 + Z1) + T)^2 */
	gf254_add(t1, t1, t0);
	gf254_sqr(w1, t1);

	/* Z3 := (W0 + Z0)^4 */
	gf254_sqr(t2, t2);
	gf254_sqr(z0, t2);

	/* Z4 := wP * T^2 */
	gf254_sqr(t0, t0);
	gf254_mul(z1, wp, t0);

	gf254_sqr(w0, t3);
}

INLINE void ladd_huff(gf_t w0, gf_t z0, gf_t w1, gf_t z1, gf_t wp, gf_t gamma) {
	gf_t t0, t1, t2, t3;
	gd_t u0, u1;
	//A := W0*Z0;
	gf254_mul(t0, w0, z0);
	//B := W0*Z1;
	gf254_mul(t1, w0, z1);
	//C := W1*Z0;
	gf254_mul(t2, w1, z0);
	//Z3 := (W0/gamma + Z0)^4;
	gf254_mul1(t3, w0, gamma);
	gf254_add(z0, z0, t3);
	gf254_sqr(z0, z0);
	gf254_sqr(z0, z0);
	//W3 := A^2;
	gf254_sqr(w0, t0);

	//W4 := (B + C)^2;
	gf254_add(w1, t1, t2);
	gf254_sqr(w1, w1);
	//Z4 := B*C + wP*W4;
	gf254_muln(u0, t1, t2);
	gf254_muln(u1, wp, w1);
	gf254_add2(u0, u0, u1);
	gf254_rdcn(z1, u0);
}

INLINE void ladd_edws_old(gf_t w0, gf_t z0, gf_t w1, gf_t z1, gf_t w,
	gf_t c1, gf_t c2, gf_t c3) {
	gf_t t0, t1, t2, t3;
	gd_t u0, u1, u2;

	/* C := W0*(W0 + Z0); */
	gf254_add(t0, w0, z0);
	gf254_mul(t0, w0, t0);

	/* D := W1*(W1 + Z1); */
	gf254_add(t1, w1, z1);
	gf254_mul(t1, t1, w1);

	/* E := (W0 + W1)*(W0 + W1 + Z0 + Z1) + C + D; */
	gf254_add(t2, w0, w1);
	gf254_add(t3, z0, z1);
	gf254_add(t3, t3, t2);
	gf254_mul(t2, t2, t3);
	gf254_add(t2, t2, t0);
	gf254_add(t2, t2, t1);

	/* V := C * D; */
	gf254_muln(u0, t0, t1);

	/* W4 := U := d1 * E^2; */
	gf254_mul1(w1, t2, c3);
	gf254_mul_1u(w1, w1);
	gf254_sqr(w1, w1);

	/* Z4 := V + U/wP; */
	gf254_muln(u1, w1, w);
	gf254_add2(u0, u0, u1);
	gf254_rdcn(z1, u0);

	/* Z3 := W3 + (d1 * Z0^4 + (d2/d1 + 1)*W0^4); */
	gf254_mul1(t2, z0, c1);
	gf254_mul_0u(t2, t2);
	gf254_mul(t1, w0, c2);
	gf254_add(t1, t1, t2);
	gf254_sqr(t1, t1);
	gf254_sqr(t1, t1);

	/* W3 := C^2; */
	gf254_sqr(w0, t0);
	gf254_add(z0, t1, w0);
}

INLINE void ladd_edws(gf_t w0, gf_t z0, gf_t w1, gf_t z1, gf_t w,
		gf_t c) {
	gf_t t0, t1, t2, t3;
	gd_t u0, u1, u2;

	gf254_mul(t0, w0, z0);
	gf254_mul(t1, w0, z1);
	gf254_mul(t2, w1, z0);

	gf254_lsh(t3, w0);
	gf254_mul1(t3, t3, c);
	gf254_add(t3, t3, z0);
	gf254_sqr(z0, t3);
	gf254_sqr(z0, z0);
	gf254_sqr(w0, t0);

	gf254_add(t0, t1, t2);
	gf254_sqr(w1, t0);
	gf254_muln(u0, t1, t2);
	gf254_muln(u1, w, w1);
	gf254_add2(u0, u0, u1);
	gf254_rdcn(z1, u0);
}
