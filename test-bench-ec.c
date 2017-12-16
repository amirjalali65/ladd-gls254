#include <stdio.h>
#include <stdint.h>

#include "test.h"
#include "bench.h"

#include "gf254.h"
#include "gf254.c"

#define INT(A, B, C, D) B, A, D, C
#define REV(A, B, C, D) D, C, B, A

int djb(void) {
	int code = 0;
	gf_t u, v, xp, yp, b;
	uint64_t kw[] = {REV(0x1FFFFFFFFFFFFFFF, 0xFFFFFFFFFFFFFFFF, 0xA6B89E49D3FECD82, 0x8CA8D66BF4B88ED4)};
	uint64_t _xw[] = {INT(0x1648F216FF30D3A8, 0x99A10A375488D5C0, 0x4A24A216ECC58834, 0xFE641E152A850623)};
	uint64_t _yw[] = {INT(0x45472CD23390FC0B, 0x178D46085D0EC0FF, 0x5BF6294548BBB36A, 0xF0D6AA893A9493C9)};
	uint64_t _u[] ={INT(0x0615FFEB31812126, 0x0D128268AF597EB2, 0x62F8B8219572E884, 0xEF70ED6C7581312F)};
	uint64_t _b[] = { 0xE2DA921E91E38DD1, 0, 0, 0 };

	gf254_inp(xp, _xw);
	gf254_inp(yp, _yw);
	gf254_inp(u, _u);
	gf254_inp(b, _b);

	TEST_ONCE("2-dim djb laddering on weis model is correct") {
		ec_mul_djb_weis(v, xp, yp, b, kw);
		TEST_ASSERT(gf254_cmp(v, u) == 0, end);
	}
	TEST_END;

	uint64_t  kh[] = {REV(0X1FFFFFFFFFFFFFFF, 0xFFFFFFFFFFFFFFFF, 0x83260E3B053F9BDA, 0xB3CDFF18DB6C3040)};
	uint64_t _xp[] = {INT(0x0863B55F5D1D2607, 0xCFB7915C977B6807, 0x44D89176116A9522, 0x9E5EAC8B9778288D)};
	uint64_t _yp[] = {INT(0x7594A667E3C170CC, 0x0BE0BFD37BFFE8FA, 0x283D5DB73F3B858F, 0x14A00C2772907470)};
	uint64_t _w[]  = {INT(0x0783868507D6BDA2, 0x6FED7C627E0591D1, 0x62EC7879E2B00DF1, 0xA891949B2B4DB590)};
	uint64_t _wp[] = {INT(0x0CFA49A060A65B8E, 0x1B9A878F1F35AF63, 0x68A94B9222BBE292, 0x4CFF29D0DD7E1A83)};
	uint64_t _wq[] = {INT(0x64530232421DB91C, 0x5765AE5FC24BB5E0, 0x68A94B9222BBE292, 0x4CFF29D0DD7E1A83)};
	uint64_t _wpq[]= {INT(0x5AB3CCE342750B91, 0xE3E4D3B015662F28, 0x63EA3316F1519BA6, 0x7295CE7B49C37CD2)};
	uint64_t _wpmq[]={INT(0x3959FFF5B3249037, 0x91711DCB5CA553FA, 0x63EA3316F1519BA6, 0x7295CE7B49C37CD2)};
	//uint64_t _gamma[] = { 0x2000200000, 0, 0, 0 };
	uint64_t _gamma[] = { 0x20000000200000, 0, 0, 0 };
	uint64_t __b[] = {INT(0x400, 0x0000040000000000, 0, 0)};

	gf_t wp, wq, wpq, wpmq, gamma;
	gf254_inp(xp, _xp);
	gf254_inp(yp, _yp);
	gf254_inp(wp, _wp);
	gf254_inp(wq, _wq);
	gf254_inp(wpq, _wpq);
	gf254_inp(wpmq, _wpmq);
	gf254_inp(gamma, _gamma);
	gf254_inp(b, __b);
	gf254_inp(u, _w);

	TEST_ONCE("2-dim djb laddering on huff model is correct") {
		ec_mul_djb_huff(v, xp, yp, gamma, b, kh);
		TEST_ASSERT(gf254_cmp(v, u) == 0, end);
	}
	TEST_END;

	code = 1;
end:
	return code;
}

/*
int ak(void) {
	int code = 0;
	gf_t u, v, xp, yp, b;
	uint64_t kw[] = {REV(0x0C4F2FBFB3CFC3A4, 0xFAD834D039EF44AA, 0xD30D848B1DB31FF1, 0x423A54266EE690C0)};
	uint64_t _xp[] = {INT(0x4880D3D964EFB8EB, 0xBCBF6A65644CF23B, 0x46DF3DC5688CB445, 0x3F5372D3598E7579)};
	uint64_t _yp[] = {INT(0x505FD79B9FD1929D, 0x31C43E6EF94FF925, 0x39AAC98051735021, 0x83244DEA20EBCE43)};
	uint64_t _u[] ={INT(0x63760666CCF5316D, 0xDEF120C9A0618559, 0x2F6A4EA81FD3710F, 0x0D30DF8384BA2173)};
	uint64_t _b[] = { 0x2000200000, 0, 0, 0 };

	gf254_inp(xp, _xp);
	gf254_inp(yp, _yp);
	gf254_inp(u, _u);
	gf254_inp(b, _b);

	TEST_ONCE("2-dim ak laddering on weis model is correct") {
		ec_mul_ak_weis(v, xp, yp, b, kw);
		TEST_ASSERT(gf254_cmp(v, u) == 0, end);
	}
	TEST_END;
        
    BENCH_BEGIN("scmul_ak_recoding") {
        int S[256], s1, s2, j;
        uint64_t k1[2], k2[2];
        gf254_rnd((__m128i*)kw);
        gls_recoding(kw, k1, k2, &s1, &s2, BETA_E1);
        scmul_ak_recoding(S, &j, k1, k2, s1, s2);
        BENCH_ADD(scmul_ak_recoding(S, &j, k1, k2, s1, s2));
    } BENCH_END;


	BENCH_BEGIN("ec_mul_ak_weis") {
		gf254_rnd((__m128i*)kw);
		ec_mul_ak_weis(v, xp, yp, b, kw);
		BENCH_ADD(ec_mul_ak_weis(v, xp, yp, b, kw));
	}
	BENCH_END;

	code = 1;
end:
	return code;
}
*/
void bench(void) {
	gf_t x0, y0, w0, z0, x1, y1, w1, z1, x, w, b, gamma, c;
	uint64_t k[4];
	uint64_t _gamma[] = { 0x20000000200000, 0, 0, 0 };
	uint64_t _b[] = { 0xE2DA921E91E38DD1, 0, 0, 0 };
	uint64_t _c[] = {INT(0x400, 0x0000040000000000, 0, 0)};

	gf254_inp(gamma, _gamma);
	gf254_inp(b, _b);
	gf254_inp(c, _c);

	BENCH_BEGIN("ec_mul_djb_weis") {
		gf254_rnd((__m128i*)k);
		gf254_rnd(x0);
		gf254_rnd(y0);
		djb_weis(w0, x0, y0, b, k);
		BENCH_ADD(djb_weis(w0, x0, y0, b, k));
	}
	BENCH_END;

	BENCH_BEGIN("ec_mul_djb_huff") {
		gf254_rnd((__m128i*)k);
		gf254_rnd(x0);
		gf254_rnd(y0);
		djb_huff(v, xp, yp, gamma, b, k);
		BENCH_ADD(djb_huff(w0, x0, y0, gamma, c, k));
	}
	BENCH_END;
}

int main(int argc, char *argv[]) {

	bench();

	return 0;
}
