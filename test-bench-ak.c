#include <stdio.h>
#include <stdint.h>

#include "test.h"
#include "bench.h"

#include "gf254.h"
#include "gf254.c"

#define INT(A, B, C, D) B, A, D, C
#define REV(A, B, C, D) D, C, B, A

int ak(void) {
	int code = 0;
	gf_t u, v, xp, yp, b;
	uint64_t kw[] = {REV(0x0C4F2FBFB3CFC3A4, 0xFAD834D039EF44AA, 0xD30D848B1DB31FF1, 0x423A54266EE690C0)};
	uint64_t _xp[] = {INT(0x4880D3D964EFB8EB, 0xBCBF6A65644CF23B, 0x46DF3DC5688CB445, 0x3F5372D3598E7579)};
	uint64_t _yp[] = {INT(0x505FD79B9FD1929D, 0x31C43E6EF94FF925, 0x39AAC98051735021, 0x83244DEA20EBCE43)};
	//uint64_t _u[] = {INT(0x63760666CCF5316D, 0xDEF120C9A0618559, 0x2F6A4EA81FD3710F, 0x0D30DF8384BA2173)};
	uint64_t _u[] = {INT(0x10A821DC6FEFB251, 0xD551FBDD5C6F9A78, 0x68E61727F919F2E8, 0x1BBD7628274C0E53)};
	uint64_t _b[] = { 0x2000200000, 0, 0, 0 };

	gf254_inp(xp, _xp);
	gf254_inp(yp, _yp);
	gf254_inp(u, _u);
	gf254_inp(b, _b);

	TEST_ONCE("2-dim ak laddering on weis model is correct") {
		ak_weis(v, xp, yp, b, kw);
		TEST_ASSERT(gf254_cmp(v, u) == 0, end);
	}
	TEST_END;
        
	code = 1;
end:
	return code;
}

void bench(void) {
	gf_t u, v, xp, yp, b;
	uint64_t k[4];

    BENCH_BEGIN("ak_recoding") {
        int S[256], s1, s2, l;
        uint64_t k1[2], k2[2];
        gf254_rnd((__m128i*)k);
        gls_recoding(k, k1, k2, &s1, &s2, BETA_E1);
        ak_recoding(S, &l, k1, k2, s1, s2);
        BENCH_ADD(ak_recoding(S, &l, k1, k2, s1, s2));
    } BENCH_END;


	BENCH_BEGIN("ak_weis") {
		gf254_rnd((__m128i*)k);
		ak_weis(v, xp, yp, b, k);
		BENCH_ADD(ak_weis(v, xp, yp, b, k));
	}
	BENCH_END;
}

int main(int argc, char *argv[]) {

	printf("\n-- Running tests for montgomery laddering:\n\n");

	ak();

	printf("\n-- Running benchmarks:\n\n");

	bench();

	return 0;
}
