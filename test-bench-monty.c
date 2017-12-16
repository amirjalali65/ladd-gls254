#include <stdio.h>
#include <stdint.h>

#include "test.h"
#include "bench.h"

#include "gf254.h"
#include "gf254.c"

#define INT(A, B, C, D) B, A, D, C
#define REV(A, B, C, D) D, C, B, A

void monty_weis(gf_t v, gf_t w, gf_t x, gf_t y, gf_t b, uint64_t k[], int recover);
void monty_huff(gf_t v, gf_t w, gf_t gamma, uint64_t k[]);
void monty_edws(gf_t v, gf_t w, gf_t c1, uint64_t k[]);

int monty(void) {
	int code = 0;
	gf_t x0, y0, w0, z0, x1, y1, w1, z1, x, w, b, gamma, c;
	uint64_t _one[] = { 1, 0, 0, 0 };

	uint64_t _x0[] = {INT(0x1648F216FF30D3A8, 0x99A10A375488D5C0, 0x4A24A216ECC58834, 0xFE641E152A850623)};
	uint64_t _y0[] = {INT(0x530FDEC4CCA02FA3, 0x8E2C4C3F0986153F, 0x11D28B53A47E3B5E, 0x0EB2B49C101195EA)};
	uint64_t _x1[] = {INT(0, 0, 0, 0)};
	uint64_t _y1[] = {INT(0, 0, 0, 0)};
	uint64_t kw[] = {REV(0x1FFFFFFFFFFFFFFF, 0xFFFFFFFFFFFFFFFF, 0xA6B89E49D3FECD82, 0x8CA8D66BF4B88ED4)};
	uint64_t _b[] = { 0xE2DA921E91E38DD1, 0, 0, 0 };

	gf254_inp(x0, _x0);
	gf254_inp(y0, _y0);
	gf254_inp(x1, _x1);
	gf254_inp(y1, _y1);
	gf254_inp(z0, _one);
	gf254_inp(z1, _one);
	gf254_inp(b, _b);

	TEST_ONCE("montgomery laddering on weis model is correct") {
		monty_weis(x1, y1, x0, y0, b, kw, 1);
		TEST_ASSERT(gf254_cmp(x0, x1) == 0, end);
		gf254_add(y0, y0, x0);
		TEST_ASSERT(gf254_cmp(y0, y1) == 0, end);
		gf254_add(y0, y0, x0);
	}
	TEST_END;

	uint64_t _w0[] = {INT(0x6AF8AEC822C08DAE, 0x1E0C4242B6647087, 0x28113AB224112618, 0x71DA56AE7D49784C)};
	uint64_t _w1[] = {INT(0, 0, 0, 0)};
	uint64_t kh[] = {REV(0X1FFFFFFFFFFFFFFF, 0xFFFFFFFFFFFFFFFF, 0x83260E3B053F9BDA, 0xB3CDFF18DB6C3040)};
	uint64_t _gamma[] = { 0x20000000200000, 0, 0, 0 };

	gf254_inp(w0, _w0);
	gf254_inp(w1, _w1);
	gf254_inp(z0, _one);
	gf254_inp(z1, _one);
	gf254_inp(gamma, _gamma);
	gf254_inp(w, _w0);

	TEST_ONCE("montgomery laddering on huff model is correct") {
		monty_huff(w1, w0, gamma, kh);
		TEST_ASSERT(gf254_cmp(w, w1) == 0, end);
	}
	TEST_END;

	uint64_t _o0[] = {INT(0x4CCAE2F3380BC7C4, 0x8412E5A4106B956B, 0x7FD33B6D6AB4B03C, 0x751F42496229BA72)};
	uint64_t _o1[] = {INT(0, 0, 0, 0)};
	uint64_t _o[] = {INT(0x4CCAE2F3380BC7C4, 0x8412E5A4106B956B, 0x7FD33B6D6AB4B03C, 0x751F42496229BA72)};
	uint64_t _c[] = {0x4000400040004, 0, 0, 0};

	gf254_inp(w0, _o0);
	gf254_inp(w1, _o1);
	gf254_inp(z0, _one);
	gf254_inp(z1, _one);
	gf254_inp(c, _c);
	gf254_inp(w, _o);

	TEST_ONCE("montgomery laddering on edws model is correct") {
		monty_edws(w1, w0, c, kh);
		TEST_ASSERT(gf254_cmp(w, w1) == 0, end);
	}
	TEST_END;

	code = 1;
end:
	return code;
}

void bench(void) {
	gf_t x0, y0, w0, z0, x1, y1, w1, z1, x, w, b, gamma, c;
	uint64_t k[4];
	uint64_t _gamma[] = { 0x20000000200000, 0, 0, 0 };
	uint64_t _b[] = { 0xE2DA921E91E38DD1, 0, 0, 0 };
	uint64_t _c[] = {0x4000400040004, 0, 0, 0};

	gf254_inp(gamma, _gamma);
	gf254_inp(b, _b);
	gf254_inp(c, _c);

	BENCH_BEGIN("monty_weis") {
		gf254_rnd((__m128i*)k);
		gf254_rnd(x0);
		gf254_rnd(y0);
		gf254_rnd(x1);
		gf254_rnd(y1);
		monty_weis(x1, y1, x0, y0, b, k, 0);
		BENCH_ADD(monty_weis(x1, y1, x0, y0, b, k, 0));
	}
	BENCH_END;

	BENCH_BEGIN("monty_weis (recover y)") {
		gf254_rnd((__m128i*)k);
		gf254_rnd(x0);
		gf254_rnd(y0);
		gf254_rnd(x1);
		gf254_rnd(y1);
		monty_weis(x1, y1, x0, y0, b, k, 1);
		BENCH_ADD(monty_weis(x1, y1, x0, y0, b, k, 1));
	}
	BENCH_END;

	BENCH_BEGIN("monty_huff") {
		gf254_rnd((__m128i*)k);
		gf254_rnd(w0);
		gf254_rnd(w1);
		monty_huff(w1, w0, gamma, k);
		BENCH_ADD(monty_huff(w1, w0, gamma, k));
	}
	BENCH_END;

	BENCH_BEGIN("monty_edws") {
		gf254_rnd((__m128i*)k);
		gf254_rnd(w0);
		gf254_rnd(w1);
		monty_edws(w1, w0, c, k);
		BENCH_ADD(monty_edws(w1, w0, c, k));
	}
	BENCH_END;
}

int main(int argc, char *argv[]) {

	printf("\n-- Running tests for montgomery laddering:\n\n");

	monty();

	printf("\n-- Running benchmarks:\n\n");

	bench();

	return 0;
}
