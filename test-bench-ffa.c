#include <stdio.h>
#include <stdint.h>

#include "test.h"
#include "bench.h"

#include "gf254.h"
#undef INLINE
#define INLINE /* */
#include "gf254.c"

int arith(void) {
	gf_t a, b, c, d, e;
	__m128i _a, _b, _c, _d;
	int code = 0;

	TEST_BEGIN("addition is commutative") {
		gf254_rnd(a);
		gf254_rnd(b);
		gf254_add(c, a, b);
		gf254_add(d, b, a);
		TEST_ASSERT(gf254_cmp(c, d) == 0, end);
	}
	TEST_END;

	TEST_BEGIN("addition is associative") {
		gf254_rnd(a);
		gf254_rnd(b);
		gf254_rnd(c);
		gf254_add(d, a, b);
		gf254_add(d, d, c);
		gf254_add(e, b, c);
		gf254_add(e, a, e);
		TEST_ASSERT(gf254_cmp(d, e) == 0, end);
	}
	TEST_END;

	TEST_BEGIN("multiplication is commutative") {
		gf254_rnd(a);
		gf254_rnd(b);
		gf254_mul(c, a, b);
		gf254_mul(d, b, a);
		TEST_ASSERT(gf254_cmp(c, d) == 0, end);
	}
	TEST_END;

	TEST_BEGIN("multiplication is associative") {
		gf254_rnd(a);
		gf254_rnd(b);
		gf254_rnd(c);
		gf254_mul(d, a, b);
		gf254_mul(d, d, c);
		gf254_mul(e, b, c);
		gf254_mul(e, a, e);
		TEST_ASSERT(gf254_cmp(d, e) == 0, end);
	}
	TEST_END;

	TEST_BEGIN("squaring is consistent") {
		gf254_rnd(a);
		gf254_sqr(b, a);
		gf254_mul(c, a, a);
		TEST_ASSERT(gf254_cmp(b, c) == 0, end);
	}
	TEST_END;

	TEST_BEGIN("inversion is consistent") {
		uint64_t _one[] = { 1, 0, 0, 0 };
		gf254_rnd(a);
		gf254_inv(b, a);
		gf254_mul(c, a, b);
		low_rdc(&c[0], &c[1]);
		gf254_inp(d, _one);
		TEST_ASSERT(gf254_cmp(c, d) == 0, end);
	}
	TEST_END;

	code = 1;
end:
	return code;
}

void bench(void) {
	gf_t a, b, c, d, e, f;

	BENCH_BEGIN("gf254_mul") {
		gf254_rnd(a);
		gf254_rnd(b);
		gf254_rnd(c);
		BENCH_ADD(gf254_mul(a, a, b));
	} BENCH_END;

	BENCH_BEGIN("gf254_mul1") {
		gf254_rnd(a);
		gf254_rnd(b);
		gf254_rnd(c);
		BENCH_ADD(gf254_mul1(c, a, b));
	} BENCH_END;

	BENCH_BEGIN("gf254_mul_gamma") {
		gf254_rnd(a);
		gf254_rnd(b);
		gf254_rnd(c);
		BENCH_ADD(gf254_mul_gamma(c, a));
	} BENCH_END;

	BENCH_BEGIN("gf254_sqr") {
		gf254_rnd(a);
		gf254_rnd(c);
		BENCH_ADD(gf254_sqr(c, a));
	} BENCH_END;

	BENCH_BEGIN("gf254_inv") {
		gf254_rnd(a);
		gf254_rnd(c);
		BENCH_ADD(gf254_inv(c, a));
	} BENCH_END;
        
	BENCH_BEGIN("gf254_cswap") {
		gf254_rnd(a);
		gf254_rnd(c);
		BENCH_ADD(gf254_cswap(c, a, 0));
	} BENCH_END;

	BENCH_BEGIN("gf254_sel") {
		gf254_rnd(a);
		gf254_rnd(b);
		gf254_rnd(c);
		BENCH_ADD(gf254_sel(c, a, b, 0));
	} BENCH_END;
}

int main(int argc, char *argv[]) {

	printf("\n-- Running tests for finite field arithmetic:\n\n");

	arith();

	printf("\n-- Running benchmarks:\n\n");

	bench();

	return 0;
}
