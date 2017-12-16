#include <stdio.h>
#include <stdint.h>

#include "test.h"
#include "bench.h"

#include "gf254.h"
#define INLINE static inline
#include "gf254.c"
#undef INLINE
#define INLINE /* */
#include "ladd.c"

void bench(void) {
	gf_t a, b, c, d, e, f;

	BENCH_BEGIN("ladd_weis") {
		gf254_rnd(a);
		gf254_rnd(b);
		gf254_rnd(c);
		gf254_rnd(d);
		gf254_rnd(e);
		gf254_rnd(f);
		BENCH_ADD(ladd_weis(a, b, c, d, e, f));
	} BENCH_END;

    BENCH_BEGIN("ladd_weis_proj") {
		gf254_rnd(a);
		gf254_rnd(b);
		gf254_rnd(c);
		gf254_rnd(d);
		gf254_rnd(e);
		gf254_rnd(f);
		BENCH_ADD(ladd_weis_proj(a, b, c, d, e, f, a));
    } BENCH_END;

    BENCH_BEGIN("ladd_huff_old") {
		gf254_rnd(a);
		gf254_rnd(b);
		gf254_rnd(c);
		gf254_rnd(d);
		gf254_rnd(e);
		gf254_rnd(f);
		BENCH_ADD(ladd_huff_old(a, b, c, d, e, f));
    } BENCH_END;
             
    BENCH_BEGIN("ladd_huff") {
		gf254_rnd(a);
		gf254_rnd(b);
		gf254_rnd(c);
		gf254_rnd(d);
		gf254_rnd(e);
		gf254_rnd(f);
		BENCH_ADD(ladd_huff(a, b, c, d, e, f));
    } BENCH_END;

    BENCH_BEGIN("ladd_edws_old") {
		gf254_rnd(a);
		gf254_rnd(b);
		gf254_rnd(c);
		gf254_rnd(d);
		gf254_rnd(e);
		gf254_rnd(f);
		BENCH_ADD(ladd_edws_old(a, b, c, d, e, f, a, b));
    } BENCH_END;

    BENCH_BEGIN("ladd_edws") {
		gf254_rnd(a);
		gf254_rnd(b);
		gf254_rnd(c);
		gf254_rnd(d);
		gf254_rnd(e);
		gf254_rnd(f);
		BENCH_ADD(ladd_edws(a, b, c, d, e, f));
    } BENCH_END;         
}

int main(int argc, char *argv[]) {

	printf("\n-- Running benchmarks for laddering steps:\n\n");

	bench();

	return 0;
}
