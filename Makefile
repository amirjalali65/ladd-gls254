# Pick compiler
ifeq ($(findstring icc, $(MAKECMDGOALS)), icc)
CC=icc
else
ifeq ($(findstring clang, $(MAKECMDGOALS)), clang)
CC=clang
else
CC=gcc
endif
endif

# Pick architecture
ifeq ($(findstring haswell, $(MAKECMDGOALS)), haswell)
ARCH=haswell
else
ifeq ($(findstring skylake, $(MAKECMDGOALS)), haswell)
ARCH=skylake
else
ARCH=native
endif
endif

CFLAGS=-O3 -funroll-loops -fomit-frame-pointer
COMMON=test.c bench.c util.c

gcc: all
icc: all
clang: all
haswell: all
skylake: all

folder:
	@if [ ! -x $(ARCH) ]; then mkdir $(ARCH); fi

all: folder test-bench-ffa bench-ladd test-bench-monty test-bench-djb test-bench-ak

test-bench-ffa:
	$(CC) $(CFLAGS) -march=$(ARCH) $(COMMON) test-bench-ffa.c  -o $(ARCH)/test-bench-ffa-$(CC)

bench-ladd:
	$(CC) $(CFLAGS) -march=$(ARCH) $(COMMON) bench-ladd.c  -o $(ARCH)/bench-ladd-$(CC)

test-bench-monty:
	$(CC) $(CFLAGS) -march=$(ARCH) $(COMMON) monty.c test-bench-monty.c  -o $(ARCH)/test-bench-monty-$(CC)

test-bench-djb:
	$(CC) $(CFLAGS) -march=$(ARCH) $(COMMON) djb.c gls.c test-bench-djb.c  -o $(ARCH)/test-bench-djb-$(CC)

test-bench-ak:
	$(CC) $(CFLAGS) -march=$(ARCH) $(COMMON) ak.c gls.c test-bench-ak.c  -o $(ARCH)/test-bench-ak-$(CC)

clean:
	rm -rf haswell skylake native
