#include <stdio.h>
#include <stdint.h>

typedef unsigned int uint128_t __attribute__((mode(TI)));

/* schoolbook multiplication (4 x 2) 64-bit words
   (optimized for h&a scalar representation) */
#define SCHBOOK_4x2_HA(h, c, a, b)\
    h = ((uint128_t) a[0]*b[0]);\
    c[0] = h; c[1] = h >> 64;\
    h = ((uint128_t) a[0]*b[1]);\
    c[2] = h; c[3] = h >> 64;\
    h = ((uint128_t) a[0]*b[2]);\
    c[4] = h; c[5] = h >> 64;\
    c[6] = a[0]*b[3];\
    h = ((uint128_t) a[1]*b[0]);\
    c[8] = h; c[9] = h >> 64;\
    h = ((uint128_t) a[1]*b[1]);\
    c[10] = h; c[11] = h >> 64;\
    c[12] = a[1]*b[2];

/* schoolbook addition (4 x 2) 64-bit words
   result on MSW c[5] | c[3] | c[1] | c[0] LSW
   (optimized for h&a scalar representation) */
#define SCHBOOK_SUM_4x2_HA(c)\
    asm ("addq %3, %0 \n\t"\
         "adcq %4, %1 \n\t"\
         "adcq %5, %2 \n\t"\
         "addq %6, %0 \n\t"\
         "adcq %7, %1 \n\t"\
         "adcq %9, %2 \n\t"\
         "addq %8, %1 \n\t"\
         "adcq %10, %2 \n\t"\
    : "+r" (c[1]), "+r" (c[3]), "+r" (c[5])\
    : "r" (c[2]), "r" (c[4]), "r" (c[6]), "r" (c[8]), "r" (c[9]), "r" (c[10]), "r" (c[11]), "r" (c[12])\
    );

/* schoolbook multiplication (4 x 1) 64-bit words */
#define SCHBOOK_4x1(h, c, a, b)\
    h = ((uint128_t) a*b[0]);\
    c[0] = h; c[1] = h >> 64;\
    h = ((uint128_t) a*b[1]);\
    c[2] = h; c[3] = h >> 64;\
    h = ((uint128_t) a*b[2]);\
    c[4] = h; c[5] = h >> 64;\
    h = ((uint128_t) a*b[3]);\
    c[6] = h; c[7] = h >> 64;

/* schoolbook addition (4 x 1) 64-bit words
   result on MSW c[7] | c[5] | c[3] | c[1] | c[0] LSW*/
#define SCHBOOK_SUM_4x1(c)\
    asm ("addq %4, %0 \n\t"\
         "adcq %5, %1 \n\t"\
         "adcq %6, %2 \n\t"\
         "adcq $0, %3 \n\t"\
    : "+r" (c[1]), "+r" (c[3]), "+r" (c[5]), "+r" (c[7])\
    : "r" (c[2]), "r" (c[4]), "r" (c[6])\
    );

/* 128-bit addition with carry */
#define SUM_128(c0, c1, a0, a1)\
    asm ("addq %2, %0 \n\t"\
         "adcq %3, %1"\
    : "+r" (c0), "+r" (c1)\
    : "r" (a0), "r" (a1) : "cc"\
    );

/* 128-bit subtraction with carry */
#define SUB_128(c0, c1, a0, a1)\
    asm ("subq %2, %0 \n\t"\
         "sbbq %3, %1"\
    : "+r" (c0), "+r" (c1)\
    : "r" (a0), "r" (a1) : "cc"\
    );

/* 192-bit addition with carry */
#define SUM_192(c0, c1, c2, a0, a1, a2)\
    asm ("addq %3, %0 \n\t"\
         "adcq %4, %1 \n\t"\
         "adcq %5, %2"\
    : "+r" (c0), "+r" (c1), "+r" (c2)\
    : "m" (a0), "m" (a1), "m" (a2) : "cc"\
    );

/* 192-bit subtraction with carry */
#define SUB_192(c0, c1, c2, a0, a1, a2)\
    asm ("subq %3, %0 \n\t"\
         "sbbq %4, %1 \n\t"\
         "sbbq %5, %2"\
    : "+r" (c0), "+r" (c1), "+r" (c2)\
    : "m" (a0), "m" (a1), "m" (a2) : "cc"\
    );

/* 256-bit addition with carry */
#define SUM_256(c0, c1, c2, c3, a0, a1, a2, a3)\
    asm ("addq %4, %0 \n\t"\
         "adcq %5, %1 \n\t"\
         "adcq %6, %2 \n\t"\
         "adcq %7, %3"\
    : "+r" (c0), "+r" (c1), "+r" (c2), "+r" (c3)\
    : "m" (a0), "m" (a1), "m" (a2), "m" (a3) : "cc"\
    );

/* 256-bit subtraction with carry */
#define SUB_256(c0, c1, c2, c3, a0, a1, a2, a3)\
    asm ("subq %4, %0 \n\t"\
         "sbbq %5, %1 \n\t"\
         "sbbq %6, %2 \n\t"\
         "sbbq %7, %3"\
    : "+r" (c0), "+r" (c1), "+r" (c2), "+r" (c3)\
    : "m" (a0), "m" (a1), "m" (a2), "m" (a3) : "cc"\
    );

/* PROTECTED DIRECT RECODING (k -> k1, k2)
   Method described in http://cacr.uwaterloo.ca/techreports/2012/cacr2012-24.pdf (Sec. 3.2) */
void gls_recoding(uint64_t k[], uint64_t k1[], uint64_t k2[], int *k1neg, int *k2neg, uint64_t beta) {
    const uint64_t BETA_22 = beta; /* "t" term of #E = t^2 - (q-1)^2 */
    const uint64_t ALL_ZERO = 0;

    uint128_t reg_128; /* 128-bit "register" */

    uint64_t tmp[8], sign;
    uint64_t result_4x1[8];
    uint64_t b1[2], b1_times_t[3], b2, b2_times_t[2];

    /* b1 (-k div 2^127) */
    b1[1] = (k[3] << 1) | (k[2] >> 63);
    b1[0] = (k[2] << 1) | (k[1] >> 63);

    /* b2 (k*BETA_22 div 2^254) */
    SCHBOOK_4x1(reg_128, result_4x1, BETA_22, k);
    SCHBOOK_SUM_4x1(result_4x1);
    b2 = (result_4x1[5] >> 62) | (result_4x1[7] << 2);

    //round
    b1[0] = b1[0] + ((k[1] >> 62) & 0x1);
    b2 = b2 + ((result_4x1[5] >> 61) & 0x1);

    /* b1*t */
    reg_128 = ((uint128_t) BETA_22*b1[0]);
    b1_times_t[0] = reg_128; b1_times_t[1] = reg_128 >> 64;
    reg_128 = ((uint128_t) BETA_22*b1[1]);
    b1_times_t[2] = reg_128 >> 64;
    SUM_128(b1_times_t[1], b1_times_t[2], (uint64_t) reg_128, ALL_ZERO);

    /* b2*t */
    reg_128 = ((uint128_t) BETA_22*b2);
    b2_times_t[0] = reg_128; b2_times_t[1] = reg_128 >> 64;

    /** k1 computation */

    /* b1 */
    tmp[0] = b1[0];
    tmp[1] = b1[1];
    tmp[2] = 0;
    tmp[3] = 0;

    /* b1 + k */
    SUM_256(tmp[0], tmp[1], tmp[2], tmp[3], k[0], k[1], k[2], k[3]);

    /* b1*q (q = 2^127) */
    tmp[4] = 0;
    tmp[5] = b1[0] << 63;
    tmp[6] = b1[0] >> 1 | b1[1] << 63;
    tmp[7] = b1[1] >> 1;

    /* b1*q + b2*t */
    SUM_256(tmp[4], tmp[5], tmp[6], tmp[7], b2_times_t[0], b2_times_t[1], ALL_ZERO, ALL_ZERO);

    /* k1 sign (0 for positive, 1 for negative) */
    sign = (tmp[6] > tmp[2]) | ((tmp[6] == tmp[2]) & (tmp[5] > tmp[1]));

    /* final subtraction */
    SUB_256(tmp[0], tmp[1], tmp[2], tmp[3], tmp[4], tmp[5], tmp[6], tmp[7]);

    /* Two's complement (if necessary) */
    tmp[0] = tmp[0] ^ (ALL_ZERO - sign);
    tmp[1] = tmp[1] ^ (ALL_ZERO - sign);
    SUM_128(tmp[0], tmp[1], sign, ALL_ZERO);

    /* output */
    *k1neg = (int) sign;
    k1[0] = tmp[0];
    k1[1] = tmp[1];

    /** k2 computation */

    /* b1t + b2 */
    tmp[0] = b1_times_t[0];
    tmp[1] = b1_times_t[1];
    tmp[2] = b1_times_t[2];

    SUM_192(tmp[0], tmp[1], tmp[2], b2, ALL_ZERO, ALL_ZERO);

    /* b2*q (q = 2^127) */
    tmp[3] = 0;
    tmp[4] = b2 << 63;
    tmp[5] = b2 >> 1;

    /* k2 sign (0 for positive, 1 for negative) */
    sign = (tmp[5] > tmp[2]) | ((tmp[5] == tmp[2]) & (tmp[4] > tmp[1]));

    /* final subtraction */
    SUB_192(tmp[0], tmp[1], tmp[2], tmp[3], tmp[4], tmp[5]);

    /* Two's complement (if necessary) */
    tmp[0] = tmp[0] ^ (ALL_ZERO - sign);
    tmp[1] = tmp[1] ^ (ALL_ZERO - sign);
    SUM_128(tmp[0], tmp[1], sign, ALL_ZERO);

    /* output */
    *k2neg = (int) sign;
    k2[0] = tmp[0];
    k2[1] = tmp[1];
}
