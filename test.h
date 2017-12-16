#include <stdio.h>
#include <stdint.h>
#include <string.h>
#include <assert.h>

#ifndef TEST_H
#define TEST_H

/** Maximum macro. */
#define min(x, y) ({                \
    typeof(x) _min1 = (x);          \
    typeof(y) _min2 = (y);          \
    (void) (&_min1 == &_min2);      \
    _min1 < _min2 ? _min1 : _min2; })

/** Red color. */
#define FAIL_COLOR              31

/** Green color. */
#define PASS_COLOR              32

/** Other colors. */
#define CMD_SET                 27
#define CMD_RESET               0
#define CMD_ATTR                1

/** Run a single test. */
#define TEST_ONCE(P)                                                    \
    printf("Testing if " P "...%*c", (64 - strlen(P)), ' ');            \

/** Run several tests. */
#define TEST_BEGIN(P)                                                   \
    printf("Testing if " P "...%*c", (64 - strlen(P)), ' ');            \
    for (int i = 0; i < TESTS; i++)                                     \

/** Check test result. */
#define TEST_ASSERT(C, LABEL)                                           \
    if (!(C)) {                                                         \
        test_fail();                                                    \
        printf("(at ");                                                 \
        printf(__FILE__);                                               \
        printf(":%d)\n", __LINE__);                                     \
        goto LABEL;                                                     \
    }

/** Compare two strings to check test result. */
#define TEST_EQUAL(A, B, LABEL)                                         \
    TEST_ASSERT(strncmp((const char *)A, (const char *)B,               \
        MAX(sizeof(A), sizeof(B)) == 0, LABEL)                          \

/** Finish test successfully. */
#define TEST_END                                                        \
    test_pass()                                                         \

/** Declare a test vector. */
#define TEST_VECTOR(VAR, STR)                                           \
    uint8_t VAR[(strlen(STR)) >> 1];                                    \
    assert(test_conv(VAR, sizeof(VAR), STR) == (strlen(STR) >> 1));     \

/** Print something in a test. */
#define TEST_PRINT(VAR)                                                 \
    test_print(#VAR, VAR, sizeof(VAR))                                  \

/** Number of tests to run. */
#define TESTS 100

/**
 * Prints a string indicating that the test failed.
 */
void test_fail(void);

/**
 * Prints a string indicating that the test passed.
 */
void test_pass(void);

/**
 * Convert a hex string to a byte buffer.
 */
size_t test_conv(uint8_t *buf, size_t len, const char *str);

/**
 * Prints a byte buffer.
 */
void test_print(const char *label, const uint8_t *buf, size_t len);

#endif /* TEST_H */
