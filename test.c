#include <stdio.h>
#include <stdint.h>

#include "test.h"
#include "bench.h"

void test_fail(void) {
        printf("[%c[%d;%dm", CMD_SET, CMD_ATTR, FAIL_COLOR);
        printf("FAIL");
        printf("%c[%dm]\n", CMD_SET, CMD_RESET);
}

void test_pass(void) {
        printf("[%c[%d;%dm", CMD_SET, CMD_ATTR, PASS_COLOR);
        printf("PASS");
        printf("%c[%dm]\n", CMD_SET, CMD_RESET);
}
