//
// Created by azhukova on 2/27/18.
//

#include <stdio.h>
#include <stdarg.h>
#include "pastml.h"

int QUIET = FALSE;

void log_info(const char* message, ...) {
    if (QUIET == FALSE) {
        va_list args;
        va_start(args, message);
        vprintf(message, args);
        va_end(args);
    }
}
