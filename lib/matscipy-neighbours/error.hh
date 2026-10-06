/*
 * matscipy-neighbours — Neighbour list for particle simulations
 * https://github.com/libAtoms/matscipy-neighbours
 *
 * SPDX-License-Identifier: MIT
 * Copyright (2014-2026) James Kermode, University of Warwick
 *                       Lars Pastewka, University of Freiburg
 *                       and others (see toplevel AUTHORS file)
 */

#ifndef MATSCIPY_ERROR_HH
#define MATSCIPY_ERROR_HH

#include <new>

#include "types.hh"

namespace matscipy {

/* Simple, Python-free error reporting. Core routines record a message here and
   return NL_ERROR (an internal/runtime failure), NL_INVALID_ARGUMENT (the
   caller passed something unusable) or NL_OUT_OF_MEMORY (an allocation, host
   or GPU, failed); the binding layer turns these into a Python RuntimeError /
   ValueError / MemoryError. The state is global (protected by the GIL
   when called from Python) and is cleared at the start of each top-level core
   call. */

constexpr int MAX_ERROR_STRING = 1024;

extern bool has_error;
extern char error_string[MAX_ERROR_STRING];

/* Record an error message. Returns NL_ERROR so callers can `return
   set_error(...)`. */
error_t set_error(const char *msg);

/* Record an invalid-argument message. Returns NL_INVALID_ARGUMENT. */
error_t set_invalid_argument(const char *msg);

/* printf-style variants. */
error_t set_errorf(const char *fmt, ...);
error_t set_invalid_argumentf(const char *fmt, ...);

/* Record an out-of-memory message. Returns NL_OUT_OF_MEMORY. */
error_t set_out_of_memory(const char *msg);

/* Reset error state. */
void clear_error();

/* Run `body` (returning error_t) and turn an allocation failure inside it,
   which the memory layer reports as std::bad_alloc, into NL_OUT_OF_MEMORY with
   `msg`. Wraps the public entry points, so the core API reports running out of
   memory as an error code rather than an exception or an abort. */
template <typename F>
error_t catch_out_of_memory(const char *msg, F &&body) {
    try {
        return body();
    } catch (const std::bad_alloc &) {
        return set_out_of_memory(msg);
    }
}

}  // namespace matscipy

#endif
