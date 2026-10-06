#include "parallel.h"

namespace pyscal {

namespace {

int hardware_threads() {
    const unsigned n = std::thread::hardware_concurrency();
    return n > 0 ? static_cast<int>(n) : 1;
}

std::atomic<int> thread_count(hardware_threads());

}  // namespace

int num_threads() { return thread_count.load(); }

void set_num_threads(int n) { thread_count.store(n < 1 ? hardware_threads() : n); }

}  // namespace pyscal
