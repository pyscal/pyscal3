/*
Threads for the per-atom loops of pyscal3 and of the vendored neighbour search.

parallel_for(n, body, grain) calls body(begin, end) on consecutive ranges of
[0, n), each at most grain long, spread over up to num_threads() threads. The
ranges are handed out one after the other, so the split does not depend on
the timing, and body must only write data that belongs to its own indices;
then the result is the same for any number of threads. Loops shorter than
two grains run on the calling thread. Threads are started for each call and
joined before it returns, so nothing is left running (safe with fork). The
first exception thrown by body is rethrown on the calling thread.
*/
#ifndef PYSCAL_PARALLEL_H
#define PYSCAL_PARALLEL_H

#include <algorithm>
#include <atomic>
#include <cstdint>
#include <exception>
#include <mutex>
#include <thread>
#include <vector>

namespace pyscal {

// number of threads used by parallel_for (at least 1)
int num_threads();
// n < 1 means all hardware threads
void set_num_threads(int n);

template <typename Body>
void parallel_for(std::int64_t n, Body &&body, std::int64_t grain = 256) {
    if (n <= 0) return;
    grain = std::max<std::int64_t>(grain, 1);
    const std::int64_t chunks = (n + grain - 1) / grain;
    const std::int64_t nthreads = std::min<std::int64_t>(num_threads(), chunks);
    if (nthreads <= 1 || chunks < 2) {
        body(std::int64_t(0), n);
        return;
    }
    std::atomic<std::int64_t> next(0);
    std::exception_ptr error;
    std::mutex error_lock;
    auto work = [&]() {
        try {
            for (;;) {
                const std::int64_t begin = next.fetch_add(grain);
                if (begin >= n) break;
                body(begin, std::min(begin + grain, n));
            }
        } catch (...) {
            std::lock_guard<std::mutex> lock(error_lock);
            if (!error) error = std::current_exception();
            next.store(n);
        }
    };
    std::vector<std::thread> threads;
    threads.reserve(nthreads - 1);
    for (std::int64_t t = 1; t < nthreads; t++) threads.emplace_back(work);
    work();
    for (auto &t : threads) t.join();
    if (error) std::rethrow_exception(error);
}

}  // namespace pyscal

#endif
