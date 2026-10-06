/*
 * matscipy-neighbours — Neighbour list for particle simulations
 * https://github.com/libAtoms/matscipy-neighbours
 *
 * SPDX-License-Identifier: MIT
 * Copyright (2014-2026) James Kermode, University of Warwick
 *                       Lars Pastewka, University of Freiburg
 *                       and others (see toplevel AUTHORS file)
 */

#ifndef MATSCIPY_MEMORY_SPACE_HH
#define MATSCIPY_MEMORY_SPACE_HH

#include <cassert>
#include <cstddef>
#include <cstring>
#include <type_traits>

namespace matscipy {

/*
 * Memory-space abstraction. A buffer lives in exactly one space;
 * `Array<T, Space>` owns it and `deep_copy` moves bytes between spaces. The
 * same algorithm code runs on host (OpenMP) or device (CUDA/HIP) by choosing
 * the Space tag. Provides the operations the neighbour-list pipeline needs:
 * own, resize, deep-copy.
 *
 * The device-type codes match DLPack (kDLCPU=1, kDLCUDA=2, kDLROCm=10) so a
 * DLPack export is a direct tag read with no translation table.
 */
enum class DeviceType : int {
    CPU = 1,    /* kDLCPU  */
    CUDA = 2,   /* kDLCUDA */
    ROCm = 10,  /* kDLROCm */
};

struct HostSpace {
    static constexpr DeviceType device = DeviceType::CPU;
    static constexpr const char *name = "host";
};
struct CudaSpace {
    static constexpr DeviceType device = DeviceType::CUDA;
    static constexpr const char *name = "cuda";
};
struct HipSpace {
    static constexpr DeviceType device = DeviceType::ROCm;
    static constexpr const char *name = "hip";
};

/* The active GPU space. `Array<T, DeviceSpace>` is the device buffer type the
   backend code uses, so the same source compiles under CUDA or HIP. */
#if defined(MATSCIPY_ENABLE_CUDA)
using DeviceSpace = CudaSpace;
#elif defined(MATSCIPY_ENABLE_HIP)
using DeviceSpace = HipSpace;
#endif

namespace detail {

/* Raw byte-level hooks, one set per backend. Host is always linked; the device
   hooks live in a backend translation unit compiled only when that backend is
   enabled. They are non-template free functions, so Array<> stays header-only
   while the actual cudaMalloc/cudaMemcpy calls sit in the GPU backend source
   (memory_space_gpu.cc, compiled by nvcc/hipcc). */
void *alloc_host(std::size_t bytes);
void free_host(void *ptr);

#ifdef MATSCIPY_ENABLE_CUDA
void *alloc_cuda(std::size_t bytes);
void free_cuda(void *ptr);
/* Handles every direction (H2D/D2H/D2D/H2H) from the space tags. */
void copy_cuda(void *dst, DeviceType dst_dev, const void *src,
               DeviceType src_dev, std::size_t bytes);
#endif

}  // namespace detail

#if defined(MATSCIPY_ENABLE_CUDA) || defined(MATSCIPY_ENABLE_HIP)
/* Id of the currently-active GPU device (for the DLPack device tuple). */
int current_device_id();

/* GPU allocations go through a per-device caching pool (memory_space_gpu.cc):
   freed blocks stay reserved for reuse by later calls. Return them to the
   driver, e.g. before handing the memory to another library. */
void empty_gpu_cache();
#endif

namespace detail {

#ifdef MATSCIPY_ENABLE_HIP
void *alloc_hip(std::size_t bytes);
void free_hip(void *ptr);
void copy_hip(void *dst, DeviceType dst_dev, const void *src,
              DeviceType src_dev, std::size_t bytes);
#endif

template <typename Space>
inline void *space_alloc(std::size_t bytes) {
    if (bytes == 0) return nullptr;
    if constexpr (std::is_same_v<Space, HostSpace>) {
        return alloc_host(bytes);
    } else if constexpr (std::is_same_v<Space, CudaSpace>) {
#ifdef MATSCIPY_ENABLE_CUDA
        return alloc_cuda(bytes);
#else
        static_assert(!std::is_same_v<Space, CudaSpace>,
                      "Array<T, CudaSpace> requires the CUDA backend "
                      "(-DENABLE_CUDA=ON).");
        return nullptr;
#endif
    } else if constexpr (std::is_same_v<Space, HipSpace>) {
#ifdef MATSCIPY_ENABLE_HIP
        return alloc_hip(bytes);
#else
        static_assert(!std::is_same_v<Space, HipSpace>,
                      "Array<T, HipSpace> requires the HIP backend "
                      "(-DENABLE_HIP=ON).");
        return nullptr;
#endif
    } else {
        static_assert(std::is_same_v<Space, HostSpace>,
                      "Unknown memory space.");
        return nullptr;
    }
}

template <typename Space>
inline void space_free(void *ptr) {
    if (!ptr) return;
    if constexpr (std::is_same_v<Space, HostSpace>) {
        free_host(ptr);
    } else if constexpr (std::is_same_v<Space, CudaSpace>) {
#ifdef MATSCIPY_ENABLE_CUDA
        free_cuda(ptr);
#endif
    } else if constexpr (std::is_same_v<Space, HipSpace>) {
#ifdef MATSCIPY_ENABLE_HIP
        free_hip(ptr);
#endif
    }
}

}  // namespace detail

/*
 * Owning, move-only buffer of `T` in `Space`. A sized allocation, not a
 * std::vector: no per-element construction (T is a trivial numeric type) and no
 * copy ctor, so deep copies are explicit via deep_copy and a device round-trip
 * is never silent.
 *
 * resize() that grows DISCARDS existing contents; shrinking keeps the buffer
 * and just lowers the logical size. An allocation that fails for lack of
 * memory throws std::bad_alloc. data() on a device array is a device
 * pointer, only valid in kernels or as a deep_copy argument, never dereferenced
 * on the host.
 */
template <typename T, typename Space = HostSpace>
class Array {
    static_assert(std::is_trivially_copyable_v<T>,
                  "Array<T> requires a trivially copyable element type.");

    T *data_ = nullptr;
    std::size_t size_ = 0;
    std::size_t capacity_ = 0;

    void steal(Array &o) noexcept {
        data_ = o.data_;
        size_ = o.size_;
        capacity_ = o.capacity_;
        o.data_ = nullptr;
        o.size_ = 0;
        o.capacity_ = 0;
    }

   public:
    using value_type = T;
    using space_type = Space;
    static constexpr DeviceType device = Space::device;

    Array() = default;
    explicit Array(std::size_t n) { resize(n); }

    Array(const Array &) = delete;
    Array &operator=(const Array &) = delete;

    Array(Array &&o) noexcept { steal(o); }
    Array &operator=(Array &&o) noexcept {
        if (this != &o) {
            detail::space_free<Space>(data_);
            steal(o);
        }
        return *this;
    }

    ~Array() { detail::space_free<Space>(data_); }

    std::size_t size() const { return size_; }
    bool empty() const { return size_ == 0; }
    T *data() { return data_; }
    const T *data() const { return data_; }

    /* Grow (reallocating, discarding old contents) or shrink (in place). */
    void resize(std::size_t n) {
        if (n <= capacity_) {
            size_ = n;
            return;
        }
        T *p = static_cast<T *>(detail::space_alloc<Space>(n * sizeof(T)));
        detail::space_free<Space>(data_);
        data_ = p;
        capacity_ = n;
        size_ = n;
    }

    void clear() { size_ = 0; }
};

/*
 * Copy `src.size()` elements from `src` to `dst` (which must already be sized
 * to match). Host<->host is a memcpy; any device endpoint dispatches to the
 * backend copy, which infers the cudaMemcpy/hipMemcpy direction from the tags.
 */
template <typename T, typename DstSpace, typename SrcSpace>
void deep_copy(Array<T, DstSpace> &dst, const Array<T, SrcSpace> &src) {
    assert(dst.size() == src.size() && "deep_copy size mismatch");
    const std::size_t bytes = src.size() * sizeof(T);
    if (bytes == 0) return;

    if constexpr (std::is_same_v<DstSpace, HostSpace> &&
                  std::is_same_v<SrcSpace, HostSpace>) {
        std::memcpy(dst.data(), src.data(), bytes);
    } else if constexpr (std::is_same_v<DstSpace, HipSpace> ||
                         std::is_same_v<SrcSpace, HipSpace>) {
#ifdef MATSCIPY_ENABLE_HIP
        detail::copy_hip(dst.data(), DstSpace::device, src.data(),
                         SrcSpace::device, bytes);
#else
        static_assert(std::is_same_v<DstSpace, HostSpace>,
                      "deep_copy to/from HipSpace requires the HIP backend.");
#endif
    } else {
#ifdef MATSCIPY_ENABLE_CUDA
        detail::copy_cuda(dst.data(), DstSpace::device, src.data(),
                          SrcSpace::device, bytes);
#else
        static_assert(std::is_same_v<DstSpace, HostSpace> &&
                          std::is_same_v<SrcSpace, HostSpace>,
                      "deep_copy to/from CudaSpace requires the CUDA backend "
                      "(-DENABLE_CUDA=ON).");
#endif
    }
}

}  // namespace matscipy

#endif
