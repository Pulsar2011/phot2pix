//
// Created by gillard on 2026-11-25.
// Centre de Physique des Particule de Marseille (CPPM)
// Licence under PolyForm Noncommercial License 1.0.0 https://polyformproject.org/licenses/noncommercial/1.0.0
//  You may share and adapt this code with attribution, Non Comercial use only.
//
// DetectorTools.h
//
/// @file DetectorTools.h
/// @brief Tools shared by the detectors and their effects: threaded loop, typed access to a FITS cube.

#ifndef P2P_DETECTOR_TOOLS_H
#define P2P_DETECTOR_TOOLS_H

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <thread>
#include <valarray>
#include <vector>

#include <DSTfits/FITSimg.h>

namespace phot2pix
{
    /// @brief Run work(begin, end) over [0, n), split into contiguous chunks on nThreads std::threads
    ///        (so the compiler can vectorise each loop).
    /// @details A single call when the job is too small to be worth threading.
    /// @param n        Number of elements.
    /// @param nThreads Number of threads.
    /// @param work     Called with the bounds of a chunk.
    template <typename Work>
    void parallelFor(size_t n, size_t nThreads, const Work& work)
    {
        if (nThreads <= 1 || n < 100 * nThreads)
        {
            work(size_t(0), n);
            return;
        }

        const size_t chunk = (n + nThreads - 1) / nThreads;
        std::vector<std::thread> workers;
        workers.reserve(nThreads);
        for (size_t b = 0; b < n; b += chunk)
            workers.emplace_back(work, b, std::min(n, b + chunk));
        for (auto& w : workers)
            w.join();
    }

    /// @brief Call f with the cube's data, as the std::valarray<T>& of the type the cube is stored in
    ///        (8- to 64-bit integers, float, double).
    /// @return False if it is none of them.
    template <typename F>
    bool withCubeData(DSL::FITScube& cube, F&& f)
    {
        if (auto* d = cube.GetData<uint8_t>())  { f(*d); return true; }
        if (auto* d = cube.GetData<uint16_t>()) { f(*d); return true; }
        if (auto* d = cube.GetData<uint32_t>()) { f(*d); return true; }
        if (auto* d = cube.GetData<uint64_t>()) { f(*d); return true; }
        if (auto* d = cube.GetData<int8_t>())   { f(*d); return true; }
        if (auto* d = cube.GetData<int16_t>())  { f(*d); return true; }
        if (auto* d = cube.GetData<int32_t>())  { f(*d); return true; }
        if (auto* d = cube.GetData<int64_t>())  { f(*d); return true; }
        if (auto* d = cube.GetData<float>())    { f(*d); return true; }
        if (auto* d = cube.GetData<double>())   { f(*d); return true; }
        return false;
    }
} // namespace phot2pix

#endif // P2P_DETECTOR_TOOLS_H
