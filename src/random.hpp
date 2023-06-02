#pragma once

#include <qc2ga/Mvec.hpp>
#include "qc2gaTools.hpp"
#include <vector>
#include <array>
#include <random>

template<typename S>
struct RNG
{
private:
    std::mt19937 gen;
    std::uniform_real_distribution<S> distrib;
public:
    void
    print_state() const
    {
        std::cout << gen << std::endl;
    }
    
    RNG(const int seed):
        gen( seed),
        distrib(S(0), S(1))
    {}
    
    /** Returns a scalar in [0, 1] */
    S
    scalar(const S min, const S max)
    {
        return distrib(gen)*(max-min) + min;
    }

    template<uintptr_t SIZE>
    std::array<S, SIZE>
    scalar_array(const S min, const S max)
    {
        std::array<S, SIZE> ar;
        for (S & val: ar)
        {
            val = scalar(min, max);
        }
        return ar;
    }
    
    qc2ga::Mvec<S>
    point()
    {
        return qc2ga::point<S>(scalar(-1, 1), scalar(-1, 1));
    }

    template<uintptr_t SIZE>
    std::array<qc2ga::Mvec<S>, SIZE>
    point_array()
    {
        std::array<qc2ga::Mvec<S>, SIZE> ar;
        for (uint32_t i = 0; i < SIZE; ++i)
        {
            ar[i] = point();
        }
        return ar;
    }
};

    
