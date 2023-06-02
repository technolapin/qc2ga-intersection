#pragma once

#include <qc2ga/Mvec.hpp>
namespace qcga_tools
{
    template<typename S>
    S
    norm_l0(const qc2ga::Mvec<S> & m)
    {
        S max = {};
        for (uint32_t i = 0; i < (1 << 8); ++i)
        {
            const S val = std::abs(m[i]);
            if (max < val) max = val;
        }
        return max;
    }
    template<typename S>
    S
    norm_l1(const qc2ga::Mvec<S> & m)
    {
        S sum = {};
        for (uint32_t i = 0; i < (1 << 8); ++i)
        {
            sum += std::abs(m[i]);
        }
        return sum;
    }

    template<typename S>
    S
    norm_l2(const qc2ga::Mvec<S> & m)
    {
        S sum = {};
        for (uint32_t i = 0; i < (1 << 8); ++i)
        {
            sum += std::pow(m[i], 2);
        }
        return std::sqrt(sum);
    }

    template<typename S>
    S
    sum_of_terms(const qc2ga::Mvec<S> & m)
    {
        S sum = {};
        for (uint32_t i = 0; i < (1 << 8); ++i)
        {
            sum += m[i];
        }
        return sum;
    }


    template<typename S>
    S
    scalar_product(const qc2ga::Mvec<S> & a,
                   const qc2ga::Mvec<S> & b)
    {
        S sum = {};
        for (uint32_t i = 0; i < (1 << 8); ++i)
        {
            sum += a[i]*b[i];
        }
        return sum;
    }

    // m ^ out = norm_l2(m)^2 I
    template<typename S>
    qc2ga::Mvec<S>
    right_complement(const qc2ga::Mvec<S> & m)
    {
        qc2ga::Mvec<S> out = {};
        for (uint32_t i = 0; i < (1 << 8); ++i)
        {
            const uint32_t j = i ^ 255;
            qc2ga::Mvec<S> a; a[i] = S(1);
            qc2ga::Mvec<S> b; b[j] = S(1);
            const auto ab = a^b;
            out[j] = m[i]*sgn(ab[255]);
        }
        return out;
    }
    
}
