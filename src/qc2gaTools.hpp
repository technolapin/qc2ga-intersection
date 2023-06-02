// Copyright (c) 2018 by University Paris-Est Marne-la-Vallee
// qc2ga2Tools.hpp
// Authors: Vincent Nozick and Stephane Breuils 
// Contact: vincent.nozick@u-pem.fr


/// \file qc2ga2Tools.hpp
/// \author Vincent Nozick, and Stephane Breuils
/// \brief some useful functions when using Quadric Conformal Geometric Algebra of R^2. Use this file if you generated the lib using the "qc2ga.conf" configuration file. Quadric conformal geometric algebra of R^2 is inspired from the following paper: Stéphane Breuils, Vincent Nozick, Akihiro Sugimoto, Eckhard Hitzer, Quadric Conformal Geometric Algebra of R9,6, Advances in Applied Clifford Algebras, Springer Verlag, 2018, 28 (2), pp.35 


// Anti-doublon
#ifndef QC2GA_TOOLS_HPP__
#define QC2GA_TOOLS_HPP__
#pragma once

// External includes
#include <string>
#include <vector>
#include <random>
#include <chrono>

// Internal Includes
#include <qc2ga/Mvec.hpp>
#include <Eigen/Dense>        // rank in display conic
#include <Eigen/Eigenvalues>  // eigen solver

/// \namespace grouping the multivectors object
namespace qc2ga{

    // default seed
    const auto defaultSeed = std::chrono::system_clock::now().time_since_epoch().count();

    // random engine
    std::default_random_engine generator(defaultSeed);

    // reset random engine seed
    void resetRandomEngineSeed(const unsigned int seed){
        generator = std::default_random_engine(seed);
    }

    // for random points
    template<typename T>
    std::uniform_real_distribution<T> uniformRealDistribution(-1.0,1.0);

    // Ii = ei1^ei2^ei3
    template<typename T>
    inline Mvec<T> Ii() { return qc2ga::ei1i2i3<T>(); }

    // I0 = e01^e02^e03
    template<typename T>
    inline Mvec<T> I0() { return qc2ga::e010203<T>(); }

    // Iit = (ei1-ei2)^ei3
    template<typename T>
    inline Mvec<T> Iit() { return (qc2ga::ei1<T>() - qc2ga::ei2<T>()) ^ qc2ga::ei3<T>(); }

    // I0t = (e01-e02)^e03
    template<typename T>
    inline Mvec<T> I0t() { return (qc2ga::e01<T>() - qc2ga::e02<T>()) ^ qc2ga::e03<T>(); }

    // ei = 1/2(ei1+ei2)
    template<typename T>
    inline Mvec<T> ei() { return 0.5*(qc2ga::ei1<T>() + qc2ga::ei2<T>()); }

    // e0 = e01+e02
    template<typename T>
    inline Mvec<T> e0() { return qc2ga::e01<T>() + qc2ga::e02<T>(); }


    /// \brief build a point from a vector
    /// \param x vector component related to e1
    /// \param y vector component related to e2
    /// \return a multivector corresponding to a point of qc2ga
    template<typename T>
    qc2ga::Mvec<T> point(const T &x, const T &y){
        qc2ga::Mvec<T> mv;
        mv[qc2ga::E1] = x;
        mv[qc2ga::E2] = y;
        mv[qc2ga::Ei1] = 0.5*(x*x);
        mv[qc2ga::Ei2] = 0.5*(y*y);
        mv[qc2ga::Ei3] = (x*y);
  		mv[qc2ga::E01] = mv[qc2ga::E02] = 1.0;
        return mv;
    }

    /// \brief build a random point with Euclidean coordinates ranging in [-1,1]
    /// \return a point 
    template<typename T>
    qc2ga::Mvec<T> randomPoint(){
        return point<T>(uniformRealDistribution<T>(generator), uniformRealDistribution<T>(generator));
    }

    template<typename T>
    qc2ga::Mvec<T> circle(const T &cx, const T &cy, const T &radius){
        qc2ga::Mvec<T> mv = qc2ga::point<T>(cx, cy) - radius * radius * (qc2ga::ei1<T>() + qc2ga::ei2<T>() ) / 6.0;
        return mv;
    }

    /// a conic has the form (a b c d e f)
    /// with ax^2 + by^2 + cxy + dxw + eyw + fw^2 = 0 
    /// \todo never tested !!!
    template<typename T>
    std::vector<T> ga2DualConic(const qc2ga::Mvec<T> &mv){

        std::vector<T> conic(10);
        conic[0] = - mv[qc2ga::E01] / 2.0;
        conic[1] = - mv[qc2ga::E02] / 2.0;
        conic[2] = - mv[qc2ga::E03] ;

        conic[3] =  mv[qc2ga::E1];
        conic[4] =  mv[qc2ga::E2];

        conic[5] = - 2 * mv[qc2ga::Ei1];

        return conic;
    }


    /// a conic has the form (a b c d e f)
    /// with ax^2 + by^2 + cxy + dxw + eyw + fw^2 = 0 
    /// \todo never tested !!!
    template<typename T>
    qc2ga::Mvec<T> dualConic2ga(const T &a, const T &b, const T &c, const T &d, const T &e, const T &f){
        qc2ga::Mvec<T> mv;
        mv[qc2ga::E01] = - 2.0 * a; // x^2
        mv[qc2ga::E02] = - 2.0 * b; // y^2
        mv[qc2ga::E03] = - c;       // xy
        mv[qc2ga::E1] = d;          // x
        mv[qc2ga::E2] = e;          // y

        mv[qc2ga::Ei1] = mv[qc2ga::Ei2] = - f / 2.0;
        
        return mv;
    }

    /// a conic has the form (a b c d e f)
    /// with ax^2 + by^2 + cxy + dxw + eyw + fw^2 = 0 
    /// \todo never tested !!!
    template<typename T>
    qc2ga::Mvec<T> lineToDualGA(const T a, const T b, const T c){
        qc2ga::Mvec<T> mv;
        mv[qc2ga::E1] = a;          // x
        mv[qc2ga::E2] = b;          // y

        mv[qc2ga::Ei1] = mv[qc2ga::Ei2] = mv[qc2ga::Ei3] = - c / 2.0; // 1
        
        return mv;
    }

    /// a conic has the form (a b c d e f)
    /// with ax^2 + by^2 + cxy + dxw + eyw + fw^2 = 0 
    /// \todo never tested !!!
    template<typename T>
    qc2ga::Mvec<T> dualConic2ga(const std::vector<T> &conic){
        return dualConic2ga<T>(conic[0], conic[1], conic[2], conic[3], conic[4], conic[5]);
    }

	
    /// a conic has the form (a b c d e f)
    /// with ax^2 + by^2 + cxy + dxw + eyw + fw^2 = 0 
    /// \todo never tested !!!
    // template<typename T>
    // std::vector<T> ga2Conic(const qc2ga::Mvec<T> &mv){

    //     std::vector<T> conic(10);
    //     conic[0] = - mv[qc2ga::E1230102i203i304i405i506i6] / 2.0; // dual(e01)
    //     conic[1] = - mv[qc2ga::E12301i10203i304i405i506i6] / 2.0; // dual(e02)
    //     conic[2] = - mv[qc2ga::E12301i102i20304i405i506i6] / 2.0; // dual(e03)
    //     conic[3] = - mv[qc2ga::E12301i102i203i30405i506i6]; // dual(e04)
    //     conic[4] = - mv[qc2ga::E12301i102i203i304i40506i6]; // dual(e05)
    //     conic[5] = - mv[qc2ga::E12301i102i203i304i405i506]; // dual(e06)
		
    //     conic[6] = -mv[qc2ga::E2301i102i203i304i405i506i6]; // dual(e1)
    //     conic[7] =  mv[qc2ga::E1301i102i203i304i405i506i6]; // dual(e2)
    //     conic[8] = -mv[qc2ga::E1201i102i203i304i405i506i6]; // dual(e3)
		
    //     conic[9] = 3 * mv[qc2ga::E123i102i203i304i405i506i6]; // dual(ei1)

    //     return conic;
    // }

    // /// a conic has the form (a b c d e f g h i j)
    // /// with ax^2 + by^2 + cz^2 + dxy + exz + fyz + gxw + hyw + izw + jw^2 = 0
    // /// defines the primal form of a conic
    // template<typename T>
    // qc2ga::Mvec<T> conic2ga(const T a, const T b, const T c, const T d, const T e, const T f, const T g, const T h, const T i, const T j){
    //     qc2ga::Mvec<T> mv;

    //     mv[qc2ga::E1230102i203i304i405i506i6] = - 2.0 * a; // dual(e01) component
    //     mv[qc2ga::E12301i10203i304i405i506i6] = - 2.0 * b; // dual(e02) component
    //     mv[qc2ga::E12301i102i20304i405i506i6] = - 2.0 * c; // dual(e03) component
    //     mv[qc2ga::E12301i102i203i30405i506i6] = - d; // dual(e04) component
    //     mv[qc2ga::E12301i102i203i304i40506i6] = - e; // dual(e05) component
    //     mv[qc2ga::E12301i102i203i304i405i506] = - f; // dual(e06) component


    //     mv[qc2ga::E2301i102i203i304i405i506i6] = g; // dual(e1) component 
    //     mv[qc2ga::E1301i102i203i304i405i506i6] = -h; // dual(e2) component
    //     mv[qc2ga::E1201i102i203i304i405i506i6] = i; // dual(e3) component

    //     mv[qc2ga::E123i102i203i304i405i506i6] = mv[qc2ga::E12301i1i203i304i405i506i6] = mv[qc2ga::E12301i102i2i304i405i506i6] = j / 3.0; // dual(eij) component 
        
    //     return mv;
    // }

    /// a conic has the form (a b c d e f)
    /// with ax^2 + by^2 + cxy + dxw + eyw + fw^2 = 0 
    /// \todo never tested !!!    
    template<typename T>
    qc2ga::Mvec<T> conic2ga(const std::vector<T> &conic){
        return qc2ga::conic2ga<T>(conic[0], conic[1], conic[2], conic[3], conic[4], conic[5]);
    }






} // namespace

#endif // projection_inclusion_guard
