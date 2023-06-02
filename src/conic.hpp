#pragma once

#include "consts.hpp"

#include <qc2ga/Mvec.hpp>

template<typename S>
struct Conic
{
    S a,b,c,d,e,f;
    Conic() = default;
    Conic(const S a,
          const S b,
          const S c,
          const S d,
          const S e,
          const S f): a(a), b(b), c(c), d(d), e(e), f(f)
    {} 

/**
       Extracts the parameters from a dual quadric !(q^Ito)
     */
    Conic(const qc2ga::Mvec<S> & dual_conic)
    {
        a = -0.5*dual_conic[qc2ga::E01];
        b = -0.5*dual_conic[qc2ga::E02];
        c = -dual_conic[qc2ga::E03];
        d = dual_conic[qc2ga::E1];
        e = dual_conic[qc2ga::E2];
        f = -2.0*dual_conic[qc2ga::Ei1]; // can use Ei2 too
    }
    
    /**
       Returns !(q^Ito)
     */
    qc2ga::Mvec<S>
    qcga() const
    {
        qc2ga::Mvec<S> m;
        m[qc2ga::E01] = - 2.0*a;
        m[qc2ga::E02] = - 2.0*b;
        m[qc2ga::E03] = -     c;
        m[qc2ga::E1 ] =       d;
        m[qc2ga::E2 ] =       e;
        m[qc2ga::Ei1] = - 0.5*f;
        m[qc2ga::Ei2] = m[qc2ga::Ei1];
        return m;
    }


    S
    evaluate(const S x,
             const S y) const
    {
        return a*x*x + b*y*y + c*x*y + d*x + e*y + f;
    }

    /**
       |a   c/2|
       |c/2 b  |
     */
    S
    delta_2() const
    {
        return a*b - 0.25*c*c;
    }
    /**
       |a   c/2 d/2|
       |c/2 b   e/2|
       |d/2 e/2 f  |
     */
    S
    delta_3() const
    {
        return (a*b - c*c/4.0)*f + (c*e*d - b*d*d - a*e*e)/4.;
    }
    Conic 
    operator+(const Conic & q) const
    {
        return {
            a + q.a,b + q.b,c + q.c,d + q.d,e + q.e,f + q.f
        };
    }
    Conic 
    operator-(const Conic & q) const
    {
        return {
            a - q.a,b - q.b,c - q.c,d - q.d,e - q.e,f - q.f
        };
    }
    Conic
    operator*(const S s) const
    {
        return {
            a*s,b*s,c*s,d*s,e*s,f*s
        };
    }
    Conic
    operator*(const Conic rhs) const
    {
        return {
            a*rhs.f+d*rhs.d+f*rhs.a,
            b*rhs.f+e*rhs.e+f*rhs.b,
            c*rhs.f+d*rhs.e+e*rhs.d+f*rhs.c,
            d*rhs.f+f*rhs.d,
            e*rhs.f+f*rhs.e,
            f*rhs.f
        };
    }
    Conic
    operator/(const S s) const
    {
        return (*this)*(1.0d/s);
    }

    S
    norm() const
    {
        return std::sqrt(a*a + b*b + c*c + d*d + e*e +f*f);
    }
    Conic
    normalize() const
    {
        return (*this)/norm();
    }


    Conic
    rotate(const S t) const
    {
        Conic transformed(*this);

        transformed.a = ( b*std::pow(std::sin(t), 2) - c*std::cos(t)*std::sin(t) + a*std::pow(std::cos(t), 2));

        transformed.b = ( a*std::pow(std::sin(t), 2) + c*std::cos(t)*std::sin(t) + b*std::pow(std::cos(t), 2));
        transformed.c = (-c*std::pow(std::sin(t), 2) - 2*b*std::cos(t)*std::sin(t) + 2*a*std::cos(t)*std::sin(t) + c*std::pow(std::cos(t), 2));
        transformed.d = (d*std::cos(t) - e*std::sin(t));
        transformed.e = (e*std::cos(t) + d*std::sin(t));
        
        transformed.f = f;
        
        return transformed;
    }

    
    std::array<S, 2>
    center() const
    {
        const S det = -S(4)*delta_2();
        return {(S(2)*b*d-c*e)/det, (S(2)*a*e - c*d)/det};
    }
    Conic
    translate(const S dx,
              const S dy) const
    {
        Conic<S> transformed(*this);
        transformed.f += a*dx*dx + b*dy*dy + c*dx*dy - (d*dx + e*dy);
        transformed.d += -(2*a*dx + c*dy);
        transformed.e += -(2*b*dy + c*dx);
        return transformed;
    }
    S
    angle() const
    {
        return std::atan2(c,a-b)/S(2);
    }

    Conic<S>
    dilate(const S sx, const S sy) const
    {
        return {
            a/std::pow(sx,2),
            b/std::pow(sy,2),
            c/(sx*sy),
            d/sx,
            e/sy,
            f
        };
    }
    /**
       Set the near-zero parameters to 0
    */
    Conic<S>
    epurate() const
    {
        return {
            (std::abs(a) < EPSILON)? S(0) : a,
            (std::abs(b) < EPSILON)? S(0) : b,
            (std::abs(c) < EPSILON)? S(0) : c,
            (std::abs(d) < EPSILON)? S(0) : d,
            (std::abs(e) < EPSILON)? S(0) : e,
            (std::abs(f) < EPSILON)? S(0) : f
        };
    }
    int32_t
    signature() const
    {
        const int32_t sgn_det2 =  sgn(delta_2());
        const int32_t sgn_det3 =  sgn(delta_3());
        const int32_t u = sgn(a+b + std::sqrt(std::pow(a+b, 2) - S(4)*delta_2()));
        const int32_t v = sgn(a+b - std::sqrt(std::pow(a+b, 2) - S(4)*delta_2()));
        const int32_t w = sgn_det2*sgn_det3;

        if ((u > 0 && v > 0 && w > 0) || (u < 0 && v < 0 && w < 0))
        {
            return 1;
        }
        else if (sgn_det3 == 0)
        {
            return 0;
        }
        else
        {
            return -1;
        }
        
    }
    Conic<S>
    normalized() const
    {
        return (*this)/norm();
    }
    
    std::pair<std::array<S, 5>, Conic<S>>
    canonical_transforms() const
    {
        /*
          rotation THEN translation THEN dilation
         */
        const int32_t sgn_det2 = sgn(delta_2());
        const int32_t sgn_det3 = sgn(delta_3());
        //   S sgn_f = sgn(f);
        S alpha = angle();
        // chirality issues
        if (sgn_det3 > 0)
        {
             alpha += M_PI/S(2);
        }
        const Conic<S> con_aligned = rotate(-alpha);
        const int32_t sig = signature();
        if (sig > 0)
        { // imaginary conic (useless)
            std::cout << "IMAGINARY CONIC!\n";
            return {};
        }
        else if (sig == 0)
        { // degenerate conic
            // NOT SUPPORTED YET (I have no use for this)
            return {};
            /*
            if (sgn_det2 != 0)
            { // intersecting lines (x^2 - y^2 = 0) or point (x^2 + y^2 = 0)
                const auto [cx,cy] = con_aligned.center();
                const Conic<S> con_centered = con_aligned.translate(-cx,-cy);
                S sx = std::sqrt(std::abs(con_centered.a));
                S sy = std::sqrt(std::abs(con_centered.b));
                const Conic<S> con_scaled = con_centered.dilate(sx,sy);
                return {
                    {alpha, -cx, -cy, sx, sy},
                    con_scaled*(-sgn_det2)
                };
            }
            else
            { // pair of lines (x^2 - 1 = 0)
                std::cout << "PARRALLEL\n";
                const S cx = -con_aligned.d / (S(2)*con_aligned.a);
                const S cy = -(S(4)*con_aligned.a*con_aligned.f - std::pow(con_aligned.d,2))/(S(4)*con_aligned.e*con_aligned.a);
                const Conic<S> con_centered = con_aligned.translate(-cx,-cy);
                S sx = std::sqrt(std::abs(con_centered.a));
                S sy = -con_centered.e;
                if (std::abs(sx) < EPSILON) sx = S(1);
                if (std::abs(sy) < EPSILON) sy = S(1);
                const Conic<S> con_scaled = con_centered.dilate(sx,sy);
                return {
                    {alpha, -cx, -cy, sx, sy},
                    con_scaled
                };
            }
            */
        }
        else
        { // real conic
            if (sgn_det2 != 0)
            { // ellipse (x^2 + y^2 - 1 = 0) or hyperbola (x^2 - y^2 + 1 = 0)
                const auto [cx,cy] = con_aligned.center();
                const Conic<S> con_centered = con_aligned.translate(-cx,-cy);
                S sx = std::sqrt(std::abs(con_centered.a/con_centered.f));
                S sy = std::sqrt(std::abs(con_centered.b/con_centered.f));
                const Conic<S> con_scaled = con_centered.dilate(sx,sy) * (S(1)/con_centered.f);
                return {
                    {alpha, -cx, -cy, sx, sy},
                    con_scaled*(-sgn_det2)
                };
            }
            else
            { // parabola (x^2 - y = 0)
                const S cx = -con_aligned.d / (S(2)*con_aligned.a);
                const S cy = -(S(4)*con_aligned.a*con_aligned.f - std::pow(con_aligned.d,2))/(S(4)*con_aligned.e*con_aligned.a);
                const Conic<S> con_centered = con_aligned.translate(-cx,-cy);
                S sx = std::sqrt(std::abs(con_centered.a));
                S sy = -con_centered.e;
                if (std::abs(sx) < EPSILON) sx = S(1);
                if (std::abs(sy) < EPSILON) sy = S(1);
                const Conic<S> con_scaled = con_centered.dilate(sx,sy);
                return {
                    {alpha, -cx, -cy, sx, sy},
                    con_scaled
                };
            }
        }
        

        
   }
};



template<typename S>
std::ostream& operator<<(std::ostream& os, const Conic<S>& c)
{
    os << "["  << c.a
       << ", " << c.b
       << ", " << c.c
       << ", " << c.d
       << ", " << c.e
       << ", " << c.f << "]";

    return os;
}



/**
   Factor a degenerate conic into two lines,
   assuming C = (u_1 x + v_1 y + w_1)(u_2 x + v_2 y + w_2)
   The lines extracted are returned as degenerate conics of the form ux+vy+w=0
*/
template<typename S>
std::array<Conic<S>, 2>
factor_pair_of_line(const Conic<S> conic)
{
    /*
      C = (u_1 x + v_1 y + w_1)(u_2 x + v_2 y + w_2)
      u_1 = cos(alpha     )    v_1 = sin(alpha     )
      u_2 = cos(alpha+beta)    v_2 = sin(alpha+beta)
     */
    const S norm = std::sqrt(std::pow(conic.a-conic.b,2) + conic.c*conic.c);
    const S an = conic.a/norm;
    const S bn = conic.b/norm;
    const S cn = conic.c/norm;
    const S dn = conic.d/norm;
    const S en = conic.e/norm;
    const S fn = conic.f/norm;

    // 2alpha + beta
    const S angle_sum = std::atan2(conic.c, conic.a-conic.b);
    const S beta = std::acos(an+bn);
    const S alpha = 0.5*(angle_sum - beta);
    // alpha + beta
    const S alpha_beta = 0.5*(angle_sum + beta);

    // rotation by -alpha
    const S cr = cn*std::cos(2*alpha) - (an-bn)*std::sin(2*alpha);
    const S er = en*std::cos(alpha) - dn*std::sin(alpha);
    const S fr = fn;
    
    const S u_1 = cos(alpha);
    const S v_1 = sin(alpha);
    const S u_2 = cos(alpha_beta);
    const S v_2 = sin(alpha_beta);
    const S w_1 = er/cr;
    const S w_2 = fr/w_1;
    return { Conic<S>(0,0,0,u_1,v_1,w_1),
             Conic<S>(0,0,0,u_2,v_2,w_2) };
}
