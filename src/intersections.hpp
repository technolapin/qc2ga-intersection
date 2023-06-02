#pragma once

#include "conic.hpp"

#include "svg_export.hpp"
		


#include<tuple>
#include "random.hpp"
#include "qcga_tools.hpp"



struct Algo
{
    inline static uint32_t COUNTER = 0;
    using Mvec = qc2ga::Mvec<double>;
    const Mvec inter;

    Conic<double> con_a, con_b; // two conics of the pencil
    Conic<double> con_deg, con_other; // two conics but one is degenerate and the other is significantly different
    double alpha, beta, w1, w2, u1, u2, v1, v2; // the parameters of the two lines inside con_deg (if it is not a point)
    double alpha1, alpha2; // the parameters of the two lines inside con_deg (if it is not a point)
    bool is_deg_a_point;
    RNG<double> & rng;
    std::vector<std::array<double, 2>> inter_pts;
    
    Algo(const Mvec & inter, RNG<double> & rng): inter(inter), rng(rng) {}
    
    void
    gen_degen_and_other()
    {
        // find two perpendicular conics
        const auto [pa] = rng.point_array<1>();
        const auto con_a_qcga = inter ^ pa;
        const auto con_b_qcga = inter ^ qcga_tools::right_complement(con_a_qcga);

        con_a = Conic(!con_a_qcga).normalize();
        con_b = Conic(!con_b_qcga).normalize();

        
        // non-normalized cubic equation
        const double aa = con_b.delta_3();
        const double bb = ((con_a+con_b).delta_3()+(con_a-con_b).delta_3())/2.0 - con_a.delta_3();
        const double cc = ((con_a+con_b).delta_3()-(con_a-con_b).delta_3())/2.0 - con_b.delta_3();
        const double dd = con_a.delta_3();
        
        // normalizing (for numerical stability)
        const double magnitude = std::sqrt(aa*aa + bb*bb + cc*cc + dd*dd);
        const double a = aa/magnitude;
        const double b = bb/magnitude;
        const double c = cc/magnitude;
        const double d = dd/magnitude;

        // discriminants
        const double D0 = b*b - 3.0*a*c;
        const double D1 = 2.0*b*b*b - 9.0*a*b*c + 27*a*a*d;
        const std::complex<double> D01 = D1*D1 - 4.0*std::pow(D0, 3);

        // sign adjustement is necessary here
        const double s0 = sgn_strict(D0);
        const double s1 = sgn_strict(D1);
        const std::complex<double> gamma_p =  s1*   std::pow((std::abs(D1)+std::sqrt(D01))/2.0, 1.0/3.0);
        const std::complex<double> gamma_m = s0*s1*std::pow(s0*((std::abs(D1)-std::sqrt(D01))/2.0), 1.0/3.0);

        // real root of the cubic equation
        const double lambda = -3.0*a;
        const double mu = b + std::real(gamma_m + gamma_p);


        con_deg   = (con_a*lambda + con_b*mu).normalize();
        con_other = (con_a*(-mu) + con_b*lambda).normalize();
        is_deg_a_point = con_deg.delta_2() > 0; // might happen
    }
    
    void
    factor_line_pair()
    {
        // normalizing for good measure
        const double K = std::sqrt(std::pow(con_deg.a-con_deg.b,2) + std::pow(con_deg.c,2))*sgn_strict(con_deg.f);
        const auto & [a,b,c,d,e,f] = con_deg/K;

        // adjoint of the hessian matrix of the conic
        std::array<std::array<double,3>,3> Hadj = {      b*f-e*e/4, (d*e)/4-(c*f)/2, (c*e)/4-(b*d)/2,
                                                   (d*e)/4-(c*f)/2, a*f-d*d/4      , (c*d)/4-(a*e)/2,
                                                   (c*e)/4-(b*d)/2,(c*d)/4-(a*e)/2,a*b-c*c/4       };
        std::array<double, 3> diag = { std::sqrt(-Hadj[0][0]),
                                       std::sqrt(-Hadj[1][1]),
                                       std::sqrt(-Hadj[2][2])};
        // finding a column with a non-nul diagonal element in the adjoint matrix
        const uint32_t i = std::distance(std::begin(diag), std::max_element(std::begin(diag), std::end(diag)));

        // this point is the intersection of the two lines
        // it is also the minors m11 m22 and m33 but with the good sign
        const std::array<double, 3> p = { Hadj[i][0]/diag[i],
                                          Hadj[i][1]/diag[i],
                                          Hadj[i][2]/diag[i]};
        
        const std::array<std::array<double,3>,3> N =
            {    a     , c/2.0-p[2], d/2.0+p[1],
             c/2.0+p[2],          b, e/2.0-p[0],
             d/2.0-p[1], e/2.0+p[0],          f};

        // some magnitudes
        const double nv0 = std::pow(N[0][0], 2) + std::pow(N[1][0], 2);
        const double nv1 = std::pow(N[0][1], 2) + std::pow(N[1][1], 2);
        const double nh0 = std::pow(N[0][0], 2) + std::pow(N[0][1], 2);
        const double nh1 = std::pow(N[1][0], 2) + std::pow(N[1][1], 2);

        // process a line parameter's to make them nicer
        // and to extract the angle
        const auto angle = [](const double & uu,
                              const double & vv,
                              const double & ww)
        {
            const double sgnw = sgn_strict(ww);
            const double norm = std::sqrt(uu*uu+vv*vv);
            double u = uu/norm*sgnw;
            double v = vv/norm*sgnw;
            double w = ww/norm*sgnw;
            const int su = sgn_strict(u);
            const int sv = sgn_strict(v);
            double alpha = std::atan2(std::abs(v),std::abs(u));
            if (su < 0) alpha = M_PI - alpha;
            if (sv < 0) alpha = - alpha;
            return std::array<double, 4>{u,v,w,alpha};
        };
            
        
        if (nv0 >= nv1)
        {
            const auto [u,v,w,a] = angle(N[0][0], N[1][0], N[2][0]);
            u1 = u;
            v1 = v;
            w1 = w;
            alpha1 = a;
        }
        else
        {
            const auto [u,v,w,a] = angle(N[0][1], N[1][1], N[2][1]);
            u1 = u;
            v1 = v;
            w1 = w;
            alpha1 = a;
        }
        if (nh0 >= nh1)
        {
            const auto [u,v,w,a] = angle(N[0][0], N[0][1], N[0][2]);
            u2 = u;
            v2 = v;
            w2 = w;
            alpha2 = a;
        }
        else
        {
            const auto [u,v,w,a] = angle(N[1][0], N[1][1], N[1][2]);
            u2 = u;
            v2 = v;
            w2 = w;
            alpha2 = a;
        }
        
    }

    void
    conic_lines_inter()
    {
        const auto conic_line_inter = [&](const double alpha, const double w)
        {
            const auto & [a,b,c,d,e,f] = con_other.rotate(-alpha);
            const double x = -w;
            const double delta = std::pow(c*x+e,2) - 4.0*b*(a*x*x+d*x+f);
            const double gamma = (c*x+e)/(2.0*b);
            const std::array<double, 2> p0 = {std::cos(alpha)*x + gamma*std::sin(alpha),
                                              std::sin(alpha)*x - gamma*std::cos(alpha)};
            const std::array<double, 2> u = {-std::sin(alpha)/(2.0*b), std::cos(alpha)/(2.0*b)};

            switch (sgn(delta))
            {
                case 1:
                    inter_pts.push_back({p0[0] + std::sqrt(delta)*u[0], p0[1] + std::sqrt(delta)*u[1]});
                    inter_pts.push_back({p0[0] - std::sqrt(delta)*u[0], p0[1] - std::sqrt(delta)*u[1]});
                    break;
                case 0:
                    inter_pts.push_back(p0);
                    break;
                default:
                    break;
            }

        };
        conic_line_inter(alpha1, w1);
        conic_line_inter(alpha2, w2);
    }

    void
    extract_pts_from_inter()
    {
        gen_degen_and_other();

        if (is_deg_a_point) 
        { // the degenerate conic is a point
            
            const auto [x, y] = con_deg.center();
            if (qcga_tools::norm_l2(inter ^ qc2ga::point<double>(x,y)) < EPSILON)
            { // the point is in the intersection
                inter_pts.push_back({x,y});
            }
        }
        else
        { // the degenerate conic is a pair of line

            factor_line_pair();
            conic_lines_inter();
        }
    }
    
    Algo &
    run()
    {
        extract_pts_from_inter();

        // for svg export
        if (false)
        {
            std::cout << "IMAGE NUMBER " << COUNTER << std::endl;
            std::stringstream ss;
            ss << "image-" << COUNTER++ << ".svg";
            svg::SVGRenderer rend(ss.str(), SOLARIZEDBASE00);
            rend.draw_conic(con_a, SOLARIZEDYELLOW);
            rend.draw_conic(con_b, SOLARIZEDORANGE);
            rend.draw_conic(con_deg, SOLARIZEDGREEN);
            rend.draw_conic(con_other, SOLARIZEDBLUE);
            for (const auto & [x,y]: inter_pts)
            {
                rend.draw_marker(x,y,SOLARIZEDBASE1);
            }
        }
        
        return *this;
    }
};


