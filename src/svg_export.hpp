#pragma once

#include "conic.hpp"
// #include "renderer.hpp"


#include <iostream>
#include <iomanip>
		
#include <fstream>
#include <cmath>
#include <functional>

namespace svg
{
    template<typename S>
    std::vector<std::vector<std::array<double, 2>>>
    sample_points_conics(const Conic<S> & con, const uint32_t n_pts, const S MIN_DIST = 0.01,
                         const S MAX_RADIUS=S(2.0))
    {
       

        std::vector<std::vector<std::array<double, 2>>> pts;
        std::vector<std::array<double, 2>> subcurve;
        const auto & [transforms, con_base] = con.normalized().canonical_transforms();
        const auto & [alpha,tx,ty,sx,sy] = transforms;
        
        const S angle_step = S(2.0*M_PI)/S(n_pts);
        const auto &[a,b,c,d,e,f] = con_base;

        const auto & foo = [&](const double t)
        {
            const S r = std::sqrt(S(2)*(-f/(a+b + (a-b)*std::cos(S(2)*t))));

            S x0 = std::cos(t)*r;
            S y0 = std::sin(t)*r;
            const S x_s = x0/sx;
            const S y_s = y0/sy;
            const S x_t = x_s - tx;
            const S y_t = y_s - ty;
            const double x = double( std::cos(-alpha)*x_t + std::sin(-alpha)*y_t);
            const double y = double(-std::sin(-alpha)*x_t + std::cos(-alpha)*y_t);
            return std::array<double, 2>{x,y};
        };
        
        for (uint32_t i = 0; i < n_pts; ++i)
        {
            const S t = angle_step * S(i);

            const auto [x,y] = foo(t);

            if (std::isnan(x) || std::isnan(y) || std::abs(x) > MAX_RADIUS || std::abs(y) > MAX_RADIUS)
            {
                if (subcurve.size())
                {
// TODO: dichotomy between the current angle and the last to find the point at the limit of the area 
                    {
                        // std::cout << "\nt = " << t << std::endl;
                        // std::cout << "LAST X Y: " << subcurve.back()[0] << " " << subcurve.back()[1] << std::endl;
                        double t_min = t - angle_step;
                        double t_max = t;
                        double t_mid = (t_min + t_max)*0.5;
//                        const auto &[x0,y0] = foo(t_min);
                        // std::cout << "XMIN YMIN = " << x0 << " " << y0 << std::endl;
                        // std::cout << "XMAX YMAX = " << x << " " << y << std::endl;
                        auto [x,y] = foo(t_mid);
                        // std::cout << "T_MID=" << t_mid << std::endl;
                        double max_coord = std::max(std::abs(x), std::abs(y));
                        uint32_t i = 0;
                        do
                        {
                            
                            // std::cout << "T Min MID Max: " << t_min << " " << t_mid << " " << t_max << std::endl;
                            ++i;
                            if (max_coord > MAX_RADIUS)
                            {
                                t_max = t_mid;
                            }
                            else
                            {
                                t_min = t_mid;
                            }
                            t_mid = (t_max + t_min)*0.5;
                            auto [x2,y2] =foo(t_mid);
                            x = x2;
                            y = y2;
                            max_coord = std::max(std::abs(x), std::abs(y));
                            if (i >= 2000) break;
                        }
                        while ( std::isnan(x) || std::isnan(y) || max_coord > MAX_RADIUS || max_coord < MAX_RADIUS-0.001);
                        if (!( std::isnan(x) || std::isnan(y) || max_coord > MAX_RADIUS || max_coord < MAX_RADIUS-0.001))
                        {
                            
                            
                            subcurve.push_back({x,y});
                        }
                        
                    }

                    
                    
                    pts.push_back(subcurve);
                    subcurve = {};
                }
                continue;
            }
            

            if (subcurve.size() > 0)
            {
                const auto & [xprev, yprev] = subcurve.back();
                if (std::sqrt(std::pow(x - xprev, 2) + std::pow(y - yprev, 2)) < MIN_DIST) continue;
            }
            subcurve.push_back({x,y});
        }
        if (subcurve.size())
        {
            const auto & [x0, y0] = subcurve[0];
            const auto & [x1, y1] = subcurve.back();
            if (std::sqrt(std::pow(x0-x1,2) + std::pow(y0-y1,2)) < S(0.1))
            {
                subcurve.push_back({x0,y0});
            }
            pts.push_back(subcurve);
        }

        
        return pts;
    }
    template<typename S>
    void
    export_to_svg(const std::string & filename,
                  const Conic<S> & con,
                  const uint32_t n_pts)
    {
        const auto  curves = sample_points_conics(con, n_pts);

        const S box_radius = S(100);
        const S stroke_width = 0.4;
        std::ofstream myfile;
        myfile.open (filename);

        myfile << "<?xml version=\"1.0\" encoding=\"UTF-8\" standalone=\"no\"?>\n";
        myfile << "<!DOCTYPE svg PUBLIC \"-//W3C//DTD SVG 1.1//EN\" \"http://www.w3.org/Graphics/SVG/1.1/DTD/svg11.dtd\">\n";
        myfile << "<svg xmlns=\"http://www.w3.org/2000/svg\" viewBox=\"" <<-S(2)*box_radius<< " " << -S(2)*box_radius << " " << S(4)*box_radius << " " <<  S(4)*box_radius << "\">\n";
        myfile << "<title>Polynomial of degree 3</title>\n";
        myfile << "<desc> description ou jsp quoi </desc>\n";

        //   myfile << "<html>\n<body>\n";

//        myfile << "<h1>" << "svg name or something" << "</h1>\n";
        //       myfile << "<svg width=\"100\" height=\"100\">\n";
        {
            //    myfile << "<circle cx=\"0\" cy=\"0\" r=\"0.5\" stroke=\"green\" stroke-width=\"0.01\" fill=\"yellow\" />\n";
        }
        for (const auto & curve: curves)
        {
            for (uint32_t i = 1; i < curve.size(); ++i)
            {
                const auto & [x0, y0] = curve[i-1];
                const auto & [x1, y1] = curve[i];
                myfile << "<line x1=\"" << x0*box_radius << "\" y1=\""<<-y0*box_radius<<"\" x2=\""<<x1*box_radius<<"\" y2=\""<<-y1*box_radius<<"\" stroke=\"black\" stroke-width=\""<<stroke_width <<"\" />\n";
            }
        }
        
        myfile << "</svg>\n";
        //myfile << "</body>\n</html>\n";
        myfile.close();
        
    }

    
    struct SVGRenderer
    {
    private:
        std::ofstream file;
        double box_radius;
        double x_min;
        double y_min;
        double w;
        double h;
        
    public:
        SVGRenderer(const std::string & filename,
                    const uint32_t background)
        {
            box_radius = 200.0;
            x_min = -box_radius;
            y_min = -box_radius;
            w = 2.0*box_radius;
            h = 2.0*box_radius;
            file.open(filename);
            file << "<?xml version=\"1.0\" encoding=\"UTF-8\" standalone=\"no\"?>\n";
            file << "<!DOCTYPE svg PUBLIC \"-//W3C//DTD SVG 1.1//EN\" \"http://www.w3.org/Graphics/SVG/1.1/DTD/svg11.dtd\">\n";
                file << "<svg xmlns=\"http://www.w3.org/2000/svg\" viewBox=\"" <<x_min<< " " << y_min << " " << w << " " << h << "\" style=\"background-color:#"<< std::hex << std::setfill('0') << std::setw(6) << background << std::dec <<"\">\n";
            file << "<title>Polynomial of degree 3</title>\n";
            file << "<desc> description ou jsp quoi </desc>\n";
        }
        
        ~SVGRenderer()
        {
            const double stroke_radius = box_radius/50.;
            draw_line(-2.,-2.,-2., 2.,SOLARIZEDBASE02, stroke_radius);
            draw_line(-2., 2., 2., 2.,SOLARIZEDBASE02, stroke_radius);
            draw_line( 2., 2., 2.,-2.,SOLARIZEDBASE02, stroke_radius);
            draw_line( 2.,-2.,-2.,-2.,SOLARIZEDBASE02, stroke_radius);
            file << "</svg>\n";
            file.close();
        }
        SVGRenderer &
        draw_line(const double x0, const double y0,
                  const double x1, const double y1,
                  const uint32_t colour, const double stroke_width)
        {
            file << "<line x1=\"" <<  x0*box_radius/2.
                 <<    "\" y1=\"" << -y0*box_radius/2.
                 <<    "\" x2=\"" <<  x1*box_radius/2.
                 <<    "\" y2=\"" << -y1*box_radius/2.
                 << "\" stroke=\"#" << std::hex <<  std::hex << std::setfill('0') << std::setw(6) <<colour << std::dec
                 << "\" stroke-width=\""<<stroke_width <<"\" />\n";
            return *this;

        }
        /** the space is actually [-2, 2] */
        template<typename S>
        SVGRenderer &
        draw_conic(const Conic<S> & con, const uint32_t colour, const double thickness = 1.0)
        {
            const double stroke_width = box_radius/100.*thickness;
            const double MAX_RADIUS = 2.0;
            const double MIN_DIST = 0.01;
            // must normalize or the result won't be precise
            if (sgn(con.normalize().delta_3()) == 0)
            {
                if (sgn(con.delta_2()) > 0)
                { // Point
                    const auto &[x,y] = con.center();
                    if (std::max(std::abs(x), std::abs(y)) < MAX_RADIUS)
                        draw_marker(x,y,colour);
                }
                else
                { // pair of lines
                    const auto &[l1,l2] = factor_pair_of_line(con);
                    const auto draw_ln = [&](const double u,
                                             const double v,
                                             const double w)
                    {
                        const std::array<double, 2> p1 = {-MAX_RADIUS, -(w-MAX_RADIUS*u)/v};
                        const std::array<double, 2> p2 = { MAX_RADIUS, -(w+MAX_RADIUS*u)/v};
                        const std::array<double, 2> p3 = {-(w-MAX_RADIUS*v)/u, -MAX_RADIUS};
                        const std::array<double, 2> p4 = {-(w+MAX_RADIUS*v)/u, MAX_RADIUS};
                        std::vector<std::array<double,2>> pts;
                        for (const auto &[x,y]: {p1,p2,p3,p4})
                        {
                            if (std::isnan(x) || std::isnan(y))
                                continue;
                            if (std::abs(std::max(std::abs(x), std::abs(y)) - MAX_RADIUS) < 0.001)
                            {
                                pts.push_back({x,y});
                            }
                        }
                        if (pts.size() < 2) return;
                        this->draw_line(pts[0][0], pts[0][1], pts[1][0], pts[1][1], colour, stroke_width);
                    };
                    
                    draw_ln(l1.d,l1.e,l1.f);
                    draw_ln(l2.d,l2.e,l2.f);
                }
            }
            else
            {
                
                const auto curves = sample_points_conics(con, 100000, MIN_DIST, MAX_RADIUS);
                for (const auto & curve: curves)
                {
                    /*
                    for (uint32_t i = 1; i < curve.size(); ++i)
                    {
                        const auto & [x0, y0] = curve[i-1];
                        const auto & [x1, y1] = curve[i];
                        draw_line(x0,y0,x1,y1,colour, stroke_width);
                    }
                    */

                    file << "<polyline points=\"";
                    
                    for (const auto & [x,y]: curve)
                    {
                        file << x*box_radius/2. << "," << -y*box_radius/2. << " ";
                    }
                    file << "\" stroke=\"#" << std::hex <<  std::hex << std::setfill('0') << std::setw(6) <<colour << std::dec
                 << "\" stroke-width=\""<<stroke_width <<"\" fill=\"none\" />\n";



                    
                }
            }
             
            return *this;
        }

        SVGRenderer &
        draw_marker(const double x, const double y, const uint32_t colour)
        {
            const double r = 0.075;
            const double stroke_width = box_radius/100.;
            draw_line(x-r, y, x+r, y, colour, stroke_width);
            draw_line(x, y-r, x, y+r, colour, stroke_width);
            return *this;
        }
        
    };


    

};
