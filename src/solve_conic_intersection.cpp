#include "conic.hpp"
#include "intersections.hpp"
#include "random.hpp"
#include <complex>
#include <fstream>
#include <sstream>
#include <stdexcept>

using S = double;

qc2ga::Mvec<S>
normalize_mvec(const qc2ga::Mvec<S> & m)
{
    S max_val = 0.0f;
    for (uint32_t i = 0; i < (1<<8); ++i)
    {
        const S val = std::abs(m[i]);
        if (max_val < val) max_val = val;
    }
    return m/max_val;
}






void
test_lc(RNG<S> & gen,
        const std::vector<qc2ga::Mvec<S>> & inters,
        const uint32_t n_inter_pts,
        const bool display = false)
{
    uint32_t success = 0;
    std::vector<uint32_t> histo_lc(5);
    for (const qc2ga::Mvec<S> & inter: inters)
    {
        const auto & [p1,p2] = gen.point_array<2>();
        qc2ga::Mvec<S> qC_1 = inter ^ p1;
        qc2ga::Mvec<S> qC_2 = inter ^ p2;
        const Conic C_1(!qC_1);
        const Conic C_2(!qC_2);
        const std::vector<std::array<S,2>> pts_from_lc = Algo(inter, gen).run().inter_pts;
        bool penalize_wrong_number = true;
        bool failed = pts_from_lc.size() < n_inter_pts && penalize_wrong_number;
        if (!failed)
        for (const auto &[x,y]: pts_from_lc)
        {
            const S eval_pt = qcga_tools::norm_l2(qc2ga::point<S>(x,y) ^ inter);
            if (eval_pt > 0.001 || std::isnan(x) || std::isnan(y))
            { // failure
                failed = true;
                break;
            }
        }
        if (failed) histo_lc[pts_from_lc.size()]++;
        else success ++;
        
    }
    std::cout << "SUCCESS RATE: " << double(success)/double(inters.size())*100. << std::endl;
    std::cout << "HISTO OF NB WRONG POINTS: (%)";
    double sum = {};
    for (const uint32_t bin: histo_lc)
    {
        const double val = double(bin)/double(inters.size())*100.;
        std::cout << " | " << val;
        sum+=val;
    } std::cout << "|" << std::endl;
    std::cout << "# TOTAL ERROR = " << sum << std::endl;
    std::cout << std::endl;
}


std::vector<qc2ga::Mvec<double>>
open_file(const std::string & filename)
{
    std::string line;
    std::ifstream file(filename);

    
    if (!file.is_open())
    {
        throw std::ios_base::failure("Could not open file");
        return {};
    }

    double coef;
    std::vector<qc2ga::Mvec<double>> inters;
    while ( getline (file,line) )
    {
        std::stringstream ss(line);

        qc2ga::Mvec<double> dinter = {};
        for (uint32_t i = 0; i < 7; ++i)
        {
            for (uint32_t j = i+1; j < 8; ++j)
            {
                const uint32_t index = (1 << i) | (1 << j);
                ss >> coef;
                dinter[index] = coef;
            }
        }
        inters.push_back(dinter.dual());
    }
    file.close();
    return inters;
}








int
main(const int argc, const char ** argv)
{
    if (argc <= 1)
    {
        std::cerr << "Wrong number of arguments!\n";
        std::cout << "Usage: solve intersections.txt\n";
        return -1;
    }
    
    const auto inters = open_file(argv[1]);


    int seed = 4317732;
    RNG<S> gen(seed);
        
        
  
    test_lc(gen, inters, 0, false);
    
    return 0;
    
}
