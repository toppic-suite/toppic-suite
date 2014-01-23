#ifndef PROT_ALGORITHM_HPP_
#define PROT_ALGORITHM_HPP_

#include <vector>

namespace prot {

bool increaseIJ(unsigned int i, unsigned int j, double deviation, 
                double tolerance, std::vector<double> ms_masses, 
                std::vector<double> theo_masses);

void compMsMassPpos(std::vector<double> &ms_masses, 
                    std::vector<double> &theo_masses, 
                    double ppo, std::vector<double> &result_ppos);

double compUniqueScore (std::vector<double> &ms_masses, 
                        std::vector<double> &theo_masses, 
                        double ppo); 
}

#endif
