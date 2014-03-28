#ifndef PROT_ALGORITHM_HPP_
#define PROT_ALGORITHM_HPP_

#include <vector>

namespace prot {

/* used in find matched mass pairs */
bool increaseIJ(unsigned int i, unsigned int j, double deviation, 
                double tolerance, std::vector<double> ms_masses, 
                std::vector<double> theo_masses);

/* compute ppos for ms_masses */
void compMsMassPpos(std::vector<double> &ms_masses, 
                    std::vector<double> &theo_masses, 
                    double ppo, std::vector<double> &result_ppos);

/* compute the number of matched theoretical masses (fragment ions) */
double compNumMatchedTheoMasses (std::vector<double> &ms_masses, 
                                 std::vector<double> &theo_masses, 
                                 double ppo); 

/* compute the position of the last residue of a proteoform based 
 * on its n term shift */
int getFirstResPos(double n_term_shift,std::vector<double> prm_masses);


/* compute the position of the last residue of a proteoform based 
 * on its c term shift */
int getLastResPos(double c_term_shift,std::vector<double> prm_masses);

}

#endif
