#ifndef PROT_ALGORITHM_HPP_
#define PROT_ALGORITHM_HPP_

#include <cstddef>
#include <vector>

namespace prot {

/* used in find matched mass pairs */
bool increaseIJ(size_t i, size_t j, double deviation, double tolerance, 
                const std::vector<double> &ms_masses, 
                const std::vector<double> &theo_masses);

/* compute ppos for ms_masses */
std::vector<double> compMsMassPpos(const std::vector<double> &ms_masses, 
                                   const std::vector<double> &theo_masses, 
                                   double ppo);

/* compute the number of matched theoretical masses (fragment ions) */
double compNumMatchedTheoMasses (const std::vector<double> &ms_masses, 
                                 const std::vector<double> &theo_masses, 
                                 double ppo); 

/* compute the position of the last residue of a proteoform based 
 * on its n term shift */
int getFirstResPos(double n_term_shift, const std::vector<double> &prm_masses);


/* compute the position of the last residue of a proteoform based 
 * on its c term shift */
int getLastResPos(double c_term_shift, const std::vector<double> &prm_masses);

}

#endif
