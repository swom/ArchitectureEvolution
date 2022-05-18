#include "netwrok_spectrum.h"
#include "utilities.h"

bool operator== (const fit_and_phen_sens_t& lhs, const fit_and_phen_sens_t& rhs)
{
    auto phenotype = are_equal_with_tolerance(lhs.m_phenotype, rhs.m_phenotype);
    auto fitness = are_equal_with_tolerance(lhs.m_fitness, rhs.m_fitness);

    return phenotype && fitness;
}

bool operator!= (const fit_and_phen_sens_t& lhs, const fit_and_phen_sens_t& rhs)
{
    return !(lhs == rhs);
}

bool operator==(const react_norm_t& lhs, const react_norm_t& rhs)
{
    auto xs = lhs.m_x == rhs.m_x;
    auto ys = are_equal_with_tolerance(lhs.m_y, rhs.m_y);
    return xs && ys;
}

#ifndef NDEBUG
void test_spectrum()
{

}
#endif
