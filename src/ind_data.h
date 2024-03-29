#ifndef IND_DATA_H
#define IND_DATA_H
#include <vector>
#include <nlohmann/json.hpp>
#include "individual.h"

template<class Ind>
struct Ind_Spectrum
{
    using Net = typename Ind::net_t;
    using Net_Spect = network_spectrum<Net>;

    Ind_Spectrum(){};
    Ind_Spectrum(Ind i, Net_Spect net_s, int gen):
        m_ind{i},
        net_spectrum{net_s},
        generation{gen}
    {}

    NLOHMANN_DEFINE_TYPE_INTRUSIVE(Ind_Spectrum,
                                   m_ind,
                                   net_spectrum,
                                   generation)
    Ind m_ind;
    Net_Spect net_spectrum;
    int generation;
};
template<class Ind>
bool operator==(const Ind_Spectrum<Ind>& lhs, const Ind_Spectrum<Ind>& rhs)
{
    bool ind = lhs.m_ind == rhs.m_ind;
    bool spec = lhs.net_spectrum == rhs.net_spectrum;
    bool gen = lhs.generation == rhs.generation;
    return ind && spec && gen;
}


///Calculates the reaction_norm of individuals' networks
/// for a given range and a given number of data points
template<class Ind>
std::vector<Ind_Spectrum<Ind>> calculate_mut_spectrums(std::vector<Ind> inds,
                                                       double mut_step,
                                                       std::mt19937_64& rng,
                                                       int n_mutations,
                                                       const range& cue_range,
                                                       const int& n_data_points,
                                                       int generation)
{

    std::vector<Ind_Spectrum<Ind>> spect_vector(inds.size());
    for(auto ind = inds.begin(); ind != inds.end(); ind++)
    {

        spect_vector.at(std::distance(inds.begin(), ind)) =
                Ind_Spectrum<Ind>{
                *ind,
                ind->calc_spectrum(mut_step,
                                   rng,
                                   n_mutations,
                                   cue_range,
                                   n_data_points
                                   ),
                generation
    };
    }
    return spect_vector;
}
#endif // IND_DATA_H
