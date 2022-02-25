#ifndef IND_DATA_H
#define IND_DATA_H
#include <vector>
#include <json.hpp>
#include "netwrok_spectrum.h"

template<class Ind>
struct Ind_Data
{
    NLOHMANN_DEFINE_TYPE_INTRUSIVE(Ind_Data,
                                   m_ind,
                                   m_reac_norm)
    Ind m_ind;
    reac_norm m_reac_norm;
};

template<class Ind>
bool operator== (const Ind_Data<Ind>& lhs, const Ind_Data<Ind>& rhs)
{
    auto ind = lhs.m_ind == rhs.m_ind;
    auto reaction_norm = lhs.m_reac_norm == rhs.m_reac_norm;
    return ind && reaction_norm;
}



template<class Ind>
struct Ind_Spectrum
{
    NLOHMANN_DEFINE_TYPE_INTRUSIVE(Ind_Spectrum,
                                   m_ind,
                                   reac_norm)
    Ind m_ind;
    reac_norm reac_norm;
    network_spectrum net_spectrum;
};

///Calculates the reaction_norm of an individual's network
/// for a given range and a given number of data points
template<class Ind>
std::vector<Ind_Data<Ind>> calculate_reaction_norms(const std::vector<Ind>& inds,
                                                     const range& cue_range,
                                                     const int& n_data_points)
{

    std::vector<Ind_Data<Ind>> inds_data(inds.size());
    for(const auto& ind : inds)
    {
        auto r_norm =calculate_reaction_norm(ind.get_net(),
                                             cue_range,
                                             n_data_points);
        inds_data.push_back({ind, r_norm});
    }
    return inds_data;
}

#endif // IND_DATA_H
