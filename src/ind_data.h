#ifndef IND_DATA_H
#define IND_DATA_H
#include <vector>
#include <nlohmann/json.hpp>

template<class Ind>
struct Ind_Data
{
    NLOHMANN_DEFINE_TYPE_INTRUSIVE(Ind_Data,
                                   m_ind,
                                   reac_norm)
    Ind m_ind;
    std::vector<std::vector<double>> reac_norm;
};

template<class Ind>
bool operator== (const Ind_Data<Ind>& lhs, const Ind_Data<Ind>& rhs)
{
    auto ind = lhs.m_ind == rhs.m_ind;
    auto reaction_norm = lhs.reac_norm == rhs.reac_norm;
    return ind && reaction_norm;
}

#endif // IND_DATA_H
