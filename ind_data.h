#ifndef IND_DATA_H
#define IND_DATA_H
#include <vector>
#include <json.hpp>
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
bool operator== (const Ind_Data<Ind>& lhs, const Ind_Data<Ind>& rhs);

#endif // IND_DATA_H
