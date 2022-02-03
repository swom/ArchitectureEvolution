#ifndef IND_DATA_H
#define IND_DATA_H
#include <vector>

template<class Ind>
struct Ind_Data
{
    Ind m_ind;
    std::vector<double> reac_norm;
};

#endif // IND_DATA_H
