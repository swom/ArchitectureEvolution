#ifndef RANGE_H
#define RANGE_H
#include <vector>
#include "json.hpp"

struct range
{
    range(double start, double end);

    NLOHMANN_DEFINE_TYPE_INTRUSIVE(range,
                                   m_start,
                                   m_end);
    double m_start;
    double m_end;
};

#endif // RANGE_H
