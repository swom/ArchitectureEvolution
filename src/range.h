#ifndef RANGE_H
#define RANGE_H
#include <vector>
#include <nlohmann/json.hpp>


struct range
{
    range(double start, double end);

    NLOHMANN_DEFINE_TYPE_INTRUSIVE(range,
                                   m_start,
                                   m_end);
    double m_start;
    double m_end;
};

bool operator==(const range& lhs, const range& rhs);
#endif // RANGE_H
