#include "netwrok_spectrum.h"
bool operator==(const react_norm_t& lhs, const react_norm_t& rhs)
{
    auto xs = lhs.x == rhs.x;
    auto ys = lhs.y == rhs.y;
    return xs && ys;
}
