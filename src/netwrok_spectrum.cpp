#include "netwrok_spectrum.h"
bool operator==(const react_norm_t& lhs, const react_norm_t& rhs)
{
    auto xs = lhs.m_x == rhs.m_x;
    auto ys = lhs.m_y == rhs.m_y;
    return xs && ys;
}

#ifndef NDEBUG
void test_spectrum()
{

}
#endif
