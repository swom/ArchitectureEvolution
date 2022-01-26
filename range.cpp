#include <iostream>
#include "range.h"


range::range(double start, double end):
    m_start{start},
    m_end{end}
{
 if(m_start < m_end)
 {
     throw std::runtime_error{"invalid cue range"};
 }
}
