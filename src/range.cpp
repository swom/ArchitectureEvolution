#include <iostream>
#include "range.h"

bool operator==(const range& lhs, const range& rhs)
{
 bool starts = lhs.m_start == rhs.m_start;
 bool ends = lhs.m_end == rhs.m_end;
 return starts && ends;
}

range::range(double start, double end):
    m_start{start},
    m_end{end}
{
 if(m_start > m_end)
 {
     throw std::runtime_error{"invalid cue range"};
 }
}
