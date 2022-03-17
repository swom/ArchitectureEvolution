#ifndef ADAPTATION_PERIOD_H
#define ADAPTATION_PERIOD_H
#include<string>

enum class adaptation_period
{
    on,
    off
};

std::string convert_adapt_perios_to_string(adaptation_period e);

#endif // ADAPTATION_PERIOD_H
