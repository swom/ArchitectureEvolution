#ifndef ADAPTATION_PERIOD_H
#define ADAPTATION_PERIOD_H
#include<string>
#include<map>

enum class adaptation_period
{
    on,
    off
};

static std::map<std::string, adaptation_period> string_to_adapt_period_map
{
    {"on", adaptation_period::on},
    {"off", adaptation_period::off},
};

std::string convert_adapt_perios_to_string(adaptation_period e);

#endif // ADAPTATION_PERIOD_H
