#include <stdexcept>

#include "adaptation_period.h"

std::string convert_adapt_periods_to_string(adaptation_period a)
{
    std::string string;

    switch (a) {
    case adaptation_period::on :
        string = "on";
        return string;
        break;

    case adaptation_period::off :
        string = "off";
        return string;
        break;

    default:
        throw std::runtime_error{"could not convert adptation_period type into string"};
        return "failed";
    }
}
