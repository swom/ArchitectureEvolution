#include "evaluation_type.h"
#include <stdexcept>

std::string convert_eval_type_to_string(evaluation_type e)
{
    std::string string;

    switch (e) {
    case evaluation_type::trial :
        string = "trial";
        return string;
        break;

    case evaluation_type::full_rn :
        string = "full_rn";
        return string;
        break;

    default:
        throw std::runtime_error{"could not convert adptation_period type into string"};
        return "failed";
    }
}
