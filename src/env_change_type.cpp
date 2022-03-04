#include <stdexcept>
#include "env_change_type.h"

std::string convert_change_symmetry_type_to_string(env_change_symmetry_type e)
{
    std::string string;

    switch (e) {
    case env_change_symmetry_type::symmetrical :
        string = "symmetrical";
        return string;
        break;

    case env_change_symmetry_type::asymmetrical :
        string = "asymmetrical";
        return string;
        break;

    default:
        throw std::runtime_error{"could not convert env_change type into string"};
        return "failed";
    }
}

std::string convert_change_freq_type_to_string(env_change_freq_type e)
{
    std::string string;

    switch (e) {
    case env_change_freq_type::regular :
        string = "regular";
        return string;
        break;

    case env_change_freq_type::stochastic :
        string = "regular";
        return string;
        break;

    default:
        throw std::runtime_error{"could not convert env_change type into string"};
        return "failed";
    }
}
