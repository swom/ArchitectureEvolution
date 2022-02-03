#include "env_change_type.h"

std::string convert_change_type_to_string(env_change_type e)
{
    std::string string;

    switch (e) {
    case env_change_type::symmetrical :
        string = "symmetrical";
        return string;
        break;

    case env_change_type::asymmetrical :
        string = "asymmetrical";
        return string;
        break;

    case env_change_type::regular :
        string = "regular";
        return string;
        break;
    default:
        throw std::runtime_error{"could not convert env_change type into string"};
    }
}
