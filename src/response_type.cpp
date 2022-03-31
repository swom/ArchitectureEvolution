#include "response_type.h"
#include <stdexcept>

std::string convert_response_type_to_string(response_type a)
{
    std::string string;

    switch (a) {
    case response_type::plastic :
        string = "plastic";
        return string;
        break;

    case response_type::constitutive :
        string = "constitutive";
        return string;
        break;

    default:
        throw std::runtime_error{"could not convert response_type type into string"};
        return "failed";
    }
}

