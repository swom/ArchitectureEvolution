#include <stdexcept>
#include "selection_type.h"

std::string convert_selection_type_to_string(selection_type s)
{
    std::string string;

    switch (s) {
    case selection_type::constant :
        string = "constant";
        return string;
        break;

    case selection_type::sporadic :
        string = "sporadic";
        return string;
        break;

    default:
        throw std::runtime_error{"could not convert selection type into string"};
    }
}
