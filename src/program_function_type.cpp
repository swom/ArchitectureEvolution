#include "program_function_type.h"
#include <stdexcept>

std::string convert_eval_type_to_string(const program_function_type& p)
{
    std::string string;

    switch (p) {
    case program_function_type::simulation :
        string = "simulation";
        return string;
        break;

    case program_function_type::mutational_spectrum_calculation :
        string = "mutational_spectrum_calculation";
        return string;
        break;

    default:
        throw std::runtime_error{"could not convert program_function type into string"};
        return "failed";
    }
}
