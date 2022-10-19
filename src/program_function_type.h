#ifndef PROGRAM_FUNCTION_TYPE_H
#define PROGRAM_FUNCTION_TYPE_H
#include<string>
#include<map>

enum class program_function_type
{
    simulation,
    mutational_spectrum_calculation
};

static std::map<std::string, program_function_type> string_to_program_function_type_map
{
    {"simulation", program_function_type::simulation},
    {"mutational_spectrum_calculation", program_function_type::mutational_spectrum_calculation},
};

std::string convert_program_function_type_to_string(program_function_type p);

#endif // PROGRAM_FUNCTION_TYPE_H
