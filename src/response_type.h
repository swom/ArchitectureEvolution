#ifndef RESPONSE_TYPE_H
#define RESPONSE_TYPE_H
#include <map>
#include<string>

enum class response_type
{
    plastic,
    constitutive,
    additive
};


static std::map<std::string, response_type> string_to_response_type_map
{
    {"constitutive", response_type::constitutive},
    {"plastic", response_type::plastic},
    {"additive", response_type::additive}
};

std::string convert_response_type_to_string(response_type a);

#endif // RESPONSE_TYPE_H
