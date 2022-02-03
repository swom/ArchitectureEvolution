#ifndef SELECTION_TYPE_H
#define SELECTION_TYPE_H
#include<map>

enum class selection_type
{
    sporadic,
    constant
};

static std::map<std::string, selection_type> string_to_sel_type_map
{
    {"constant", selection_type::constant},
    {"sporadic", selection_type::sporadic}
};

std::string convert_selection_type_to_string(selection_type s);


#endif // SELECTION_TYPE_H
