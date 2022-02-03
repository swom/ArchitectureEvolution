#ifndef SELECTION_TYPE_H
#define SELECTION_TYPE_H
#include<map>

enum class selection_type
{
    sporadic,
    constant
};

std::string convert_selection_type_to_string(selection_type s);


#endif // SELECTION_TYPE_H
