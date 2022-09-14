#ifndef EVALUATION_TYPE_H
#define EVALUATION_TYPE_H
#include<string>
#include<map>

enum class evaluation_type
{
    trial,
    full_rn
};

static std::map<std::string, evaluation_type> string_to_eval_type_map
{
    {"trial", evaluation_type::trial},
    {"full_rn", evaluation_type::full_rn},
};

std::string convert_eval_type_to_string(evaluation_type e);

#endif // EVALUATION_TYPE_H
