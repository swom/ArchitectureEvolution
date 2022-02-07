#ifndef ENV_CHANGE_TYPE_H
#define ENV_CHANGE_TYPE_H
#include <map>
#include<string>

enum class env_change_type
{
    symmetrical,
    asymmetrical,
    regular
};

static std::map<std::string, env_change_type> string_to_env_change_map
{
    {"symmetrical", env_change_type::symmetrical},
    {"asymmetrical", env_change_type::asymmetrical},
    {"regular", env_change_type::regular}
};

std::string convert_change_type_to_string(env_change_type e);


#endif // ENV_CHANGE_TYPE_H
