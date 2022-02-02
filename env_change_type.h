#ifndef ENV_CHANGE_TYPE_H
#define ENV_CHANGE_TYPE_H
#include <map>

enum class env_change_type
{
    symmetrical,
    asymmetrical
};

static std::map<std::string, env_change_type> string_to_env_change_map
{
    {"symmetrical", env_change_type::symmetrical},
    {"asymmetrical", env_change_type::asymmetrical}
};

#endif // ENV_CHANGE_TYPE_H
