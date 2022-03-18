#ifndef ENV_CHANGE_TYPE_H
#define ENV_CHANGE_TYPE_H
#include <map>
#include<string>
enum class env_change_symmetry_type
{
    symmetrical,
    asymmetrical,
};

enum class env_change_freq_type
{
    stochastic,
    regular
};

static std::map<std::string, env_change_symmetry_type> string_to_env_change_symmetry_type_map
{
    {"symmetrical", env_change_symmetry_type::symmetrical},
    {"asymmetrical", env_change_symmetry_type::asymmetrical},
};

static std::map<std::string, env_change_freq_type> string_to_env_change_freq_type_map
{
    {"regular", env_change_freq_type::regular},
    {"stochastic", env_change_freq_type::stochastic},
};

std::string convert_change_symmetry_type_to_string(env_change_symmetry_type e);
std::string convert_change_freq_type_to_string(env_change_freq_type e);


#endif // ENV_CHANGE_TYPE_H
