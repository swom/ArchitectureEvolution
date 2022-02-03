#ifndef LAUNCH_ENV_CHANGE_H
#define LAUNCH_ENV_CHANGE_H
#include"launch_mut_type.h"
template<selection_type S>
void run_simulation_given_env_change_type(const cxxopts::ParseResult& results)
{
    auto env_change_type = convert_sim_args(results).change_type ;
    if(env_change_type == env_change_type::asymmetrical)
    {
        run_simulation_given_mut_type<env_change_type::asymmetrical, S>(results);
    }
    else if(env_change_type == env_change_type::symmetrical)
    {
        run_simulation_given_mut_type<env_change_type::symmetrical, S>(results);
    }
    else if(env_change_type == env_change_type::regular)
    {
        run_simulation_given_mut_type<env_change_type::regular, S>(results);
    }
    else{
        throw std::runtime_error{"unknown change type"};
    }
}

#endif // LAUNCH_ENV_CHANGE_H
