#ifndef LAUNCHER_FUNCTIONS_H
#define LAUNCHER_FUNCTIONS_H
#include"launch_env_change.h"
void run_simulation_given_sel_type(const cxxopts::ParseResult& results)
{
    auto sel_type = convert_sim_args(results).selection_type ;
    if(sel_type == selection_type::constant)
    {
        run_simulation_given_env_change_type<selection_type::constant>(results);
    }
    else if(sel_type == selection_type::sporadic)
    {
        run_simulation_given_env_change_type<selection_type::sporadic>(results);
    }
    else{
        throw std::runtime_error{"unknown change type"};
    }
}

#endif // LAUNCHER_FUNCTIONS_H
