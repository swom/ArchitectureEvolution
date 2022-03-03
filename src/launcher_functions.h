#ifndef LAUNCHER_FUNCTIONS_H
#define LAUNCHER_FUNCTIONS_H
#include"parser.h"


template <env_change_symmetry_type env_change_sym_type, env_change_freq_type env_change_freq_type, selection_type selection_type>
void run_simulation_given_mut_type(const cxxopts::ParseResult& results)
{
    auto mut_type = convert_ind_args(results).m_mutation_type;

    if(mut_type == mutation_type::weights)
    {
        using net_t = network<mutation_type::weights>;
        using ind_t = individual<net_t>;
        using pop_t = population<ind_t>;
        using sim_t = simulation<pop_t, env_change_sym_type, env_change_freq_type, selection_type>;
        observer<sim_t> o;

        auto s = create_simulation<pop_t, env_change_sym_type, env_change_freq_type, selection_type>(results);
        exec<sim_t>(s, o) ;
        save_json(o,
                  create_save_name_from_observer_data(o));
    }
    else if (mut_type == mutation_type::activation) {

        using net_t = network<mutation_type::activation>;
        using ind_t = individual<net_t>;
        using pop_t = population<ind_t>;
        using sim_t = simulation<pop_t, env_change_sym_type, env_change_freq_type, selection_type>;

        observer<sim_t> o;
        auto s = create_simulation<pop_t, env_change_sym_type, env_change_freq_type, selection_type>(results);
        exec<sim_t>(s, o) ;
        save_json(o,
                  create_save_name_from_observer_data(o));
    }
    else if (mut_type == mutation_type::weights_and_activation) {

        using net_t = network<mutation_type::weights_and_activation>;
        using ind_t = individual<net_t>;
        using pop_t = population<ind_t>;
        using sim_t = simulation<pop_t, env_change_sym_type, env_change_freq_type, selection_type>;

        observer<sim_t> o;
        auto s = create_simulation<pop_t, env_change_sym_type, env_change_freq_type, selection_type>(results);
        exec<sim_t>(s, o) ;
        save_json(o,
                  create_save_name_from_observer_data(o));
    }
    else if (mut_type == mutation_type::duplication) {

        using net_t = network<mutation_type::duplication>;
        using ind_t = individual<net_t>;
        using pop_t = population<ind_t>;
        using sim_t = simulation<pop_t, env_change_sym_type, env_change_freq_type, selection_type>;

        observer<sim_t> o;
        auto s = create_simulation<pop_t, env_change_sym_type, env_change_freq_type, selection_type>(results);
        exec<sim_t>(s, o) ;
        save_json(o,
                  create_save_name_from_observer_data(o));
    }
    else if (mut_type == mutation_type::NRduplication) {

        using net_t = network<mutation_type::NRduplication>;
        using ind_t = individual<net_t>;
        using pop_t = population<ind_t>;
        using sim_t = simulation<pop_t, env_change_sym_type, env_change_freq_type, selection_type>;

        observer<sim_t> o;
        auto s = create_simulation<pop_t, env_change_sym_type, env_change_freq_type, selection_type>(results);
        exec<sim_t>(s, o) ;
        save_json(o,
                  create_save_name_from_observer_data(o));
    }
    else if (mut_type == mutation_type::addition) {

        using net_t = network<mutation_type::addition>;
        using ind_t = individual<net_t>;
        using pop_t = population<ind_t>;
        using sim_t = simulation<pop_t, env_change_sym_type, env_change_freq_type, selection_type>;

        observer<sim_t> o;
        auto s = create_simulation<pop_t, env_change_sym_type, env_change_freq_type, selection_type>(results);
        exec<sim_t>(s, o) ;
        save_json(o,
                  create_save_name_from_observer_data(o));
    }
    else if (mut_type == mutation_type::NRaddition) {

        using net_t = network<mutation_type::NRaddition>;
        using ind_t = individual<net_t>;
        using pop_t = population<ind_t>;
        using sim_t = simulation<pop_t, env_change_sym_type, env_change_freq_type, selection_type>;

        observer<sim_t> o;
        auto s = create_simulation<pop_t, env_change_sym_type, env_change_freq_type, selection_type>(results);
        exec<sim_t>(s, o) ;
        save_json(o,
                  create_save_name_from_observer_data(o));
    }
    else
    {
        throw std::runtime_error{"unknown mutation type"};
    }
}

template<env_change_freq_type E, selection_type S>
void run_simulation_given_env_change_type(const cxxopts::ParseResult& results)
{
    auto env_change_symmetry_type = convert_sim_args(results).change_sym_type ;
    if(env_change_symmetry_type == env_change_symmetry_type::asymmetrical)
    {
        run_simulation_given_mut_type<env_change_symmetry_type::asymmetrical, E, S>(results);
    }
    else if(env_change_symmetry_type == env_change_symmetry_type::symmetrical)
    {
        run_simulation_given_mut_type<env_change_symmetry_type::symmetrical, E, S>(results);
    }
    else{
        throw std::runtime_error{"unknown change type"};
    }
}

template<env_change_freq_type E>
void run_simulation_given_sel_type(const cxxopts::ParseResult& results)
{
    auto sel_type = convert_sim_args(results).sel_type ;
    if(sel_type == selection_type::constant)
    {
        run_simulation_given_env_change_type<E,selection_type::constant>(results);
    }
    else if(sel_type == selection_type::sporadic)
    {
        run_simulation_given_env_change_type<E, selection_type::sporadic>(results);
    }
    else{
        throw std::runtime_error{"unknown change type"};
    }
}

void run_simulation_given_env_freq_type(const cxxopts::ParseResult& results)
{
    auto env_change_freq_type = convert_sim_args(results).change_freq_type ;
    if(env_change_freq_type == env_change_freq_type::regular)
    {
        run_simulation_given_sel_type<env_change_freq_type::regular>(results);
    }
    else if(env_change_freq_type == env_change_freq_type::stochastic)
    {
        run_simulation_given_sel_type<env_change_freq_type::stochastic>(results);
    }
    else{
        throw std::runtime_error{"unknown change type"};
    }
}


#endif // LAUNCHER_FUNCTIONS_H
