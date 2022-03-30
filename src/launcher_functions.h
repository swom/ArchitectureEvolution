#ifndef LAUNCHER_FUNCTIONS_H
#define LAUNCHER_FUNCTIONS_H
#include"parser.h"


template <env_change_symmetry_type S,
          env_change_freq_type F,
          selection_type Sel,
          adaptation_period A>
void run_simulation_given_mut_type(const cxxopts::ParseResult& results)
{
    auto mut_type = convert_ind_args(results).m_mutation_type;

    if(mut_type == mutation_type::weights)
    {
        using net_t = network<mutation_type::weights>;
        using ind_t = individual<net_t>;
        using pop_t = population<ind_t>;
        using sim_t = simulation<pop_t, S, F, Sel, A>;

        auto s = create_simulation<pop_t, S, F, Sel, A>(results);
        observer<sim_t> o{convert_obs_args(results), s.get_params()};
        exec<sim_t>(s, o) ;
        save_json(o,
                  create_save_name_from_params(o.get_params()));
    }
    else if (mut_type == mutation_type::activation) {

        using net_t = network<mutation_type::activation>;
        using ind_t = individual<net_t>;
        using pop_t = population<ind_t>;
        using sim_t = simulation<pop_t, S, F, Sel, A>;

        auto s = create_simulation<pop_t, S, F, Sel, A>(results);
        observer<sim_t> o{convert_obs_args(results), s.get_params()};
        exec<sim_t>(s, o) ;
        save_json(o,
                  create_save_name_from_params(o.get_params()));
    }
    else if (mut_type == mutation_type::weights_and_activation) {

        using net_t = network<mutation_type::weights_and_activation>;
        using ind_t = individual<net_t>;
        using pop_t = population<ind_t>;
        using sim_t = simulation<pop_t, S, F, Sel, A>;

        auto s = create_simulation<pop_t, S, F, Sel, A>(results);
        observer<sim_t> o{convert_obs_args(results), s.get_params()};
        exec<sim_t>(s, o) ;
        save_json(o,
                  create_save_name_from_params(o.get_params()));
    }
    else if (mut_type == mutation_type::duplication) {

        using net_t = network<mutation_type::duplication>;
        using ind_t = individual<net_t>;
        using pop_t = population<ind_t>;
        using sim_t = simulation<pop_t, S, F, Sel, A>;

        auto s = create_simulation<pop_t, S, F, Sel, A>(results);
        observer<sim_t> o{convert_obs_args(results), s.get_params()};
        exec<sim_t>(s, o) ;
        save_json(o,
                  create_save_name_from_params(o.get_params()));
    }
    else if (mut_type == mutation_type::NRduplication) {

        using net_t = network<mutation_type::NRduplication>;
        using ind_t = individual<net_t>;
        using pop_t = population<ind_t>;
        using sim_t = simulation<pop_t, S, F, Sel, A>;

        auto s = create_simulation<pop_t, S, F, Sel, A>(results);
        observer<sim_t> o{convert_obs_args(results), s.get_params()};
        exec<sim_t>(s, o) ;
        save_json(o,
                  create_save_name_from_params(o.get_params()));
    }
    else if (mut_type == mutation_type::addition) {

        using net_t = network<mutation_type::addition>;
        using ind_t = individual<net_t>;
        using pop_t = population<ind_t>;
        using sim_t = simulation<pop_t, S, F, Sel, A>;

        auto s = create_simulation<pop_t, S, F, Sel, A>(results);
        observer<sim_t> o{convert_obs_args(results), s.get_params()};
        exec<sim_t>(s, o) ;
        save_json(o,
                  create_save_name_from_params(o.get_params()));
    }
    else if (mut_type == mutation_type::NRaddition) {

        using net_t = network<mutation_type::NRaddition>;
        using ind_t = individual<net_t>;
        using pop_t = population<ind_t>;
        using sim_t = simulation<pop_t, S, F, Sel, A>;

        auto s = create_simulation<pop_t, S, F, Sel, A>(results);
        observer<sim_t> o{convert_obs_args(results), s.get_params()};
        exec<sim_t>(s, o) ;
        save_json(o,
                  create_save_name_from_params(o.get_params()));
    }
    else
    {
        throw std::runtime_error{"unknown mutation type"};
    }
}

template<env_change_freq_type E,
         selection_type S,
         adaptation_period A>
void run_simulation_given_eptype(const cxxopts::ParseResult& results)
{
    auto env_change_symmetry_type = convert_sim_args(results).change_sym_type ;
    if(env_change_symmetry_type == env_change_symmetry_type::asymmetrical)
    {
        run_simulation_given_mut_type<env_change_symmetry_type::asymmetrical, E, S, A>(results);
    }
    else if(env_change_symmetry_type == env_change_symmetry_type::symmetrical)
    {
        run_simulation_given_mut_type<env_change_symmetry_type::symmetrical, E, S, A>(results);
    }
    else{
        throw std::runtime_error{"unknown change type"};
    }
}

template<env_change_freq_type E,
         selection_type S,
         adaptation_period A>
void run_simulation_given_env_change_type(const cxxopts::ParseResult& results)
{
    auto env_change_symmetry_type = convert_sim_args(results).change_sym_type ;
    if(env_change_symmetry_type == env_change_symmetry_type::asymmetrical)
    {
        run_simulation_given_mut_type<env_change_symmetry_type::asymmetrical, E, S, A>(results);
    }
    else if(env_change_symmetry_type == env_change_symmetry_type::symmetrical)
    {
        run_simulation_given_mut_type<env_change_symmetry_type::symmetrical, E, S, A>(results);
    }
    else{
        throw std::runtime_error{"unknown change type"};
    }
}

template<env_change_freq_type E,
         adaptation_period A>
void run_simulation_given_sel_type(const cxxopts::ParseResult& results)
{
    auto sel_type = convert_sim_args(results).sel_type ;
    if(sel_type == selection_type::constant)
    {
        run_simulation_given_env_change_type<E,selection_type::constant,A>(results);
    }
    else if(sel_type == selection_type::sporadic)
    {
        run_simulation_given_env_change_type<E, selection_type::sporadic, A>(results);
    }
    else{
        throw std::runtime_error{"unknown change type"};
    }
}

template<adaptation_period A>
void run_simulation_given_env_freq_type(const cxxopts::ParseResult& results)
{
    auto env_change_freq_type = convert_sim_args(results).change_freq_type ;
    if(env_change_freq_type == env_change_freq_type::regular)
    {
        run_simulation_given_sel_type<env_change_freq_type::regular,A>(results);
    }
    else if(env_change_freq_type == env_change_freq_type::stochastic)
    {
        run_simulation_given_sel_type<env_change_freq_type::stochastic, A>(results);
    }
    else{
        throw std::runtime_error{"unknown change type"};
    }
}

void run_simulation_given_adaptation_period(const cxxopts::ParseResult& results)
{
    auto adaptation_period = convert_sim_args(results).adaptation_per ;
    if(adaptation_period == adaptation_period::on)
    {
        run_simulation_given_env_freq_type<adaptation_period::on>(results);
    }
    else if(adaptation_period == adaptation_period::off)
    {
        run_simulation_given_env_freq_type<adaptation_period::off>(results);
    }
    else{
        throw std::runtime_error{"unknown adaptation period command"};
    }
}
#endif // LAUNCHER_FUNCTIONS_H
