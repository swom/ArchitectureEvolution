#ifndef LAUNCHER_FUNCTIONS_H
#define LAUNCHER_FUNCTIONS_H
#include"parser.h"


template < response_type R,
           env_change_symmetry_type S,
           selection_type Sel,
           env_change_freq_type E,
           adaptation_period A,
           evaluation_type Ev>
void run_simulation_given_mut_type(const cxxopts::ParseResult& results)
{
    auto mut_type = convert_ind_args(results).m_mutation_type;

    if(mut_type == mutation_type::weights)
    {
        using net_t = network<mutation_type::weights, R>;
        using ind_t = individual<net_t>;
        using pop_t = population<ind_t>;
        using sim_t = simulation<pop_t, S, E, Sel, A, Ev>;

        auto s = create_simulation<pop_t, S, E, Sel, A, Ev>(results);
        observer<sim_t> o{convert_obs_args(results), s.get_params()};
        exec<sim_t>(s, o) ;

        std::cout <<"sel_str = " << o.get_params().s_p.selection_strength << std::endl;

        save_json(o,
                  create_save_name_from_params(o.get_params()));
    }
    else if (mut_type == mutation_type::activation) {

        using net_t = network<mutation_type::activation, R>;
        using ind_t = individual<net_t>;
        using pop_t = population<ind_t>;
        using sim_t = simulation<pop_t, S, E, Sel, A, Ev>;

        auto s = create_simulation<pop_t, S, E, Sel, A, Ev>(results);
        observer<sim_t> o{convert_obs_args(results), s.get_params()};
        exec<sim_t>(s, o) ;
        save_json(o,
                  create_save_name_from_params(o.get_params()));
    }
    else if (mut_type == mutation_type::weights_and_activation) {

        using net_t = network<mutation_type::weights_and_activation, R>;
        using ind_t = individual<net_t>;
        using pop_t = population<ind_t>;
        using sim_t = simulation<pop_t, S, E, Sel, A, Ev>;

        auto s = create_simulation<pop_t, S, E, Sel, A, Ev>(results);
        observer<sim_t> o{convert_obs_args(results), s.get_params()};
        exec<sim_t>(s, o) ;
        save_json(o,
                  create_save_name_from_params(o.get_params()));
    }
    else if (mut_type == mutation_type::duplication) {

        using net_t = network<mutation_type::duplication, R>;
        using ind_t = individual<net_t>;
        using pop_t = population<ind_t>;
        using sim_t = simulation<pop_t, S, E, Sel, A, Ev>;

        auto s = create_simulation<pop_t, S, E, Sel, A, Ev>(results);
        observer<sim_t> o{convert_obs_args(results), s.get_params()};
        exec<sim_t>(s, o) ;
        save_json(o,
                  create_save_name_from_params(o.get_params()));
    }
    else if (mut_type == mutation_type::NRduplication) {

        using net_t = network<mutation_type::NRduplication, R>;
        using ind_t = individual<net_t>;
        using pop_t = population<ind_t>;
        using sim_t = simulation<pop_t, S, E, Sel, A, Ev>;

        auto s = create_simulation<pop_t, S, E, Sel, A, Ev>(results);
        observer<sim_t> o{convert_obs_args(results), s.get_params()};
        exec<sim_t>(s, o) ;
        save_json(o,
                  create_save_name_from_params(o.get_params()));
    }
    else if (mut_type == mutation_type::addition) {

        using net_t = network<mutation_type::addition, R>;
        using ind_t = individual<net_t>;
        using pop_t = population<ind_t>;
        using sim_t = simulation<pop_t, S, E, Sel, A, Ev>;

        auto s = create_simulation<pop_t, S, E, Sel, A, Ev>(results);
        observer<sim_t> o{convert_obs_args(results), s.get_params()};
        exec<sim_t>(s, o) ;
        save_json(o,
                  create_save_name_from_params(o.get_params()));
    }
    else if (mut_type == mutation_type::NRaddition) {

        using net_t = network<mutation_type::NRaddition, R>;
        using ind_t = individual<net_t>;
        using pop_t = population<ind_t>;
        using sim_t = simulation<pop_t, S, E, Sel, A, Ev>;

        auto s = create_simulation<pop_t, S, E, Sel, A, Ev>(results);
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

template <env_change_symmetry_type S,
          selection_type Sel,
          env_change_freq_type E,
          adaptation_period A,
          evaluation_type Ev>
void run_simulation_given_resp_type(const cxxopts::ParseResult& results)
{
    auto response_type = convert_net_args(results).resp_type;
    if(response_type == response_type::constitutive)
    {
        run_simulation_given_mut_type< response_type::constitutive, S, Sel, E, A, Ev>(results);
    }
    else if(response_type == response_type::plastic)
    {
        run_simulation_given_mut_type< response_type::plastic, S, Sel, E, A, Ev>(results);
    }
    else{
        throw std::runtime_error{"unknown change type"};
    }
}

template<selection_type S,
         env_change_freq_type E,
         adaptation_period A,
         evaluation_type Ev
         >
void run_simulation_given_env_change_type(const cxxopts::ParseResult& results)
{
    auto env_change_symmetry_type = convert_sim_args(results).change_sym_type ;
    if(env_change_symmetry_type == env_change_symmetry_type::asymmetrical)
    {
        run_simulation_given_resp_type<env_change_symmetry_type::asymmetrical, S, E, A, Ev>(results);
    }
    else if(env_change_symmetry_type == env_change_symmetry_type::symmetrical)
    {
        run_simulation_given_resp_type<env_change_symmetry_type::symmetrical, S, E, A, Ev>(results);
    }
    else{
        throw std::runtime_error{"unknown change type"};
    }
}

template< env_change_freq_type E,
          adaptation_period A,
          evaluation_type Ev
          >
void run_simulation_given_sel_type(const cxxopts::ParseResult& results)
{
    auto sel_type = convert_sim_args(results).sel_type ;
    if(sel_type == selection_type::constant)
    {
        run_simulation_given_env_change_type<selection_type::constant, E, A, Ev>(results);
    }
    else if(sel_type == selection_type::sporadic)
    {
        run_simulation_given_env_change_type< selection_type::sporadic, E, A, Ev>(results);
    }
    else{
        throw std::runtime_error{"unknown change type"};
    }
}

template<adaptation_period A,
         evaluation_type Ev>
void run_simulation_given_env_freq_type(const cxxopts::ParseResult& results)
{
    auto env_change_freq_type = convert_sim_args(results).change_freq_type ;
    if(env_change_freq_type == env_change_freq_type::regular)
    {
        run_simulation_given_sel_type<env_change_freq_type::regular,A, Ev>(results);
    }
    else if(env_change_freq_type == env_change_freq_type::stochastic)
    {
        run_simulation_given_sel_type<env_change_freq_type::stochastic, A, Ev>(results);
    }
    else{
        throw std::runtime_error{"unknown change type"};
    }
}

template<evaluation_type Ev>
void run_simulation_given_adaptation_period(const cxxopts::ParseResult& results)
{
    auto adaptation_period = convert_sim_args(results).adaptation_per ;
    if(adaptation_period == adaptation_period::on)
    {
        run_simulation_given_env_freq_type<adaptation_period::on, Ev>(results);
    }
    else if(adaptation_period == adaptation_period::off)
    {
        run_simulation_given_env_freq_type<adaptation_period::off, Ev>(results);
    }
    else{
        throw std::runtime_error{"unknown adaptation period command"};
    }
}

void run_simulation_given_evaluation_type(const cxxopts::ParseResult& results)
{
    auto evaluation_type = convert_sim_args(results).evalu_type;
    if(evaluation_type == evaluation_type::trial)
    {
        run_simulation_given_adaptation_period<evaluation_type::trial>(results);
    }
    else if(evaluation_type == evaluation_type::full_rn)
    {
        run_simulation_given_adaptation_period<evaluation_type::full_rn>(results);
    }
    else{
        throw std::runtime_error{"unknown adaptation period command"};
    }

}
#endif // LAUNCHER_FUNCTIONS_H
