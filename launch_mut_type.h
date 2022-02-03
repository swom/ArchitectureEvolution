#ifndef LAUNCH_MUT_TYPE_H
#define LAUNCH_MUT_TYPE_H
#include "parser.h"

template <env_change_type env_change_type, selection_type selection_type>
void run_simulation_given_mut_type(const cxxopts::ParseResult& results)
{
    auto mut_type = convert_ind_args(results).m_mutation_type;

    if(mut_type == mutation_type::weights)
    {
        using net_t = network<mutation_type::weights>;
        using ind_t = individual<net_t>;
        using pop_t = population<ind_t>;
        using sim_t = simulation<pop_t, env_change_type, selection_type>;
        observer<sim_t> o;

        auto s = create_simulation<pop_t, env_change_type, selection_type>(results);
        exec<sim_t>(s, o) ;
        save_json(o,
                  create_save_name_from_observer_data(o));
    }
    else if (mut_type == mutation_type::activation) {

        using net_t = network<mutation_type::activation>;
        using ind_t = individual<net_t>;
        using pop_t = population<ind_t>;
        using sim_t = simulation<pop_t, env_change_type, selection_type>;

        observer<sim_t> o;
        auto s = create_simulation<pop_t, env_change_type, selection_type>(results);
        exec<sim_t>(s, o) ;
        save_json(o,
                  create_save_name_from_observer_data(o));
    }
    else if (mut_type == mutation_type::weights_and_activation) {

        using net_t = network<mutation_type::weights_and_activation>;
        using ind_t = individual<net_t>;
        using pop_t = population<ind_t>;
        using sim_t = simulation<pop_t, env_change_type, selection_type>;

        observer<sim_t> o;
        auto s = create_simulation<pop_t, env_change_type, selection_type>(results);
        exec<sim_t>(s, o) ;
        save_json(o,
                  create_save_name_from_observer_data(o));
    }
    else if (mut_type == mutation_type::duplication) {

        using net_t = network<mutation_type::duplication>;
        using ind_t = individual<net_t>;
        using pop_t = population<ind_t>;
        using sim_t = simulation<pop_t, env_change_type, selection_type>;

        observer<sim_t> o;
        auto s = create_simulation<pop_t, env_change_type, selection_type>(results);
        exec<sim_t>(s, o) ;
        save_json(o,
                  create_save_name_from_observer_data(o));
    }
    else if (mut_type == mutation_type::NRduplication) {

        using net_t = network<mutation_type::NRduplication>;
        using ind_t = individual<net_t>;
        using pop_t = population<ind_t>;
        using sim_t = simulation<pop_t, env_change_type, selection_type>;

        observer<sim_t> o;
        auto s = create_simulation<pop_t, env_change_type, selection_type>(results);
        exec<sim_t>(s, o) ;
        save_json(o,
                  create_save_name_from_observer_data(o));
    }
    else if (mut_type == mutation_type::addition) {

        using net_t = network<mutation_type::addition>;
        using ind_t = individual<net_t>;
        using pop_t = population<ind_t>;
        using sim_t = simulation<pop_t, env_change_type, selection_type>;

        observer<sim_t> o;
        auto s = create_simulation<pop_t, env_change_type, selection_type>(results);
        exec<sim_t>(s, o) ;
        save_json(o,
                  create_save_name_from_observer_data(o));
    }
    else if (mut_type == mutation_type::NRaddition) {

        using net_t = network<mutation_type::NRaddition>;
        using ind_t = individual<net_t>;
        using pop_t = population<ind_t>;
        using sim_t = simulation<pop_t, env_change_type, selection_type>;

        observer<sim_t> o;
        auto s = create_simulation<pop_t, env_change_type, selection_type>(results);
        exec<sim_t>(s, o) ;
        save_json(o,
                  create_save_name_from_observer_data(o));
    }
    else
    {
        throw std::runtime_error{"unknown mutation type"};
    }
}

#endif // LAUNCH_MUT_TYPE_H
