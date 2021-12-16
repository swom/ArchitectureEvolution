#include "parser.h"
#include <cassert>
#include <string>
#include <vector>


#ifndef NDEBUG
void test() {
        test_environment();
        test_individual();
        test_mutation_type();
        test_network();
        test_observer();
        test_population();
        test_simulation();
        test_weight();
        test_node();
}
#endif

void run_simulation_given_arguments(const cxxopts::ParseResult& results)
{
    auto mut_type = convert_ind_args(results).m_mutation_type;

    if(mut_type == mutation_type::weights)
    {

        using net_t = network<mutation_type::weights>;
        using ind_t = individual<net_t>;
        using pop_t = population<ind_t>;
        using sim_t = simulation<pop_t>;

        observer<sim_t> o;
        auto s = create_simulation<pop_t>(results);
        exec<sim_t>(s, o) ;
        save_json(o,
                  convert_mut_type_to_string(o.get_params().i_p.m_mutation_type) +
                  "_" + convert_arc_to_string(o.get_params().i_p.net_par.net_arc) +
                  "_Changefreq_" + std::to_string(o.get_params().s_p.change_freq).substr(0, 5) +
                  "_" + std::to_string(o.get_params().s_p.seed) + ".json");
    }
    else if (mut_type == mutation_type::activation) {

        using net_t = network<mutation_type::activation>;
        using ind_t = individual<net_t>;
        using pop_t = population<ind_t>;
        using sim_t = simulation<pop_t>;

        observer<sim_t> o;
        auto s = create_simulation<pop_t>(results);
        exec<sim_t>(s, o) ;
        save_json(o,
                  convert_mut_type_to_string(o.get_params().i_p.m_mutation_type) +
                  "_" + convert_arc_to_string(o.get_params().i_p.net_par.net_arc) +
                  "_Changefreq_" + std::to_string(o.get_params().s_p.change_freq).substr(0, 5) +
                  "_" + std::to_string(o.get_params().s_p.seed) + ".json");
    }
    else if (mut_type == mutation_type::weights_and_activation) {

        using net_t = network<mutation_type::weights_and_activation>;
        using ind_t = individual<net_t>;
        using pop_t = population<ind_t>;
        using sim_t = simulation<pop_t>;

        observer<sim_t> o;
        auto s = create_simulation<pop_t>(results);
        exec<sim_t>(s, o) ;
        save_json(o,
                  convert_mut_type_to_string(o.get_params().i_p.m_mutation_type) +
                  "_" + convert_arc_to_string(o.get_params().i_p.net_par.net_arc) +
                  "_Changefreq_" + std::to_string(o.get_params().s_p.change_freq).substr(0, 5) +
                  "_" + std::to_string(o.get_params().s_p.seed) + ".json");
    }
    else if (mut_type == mutation_type::duplication) {

        using net_t = network<mutation_type::duplication>;
        using ind_t = individual<net_t>;
        using pop_t = population<ind_t>;
        using sim_t = simulation<pop_t>;

        observer<sim_t> o;
        auto s = create_simulation<pop_t>(results);
        exec<sim_t>(s, o) ;
        save_json(o,
                  convert_mut_type_to_string(o.get_params().i_p.m_mutation_type) +
                  "_" + convert_arc_to_string(o.get_params().i_p.net_par.net_arc) +
                  "_Changefreq_" + std::to_string(o.get_params().s_p.change_freq).substr(0, 5) +
                  "_Duprate_" + std::to_string(o.get_params().p_p.mut_rate_duplication).substr(0, 5)+
                  "_Maxarch_" + convert_arc_to_string(o.get_params().i_p.net_par.max_arc) +
                  "_" + std::to_string(o.get_params().s_p.seed) + ".json");
    }
    else
    {
        throw std::runtime_error{"unknown mutation type"};
    }
}

int main(int argc, char ** argv) //!OCLINT tests may be long
{

    auto results = create_parser().parse(argc,argv);

    try {
#ifndef NDEBUG
    if (results.count("test"))
    {
        test();
        // We've already tested, so the program is done
        return 0;
    }
#else
    // In release mode, all asserts are removed from the code
    assert(1 == 2);
#endif

    run_simulation_given_arguments(results);

    }  catch (int exc) {
        if(exc==1)
            throw std::runtime_error("The current and maximum architectures are not compatible");
        if(exc==2)
            throw std::runtime_error("A wrong mutation type has been entered");

    }

    return 0;
}
