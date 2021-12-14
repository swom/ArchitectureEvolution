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

        observer<mutation_type::weights> o;
        auto s = create_simulation<mutation_type::weights>(results);
        exec<mutation_type::weights>(s, o) ;
        save_json(o,
                  convert_arc_to_string(o.get_params().i_p.net_par.net_arc) +
                  "_" + std::to_string(o.get_params().s_p.seed) + ".json");
    }
    else if (mut_type == mutation_type::activation) {

        observer<mutation_type::activation> o;
        auto s = create_simulation<mutation_type::activation>(results);
        exec<mutation_type::activation>(s, o) ;
        save_json(o,
                  convert_arc_to_string(o.get_params().i_p.net_par.net_arc) +
                  "_" + std::to_string(o.get_params().s_p.seed) + ".json");
    }
    else if (mut_type == mutation_type::weights_and_activation) {
        observer<mutation_type::weights_and_activation> o;
        auto s = create_simulation<mutation_type::weights_and_activation>(results);
        exec<mutation_type::weights_and_activation>(s, o) ;
        save_json(o,
                  convert_arc_to_string(o.get_params().i_p.net_par.net_arc) +
                  "_" + std::to_string(o.get_params().s_p.seed) + ".json");
    }
        else if (mut_type == mutation_type::duplication) {
            observer<mutation_type::duplication> o;
            auto s = create_simulation<mutation_type::duplication>(results);
            exec<mutation_type::duplication>(s, o) ;
            save_json(o,
                      convert_arc_to_string(o.get_params().i_p.net_par.net_arc) +
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

    }

    return 0;
}
