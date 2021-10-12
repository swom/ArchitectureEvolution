#include "parser.h"
#include <cassert>
#include <string>
#include <vector>


#ifndef NDEBUG
void test() {
    test_environment();
    test_individual();
    test_network();
    test_observer();
    test_parser();
    test_population();
    test_simulation();
}
#endif

int main(int argc, char ** argv) //!OCLINT tests may be long
{
    const std::vector<std::string> args(argv, argv + argc);
#ifndef NDEBUG
    if (args.size() > 1 && args[1] == "--test")
    {
        test();
        // We've already tested, so the program is done
        return 0;
    }
#else
    // In release mode, all asserts are removed from the code
    assert(1 == 2);
#endif


    auto results = create_parser().parse(argc,argv);
    all_params params{
        convert_env_args(results),
                convert_ind_args(results),
                convert_pop_args(results),
                convert_sim_args(results)
    };



    simulation s{params};
    observer o;
    exec(s, o);

    save_json(o, convert_arc_to_string(params.i_p.net_par.net_arc)+ "_" + std::to_string(params.s_p.seed) + ".json");

    return 0;
}
