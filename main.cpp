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
    test_population();
    test_simulation();
    test_weight();
}
#endif

int main(int argc, char ** argv) //!OCLINT tests may be long
{

    auto results = create_parser().parse(argc,argv);

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

    all_params params{
        convert_env_args(results),
                convert_ind_args(results),
                convert_pop_args(results),
                convert_sim_args(results)
    };

    simulation s{params};
    observer o;
    //exec(s, o);

    save_json(o,
              convert_arc_to_string(params.i_p.net_par.net_arc) +
              "_" + std::to_string(params.s_p.seed) + ".json");

    return 0;
}
