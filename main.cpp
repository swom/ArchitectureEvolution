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

    simulation s = create_simulation(results);
    observer o;
    exec(s, o);

    save_json(o,
              convert_arc_to_string(s.get_params().i_p.net_par.net_arc) +
              "_" + std::to_string(s.get_params().s_p.seed) + ".json");

    return 0;
}
