#include "launcher_functions.h"
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

     //   run_simulation_given_sel_type(results);

    }  catch (int exc) {
        if(exc==1)
            std::cerr << "The current and maximum architectures are not compatible";
        if(exc==2)
            std::cerr << "A wrong mutation type has been entered";
    }

    return 0;
}
