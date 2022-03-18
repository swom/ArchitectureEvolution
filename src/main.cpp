#include "launcher_functions.h"
#include <cassert>
#include <string>
#include <vector>


#ifndef NDEBUG
void test() {
    std::cout << "testing environment" << std::endl;
    test_environment();
    std::cout << "testing individual" << std::endl;
    test_individual();
    std::cout << "testing mutation_type" << std::endl;
    test_mutation_type();
    std::cout << "testing network" << std::endl;
    test_network();
    std::cout << "testing observer" << std::endl;
    test_observer();
    std::cout << "testing population" << std::endl;
    test_population();
    std::cout << "testing simulation" << std::endl;
    test_simulation();
    std::cout << "testing weight" << std::endl;
    test_weight();
    std::cout << "testing node" << std::endl;
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

       run_simulation_given_adaptation_period(results);

    }  catch (int exc) {
        if(exc==2)
            std::cerr << "A wrong mutation type has been entered";
    }
    catch( const std::runtime_error& e ) {
        std::cerr << e.what() << std::endl;
    }
    catch( const std::invalid_argument& e ) {
        std::cerr << e.what() << std::endl;
    };

    return 0;
}
