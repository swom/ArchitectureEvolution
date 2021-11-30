#include "mutation_type.h"
#include "cassert"
void test_mutation_type()
{
    assert(mutation_type::activation != mutation_type::weights &&
            mutation_type::weights != mutation_type::weights_and_activation &&
            mutation_type::weights_and_activation != mutation_type::activation);
}

