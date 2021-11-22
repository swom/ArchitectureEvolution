#include "mutation_type.h"
#include "cassert"




#ifndef NDEBUG
void test_mutation_type()
{

#define FIX_ISSUE_108
#ifdef FIX_ISSUE_108
    {
        assert(mutation_type::activation != mutation_type::weights &&
                mutation_type::weights != mutation_type::weights_and_activation &&
                mutation_type::weights_and_activation != mutation_type::activation);
    }
#endif

}
#endif
