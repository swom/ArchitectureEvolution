#include "mutation_type.h"
#include "cassert"

std::string convert_mut_type_to_string(mutation_type m)
{
    std::string string;

    switch (m) {
    case mutation_type::weights :
        string = "weights";
        return string;
        break;

    case mutation_type::weights_and_activation :
        string = "weights_and_activation";
        return string;
        break;

    case mutation_type::activation :
        string = "activation";
        return string;
        break;

    case mutation_type::duplication :
        string = "duplication";
        return string;
        break;

      case mutation_type::addition :
          string = "addition";
          return string;
          break;

      case mutation_type::NRduplication :
          string = "NRduplication";
          return string;
          break;

      case mutation_type::NRaddition :
          string = "NRaddition";
          return string;
          break;

    default:
        throw 2;
    }
}

void test_mutation_type()
{
    assert(mutation_type::activation != mutation_type::weights &&
            mutation_type::weights != mutation_type::weights_and_activation &&
            mutation_type::weights_and_activation != mutation_type::activation &&
            mutation_type::duplication != mutation_type::activation &&
            mutation_type::duplication != mutation_type::weights &&
            mutation_type::duplication != mutation_type::weights_and_activation &&
            mutation_type::addition != mutation_type::weights &&
            mutation_type::addition != mutation_type::weights_and_activation &&
            mutation_type::addition != mutation_type::activation &&
            mutation_type::addition != mutation_type::duplication);
}

