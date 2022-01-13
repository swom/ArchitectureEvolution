#ifndef MUTATION_TYPE_H
#define MUTATION_TYPE_H
#include <iostream>


enum class mutation_type
{
 weights,
    activation,
    weights_and_activation,
    duplication,
    addition,
  NRduplication,
  NRaddition
};

std::string convert_mut_type_to_string(mutation_type m);

void test_mutation_type();

#endif // MUTATION_TYPE_H
