#ifndef MUTATION_TYPE_H
#define MUTATION_TYPE_H


enum class mutation_type
{
 weights,
    activation,
    weights_and_activation,
    duplication
};

void test_mutation_type();

#endif // MUTATION_TYPE_H
