#ifndef NETWORK_H
#define NETWORK_H
#include "utilities.h"
#include <iostream>
#include <random>
#include "json.hpp"
#include "weight.h"
#include "mutation_type.h"

double sigmoid(double x);
double linear(double x);

static std::map<std::string, std::function<double(double)>> string_to_act_func_map
{
{"linear", linear},
{"sigmoid", sigmoid}
};

struct net_param
{
    NLOHMANN_DEFINE_TYPE_INTRUSIVE(net_param,
                                   net_arc
                                   )
    net_param(const std::vector<int>& net_arch = {1,2,1},
              std::function<double(double)> func = linear):
        net_arc{net_arch},
        function{func}
    {};

    std::vector<int> net_arc;
    std::function<double(double)> function;
};


class network
{
public:
    network(std::vector<int> nodes_per_layer, std::function<double(double)> activation_function = &linear);
    network (net_param n_p);

    virtual ~network() {}
    void mutate(const double& ,
                        const double& ,
                        std::mt19937_64& )
    {};

    NLOHMANN_DEFINE_TYPE_INTRUSIVE(network,
                                   m_input_size, 
                                   m_network_weights
                                   );

    ///Returns the activation function
    std::function<double(double)> get_activation_function() const noexcept{return m_activation_function;}

    ///Returns the const ref to the node biases
    const std::vector<std::vector<double>>& get_biases() const noexcept{return m_nodes_biases;}


    ///Returns the input size
    size_t get_input_size() const noexcept {return static_cast<size_t>(m_input_size);}

    ///Mutates the weights of the network
    void mutate_weights(const double& mut_rate, const double& mut_step, std::mt19937_64 &rng);

    ///Mutates the activation of the weights of the network - they get switched on and off
    void mutate_activation(const double &mut_rate, std::mt19937_64 &rng);

    double operator ()(double n) const {return m_activation_function(n);}

    ///Returns const ref to vector of weights
    const std::vector<std::vector<std::vector<weight>>>& get_net_weights() const noexcept{return m_network_weights;}

    ///Returns not constant ref to vector of weights
    std::vector<std::vector<std::vector<weight>>>& get_net_weights() noexcept{return m_network_weights;}
private:
    ///Vector of of vectors, representing the weights coming into each node
    std::vector<std::vector<std::vector<weight>>> m_network_weights;

    ///Vector of vectors containing the nodes biases stored per layer per node
    std::vector<std::vector<double>> m_nodes_biases;

    ///The size of the input vector the network will receive
    int m_input_size;

    ///The activation function of the nodes
    std::function<double(double)> m_activation_function;
};


template <mutation_type mutation_type>
class mutator_network : public network
{
public:
    mutator_network(const net_param& p) : network{p} {};
    void mutate(const double& mut_rate,
                const double& mut_step,
                std::mt19937_64& rng);
};

bool operator==(const network& lhs, const network& rhs);

bool operator!=(const network& lhs, const network& rhs);

network change_all_weights(network n, double new_weight);
network change_all_weights(network n, weight new_weight);

///Mutates a network n times with given mutation
/// rate and step and returns vector all mutated weights
/// of network in all times it was mutated
std::vector<weight> register_n_weight_mutations(network n, double mut_rate, double mut_step, std::mt19937_64 &rng, int repeats);

std::vector<weight> register_n_activation_mutations(network n, double mut_rate, std::mt19937_64 &rng, int repeats);

template <typename Fun>
inline std::vector<double> response(const network& n, std::vector<double> input, Fun fun = &linear)
{
    assert(input.size() == n.get_input_size());

    for(size_t layer = 0; layer != n.get_net_weights().size(); layer++)
    {
        auto output = std::vector<double>(n.get_net_weights()[layer].size());

        for(size_t node = 0; node != n.get_net_weights()[layer].size(); node++)
        {
            std::vector<double> w = convert_to_double_or_zero(n.get_net_weights()[layer][node]);
            double node_value = n.get_biases()[layer][node] +
                    std::inner_product(input.begin(),
                                       input.end(),
                                       w.begin(),
                                       0.0);

            output[node] = fun(node_value);
        }
        input = std::move(output);
    }

    return input;
}


std::vector<double> response(const network& n, std::vector<double> input);

///Checks if a network and a function return the same output
bool net_behaves_like_the_function(const network &n, const std::function<double(std::vector<double>)> &f, int n_repeats = 1000);

///Checks whether all connections of the network are active
bool all_weigths_are_active(const network &n);

///Checks that all weights have a certain value
bool all_weigths_have_value(const network &n, double value);

///Checks that the registered_mutations correspond to the given mutation rate
bool on_average_an_nth_of_the_weights_are_inactive(const network &n, const std::vector<weight>&registered_mutations,
                                                      const double &proportion, int repeats);

///Returns the total number of connections in the network
int get_number_weights(const network &n);



void test_network();

#endif // NETWORK_H
