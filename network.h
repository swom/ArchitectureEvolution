#ifndef NETWORK_H
#define NETWORK_H
#include "utilities.h"
#include <iostream>
#include <vector>
#include <random>
#include "json.hpp"


double sigmoid(double x);
double linear(double x);

static std::map<std::string, std::function<double(double)>> string_to_act_func_map
{
{"linear", linear},
{"sigmoid", sigmoid}
};

//static std::map<std::function<double(double)>, std::string> act_funct_to_string_map
//{
//{linear, "linear"},
//{sigmoid, "sigmoid"}
//};

struct net_param
{
    NLOHMANN_DEFINE_TYPE_INTRUSIVE(net_param,
                                   net_arc
                                   )

    std::vector<int> net_arc {1,2,1};
    std::function<double(double)> function;
//    std::string str_func = act_funct_to_string_map.find(function)->second;
};


class network
{
public:
    network(std::vector<int> nodes_per_layer, std::function<double(double)> activation_function = &linear);
    network (net_param n_p);

    NLOHMANN_DEFINE_TYPE_INTRUSIVE(network,
                                   m_input_size,
                                   m_network_weights);

    ///Returns the activation function
    std::function<double(double)> get_activation_function() const noexcept{return m_activation_function;}

    ///Returns the const ref to the node biases
    const std::vector<std::vector<double>>& get_biases() const noexcept{return m_nodes_biases;}

    ///Returns const ref to vector of weights
    const std::vector<std::vector<std::vector<double>>>& get_net_weights() const noexcept{return m_network_weights;}

    ///Returns not constant ref to vector of weights
    std::vector<std::vector<std::vector<double>>>& get_net_weights() noexcept{return m_network_weights;}

    ///Returns the input size
    size_t get_input_size() const noexcept {return static_cast<size_t>(m_input_size);}

    ///Mutates the network
    void mutate(const double& mut_rate, const double& mut_step, std::mt19937_64 &rng);

    double operator ()(double n) const {return m_activation_function(n);}

private:
    ///Vector of of vectors, representing the weights coming into each node
    std::vector<std::vector<std::vector<double>>> m_network_weights;

    ///Vector of vectors containing the nodes biases stored per layer per node
    std::vector<std::vector<double>> m_nodes_biases;

    ///The size of the input vector the network will receive
    int m_input_size;

    ///The activation function of the nodes
    std::function<double(double)> m_activation_function;
};

bool operator==(const network& lhs, const network& rhs);

bool operator!=(const network& lhs, const network& rhs);

network change_all_weights(network n, double new_weight);

///Mutates a network n times with given mutation
/// rate and step and returns vector all mutated weights
/// of network in all times it was mutated
std::vector<double> register_n_mutations(network n, double mut_rate, double mut_step, std::mt19937_64 &rng, int repeats);

template <typename Fun>
inline std::vector<double> response(const network& n, std::vector<double> input, Fun fun = &linear)
{
    assert(input.size() == n.get_input_size());

    for(size_t layer = 0; layer != n.get_net_weights().size(); layer++)
    {
        auto output = std::vector<double>(n.get_net_weights()[layer].size());

        for(size_t node = 0; node != n.get_net_weights()[layer].size(); node++)
        {
            double node_value = n.get_biases()[layer][node] +
                    std::inner_product(input.begin(),
                                       input.end(),
                                       n.get_net_weights()[layer][node].begin(),
                                       0.0);

            output[node] = fun(node_value);
        }
        input = std::move(output);
    }

    return input;
}

std::vector<double> response(const network& n, std::vector<double> input);


void test_network();

#endif // NETWORK_H
