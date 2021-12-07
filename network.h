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

template<mutation_type M = mutation_type::weights>
class network
{
public:
    network(std::vector<int> nodes_per_layer,
            std::function<double(double)> activation_function = &linear);
    network (const net_param &n_p);

    virtual ~network() {}
    virtual void mutate(const double& ,
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

    ///Sets the value of the nodes to the given value
    void change_biases(std::vector<std::vector<double>> new_biases) {m_nodes_biases = new_biases;}

    ///Returns the input size
    size_t get_input_size() const noexcept {return static_cast<size_t>(m_input_size);}

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


///Takes a vector of nodes, goes through them & mutates their biases, returns the mutated vector
/// Not quite sure why this would have to be here, but it makes it run happily so let's go for it
std::vector<std::vector<double>> mutate_biases(const double& mut_rate,
                                               const double& mut_step,
                                               std::mt19937_64& rng,
                                               const std::vector<std::vector<double>>& biases);

template <mutation_type M>
class mutator_network : public network<M>
{
public:
    mutator_network(const net_param& p) : network<M>{p} {};

    virtual void mutate(const double& mut_rate,
                        const double& mut_step,
                        std::mt19937_64& rng) override
    {

        if constexpr (M == mutation_type::activation)
        {
            mutate_activation(*this, mut_rate, rng);
        }

        else if constexpr(M == mutation_type::weights)
        {
            mutate_weights(*this, mut_rate, mut_step, rng);
        }

        else if constexpr(M == mutation_type::weights_and_activation)
        {
            mutate_activation(*this, mut_rate, rng);
            mutate_weights(*this, mut_rate, mut_step, rng);
        }
        this->change_biases(mutate_biases(mut_rate, mut_step, rng, this->get_biases()));
    };
};

template<mutation_type M>
bool operator==(const network<M>& lhs, const network<M> &rhs)
{
    return lhs.get_input_size() == rhs.get_input_size() &&
            lhs.get_net_weights() == rhs.get_net_weights();
}

template<mutation_type M>
bool operator!=(const network<M>& lhs, const network<M>& rhs)
{
    return !(lhs == rhs);
}

template<mutation_type M_lhs, mutation_type M_rhs>
bool are_equal_except_mutation_type(const network<M_lhs>& lhs, const network<M_rhs> &rhs)
{
    return lhs.get_input_size() == rhs.get_input_size() &&
            lhs.get_net_weights() == rhs.get_net_weights();
}

template<class Net>
Net change_all_weights_values(Net n, double new_weight)
{
    for(auto& layer : n.get_net_weights())
        for(auto& node : layer)
            for(auto& weight : node)
            {
                weight.change_weight(new_weight);
            }
    return n;
}

template<class Net>
Net change_all_weights_values_and_activations(Net n, weight new_weight)
{
    for(auto& layer : n.get_net_weights())
        for(auto& node : layer)
            for(auto& weight : node)
            {
                weight.change_weight(new_weight.get_weight());
                weight.change_activation(new_weight.is_active());
            }
    return n;
}


///Mutates a network n times with given mutation
/// rate and step and returns vector all mutated weights
/// of network in all times it was mutated
template<class Net>
std::vector<weight> register_n_weight_mutations(Net n,
                                                double mut_rate,
                                                double mut_step,
                                                std::mt19937_64 &rng,
                                                int repeats);
template<class Net>
std::vector<weight> register_n_activation_mutations(Net n,
                                                    double mut_rate,
                                                    std::mt19937_64 &rng,
                                                    int repeats);


///Checks if a network and a function return the same output
template<class Net>
bool net_behaves_like_the_function(const Net &n,
                                   const std::function<double(std::vector<double>)> &f,
                                   int n_repeats = 1000)
{
    std::vector<std::vector<double>> input_series;
    size_t input_size = n.get_input_size();
    std::vector<double> n_output;
    std::vector<double> f_output;

    for(int i = 0; i != n_repeats; ++i)
    {
        std::vector<double> input;
        for(size_t j = 0; j != input_size; ++j){
            input.push_back(i + j);
        }
        assert(response(n, input).size() == 1);
        if(response(n, input) [0] != f(input))
        {
            return false;
        };
    }

    return true;
}

///Checks whether all connections of the network are active
template<class Net>
bool all_weigths_are_active(const Net &n);

///Checks that all weights have a certain value
template<class Net>
bool all_weigths_have_value(const Net &n, double value);

///Checks that the registered_mutations correspond to the given mutation rate
template<class Net>
bool on_average_an_nth_of_the_weights_are_inactive(const Net &n,
                                                   const std::vector<weight>& registered_mutations,
                                                   const double &proportion,
                                                   int repeats);

///Returns the total number of connections in the network
template<class Net>
int get_number_weights(const Net &n);

template <typename Fun, class Net>
inline std::vector<double> response(const Net& n, std::vector<double> input, Fun fun = &linear)
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


template<class Net>
std::vector<double> response(const Net& n, std::vector<double> input);

///Checks whether all connections of the network are active
template<class Net>
bool all_weigths_are_active(const Net &n);

///Checks that all weights have a certain value
template<class Net>
bool all_weigths_have_value(const Net &n, double value);

///Checks that the registered_mutations correspond to the given mutation rate
template<class Net>
bool on_average_an_nth_of_the_weights_are_inactive(const Net &n, const std::vector<weight>&registered_mutations,
                                                   const double &proportion, int repeats);

///Returns the total number of connections in the network
template<class Net>
int get_number_weights(const Net &n);

///Checks that both networks are mutator_networks that have the same mutation function
///In addition to checking the normal equality
template<class Net_lhs, class Net_rhs>
bool is_same_mutator_network(const Net_lhs &lhs, const Net_rhs &rhs)
{
    if(!are_equal_except_mutation_type(lhs, rhs)){
        return false;
    }

    if(typeid(lhs) == typeid(rhs)){
        return true;
    }
    else{
        return false;
    }
}

///Mutates the weights of a network
template<class Net>
void mutate_weights(Net &n, const double& mut_rate, const double& mut_step, std::mt19937_64 &rng);

///Mutates the activation of the weights of the network - they get switched on and off
template<class Net>
void mutate_activation(Net &n, const double &mut_rate, std::mt19937_64 &rng);


void test_network();

#endif // NETWORK_H
