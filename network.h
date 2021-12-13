#ifndef NETWORK_H
#define NETWORK_H
#include "utilities.h"
#include <iostream>
#include <random>
#include "json.hpp"
#include "node.h"
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
                                   net_arc,
                                   max_arc
                                   )

    net_param(const std::vector<int>& net_arch = {1,2,1},
              std::function<double(double)> func = linear,
              const std::vector<int>& max_arch = {1,8,1}):
        net_arc{net_arch},
        function{func},
        max_arc{max_arch}
    {};

    std::vector<int> net_arc;
    std::function<double(double)> function;
    std::vector<int> max_arc;
};


///Mutates the weights of a network
template<class Net>
void mutate_weights(Net& n, const double& mut_rate,
                    const double& mut_step,
                    std::mt19937_64& rng)
{

    std::bernoulli_distribution mut_p{mut_rate};
    std::normal_distribution<double> mut_st{0,mut_step};

    for(auto& layer : n.get_net_weights())
        for(auto& node : layer)
            for(size_t i = 0; i != node.get_vec_weights().size(); ++i)
            {
                if(mut_p(rng)){
                    const weight &current_weight = node.get_vec_weights()[i];
                    weight mutated_weight(current_weight.get_weight() + mut_st(rng),
                                          current_weight.is_active());
                    node.change_nth_weight(mutated_weight, i);
                }
            }

}

/////Mutates the weights of a network
//template<class Net>
//void mut_dupl_node(Net& n,
//                  const double& mut_rate,
//                    std::mt19937_64& rng)
//{

//    std::bernoulli_distribution mut_p{mut_rate};

//    for(size_t layer = 0; layer != n.get_current_arc() - 1; layer++)
//    {
//        auto& current_layer = n.get_connection_weights()[layer];

//        for(size_t node = 0; node != current_layer.size(); node++)
//        {

//            const auto& current_node = current_layer[node];

//            if(current_node.m_active && mut_p(rng))
//            {
//               //this returns an iterator if you use std::find()
//                auto free_node = find_first_free_nodes(current_layer);

//                if(free_node != current_layer.end())
//                {
//                        *free_node = current_node;

//                  for(auto& next_node : n.get_connection_weights()[layer + 1])
//                  {
//                      next_node.m_weights[free_node] = next_node.m_weights[node];
//                  }
//                }
//            }
//        }
//    }

//}

///Mutates the activation of the weights of the network - they get switched on and off
template<class Net>
void mutate_activation(Net &n, const double &mut_rate, std::mt19937_64 &rng)
{
    std::bernoulli_distribution mut_p{mut_rate};

    for(auto& layer : n.get_net_weights())
        for(auto& node : layer)
            for(size_t i = 0; i != node.get_vec_weights().size(); ++i)
            {
                if(mut_p(rng)){
                    const weight &current_weight = node.get_vec_weights()[i];
                    weight mutated_weight(current_weight.get_weight(),
                                          !current_weight.is_active());
                    node.change_nth_weight(mutated_weight, i);
                }
            }
}

///Takes a vector of nodes, goes through them & mutates their biases, returns the mutated vector
/// Not quite sure why this would have to be here, but it makes it run happily so let's go for it
std::vector<std::vector<double>> mutate_biases(const double& mut_rate,
                                               const double& mut_step,
                                               std::mt19937_64& rng,
                                               const std::vector<std::vector<double>>& biases);

template<mutation_type M = mutation_type::weights>
class network
{
public:
    network(std::vector<int> nodes_per_layer,
            std::function<double(double)> activation_function = &linear);


        network(const net_param &n_p):
        m_input_size{n_p.net_arc[0]},
        m_activation_function{n_p.function},
        m_current_arc{n_p.net_arc},
        m_max_arc{n_p.max_arc}
    {

        for (size_t i = 1; i != n_p.net_arc.size(); i++ )
        {
            std::vector<node>temp_layer_vector;
            size_t n_nodes_prev_layer = n_p.max_arc[i-1];
            for(int j = 0; j != n_p.max_arc[i]; j++)
            {
                std::vector<weight> temp_weights(n_nodes_prev_layer);

                node temp_node(temp_weights);
                if((j+1) <= n_p.net_arc[i]){
                    temp_node.activate();
               }
                temp_layer_vector.push_back(temp_node);
            }

            //A vector of the size of the number of connections is pushed back in the weight matrix
            m_network_weights.push_back(temp_layer_vector);

            //A vector of the size of the nodes in the layer is pushed back;
            m_nodes_biases.push_back(std::vector<double>(n_p.net_arc[i],0));
        }
        if(!net_arc_and_max_arc_are_compatible(m_current_arc, m_max_arc)){
            throw 1;
          }
    }

    void mutate(const double& mut_rate_weight,
                           const double& mut_step,
                           std::mt19937_64& rng,
                const double& mut_rate_act)
       {

           if constexpr (M == mutation_type::activation)
           {
               mutate_activation(*this, mut_rate_act, rng);
           }

           else if constexpr(M == mutation_type::weights)
           {
               mutate_weights(*this, mut_rate_weight, mut_step, rng);
           }

           else if constexpr(M == mutation_type::weights_and_activation)
           {
               mutate_activation(*this, mut_rate_act, rng);
               mutate_weights(*this, mut_rate_weight, mut_step, rng);
           }
           this->change_biases(mutate_biases(mut_rate_weight, mut_step, rng, this->get_biases()));
       };

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
    const std::vector<std::vector<node>>& get_net_weights() const noexcept{return m_network_weights;}

    ///Returns not constant ref to vector of weights
    std::vector<std::vector<node>>& get_net_weights() noexcept{return m_network_weights;}

    ///Returns the current architecture
    const std::vector<int>& get_current_arc() const noexcept{return m_current_arc;}

    ///Returns the maximum architecture
    const std::vector<int>& get_max_arc() const noexcept{return m_max_arc;}

    void change_network_arc(std::vector<int> new_arc);

    std::vector<node>::iterator get_empty_node_in_layer(size_t l);

    void duplicate_node(const node &to_duplicate, size_t layer, size_t index_to_duplicate,
                        const std::vector<node>::iterator &empty_node_iterator);

private:

    ///Vector of of vectors, representing the weights coming into each node
    std::vector<std::vector<node>> m_network_weights;

    ///Vector of vectors containing the nodes biases stored per layer per node
    std::vector<std::vector<double>> m_nodes_biases;

    ///The size of the input vector the network will receive
    int m_input_size;

    ///The activation function of the nodes
    std::function<double(double)> m_activation_function;

    ///The current architecture
    std::vector<int> m_current_arc;

    ///The maximum architecture
    std::vector<int> m_max_arc;

    inline bool net_arc_and_max_arc_are_compatible(const std::vector<int> &net_arc, const std::vector<int> &max_arc)
    {
      if(net_arc.size() != max_arc.size()){
       return false;
        }

      bool there_is_a_wrong_number_in_net_arc = false;
      for(auto & layer : net_arc){
          if(layer <= 0){
            there_is_a_wrong_number_in_net_arc = true;
            }
        }

      bool there_is_a_wrong_number_in_max_arc = false;
      for(auto & layer : max_arc){
          if(layer <= 0){
            there_is_a_wrong_number_in_max_arc = true;
            }
        }

      if(there_is_a_wrong_number_in_max_arc || there_is_a_wrong_number_in_net_arc){
          return false;
        }

      bool net_arc_smaller_than_max_arc = true;
      for(size_t i = 0; i != net_arc.size(); ++i){
          if(net_arc[i] > max_arc[i]){
              net_arc_smaller_than_max_arc = false;
            }
        }

      if(!net_arc_smaller_than_max_arc){
          return false;
        }

      if(net_arc[0] != max_arc[0] || net_arc.back() != max_arc.back()){
          return false;
        }
      else
        return true;
    }
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
            for(size_t i = 0; i != node.get_vec_weights().size(); ++i)
            {
                const weight &current_weight = node.get_vec_weights()[i];
                weight changed_weight(new_weight, current_weight.is_active());
                node.change_nth_weight(changed_weight, i);
            }

    return n;
}

template<class Net>
Net change_all_weights_values_and_activations(Net n, weight new_weight)
{
    for(auto& layer : n.get_net_weights())
        for(auto& node : layer)
            for(size_t i = 0; i != node.get_vec_weights().size(); ++i)
            {
                node.change_nth_weight(new_weight, i);
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
        assert(output(n, input).size() == 1);
        if(output(n, input) [0] != f(input))
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

///Returns the output of a network for a given input
template<class Net>
std::vector<double> output(const Net& n, std::vector<double> input)
{
    assert(input.size() == n.get_input_size());

    for(size_t layer = 0; layer != n.get_net_weights().size(); layer++)
    {
        auto output = std::vector<double>(n.get_net_weights()[layer].size());

        for(size_t node = 0; node != n.get_net_weights()[layer].size(); node++)
        {
            const class node &current_node = n.get_net_weights()[layer][node];
            std::vector<weight> vec_w = current_node.get_vec_weights();
            std::vector<double> w = convert_to_double_or_zero(vec_w);
            double node_value = n.get_biases()[layer][node] +
                    std::inner_product(input.begin(),
                                       input.end(),
                                       w.begin(),
                                       0.0);

            output[node] = n(node_value);
        }

        input = std::move(output);
    }

    return input;
}

///Returns the output of a network for a given input and a given activation function
template <typename Fun, mutation_type M>
inline std::vector<double> output(const network<M>& n, std::vector<double> input, Fun fun)
{
    assert(input.size() == n.get_input_size());

    for(size_t layer = 0; layer != n.get_net_weights().size(); layer++)
    {
        auto output = std::vector<double>(n.get_net_weights()[layer].size());

        for(size_t node = 0; node != n.get_net_weights()[layer].size(); node++)
        {
            const class node &current_node = n.get_net_weights()[layer][node];
            std::vector<weight> vec_w = current_node.get_vec_weights();
            std::vector<double> w = convert_to_double_or_zero(vec_w);
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
    if(typeid(lhs) != typeid(rhs))
    {
        return false;
    }

    if(!are_equal_except_mutation_type(lhs, rhs))
    {
        return false;
    }

    return true;
}



void test_network();

#endif // NETWORK_H
