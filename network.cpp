#include "network.h"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <numeric>


network::network(net_param n_p):
    m_input_size{n_p.net_arc[0]},
    m_activation_function{n_p.function}
{
    for (size_t i = 1; i != n_p.net_arc.size(); i++ )
    {
        std::vector<std::vector<weight>>temp_layer_vector;
        size_t n_nodes_prev_layer = n_p.net_arc[i-1];
        for(int j = 0; j != n_p.net_arc[i]; j++)
        {
            std::vector<weight> temp_weights(n_nodes_prev_layer);
            temp_layer_vector.push_back(temp_weights);
        }

        //A vector of the size of the number of connections is pushed back in the weight matrix
        m_network_weights.push_back(temp_layer_vector);

        //A vector of the size of the nodes in the layer is pushed back;
        m_nodes_biases.push_back(std::vector<double>(n_p.net_arc[i],0));
    }
}


network::network(std::vector<int> nodes_per_layer, std::function<double(double)> activation_function):
    m_input_size{nodes_per_layer[0]},
    m_activation_function{activation_function}
{
    for (size_t i = 1; i != nodes_per_layer.size(); i++ )
    {
        std::vector<std::vector<weight>>temp_layer_vector;
        size_t n_nodes_prev_layer = nodes_per_layer[i-1];
        for(int j = 0; j != nodes_per_layer[i]; j++)
        {
            std::vector<weight> temp_weights(n_nodes_prev_layer);
            temp_layer_vector.push_back(temp_weights);
        }

        //A vector of the size of the number of connections is pushed back in the weight matrix
        m_network_weights.push_back(temp_layer_vector);

        //A vector of the size of the nodes in the layer is pushed back;
        m_nodes_biases.push_back(std::vector<double>(nodes_per_layer[i],0));
    }
}

bool operator==(const network& lhs, const network& rhs)
{
    return lhs.get_input_size() == rhs.get_input_size() &&
            lhs.get_net_weights() == rhs.get_net_weights();
}

bool operator!=(const network& lhs, const network& rhs)
{
    return !(lhs == rhs);
}

network change_all_weights(network n, double new_weight)
{
    for(auto& layer : n.get_net_weights())
        for(auto& node : layer)
            for(auto& weight : node)
            {
                weight.change_weight(new_weight);
            }
    return n;
}

network change_all_weights(network n, weight new_weight)
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

std::vector<weight> register_n_weight_mutations(network n, double mut_rate, double mut_step, std::mt19937_64& rng, int repeats)
{
    std::vector<weight> networks_weights;
    for(int i = 0; i != repeats; i++)
    {

        auto n_new = n;
        n_new.mutate_weights(mut_rate, mut_step, rng);

        for(auto& layer : n_new.get_net_weights())
            for(auto& node : layer)
                for(auto& weight : node)
                {
                    if(weight.get_weight() < -0.0001 || weight.get_weight() > 0.0001)
                        networks_weights.push_back(weight);
                }
    }
    return  networks_weights;
}

std::vector<weight> register_n_activation_mutations(network n, double mut_rate, std::mt19937 &rng, int repeats)
{
    std::vector<weight> networks_weights;
    for(int i = 0; i != repeats; i++)
    {
       auto n_new = n;
        n_new.mutate_activation(mut_rate, rng);
        auto weights = n_new.get_net_weights();

        for(size_t j=0; j != weights.size(); ++j)
            for(size_t k=0; k != weights[j].size(); ++k)
                for(size_t l=0; l != weights[j][k].size(); ++l)
                {
                    if(!weights[j][k][l].is_active())
                        networks_weights.push_back(weights[j][k][l]);
                }
    }
    return  networks_weights;
}


double sigmoid(double x)
{
    return x / (1 + std::abs(x));
}

double linear(double x)
{
    return x;
}

void network::mutate_weights(const double& mut_rate,
                     const double& mut_step,
                     std::mt19937_64& rng)
{

    std::bernoulli_distribution mut_p{mut_rate};
    std::normal_distribution<double> mut_st{0,mut_step};

    for(auto& layer : m_network_weights)
        for(auto& node : layer)
            for(auto& weight : node)
            {
                if(mut_p(rng))
                {weight.change_weight(weight.get_weight() + mut_st(rng));}
            }

    for(auto& layer : m_nodes_biases)
        for(auto& bias : layer)
        {
            if(mut_p(rng))
            {bias += mut_st(rng);}
        }
}

void network::mutate_activation(const double &mut_rate, std::mt19937 &rng)
{
  std::bernoulli_distribution mut_p{mut_rate};

  for(auto& layer : m_network_weights)
      for(auto& node : layer)
          for(auto& weight : node)
          {
              if(mut_p(rng))
              {weight.change_activation(!weight.is_active());}
          }
}



std::vector<double> response(const network& n, std::vector<double> input)
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

            output[node] = n(node_value);
        }

        input = std::move(output);
    }

    return input;
}

bool net_behaves_like_the_function(const network &n, const std::function<double(std::vector<double>)> &f, int n_repeats)
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

bool all_weigths_are_active(const network &n)
{
 auto weights = n.get_net_weights();

  for(auto &layer : weights ){
      for(auto &node : layer){
          for (auto &weight : node){
              if(!weight.is_active()){
                  return false;
                }
            }
        }
    }
  return true;
}

bool on_average_an_nth_of_the_activations_are_mutated(const network &n, const std::vector<weight>&registered_mutations,
                                                      const double &mut_rate, int repeats)
{
  int number_of_weights = get_number_weights(n);
  size_t number_mutations = registered_mutations.size();
  double error = 0.1 * repeats;
  return !(number_mutations < mut_rate * repeats * number_of_weights - error ||
           number_mutations > mut_rate * repeats * number_of_weights + error);
}

int get_number_weights(const network &n)
{
  size_t number_weights = 0;
  for(const auto &layer : n.get_net_weights() ){
      for(const auto &node : layer){
          number_weights += node.size();
            }
        }
    return (int) number_weights;
}



#ifndef NDEBUG
void test_network() //!OCLINT
{
    ///A network can be initialized
    /// by a vector of int that is going to
    /// define the number of layers(size of the vector)
    /// and nodes in those layers(values of the vector)
    {
        std::vector<int> testvec{1, 2, 3, 1} ;
        network n {testvec} ;

        //The number of rows of the weight matrix
        //should be the number elements of the vector
        //provided to the constructor
        // minus 1 (the input layer does not have connections to previous layers)
        assert(n.get_net_weights().size() == testvec.size() -1);


        //The size of the input is stored even though is not saved in n_connection_weights
        //
        for(size_t  i = 0; i != n.get_net_weights().size(); i++)
        {
            assert (static_cast<int>(n.get_net_weights()[i].size()) == testvec[i + 1]);
        }

    }
    ///A network can be initialized with a specific activation function
    {
        network n{{1,2,1}, linear};
        n = change_all_weights(n,1);
        std::vector<double> input{1};
        auto using_member_function = response(n,{input});
        auto using_given_function = response(n,input, &linear);
        auto using_different_given_function = response(n, input, &sigmoid);
        assert(using_given_function == using_member_function);
        assert(using_different_given_function != using_member_function);
    }

    //#ifdef FIX_ISSUE_46
    /// A network can be initilized with a parameter struct net_param
    {
        std::vector<int> net_arc{1, 2, 3, 1} ;
        std::function<double(double)> function = linear;
        net_param n_p{net_arc, function};
        network n{n_p};
        //Set weigths to one
        n = change_all_weights(n, 1);
        //Check architecture
        assert(n.get_net_weights().size() == net_arc.size() -1);
        for(size_t  i = 0; i != n.get_net_weights().size(); i++)
        {
            assert (static_cast<int>(n.get_net_weights()[i].size()) == net_arc[i + 1]);
        }
        //Check activation func
        std::vector<double> input{1};
        auto using_member_function = response(n,{input});
        auto using_given_function = response(n,input, &linear);
        auto using_different_given_function = response(n, input, &sigmoid);
        assert(using_given_function == using_member_function);
        assert(using_different_given_function != using_member_function);
    }
    //#endif

    ///The function resposne returns the output of the network
    {
        auto very_simple_nodes = std::vector<int>{1,2,1};
        auto input = std::vector<double>{1};

        //For simple net with weights == 0
        auto expected_output = std::vector<double>{0};
        network n{very_simple_nodes};
        auto output = response(n, input, &linear);
        assert(output == expected_output);

        //Let's change the weights of the network to something else than 0(e.g 1)
        double new_weight_value = 1.0;
        n = change_all_weights(n,new_weight_value);
        expected_output = std::vector<double>{2};
        output = response(n,input, &linear);
        assert(output == expected_output);

        //Testing a more complex arhitecture
        std::vector<int> not_too_simple_nodes{1,3,3,1};
        network not_too_simple{not_too_simple_nodes};
        not_too_simple = change_all_weights(not_too_simple, new_weight_value);
        expected_output = {1 * 3 * 3};
        output = response(not_too_simple, input, &linear);
        assert(output == expected_output);

    }

    ///Network weights mutate following a normal distribution
    {
        double mut_rate = 0.01;
        double mut_step = 0.1;
        std::mt19937_64 rng;

        auto expected_mean_value  = 0;
        auto expected_stdev = mut_step;

        int repeats = 100000;

        auto very_simple_nodes = std::vector<int>{1,2,1};
        network n{very_simple_nodes};

        std::vector<weight> networks_weights = register_n_weight_mutations(n,
                                                                    mut_rate,
                                                                    mut_step,
                                                                    rng,
                                                                    repeats);

        std::vector<double> weights_as_doubles = convert_to_double(networks_weights);

        double mean = calc_mean(weights_as_doubles);
        double stdev = calc_stdev(weights_as_doubles);

        assert(mean - expected_mean_value < 0.01 && mean - expected_mean_value > -0.01);
        assert(stdev - expected_stdev < 0.01 && stdev - expected_stdev > -0.01);

    }

    ///A network can be initialized with one bias per node
    /// stored by layer and by node in a vector of vectors
    {
        std::vector<int> net_arch{1,2,2,3,1};
        network n{net_arch};
        for(size_t i = 1; i != net_arch.size(); i++)
        {
            auto n_biases = n.get_biases()[i - 1].size();
            auto n_nodes = static_cast<size_t>(net_arch[i]);
            assert(n_biases == n_nodes);
        }

    }

#define FIX_ISSUE_98
#ifdef FIX_ISSUE_98
    {
        std::mt19937 rng;
        network n{net_param{}};
        double mut_rate = 1;

        assert(all_weigths_are_active(n));
        n.mutate_activation(mut_rate, rng);
        assert(!all_weigths_are_active(n));
    }
#endif


#define FIX_ISSUE_100
#ifdef FIX_ISSUE_100
    {
        std::mt19937 rng;

        std::vector<int> net_arch{5,5,5,5,5};
        network n{net_param{net_arch, linear}};

        std::vector<double> mut_rates{0.5, 0.33, 0.25, 0.56789};
        int repeats = 10000;

        for(const auto& mut_rate : mut_rates)
        {
            assert(all_weigths_are_active(n));
            auto active_connections = register_n_activation_mutations(n, mut_rate, rng, repeats);
            assert(on_average_an_nth_of_the_activations_are_mutated(n, active_connections, mut_rate, repeats));
        }
    }
#endif

#define FIX_ISSUE_87
#ifdef FIX_ISSUE_87
    /// A network contains a vector of vectors of vectors of weight objects
    {
        net_param n_p{};
        network n{n_p};

        std::vector<std::vector<std::vector<weight>>> weights = n.get_net_weights();
    }
#endif

#define FIX_ISSUE_96
#ifdef FIX_ISSUE_96
    /// The response function works as intended when some connections are switched off
    {
        auto very_simple_nodes = std::vector<int>{1,2,1};
        auto input = std::vector<double>{1};

        //Comparing a default network with weights = 0 to a network with other weights, all switched off
        network n_default{very_simple_nodes};
        auto output_default = response(n_default, input, &linear);

        weight new_weight_value{1.0, false};
        network n_changed = change_all_weights(n_default, new_weight_value);
        auto output_changed = response(n_changed,input, &linear);

        assert(output_default == output_changed);
    }
#endif

}
#endif
