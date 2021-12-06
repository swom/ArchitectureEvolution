#include "network.h"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <numeric>


network::network(const net_param &n_p):
    m_input_size{n_p.net_arc[0]},
    m_activation_function{n_p.function},
    m_current_arc{n_p.net_arc},
    m_max_arc{n_p.max_arc}
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
    if(!net_arc_and_max_arc_are_compatible(m_current_arc, m_max_arc)){
        throw 1;
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
        mutate_weights(n_new, mut_rate, mut_step, rng);

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

std::vector<weight> register_n_activation_mutations(network n, double mut_rate, std::mt19937_64 &rng, int repeats)
{
    std::vector<weight> networks_weights;
    for(int i = 0; i != repeats; i++)
    {
       auto n_new = n;
        mutate_activation(n_new, mut_rate, rng);
        auto weights = n_new.get_net_weights();

        for(size_t j=0; j != weights.size(); ++j)
            for(size_t k=0; k != weights[j].size(); ++k)
                for(size_t l=0; l != weights[j][k].size(); ++l)
                {
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

bool all_weigths_have_value(const network &n, double value)
{
 auto weights = n.get_net_weights();

  for(auto &layer : weights ){
      for(auto &node : layer){
          for (auto &weight : node){
              if(weight.get_weight() != value)
              {
                  return false;
              }
            }
        }
    }
  return true;
}

bool on_average_an_nth_of_the_weights_are_inactive(const network &n, const std::vector<weight>&registered_mutations,
                                                      const double &proportion, int repeats)
{
  int number_of_weights = get_number_weights(n);
  int inactive_weights = 0;

  for (auto &weight : registered_mutations){
      if(!weight.is_active())
        ++inactive_weights;
    }

  double error = 0.1 * repeats;
  return (inactive_weights > proportion * repeats * number_of_weights - error &&
           inactive_weights < proportion * repeats * number_of_weights + error);
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

bool is_same_mutator_network(const network &lhs, const network &rhs)
{
  if(lhs != rhs){
    return false;
    }

  if(typeid(lhs) == typeid(rhs)){
    return true;
    }
  else{
    return false;
    }
}

void mutate_weights(network& n, const double& mut_rate,
                     const double& mut_step,
                     std::mt19937_64& rng)
{

    std::bernoulli_distribution mut_p{mut_rate};
    std::normal_distribution<double> mut_st{0,mut_step};

    for(auto& layer : n.get_net_weights())
        for(auto& node : layer)
            for(auto& weight : node)
            {
                if(mut_p(rng))
                {weight.change_weight(weight.get_weight() + mut_st(rng));}
            }

}

void mutate_activation(network &n, const double &mut_rate, std::mt19937_64 &rng)
{
  std::bernoulli_distribution mut_p{mut_rate};

  for(auto& layer : n.get_net_weights())
      for(auto& node : layer)
          for(auto& weight : node)
          {
              if(mut_p(rng))
              {weight.change_activation(!weight.is_active());}
          }
}

std::vector<std::vector<double>> mutate_biases(const double& mut_rate,
                                               const double& mut_step,
                                               std::mt19937_64& rng,
                                               const std::vector<std::vector<double>>& biases)
{
  std::vector<std::vector<double>> new_biases;
  std::bernoulli_distribution mut_p{mut_rate};
  std::normal_distribution<double> mut_st{0,mut_step};

  for(auto& layer : biases){
    std::vector<double> new_layer;
      for(auto& bias : layer)
      {

          if(mut_p(rng)){
              new_layer.push_back(bias + mut_st(rng));
            }
          else{
              new_layer.push_back(bias);
            }
      }
     new_biases.push_back(new_layer);
    }
  return(new_biases);
}


bool net_arc_and_max_arc_are_compatible(const std::vector<int> &net_arc, const std::vector<int> &max_arc)
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
        std::mt19937_64 rng;
        network n{net_param{}};
        double mut_rate = 1;

        assert(all_weigths_are_active(n));
        mutate_activation(n, mut_rate, rng);
        assert(!all_weigths_are_active(n));
    }
#endif


#define FIX_ISSUE_100
#ifdef FIX_ISSUE_100
    {
        std::mt19937_64 rng;

        std::vector<int> net_arch{5,5,5,5,5};
        network n{net_param{net_arch, linear}};

        std::vector<double> mut_rates{0.5, 0.33, 0.25, 0.56789};
        int repeats = 10000;

        for(const auto& mut_rate : mut_rates)
        {
            assert(all_weigths_are_active(n));
            auto active_connections = register_n_activation_mutations(n, mut_rate, rng, repeats);
            assert(on_average_an_nth_of_the_weights_are_inactive(n, active_connections, mut_rate, repeats));
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

#define FIX_ISSUE_112
#ifdef FIX_ISSUE_112
    {
        mutator_network<mutation_type::activation> n_activation{net_param()};
        assert(all_weigths_are_active(n_activation));

        mutator_network<mutation_type::weights> n_weights{net_param()};
        assert(all_weigths_have_value(n_weights, 0));

        auto mutation_rate = 1;
        auto mutation_step = 1;
        std::mt19937_64 rng;
        auto rng_copy = rng;

        auto before_mutation = n_weights;
        n_weights.mutate(mutation_rate, mutation_step, rng_copy);
        n_activation.mutate(mutation_rate, mutation_step, rng);

        assert(n_activation.get_net_weights() != n_weights.get_net_weights());
        assert(!all_weigths_are_active(n_activation));

        assert(all_weigths_are_active(n_weights));
        assert(before_mutation != n_weights);
    }
#endif

#define FIX_ISSUE_126
#ifdef FIX_ISSUE_126
    {
        auto pars = net_param();
        mutator_network<mutation_type::activation> n_activation{pars};
        mutator_network<mutation_type::weights> n_weights{pars};
        network n{pars};

        assert(n == n_activation);
        assert(n == n_weights);
    }
#endif

#define FIX_ISSUE_159
#ifdef FIX_ISSUE_159
  ///Network has a (current) network architecture *and* a maximum architecture
    {
        std::vector<int> start_arc{1,2,2,1};
        std::vector<int> max_arc_that_works{1,8,8,1};
        std::vector<int> max_arc_too_few_nodes{1,1,1,1};
        std::vector<int> max_arc_too_many_layers{1,8,8,8,1};
        std::vector<int> max_arc_too_few_layers{1,8,1};
        std::vector<int> max_arc_wrong_input{2,8,8,1};
        std::vector<int> max_arc_wrong_output{1,8,8,2};

        auto pars = net_param();
        pars.net_arc = start_arc;
        pars.max_arc = max_arc_that_works;

        bool exception_thrown = false;

        network n{net_param{}};

        try{
        n = network{pars};
        }
        catch(int exc){
          exception_thrown = true;
        }

        assert(exception_thrown == false);
        assert(n.get_current_arc() == start_arc);
        assert(n.get_max_arc() == max_arc_that_works);

        exception_thrown = false;
        pars.max_arc = max_arc_too_few_nodes;
        try{
        n = network{pars};
        }
        catch(int exc){
          exception_thrown = true;
        }
        assert(exception_thrown == true);

        exception_thrown = false;
        pars.max_arc = max_arc_too_many_layers;
        try{
        n = network{pars};
        }
        catch(int exc){
          exception_thrown = true;
        }
        assert(exception_thrown == true);

        exception_thrown = false;
        pars.max_arc = max_arc_too_few_layers;
        try{
        n = network{pars};
        }
        catch(int exc){
          exception_thrown = true;
        }
        assert(exception_thrown == true);

        exception_thrown = false;
        pars.max_arc = max_arc_wrong_input;
        try{
        n = network{pars};
        }
        catch(int exc){
          exception_thrown = true;
        }
        assert(exception_thrown == true);

        exception_thrown = false;
        pars.max_arc = max_arc_wrong_output;
        try{
        n = network{pars};
        }
        catch(int exc){
          exception_thrown = true;
        }
        assert(exception_thrown == true);
    }
#endif

}
#endif
