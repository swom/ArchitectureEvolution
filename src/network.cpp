#include "network.h"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <numeric>
bool operator==(const net_param& lhs, const net_param& rhs)
{
    bool net_arcs = lhs.net_arc == rhs.net_arc;
    bool max_arcs = lhs.max_arc == rhs.max_arc;
    return max_arcs && net_arcs;
}

template<class Net>
std::vector<weight> register_n_weight_mutations(Net n, double mut_rate, double mut_step, std::mt19937_64& rng, int repeats)
{
    std::vector<weight> networks_weights;
    for(int i = 0; i != repeats; i++)
    {

        auto n_new = n;
        mutate_weights(n_new, mut_rate, mut_step, rng);

        for(auto& layer : n_new.get_net_weights())
            for(auto& node : layer)
                for(size_t j = 0; j != node.get_vec_weights().size(); ++j)
                {
                    if(node.get_vec_weights()[j].get_weight() < -0.0001 || node.get_vec_weights()[j].get_weight() > 0.0001)
                        networks_weights.push_back(node.get_vec_weights()[j]);
                }
    }
    return  networks_weights;
}

template<class Net>
std::vector<weight> register_n_activation_mutations(Net n, double mut_rate, std::mt19937_64 &rng, int repeats)
{
    std::vector<weight> networks_weights;
    for(int i = 0; i != repeats; i++)
    {
        auto n_new = n;
        mutate_activation(n_new, mut_rate, rng);
        auto weights = n_new.get_net_weights();

        for(auto& layer : n_new.get_net_weights())
            for(auto& node : layer)
                for(size_t j = 0; j != node.get_vec_weights().size(); ++j)
                {
                    networks_weights.push_back(node.get_vec_weights()[j]);
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

template<class Net>
bool all_weigths_are_active(const Net &n)
{
    auto weights = n.get_net_weights();

    for(auto &layer : weights ){
        for(auto &node : layer){
            for (size_t i = 0; i != node.get_vec_weights().size(); ++i){
                weight current_weight = node.get_vec_weights()[i];
                if(!current_weight.is_active()){
                    return false;
                }
            }
        }
    }
    return true;
}

template<class Net>
bool all_weigths_have_value(const Net &n, double value)
{
    auto weights = n.get_net_weights();

    for(auto &layer : weights ){
        for(auto &node : layer){
            for (size_t i = 0; i != node.get_vec_weights().size(); ++i){
                weight current_weight = node.get_vec_weights()[i];
                if(current_weight.get_weight() != value)
                {
                    return false;
                }
            }
        }
    }
    return true;
}

template<class Net>
bool on_average_an_nth_of_the_weights_are_inactive(const Net &n,
                                                   const std::vector<weight>&registered_mutations,
                                                   const double &proportion,
                                                   int repeats)
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

template<class Net>
int get_number_weights(const Net &n)
{
    size_t number_weights = 0;
    for(const auto &layer : n.get_net_weights() ){
        for(const auto &node : layer){
            number_weights += node.get_vec_weights().size();
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
        net_param n_p{{1,2,1}, linear};
        network n{n_p};
        n = change_all_weights_values_and_activations(n,1);
        std::vector<double> input{1};
        auto using_member_function = output(n,{input});
        auto using_given_function = output(n,input, &linear);
        auto using_different_given_function = output(n, input, &sigmoid);
        assert(using_given_function == using_member_function);
        assert(using_different_given_function != using_member_function);
    }

    //#ifdef FIX_ISSUE_46
    /// A network can be initilized with a parameter struct net_param
    {
        std::vector<int> net_arc{1, 2, 3, 1} ;
        std::function<double(double)> function = linear;
        net_param n_p{net_arc, function, net_arc};
        network n{n_p};
        //Set weigths to one
        n = change_all_weights_values_and_activations(n, 1);
        //Check architecture
        assert(n.get_net_weights().size() == net_arc.size() -1);
        for(size_t  i = 0; i != n.get_net_weights().size(); i++)
        {
            assert (static_cast<int>(n.get_net_weights()[i].size()) == net_arc[i + 1]);
        }
        //Check activation func
        std::vector<double> input{1};
        auto using_member_function = output(n,{input});
        auto using_given_function = output(n,input, &linear);
        auto using_different_given_function = output(n, input, &sigmoid);
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
        net_param n_p{very_simple_nodes};
        network n{n_p};
        auto net_output = output(n, input, &linear);
        assert(net_output == expected_output);

        //Let's change the weights of the network to something else than 0(e.g 1)
        double new_weight_value = 1.0;
        n = change_all_weights_values_and_activations(n,new_weight_value);
        expected_output = std::vector<double>{2};
        net_output = output(n,input, &linear);
        assert(net_output == expected_output);

        //Testing a more complex arhitecture
        std::vector<int> not_too_simple_nodes{1,3,3,1};
        net_param n_p_complex{not_too_simple_nodes, linear, not_too_simple_nodes};
        network not_too_simple{n_p_complex};
        not_too_simple = change_all_weights_values_and_activations(not_too_simple, new_weight_value);
        expected_output = {1 * 3 * 3};
        net_output = output(not_too_simple, input, &linear);
        assert(net_output == expected_output);

    }

    ///Network weights mutate following a normal distribution
    {
        double mut_rate = 1;
        double mut_step = 0.1;
        std::mt19937_64 rng{};

        auto expected_mean_value  = 0;
        auto expected_stdev = mut_step;

        int repeats = 1000;

        auto very_simple_nodes = std::vector<int>{1,1};
        net_param n_p{very_simple_nodes, linear, very_simple_nodes};
        network n{n_p};

        std::vector<weight> networks_weights = register_n_weight_mutations(n,
                                                                           mut_rate,
                                                                           mut_step,
                                                                           rng,
                                                                           repeats);

        std::vector<double> weights_as_doubles = convert_to_double(networks_weights);

        double mean = calc_mean(weights_as_doubles);
        double stdev = calc_stdev(weights_as_doubles);

        auto delta_mean = mean - expected_mean_value;
        auto delta_stdev = stdev - expected_stdev;
        assert(delta_mean < 0.01 && delta_mean > -0.01);
        assert(delta_stdev < 0.01 && delta_stdev > -0.01);

    }

//    ///A network can be initialized with one bias per node
//    /// stored by layer and by node in a vector of vectors
//    {
//        std::vector<int> net_arch{1,2,2,3,1};
//        network n{net_arch};
//        for(size_t i = 1; i != net_arch.size(); i++)
//        {
//            auto n_biases = n.get_biases()[i - 1].size();
//            auto n_nodes = static_cast<size_t>(net_arch[i]);
//            assert(n_biases == n_nodes);
//        }

//    }
    //This is no longer relevant as biases are stored in nodes

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
        network n{net_param{net_arch, linear, net_arch}};

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
    /// A network contains a vector of vectors of nodes
    {
        net_param n_p{};
        network n{n_p};

        std::vector<std::vector<node>> nodes = n.get_net_weights();
    }
#endif

#define FIX_ISSUE_96
#ifdef FIX_ISSUE_96
    /// The output function works as intended when some connections are switched off
    {
        auto very_simple_nodes = std::vector<int>{1,2,1};
        auto input = std::vector<double>{1};

        //Comparing a default network with weights = 0 to a network with other weights, all switched off
        network n_default{very_simple_nodes};
        auto output_default = output(n_default, input, &linear);

        weight new_weight_value{1.0, false};
        network n_changed = change_all_weights_values_and_activations(n_default, new_weight_value);
        auto output_changed = output(n_changed,input, &linear);

        assert(output_default == output_changed);
    }
#endif

#define FIX_ISSUE_112
#ifdef FIX_ISSUE_112
    {
        network<mutation_type::activation> n_activation{net_param()};
        assert(all_weigths_are_active(n_activation));

        network<mutation_type::weights> n_weights{net_param()};
        assert(all_weigths_have_value(n_weights, 0));

        auto mutation_rate = 1;
        auto mutation_step = 1;
        std::mt19937_64 rng;
        auto rng_copy = rng;

        auto before_mutation = n_weights;
        n_weights.mutate(mutation_rate, mutation_step, rng_copy);
        n_activation.mutate(mutation_rate, mutation_step, rng, mutation_rate);

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
        network<mutation_type::activation> n_activation{pars};
        network<mutation_type::weights> n_weights{pars};
        network n{pars};

        assert(are_equal_except_mutation_type(n, n_activation));
        assert(n == n_weights);
    }
#endif

#define FIX_ISSUE_159
#ifdef FIX_ISSUE_159
  ///Network has a (current) network architecture *and* a maximum architecture
    {
        std::vector<int> start_arc{1,2,2,1};
        std::vector<int> max_arc{1,8,8,1};

        auto pars = net_param();
        pars.net_arc = start_arc;
        pars.max_arc = max_arc;

        network n{pars};


        assert(n.get_current_arc() == start_arc);
        assert(n.get_max_arc() == max_arc);
    }
#endif


#define FIX_ISSUE_187
#ifdef FIX_ISSUE_187
  ///There is a way to change the current architecture to another architecture, as long as it is compatible with the max_architecture
    {
        network n{net_param{}};
        bool exception_thrown = false;
        std::vector<int> new_arc{1,5,1};
        std::vector<int> invalid_new_arc{1,10,1}; //default max_arc is 1-8-1

        try{
        n.change_network_arc(new_arc);
        }
        catch(const std::invalid_argument& e){
            std::cerr << e.what() << std::endl;
            exception_thrown = true;
          }
        assert(exception_thrown == false);
        assert(n.get_current_arc() == new_arc);

        try{
        n.change_network_arc(invalid_new_arc);
        }
        catch(const std::invalid_argument& e){
            std::cerr << e.what() << std::endl;
            exception_thrown = true;
          }
        assert(exception_thrown == true);
    }
#endif

#define FIX_ISSUE_198
#ifdef FIX_ISSUE_198
  ///A node can be duplicated into an inactive node
    {
        net_param n_p{};
        n_p.net_arc = {1,2,1};
        n_p.max_arc = {1,3,1};

        network n{n_p};

        std::mt19937_64 rng;

        mutate_weights(n, 1, 0.1, rng); //so the weights will be different

        network n_before = n;

        node &to_duplicate = n.get_net_weights()[0][0];

        assert(*n.get_empty_node_in_layer(0) == n.get_net_weights()[0][2]); //This should be in third position (index 2)
        auto empty_node_iterator = n.get_empty_node_in_layer(0);

        assert(to_duplicate != *empty_node_iterator);

        n.duplicate_node(to_duplicate, 0, 0, empty_node_iterator);

        ///Checking that the node has been copied - but not everywhere!
        assert(n != n_before);
        assert(to_duplicate == *empty_node_iterator);
        assert(to_duplicate != n.get_net_weights()[0][1]);

        ///Checking that outgoing weights have been copied as well
        const node &second_l_node = n.get_net_weights()[1][0];
        assert(second_l_node.get_vec_weights()[0] == second_l_node.get_vec_weights()[2]);
        assert(second_l_node.get_vec_weights()[0] != second_l_node.get_vec_weights()[1]);

        ///Current architecture has been updated; since it is the max arc, further duplication does nothing
        assert(n.get_current_arc() == n.get_max_arc());
        network n_after_one = n;
        n.duplicate_node(to_duplicate, 0, 0, n.get_empty_node_in_layer(0));
        assert(n == n_after_one);
    }
#endif

#define FIX_ISSUE_201
#ifdef FIX_ISSUE_201
    ///There is a new mutation mode where nodes can be duplicated by mutation
    {
        network<mutation_type::weights_and_activation> n_simple{net_param()};
        std::vector<int> basic_arc{1,2,1};
        assert(n_simple.get_current_arc() == basic_arc);

        network<mutation_type::duplication> n_dup{net_param()};
        assert(n_dup.get_current_arc() == basic_arc);

        auto mutation_rate = 1;
        auto mutation_step = 1;
        std::mt19937_64 rng;
        auto rng_copy = rng;

        n_simple.mutate(mutation_rate, mutation_step, rng_copy, mutation_rate, mutation_rate);
        n_dup.mutate(mutation_rate, mutation_step, rng, mutation_rate, mutation_rate);

        for(size_t i = 0; i != 2; ++i){
            for(int j = 0; j != basic_arc[i+1]; ++j){
                for(int k = 0; k != basic_arc[i]; ++k){
                 assert(n_simple.get_net_weights()[i][j].get_vec_weights()[k] == n_dup.get_net_weights()[i][j].get_vec_weights()[k]);
                }
            }
        }
        assert(n_simple.get_net_weights() != n_dup.get_net_weights());
        assert(n_simple.get_current_arc() != n_dup.get_current_arc());

        assert(n_dup.get_net_weights()[0][2] == n_dup.get_net_weights()[0][1]);
        assert(n_dup.get_net_weights()[0][3] == n_dup.get_net_weights()[0][0]);

        std::vector<int> arc_theo{1,4,1};

        assert(n_dup.get_current_arc() == arc_theo);
    }
#endif

#define FIX_ISSUE_202
#ifdef FIX_ISSUE_202
    ///Inactive nodes are not counted in the response function
    {
        network n1{net_param({1,2,1}, linear, {1,8,1})};
        network n2{net_param({1,2,1}, linear, {1,2,1})};

        n1 = change_all_weights_values(n1, 0.5);
        n2 = change_all_weights_values(n2, 0.5);

        std::vector<double> output1 = output(n1, {1});
        std::vector<double> output2 = output(n2, {1});
        assert(output1 == output2);
    }
#endif

#define FIX_ISSUE_205
#ifdef FIX_ISSUE_205
    ///Connections linked to inactive nodes don't mutate
    {
        network n_before{net_param({1,1,1}, linear, {1,2,1})};
        std::mt19937_64 rng;

        assert(!n_before.get_net_weights()[0][1].is_active());

        network n_w = n_before;
        mutate_weights(n_w, 1, 0.1, rng);

        network n_a = n_before;
        mutate_activation(n_a, 1, rng);

        network n_b = n_before;
        mutate_biases(n_b, 1, 0.1, rng);


        assert(n_w.get_net_weights()[0][0] != n_before.get_net_weights()[0][0]);
        assert(n_w.get_net_weights()[0][1] == n_before.get_net_weights()[0][1]);
        assert(n_w.get_net_weights()[1][0].get_vec_weights()[0] !=
               n_before.get_net_weights()[1][0].get_vec_weights()[0]);
        assert(n_w.get_net_weights()[1][0].get_vec_weights()[1] ==
               n_before.get_net_weights()[1][0].get_vec_weights()[1]);

        assert(n_b.get_net_weights()[0][0] != n_before.get_net_weights()[0][0]);
        assert(n_b.get_net_weights()[0][1] == n_before.get_net_weights()[0][1]);

        assert(n_a.get_net_weights()[0][0] != n_before.get_net_weights()[0][0]);
        assert(n_a.get_net_weights()[0][1] == n_before.get_net_weights()[0][1]);
        assert(n_a.get_net_weights()[1][0].get_vec_weights()[0] !=
               n_before.get_net_weights()[1][0].get_vec_weights()[0]);
        assert(n_a.get_net_weights()[1][0].get_vec_weights()[1] ==
               n_before.get_net_weights()[1][0].get_vec_weights()[1]);

    }
#endif


#define FIX_ISSUE_209
#ifdef FIX_ISSUE_209
    ///There is a function that randomy adds a node to the network with the correct number of connections
    /// With a number of edges depending on the average degree
    {
    net_param n_p{};
    n_p.net_arc = {1,2,1};
    n_p.max_arc = {1,3,1};

    network n{n_p};

    std::mt19937_64 rng;

    assert(*n.get_empty_node_in_layer(0) == n.get_net_weights()[0][2]); //This should be in third position (index 2)
    auto empty_node_iterator = n.get_empty_node_in_layer(0);

    n.add_node(0, empty_node_iterator, rng);

    auto added_node = *empty_node_iterator;

    ///Checking that the node is now active
    assert(added_node.is_active());

    ///Checking that it has the right number of active incoming connections
    size_t n_in_con = 0;
    for(const auto &con : added_node.get_vec_weights())
      if(con.is_active()) ++n_in_con;

    assert(n_in_con == std::round(average_number_incoming_weights(n, 0)));

    ///Checking that it has the right number of active outgoing connections
    size_t n_out_con = 0;
    for(const auto &node : n.get_net_weights()[1])
      if(node.get_vec_weights()[2].is_active()) ++n_out_con;

    assert(n_out_con == std::round(average_number_outgoing_weights(n, 0))); //0 corresponds to the layer

    }
#endif

#define FIX_ISSUE_210
#ifdef FIX_ISSUE_210
    ///When adding a node to a network randomly, which nodes it is connected to is random
    {
    net_param n_p{};
    n_p.net_arc = {1,8,1,8,1};
    n_p.max_arc = {1,8,2,8,1};

    network n1{n_p};
    std::mt19937_64 rng1(0);
    std::mt19937_64 rng2(1);

    for(int i=0; i!= 100; ++i){
        mutate_activation(n1, 0.2, rng1);
      }


    network n2 = n1;

    auto empty_node_iterator1 = n1.get_empty_node_in_layer(1);
    auto empty_node_iterator2 = n2.get_empty_node_in_layer(1);

    n1.add_node(1, empty_node_iterator1, rng1);
    n2.add_node(1, empty_node_iterator2, rng2);

    auto added_node1 = *empty_node_iterator1;
    auto added_node2 = *empty_node_iterator2;

    ///Checking that the node is now active in both cases
    assert(added_node1.is_active() && added_node2.is_active());

    ///Checking that the incoming connections are different in activation
    int dif = 0;
    for(size_t i=0; i != added_node1.get_vec_weights().size(); ++i){
        weight w1 = added_node1.get_vec_weights()[i];
        weight w2 = added_node2.get_vec_weights()[i];
        if(w1.is_active() != w2.is_active()){
            ++dif;
          }
      }
    assert(dif != 0);

    ///Checking that the outgoing connections are different in activation

    dif = 0;
    for(size_t i=0; i != n1.get_net_weights()[2].size(); ++i){
        node node1 = n1.get_net_weights()[2][i];
        node node2 = n2.get_net_weights()[2][i];

        if(node1.get_vec_weights()[1].is_active() != node2.get_vec_weights()[1].is_active()){
            ++dif;
          }
      }
    assert(dif != 0);

    }
#endif

#define FIX_ISSUE_211
#ifdef FIX_ISSUE_211
    ///When adding a new node to a network randomly, the value of the weights and of its bias are chosen randomly
    {
    net_param n_p{};
    n_p.net_arc = {1,2,1};
    n_p.max_arc = {1,3,1};
    std::mt19937_64 rng;

    network n1{n_p};
    mutate_biases(n1, 1, 0.5, rng);
    mutate_weights(n1, 1, 0.5, rng);

    network n2 = n1;
    std::mt19937_64 rng1(0);
    std::mt19937_64 rng2(1);

    auto empty_node_iterator1 = n1.get_empty_node_in_layer(0);
    auto empty_node_iterator2 = n2.get_empty_node_in_layer(0);

    n1.add_node(0, empty_node_iterator1, rng1);
    n2.add_node(0, empty_node_iterator2, rng2);

    auto added_node1 = *empty_node_iterator1;
    auto added_node2 = *empty_node_iterator2;

    ///Checking that the node is now active in both cases
    assert(added_node1.is_active() && added_node1.is_active());

    ///Checking that the incoming connections are different in weight value
    weight w1 = added_node1.get_vec_weights()[0];
    weight w2 = added_node2.get_vec_weights()[0];
    assert(w1.get_weight() != w2.get_weight());

    ///Checking that the outgoing connections are different in weight value
    node node1 = n1.get_net_weights()[1][0];
    node node2 = n2.get_net_weights()[1][0];
    assert(node1.get_vec_weights()[2].get_weight() != node2.get_vec_weights()[2].get_weight());

    ///Checking that the bias value of the nodes is different
    assert(added_node1.get_bias() != added_node2.get_bias());
    }
  #endif

#define FIX_ISSUE_212
#ifdef FIX_ISSUE_212
    ///There is a new mutation mode where nodes can randomly added by mutation
    {
        net_param n_p{{1,1,1}};
        network<mutation_type::addition> n_mut{n_p};
        network<mutation_type::duplication> n_dup{n_p};
        network n_add{n_p};

        auto mutation_rate = 1;
        std::mt19937_64 rng;
        auto rng_copy = rng;
        auto rng_copy_2 = rng;

        n_mut.get_net_weights()[0][2].change_bias(1); //to have some randomness in the values generated
        n_dup.get_net_weights()[0][2].change_bias(1);
        n_add.get_net_weights()[0][2].change_bias(1);//so that there is difference between dup and add

        rng.discard(1); //to compensate for the rng being called to know if there is mutation in mutate

        n_mut.mutate(0, 0, rng_copy, 0, mutation_rate);
        n_add.add_node(0, n_add.get_empty_node_in_layer(0), rng);
        n_dup.mutate(0, 0, rng_copy_2, 0, mutation_rate);

        assert(are_equal_except_mutation_type(n_mut, n_add));
        assert(!are_equal_except_mutation_type(n_mut, n_dup));
    }
#endif

#define FIX_ISSUE_231
#ifdef FIX_ISSUE_231
  ///There is a function that deletes a given node
  {
  net_param n_p{};
  n_p.net_arc = {1,2,1};

  network n{n_p};

  n.delete_node(0, 0); //Deleting the first node of the middle layer

  auto deleted_node = n.get_net_weights()[0][0];

  ///Checking that the node has been inactivated, its bias put back to 0 and all its weights reset
  assert(!deleted_node.is_active());
  assert(deleted_node.get_bias() == 0);
  weight w_theo{};
  for(const auto &weight: deleted_node.get_vec_weights()){
      assert(weight == w_theo);
    }


  ///Checking that outgoing weights have also been reset
  for(const auto &node : n.get_net_weights()[1]){
      assert(node.get_vec_weights()[0] == w_theo);
    }
  }
#endif

#define FIX_ISSUE_234
#ifdef FIX_ISSUE_234
    ///There is a new mutation function where nodes can randomly be deleted by mutation
    {
        net_param n_p{{1,2,1}};
        network n_mut{n_p};
        network n_del{n_p};

        auto mutation_rate = 1;
        std::mt19937_64 rng;


        mut_del(n_mut, mutation_rate, rng);
        n_del.delete_node(0, 0);

        assert(n_mut == n_del);
    }
#endif

#define FIX_ISSUE_235
#ifdef FIX_ISSUE_235
    ///There is a non-ratchet duplication mutation mode
    /// ///There is a non-ratchet addition mutation mode
    {
        net_param n_p{{1,3,1}};

        network<mutation_type::NRduplication> n_nrd_mut{n_p};
        network n_dupdel{n_p};

        auto mutation_rate = 0.5;
        std::mt19937_64 rng;
        std::mt19937_64 rng_copy = rng;

        n_nrd_mut.mutate(mutation_rate, 0.1, rng, mutation_rate, mutation_rate);

        mutate_biases(n_dupdel, mutation_rate, 0.1, rng_copy);
        mutate_activation(n_dupdel, mutation_rate, rng_copy);
        mutate_weights(n_dupdel, mutation_rate, 0.1, rng_copy);
        mut_dupl_node(n_dupdel, mutation_rate, rng_copy);
        mut_del(n_dupdel, mutation_rate, rng_copy);

        assert(are_equal_except_mutation_type(n_nrd_mut, n_dupdel));

        network<mutation_type::NRaddition> n_nra_mut{n_p};
        network n_addel{n_p};

        n_nra_mut.mutate(mutation_rate, 0.1, rng, mutation_rate, mutation_rate);

        mutate_biases(n_addel, mutation_rate, 0.1, rng_copy);
        mutate_activation(n_addel, mutation_rate, rng_copy);
        mutate_weights(n_addel, mutation_rate, 0.1, rng_copy);
        mut_add_node(n_addel, mutation_rate, rng_copy);
        mut_del(n_addel, mutation_rate, rng_copy);
    }
#endif

#define FIX_ISSUE_226
#ifdef FIX_ISSUE_226
    ///The range from which new nodes bias and weights are drawn during random addition
    /// depend on min and max of existing values
    {
    net_param n_p{};
    n_p.net_arc = {1,7,1};
    network n{n_p};
    std::mt19937_64 rng;
    mutate_biases(n, 1, 0.5, rng);
    mutate_weights(n, 1, 0.5, rng);

    for(int i = 0; i != 100; ++i){
        network n_copy = n;
        std::mt19937_64 rng_2(i);

        auto empty_node_iterator = n_copy.get_empty_node_in_layer(0);
        n_copy.add_node(0, empty_node_iterator, rng_2);

        const auto &node = n_copy.get_net_weights()[0][7];
        assert(node.get_bias() > min_bias_in_layer(n_copy, 0));
        assert(node.get_bias() < max_bias_in_layer(n_copy, 0));

        for(const auto &weight : node.get_vec_weights()){
            assert(weight.get_weight() > min_weight_in_layer(n_copy, 0));
            assert(weight.get_weight() < max_weight_in_layer(n_copy, 0));
          }

        const auto &node_next_l = n_copy.get_net_weights()[1][0];
        const auto &outgoing_weight = node_next_l.get_vec_weights()[7];
        assert(outgoing_weight.get_weight() > min_weight_in_layer(n_copy, 1));
        assert(outgoing_weight.get_weight() < max_weight_in_layer(n_copy, 1));
      }
    }
#endif

#define FIX_ISSUE_239
#ifdef FIX_ISSUE_239
    ///Stochastic duplication should only change connections to and from active nodes
    /// when calculating the average number of active ingoing/ outgoing connections,
    /// this should only take into account those to and from active nodes.
    {
    net_param n_p{};
    n_p.net_arc = {1,2,2,2,1};
    n_p.max_arc = {1,3,3,3,1};
    network n{n_p};
    std::mt19937_64 rng;

    //Set up so that the active nodes each have one deactivated connection

    for(size_t i = 1; i != 3; ++i){
        for(size_t j = 0; j !=2; ++j){
            auto &node = n.get_net_weights()[i][j];
            weight w(0, false);
            node.change_nth_weight(w, 0);
          }
      }

    for(int i = 0; i != 100; ++i){
    auto n_copy = n;
    auto empty_node_iterator = n_copy.get_empty_node_in_layer(1);
    n_copy.add_node(1, empty_node_iterator, rng);
    weight w{};

    auto added_node = n_copy.get_net_weights()[1][2];
    assert(added_node.get_vec_weights()[2] == w);

    auto empty_node = n_copy.get_net_weights()[2][2];
    assert(empty_node.get_vec_weights()[2] == w);
      }
    }
#endif

    ///Networks can be templated to be plastic and have thus one extra input
    {
        int n_inputs = 1;
        int n_outputs = 1;

        net_param n_p{};
        n_p.net_arc = {n_inputs,n_outputs};
        n_p.max_arc = {n_inputs,n_outputs};

        network<mutation_type::weights, response_type::plastic> n_plastic{n_p};
        network<mutation_type::weights, response_type::constitutive> n_constitutive{n_p};

        assert(n_plastic.get_max_arc()[0] == n_constitutive.get_max_arc()[0] + 1);
        assert(n_plastic.get_current_arc()[0] == n_constitutive.get_current_arc()[0] + 1);

        //check that nodes in first layer have one more weight
        ///that connects to extra input in plastic network

        auto  nodes_in_constitutive = n_constitutive.get_net_weights()[0];
        auto  nodes_in_plastic = n_plastic.get_net_weights()[0];

        assert(nodes_in_plastic.size() == nodes_in_constitutive.size());
        for(int node = 0; node !=  n_plastic.get_net_weights()[0].size(); node++)
        {
            assert(nodes_in_constitutive[node].get_vec_weights().size() + 1 == nodes_in_plastic[node].get_vec_weights().size());
        }
    }

}
#endif
