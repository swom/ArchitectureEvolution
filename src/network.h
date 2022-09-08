#ifndef NETWORK_H
#define NETWORK_H
#include "utilities.h"
#include <iostream>
#include "node.h"
#include "mutation_type.h"
#include "rndutils.hpp"
#include "response_type.h"
#include "netwrok_spectrum.h"
#include <random>
#include <ranges>
#include <mutex>
#include <nlohmann/json.hpp>

double sigmoid(double x);
double linear(double x);
double constant_one(const std::vector<double>&);
double constant_zero(const std::vector<double>&);


static std::map<std::string, std::function<double(double)>> string_to_act_func_map
{
{"linear", linear},
{"sigmoid", sigmoid}
                                                            };

struct net_param
{
    NLOHMANN_DEFINE_TYPE_INTRUSIVE(net_param,
                                   net_arc,
                                   max_arc,
                                   resp_type,
                                   input_range,
                                   n_sampled_inputs
                                   )

    net_param(const std::vector<int>& net_arch = {1,2,1},
              std::function<double(double)> func = linear,
              const std::vector<int>& max_arch = {1,8,1},
              response_type response_type = response_type::constitutive,
              const std::vector<double>& inp_range = {0,0},
              const int& n_sampl_inps = 0):
        net_arc{net_arch},
        function{func},
        max_arc{max_arch},
        resp_type{response_type},
        input_range{inp_range.front(), inp_range.back()},
        n_sampled_inputs{n_sampl_inps}
    {};

    std::vector<int> net_arc;
    std::function<double(double)> function;
    std::vector<int> max_arc;
    response_type resp_type;
    range input_range;
    int n_sampled_inputs;
};

bool operator==(const net_param& lhs, const net_param& rhs);

///Mutates the weights of a network
template<class Net>
void mutate_weights(Net& n, const double& mut_rate,
                    const double& mut_step,
                    std::mt19937_64& rng)
{
    if(mut_rate != 0)
    {
        std::bernoulli_distribution mut_p{mut_rate};
        std::normal_distribution<double> mut_st{0,mut_step};

        for(size_t layer = 0; layer != n.get_net_weights().size(); ++layer)
        {
            auto& current_layer = n.get_net_weights()[layer];
            for(size_t node = 0; node != current_layer.size(); ++node)
            {
                auto &current_node = current_layer[node];
                if(current_node.is_active())
                {
                    for(size_t weights = 0; weights != current_node.get_vec_weights().size(); ++weights)
                    {
                        if(mut_p(rng) &&
                                (layer == 0 ? true : n.get_net_weights()[layer -1][weights].is_active()))
                        {
                            const auto &current_weight = current_node.get_vec_weights().at(weights);
                            weight mutated_weight(current_weight.get_weight() + mut_st(rng),
                                                  current_weight.is_active());
                            current_node.change_nth_weight(mutated_weight, weights);
                        }
                    }
                }
            }
        }
    }
}

///Mutates the nodes of a network via duplication
template<class Net>
void mut_dupl_node(Net& n,
                   const double& mut_rate,
                   std::mt19937_64& rng)
{

    std::bernoulli_distribution mut_p{mut_rate};

    for(size_t layer = 0; layer != n.get_current_arc().size() - 1; layer++)
    {
        auto& current_layer = n.get_net_weights().at(layer);

        // to avoid chain duplications
        for(int node = int(current_layer.size() - 1); node >= 0; node--)
        {
            const auto& current_node = current_layer.at(node);

            if(current_node.is_active() && mut_p(rng))
            {
                //this returns an iterator if you use std::find()
                auto free_node = n.get_empty_node_in_layer(layer);

                if(free_node != current_layer.end()){
                    n.duplicate_node(current_node, layer, node, free_node);
                }
            }
        }
    }

}

///Mutates a network via random addition of nodes
template<class Net>
void mut_add_node(Net& n,
                  const double& mut_rate,
                  std::mt19937_64& rng)
{
    if(mut_rate){
        std::bernoulli_distribution mut_p{mut_rate};

        for(size_t layer = 0; layer != n.get_current_arc().size() - 1; layer++)
        {
            auto& current_layer = n.get_net_weights().at(layer);

            // to avoid chain duplications
            for(int node = int(current_layer.size() - 1); node >= 0; node--)
            {

                const auto& current_node = current_layer.at(node);

                ///To keep it comparable with duplication, this is also per active node
                /// Networks with more active nodes will have more mutations

                if(current_node.is_active() && mut_p(rng))
                {
                    //this returns an iterator if you use std::find()
                    auto free_node = n.get_empty_node_in_layer(layer);

                    if(free_node != current_layer.end()){
                        n.add_node(layer, free_node, rng);
                    }
                }
            }
        }
    }
}

///Mutates a network by deleting some nodes
template<class Net>
void mut_del(Net& n,
             const double& mut_rate,
             std::mt19937_64& rng)
{

    std::bernoulli_distribution mut_p{mut_rate};

    for(size_t layer = 0; layer != n.get_net_weights().size() - 1; layer++)
    {
        auto& current_layer = n.get_net_weights().at(layer);

        for(size_t node_index = 0; node_index != current_layer.size(); node_index++)
        {
            const auto& current_node = current_layer.at(node_index);

            if(current_node.is_active() && mut_p(rng))
            {
                n.delete_node(layer, node_index);
            }
        }
    }
}

///Mutates the activation of the weights of the network - they get switched on and off
template<class Net>
void mutate_activation(Net &n, const double &mut_rate, std::mt19937_64 &rng)
{
    if(mut_rate != 0){
        std::bernoulli_distribution mut_p{mut_rate};

        for(size_t i = 0; i != n.get_net_weights().size(); i++)
            for(size_t j = 0; j != n.get_net_weights()[i].size(); j++){
                auto &current_node = n.get_net_weights()[i][j];
                if(current_node.is_active()){
                    for(size_t k = 0; k != current_node.get_vec_weights().size(); k++)
                    {
                        if(mut_p(rng) && (i == 0 ? true : n.get_net_weights()[i -1][k].is_active())){
                            const weight &current_weight = current_node.get_vec_weights()[k];
                            weight mutated_weight(current_weight.get_weight(),
                                                  !current_weight.is_active());
                            current_node.change_nth_weight(mutated_weight, k);
                        }
                    }
                }
            }
    }
}

///Mutates the biases of the nodes
template<class Net>
void mutate_biases(Net& n, const double& mut_rate,
                   const double& mut_step,
                   std::mt19937_64& rng)
{
    if(mut_rate != 0){
        std::bernoulli_distribution mut_p{mut_rate};
        std::normal_distribution<double> mut_st{0,mut_step};

        auto& vector = n.get_net_weights();
        for(auto& layer : vector){
            for(auto& node : layer)
            {
                if(mut_p(rng) && node.is_active()){
                    node.change_bias(node.get_bias() + mut_st(rng));
                }
            }
        }
    }
}

///Calculates the reaction_norm that is produced by an optimal function
/// for a given range and a given number of data points
template<class Func>
reac_norm calculate_reaction_norm_from_function(const Func& func,
                                                const range& cue_range,
                                                const int& n_data_points);

template<mutation_type M = mutation_type::weights,
         response_type R = response_type::constitutive>
class network
{
public:
    static constexpr mutation_type mutation_t = M;
    static constexpr response_type response_t = R;

private:

    ///Calculates the reaction norm fo the network after a mutation on a weight
    //make a private function of net
    reac_norm calc_alternative_reac_norm(reac_norm rn,
                                         double new_weight,
                                         const range& input_range,
                                         int n_inputs,
                                         size_t layer_index,
                                         size_t node_index,
                                         size_t weight_index
                                         )
    {
        rn = calculate_reaction_norm_with_modified_weight(*this,
                                                          input_range,
                                                          n_inputs,
                                                          layer_index,
                                                          node_index,
                                                          weight_index,
                                                          new_weight);
        return rn;
    }
public:

    network(std::vector<int> nodes_per_layer,
            std::function<double(double)> activation_function = &linear):
        m_network_weights{},
        m_input_size{nodes_per_layer[0]},
        m_activation_function{activation_function},
        m_current_arc{nodes_per_layer},
        m_max_arc{nodes_per_layer}
    {

        for (size_t i = 1; i != nodes_per_layer.size(); i++ )
        {
            std::vector<node>temp_layer_vector;
            size_t n_nodes_prev_layer = nodes_per_layer[i-1];
            for(int j = 0; j != nodes_per_layer[i]; j++)
            {
                std::vector<weight> temp_weights(n_nodes_prev_layer);
                node temp_node(temp_weights);
                temp_layer_vector.push_back(temp_node);
            }


            //A vector of the size of the number of connections is pushed back in the weight matrix
            m_network_weights.push_back(temp_layer_vector);

        }
    }


    network(const net_param &n_p):
        m_network_weights{},
        m_input_size{n_p.net_arc[0]},
        m_activation_function{n_p.function},
        m_current_arc{n_p.net_arc},
        m_max_arc{n_p.max_arc}
    {

        if constexpr (R == response_type::additive)
        {
            m_additive_genes = calculate_reaction_norm_from_function(constant_zero,
                                                                     n_p.input_range,
                                                                     n_p.n_sampled_inputs);
            return;
        }
        ///Change architecture by adding one extra input
        ///for environment if response is plastic
        if constexpr (R == response_type::plastic)
        {
            m_input_size ++;
            m_current_arc[0] ++;
            m_max_arc[0] ++;
        }

        if(!net_arc_and_max_arc_are_compatible(m_current_arc, m_max_arc)){
            throw std::runtime_error{"starting and maximum architecture are not compatible"};
        }

        for (size_t layer = 1; layer != m_current_arc.size(); layer++ )
        {
            std::vector<node>temp_layer_vector;
            size_t n_nodes_prev_layer = m_max_arc[layer-1];
            for(int node = 0; node != m_max_arc[layer]; node++)
            {
                std::vector<weight> temp_weights(n_nodes_prev_layer);

                class node temp_node(temp_weights);
                if(node < m_current_arc[layer])
                {
                    temp_node.activate();
                }

                temp_layer_vector.push_back(temp_node);
            }

            //A vector of the size of the number of connections is pushed back in the weight matrix
            m_network_weights.push_back(temp_layer_vector);
        }
    }


    ///Checks if hhe layer has no nodes
    inline bool any_layer_has_no_nodes() const noexcept{ return std::any_of(
                    m_network_weights.begin(),
                    m_network_weights.end(),
                    [] (const auto& i) {return i.size() == 0;});}

    ///Changes the value of all weights and biases in the network to a certain value
    void change_all_weights_and_biases_values(double new_value)
    {
        for(auto& layer : m_network_weights)
            for(auto& node : layer)
            {
                node.change_bias(new_value);
                for(size_t i = 0; i != node.get_vec_weights().size(); ++i)
                {
                    const weight &current_weight = node.get_vec_weights()[i];
                    weight changed_weight(new_value, current_weight.is_active());
                    node.change_nth_weight(changed_weight, i);
                }
            }
    }

    ///Changes the value of all weights in the network to a certain value
    void change_all_weights_values(double new_weight)
    {
        for(auto& layer : m_network_weights)
            for(auto& node : layer)
                for(size_t i = 0; i != node.get_vec_weights().size(); ++i)
                {
                    const weight &current_weight = node.get_vec_weights()[i];
                    weight changed_weight(new_weight, current_weight.is_active());
                    node.change_nth_weight(changed_weight, i);
                }
    }

    void mutate_genes(const double& mut_rate_weight,
                      const double& mut_step,
                      std::mt19937_64& rng)
    {
        std::bernoulli_distribution p_mut{mut_rate_weight};
        std::normal_distribution p_step{mut_step};
        for(auto& gene : m_additive_genes)
        {
            gene.m_y += p_step(rng) * p_mut(rng);
        }
    }

    void mutate(const double& mut_rate_weight,
                const double& mut_step,
                std::mt19937_64& rng,
                const double& mut_rate_act = 0.01,
                const double& mut_rate_dup = 0.01
            )
    {
        if constexpr(R == response_type::additive)
        {
            mutate_genes(mut_rate_weight, mut_step, rng);
        }

        if constexpr (M == mutation_type::activation)
        {
            mutate_biases(*this, mut_rate_weight, mut_step, rng);
            mutate_activation(*this, mut_rate_act, rng);
        }

        else if constexpr(M == mutation_type::weights)
        {
            mutate_biases(*this, mut_rate_weight, mut_step, rng);
            mutate_weights(*this, mut_rate_weight, mut_step, rng);
        }

        else if constexpr(M == mutation_type::weights_and_activation)
        {
            mutate_biases(*this, mut_rate_weight, mut_step, rng);
            mutate_activation(*this, mut_rate_act, rng);
            mutate_weights(*this, mut_rate_weight, mut_step, rng);
        }

        else if constexpr(M == mutation_type::duplication)
        {
            mutate_biases(*this, mut_rate_weight, mut_step, rng);
            mutate_activation(*this, mut_rate_act, rng);
            mutate_weights(*this, mut_rate_weight, mut_step, rng);
            mut_dupl_node(*this, mut_rate_dup, rng);
        }

        else if constexpr(M == mutation_type::addition)
        {
            mutate_biases(*this, mut_rate_weight, mut_step, rng);
            mutate_activation(*this, mut_rate_act, rng);
            mutate_weights(*this, mut_rate_weight, mut_step, rng);
            mut_add_node(*this, mut_rate_dup, rng);
        }

        else if constexpr(M == mutation_type::NRduplication)
        {
            mutate_biases(*this, mut_rate_weight, mut_step, rng);
            mutate_activation(*this, mut_rate_act, rng);
            mutate_weights(*this, mut_rate_weight, mut_step, rng);
            mut_dupl_node(*this, mut_rate_dup, rng);
            mut_del(*this, mut_rate_dup, rng);
        }

        else if constexpr(M == mutation_type::NRaddition)
        {
            mutate_biases(*this, mut_rate_weight, mut_step, rng);
            mutate_activation(*this, mut_rate_act, rng);
            mutate_weights(*this, mut_rate_weight, mut_step, rng);
            mut_add_node(*this, mut_rate_dup, rng);
            mut_del(*this, mut_rate_dup, rng);
        }

    };

    ///Returns reference to the first node of the first layer
    node& get_first_node(){return m_network_weights[0][0];}


    ///Returns the index of the first inactive node in the layer l
    inline std::vector<node>::iterator get_empty_node_in_layer(size_t l)
    {
        std::vector<node> &layer = get_net_weights()[l];
        return std::find_if(layer.begin(), layer.end(), is_inactive);
    }


    inline void duplicate_node(const node &node_to_duplicate,
                               size_t layer,
                               size_t index_to_duplicate,
                               const std::vector<node>::iterator &empty_node_iterator)
    {
        if (layer == (m_current_arc.size() - 1))
        {
            throw std::runtime_error{"help :(, trying to duplicate one node in output"};
        }

        //This is bad, very confusing
        //Check if layer has free nodes, if not just return
        //Should be already checked in mutate_duplication
        if(m_current_arc.at(layer + 1) == m_max_arc.at(layer + 1))
        {
            return;
        }

        size_t index = empty_node_iterator - get_net_weights().at(layer).begin() ;
        m_network_weights.at(layer).at(index) = node_to_duplicate;

        //This also should trigger if first check fails
        for(auto &node : m_network_weights.at(layer+1)){
            weight weight_to_duplicate = node.get_vec_weights().at(index_to_duplicate);
            node.change_nth_weight(weight_to_duplicate, index);
        }
        ++m_current_arc.at(layer+1);

        if(m_current_arc.at(layer + 1) > m_max_arc.at(layer + 1))
        {
            std::runtime_error{"One node too much has been added"};
        }
    }

    inline void add_node(size_t layer,
                         const std::vector<node>::iterator &empty_node_iterator,
                         std::mt19937_64 &rng)
    {
        if(m_current_arc[layer + 1] == m_max_arc[layer + 1])
        {
            return;
        }

        //Preparing the indexes of which connections to inactivate
        node &node_to_add = *empty_node_iterator;
        std::vector<size_t> vec_indexes(node_to_add.get_vec_weights().size());
        std::iota(std::begin(vec_indexes), std::end(vec_indexes), 0);
        if(layer != 0){
            for(size_t i = 0; i != m_network_weights[layer-1].size(); ++i){
                if(!m_network_weights[layer-1][i].is_active()){
                    vec_indexes[i] = 999;
                }
            }
            vec_indexes.erase(std::remove(vec_indexes.begin(), vec_indexes.end(), 999), vec_indexes.end());
        }

        std::vector<size_t> indexes_to_activate;
        int nb_incoming_weights = int(std::round(average_number_incoming_weights(*this, layer)));
        std::sample(vec_indexes.begin(), vec_indexes.end(), std::back_inserter(indexes_to_activate), nb_incoming_weights, rng);

        std::uniform_real_distribution<double> dist_in(min_weight_in_layer(*this, layer),
                                                       max_weight_in_layer(*this, layer));


        std::vector<size_t> vec_indexes_out(m_network_weights[layer + 1].size());
        std::iota(std::begin(vec_indexes_out), std::end(vec_indexes_out), 0);

        for(size_t i = 0; i != m_network_weights[layer+1].size(); ++i){
            if(!m_network_weights[layer+1][i].is_active()){
                vec_indexes_out[i] = 999;
            }
            vec_indexes_out.erase(std::remove(vec_indexes_out.begin(), vec_indexes_out.end(), 999), vec_indexes_out.end());

        }

        std::vector<size_t> indexes_to_activate_out;
        int nb_outgoing_weights = int(std::round(average_number_outgoing_weights(*this, layer)));
        std::sample(vec_indexes_out.begin(), vec_indexes_out.end(), std::back_inserter(indexes_to_activate_out), nb_outgoing_weights, rng);

        std::uniform_real_distribution<double> dist_out(min_weight_in_layer(*this, layer + 1),
                                                        max_weight_in_layer(*this, layer + 1));

        //activating
        size_t index = empty_node_iterator - get_net_weights()[layer].begin() ;
        node &added_node = m_network_weights[layer][index];
        added_node.activate();

        //Making the bias random within a given range
        std::uniform_real_distribution<double> dist_bias(min_bias_in_layer(*this, layer),
                                                         max_bias_in_layer(*this, layer));
        added_node.change_bias(dist_bias(rng));

        //adding incoming connections
        for(size_t i=0; i!= added_node.get_vec_weights().size(); ++i){
            weight w = added_node.get_vec_weights()[i];
            if(layer != 0){
                if((m_network_weights[layer-1][i].is_active())){
                    if(std::count(indexes_to_activate.begin(), indexes_to_activate.end(), i)){
                        w.change_activation(true);
                        w.change_weight(dist_in(rng));
                        added_node.change_nth_weight(w,i);
                    }
                    else{
                        w.change_activation(false);
                        added_node.change_nth_weight(w,i);
                    }
                }
            }
            else{
                if(std::count(indexes_to_activate.begin(),
                              indexes_to_activate.end(), i))
                {
                    w.change_activation(true);
                    w.change_weight(dist_in(rng));
                    added_node.change_nth_weight(w,i);
                }
                else{
                    w.change_activation(false);
                    added_node.change_nth_weight(w,i);
                }
            }
        }

        //adding outgoing connections
        for(size_t i=0; i!= m_network_weights[layer + 1].size(); ++i){
            weight w = m_network_weights[layer + 1][i].get_vec_weights()[index];
            if((m_network_weights[layer+1][i].is_active())){
                if(std::count(indexes_to_activate_out.begin(), indexes_to_activate_out.end(), i)){
                    w.change_activation(true);
                    w.change_weight(dist_out(rng));
                    m_network_weights[layer + 1][i].change_nth_weight(w, index);
                }
                else{
                    w.change_activation(false);
                    m_network_weights[layer + 1][i].change_nth_weight(w, index);
                }
            }
        }

        ++m_current_arc[layer + 1];
    }

    inline void delete_node(size_t layer, size_t index)
    {
        if(m_current_arc[layer + 1] == 1)
            return;

        //de-activating
        node &deleted_node = m_network_weights.at(layer).at(index);
        deleted_node.deactivate();

        //Resetting bias
        deleted_node.change_bias(0);

        //Resetting incoming connections
        weight default_w{};
        for(size_t i = 0; i != deleted_node.get_vec_weights().size(); i++){
            deleted_node.change_nth_weight(default_w,i);
        }

        //Resetting outgoing connections
        for(size_t i=0; i!= m_network_weights.at(layer+1).size(); i++){
            m_network_weights.at(layer+1).at(i).change_nth_weight(default_w, index);
        }

        --m_current_arc.at(layer+1);

        if(any_layer_has_no_nodes())
        {
            throw std::runtime_error{"One layer got has 0 nodes"};
        }
    }


    NLOHMANN_DEFINE_TYPE_INTRUSIVE(network,
                                   m_input_size,
                                   m_network_weights
                                   );

    ///Returns the activation function
    std::function<double(double)> get_activation_function() const noexcept{return m_activation_function;}

    ///Returns the input size
    size_t get_input_size() const noexcept {return static_cast<size_t>(m_input_size);}

    double operator ()(double n) const {return m_activation_function(n);}

    ///Returns the additive genes
    const reac_norm& get_genes() const noexcept {return m_additive_genes;}

    ///Returns const ref to vector of weights
    const std::vector<std::vector<node>>& get_net_weights() const noexcept{return m_network_weights;}

    ///Returns not constant ref to vector of weights
    std::vector<std::vector<node>>& get_net_weights() noexcept{return m_network_weights;}

    ///Returns the current architecture
    const std::vector<int>& get_current_arc() const noexcept{return m_current_arc;}

    ///Returns the maximum architecture
    const std::vector<int>& get_max_arc() const noexcept{return m_max_arc;}

    ///Changes the current arc of a network to a new one
    void change_network_arc(std::vector<int> new_arc){


        if(net_arc_and_max_arc_are_compatible(new_arc, m_max_arc)){
            m_current_arc = new_arc;
        }
        else
        {
            throw std::invalid_argument{"The current and maximum architectures are not compatible"};
        }


    }


    //make a public function of net
    net_weight_mut_spectrum calc_spectrum_weights_for_weights_mutation(double mut_step,
                                                                       std::mt19937_64& rng,
                                                                       int n_mutations,
                                                                       range input_range,
                                                                       int n_inputs
                                                                       )

    {
        auto mutable_net = *this;
        std::normal_distribution<double> mut_dist(0, mut_step);

        net_weight_mut_spectrum network_weights_spectrum(mutable_net.get_net_weights().size());
        layer_weight_mut_spectrum layer_spectrum;
        node_weight_mut_spectrum node_spectrum;
        weight_mut_spectrum weight_spectrum(n_mutations);
        reac_norm rn;
        rn.reserve(n_inputs);

        std::vector<double> mutations(n_mutations);
        std::generate(mutations.begin(), mutations.end(),
                      [&rng, &mut_dist]{return mut_dist(rng);});

        for(size_t layer_index = 0; layer_index != mutable_net.get_net_weights().size(); layer_index++)
        {
            auto current_layer = mutable_net.get_net_weights()[layer_index];
            layer_spectrum.resize(current_layer.size());

            for(size_t node_index = 0; node_index != current_layer.size(); ++node_index)
            {
                if(!current_layer.at(node_index).is_active()) continue;

                node_spectrum.resize(current_layer[node_index].get_vec_weights().size());

                for(size_t weight_index = 0;
                    weight_index != current_layer.at(node_index).get_vec_weights().size(); //Demetra rule!!!
                    ++weight_index)
                {
                    if(!current_layer[node_index].get_vec_weights()[weight_index].is_active()) continue;
#pragma omp parallel for
                    for (int mut = 0; mut < int(mutations.size()); mut++)
                    {
                        auto new_weight = current_layer[node_index].get_vec_weights()[weight_index].get_weight() + mutations[mut];
                        weight_spectrum[mut] = calc_alternative_reac_norm(rn,
                                                                          new_weight,
                                                                          input_range,
                                                                          n_inputs,
                                                                          layer_index,
                                                                          node_index,
                                                                          weight_index
                                                                          );
                    }
                    node_spectrum[weight_index] = weight_spectrum;
                }
                layer_spectrum[node_index] = node_spectrum;
            }
            network_weights_spectrum[layer_index] = layer_spectrum;
        }
        return network_weights_spectrum;
    }


private:

    ///Additivegnes used if network is templated as response_type::additive
    reac_norm m_additive_genes;

    ///Vector of of vectors, representing the weights coming into each node
    std::vector<std::vector<node>> m_network_weights;

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

///Creates N mutation values given a mutation step S
std::vector<double> create_mutations(int n_mutations,
                                     double mutation_step, std::mt19937_64 &rng);

///Calculates the distance between the y elements of 2 reaction norms
double rn_distance(const reac_norm& lhs, const reac_norm& rhs);

///Calculates the reaction_norm of an individual's network
/// for a given range and a given number of data points
template<class Net>
double calculate_net_distance_from_reaction_norm(const Net& net,
                                                 const reac_norm& r_norm
                                                 )
{
    double sqr_distance = 0;
    std::vector<double> input;
    std::vector<double> output;
    for (int i = 0; i != r_norm.size(); i++)
    {
        input.assign({r_norm[i].m_x});
        output_ugly_but_fast(net, input, output);
        sqr_distance += (r_norm[i].m_y - output[0]) *
                (r_norm[i].m_y - output[0]);
    }
    return std::sqrt(sqr_distance);
}

///Caluclates the distance from a base reaction norm
/// of a netwrok when one of its eights undergoes a mutation
template<class Net>
double calc_rn_distance_for_weight_mut(Net& net,
                                       weight& current_weight,
                                       const reac_norm& base_reac_norm,
                                       const double& mutation)
{

    current_weight.mutate_weight(mutation);
    auto dist = calculate_net_distance_from_reaction_norm(net,
                                                          base_reac_norm);
    current_weight.reverse_mutate_weight(mutation);
    return dist;
}

///Calculates the distances from a base reaction norm
/// of a netwrok when the weights of one of its
/// nodes undergo a mutation
template<class Net>
void  calc_rn_distance_for_weights_mut(Net& net,
                                       node& node,
                                       const reac_norm& reac_norm,
                                       const double& mutation,
                                       std::vector<double>& distance)
{
    for(auto& current_weight : node.get_vec_mutable_weights())
    {
        distance.emplace_back(calc_rn_distance_for_weight_mut(net,
                                                              current_weight,
                                                              reac_norm,
                                                              mutation));
    }
}

///Calculates the distances from a base reaction norm
/// of a netwrok when the weights of one of its
/// nodes undergo a mutation
template<class Net>
void  calc_delta_distances_for_weights_mut(Net& net,
                                           node& node,
                                           const reac_norm& reac_norm,
                                           const double& mutation,
                                           std::vector<double>& distance,
                                           double distance_from_base_norm)
{
    for(auto& current_weight : node.get_vec_mutable_weights())
    {
        distance.emplace_back(distance_from_base_norm - calc_rn_distance_for_weight_mut(net,
                                                                                        current_weight,
                                                                                        reac_norm,
                                                                                        mutation));
    }
}

///Calculates the distance from a base reaction norm
/// of a network whose node's bias is mutated
template<class Net>
double calculate_rn_distance_for_bias_mut(Net& net,
                                          node& node,
                                          const reac_norm& base_reac_norm,
                                          const double& mutation)
{
    node.mutate_bias(mutation);
    auto dist = calculate_net_distance_from_reaction_norm(net, base_reac_norm);
    node.reverse_mutate_bias(mutation);
    return dist;
}

////Calculates the reaction norm of a network with a mutated bias
/// using a scratch reaction norm
template<class Net>
void calculate_rn_for_bias_mut(Net& net,
                               node& node,
                               reac_norm& scratch_reaction_norm,
                               const double& mutation)
{
    node.mutate_bias(mutation);
    recalculate_scratch_reaction_norm(net, scratch_reaction_norm);
    node.reverse_mutate_bias(mutation);
}

////Calculates the reaction norm of a network with a mutated weights
/// using a scratch reaction norm
template<class Net>
void calculate_rn_for_weights_mut(Net& net,
                                  weight& weight,
                                  reac_norm& scratch_reaction_norm,
                                  const double& mutation)
{
    weight.mutate_weight(mutation);
    recalculate_scratch_reaction_norm(net, scratch_reaction_norm);
    weight.reverse_mutate_weight(mutation);
}

///Calculates the robustness of the fitness
/// of a network to a series of mutations
/// on all its loci
template<typename Func, class Net>
double calc_fitness_mutational_sensibility(Net& net,
                                           const std::vector<double> mutations,
                                           Func optimal_function,
                                           const range& input_range = {-1,1},
                                           int n_points = 100)
{

    std::vector<double> distances_differences;

    auto optimal_reac_norm = calculate_reaction_norm_from_function(optimal_function, input_range, n_points);
    auto base_reac_norm = calculate_reaction_norm(net, input_range, n_points);

    auto distance_base_rn_from_optimal_rn = rn_distance(optimal_reac_norm,base_reac_norm);

    for(const auto& mutation : mutations)
        for(auto& layer : net.get_net_weights())
            for(auto& node : layer)
            {
                auto distance_from_optimal_after_bias_mut = calculate_rn_distance_for_bias_mut(net, node, optimal_reac_norm, mutation);
                distances_differences.emplace_back(distance_base_rn_from_optimal_rn - distance_from_optimal_after_bias_mut);

                for(auto& current_weight : node.get_vec_mutable_weights())
                {
                    auto distance_from_optimal_after_weight_mut =  calc_rn_distance_for_weight_mut(net,
                                                                                                   current_weight,
                                                                                                   optimal_reac_norm,
                                                                                                   mutation);

                    distances_differences.emplace_back(distance_base_rn_from_optimal_rn - distance_from_optimal_after_weight_mut);
                }
            }

    return distances_differences.empty() ? 0 : calc_mean(distances_differences);
}

///Calculates the robustness of the phenotype produced by a network
/// to a series of mutations on all its loci
template<class Net>
double calc_phenotype_mutational_sensibility(Net& net,
                                             const std::vector<double> mutations,
                                             const range& input_range = {-1,1},
                                             int n_points = 100)
{

    std::vector<double> distances;
    auto base_reac_norm = calculate_reaction_norm(net, input_range, n_points);

    for(const auto& mutation : mutations)
        for(auto& layer : net.get_net_weights())
            for(auto& node : layer)
            {
                distances.emplace_back(calculate_rn_distance_for_bias_mut(net, node, base_reac_norm, mutation));

                for(auto& current_weight : node.get_vec_mutable_weights())
                {
                    distances.emplace_back(calc_rn_distance_for_weight_mut(net,
                                                                           current_weight,
                                                                           base_reac_norm,
                                                                           mutation));
                }
            }
    return distances.empty() ? 0 : calc_mean(distances);
}

///Calculates the robustness of the phenotype produced by a network
/// and its fitness
/// to a series of mutations on all its loci
template<class Net, typename Func>
fit_and_phen_sens_t calc_phen_and_fit_mut_sensibility(Net& net,
                                                      const std::vector<double>& mutations,
                                                      Func optimal_function,
                                                      const range& input_range = {-1,1},
                                                      int n_points = 100)
{

    std::vector<double> fitness_distances;
    std::vector<double> phenotype_distances;

    reac_norm optimal_reac_norm = calculate_reaction_norm_from_function(optimal_function, input_range, n_points);
    reac_norm base_reac_norm = calculate_reaction_norm(net, input_range, n_points);
    reac_norm scratch_reac_norm = base_reac_norm;

    double distance_base_rn_from_optimal_rn = rn_distance(optimal_reac_norm,base_reac_norm);

    for(const auto& mutation : mutations)
        for(auto& layer : net.get_net_weights())
        {
            auto active_nodes = layer | std::views::filter([](node& n){return is_active(n);});
            for(auto& node : active_nodes)
            {
                calculate_rn_for_bias_mut(net, node, scratch_reac_norm, mutation);
                fitness_distances.emplace_back(distance_base_rn_from_optimal_rn - rn_distance(optimal_reac_norm, scratch_reac_norm));
                phenotype_distances.emplace_back(rn_distance(base_reac_norm, scratch_reac_norm));

                auto active_weights = node.get_vec_mutable_weights() | std::views::filter([](weight& w){return is_active(w);});
                for(auto& current_weight : active_weights)
                {
                    calculate_rn_for_weights_mut(net, current_weight, scratch_reac_norm, mutation);
                    fitness_distances.emplace_back(distance_base_rn_from_optimal_rn - rn_distance(optimal_reac_norm, scratch_reac_norm));
                    phenotype_distances.emplace_back(rn_distance(base_reac_norm, scratch_reac_norm));
                }
            }
        }

    return fit_and_phen_sens_t{calc_mean(fitness_distances), calc_mean(phenotype_distances)};
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

///Counts the number of biases
/// in a network
template<class Net>
int count_biases(const Net& n)
{
    int count = 0;
    for(const auto& layer : n.get_net_weights())
        for(const auto& node : layer)
        {
            count++;
        }
    return count;
}

///Counts the number of nodes
/// in a network
template<class Net>
int count_nodes(const Net& n)
{
    int count = 0;
    for(const auto& layer : n.get_net_weights())
        for(const auto& node : layer)
        {
            count += node.get_vec_weights().size();
        }
    return count;
}

///Counts the number of weights and biases (mutable loci by standard mode of mutation)
/// in a network
template<class Net>
int count_weights_and_biases(const Net& n)
{
    int count = 0;
    for(const auto& layer : n.get_net_weights())
        for(const auto& node : layer)
        {
            count++;
            count += node.get_vec_weights().size();
        }
    return count;
}


///Calculates the reaction_norm that is produced by an optimal function
/// for a given range and a given number of data points
template<class Func>
reac_norm calculate_reaction_norm_from_function(const Func& func,
                                                const range& cue_range,
                                                const int& n_data_points)
{
    reac_norm r_norm;
    r_norm.reserve(n_data_points);
    std::vector<double> input;
    double step_size = (cue_range.m_end - cue_range.m_start) / n_data_points;
    if(step_size == 0)
    {
        input.assign({ cue_range.m_start });
        r_norm.emplace_back(cue_range.m_start, func(input));
        return r_norm;
    }

    for (double i = cue_range.m_start; i < cue_range.m_end; i += step_size)
    {
        input.assign({ i });
        r_norm.emplace_back(i, func(input));
    }
    return r_norm;
}

///Checks if the distance of the reaction norm of a network
/// from an optimal function is 0
template<typename Fun, class Net>
bool rn_is_equal_to_optimal_rn(const Net& net,
                               Fun function,
                               range input_range,
                               int n_points)
{
    auto net_rn = calculate_reaction_norm(net, input_range, n_points);
    auto optimal_rn = calculate_reaction_norm_from_function(function, input_range, n_points);

    return net_rn == optimal_rn;
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

    for(size_t i = 0; i != size_t(n_repeats); ++i)
    {
        std::vector<double> input;
        for(size_t j = 0; j != input_size; ++j){
            input.push_back(int(i + j));
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
bool all_weigths_have_value(const Net &n, double value)
{
    auto weights = n.get_net_weights();

    for(auto &layer : weights ){
        for(auto &node : layer){
            for (size_t i = 0; i != node.get_vec_weights().size(); ++i){
                const class weight& current_weight = node.get_vec_weights()[i];
                if(current_weight.get_weight() != value)
                {
                    return false;
                }
            }
        }
    }
    return true;
}

///Checks that all weights have an absolute value higher than a certain amount
template<class Net>
bool all_weigths_have_higher_abs_value(const Net &n, double value);
//{
//    auto weights = n.get_net_weights();

//    for(const auto &layer : weights )
//    {
//        for(const auto &node : layer)
//        {
//            if(is_active(node))
//            {
//                for(const auto& current_weight : node.get_vec_weights())
//                {
//                    if(is_active(current_weight) &&
//                            std::abs(current_weight.get_weight()) < value)
//                    {
//                        return false;
//                    }
//                }
//            }
//        }
//    }
//    return true;
//}

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

    if(n.any_layer_has_no_nodes())
    {
        throw std::runtime_error{"One layer has 0 nodes"};
    }
    std::vector<double> w{0};
    std::vector<double> output;

    for(size_t layer = 0; layer != n.get_net_weights().size(); layer++)
    {
        output.resize(n.get_net_weights().at(layer).size());

        for(size_t node = 0; node != n.get_net_weights().at(layer).size(); node++)
        {
            if(n.get_net_weights().at(layer).size() != output.size())
            {
                throw std::runtime_error{"A node got lost somewhere in output()"};
            }

            const auto &current_node = n.get_net_weights().at(layer).at(node);

            double node_value = 0;

            if(current_node.is_active())
            {
                //std::vector<weight> vec_w = current_node.get_vec_weights();
                //convert_to_double_or_zero(vec_w).swap(w);

                //if(input.size() != w.size())
                //{
                //    throw std::runtime_error{"incoming weights and incoming inputs are not the same number"};
                //}

                //node_value = current_node.get_bias() +
                //        std::inner_product(input.begin(),
                //                           input.end(),
                //                           w.begin(),
                //                           0.0);

                const std::vector<weight>& vec_w = current_node.get_vec_weights();
                node_value = current_node.get_bias() +
                        std::inner_product(input.begin(),
                                           input.end(),
                                           vec_w.begin(),
                                           0.0,
                                           [](double a, double b) { return a + b; },
                [](double a, const auto& b) { return a * (b.get_weight() * b.is_active()); }
                );
            }

            output.at(node) = n(node_value);
        }

        output.swap(input);
    }

    return input;
}


// populates 'output' with the output of a network for a given input
// 'input' is used as scratch-memory
template<class Net>
void output_ugly_but_fast(const Net& n, std::vector<double>& input, std::vector<double>& output)
{
    if constexpr(Net::response_t == response_type::additive)
    {
        output.resize(1);
        const react_norm_t& gene =  *std::find_if(n.get_genes().begin(), n.get_genes().end(),
                                                  [&](const auto& gene){return gene.m_x == input[0];});
        output[0] = gene.m_y;
        return;
    }
    assert(input.size() == n.get_input_size());

    if (n.any_layer_has_no_nodes())
    {
        throw std::runtime_error{ "One layer has 0 nodes" };
    }

    auto n_layers = n.get_net_weights().size();
    for (size_t layer = 0; layer != n_layers; layer++)
    {
        auto n_nodes = n.get_net_weights()[layer].size();
        output.resize(n_nodes);
        for (size_t node = 0; node != n_nodes; node++)
        {
            const auto& current_node = n.get_net_weights()[layer][node];
            if (current_node.is_active())
            {
                const std::vector<weight>& vec_w = current_node.get_vec_weights();
                output[node] = n(current_node.get_bias() +
                                 std::inner_product(input.begin(),
                                                    input.end(),
                                                    vec_w.begin(),
                                                    0.0,
                                                    [](double a, double b) { return a + b; },
                [](double a, const auto& b) { return a * (b.get_weight() * b.is_active()); }
                ));
            }
            else
            {
                output[node] = 0;
            }
        }
        output.swap(input);
    }
    output.swap(input); //undo swap in last iteration
}

// populates 'output' with the output of a network for a given input
// 'input' is used as scratch-memory
template<class Net>
void output_with_modified_weight(const Net& n,
                                 std::vector<double>& input,
                                 std::vector<double>& output,
                                 size_t modified_layer_index,
                                 size_t modified_node_index,
                                 size_t modified_weight_index,
                                 double modified_weight)
{
    assert(input.size() == n.get_input_size());

    if (n.any_layer_has_no_nodes())
    {
        throw std::runtime_error{ "One layer has 0 nodes" };
    }

    bool is_modified_layer = false;
    bool is_modified_node = false;

    for (size_t layer = 0; layer != n.get_net_weights().size(); layer++)
    {
        is_modified_layer = layer == modified_layer_index;

        output.clear();
        for (size_t node = 0; node != n.get_net_weights().at(layer).size(); node++)
        {

            is_modified_node = node == modified_node_index;

            const auto& current_node = n.get_net_weights().at(layer).at(node);

            double node_value = 0;
            if (current_node.is_active())
            {
                if(is_modified_layer && is_modified_node)
                {
                    std::vector<weight> changed_weights = current_node.get_vec_weights();
                    changed_weights[modified_weight_index].change_weight(modified_weight);

                    node_value = current_node.get_bias() +
                            std::inner_product(input.begin(),
                                               input.end(),
                                               changed_weights.begin(),
                                               0.0,
                                               [](double a, double b) { return a + b; },
                    [](double a, const auto& b) { return a * (b.get_weight() * b.is_active()); }
                    );
                }
                else
                {
                    const std::vector<weight>& vec_w = current_node.get_vec_weights();
                    node_value = current_node.get_bias() +
                            std::inner_product(input.begin(),
                                               input.end(),
                                               vec_w.begin(),
                                               0.0,
                                               [](double a, double b) { return a + b; },
                    [](double a, const auto& b) { return a * (b.get_weight() * b.is_active()); }
                    );
                }
            }
            output.push_back(n(node_value));
        }
        output.swap(input);
    }
    output.swap(input); //undo swap in last iteration
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

            double node_value = 0;

            if(current_node.is_active()){
                node_value = current_node.get_bias() +
                        std::inner_product(input.begin(),
                                           input.end(),
                                           current_node.get_vec_weights().begin(),
                                           0.0,
                                           std::plus<>(),
                                           [](const double& lhs, const weight& rhs)
                {return lhs * rhs.is_active() * rhs.get_weight();});
            }
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

///Calculates the average number of incoming weights in a layer
///Only counting the active nodes and connections coming from active nodes
template<class Net>
inline double average_number_incoming_weights(const Net &n, size_t layer_index){
    std::vector<node> layer = n.get_net_weights()[layer_index];
    double total = 0;
    size_t layer_size_active = 0;

    for(const auto &node : layer){
        if(node.is_active()){
            ++layer_size_active;
            for(size_t i = 0; i != node.get_vec_weights().size(); ++i){
                if(node.get_vec_weights()[i].is_active() &&
                        (layer_index == 0 ? true : n.get_net_weights()[layer_index - 1][i].is_active())){
                    ++total;
                }
            }}
    }
    return total / layer_size_active;
}

///Calculates the average number of weights going out of a layer
template<class Net>
inline double average_number_outgoing_weights(const Net &n, size_t layer_index){
    size_t layer_size_active = 0;
    size_t next_layer_size_active = 0;
    std::vector<node> layer = n.get_net_weights()[layer_index];
    std::vector<node> next_layer = n.get_net_weights()[layer_index + 1];

    for(const auto &node : layer){
        if(node.is_active()){
            ++layer_size_active;
        }
    }
    for(const auto &node : next_layer){
        if(node.is_active()){
            ++next_layer_size_active;
        }
    }

    return (average_number_incoming_weights(n, layer_index + 1) * next_layer_size_active / layer_size_active);
}

///Returns the minimum bias of all the nodes in a layer
template<class Net>
inline double min_bias_in_layer(const Net &n, size_t layer){
    std::vector<double> biases;
    for(const node &node : n.get_net_weights()[layer]){
        biases.push_back(node.get_bias());
    }
    return *std::min_element(biases.begin(), biases.end());
}

///Returns the maximum bias of all the nodes in a layer
template<class Net>
inline double max_bias_in_layer(const Net &n, size_t layer){
    std::vector<double> biases;
    for(const node &node : n.get_net_weights()[layer]){
        biases.push_back(node.get_bias());
    }
    return *std::max_element(biases.begin(), biases.end());
}

///Returns the minimum weight of all the connections of all nodes in a layer
template<class Net>
inline double min_weight_in_layer(const Net &n, size_t layer){
    std::vector<double> weights;
    for(const node &node : n.get_net_weights()[layer]){
        for(const weight &weight : node.get_vec_weights()){
            weights.push_back(weight.get_weight());
        }
    }
    return *std::min_element(weights.begin(), weights.end());
}

///Returns the maximum weight of all the connections of all nodes in a layer
template<class Net>
inline double max_weight_in_layer(const Net &n, size_t layer){
    std::vector<double> weights;
    for(const node &node : n.get_net_weights()[layer]){
        for(const weight &weight : node.get_vec_weights()){
            weights.push_back(weight.get_weight());
        }
    }
    return *std::max_element(weights.begin(), weights.end());
}


///Calculates the reaction_norm of an individual's network
/// for a given range and a given number of data points
template<class Net>
reac_norm calculate_reaction_norm(const Net& net,
                                  const range& cue_range,
                                  const int& n_data_points)
{
    reac_norm r_norm;
    r_norm.reserve(n_data_points);
    std::vector<double> input;
    std::vector<double> output;
    double step_size = (cue_range.m_end - cue_range.m_start) / n_data_points;
    if(step_size == 0)
    {
        input.assign({ cue_range.m_start });
        output_ugly_but_fast(net, input, output);
        r_norm.emplace_back(cue_range.m_start, output[0]);
        return r_norm;
    }

    for (double i = cue_range.m_start; i < cue_range.m_end; i += step_size)
    {
        input.assign({ i });
        output_ugly_but_fast(net, input, output);
        r_norm.emplace_back(i, output[0]);
    }
    return r_norm;
}

///Recalculates the reaction_norm of an individual's network
///based on the inputs (.m_x) of a given reaction norm
template<class Net>
void recalculate_scratch_reaction_norm(const Net& net,
                                       reac_norm& r_norm)
{

    std::vector<double> input;
    std::vector<double> output;

    for (int i = 0; i < r_norm.size(); i++)
    {
        input.assign({ r_norm[i].m_x });
        output_ugly_but_fast(net, input, output);
        r_norm[i].m_y = output[0];
    }
}


///Calculates the reaction_norm of an individual's network
/// for a given range and a given number of data points
template<class Net>
reac_norm calculate_reaction_norm_with_modified_weight(const Net& net,
                                                       const range& cue_range,
                                                       const int& n_data_points,
                                                       size_t layer_index,
                                                       size_t node_index,
                                                       size_t weight_index,
                                                       double modified_weight)
{
    reac_norm r_norm;
    r_norm.reserve(n_data_points);
    std::vector<double> input;
    std::vector<double> output;
    double step_size = (cue_range.m_end - cue_range.m_start)/n_data_points;

    for(double i = cue_range.m_start; i < cue_range.m_end; i += step_size)
    {
        input.assign({ i });
        output_with_modified_weight(net,
                                    input,
                                    output,
                                    layer_index,
                                    node_index,
                                    weight_index,
                                    modified_weight);

        r_norm.emplace_back(i, output[0]);
    }

    return r_norm;
}

template<class N>
double weights_sum(const N& network)
{
    double sum = 0;
    for(const auto& layer : network.get_net_weights())
    {
        for(const auto& node : layer)
        {
            for(const auto& weight : node.get_vec_weights())
            {
                sum += weight.get_weight();
            }
        }
    }
    return sum;
}
void test_network();

#endif // NETWORK_H
