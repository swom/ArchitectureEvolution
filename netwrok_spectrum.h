#ifndef NETWROK_SPECTRUM_H
#define NETWROK_SPECTRUM_H
#include "network.h"
#include "range.h"

using reac_norm = std::unordered_map<double,double>;

//typedefs for spectrum of weigths
using weight_mut_spectrum = std::vector<reac_norm>;
using node_weight_mut_spectrum = std::vector<weight_mut_spectrum>;
using layer_weight_mut_spectrum = std::vector<node_weight_mut_spectrum>;
using net_weight_mut_spectrum = std::vector<layer_weight_mut_spectrum>;

//typedefs for spectrum of biases
using bias_weight_mut_spectrum = std::vector<reac_norm>;
using layer_b_weight_mut_spectrum = std::vector<bias_weight_mut_spectrum>;
using net_b_weight_mut_spectrum = std::vector<layer_b_weight_mut_spectrum>;

//typedefs for spectrum of nodes active connections
using node_act_mut_spectrum = std::vector<reac_norm>;
using layer_act_mut_spectrum = std::vector<node_act_mut_spectrum>;
using net_act_mut_spectrum = std::vector<layer_act_mut_spectrum>;

//typedefs for spectrum of layers' nodes
using layer_mut_spectrum = std::vector<reac_norm>;
using net_node_level_mutation_mut_spectrum = std::vector<layer_mut_spectrum>;

struct network_spectrum
{
public:
    network_spectrum();
    ///Reaction norm spectrum for each WEIGHT when WEIGHT MUTATION
    net_weight_mut_spectrum m_net_spectrum_weights_for_weights_mutation;
    ///Reaction norm spectrum for each BIAS when WEIGHT MUTATION
    net_b_weight_mut_spectrum m_net_spectrum_biases_for_weights_mutation;
    ///Reaction norm spectrum for each NODE when ACTIVATION MUTATION
    net_act_mut_spectrum m_net_spectrum_for_act_mutation;
    /// Reaction norm spectrum for each LAYER when DUPLICATION MUTATION
    net_node_level_mutation_mut_spectrum m_net_spectrum_for_dup;
    /// Reaction norm spectrum for each LAYER when NRADDITION MUTATION
    net_node_level_mutation_mut_spectrum m_net_spectrum_for_nradd;
    /// Reaction norm spectrum for each LAYER when DELETION MUTATION
    net_node_level_mutation_mut_spectrum m_net_spectrum_for_del;
};

///Calculates the reaction_norm of an individual's network
/// for a given range and a given number of data points
template<class Net>
reac_norm calculate_reaction_norm(const Net& net,
                                  const range& cue_range,
                                  const int& n_data_points)
{
    reac_norm r_norm;
    r_norm.reserve(n_data_points);
    double step_size = (cue_range.m_end - cue_range.m_start)/n_data_points;
    for(double i = cue_range.m_start; i < cue_range.m_end; i += step_size)
    {
        r_norm.insert({i, output(net, std::vector<double>{i}).at(0)});
    }
    return r_norm;
}

template<class Network>
net_weight_mut_spectrum calc_spectrum_weights_for_weights_mutation(const Network& n,
                                                                   double mut_step,
                                                                   std::mt19937_64& rng,
                                                                   int n_mutations,
                                                                   int n_inputs,
                                                                   range input_range)

{
    auto mutable_net = n;
    std::normal_distribution<double> mut_dist(0, mut_step);

    net_weight_mut_spectrum network_weights_spectrum(mutable_net.get_net_weights().size());
    layer_weight_mut_spectrum layer_spectrum;
    node_weight_mut_spectrum node_spectrum;
    weight_mut_spectrum weight_spectrum;
    reac_norm rn;
    rn.reserve(n_inputs);

    for(auto layer_it = mutable_net.get_net_weights().begin(); layer_it != mutable_net.get_net_weights().end(); layer_it++)
    {
        layer_spectrum.resize(layer_it->size());
        for(auto node_it = layer_it->begin(); node_it != layer_it->end(); node_it++)
        {
            node_spectrum.resize(node_it->get_vec_weights().size());
            for(auto weight = 0; weight != node_it->get_vec_weights().size(); weight++)
            {
                auto original_weight = node_it->get_vec_weights()[weight];
                for(int i = 0; i != n_mutations; i++)
                {
                    node_it->change_nth_weight(original_weight +  mut_dist(rng), weight);
                    rn = calculate_reaction_norm(n, input_range, n_inputs);
                    node_it->change_nth_weight(original_weight, weight);
                }

                auto index_weight = std::distance(node_it->begin(), weight);
                node_spectrum[index_weight] = rn;
                rn.clear();
                rn.reserve(n_inputs);
            }
            auto index_node = std::distance(layer_it->begin(), node_it);
            layer_spectrum[index_node] = node_spectrum;
        }
        auto index_layer = std::distance( mutable_net.get_net_weights().begin(), layer_it);
        network_weights_spectrum[index_layer] = layer_spectrum;
    }
    return network_weights_spectrum;

}
#endif // NETWROK_SPECTRUM_H
