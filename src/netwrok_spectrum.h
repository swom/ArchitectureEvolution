#ifndef NETWROK_SPECTRUM_H
#define NETWROK_SPECTRUM_H
#include <nlohmann/json.hpp>
#include<random>
#include <vector>
#include "range.h"

struct react_norm_t { 
  NLOHMANN_DEFINE_TYPE_INTRUSIVE(react_norm_t, x, y); 
  friend bool operator==(const react_norm_t& a, const react_norm_t& b) = default;

  double x, y; 
};

using reac_norm = std::vector<react_norm_t>;

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

template<class Network>
class network_spectrum
{
public:
    network_spectrum(){};
    network_spectrum(Network& net,
                     double mut_step,
                     std::mt19937_64 rng,
                     int n_mutations,
                     range input_range,
                     int input_number):
        m_current_reac_norm{calculate_reaction_norm(net, input_range, input_number)},
        m_net_spectrum_weights_for_weights_mutation{net.calc_spectrum_weights_for_weights_mutation( mut_step,
                                                                                                    rng,
                                                                                                    n_mutations,
                                                                                                    input_range,
                                                                                                    input_number
                                                                                                    )}
    {};

    NLOHMANN_DEFINE_TYPE_NON_INTRUSIVE(network_spectrum,
                                       m_net_spectrum_weights_for_weights_mutation,
                                       m_net_spectrum_biases_for_weights_mutation,
                                       m_net_spectrum_for_act_mutation,
                                       m_net_spectrum_for_dup,
                                       m_net_spectrum_for_nradd,
                                       m_net_spectrum_for_del
                                       )
    private :

        ///The current reaction norm without any mutations
        reac_norm m_current_reac_norm;
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

#endif // NETWROK_SPECTRUM_H
