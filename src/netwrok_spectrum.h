#ifndef NETWROK_SPECTRUM_H
#define NETWROK_SPECTRUM_H
#include <nlohmann/json.hpp>
#include<random>
#include <vector>
#include "range.h"

///Stores the phenpotype and fitness sensibilities to mutations of a network
struct fit_and_phen_sens_t
{

    fit_and_phen_sens_t(double fitness_sens = 123456,
                        double phenotype_sens = 123456,
                        int rank = -13456,
                        double fitness = -11111):
        m_fitness_sens{fitness_sens},
        m_phenotype_sens{phenotype_sens},
        m_rank{rank},
        m_fitness{fitness}
    {if(m_phenotype_sens < 0)
        {
            throw std::runtime_error{"During construction of fit_and_phen_t:"
                                  "sensibility constructed with phenotipyc sensibility < 0,"
                                  "phenotypic sensibility cannot be < 0"};}
    }

    NLOHMANN_DEFINE_TYPE_INTRUSIVE(fit_and_phen_sens_t,
                                   m_fitness_sens,
                                   m_phenotype_sens,
                                   m_fitness,
                                   m_rank,
                                   m_ID,
                                   m_ancestor_ID)

    ///Returns the rank of the idnidvidual (whihc works also as ID)
    int get_rank() const noexcept {return m_rank;}

    ///Returns the rank of the idnidvidual (whihc works also as ID)
    int get_ID() const noexcept {return m_ID;}

    ///Returns the rank of the idnidvidual (whihc works also as ID)
    int get_ancestor_ID() const noexcept {return m_ancestor_ID;}

    ///The sensitibility of fitness to mutations,
    /// negative values signify that mutations on average move
    /// the network away from the fitness optima
    /// positive values signify that mutations move it closer
    double m_fitness_sens;

    ///The sensitibility of the phenotype to mutations,
    /// higher valuessignify that mutation on average
    /// move the phenotype further away from its original state
    double m_phenotype_sens;

    ///The rank of the individual from which the sensibilities are calculated
    int m_rank;

    ///The ID of the individual from which the sensibilities are calculated
    int m_ID;

    ///The ID of the ancestir if individual from which the sensibilities are calculated
    int m_ancestor_ID;

    ///The fitness of the individual
   double m_fitness;
};

bool operator== (const fit_and_phen_sens_t& lhs, const fit_and_phen_sens_t& rhs);
bool operator!= (const fit_and_phen_sens_t& lhs, const fit_and_phen_sens_t& rhs);

struct react_norm_t { 
    react_norm_t(){};
    react_norm_t(double x, double y):
        m_x{x},
        m_y{y}
    {};
    NLOHMANN_DEFINE_TYPE_INTRUSIVE(react_norm_t, m_x, m_y);

    double m_x;
    double m_y;
};
bool operator==(const react_norm_t& lhs, const react_norm_t& rhs);

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

    NLOHMANN_DEFINE_TYPE_INTRUSIVE(network_spectrum,
                                   m_current_reac_norm,
                                   m_net_spectrum_weights_for_weights_mutation
                                   //                                       ,
                                   //                                       m_net_spectrum_biases_for_weights_mutation,
                                   //                                       m_net_spectrum_for_act_mutation,
                                   //                                       m_net_spectrum_for_dup,
                                   //                                       m_net_spectrum_for_nradd,
                                   //                                       m_net_spectrum_for_del
                                   )

    ///
    ////// \brief get_reac_norm
    /// Returns the reaction norm stored
    ////// \return
    const reac_norm &get_reaction_norm() const noexcept{return m_current_reac_norm;}
    ///
    ////// \brief get_weights_mut_spec
    /// Returns the mutational spectrum when weights are mutated
    ////// \return
    const net_weight_mut_spectrum &get_weights_mut_spec() const noexcept{return m_net_spectrum_weights_for_weights_mutation;}


private :

    ///The current reaction norm without any mutations
    reac_norm m_current_reac_norm;

    ///Reaction norm spectrum for each WEIGHT when WEIGHT MUTATION
    net_weight_mut_spectrum m_net_spectrum_weights_for_weights_mutation;

    ///Reaction norm spectrum for each BIAS when WEIGHT MUTATION
    //net_b_weight_mut_spectrum m_net_spectrum_biases_for_weights_mutation;

    ///Reaction norm spectrum for each NODE when ACTIVATION MUTATION
    //net_act_mut_spectrum m_net_spectrum_for_act_mutation;

    /// Reaction norm spectrum for each LAYER when DUPLICATION MUTATION
    //net_node_level_mutation_mut_spectrum m_net_spectrum_for_dup;

    /// Reaction norm spectrum for each LAYER when NRADDITION MUTATION
    //net_node_level_mutation_mut_spectrum m_net_spectrum_for_nradd;

    /// Reaction norm spectrum for each LAYER when DELETION MUTATION
    //net_node_level_mutation_mut_spectrum m_net_spectrum_for_del;
};

template<class Net>
bool operator==(const network_spectrum<Net>& lhs, const network_spectrum<Net>& rhs)
{
    bool react_norm = lhs.get_reaction_norm() == rhs.get_reaction_norm();
    bool mut_spet_weights = lhs.get_weights_mut_spec() == rhs.get_weights_mut_spec();
    return react_norm && mut_spet_weights;
}

#endif // NETWROK_SPECTRUM_H
