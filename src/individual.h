#ifndef INDIVIDUAL_H
#define INDIVIDUAL_H

#include "network.h"




static std::map<std::string, mutation_type> string_to_mut_type_map
{
    {"weights", mutation_type::weights},
    {"activation", mutation_type::activation},
    {"weights_and_activation", mutation_type::weights_and_activation},
    {"duplication", mutation_type::duplication},
    {"NRduplication", mutation_type::NRduplication},
    {"addition", mutation_type::addition},
    {"NRaddition", mutation_type::NRaddition}
};

struct ind_param
{
    ind_param(net_param net_pars = net_param(),
              enum mutation_type mut = mutation_type::weights):
        net_par{net_pars},
        m_mutation_type{mut}
    {}

    NLOHMANN_DEFINE_TYPE_INTRUSIVE(ind_param,
                                   net_par,
                                   m_mutation_type)
    net_param net_par;
    enum mutation_type m_mutation_type;


};
bool operator==(const ind_param& lhs, const ind_param& rhs);

template<class Net = network<>>
class individual
{
public:
    using net_t = Net;

    individual(const ind_param &i_p = ind_param{}) :
        ///!!!!Attention!!!! input values are for now a fixed amount
        m_input_values(i_p.net_par.net_arc[0], 1.0),
        m_network{i_p.net_par}
    {}


    NLOHMANN_DEFINE_TYPE_INTRUSIVE(individual,
                                   m_fitness,
                                   m_input_values,
                                   m_network);

    ///Calculates the mutational spectrum of the individual's network
    network_spectrum<Net> calc_spectrum(double mut_step,
                                        std::mt19937_64& rng,
                                        int n_mutations,
                                        const range& cue_range,
                                        const int& n_data_points
                                        )
    {

        return network_spectrum<Net> (m_network,
                                      mut_step,
                                      rng,
                                      n_mutations,
                                      cue_range,
                                      n_data_points
                                      );
    }
    ///Changes the network of an individual with another network
    void change_net(const Net& n)
    {
        m_network = n;
    }

    ///Changes all weights of a network to a given value
    void change_all_weights(double new_weight)
    {
        m_network.change_all_weights_values(new_weight);
    }

    ///Returns copy of fitness
    const double& get_fitness() const noexcept {return m_fitness;}

    ///Return const referernce to vector of fixed input values
    const std::vector<double>& get_input_values() const noexcept {return m_input_values;}

    ///Returns const ref to network
    const Net& get_net() const noexcept {return m_network;}

    ///Returns const ref to network
    Net& get_mutable_net() noexcept {return m_network;}

    ///Returns ref to fitness USED FOR JSON SAVING
    double& get_to_fitness() noexcept {return m_fitness;}

    ///Returns ref to inputs
    std::vector<double>& get_to_input_values() noexcept {return m_input_values;}

    ///Returns ref to network USED FOR JSON SAVING
    Net& get_to_net() noexcept {return m_network;}

    ///Returns the rank of the individual
    const int& get_rank() const noexcept
    {
        return m_rank;
    }

    ///Returns the rank of the ancestor of the individual
    const std::string& get_ID() const noexcept
    {
        return m_ID;
    }

    ///Returns the rank of the ancestor of the individual
    const std::string& get_ancestor_ID() const noexcept
    {
        return m_ancestor_ID;
    }

    ///Sets the rank of the individual
    void set_ID(const std::string& ID) noexcept {m_ID =  ID;}

    ///Sets the rank of the individual
    void set_rank(int rank) noexcept {m_rank =  rank;}

    ///Makes the actual ID of the individual the ancestor ID
    ///to be used when saving to set a new lineage
    void make_ID_ancestor_ID() {m_ancestor_ID = m_ID;};

    ///Mutates the network of an individual
    void mutate(double mut_rate_w, double mut_step, std::mt19937_64 &rng, double mut_rate_a, double mut_rate_d)
    {
        m_network.mutate(mut_rate_w, mut_step, rng, mut_rate_a, mut_rate_d);
    }

    ///Sets the fitness of an ind
    void set_fitness(double fitness) {m_fitness = fitness;}

    ///Resets fitness to  0
    void reset_fitness() {m_fitness = 0;}

    ///Set the input values of an individual
    void assign_input(const std::vector<double> &input) {m_input_values = input;}

private:

    ///The fitness of an individual
    double m_fitness = 0;

    ///The vector of fixed input values that will be given to the network
    std::vector<double> m_input_values;

    ///The network of an individual
    Net m_network;

    ///The rank in terms of fitness of the individual in the population
    int m_rank = 0;

    ///The rank in terms of fitness of the individual in the population
    std::string m_ID = "0";

    ///The rank of the ancestor
     std::string m_ancestor_ID = "0";
};

/// Checks if 2 individuals are the same
template<class Net>
bool operator== (const individual<Net>& lhs, const individual<Net>& rhs)
{
    bool fitness = are_equal_with_tolerance(lhs.get_fitness(), rhs.get_fitness());
    bool network = lhs.get_net() == rhs.get_net();
    bool inputs = lhs.get_input_values() == rhs.get_input_values();

    return fitness && network && inputs;
}

namespace ind {

///Lets a network send out an ouput signal
///!!!!Attention!!! for now no input is provided
template<class Ind>
std::vector<double> response(const Ind& ind,
                             const std::vector<double>& input)
{
    return output(ind.get_net(),input);
}

///Lets a network send out an ouput signal
/// using scratch memory vectors
template<class Ind>
void response_scratch(const Ind& ind,
                      std::vector<double>& input,
                      std::vector<double>& ouput
                      )
{
    return output_ugly_but_fast(ind.get_net(),input, ouput);
}

///Calculates the distance of a response of a network
/// and a given value
template<class Ind>
double calc_sqr_distance(const Ind &i,
                         double env_value,
                         const std::vector<double>& input)
{
    auto output = response(i, input);
    return (output[0] - env_value) * (output[0] - env_value);
}

///Calculates the distance of a response of a network
/// and a given value using scratch memory vectors
template<class Ind>
double calc_sqr_distance_scratch(const Ind &i,
                                 double env_value,
                                 std::vector<double>& input,
                                 std::vector<double>& output
                                 )
{
    response_scratch(i, input, output);
    return std::fabs(output[0] - env_value);
}

}


void test_individual();

#endif // INDIVIDUAL_H
