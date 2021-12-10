#ifndef INDIVIDUAL_H
#define INDIVIDUAL_H

#include "network.h"

static std::map<std::string, mutation_type> string_to_mut_type_map
{
    {"weights", mutation_type::weights},
    {"activation", mutation_type::activation},
    {"weights_and_activation", mutation_type::weights_and_activation}
};

struct ind_param
{
    NLOHMANN_DEFINE_TYPE_INTRUSIVE(ind_param,
                                   net_par,
                                   m_mutation_type)
    net_param net_par;
    enum mutation_type m_mutation_type;

    ind_param(net_param net_pars = net_param(),
              enum mutation_type mut = mutation_type::weights):
        net_par{net_pars},
        m_mutation_type{mut}
    {}
};


template<class Net = network<>>
class individual
{
public:
    using network_type = Net;

    individual(const ind_param &i_p = ind_param{}) :
        ///!!!!Attention!!!! input values are for now a fixed amount
        m_input_values(i_p.net_par.net_arc[0], 1.0),
        m_network{i_p.net_par}
      {}

    ///Changes the network of an individual with another network
    void change_net(const Net& n)
    {
        m_network = n;
    }

    ///Returns copy of fitness
    const double& get_fitness() const noexcept {return m_fitness;}

    ///Return const referernce to vector of fixed input values
    const std::vector<double>& get_input_values() const noexcept {return m_input_values;}

    ///Returns const ref to network
    const Net& get_net() const noexcept {return m_network;}

    ///Returns ref to the pointer to network
    Net& get_net_ptr() {return m_network;}

    ///Returns ref to fitness USED FOR JSON SAVING
    double& get_to_fitness() noexcept {return m_fitness;}

    ///Returns ref to inputs
    std::vector<double>& get_to_input_values() noexcept {return m_input_values;}

    ///Returns ref to network USED FOR JSON SAVING
    Net& get_to_net() noexcept {return m_network;}

    ///Mutates the network of an individual
    void mutate(double mut_rate, double mut_step, std::mt19937_64 &rng)
    {
        m_network.mutate(mut_rate, mut_step, rng);
    }

    ///Sets the fitness of an ind
    void set_fitness(double fitness) {m_fitness = fitness;}

    ///Set the input values of an individual
    void assign_input(const std::vector<double> &input) {m_input_values = input;}

private:

    ///The fitness of an individual
    double m_fitness = 0;

    ///The vector of fixed input values that will be given to the network
    std::vector<double> m_input_values;

    ///The network of an individual
    Net m_network;
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


///Functions required to save to json format
using json = nlohmann::json;

template<class Ind>
void to_json(json& j, const Ind& ind)
{
    j = json{
    {"fitness", ind.get_fitness()},
    {"input_values", ind.get_input_values()},
    {"network", ind.get_net()}
};
}

template<class Ind>
void from_json(const json& j, Ind& ind) {
    j.at("fitness").get_to(ind.get_to_fitness());
    j.at("input_values").get_to(ind.get_to_input_values());
    j.at("network").get_to(ind.get_to_net());
}



///Lets a network send out an ouput signal
///!!!!Attention!!! for now no input is provided
template<class Ind>
std::vector<double> response(const Ind& ind)
{
    return output(ind.get_net(),ind.get_input_values());
}

///Calculates the distance of a response of a network
/// and a given value
template<class Ind>
double calc_sqr_distance(const Ind &i, double env_value)
{
    auto output = response(i);
    return (output[0] - env_value) * (output[0] - env_value);
}

void test_individual();
}
#endif // INDIVIDUAL_H
