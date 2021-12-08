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
              enum mutation_type mut = mutation_type::activation):
        net_par{net_pars},
        m_mutation_type{mut}
    {}
};


template<mutation_type M = mutation_type::weights>
class individual
{
public:

    individual(const ind_param &i_p = ind_param{}) :
        ///!!!!Attention!!!! input values are for now a fixed amount
        m_input_values(i_p.net_par.net_arc[0], 1.0)
      {
           m_network = std::make_unique<mutator_network<M>>(i_p.net_par);
      }

    individual(individual<M>&&) = default;
    individual(const individual<M> &i) noexcept:
        m_fitness{i.get_fitness()},
        m_input_values{i.get_input_values()},
        m_network{std::make_unique<network<M>>(i.get_net())}
    {}

    ///Overload of copy assignment where pointer to netowrk is not copied
    /// but the pointee is copied and assigned
    //copy assignment
    individual& operator=(const individual& other)
    {
        // Guard self assignment
        if (this == &other)
        {
            return *this;
        }

        m_fitness = other.get_fitness();
        m_input_values = other.get_input_values();
        if(m_network == nullptr || get_net() != other.get_net())
        {
            *m_network = other.get_net();
        }

        return *this;
    };

    // move assignment
    individual& operator=(individual&& other) noexcept
    {
        // Guard self assignment
        if (this == &other)
            return *this;

        m_input_values = std::exchange(other.get_to_input_values(),
                                       std::vector<double>{});
        m_fitness = std::exchange(other.get_to_fitness(),0);

        m_network = std::move(other.get_net_ptr());
        other.get_net_ptr().reset();

        return *this;
    }

    ///Changes the network of an individual with another network
    void change_net(const network<M>& n)
    {
        *m_network = n;
    }

    ///Returns copy of fitness
    const double& get_fitness() const noexcept {return m_fitness;}

    ///Return const referernce to vector of fixed input values
    const std::vector<double>& get_input_values() const noexcept {return m_input_values;}

    ///Returns const ref to network
    const network<M>& get_net() const noexcept {return *m_network;}

    ///Returns ref to the pointer to network
    std::unique_ptr<network<M>>& get_net_ptr() {return m_network;}


    ///Returns ref to fitness USED FOR JSON SAVING
    double& get_to_fitness() noexcept {return m_fitness;}

    ///Returns ref to inputs
    std::vector<double>& get_to_input_values() noexcept {return m_input_values;}

    ///Returns ref to network USED FOR JSON SAVING
    network<M>& get_to_net() noexcept {return *m_network;}

    ///Mutates the network of an individual
    void mutate(double mut_rate, double mut_step, std::mt19937_64 &rng)
    {
        m_network->mutate(mut_rate, mut_step, rng);
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
    std::unique_ptr<network<M>> m_network;
};


///Functions required to save to json format
using json = nlohmann::json;

template<mutation_type M>
void to_json(json& j, const individual<M>& ind)
{
    j = json{
    {"fitness", ind.get_fitness()},
    {"input_values", ind.get_input_values()},
    {"network", ind.get_net()}
};
}

template<mutation_type M>
void from_json(const json& j, individual<M>& ind) {
    j.at("fitness").get_to(ind.get_to_fitness());
    j.at("input_values").get_to(ind.get_to_input_values());
    j.at("network").get_to(ind.get_to_net());
}



/// Checks if 2 individuals are the same
template<mutation_type M>
bool operator== (const individual<M>& lhs, const individual<M>& rhs)
{
  bool fitness = are_equal_with_tolerance(lhs.get_fitness(), rhs.get_fitness());
  bool network = lhs.get_net() == rhs.get_net();
  bool inputs = lhs.get_input_values() == rhs.get_input_values();

  return fitness && network && inputs;
}

///Lets a network send out an ouput signal
///!!!!Attention!!! for now no input is provided
template<mutation_type M>
std::vector<double> response(const individual<M>& ind)
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
#endif // INDIVIDUAL_H
