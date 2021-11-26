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
                                   mutation_type)
    net_param net_par;
    enum mutation_type mutation_type;

    ind_param(net_param net_pars = net_param(),
              enum mutation_type mut = mutation_type::activation):
        net_par{net_pars},
        mutation_type{mut}
    {}
};



class individual
{
public:

    individual(ind_param i_p = {});
    individual(individual&&) = default;
    individual(const individual& i) noexcept;

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
               m_network = std::make_unique<network>(other.get_net());
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
    void change_net(const network& n);

    ///Returns copy of fitness
    const double& get_fitness() const noexcept {return m_fitness;}

    ///Return const referernce to vector of fixed input values
    const std::vector<double>& get_input_values() const noexcept {return m_input_values;}

    ///Returns const ref to network
    const network& get_net() const noexcept {return *m_network;}

    std::unique_ptr<network>& get_net_ptr() {return m_network;}

    ///Returns ref to fitness USED FOR JSON SAVING
    double& get_to_fitness() noexcept {return m_fitness;}

    ///Returns ref to inputs
    std::vector<double>& get_to_input_values() noexcept {return m_input_values;}

    ///Returns ref to network USED FOR JSON SAVING
    network& get_to_net() noexcept {return *m_network;}

    ///Mutates the network of an individual
    void mutate(double mut_rate, double mut_step, std::mt19937_64 &rng);

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
    std::unique_ptr<network> m_network;
};


///Functions required to save to json format
using json = nlohmann::json;

void to_json(json& j, const individual& ind);

void from_json(const json& j, individual& ind);


/// Checks if 2 individuals are the same
bool operator== (const individual& lhs, const individual& rhs);

///Calculates the distance of a response of a network
/// and a given value
double calc_sqr_distance(const individual& i, double env_value);

///Lets a network send out an ouput signal
///!!!!Attention!!! for now no input is provided
std::vector<double> response(const individual& ind);

void test_individual();
#endif // INDIVIDUAL_H
