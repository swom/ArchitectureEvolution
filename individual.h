#ifndef INDIVIDUAL_H
#define INDIVIDUAL_H

#include "network.h"

struct ind_param
{
    NLOHMANN_DEFINE_TYPE_INTRUSIVE(ind_param,
                                   net_par,
                                   age)
net_param net_par;
int age;
};



class individual
{
public:
  individual( std::vector<int> net_arch = std::vector<int>{1,2,1}, int age = 0);
  individual(ind_param i_p);

  NLOHMANN_DEFINE_TYPE_INTRUSIVE(individual,
                                 m_fitness,
                                 m_input_values,
                                 m_network);

  ///Returns the age of the individual
  int get_age() const noexcept {return m_age;}

  ///Return const referernce to vector of fixed input values
  const std::vector<double>& get_input_values() const noexcept {return m_input_values;}

  ///Returns the fitness
  double get_fitness() const noexcept {return m_fitness;}

  ///Returns const ref to network
  const network& get_net() const noexcept {return m_network;}

  ///Returns ref to network
  network& get_net() noexcept {return m_network;}

  ///Mutates the network of an individual
  void mutate(double mut_rate, double mut_step, std::mt19937_64 &rng);

  ///Sets the fitness of an ind
  void set_fitness(double fitness) {m_fitness = fitness;}

private:

  ///The age of the individual
  int m_age;

  ///The fitness of an individual
  double m_fitness = 0;

  ///The vector of fixed input values that will be given to the network
  std::vector<double> m_input_values;

  ///The network of an individual
  network m_network;
};

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
