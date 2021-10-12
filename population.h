#ifndef POPULATION_H
#define POPULATION_H

#include "individual.h"
#include "rndutils.hpp"
#include <vector>


struct pop_param
{
    NLOHMANN_DEFINE_TYPE_INTRUSIVE(pop_param,
                                   number_of_inds,
                                   mut_rate,
                                   mut_step)
int number_of_inds;
double mut_rate;
double mut_step;
};

class population
{
public:
  population(int init_nr_indiv = 1,
             double mut_rate = 0.01,
             double mut_step = 0.1
             , std::vector<int> net_arch = {1,2,1});
   population(pop_param p_p, ind_param i_p);

  NLOHMANN_DEFINE_TYPE_INTRUSIVE(population,
                               m_vec_indiv,
                               m_mut_rate,
                               m_mut_step);

  ///Get const ref to vector of individuals
  const std::vector<individual>& get_inds() const noexcept{return m_vec_indiv;}

  ///Get ref to vector of individuals
  std::vector<individual>& get_inds() noexcept{return m_vec_indiv;}

  ///Returns the ref tot the mutable fitness distribution
  rndutils::mutable_discrete_distribution<>& get_fitness_dist() noexcept{return m_fitness_dist;}
  ///Get const ref to vector of individuals
  const std::vector<individual>& get_new_inds() const noexcept{return m_vec_new_indiv;}

  ///Get ref to vector of individuals
  std::vector<individual>& get_new_inds() noexcept{return m_vec_new_indiv;}

  ///Return mutation rate
  double get_mut_rate() const noexcept {return m_mut_rate;}

  ///Return mutation step
  double get_mut_step() const noexcept {return m_mut_step;}

private:

  std::vector<individual> m_vec_indiv;
  std::vector<individual> m_vec_new_indiv;
  double m_mut_rate;
  double m_mut_step;
  rndutils::mutable_discrete_distribution<> m_fitness_dist;

};
///Checks that 2 populations are equal
bool operator== (const population& lhs, const population& rhs);

///Calculates the avg_fitness of the population
double avg_fitness(const population& p);

///Calculates the fitness of inds in pop given a target env_value
population calc_fitness(population p, const double &env_value, const double &sel_str);

///changes the net of the nth individual to a given net
void change_nth_ind_net(population& p, size_t ind_index, network n);

///Creates a mutable distribution from whihc to draw inds based on fitness
rndutils::mutable_discrete_distribution<>  create_mut_dist_fit(population& p);

///Extracts a vector of the fitnesses of individuals into a double vectors
std::vector<double> extract_fitnesses(const std::vector<individual>& inds);

///Gets the best n individuals in a pop
std::vector<individual> get_best_n_inds(const population& p, int nth);

const individual& get_nth_ind(const population& p, size_t ind_index);
individual& get_nth_ind(population& p, size_t ind_index);

///Returns the fitness of the nth individual
double get_nth_ind_fitness(const population& p, const size_t& ind_index);

const network& get_nth_ind_net(const population& p, size_t ind_index);
network& get_nth_ind_net( population& p, size_t ind_index);

///Rescales the distance fro the target of an ind
///to a fitness value between 0  and 1
std::vector<double> rescale_dist_to_fit(std::vector<double> distance_from_target,
                                        double selection_strength);

///Reproduces inds with a probability proportional to their fitness
void reproduce(population& p, std::mt19937_64& rng);

///Select inds for new pop from old pop based on mutable dist
/// and mutates them
void select_new_pop(population& p,
                    const rndutils::mutable_discrete_distribution<>& mut_dist,
                    std::mt19937_64 &rng);

///Sets the fitness of the individuals to the one contained in the fitness vector
void set_fitness_inds(population& p, const std::vector<double>& fitness_vector);

///Swaps a vector of new_inds with the vector of old inds
void swap_new_with_old_pop(population& p);

///Sets the fitness of the nth ind in the population
void set_nth_ind_fitness (population& p, size_t ind_index, double fitness);

///Calculates the standard deviation
double var_fitness(const population &p);

void test_population() noexcept;

#endif // POPULATION_H
