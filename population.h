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

template <mutation_type M = mutation_type::weights>
class population
{
public:
  population(int init_nr_indiv = 1,
             double mut_rate = 0.01,
             double mut_step = 0.1,
             std::vector<int> net_arch = {1,2,1});
   population(const pop_param &p_p, const ind_param &i_p);

  NLOHMANN_DEFINE_TYPE_INTRUSIVE(population,
                               m_vec_indiv,
                               m_mut_rate,
                               m_mut_step);

  ///Changes the network of the nth individual to a given network
  void change_nth_ind_net(size_t ind_index, const network& n){
      m_vec_indiv[ind_index].change_net(n);
  }


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
template< mutation_type M>
bool operator== (const population<M>& lhs, const population<M>& rhs);

///Calculates the avg_fitness of the population
template< mutation_type M>
double avg_fitness(const population<M>& p);

///Calculates the fitness of inds in pop given a target env_value
template< mutation_type M>
population<M>& calc_fitness(population<M> &p, const double &env_value, const double &sel_str);

///changes the net of the nth individual to a given net
template< mutation_type M>
void change_nth_ind_net(population<M>& p, size_t ind_index, network n);

///Creates a mutable distribution from whihc to draw inds based on fitness
template< mutation_type M>
rndutils::mutable_discrete_distribution<>  create_mut_dist_fit(population<M>& p);

///Extracts a vector of the fitnesses of individuals into a double vectors
std::vector<double> extract_fitnesses(const std::vector<individual>& inds);

///Gets the best n individuals in a pop
template< mutation_type M>
std::vector<individual> get_best_n_inds(const population<M>& p, int nth);

template< mutation_type M>
const individual& get_nth_ind(const population<M>& p, size_t ind_index);

template< mutation_type M>
individual& get_nth_ind(population<M>& p, size_t ind_index);

///Returns the fitness of the nth individual
template< mutation_type M>
double get_nth_ind_fitness(const population<M>& p, const size_t& ind_index);

template<mutation_type M>
const network& get_nth_ind_net(const population<M>& p, size_t ind_index);

///Rescales the distance fro the target of an ind
///to a fitness value between 0  and 1
std::vector<double> rescale_dist_to_fit(std::vector<double> distance_from_target,
                                        double selection_strength);

///Reproduces inds with a probability proportional to their fitness
template<mutation_type M>
void reproduce(population<M>& p, std::mt19937_64& rng);

///Select inds for new pop from old pop based on mutable dist
/// and mutates them
template<mutation_type M>
void select_new_pop(population<M>& p,
                    const rndutils::mutable_discrete_distribution<>& mut_dist,
                    std::mt19937_64 &rng);

///Sets the fitness of the individuals to the one contained in the fitness vector
template<mutation_type M>
void set_fitness_inds(population<M>& p, const std::vector<double>& fitness_vector);

///Swaps a vector of new_inds with the vector of old inds
template<mutation_type M>
void swap_new_with_old_pop(population<M> &p);

///Sets the fitness of the nth ind in the population
template<mutation_type M>
void set_nth_ind_fitness (population<M>& p, size_t ind_index, double fitness);

///Calculates the standard deviation
template<mutation_type M>
double var_fitness(const population<M> &p);

///Checks that all individuals in the pop have the same input
template<mutation_type M>
bool all_individuals_have_same_input(const population<M> &p);

///Returns the input of the nth individual
template<mutation_type M>
const std::vector<double> &get_nth_individual_input(const population<M> &p, const int n);

void test_population() noexcept;

#endif // POPULATION_H
