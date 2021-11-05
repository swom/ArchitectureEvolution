#ifndef SIMULATION_H
#define SIMULATION_H

#include "environment.h"
#include "population.h"
#include <vector>

struct sim_param
{
    NLOHMANN_DEFINE_TYPE_INTRUSIVE(sim_param,
                                   seed,
                                   change_freq,
                                   selection_strength,
                                   n_generations)
 int seed;
 double change_freq;
 double selection_strength;
 int n_generations;

};

struct all_params
{
    NLOHMANN_DEFINE_TYPE_INTRUSIVE(all_params,
                                   e_p,
                                   i_p,
                                   p_p,
                                   s_p)
 env_param e_p;
 ind_param i_p;
 pop_param p_p;
 sim_param s_p;

};

class simulation
{
public:

  simulation(double targetA = 0,
             double targetB = 0.1,
             int init_pop_size = 1,
             int seed = 0,
             double t_change_interval = 0.1,
             std::vector<int> net_arch = {1,2,1},
             double sel_str = 2,
             int number_of_generations = 1000
          );
  simulation (all_params params);

  NLOHMANN_DEFINE_TYPE_INTRUSIVE(simulation,
                                 m_environment,
                                 m_population,
                                 m_time,
                                 m_change_freq,
                                 m_sel_str,
                                 m_seed)

  ///Returns const ref ot population memeber
  const population& get_pop() const noexcept {return m_population;}

  ///Returns const ref ot population memeber
  population& get_pop() noexcept {return m_population;}

  ///Returns ref to rng
  std::mt19937_64& get_rng() noexcept {return m_rng;}

  ///Returns const ref to env_member
  const environment& get_env() const noexcept {return m_environment;}

  ///Returns const ref to env_member
  environment& get_env() noexcept {return m_environment;}

  ///Returns the number of generatiosn for which the simualtion has to run
  const int& get_n_gen() const noexcept {return m_n_generations;}

  ///returns const ref to
  const std::bernoulli_distribution& get_t_change_env_distr() const noexcept {return m_t_change_env_distr;}
  std::bernoulli_distribution& get_t_change_env_distr() noexcept {return m_t_change_env_distr;}
  const int& get_time() const noexcept {return m_time;}
  void increase_time() {++m_time;}

  ///Returns the strength of selection
  double get_sel_str() const noexcept {return m_sel_str;}

  ///Returns change frequency
  double get_change_freq() const noexcept {return m_change_freq;}

  ///Returns seed
  int get_seed() const noexcept {return m_seed;}

  ///Returns a reference to the vector of individuals
  const std::vector<individual> &get_inds() const;

  ///Returns the current inputs in the simulation
  const std::vector<double> &get_input() const noexcept {return m_input;}

  ///Returns the input of the individuals
  const std::vector<double> &get_inds_input() const;

  ///Returns the size of the inputs of the individuals
  size_t get_inds_input_size() const{return get_inds_input().size();}

  ///Returns the current optimal output
  const double &get_optimal() const noexcept {return m_optimal_output;}

  ///Returns the function A of the environment
  const std::function<double(std::vector<double>)> &get_env_function_A() const noexcept
    {return get_env().get_env_function_A();}

  ///Updates the optimal to the given value
  void update_optimal(double new_optimal) {m_optimal_output = new_optimal;}

  ///Updates the inputs of the simulation with new calculated inputs
  void update_inputs();


  const all_params& get_params() const noexcept {return m_params;}

  private:

   environment m_environment;
   population m_population;
   int m_n_generations;
   std::mt19937_64 m_rng;
   int m_seed;
   std::bernoulli_distribution m_t_change_env_distr;
   int m_time = 0;
   double m_sel_str;
   double m_change_freq;
   all_params m_params;

   ///The current inputs that the networks of individuals will recieve
   std::vector<double> m_input;

   ///The optimal output at a given moment; depends on inputs and environmental function
   double m_optimal_output;

};
///Checks if 2 simulations are equal
bool operator ==(const simulation& lhs, const simulation& rhs);

///Calculates the avg_fitness of the population
double avg_fitness(const simulation& s);

///Calculates fitness of inds in pop given current env values
void calc_fitness(simulation& s);

///Changes all the weights of a given individual to a given value
void change_all_weights_nth_ind(simulation& s, size_t ind_index, double new_weight);

///Changes the network of the nth individual for a given network
void change_nth_ind_net(simulation& s, size_t ind_index, const network& n);

///Gets the best n individuals in a pop
std::vector<individual> get_best_n_inds(const simulation& s, int n);

///Returns the current target value of the environment in the simualtion
double get_current_env_value(const simulation&s);

///Returns the current target value of the environment in the simualtion
double get_current_env_value(simulation&s);

///Returns the individuals in the simualtion
const std::vector<individual>& get_inds(const simulation&s);

///Returns the fitness of the nth ind in pop
double get_nth_ind_fitness(const simulation& s, const size_t ind_index);

///Returns const or non-onst ref to the network of the nth individual in the
/// popoulation member of a simulation
const network& get_nth_ind_net(const simulation& s, size_t ind_index);
network& get_nth_ind_net( simulation& s, size_t ind_index);

///Calculates fitness and selects a new population based on fitness
void select_inds(simulation& s);

///Ticks time one generation into the future
void tick(simulation &s);

///Saves the enitre GODDDAM SIMULATIONNNN!!!!!!! WHOO NEEDS MEMORRYYYY
void save_json(const simulation& s, const std::string& filename);

///Calculates the standard devaition of the population fitness
double var_fitness(const simulation&s);

///Assign inputs to a population
void assign_new_inputs_to_inds(population &p, const std::vector<double> &inputs);

///Checks if all the individuals in a simulated population have the same input
bool all_individuals_have_same_input(const simulation &s);

///Get the inputs of the individuals in the simulation. Requires all individuals to have the same input.
const std::vector<double> &get_current_input(const simulation &s);

///Returns the input of the nth individual in the population
const std::vector<double> &get_nth_individual_input(const simulation &s, const int n);

///Changes the inputs in the environment of the simulation
std::vector<double> create_inputs(simulation &s);

///Calculates the optimal output
double calculate_optimal(const simulation &s);

///Assigns the given new input to each individual in the simulation
void assign_new_inputs_to_inds(simulation &s, std::vector<double> new_input);


void test_simulation() noexcept;

#endif // SIMULATION_H
