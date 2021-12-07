#ifndef SIMULATION_H
#define SIMULATION_H

#include "environment.h"
#include "population.h"
#include <vector>

double identity_first_element(const std::vector<double>& vector);

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
                                   i_p,
                                   p_p,
                                   s_p)
 env_param e_p;
 ind_param i_p;
 pop_param p_p;
 sim_param s_p;


};

template<mutation_type M= mutation_type::weights>
class simulation
{
public:

  simulation(int init_pop_size = 1,
             int seed = 0,
             double t_change_interval = 0.1,
             std::vector<int> net_arch = {1,2,1},
             double sel_str = 2,
             int number_of_generations = 1000
          );
  simulation (const all_params& params);

  NLOHMANN_DEFINE_TYPE_INTRUSIVE(simulation,

                                 m_population,
                                 m_time,
                                 m_change_freq,
                                 m_sel_str,
                                 m_seed)

  ///Returns const ref ot population memeber
  const population<M>& get_pop() const noexcept {return m_population;}

  ///Returns const ref ot population memeber
  population<M>& get_pop() noexcept {return m_population;}

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
  const std::vector<individual<M>> &get_inds() const;

  ///Returns the current inputs in the simulation
  const std::vector<double> &get_input() const noexcept {return m_input;}

  ///Returns the current optimal output
  const double &get_optimal() const noexcept {return m_optimal_output;}

  ///Returns the function A of the environment
  const std::function<double(std::vector<double>)> &get_env_function_A() const noexcept
    {return get_env().get_env_function_A();}

  ///Updates the optimal to the given value
  void update_optimal(double new_optimal) {m_optimal_output = new_optimal;}

  ///Updates the inputs of the simulation with new calculated inputs
  void update_inputs(std::vector<double> new_inputs){m_input = new_inputs;}



  const all_params& get_params() const noexcept {return m_params;}

  private:

   environment m_environment;
   population<M> m_population;
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
template<mutation_type M>
bool operator ==(const simulation<M>& lhs, const simulation<M>& rhs);

///Checks if all the individuals in a simulated population have the same input
template<mutation_type M>
bool all_individuals_have_same_input(const simulation<M> &s);

///Assigns the given new input to each individual in the simulation
template<mutation_type M>
void assign_new_inputs_to_inds(simulation<M> &s, std::vector<double> new_input);

///Assigns the input in simulation<M> to individuals
template<mutation_type M>
void assign_inputs(simulation<M> &s);

///Assign inputs to a population
template<mutation_type M>
void assign_new_inputs_to_inds(population<M> &p, const std::vector<double> &inputs);

///Updates the inputs in simulation and assigns them to individuals
template<mutation_type M>
void assign_new_inputs(simulation<M> &s);

///Calculates the avg_fitness of the population
template<mutation_type M>
double avg_fitness(const simulation<M>& s)
{
    return avg_fitness(s.get_pop());
}

///Calculates fitness of inds in pop given current env values
template<mutation_type M>
void calc_fitness(simulation<M>& s);

///Changes all the weights of a given individual to a given value
template<mutation_type M>
void change_all_weights_nth_ind(simulation<M>& s, size_t ind_index, double new_weight);

///Changes the network of the nth individual for a given network
template<mutation_type M>
void change_nth_ind_net(simulation<M>& s, size_t ind_index, const network<M>& n);

///Gets the best n individuals in a pop
template<mutation_type M>
std::vector<individual<M>> get_best_n_inds(const simulation<M>& s, int n)
{
    return get_best_n_inds(s.get_pop(), n);
}

///Returns the input of the individuals
template<mutation_type M>
std::vector<double> get_inds_input(const simulation<M> &s);

///Returns the size of the inputs of the individuals
template<mutation_type M>
size_t get_inds_input_size(const simulation<M> &s);

///Returns the current optimal function of the environment
template<mutation_type M>
std::function<double(std::vector<double>)> get_current_env_function(const simulation<M> &s);

///Gets the name of the current environmental function
template<class S>
char get_name_current_function(const S& s) noexcept
{
    return s.get_env().get_name_current_function();
}

///Returns the individuals in the simualtion
template<mutation_type M>
const std::vector<individual<M>>& get_inds(const simulation<M>&s);

///Returns the fitness of the nth ind in pop
template<mutation_type M>
double get_nth_ind_fitness(const simulation<M>& s, const size_t ind_index);

///Returns const or non-onst ref to the network of the nth individual in the
/// popoulation member of a simulation
template<mutation_type M>
const network<M>& get_nth_ind_net(const simulation<M>& s, size_t ind_index);

///Saves the enitre GODDDAM SIMULATIONNNN!!!!!!! WHOO NEEDS MEMORRYYYY
template<mutation_type M>
void save_json(const simulation<M>& s, const std::string& filename);

///Calculates fitness and selects a new population based on fitness
template<mutation_type M>
void select_inds(simulation<M>& s);

///Ticks time one generation into the future
template<mutation_type M>
void tick(simulation<M> &s);


///Calculates the standard devaition of the population fitness
template<mutation_type M>
double var_fitness(const simulation<M>&s)
{
    return var_fitness(s.get_pop());
}


///Get the inputs of the individuals in the simulation. Requires all individuals to have the same input.
template<mutation_type M>
const std::vector<double> &get_current_input(const simulation<M> &s);

///Returns the input of the nth individual in the population
template<mutation_type M>
const std::vector<double> &get_nth_individual_input(const simulation<M> &s, const int n);

///Changes the inputs in the environment of the simulation
template<mutation_type M>
std::vector<double> create_inputs(simulation<M> s);

///Calculates the optimal output
template<mutation_type M>
double calculate_optimal(const simulation<M> &s);


///Switches the function of the environment used to calculate the optimal output
template<mutation_type M>
void switch_optimal_function(simulation<M> &s);


void test_simulation() noexcept;

#endif // SIMULATION_H
