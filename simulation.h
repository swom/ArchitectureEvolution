#ifndef SIMULATION_H
#define SIMULATION_H

#include "environment.h"
#include "population.h"

#include <fstream>
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


template<class Pop = population<>>
class simulation
{
public:

    using pop_t = Pop;

    simulation(int init_pop_size = 1,
               int seed = 0,
               double t_change_interval = 0.1,
               std::vector<int> net_arch = {1,2,1},
               double sel_str = 2,
               int number_of_generations = 1000);

    simulation(const all_params& params):
        m_environment{params.e_p},
        m_population{params.p_p, params.i_p},
        m_n_generations{params.s_p.n_generations},
        m_seed{params.s_p.seed},
        m_t_change_env_distr{static_cast<double>(params.s_p.change_freq)},
        m_sel_str{params.s_p.selection_strength},
        m_change_freq {static_cast<double>(params.s_p.change_freq)},
        m_params {params},
        m_input(params.i_p.net_par.net_arc[0], 1),
        m_optimal_output{1}
    {
        m_rng.seed(m_seed);
    }

    NLOHMANN_DEFINE_TYPE_INTRUSIVE(simulation,
                                  m_population,
                                   m_time,
                                   m_change_freq,
                                   m_sel_str,
                                   m_seed)

    ///Returns const ref ot population memeber
    const Pop& get_pop() const noexcept {return m_population;}

    ///Returns const ref ot population memeber
    Pop& get_pop() noexcept {return m_population;}

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
    const std::vector<typename Pop::ind_t> &get_inds() const;

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

    ///Changes the last input (env function indicator) from 1 to -1 or vice versa
    void switch_env_indicator()
    {
        if(get_input().size() > 1){
            m_input.back() = -m_input.back();
        }
    }

    const all_params& get_params() const noexcept {return m_params;}

private:

    environment m_environment;
    Pop m_population;
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

///Loads a sim object from json
template<class Class>
Class load_json(
        const std::string& filename
        )
{
    std::ifstream f(filename);
    nlohmann::json json_in;
    Class s;
    f >> json_in;
    return s = json_in;
}

///Saves the enitre GODDDAM SIMULATIONNNN!!!!!!! WHOO NEEDS MEMORRYYYY
template<class Class>
void save_json(const Class& s, const std::string& filename)
{
    std::ofstream  f(filename);
    nlohmann::json json_out;
    json_out = s;
    f << json_out;
}

namespace sim {


///Checks if 2 simulations are equal
template<class Pop>
bool operator ==(const simulation<Pop>& lhs, const simulation<Pop>& rhs);

///Checks if all the individuals in a simulated population have the same input
template<class Sim>
bool all_individuals_have_same_input(const Sim &s)
{
    auto p = s.get_pop();

    return pop::all_individuals_have_same_input(p);
}

///Assigns the given new input to each individual in the simulation
template<class Sim>
void assign_new_inputs_to_inds(Sim &s, std::vector<double> new_input)
{
    pop::assign_new_inputs_to_inds(s.get_pop(), new_input);
}

///Assigns the input in simulation<M> to individuals
template<class Sim>
void assign_inputs(Sim &s)
{
    pop::assign_new_inputs_to_inds(s.get_pop(), s.get_input());
}

///Returns the individuals in the simualtion
template<class Sim>
const std::vector<typename Sim::pop_t::ind_t> &get_inds(const Sim&s)
{
    return s.get_pop().get_inds();
}

///Returns the input of the individuals
template<class Sim>
std::vector<double> get_inds_input(const Sim &s)
{
    assert(all_individuals_have_same_input(s));
    return get_inds(s)[0].get_input_values();
}

///Returns the size of the inputs of the individuals
template<class Sim>
size_t get_inds_input_size(const Sim &s)

{
    return get_inds_input(s).size();
}

///Changes the inputs in the environment of the simulation
template<class Sim>
std::vector<double> create_inputs(Sim s)
{
    auto &e = s.get_env();
    return(env::create_n_inputs(e, get_inds_input_size(s), s.get_rng() ));
}

///Updates the inputs in simulation and assigns them to individuals
template<class Sim>
void assign_new_inputs(Sim &s)
{
    std::vector<double> new_inputs = create_inputs(s);

    if(s.get_input().size() > 1){
        new_inputs.back() = s.get_input().back();
    }

    s.update_inputs(new_inputs);
    assign_inputs(s);
}

///Calculates the optimal output
template<class Sim>
double calculate_optimal(const Sim &s)
{
    return(env::calculate_optimal(s.get_env(), s.get_input()));
}

///Calculates the avg_fitness of the population
template<class Sim>
double avg_fitness(const Sim& s)
{
    return pop::avg_fitness(s.get_pop());
}
///Calculates fitness of inds in pop given current env values
template<class Sim>
void calc_fitness(Sim &s)
{
    s.update_optimal(calculate_optimal(s));
    s.get_pop() = pop::calc_fitness(s.get_pop(),
                                    s.get_optimal(),
                                    s.get_sel_str());
}

///Changes all the weights of a given individual to a given value
template<class Sim>
void change_all_weights_nth_ind(Sim& s, size_t ind_index, double new_weight);

///Changes the network of the nth individual for a given network
template<class Pop>
void change_nth_ind_net(simulation<Pop>& s, size_t ind_index, const typename Pop::ind_t::net_t &n)
{
    pop::change_nth_ind_net(s.get_pop(), ind_index, n) ;
}

///Gets the best n individuals in a pop
template<class Sim>
std::vector<typename Sim::pop_t::ind_t> get_best_n_inds(const Sim& s, int n)
{
    return pop::get_best_n_inds(s.get_pop(), n);
}

///Returns the current optimal function of the environment
template<class Sim>
std::function<double(std::vector<double>)> get_current_env_function(const Sim &s);

///Gets the name of the current environmental function
template<class Sim>
char get_name_current_function(const Sim& s) noexcept
{
    return s.get_env().get_name_current_function();
}

///Returns the fitness of the nth ind in pop
template<class Sim>
double get_nth_ind_fitness(const Sim& s, const size_t ind_index);

///Returns const or non-onst ref to the network of the nth individual in the
/// popoulation member of a simulation
template<class Sim>
const typename Sim::pop_t::ind_t::net_t& get_nth_ind_net(const Sim& s, size_t ind_index);

///Reproduces inds to next gen based on their fitness
template<class Sim>
void reproduce(Sim& s)
{
    pop::reproduce(s.get_pop(), s.get_rng());
}

///Calculates fitness and selects a new population based on fitness
template<class Sim>
void select_inds(Sim& s)
{
    calc_fitness(s);
    reproduce(s);
}

///Checks if environment should change
template<class Sim>
bool is_environment_changing(Sim &s) {
    std::bernoulli_distribution distro = s.get_t_change_env_distr();
    return distro (s.get_rng());
}

///Switches the function of the environment used to calculate the optimal output
template<class Sim>
void switch_optimal_function(Sim &s)
{
    env::switch_env_function(s.get_env());
}

///Wrapper function; does everything that needs doing when the environment changes
template<class Sim>
void perform_environment_change(Sim &s)
{
    switch_optimal_function(s);
    s.switch_env_indicator();
}

///Ticks time one generation into the future
template<class Sim>
void tick(Sim &s)
{
    s.increase_time();

    if(is_environment_changing(s)){

        perform_environment_change(s);
    }

    if(get_inds(s).size()){

        assign_new_inputs(s);

    }

    select_inds(s);
}

///Calculates the standard devaition of the population fitness
template<class Sim>
double var_fitness(const Sim&s)
{
    return pop::var_fitness(s.get_pop());
}

///Get the inputs of the individuals in the simulation. Requires all individuals to have the same input.
template<class Sim>
const std::vector<double> &get_current_input(const Sim &s);

///Returns the input of the nth individual in the population
template<class Sim>
const std::vector<double> &get_nth_individual_input(const Sim &s, const int n);

///Updates the input with the current environmental indicator
template<class Sim>
void update_env_indicator(Sim &s);

}

void test_simulation() noexcept;

#endif // SIMULATION_H
