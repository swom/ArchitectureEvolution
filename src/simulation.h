#ifndef SIMULATION_H
#define SIMULATION_H

#include "selection_type.h"
#include "env_change_type.h"
#include "environment.h"
#include "population.h"
//#include <omp.h>

#include <fstream>
#include <vector>


double identity_first_element(const std::vector<double>& vector);

struct sim_param
{
    NLOHMANN_DEFINE_TYPE_INTRUSIVE(sim_param,
                                   seed,
                                   change_freq_A,
                                   change_freq_B,
                                   selection_strength,
                                   n_generations,
                                   selection_freq)

    sim_param(int seed_n = 0,
              double change_frequency_A = 0.1,
              double change_frequency_B = 0.01,
              double sel_strength = 1,
              int generations = 100,
              int selection_frequency = 1,
              env_change_symmetry_type env_change_symmetry_type = env_change_symmetry_type::symmetrical,
              env_change_freq_type env_change_freq_type = env_change_freq_type::stochastic,
              selection_type selec_type = selection_type::constant):
        seed{seed_n},
        change_freq_A{change_frequency_A},
        change_freq_B{change_frequency_B},
        selection_strength{sel_strength},
        n_generations{generations},
        selection_freq{selection_frequency},
        change_sym_type{env_change_symmetry_type},
        change_freq_type{env_change_freq_type},
        sel_type{selec_type}
    {}

    int seed;
    double change_freq_A;
    double change_freq_B;
    double selection_strength;
    int n_generations;
    int selection_freq;
    env_change_symmetry_type change_sym_type;
    env_change_freq_type change_freq_type;
    selection_type sel_type;

};

bool operator==(const sim_param& lhs, const sim_param& rhs);

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
    //assert(all_individuals_have_same_input(s));
    return get_inds(s)[0].get_input_values();
}

///Returns the size of the inputs of the individuals
template<class Sim>
size_t get_inds_input_size(const Sim &s)

{
    return get_inds_input(s).size();
}

///Updates the inputs in simulation and assigns them to individuals
template<class Sim>
void assign_new_inputs(Sim &s)
{
    auto new_inputs = s.create_inputs();

    if(s.get_input().size() > 1){
        new_inputs.back() = s.get_input().back();
    }

    s.update_inputs(new_inputs);
    assign_inputs(s);
}

template<class Pop = population<>,
         enum env_change_symmetry_type Env_change_sym = env_change_symmetry_type::symmetrical,
         enum env_change_freq_type Env_change_freq = env_change_freq_type::stochastic,
         enum selection_type Sel_Type = selection_type::constant>
class simulation
{
public:

    using pop_t = Pop;
    using env_ch_s_t = env_change_symmetry_type;
    using env_ch_f_t = env_change_freq_type;

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
        m_t_change_env_distr_A{static_cast<double>(params.s_p.change_freq_A)},
        m_t_change_env_distr_B{static_cast<double>(params.s_p.change_freq_B)},
        m_sel_str{params.s_p.selection_strength},
        m_change_freq_A {static_cast<double>(params.s_p.change_freq_A)},
        m_change_freq_B {static_cast<double>(params.s_p.change_freq_B)},
        m_selection_frequency{params.s_p.selection_freq},
        m_selection_duration{params.s_p.selection_freq / 10},
        m_params {params},
        m_input(params.i_p.net_par.net_arc[0], 1), //BAD!!! implementation of env function input
        m_optimal_output{1}
    {
        m_rng.seed(m_seed);
    }

    NLOHMANN_DEFINE_TYPE_INTRUSIVE(simulation,
                                   m_environment,
                                   m_population,
                                   m_time,
                                   m_change_freq_A,
                                   m_change_freq_B,
                                   m_sel_str,
                                   m_seed)

    ///Returns const ref ot population memeber
    const Pop& get_pop() const noexcept {return m_population;}

    ///Returns const ref ot population memeber
    Pop& get_pop() noexcept {return m_population;}

    ///Returns ref to rng
    std::mt19937_64& get_rng() noexcept {return m_rng;}

    ///Returns ref to environmental rng
    std::mt19937_64 &get_env_rng() noexcept {return m_environment.get_rng();}

    ///Returns const ref to env_member
    const environment& get_env() const noexcept {return m_environment;}

    ///Returns const ref to env_member
    environment& get_env() noexcept {return m_environment;}

    ///Returns the range of the inputs provided by environment
    const range& get_env_cue_range() const noexcept {return m_environment.get_cue_range();}

    ///Returns the number of generatiosn for which the simualtion has to run
    const double& get_mut_step() const noexcept {return m_params.p_p.mut_step;}

    ///Returns the number of generatiosn for which the simualtion has to run
    const int& get_n_gen() const noexcept {return m_n_generations;}

    ///Returns the number of generatiosn for which the simualtion has to run
    int get_n_trials() const noexcept {return m_population.get_n_trials();}

    ///returns const ref to Bernoulli distribution for change freq of A
    const std::bernoulli_distribution& get_t_change_env_distr_A() const noexcept {return m_t_change_env_distr_A;}

    ///returns const ref to Bernoulli distribution for change freq of B
    const std::bernoulli_distribution& get_t_change_env_distr_B() const noexcept {return m_t_change_env_distr_B;}

    ///returns the number of generations the simualtion has run for
    const int& get_time() const noexcept {return m_time;}

    ///increases the number of genration the simulations has run for
    void increase_time() {++m_time;}

    ///Returns the strength of selection
    double get_sel_str() const noexcept {return m_sel_str;}

    ///Returns the number of generations after which
    ///selection takes place
    int get_sel_freq() const noexcept {return m_selection_frequency;}

    ///Returns the number of generations for which
    ///selection takes place when selection is 'sporadic'
    int get_sel_duration() const noexcept {return m_selection_duration;}

    ///Returns change frequency of environment/function A
    double get_change_freq_A() const noexcept {return m_change_freq_A;}

    ///Returns change frequency of environment/function B
    double get_change_freq_B() const noexcept {return m_change_freq_B;}

    ///Returns seed
    int get_seed() const noexcept {return m_seed;}

    ///Returns a reference to the vector of individuals
    const std::vector<typename Pop::ind_t> &get_inds() const {return m_population.get_inds();};

    ///Returns the current inputs in the simulation
    const std::vector<double> &get_input() const noexcept {return m_input;}

    ///Returns the current optimal output
    const double &get_optimal() const noexcept {return m_optimal_output;}

    ///Checks if environment needs to change
    bool is_environment_changing(){
        if constexpr( Env_change_freq == env_change_freq_type::regular)
        {
            return std::fmod(get_time(), 1.0/m_change_freq_A)  == 0;
        }
        else if( Env_change_freq == env_change_freq_type::stochastic)
        {
            if( m_environment.get_name_current_function() == 'A' )
            {
                std::bernoulli_distribution distro = get_t_change_env_distr_A();
                return distro (get_env_rng());
            }
            else if (m_environment.get_name_current_function() == 'B')
            {
                std::bernoulli_distribution distro;
                if constexpr( Env_change_sym == env_change_symmetry_type::asymmetrical)
                {
                    distro = get_t_change_env_distr_B();
                }
                else if(Env_change_sym == env_change_symmetry_type::symmetrical)
                {
                    distro = get_t_change_env_distr_A();
                }
                return distro (get_env_rng());
            }
            else
                throw std::runtime_error{"invalid current function name"};
        }
    }

    ///Returns the function A of the environment
    const std::function<double(std::vector<double>)> &get_env_function_A() const noexcept
    {return get_env().get_env_function_A();}

    ///Updates the optimal to the given value
    void update_optimal(double new_optimal) {m_optimal_output = new_optimal;}

    ///Updates the inputs of the simulation with new calculated inputs
    void update_inputs(std::vector<double> new_inputs){m_input = new_inputs;}

    ///Evaluates the operformance of all indiivduals in a population
    std::vector<double> evaluate_inds(){

        std::vector<double> cumulative_performance(get_inds().size(), 0);

        std::vector<std::vector<double>> inputs(m_population.get_n_trials());
        std::vector<double> optimals(inputs.size());

        for(int i = 0; i != m_population.get_n_trials(); i++)
        {
            inputs[i] = create_inputs();
            optimals[i] = env::calculate_optimal(m_environment, inputs[i]);
        }

        ///BAD!!! Temporary solutions to pass test
        update_inputs(inputs[0]);
        update_optimal(optimals[0]);

#pragma omp parallel for
        for(int i = 0; i < m_population.get_n_trials(); i++)
        {
//            std::cout << "the number of threads used is "<< omp_get_num_threads() << std::endl;
            auto performance = pop::calc_dist_from_target(get_inds(),
                                                          optimals[i],
                                                          inputs[i]);
#pragma omp critical
            {
                std::transform(cumulative_performance.begin(),
                               cumulative_performance.end(),
                               performance.begin(),
                               cumulative_performance.begin(),
                               std::plus<double>());
            }
        }
        return cumulative_performance;
    }
    ///Changes the inputs in the environment of the simulation
    std::vector<double> create_inputs()
    {
        return(env::create_n_inputs(get_env(),
                                    get_inds_input_size(*this),
                                    get_rng()
                                    )
               );
    }

    ///Calculates fitness of inds in pop given current env values
    const simulation<Pop, Env_change_sym, Env_change_freq, Sel_Type>&
    calc_fitness()
    {
        auto cumulative_performance = evaluate_inds();

        auto fitness_vector = pop::rescale_dist_to_fit(cumulative_performance, get_sel_str());

        pop::set_fitness_inds(get_pop(), fitness_vector);

        return *this;
    }
    ///Reproduces inds to next gen based on their fitness
    void reproduce()
    {
        pop::reproduce(get_pop(), get_rng());
    }

    ///Reproduces inds to next gen randomly
    void reproduce_randomly()
    {
        pop::reproduce_random(get_pop(), get_rng());
    }
    ///Calculates fitness and selects a new population based on fitness
    void select_inds()
    {
        if constexpr(Sel_Type == selection_type::sporadic)
        {
            if(m_time % m_selection_frequency >= 0 &&
                    m_time % m_selection_frequency < m_selection_duration)
            {
                calc_fitness();
                reproduce();
            }
            else
            {
                calc_fitness();
                reproduce_randomly();
            }
        }
        else if constexpr(Sel_Type == selection_type::constant)
        {
            calc_fitness();
            reproduce();
        }
        else
        {
            throw std::runtime_error{"wrong type of selection"};
        }
    }

    ///Changes the last input (env function indicator) from 1 to -1 or vice versa
    ///the way this is implemented is BAD!!!
    void switch_env_indicator()
    {
        if(get_input().size() > 1){
            m_input.back() = -m_input.back();
        }
    }

    ///Resets the fitness of the population to 0
    void reset_fit_pop()
    {
        m_population.reset_fitness();
    }

    const all_params& get_params() const noexcept {return m_params;}

private:

    environment m_environment;
    Pop m_population;
    int m_n_generations;
    std::mt19937_64 m_rng;
    int m_seed;
    std::bernoulli_distribution m_t_change_env_distr_A;
    std::bernoulli_distribution m_t_change_env_distr_B;
    int m_time = 0;
    double m_sel_str;
    double m_change_freq_A;
    double m_change_freq_B;

    ///Every how many generations indiivdual are selected
    int m_selection_frequency;
    //For how many generations individuals are selected
    //A tenth of the selection frequency
    int m_selection_duration;

    all_params m_params;

    ///The current inputs that the networks of individuals will recieve
    std::vector<double> m_input;

    ///The optimal output at a given moment;
    /// depends on inputs and environmental function
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
    f >> json_in;
    Class s;
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


///Calculates the optimal output
template<class Sim>
double calculate_optimal(const Sim &s)
{
    return(env::calculate_optimal(s.get_env(), s.get_input()));
}

///Returns a population whose fitness has been calculated
template<class Sim>
typename Sim::pop_t calc_fitness_of_pop(Sim s)
{

    s.update_optimal(env::calculate_optimal(s.get_env(), s.get_input()));
    return pop::calc_fitness(s.get_pop(),
                             s.get_optimal(),
                             s.get_sel_str(),
                             s.get_input());
}

///Calculates the avg_fitness of the population
template<class Sim>
double avg_fitness(const Sim& s)
{
    return pop::avg_fitness(s.get_pop());
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

///Gets const ref the best n individuals in a pop
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
///checks if the individuals in the populations from 2 different simulations
///have exactly the same fitness values
template<class Sim>
bool pops_have_same_fitness(const Sim& lhs, const Sim& rhs)
{
    return pop::extract_fitnesses(lhs.get_inds()) == pop::extract_fitnesses(rhs.get_inds());
}

///sums the fitness of all individuals of a simulation toghether
template<class Sim>
double sum_of_fitnesses(const Sim& s)
{
    auto fitnesses = pop::extract_fitnesses(s.get_inds());
    return std::accumulate(fitnesses.begin(),
                           fitnesses.end(),
                           0.0);
}

///Ticks time one generation into the future
template<class Sim>
void tick(Sim &s)
{
    s.increase_time();

    if(s.is_environment_changing()){
        perform_environment_change(s);
    }

    s.select_inds();
}

///Calculates the standard devaition of the population fitness
template<class Sim>
double var_fitness(const Sim&s)
{
    return pop::stdev_fitness(s.get_pop());
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
