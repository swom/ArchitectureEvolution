#ifndef SIMULATION_H
#define SIMULATION_H
#include <filesystem>
#include <fstream>
#include <vector>

#include "selection_type.h"
#include "env_change_type.h"
#include "evaluation_type.h"
#include "adaptation_period.h"
#include "environment.h"
#include "population.h"
//#include <omp.h>

double identity_first_element(const std::vector<double>& vector);

struct sim_param
{
    NLOHMANN_DEFINE_TYPE_INTRUSIVE(sim_param,
                                   seed,
                                   change_freq_A,
                                   change_freq_B,
                                   selection_strength,
                                   selection_duration,
                                   n_generations,
                                   selection_freq,
                                   selection_duration_prop_to_freq,
                                   change_sym_type,
                                   change_freq_type,
                                   sel_type,
                                   adaptation_per,
                                   evalu_type,
                                   m_reac_norm_n_points)


    sim_param(int seed_n = 0,
              double change_frequency_A = 0.1,
              double change_frequency_B = 0.01,
              double sel_strength = 1,
              int generations = 100,
              int selection_frequency = 1,
              int selec_duration_prop_to_freq = 1,
              int reaction_norm_n_points = 40,
              int adaptation_period_proportion = 10,
              env_change_symmetry_type env_change_symmetry_type = env_change_symmetry_type::symmetrical,
              env_change_freq_type env_change_freq_type = env_change_freq_type::stochastic,
              selection_type selec_type = selection_type::constant,
              adaptation_period adapt_per = adaptation_period::off,
              evaluation_type eval_type = evaluation_type::full_rn):
        seed{seed_n},
        change_freq_A{change_frequency_A},
        change_freq_B{change_frequency_B},
        selection_strength{sel_strength},
        n_generations{generations},
        selection_freq{selection_frequency},
        selection_duration_prop_to_freq{selec_duration_prop_to_freq},
        selection_duration{selection_freq == 0 ? 0 : selection_freq / selec_duration_prop_to_freq},
                           m_reac_norm_n_points{reaction_norm_n_points},
                           adaptation_period_proportion{adaptation_period_proportion},
                           change_sym_type{env_change_symmetry_type},
                           change_freq_type{env_change_freq_type},
                           sel_type{selec_type},
                           adaptation_per{adapt_per},
                           evalu_type{eval_type}
    {
                           if(selection_duration &&
                              (static_cast<double>(selection_freq) / static_cast<double>(selec_duration_prop_to_freq)) < 1)
    {
                           throw std::invalid_argument{"Simulation parameters:"
                                                       "the numbers provided for the seleciton frequency "
                                                       "and the proportion of selection time between seleciton events are incorrect,"
                                                       " as they selection period would be shorter than 1 generation"};
}
}

                           int seed;
    double change_freq_A;
    double change_freq_B;
    double selection_strength;
    int n_generations;
    int selection_freq;
    int selection_duration_prop_to_freq;
    int selection_duration;

    ///The number of data points on which to calculate the reaction norm of an individual
    int m_reac_norm_n_points;

    ///The proportion of time n/sim_time in the simulation were selection will be constant (only if adaptive period is on)
    int adaptation_period_proportion;

    env_change_symmetry_type change_sym_type;
    env_change_freq_type change_freq_type;
    selection_type sel_type;
    adaptation_period adaptation_per;
    evaluation_type evalu_type;
};

bool operator==(const sim_param& lhs, const sim_param& rhs);

struct all_params
{
    all_params(const env_param& e_p = env_param{},
               const ind_param& i_p = ind_param{},
               const pop_param& p_p = pop_param{},
               const sim_param& s_p = sim_param{}):
        e_p{e_p},
        i_p{i_p},
        p_p{p_p},
        s_p{s_p}
    {}

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


///Assigns the given new input to each individual in the simulation
template<class Sim>
void assign_new_inputs_to_inds(Sim &s, std::vector<double> new_input)
{
    pop::assign_new_inputs_to_inds(s.get_pop_non_const(), new_input);
}

///Assigns the input in simulation<M> to individuals
template<class Sim>
void assign_inputs(Sim &s)
{
    pop::assign_new_inputs_to_inds(s.get_pop_non_const(), s.create_inputs());
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

///Returns the additive genes of the first individual
/// in the population
template<class Sim>
const reac_norm& get_first_ind_genes(const Sim &s)
{
    return get_inds(s)[0].get_net().get_genes();
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
         enum selection_type Sel_Type = selection_type::constant,
         enum adaptation_period Adapt_per = adaptation_period::off,
         enum evaluation_type Eval_type = evaluation_type::full_rn>
class simulation
{
public:

    using pop_t = Pop;
    using env_ch_s_t = env_change_symmetry_type;
    using env_ch_f_t = env_change_freq_type;
    using sel_t = selection_type;
    using adapt_p = adaptation_period;
    static constexpr response_type Resp_type = pop_t::ind_t::net_t::response_t;

    simulation(int init_pop_size = 1,
               int seed = 0,
               double t_change_interval = 0.1,
               std::vector<int> net_arch = {1,2,1},
               double sel_str = 2,
               int number_of_generations = 1000):
        m_environment{},
        m_population{init_pop_size},
        m_n_generations{number_of_generations},
        m_seed{seed},
        m_t_change_env_distr_A{static_cast<double>(t_change_interval)},
        m_t_change_env_distr_B{static_cast<double>(t_change_interval)},
        m_sel_str{sel_str},
        m_change_freq_A{static_cast<double>(t_change_interval)},
        m_change_freq_B{static_cast<double>(t_change_interval)},
        m_input(net_arch[0], 1),
        m_optimal_output{1}
    {
        m_rng.seed(m_seed);
        for(auto& ind : m_population.get_inds_nonconst())
        {
            ind = individual{net_param{net_arch, linear, net_arch}};
        }
    }

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
        m_selection_duration{params.s_p.selection_duration},
        m_params {params},
        m_input(params.i_p.net_par.net_arc[0], 1), //BAD!!! implementation of env function input
        m_optimal_output{1}
    {
        m_rng.seed(m_seed);
        if constexpr(Pop::ind_t::net_t::response_t == response_type::additive)
        {
            if(params.s_p.m_reac_norm_n_points != params.i_p.net_par.n_sampled_inputs)
            {
                throw std::invalid_argument{"simulation on construction: Additive response selected,"
                                            " but number of genes for inds is not the number of points in reaction norm"};
            }
            else if(get_first_ind_genes(*this) != calculate_reaction_norm_from_function(constant_zero,
                                                                                        params.i_p.net_par.input_range,
                                                                                        params.i_p.net_par.n_sampled_inputs))
            {
                throw std::invalid_argument{"simulation on construction: Additive response selected,"
                                            " but genes x values for inds will not correspond to inputs"};
            }
        }
    }

    NLOHMANN_DEFINE_TYPE_INTRUSIVE(simulation,
                                   m_environment,
                                   m_population,
                                   m_time,
                                   m_change_freq_A,
                                   m_change_freq_B,
                                   m_sel_str,
                                   m_seed)

    ///Assigns a unique ID to indivuduals
    Pop assign_ID_to_inds() noexcept {return m_population.assign_ID_to_inds(m_time);};

    ///Calculates the mutational sensibility to fitness and phenotype of all individuals in the population
    std::vector<fit_and_phen_sens_t> calculate_fit_phen_mut_sens_for_all_inds(int n_mutations, int n_points)
    {
        return m_population.calculate_fit_phen_mut_sens_for_all_inds(n_mutations,
                                                                     m_rng,
                                                                     m_environment.get_current_function(),
                                                                     m_params.e_p.cue_range,
                                                                     n_points);
    }

    ///Changes all the weights of a given individual to a given value
    void change_all_weights_nth_ind(size_t ind_index, double new_weight)
    {
        auto new_net = change_all_weights_values_and_activations(pop::get_nth_ind_net(m_population,
                                                                                      ind_index),
                                                                 new_weight);
        pop::change_nth_ind_net(m_population,
                                ind_index,
                                new_net);
    }

    ///Changes weights of all individuals
    void changel_all_inds_weights(double new_weight)
    {
        for(int ind = 0; ind != get_inds().size(); ind++)
            change_all_weights_nth_ind(ind, new_weight);
    }

    ///Returns const ref ot population memeber
    const Pop& get_pop() const noexcept {return m_population;}

    ///Returns const ref ot population memeber
    Pop& get_pop_non_const() noexcept {return m_population;}

    ///Gets the size of the population
    int get_pop_size() const noexcept {return m_population.get_inds().size();}

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

    ///Returns the type of selection used in simulation
    auto get_sel_type() const noexcept{return m_params.s_p.sel_type;}

    ///Returns the type of enviromental change frequency used in simulation
    auto get_e_change_f_type() const noexcept{return m_params.s_p.change_freq_type;}

    ///Returns the strength of selection
    double get_sel_str() const noexcept {return m_sel_str;}

    ///Returns the number of generations after which
    ///selection takes place
    int get_sel_freq() const noexcept {return m_selection_frequency;}

    ///Returns the number of generations for which
    ///selection takes place when selection is 'sporadic'
    int get_selection_duration() const noexcept {return m_selection_duration;}

    ///Returns change frequency of environment/function A
    double get_change_freq_A() const noexcept {return m_change_freq_A;}

    ///Returns change frequency of environment/function B
    double get_change_freq_B() const noexcept {return m_change_freq_B;}

    ///Returns seed
    int get_seed() const noexcept {return m_seed;}

    ///Returns a reference to the vector of individuals
    const std::vector<typename Pop::ind_t> &get_inds() const {return m_population.get_inds();};

    ///Returns a const reference to the vector of individuals
    const std::vector<typename Pop::ind_t> &get_new_inds() const noexcept {return m_population.get_new_inds();};

    ///Returns a reference to the vector of individuals
    std::vector<typename Pop::ind_t> &get_inds_non_const() {return m_population.get_inds_nonconst();};

    ///Returns a reference to the vector of individuals
    std::vector<typename Pop::ind_t> &get_new_inds_non_const() {return m_population.get_new_inds_nonconst();};

    ///Returns the current inputs in the simulation for the current or last trial
    const std::vector<double> &get_input() const noexcept {return m_input;}

    //Returns a const reference to the inputs given to individuals in this generation
    const std::vector<std::vector<double>>& get_stored_inputs() const noexcept {return m_stored_inputs;}

    ///Returns the current optimal output
    const double &get_optimal() const noexcept {return m_optimal_output;}

    ///Returns the selection duration to frequency of seleciton proportion
    /// if selection frequency / selection_duration_prop_to_freq = selection period duration;
    const int& get_selection_duration_prop_to_freq() const noexcept {return m_params.s_p.selection_duration_prop_to_freq;}

    //Returns a const reference to the inputs given to individuals in this generation
    const std::vector<double>& get_stored_optimals() const noexcept {return m_stored_optimal_output;}

    ///Checks if environment needs to change
    bool is_environment_changing(){
        if(m_time == 0) return false;

        if constexpr (Adapt_per == adaptation_period::on)
        {
            if(m_time <= (m_n_generations / m_params.s_p.adaptation_period_proportion))
                return false;
        }

        if constexpr( Env_change_freq == env_change_freq_type::regular)
        {
            if( m_environment.get_name_current_function() == 'A' )
            {
                return std::fmod(get_time(), 1.0/m_change_freq_A) == 0;
            }
            else if (m_environment.get_name_current_function() == 'B')
            {
                bool change;
                if constexpr( Env_change_sym == env_change_symmetry_type::asymmetrical)
                {
                    change = std::fmod(get_time(), 1.0/m_change_freq_B) == 0;
                }
                else if(Env_change_sym == env_change_symmetry_type::symmetrical)
                {
                    change = std::fmod(get_time(), 1.0/m_change_freq_A) == 0;
                }
                return change;
            }
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
        return false;
    }

    ///Returns the function A of the environment
    const std::function<double(std::vector<double>)> &get_env_function_A() const noexcept
    {return get_env().get_env_function_A();}

    ///Returns the number corresponding to the current environmental function
    ///0 for env_function 'A'
    ///1 for env_function 'B'
    int get_number_for_current_env_function() const noexcept {return m_environment.get_name_current_function() - 'A';}

    ///Updates the optimal to the given value
    void update_optimal(double new_optimal) {m_optimal_output = new_optimal;}

    ///Updates the inputs of the simulation with new calculated inputs
    void update_inputs(std::vector<double> new_inputs){m_input = new_inputs;}

    ///Updates the optimal to the given value
    void store_optimals(std::vector<double> new_optimal) {m_stored_optimal_output = new_optimal;}

    ///Updates the inputs of the simulation with new calculated inputs
    void store_inputs(std::vector<std::vector<double>> new_inputs){m_stored_inputs = new_inputs;}

    ///Creates the inputs for a simulation where networks
    /// use an additive gene response mechanism
    void create_inputs_additive_response( std::vector<std::vector<double>>& inputs)
    {
        auto genes = get_first_ind_genes(*this);
        std::shuffle(genes.begin(), genes.end(), m_rng);
        genes.resize(inputs.size());
        std::transform(genes.begin(), genes.end(),
                       inputs.begin(),
                       [](const auto& rn_t){return std::vector<double>{rn_t.m_x};});
    }

    ///Creates the inputs and optimals for a simulation where networks
    /// use a additive gene response mechanism
    void create_inputs_optimals_additive_response(std::vector<std::vector<double>>& inputs, std::vector<double>& optimals)
    {
        create_inputs_additive_response(inputs);
        for(int i = 0; i < inputs.size(); i++)
        {
            optimals[i] = env::calculate_optimal(m_environment, inputs[i]);
        }
    }

    ///Creates the inputs and optimals for a simulation where networks
    /// use a network response mechanism
    void create_inputs_optimals_network_response(std::vector<std::vector<double>>& inputs, std::vector<double>& optimals)
    {
        for(int i = 0; i < m_population.get_n_trials(); i++)
        {
            inputs[i] = create_inputs();
            optimals[i] = env::calculate_optimal(m_environment, inputs[i]);
        }
    }

    ///Gets inputs bsaed on the environment of the simulation
    /// and updates the input stored in simulation
    std::vector<double> create_inputs()
    {
        std::vector<double> inputs;

        inputs = env::create_n_inputs(get_env(),
                                      get_inds_input_size(*this),
                                      get_rng()
                                      );

        if constexpr(Resp_type == response_type::plastic)
        {
            inputs.push_back(get_number_for_current_env_function());
        }

        m_input = inputs;
        return inputs;
    }

    ///Creates and stores the ininputs and optimal values for those
    /// inputs to be then used inthe evaluation of individuals
    void create_and_store_inputs_and_optimals(std::vector<std::vector<double>>& inputs,
                                              std::vector<double>& optimals)
    {
        for(int i = 0; i != m_population.get_n_trials(); i++)
        {
            inputs[i] = create_inputs();
            optimals[i] = env::calculate_optimal(m_environment, inputs[i]);
        }

        store_inputs(inputs);
        store_optimals(optimals);
    }


    ///Creates the inputs and optimals for trial evaluations
    void create_input_optimals_trial(std::vector<std::vector<double>>& inputs,
                                     std::vector<double>& optimals)
    {
        auto trials = m_population.get_n_trials();
        inputs.resize(trials);
        optimals.resize(inputs.size());
        if constexpr(Resp_type == response_type::additive)
        {
            create_inputs_optimals_additive_response(inputs, optimals);
        }
        else
        {
            create_inputs_optimals_network_response(inputs, optimals);

        }
        store_inputs(inputs);
        store_optimals(optimals);
    }

    ///Creates the inputs and optimals for full_rn evalutaion
    void create_inputs_optimal_full_rn(std::vector<std::vector<double>>& inputs,
                                       std::vector<double>& optimals)
    {
        auto optimal_rn = calculate_reaction_norm_from_function(m_environment.get_current_function(),
                                                                m_environment.get_cue_range(),
                                                                m_params.s_p.m_reac_norm_n_points);

        inputs.resize(optimal_rn.size());
        optimals.resize(inputs.size());
        for(int i = 0; i != optimal_rn.size(); i++)
        {
            inputs[i] = {optimal_rn[i].m_x};
            optimals[i] = {optimal_rn[i].m_y};
        }
    }

    ///Evaluates the operformance of all indiivduals in a population
    std::vector<double> calculate_cumulative_performance_inds(){

        std::vector<double> cumulative_performance(get_inds().size(), 0);
        std::vector<std::vector<double>> performances(0, cumulative_performance);;

        std::vector<std::vector<double>> inputs;
        std::vector<double> optimals;


        if constexpr(Eval_type == evaluation_type::trial)
        {
            create_input_optimals_trial(inputs, optimals);
        }
        if constexpr(Eval_type == evaluation_type::full_rn)
        {
            create_inputs_optimal_full_rn(inputs, optimals);
        }

        performances.resize(inputs.size());

#pragma omp parallel for
        for(int i = 0; i < inputs.size(); i++)
        {
            performances[i] = pop::calc_dist_from_target(get_inds(),
                                                         optimals[i],
                                                         inputs[i]);
        }

        //Check this out!
        for(auto& performance : performances)
        {
            std::transform(cumulative_performance.begin(),
                           cumulative_performance.end(),
                           performance.begin(),
                           cumulative_performance.begin(),
                           std::plus<double>());
        }

        return cumulative_performance;
    }

    ///Calculate the performance as the cumulative performance (sum of all trials distances)
    /// divided by the number of trials (the performance is the mean distance per trial)
    std::vector<double> calculate_performances_inds(std::vector<double> cumulative_performances)
    {
        if constexpr(Eval_type == evaluation_type::trial)
        {
            for(auto& cumulative_performance : cumulative_performances)
            {
                cumulative_performance = cumulative_performance / m_population.get_n_trials();
                cumulative_performance = std::sqrt(cumulative_performance);
            }
        }
        if constexpr(Eval_type == evaluation_type::full_rn)
        {
            for(auto& cumulative_performance : cumulative_performances)
            {
                cumulative_performance = cumulative_performance / m_params.s_p.m_reac_norm_n_points;
                cumulative_performance = std::sqrt(cumulative_performance);
            }
        }
        return cumulative_performances;
    }

    ///Calculates fitness of inds in pop given current env values
    const simulation<Pop,
    Env_change_sym,
    Env_change_freq,
    Sel_Type,
    Adapt_per,
    Eval_type>&
    calc_fitness()
    {
        auto cumulative_performance = calculate_cumulative_performance_inds();

        auto performance = calculate_performances_inds(cumulative_performance);

        auto fitness_vector = pop::rescale_dist_to_fit(performance, get_sel_str());

        pop::set_fitness_inds(get_pop_non_const(), fitness_vector);

        return *this;
    }
    ///Reproduces inds to next gen based on their fitness
    void reproduce()
    {
        pop::reproduce(get_pop_non_const(), get_rng());
    }

    ///Reproduces inds to next gen randomly
    void reproduce_randomly()
    {
        pop::reproduce_random(get_pop_non_const(), get_rng());
    }

    ///Sorts indivudals in the vector of popoulation by fitness and assigns thema rank based on their position
    Pop sort_and_assign_ranks_by_fitness()
    {return m_population.sort_and_assign_ranks_by_fitness();}

    ///Makes the rank of the indivdual become the ancestor rank
    Pop ind_ID_becomes_ancestor_ID(std::vector<typename Pop::ind_t>& inds)
    {return m_population.ind_ID_becomes_ancestor_ID(inds);}

    ///Changes the ancestor rank for the current rank of the individual
    /// and then sorts the individuals by fitness and assigns new ranks
    Pop update_ancestor_rank_and_sort_by_new_rank()
    {
        sort_and_assign_ranks_by_fitness();
        ind_ID_becomes_ancestor_ID(get_new_inds_non_const());

        return m_population;
    }

    ///Calculates fitness and selects a new population based on fitness
    void select_inds()
    {
        if constexpr(Sel_Type == selection_type::sporadic)
        {
            if constexpr(Adapt_per == adaptation_period::on)
            {
                if(m_time < (m_n_generations /m_params.s_p.adaptation_period_proportion))
                {
                    calc_fitness();
                    reproduce();
                    return;
                }
            }

            if( m_selection_frequency != 0 &&
                    m_time % m_selection_frequency >= 0 &&
                    m_time % m_selection_frequency < m_selection_duration
                    )
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

    ///The current inputs that the networks of individuals will recieve in a given trial
    std::vector<double> m_input;

    ///The optimal output at a given trial;
    /// depends on inputs and environmental function
    double m_optimal_output;

    ///The series of inputs individuals will receive during one update of the simulation
    std::vector<std::vector<double>> m_stored_inputs;

    ///The optimal phenotypic outputs for one simulation update
    std::vector<double> m_stored_optimal_output;

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
    auto path = std::filesystem::current_path();
    path /= filename;
    if(std::filesystem::exists(path))
    {
        std::cout << "overriding previous results" << std::endl;
    }

    std::ofstream f(filename);
    if(f)
    {
        nlohmann::json json_out;
        json_out = s;
        f << json_out;
    }
    else
    {
        throw std::runtime_error{"could not open output stream for saving!"};
    }
    f.close();
}

namespace sim {


///Checks if 2 simulations are equal
template<class Pop>
bool operator ==(const simulation<Pop>& lhs, const simulation<Pop>& rhs);

///Check that all individuals have the same network
template<class Sim>
bool all_inds_have_same_net(const Sim& s)
{
    auto inds = s.get_pop().get_inds();
    return pop::all_inds_have_same_net(inds);
}

///Checks if all the individuals in a simulated population have the same input
template<class Sim>
bool all_individuals_have_same_input(const Sim &s)
{
    auto p = s.get_pop();

    return pop::all_individuals_have_same_input(p);
}

///Checks that all fitnesses in the population are not equal
template<class Sim>
bool all_fitnesses_are_not_equal(const Sim& s)
{
    return pop::all_fitnesses_are_not_equal(s.get_pop().get_inds());
}

///Checks that individuals all have a specific fitness value
template<class Sim>
bool all_inds_have_fitness(double value, const Sim& s)
{
    return pop::all_fitnesses_are(value, s.get_pop());
}

///Checks that the population in this simulation
/// has individuals sorted by rank
template<class Sim>
bool all_ranks_are_equal(const Sim& s)
{
    return pop::all_ranks_are_equal(s.get_pop().get_inds());
}

///Assign random ranks/Ids to indiviudal in a population
simulation<> assign_random_IDs_to_inds(simulation<> s, rndutils::xorshift128 &rng);

///Calculates the optimal output
template<class Sim>
double calculate_optimal(const Sim &s)
{
    return(env::calculate_optimal(s.get_env(), s.get_input()));
}

///Calculates the time to add to the seleciton frequency to record
/// data at the end of a selection period
template<class S>
int calculate_selection_duration(const S& o)
{
    int rec_freq_shift = 0;

    if(o.get_sel_type() == selection_type::sporadic &&
            o.get_e_change_f_type() == env_change_freq_type::regular)
    {
        rec_freq_shift += o.get_selection_duration();
    }

    return rec_freq_shift;
}


///Returns a population whose fitness has been calculated
template<class Sim>
typename Sim::pop_t calc_fitness_of_pop(Sim s)
{

    s.update_optimal(env::calculate_optimal(s.get_env(), s.get_input()));
    return pop::calc_fitness(s.get_pop_non_const(),
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

///Changes the network of the nth individual for a given network
template<class Pop>
void change_nth_ind_net(simulation<Pop>& s, size_t ind_index, const typename Pop::ind_t::net_t &n)
{
    pop::change_nth_ind_net(s.get_pop_non_const(), ind_index, n) ;
}

///Gets const ref the best n individuals in a pop
template<class Sim>
std::vector<typename Sim::pop_t::ind_t> get_best_n_inds(const Sim& s, int n)
{
    return pop::get_best_n_inds(s.get_pop(), n);
}

///Returns the current optimal function of the environment
template<class Sim>
std::function<double(std::vector<double>)> get_current_env_function(const Sim &s)
{
    auto e = s.get_env();
    return e.get_current_function();
}

///Gets the name of the current environmental function
template<class Sim>
char get_name_current_function(const Sim& s) noexcept
{
    return s.get_env().get_name_current_function();
}

///Returns the fitness of the nth ind in pop
template<class Sim>
double get_nth_ind_fitness(const Sim& s, const size_t ind_index);

///Returns the net of the nth individual in the population
template<class Sim>
const typename Sim::pop_t::ind_t::net_t & get_nth_ind_net(const Sim& s, size_t ind_index)
{
    return pop::get_nth_ind_net(s.get_pop(), ind_index);
}

///Gets the vector of individuals whihc were the parents of the current generation
template<class S>
const std::vector<typename S::pop_t::ind_t>& get_parents(const S& s)
{
    return s.get_pop().get_new_inds();
}

///Gets the vector of individuals whihc were the parents of the current generation
template<class S>
const std::vector<typename S::pop_t::ind_t>& get_inds(const S& s)
{
    return s.get_pop().get_inds();
}

///retruns the Ids of the parent population (m_vec_new_indiv)
template<class Ind>
std::vector<std::string> pop_IDs(const std::vector<Ind>& pop)
{
    std::vector<std::string> IDs(pop.size());
    std::transform(pop.begin(), pop.end(), IDs.begin(),
                   [](const Ind& ind){return ind.get_rank();});
    return IDs;}

///retruns the ancestor Ids of the population
template<class Ind>
std::vector<std::string> ancestor_IDs(const std::vector<Ind>& pop)
{
    std::vector<std::string> IDs(pop.size());
    std::transform(pop.begin(), pop.end(), IDs.begin(),
                   [](const Ind& ind){return ind.get_ancestor_ID();});
    return IDs;
}

///Checks that the ancestor IDs of inds in a population
/// are equal to their parents IDs (inds store in new_inds_vec
bool ancestor_ID_is_parent_ID(const simulation<>& s);

///Checks that the population in this simulation
/// has individuals sorted by fitness
template<class Sim>
bool is_sorted_by_fitness(const Sim& s)
{
    return pop::is_sorted_by_fitness(s.get_pop().get_inds());
}

///Checks that the population in this simulation
/// has individuals sorted by rank
template<class Sim>
bool is_sorted_by_rank(const Sim& s)
{
    return pop::is_sorted_by_rank(s.get_pop().get_inds());
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
