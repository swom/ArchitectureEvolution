#ifndef OBSERVER_H
#define OBSERVER_H
#include "simulation.h"
#include "Stopwatch.hpp"

static int sample_ind_record_freq_to_sens = 10;

struct sensibilities_to_mut
{
    sensibilities_to_mut(int gen = -1, std::vector<fit_and_phen_sens_t> sensibilities = {}):
        m_generation(gen),
        m_sensibilities(sensibilities)
    {}

    NLOHMANN_DEFINE_TYPE_INTRUSIVE(sensibilities_to_mut,
                                   m_generation,
                                   m_sensibilities)
    int m_generation;
    std::vector<fit_and_phen_sens_t> m_sensibilities;
    size_t size() const noexcept {return m_sensibilities.size();}
};

bool operator== (const sensibilities_to_mut& lhs, const sensibilities_to_mut& rhs);
bool operator!= (const sensibilities_to_mut& lhs, const sensibilities_to_mut& rhs);

struct inputs_optimals{
    inputs_optimals(std::vector<std::vector<double>> inputs = {},
                    std::vector<double> optimals = {},
                    int gen = 0):
        m_inputs{inputs},
        m_optimals{optimals},
        m_gen{gen}
    {}

    NLOHMANN_DEFINE_TYPE_INTRUSIVE(inputs_optimals,
                                   m_inputs,
                                   m_optimals,
                                   m_gen)

    ///The inputs given at a certain generation
    std::vector<std::vector<double>> m_inputs;
    ///The optimal values corresponding to the inputs given at a certain generation
    std::vector<double> m_optimals;
    ///The generation at whihc the data is collected
    int m_gen;
};

bool operator== (const inputs_optimals& lhs, const inputs_optimals& rhs);
bool operator!= (const inputs_optimals& lhs, const inputs_optimals& rhs);
///Calculates the average robustness of a population in a simulation
template<class Sim>
double calc_avg_mutation_sensibility(Sim& s, int n_mutations)
{
    return calc_mean(pop::calc_mutation_sensibility_all_inds(s.get_pop(), n_mutations, s.get_rng()));
}


struct obs_param{
    obs_param(int top_prop = 1,
              int top_ind_reg_freq = 1,
              int spectrum_reg_freq = 0,
              int n_data_points_for_reac_norm = 100,
              int n_mutations_for_mutational_spectrum = 1):
        m_top_proportion{top_prop},
        m_top_ind_reg_freq{top_ind_reg_freq},
        m_spectrum_reg_freq{spectrum_reg_freq},
        m_reac_norm_n_points{n_data_points_for_reac_norm},
        m_n_mutations_per_locus{n_mutations_for_mutational_spectrum}
    {
        if(top_prop == 0)
            throw std::invalid_argument{"the number of indidivuduals recorderd cannot be 0"};
    }

    NLOHMANN_DEFINE_TYPE_INTRUSIVE(obs_param,
                                   m_top_proportion,
                                   m_top_ind_reg_freq,
                                   m_spectrum_reg_freq,
                                   m_reac_norm_n_points,
                                   m_n_mutations_per_locus
                                   )

    ///The top n idividuals stored in the observed
    int m_top_proportion;
    ///The number of generations after which individuals are stored
    int m_top_ind_reg_freq;
    ///The number of generations after which
    /// the mutational spectrum of the top individuals is stored
    int m_spectrum_reg_freq;
    ///The number of data points on which to calculate the reaction norm of an individual
    int m_reac_norm_n_points;
    ///The number of mutations executed on each locus
    /// when calculating the mutational spectrum of an individual
    int m_n_mutations_per_locus;
};

template<class Ind>
struct Ind_Data
{
    Ind_Data(){};
    Ind_Data(Ind ind, reac_norm rn,  fit_and_phen_sens_t sensibilities, int gen):
        m_ind{ind},
        m_reac_norm{rn},
        m_sensibilities{sensibilities},
        generation{gen}
    {}

    NLOHMANN_DEFINE_TYPE_INTRUSIVE(Ind_Data,
                                   m_ind,
                                   m_reac_norm,
                                   m_sensibilities,
                                   generation)

    Ind m_ind;
    reac_norm m_reac_norm;
    fit_and_phen_sens_t m_sensibilities;
    int generation;
};

template<class Ind>
bool operator== (const Ind_Data<Ind>& lhs, const Ind_Data<Ind>& rhs)
{
    auto ind = lhs.m_ind == rhs.m_ind;
    auto reaction_norm = lhs.m_reac_norm == rhs.m_reac_norm;
    return ind && reaction_norm;
}

///Finds the corresponding sensibility of an individual
/// based on its fitness ranking
template<class Ind>
const fit_and_phen_sens_t& find_sensibility_for_ind(const Ind& ind,
                                                    const sensibilities_to_mut& sensibilities)
{
    auto rank = ind.get_rank();

    const fit_and_phen_sens_t& sens = *std::find_if(
                sensibilities.m_sensibilities.begin(),
                sensibilities.m_sensibilities.end(),
                [rank](const fit_and_phen_sens_t& sensibility){return sensibility.m_rank == rank;});

    return sens;
}

///Finds the corresponding individual for its sensibilities
/// based on its fitness ranking
template<class Ind>
const Ind& find_ind_for_sensibility(const std::vector<Ind>& inds,
                                    const fit_and_phen_sens_t& sens)
{
    auto rank = sens.m_rank;

    const Ind& ind = *std::find_if(
                inds.begin(),
                inds.end(),
                [rank](const Ind& ind){return ind.get_rank() == rank;});

    return ind;
}
///Finds the best possible combination of fitness and phenotypic sensibility
/// in a vector of sensibilities
/// the best combination possible having fitness sensibility = the max fitness sensibility
/// and the phenotypic sensibility = 0
fit_and_phen_sens_t find_best_fit_phen_combination(const sensibilities_to_mut& record);

///Calculates the reaction_norm of individuals' networks
/// for a given range and a given number of data points
template<class Ind, typename Func>
std::vector<Ind_Data<Ind>> create_inds_data(const std::vector<Ind>& inds,
                                            const sensibilities_to_mut& sensibilities,
                                            const obs_param& o_params,
                                            const all_params& a_params,
                                            const int& generation,
                                            Func optimal_function)
{

    std::vector<Ind_Data<Ind>> inds_data(inds.size());
    for(auto i = 0; i != inds.size(); i++)
    {
        inds_data[i].m_ind = inds[i];

        inds_data[i].m_reac_norm = calculate_reaction_norm(inds[i].get_net(),
                                                           a_params.e_p.cue_range,
                                                           o_params.m_reac_norm_n_points);

        inds_data[i].m_sensibilities = find_sensibility_for_ind(inds[i], sensibilities);

        inds_data[i].generation = generation;
    }
    return inds_data;
}

class base_observer{
public:

};

template<class Sim = simulation<>>
class observer : public base_observer
{
public:
    using Sim_t = Sim;
    using Pop = typename Sim_t::pop_t;
    using Ind = typename Pop::ind_t;
    using Net = typename Ind::net_t;
    using Net_Spect = network_spectrum<Net>;


private:

    std::vector<double> m_avg_fitnesses;
    std::vector<double> m_avg_mutation_sensibility;
    std::vector<sensibilities_to_mut> m_fit_phen_mut_sensibility;
    std::vector<double> m_var_fitnesses;
    std::vector<char> m_env_functions;
    std::vector<std::vector<Ind_Data<Ind>>> m_sampled_inds;
    std::vector<std::vector<Ind_Data<Ind>>> m_top_inds;
    std::vector<std::vector<Ind_Spectrum<Ind>>> m_top_spectrums;
    obs_param m_obs_param;
    all_params m_params;
    std::vector<inputs_optimals> m_inputs_optimals;

public:

    observer(
            obs_param params = obs_param{},
            all_params sim_params = all_params{}
            ):
        m_obs_param(params),
        m_params{sim_params}
    {

    }


    NLOHMANN_DEFINE_TYPE_INTRUSIVE(observer,
                                   m_avg_fitnesses,
                                   m_var_fitnesses,
                                   m_top_inds,
                                   m_top_spectrums,
                                   m_env_functions,
                                   m_params,
                                   m_obs_param,
                                   m_inputs_optimals,
                                   m_avg_mutation_sensibility,
                                   m_fit_phen_mut_sensibility,
                                   m_sampled_inds)

    ///adds a network spectrum to the vector of network spectrums
    void add_spectrum(const std::vector<Ind_Spectrum<Ind>>& spectrum){ m_top_spectrums.push_back(spectrum);}

    ///returns const ref to m_avg_fitness
    const std::vector<double>& get_avg_fitness()  const noexcept{return m_avg_fitnesses;}

    ///Returns the vector containing the avg robustness of the population for every generation
    const std::vector<double>& get_avg_mutation_sensibility() const noexcept {return m_avg_mutation_sensibility;}

    ///Returns the vector containing the fitness mutational sensibility of the population for every generation
    const std::vector<sensibilities_to_mut>& get_fit_phen_mut_sensibility() const noexcept {return m_fit_phen_mut_sensibility;}

    ///Returns the vector containing the fitness mutational sensibility of the population for every generation
    std::vector<sensibilities_to_mut>& get_fit_phen_mut_sensibility_non_const()  noexcept {return m_fit_phen_mut_sensibility;}

    ///Returns the vector containing the fitness mutational sensibility of the population for every generation
    const sensibilities_to_mut& get_first_fit_phen_mut_sens() const noexcept {return m_fit_phen_mut_sensibility.at(0);}

    ///returns const ref to vector of env_functions' names
    const std::vector<char>& get_env_funcs() const noexcept {return m_env_functions;}

    ///Returns the type of enviromental change frequency used in simulation
    auto get_e_change_f_type() const noexcept{return m_params.s_p.change_freq_type;}

    ///Returns the selection_duration used in simulation
    int get_selection_duration() const noexcept {return m_params.s_p.selection_duration;}

    ///Returns the type of selection used in simulation
    auto get_sel_type() const noexcept{return m_params.s_p.sel_type;}

    ///returns the inputs and optimals
    const std::vector<inputs_optimals>& get_inputs_and_optimals() const noexcept {return m_inputs_optimals;}

    ///Returns const ref to record frequency
    /// of top individuals
    const int& get_record_freq_top_inds() const noexcept {return m_obs_param.m_top_ind_reg_freq;}

    ///Returns const ref to record frequency
    /// of top indidivudals mutational spectrums
    const int& get_record_freq_spectrum() const noexcept {return m_obs_param.m_spectrum_reg_freq;}

    ///Returns the data about the last recorded population of individuals
    const std::vector<fit_and_phen_sens_t>& get_last_recorded_pop() const noexcept
    {return get_fit_phen_mut_sensibility().back().m_sensibilities;}

    /// Returns const ref to the number of points to be recorded for reaction norm
    const int& get_n_points_reac_norm() const noexcept {return m_obs_param.m_reac_norm_n_points;}

    /// Returns const ref to the number of mutations
    /// used for calcualting the mutational spectrum
    const int& get_n_mut_mutational_spectrum() const noexcept {return m_obs_param.m_n_mutations_per_locus;}

    ///returns const ref to m_var_fitnesses
    const std::vector<double>& get_var_fitness() const noexcept{return m_var_fitnesses;}

    ///returns const ref to the number of top individuals to be recorded
    const int& get_top_inds_proportion() const noexcept{return m_obs_param.m_top_proportion;}

    ///Returns the vector whos generation element is of a certain value
    template<class ind_data_structure>
    const std::vector<ind_data_structure>& get_generation(const std::vector<std::vector<ind_data_structure>>& data, int generation)
    {
        auto data_vec = std::find_if(data.begin(),
                                     data.end(),
                                     [&generation](const auto& inds_vec){return inds_vec.at(0).generation == generation;});
        if(data_vec != data.end())
            return *data_vec;
        else
        {
            std::string error = {"No record for the given generation for the given data_type: "};
            error += typeid(ind_data_structure).name();
            throw std::invalid_argument{error};
        }
    }

    ///Retruns a vector of individuals from a vector of individual data structures
    template<class ind_data_structure>
    std::vector<Ind> extract_inds(std::vector<ind_data_structure> data_v) noexcept
    {
        std::vector<Ind> inds(data_v.size());
        std::transform(data_v.begin(), data_v.end(), inds.begin(),
                       [](const auto& data_structure){return data_structure.m_ind;});
        return inds;
    }

    ///Returns the vector of sampled individuals
    const std::vector<std::vector<Ind_Data<Ind>>>& get_sampled_inds() const noexcept {return m_sampled_inds;}

    ///Returns a vector containing the stored top idnividuals of a given gen
    std::vector<Ind> get_top_inds_gen(int generation) noexcept
    {
        return extract_inds(get_generation(m_top_inds, generation));
    }

    ///returns const ref to best_ind vector
    const std::vector<std::vector<Ind_Data<Ind>>>& get_top_inds() const noexcept{return m_top_inds;}

    /// Returns the first top_individual from the first recorded batch
    const Ind_Data<Ind>& get_first_top_ind_of_first_gen() const noexcept {return m_top_inds.at(0).at(0);}

    /// of the best individuals in various generations
    ///returns const ref to the vector of mutational spectrum
    const std::vector<std::vector<Ind_Spectrum<Ind>>>& get_top_spectrums() const noexcept{return m_top_spectrums;}

    ///returns the time of the simulation taking in consioderation that tickhas just been run
    /// and that therefore the time count is up by one
    int get_time_before_tick(const Sim& s) const noexcept {return s.get_time() - 1;}

    ///returns const ref to the vector of mutational spectrum
    /// from individuals of a particular generation
    const std::vector<Ind_Spectrum<Ind>>& get_top_spectrums_gen(int generation) noexcept
    {
        return get_generation(m_top_spectrums, generation);
    }

    ///Returns a const referernce to the paramteres used to initialize the simulation
    const all_params& get_params() const noexcept {return m_params;};

    //    //Returns a const reference to the input vector given to individuals in the current or last trial
    //    const std::vector<std::vector<std::vector<double>>>& get_input() const noexcept {return m_input;}

    //    ////returns a constant reference to the otpimal output value given to individuals that generation
    //    const std::vector<std::vector<double>>& get_optimal() const noexcept {return m_optimal;}

    void store_avg_mut_sensibility(Sim& s) noexcept
    {
        m_avg_mutation_sensibility.push_back(calc_avg_mutation_sensibility(s, m_obs_param.m_n_mutations_per_locus));
    }


    ///Stores the data from simulation at the coorect times
    void store_data(Sim& s)
    {
        store_env_func(s);
        store_var_fit(s);
        store_ind_data(s);
    }

    ///Stores data related to the individuals that just went through tick
    /// for this reason it swaps the individual vectors in population back
    /// as they were after calc_fitness() and before reproduce()
    /// in simulation tick()
    void store_ind_data(Sim& s)
    {
        int selection_duration = sim::calculate_selection_duration(s);

        ///To look at generation that went through tick
        /// we need swap back the individuals vectors that were swapped during reproduce()

        pop::swap_new_with_old_pop(s.get_pop_non_const()); ///Swap in

        store_avg_fit(s);

        if(is_time_to_record_inds(*this,s))
        {
            store_data_based_on_sensibilities(s);
            store_inputs_and_optimals(s);
        }

        if(is_time_to_record_best_inds_mut_spectrum(*this, s, selection_duration))
        {
            store_network_spectrum_n_best(s);
        }

        pop::swap_new_with_old_pop(s.get_pop_non_const()); ///Swap out
    }

    ///Saves the avg fitness
    void store_avg_fit(const Sim &s)
    {
        m_avg_fitnesses.push_back(sim::avg_fitness(s));
    }

    ///Stores the inputs and optimals provided to indivudals and records also in which genration
    void store_inputs_and_optimals(const Sim& s){m_inputs_optimals.push_back(inputs_optimals{s.get_stored_inputs(),
                                                                                             s.get_stored_optimals(),
                                                                                             s.get_time()});}
    ///Saves the variance of the fitness
    void store_var_fit(const Sim& s)
    {
        m_var_fitnesses.push_back(sim::var_fitness(s));
    }

    ///Stores the sensibilities and the top individuals toghether
    /// the sensibilities are stored first so to eanble the assignment of
    /// the correct sensibility to the correct top individual
    void store_data_based_on_sensibilities(Sim& s)
    {
        if(s.get_inds().empty()) return;


        store_fit_phen_mut_sensibility(s);
        store_top_n_inds(s);
        if(m_fit_phen_mut_sensibility.size() % sample_ind_record_freq_to_sens == 0)
        {
            store_top_mid_low_sens_inds(s);
        }
    }

    ///Stores the network spectrum of the top n best individuals
    void store_network_spectrum_n_best(Sim& s)
    {

        m_top_spectrums.emplace_back(
                    calculate_mut_spectrums(sim::get_best_n_inds(s, get_top_inds_proportion()),
                                            s.get_mut_step(),
                                            s.get_rng(),
                                            get_n_mut_mutational_spectrum(),
                                            s.get_env_cue_range(),
                                            get_n_points_reac_norm(),
                                            get_time_before_tick(s) //when this function is called timer is ticked,
                                            //but we are still in the previous generation
                                            )
                    );
    }

    ///Calculates the mutational spectrums of individuals
    ///stored in top inds vector, given a certain generation
    /// and adds it to the vector of top individuals
    std::vector<Ind_Spectrum<Ind>> calculate_mut_spectrums_for_gen(int generation)
    {
        std::mt19937_64 rng;
        std::vector<Ind_Spectrum<Ind>> spectrums = calculate_mut_spectrums(get_top_inds_gen(generation),
                                                                           m_params.p_p.mut_step,
                                                                           rng,
                                                                           get_n_mut_mutational_spectrum(),
                                                                           m_params.e_p.cue_range,
                                                                           get_n_points_reac_norm(),
                                                                           generation);
        return spectrums;
    }

    void store_env_func (const Sim& s) noexcept {m_env_functions.push_back(sim::get_name_current_function(s));}

    ///Samples 3 individuals based on their mutational sensibilities
    /// The one with the highest, lowest and most median
    void store_top_mid_low_sens_inds(const Sim& s)
    {
        ///Make sure that the sensibilities have already been recorded
        if(get_fit_phen_mut_sensibility().size() < m_sampled_inds.size() + 1)
        {
            throw std::runtime_error{"the sensibilities have not yet been recorded, "
                                     "therefore the top individuals cannot be assigned "
                                     "a sensibility"};
        }
        m_sampled_inds.push_back(create_inds_data(sample_top_mid_low_sens_inds(s, *this),
                                                  m_fit_phen_mut_sensibility.back(),
                                                  m_obs_param,
                                                  m_params,
                                                  get_time_before_tick(s), //when this function is called timer is ticked, but we are still in the previous generation
                                                  sim::get_current_env_function(s)));
    }
private :
    ///Stores the mutational sensibilities to fitness and phenotype of all individuals in the population
    void store_fit_phen_mut_sensibility(Sim& s) noexcept
    {
        m_fit_phen_mut_sensibility.push_back(
                    {get_time_before_tick(s),
                     s.calculate_fit_phen_mut_sens_for_all_inds(
                     m_obs_param.m_n_mutations_per_locus,
                     m_obs_param.m_reac_norm_n_points)}
                    );
    }

    ///Saves the top_proportion nth best individuals in the population
    void store_top_n_inds(Sim& s)
    {
        if(get_fit_phen_mut_sensibility().size() < m_top_inds.size() + 1)
        {
            throw std::runtime_error{"the sensibilities have not yet been recorded, "
                                     "therefore the top individuals cannot be assigned "
                                     "a sensibility"};
        }

        m_top_inds.push_back(create_inds_data(sim::get_best_n_inds(s, get_top_inds_proportion()),
                                              m_fit_phen_mut_sensibility.back(),
                                              m_obs_param,
                                              m_params,
                                              get_time_before_tick(s), //when this function is called timer is ticked, but we are still in the previous generation
                                              sim::get_current_env_function(s)
                                              )
                             );
    }
};

template<class Ind>
bool operator==(const observer<Ind>& lhs, const observer<Ind>& rhs);

bool operator==(const all_params& lhs, const all_params& rhs);

bool operator!=(const all_params& lhs, const all_params& rhs);

///Creates a very simple simulation with 1 individual
/// that has a sigmoid transformation function
/// that has 1 connection and one bias
/// with input range == 1
/// with env function == y = 1
simulation<> create_simple_simulation(int n_gen = 1);

///Creates a unique saving name based on the parameters
std::string create_save_name_from_params(const all_params& p);

///Calculates the euclidean distance from the best combination of sensibilities in a vector
/// to a given set of sesnsibilites s
double squared_distance_from_best_combination(const fit_and_phen_sens_t& ind,
                                              const fit_and_phen_sens_t& best);

///Calculates the distance of a given individual sensibilities
/// from the best possible sensibilities combination of the last recorded generation
template<class Ind>
double distance_from_best_sens_comb(const Ind& i,
                                    const sensibilities_to_mut& last_gen)
{
    auto ind_sens =find_sensibility_for_ind(i, last_gen);
    auto best = find_best_fit_phen_combination(last_gen);

    return squared_distance_from_best_combination(ind_sens, best);
}

///Check if it is the time at the end of a selection period
template<class O, class S>
bool is_end_of_selection_period_and_time_to_record(const O& o, const S& s)
{
    return s.get_sel_type() == selection_type::sporadic &&
            o.get_record_freq_top_inds() != 0 &&
            (s.get_time() - s.get_selection_duration()) %  o.get_record_freq_top_inds() == 0;
}

///Check if it is the time at the start of a selection period
template<class O, class S>
bool is_before_start_of_selection_period_and_time_to_record(const O& o, const S& s)
{
    return o.get_record_freq_top_inds() != 0 &&
            (s.get_time()) %  o.get_record_freq_top_inds() == 0;
}

///Check if it is time to record
template<class O, class S>
bool is_time_to_record_inds(const O& o, const S& s)
{
    return is_end_of_selection_period_and_time_to_record(o, s) |
            is_before_start_of_selection_period_and_time_to_record(o, s);
}

///Check if it is the time to record the best individuals' mutation spectrum
template<class O, class S>
bool is_time_to_record_best_inds_mut_spectrum(const O& o, const S& s, int rec_freq_shift)
{
    return o.get_record_freq_spectrum() != 0 &&
            (s.get_time() - rec_freq_shift) % o.get_record_freq_spectrum() == 0;
}

///?Checks that an observer and a simulation haev same parameters and if not throws an exception
template<class O, class S>
bool throw_if_obs_and_sim_do_not_have_same_param(const O& o, const S& s)
{

    try
    {
        if(s.get_params() != o.get_params())
        {
            throw std::runtime_error{"During exec(): /n "
                                     "Observer was not initialized "
                                     "correctly with simulation parameters /n"};
        }
        return false;
    }
    catch (std::runtime_error& e)
    {
        std::cerr << e.what();
    }
    return true;
}

template<class Sim>
void print_elapsed_time_every_n_gen(const Sim& s, int n, stopwatch::Stopwatch& my_watch)
{
    if(s.get_time() % n == 0)
    {
        auto lap_ms = my_watch.lap<stopwatch::ms>();
        std::cout << "Cycle " << s.get_time() << " --Lap time in ms: " << lap_ms << std::endl;
    }
}


///Executes a simulation for n generations
template<class Sim>
void exec(Sim& s , observer<Sim>& o)
{
    namespace sw = stopwatch;
    sw::Stopwatch my_watch;

    throw_if_obs_and_sim_do_not_have_same_param(o,s);

    while(s.get_time() !=  s.get_n_gen())
    {
        sim::tick(s);

        o.store_data(s);

        print_elapsed_time_every_n_gen(s, 1000, my_watch);
    }
}

///Retruns a vector of the generations in which individuals were recorded
template<class ind_data_structure>
std::vector<int> extract_gens(std::vector<ind_data_structure> data_v) noexcept
{
    std::vector<int> gens(data_v.size());
    std::transform(data_v.begin(), data_v.end(), gens.begin(),
                   [](const auto& data_structure){return data_structure.at(0).generation;});
    return gens;
}
///Returns the sensibilities recorded from the individuals
/// stored in the first time the sensibilites are recorded
sensibilities_to_mut get_inds_sensibilities_of_first_record(const observer<>& o);

///finds the sensibilities of the hihgest ranking individual in a given generation
fit_and_phen_sens_t find_sensibilities_from_highest_ranking_ind(const sensibilities_to_mut& record);

///Returns the first top individual
/// recorded in the first time top individuals are stores
template<class O>
Ind_Data<typename O::Ind> get_first_top_ind_of_first_record(const O& obs)
{
    return obs.get_top_inds().at(0).at(0);
}

///Gets the inputs  of the nth generation
template<class O>
const std::vector<std::vector<double>>& get_nth_gen_inputs(const O& o, int gen)
{
    auto generation = std::find_if(o.get_inputs_and_optimals().begin(),
                                   o.get_inputs_and_optimals().end(),
                                   [&gen] (const inputs_optimals& io) {return io.m_gen == gen;});
    if( generation == o.get_inputs_and_optimals().end())
    {
        throw std::invalid_argument{"In observer, requested to acces inputs "
                                    "from a generation that has no recorded"
                                    " inputs or that does not exist"};
    }

    return generation->m_inputs;
}

///Gets the optimals of the nth generation
template<class O>
const std::vector<double>& get_nth_gen_optimals(const O& o, int gen)
{

    auto generation = std::find_if(o.get_inputs_and_optimals().begin(),
                                   o.get_inputs_and_optimals().end(),
                                   [&gen] (const inputs_optimals& io) {return io.m_gen == gen;});

    if(generation == o.get_inputs_and_optimals().end())
    {
        throw std::invalid_argument{"In observer, requested to acces optimals "
                                    "from a generation that has no recorded"
                                    " inputs or that does not exist"};
    }

    return generation->m_optimals;
}

///Checks if one observer has
/// as lower phenotypic mutation sensibility values stored
/// than another
template<class Obs>
bool lhs_has_lower_phen_mutation_sensibility_than_rhs(const Obs& lhs, const Obs& rhs)
{
    assert(lhs.get_fit_phen_mut_sensibility().size() == rhs.get_fit_phen_mut_sensibility().size());

    for(int i = 0; i != lhs.get_fit_phen_mut_sensibility().size(); i++)
    {
        assert(lhs.get_fit_phen_mut_sensibility()[i].size() == rhs.get_fit_phen_mut_sensibility()[i].size());

        for(int j = 0; j != lhs.get_fit_phen_mut_sensibility()[i].size(); j++)
        {
            if(lhs.get_fit_phen_mut_sensibility()[i].m_sensibilities[j].m_phenotype_sens >=
                    rhs.get_fit_phen_mut_sensibility()[i].m_sensibilities[j].m_phenotype_sens)
            {
                return false;
            }
        }
    }
    return true;
}

///Checks that all mutation fitness sensibilities
/// are negative
template<class Obs>
bool all_fit_mut_sens_are_negative(const Obs& observer)
{
    for(const auto& sensibilities : observer.get_fit_phen_mut_sensibility())
    {
        for(const auto& sensibility : sensibilities.m_sensibilities)
        {
            if(sensibility.m_fitness_sens >= 0)
            {
                return false;
            }
        }
    }
    return true;
}

///Checks if one observer has
/// as higher fitness mutational sensibility values stored
/// than another
template<class Obs>
bool lhs_is_more_sensible_to_mutation_effects_on_fitness_than_rhs(const Obs& lhs, const Obs& rhs)
{
    assert(lhs.get_fit_phen_mut_sensibility().size() == rhs.get_fit_phen_mut_sensibility().size());

    for(int i = 0; i != lhs.get_fit_phen_mut_sensibility().size(); i++)
    {
        assert(lhs.get_fit_phen_mut_sensibility()[i].size() == rhs.get_fit_phen_mut_sensibility()[i].size());

        for(int j = 0; j != lhs.get_fit_phen_mut_sensibility()[i].size(); j++)
        {
            if(lhs.get_fit_phen_mut_sensibility()[i].m_sensibilities[j].m_fitness_sens >=
                    rhs.get_fit_phen_mut_sensibility()[i].m_sensibilities[j].m_fitness_sens)
            {
                return false;
            }
        }
    }
    return true;
}

///load an observer of correct type based on selection frequency type
template<mutation_type M, selection_type S, env_change_freq_type F>
std::unique_ptr<base_observer> load_observer_json_change_symmetry_type(const all_params& pars)
{
    auto f_t = pars.s_p.change_sym_type;
    if (f_t == env_change_symmetry_type::symmetrical ) {
        using net_t = network<M>;
        using ind_t = individual<net_t>;
        using pop_t = population<ind_t>;
        using sim_t = simulation<pop_t, env_change_symmetry_type::symmetrical, F, S>;
        return std::make_unique<base_observer>(observer<sim_t>());
    }
    else if(f_t == env_change_symmetry_type::symmetrical)
    {
        using net_t = network<M>;
        using ind_t = individual<net_t>;
        using pop_t = population<ind_t>;
        using sim_t = simulation<pop_t, env_change_symmetry_type::asymmetrical, F, S>;
        return std::make_unique<base_observer>(observer<sim_t>());
    }
    else
        throw std::invalid_argument{"invalid symmetry type when loading observer"};


}

///load an observer of correct type based on change frequency type
template<mutation_type M, selection_type S>
std::unique_ptr<base_observer> load_observer_json_change_freq_type(const all_params& pars)
{
    auto f_t = pars.s_p.change_freq_type;
    switch (f_t) {
    case env_change_freq_type::regular :
        return load_observer_json_change_symmetry_type<M, S, env_change_freq_type::regular>(pars);
        break;
    case env_change_freq_type::stochastic :
        return load_observer_json_change_symmetry_type<M, S, env_change_freq_type::stochastic>(pars);
        break;
    default:
        throw std::invalid_argument{"invalid frequency type when loading observer"};

    }
}

///load an observer of correct type based on selection type
template<mutation_type M>
std::unique_ptr<base_observer> load_observer_selection_type(const all_params& pars)
{
    auto s_t = pars.s_p.sel_type;
    switch (s_t) {
    case selection_type::constant :
        return load_observer_json_change_freq_type<M, selection_type::constant>(pars);
        break;
    case selection_type::sporadic :
        return load_observer_json_change_freq_type<M, selection_type::sporadic>(pars);
        break;
    default:
        throw std::invalid_argument{"invalid selection_type when loading observer"};

    }
}
///load an observer of correct type based on parameter
std::unique_ptr<base_observer> load_observer_json_of_correct_type(const all_params &pars);
///loads an observer based on filename
std::unique_ptr<base_observer> load_observer_json(const std::string& filename);

///Samples three individual
/// the closest to the best combination of sesnibilities achievable in that generation
/// the median individual in the population
/// and the farthest from the best combination of sesnibilities achievable in that generation
template<class O>
std::vector<typename O::Ind> sample_top_mid_low_sens_inds(const typename O::Sim_t& s, O &o)
{

    auto& last_gen_record = o.get_fit_phen_mut_sensibility_non_const().back();
    auto best_comb = find_best_fit_phen_combination(last_gen_record);
    auto& last_gen_sens = last_gen_record.m_sensibilities;

    std::sort(last_gen_sens.begin(), last_gen_sens.end(),
              [best_comb]
              (const fit_and_phen_sens_t& lhs,
              const fit_and_phen_sens_t& rhs)
    {return squared_distance_from_best_combination(lhs, best_comb) <
                squared_distance_from_best_combination(rhs, best_comb);}
    );

    auto top_ind = find_ind_for_sensibility(s.get_inds(), *last_gen_sens.begin());
    auto low_ind = find_ind_for_sensibility(s.get_inds(), last_gen_sens.back());
    auto mid_ind = find_ind_for_sensibility(s.get_inds(), last_gen_sens[last_gen_sens.size() / 2]);

    return {top_ind, mid_ind, low_ind};
}

void test_observer();

#endif // OBSERVER_H
