#ifndef OBSERVER_H
#define OBSERVER_H
#include "simulation.h"
#include "Stopwatch.hpp"

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

struct obs_param{
    obs_param(int top_prop = 1,
              int top_ind_reg_freq = 1,
              int spectrum_reg_freq = 0,
              int n_data_points_for_reac_norm = 100,
              int n_mutations_for_mutational_spectrum = 1000):
        m_top_proportion{top_prop},
        m_top_ind_reg_freq{top_ind_reg_freq},
        m_spectrum_reg_freq{spectrum_reg_freq},
        m_reac_norm_n_points{n_data_points_for_reac_norm},
        m_n_mutations_per_locus{n_mutations_for_mutational_spectrum}
    {
        if(top_ind_reg_freq == 0)
            throw std::invalid_argument{"the number of generations after which a recording should happen cannot be 0"};
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

class base_observer{
public:

};

template<class Sim = simulation<>>
class observer : public base_observer
{
    using Sim_t = Sim;
    using Pop = typename Sim_t::pop_t;
    using Ind = typename Pop::ind_t;
    using Net = typename Ind::net_t;
    using Net_Spect = network_spectrum<Net>;


private:

    std::vector<double> m_avg_fitnesses;
    std::vector<double> m_var_fitnesses;
    std::vector<char> m_env_functions;
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
                                   m_inputs_optimals)

    ///adds a network spectrum to the vector of network spectrums
    void add_spectrum(const std::vector<Ind_Spectrum<Ind>>& spectrum){ m_top_spectrums.push_back(spectrum);}

    ///returns const ref to m_avg_fitness
    const std::vector<double>& get_avg_fitness()  const noexcept{return m_avg_fitnesses;}

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

    /// Returns const ref to the number of points to be recorded for reaction norm
    const int& get_n_points_reac_norm() const noexcept {return m_obs_param.m_reac_norm_n_points;}

    /// Returns const ref to the number of mutations
    /// used for calcualting the mutational spectrum
    const int& get_n_mut_mutational_spectrum() const noexcept {return m_obs_param.m_n_mutations_per_locus;}

    ///returns const ref to m_var_fitnesses
    const std::vector<double>& get_var_fitness() const noexcept{return m_var_fitnesses;}

    ///returns const ref to the number of top individuals to be recorded
    const int& get_top_inds_proportion() const noexcept{return m_obs_param.m_top_proportion;}

    ///Returns the vector whosgeneration element is of a certain value
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


    ///Returns a vector containing the stored top idnividuals of a given gen
    std::vector<Ind> get_top_inds_gen(int generation) noexcept
    {
        return extract_inds(get_generation(m_top_inds, generation));
    }
    ///returns const ref to best_ind vector
    const std::vector<std::vector<Ind_Data<Ind>>>& get_top_inds() const noexcept{return m_top_inds;}

    /// of the best individuals in various generations
    ///returns const ref to the vector of mutational spectrum
    const std::vector<std::vector<Ind_Spectrum<Ind>>>& get_top_spectrums() const noexcept{return m_top_spectrums;}

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

    ///Saves the top_proportion nth best individuals in the population
    void store_top_n_inds(const Sim& s)
    {
        m_top_inds.push_back(calculate_reaction_norms(sim::get_best_n_inds(s, get_top_inds_proportion()),
                                                      s.get_env_cue_range(),
                                                      get_n_points_reac_norm(),
                                                      s.get_time() - 1 //when this function is called timer is ticked,
                                                      //but we are still in the previous generation
                                                      )
                             );
    }

    ///Saves the nth best individuals in the population
    void store_top_n_inds(const Sim& s, int proportion)
    {
        m_top_inds.push_back(calculate_reaction_norms(
                                 sim::get_best_n_inds(s, proportion),
                                 s.get_env_cue_range(),
                                 get_n_points_reac_norm(),
                                 s.get_time() - 1 //when this function is called timer is ticked,
                                 //but we are still in the previous generation
                                 )
                             );
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
                                            s.get_time() - 1 //when this function is called timer is ticked,
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

    //    void store_input(const Sim& s) noexcept {m_input.push_back(s.get_stored_inputs());}

    //    void store_optimal(const Sim& s) noexcept {m_optimal.push_back(s.get_stored_optimals());}
};

template<class Ind>
bool operator==(const observer<Ind>& lhs, const observer<Ind>& rhs);

bool operator==(const all_params& lhs, const all_params& rhs);

bool operator!=(const all_params& lhs, const all_params& rhs);

///Creates a unique saving name based on the parameters
std::string create_save_name_from_params(const all_params& p);

///Executes a simulation for n generations
template<class Sim>
void exec(Sim& s , observer<Sim>& o)
{
    if(s.get_params() != o.get_params())
    {
        throw std::runtime_error{"During exec(): Observer was not initialized correctly with simulation parameters"};
    }

    int rec_freq_shift = 0;
    if(o.get_sel_type() == selection_type::sporadic &&
            o.get_e_change_f_type() == env_change_freq_type::regular)
    {
        rec_freq_shift += o.get_selection_duration();
    }

    namespace sw = stopwatch;
    sw::Stopwatch my_watch;

    while(s.get_time() !=  s.get_n_gen())
    {
        sim::tick(s);

        o.store_env_func(s);
        o.store_var_fit(s);
        o.store_avg_fit(s);

        if(o.get_record_freq_top_inds() != 0 &&
                (s.get_time() - rec_freq_shift) %  o.get_record_freq_top_inds() == 0)
        {
            o.store_top_n_inds(s);
            o.store_inputs_and_optimals(s);
        }
        if( o.get_record_freq_spectrum() != 0 &&
                (s.get_time() - rec_freq_shift) % o.get_record_freq_spectrum() == 0)
        {
            o.store_network_spectrum_n_best(s); //take out, unrealistic amount of data, convert to summary x weight
        }
        if(s.get_time() % 1000 == 0)
        {
            auto lap_ms = my_watch.lap<sw::ms>();
            std::cout << "Cycle " << s.get_time() << " --Lap time in ms: " << lap_ms << std::endl;
        }

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

///Gets the inputs  of the nth generation
template<class O>
const std::vector<std::vector<double>>& get_nth_gen_inputs(const O& o, int gen)
{
    auto generation = std::find_if(o.get_inputs_and_optimals().begin(),
                                          o.get_inputs_and_optimals().end(),
                                          [&gen] (const inputs_optimals& io) {return io.m_gen == gen;});
    if( generation != o.get_inputs_and_optimals().end())
    {
        return generation->m_inputs;
    }
    else
    {
        return std::vector<std::vector<double>>{};
    }
}

///Gets the optimals of the nth generation
template<class O>
const std::vector<double>& get_nth_gen_optimals(const O& o, int gen)
{

    auto generation = std::find_if(o.get_inputs_and_optimals().begin(),
                                          o.get_inputs_and_optimals().end(),
                                          [&gen] (const inputs_optimals& io) {return io.m_gen == gen;});

    if(generation != o.get_inputs_and_optimals().end())
    {
        return generation->m_optimals;
    }
    else
    {
        return std::vector<double>{};
    }
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
void test_observer();

#endif // OBSERVER_H
