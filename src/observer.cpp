#include "observer.h"
#include <fstream>
bool operator == (const Rn_Data& lhs, const Rn_Data& rhs)
{
    bool gen = lhs.generation == rhs.generation;
    bool rns = lhs.m_reac_norm == rhs.m_reac_norm;
    return gen && rns;
}

bool operator== (const sensibilities_to_mut& lhs, const sensibilities_to_mut& rhs)
{
    auto generation = lhs.m_generation == rhs.m_generation;
    auto sensibilities = lhs.m_sensibilities == rhs.m_sensibilities;

    return generation & sensibilities;
}

bool operator!= (const sensibilities_to_mut& lhs, const sensibilities_to_mut& rhs)
{
    return !(lhs == rhs);
}

bool operator== (const inputs_optimals& lhs, const inputs_optimals& rhs)
{
    bool inputs = lhs.m_inputs == rhs.m_inputs;
    bool optimals = lhs.m_optimals == rhs.m_optimals;
    bool gen = lhs.m_gen == rhs.m_gen;
    return inputs && optimals && gen;
}

bool operator!= (const inputs_optimals& lhs, const inputs_optimals& rhs)
{
    return !(lhs == rhs);
}

bool operator==(const all_params& lhs, const all_params& rhs)
{
    bool ind_par = lhs.i_p == rhs.i_p;
    bool env_pars = lhs.e_p == rhs.e_p;
    bool pop_pars = lhs.p_p == rhs.p_p;
    bool sim_pars = lhs.s_p == rhs.s_p;

    return ind_par && env_pars && pop_pars && sim_pars;
}

template<class Sim>
bool operator==(const observer<Sim>& lhs, const observer<Sim>& rhs)
{
    auto same_par = lhs.get_params() ==  rhs.get_params();
    auto same_inputs_and_optimals = lhs.get_inputs_and_optimals() == rhs.get_inputs_and_optimals();
    auto same_avg_fitness = lhs.get_avg_fitness() == rhs.get_avg_fitness();
    auto same_fit_var = lhs.get_var_fitness() == rhs.get_var_fitness();
    auto same_top_inds = lhs.get_top_inds() == rhs.get_top_inds();
    auto same_env_func = lhs.get_env_funcs() == rhs.get_env_funcs();

    return same_par &&
            same_inputs_and_optimals &&
            same_env_func &&
            same_avg_fitness &&
            same_fit_var &&
            same_top_inds;

}

template<class Sim>
bool operator!=(const observer<Sim>& lhs, const observer<Sim>& rhs)
{
    return !(lhs == rhs);
}

bool operator!=(const all_params& lhs, const all_params& rhs)
{
    return !(lhs == rhs);
}

bool all_inds_rns_are_equal_to(const std::vector<reac_norm>& all_inds_rns,
                               const reac_norm& rn)
{
    return std::all_of(all_inds_rns.begin(), all_inds_rns.end(),
                       [&](const reac_norm& rn_ind){return rn_ind == rn;});
}

bool ancestor_ID_is_last_recorded_ind_ID(const observer<>& o, const simulation<>& s)
{
    return sim::ancestor_IDs(sim::get_inds(s)) == sim::pop_IDs(o.get_last_recorded_pop());
}


std::string create_mut_spec_save_name(const all_params &p)
{
    std::string name = "MSP_" +
            create_save_name_from_params(p);

    return name;;
}

fit_and_phen_sens_t find_best_fit_phen_combination(const sensibilities_to_mut &record)
{
    auto max_fit_sens = std::max_element(record.m_sensibilities.begin(),
                                         record.m_sensibilities.end(),
                                         [](const fit_and_phen_sens_t& lhs, const fit_and_phen_sens_t& rhs)
    {return lhs.m_fitness_sens < rhs.m_fitness_sens;});

    return {max_fit_sens->m_fitness_sens, 0};
}

fit_and_phen_sens_t find_sensibilities_from_highest_ranking_ind(const sensibilities_to_mut& record)
{
    auto best = std::max_element(record.m_sensibilities.begin(), record.m_sensibilities.end(),
                                 [](const fit_and_phen_sens_t& lhs, const fit_and_phen_sens_t& rhs){return lhs.m_rank < rhs.m_rank;});
    return *best;
}

sensibilities_to_mut get_inds_sensibilities_of_first_record(const observer<>& o)
{
    return o.get_first_fit_phen_mut_sens();
}

///load an observer of correct type based on mutation_type
std::unique_ptr<base_observer> load_observer_json_mutation_type(const all_params &pars)
{
    auto m_t = pars.i_p.m_mutation_type;
    switch (m_t) {
    case mutation_type::NRaddition :
        return load_observer_selection_type<mutation_type::NRaddition>(pars);
        break;
    case mutation_type::NRduplication :
        return load_observer_selection_type<mutation_type::NRduplication>(pars);
        break;
    case mutation_type::activation :
        return load_observer_selection_type<mutation_type::activation>(pars);
        break;
    case mutation_type::addition :
        return load_observer_selection_type<mutation_type::addition>(pars);
        break;
    case mutation_type::duplication :
        return load_observer_selection_type<mutation_type::duplication>(pars);
        break;
    case mutation_type::weights :
        return load_observer_selection_type<mutation_type::weights>(pars);
        break;
    case mutation_type::weights_and_activation :
        return load_observer_selection_type<mutation_type::weights_and_activation>(pars);
        break;
    default:
        throw std::invalid_argument{"invalid mutation_type when loading observer"};

    }
}

///load an observer of correct type based on parameter
std::unique_ptr<base_observer> load_observer_json_of_correct_type(const all_params &pars)
{
    return load_observer_json_mutation_type(pars);
}

std::unique_ptr<base_observer> load_observer_json(const std::string &filename)
{

    std::ifstream f(filename);
    nlohmann::json json_in;
    f >> json_in;
    auto pars = json_in.get<all_params>();

    auto obs_ptr = load_observer_json_of_correct_type(pars);

    return obs_ptr;
}

observer<> load_default_observer_json(const std::string &filename)
{
    return  load_json<observer<>>(filename);
}

void calculate_mut_spec_from_obs(observer<>& o,
                                 int n_gens = 0)
{

    auto gens = extract_gens(o.get_top_inds(), n_gens);

    namespace sw = stopwatch;
    sw::Stopwatch my_watch;

    for (int i = 0 ; i < int(gens.size()); i++)
    {
        std::cout << "Generation: " << i << std::endl;
        auto spectrums = o.calculate_mut_spectrums_for_gen(gens[i]);
        o.add_spectrum(spectrums);
    }
    std::cout << "All spectrums calculated. Time: " << my_watch.lap<stopwatch::ms>() << " ms" << std::endl;
}

observer<> calculate_mut_spec_from_loaded_observer_data(const all_params& params)
{
    auto o = load_json<observer<>>(create_save_name_from_params(params));

    calculate_mut_spec_from_obs(o, gens_intervals);
    return o;
}

std::string create_save_name_from_params(const all_params& p)
{

    std::string name = "mut_" + convert_mut_type_to_string(p.i_p.m_mutation_type).substr(0, 3) +
            "_sel_" + convert_selection_type_to_string(p.s_p.sel_type).substr(0, 3) +
            "_sym_" + convert_change_symmetry_type_to_string(p.s_p.change_sym_type).substr(0, 3) +
            "_fr_" + convert_change_freq_type_to_string(p.s_p.change_freq_type).substr(0, 3) +
            "_ap_" + convert_adapt_periods_to_string(p.s_p.adaptation_per).substr(0, 3) +
            "_r_" + convert_response_type_to_string(p.i_p.net_par.resp_type).substr(0, 3) +
            "_e_" + convert_eval_type_to_string(p.s_p.evalu_type).substr(0,3) +
            "_arc_" + convert_arc_to_string(p.i_p.net_par.net_arc) +
            "_marc_" + convert_arc_to_string(p.i_p.net_par.max_arc) +
//            "_wr_" + std::to_string(p.p_p.mut_rate_weight).substr(0,3) +
//            "_ar_" + std::to_string(p.p_p.mut_rate_activation).substr(0, 3) +
            "_dup_" + std::to_string(p.p_p.mut_rate_duplication).substr(0, 3) +
            "_cA_" + std::to_string(p.s_p.change_freq_A).substr(0, 3) +
            "_cB_" + std::to_string(p.s_p.change_freq_B).substr(0, 3) +
            "_st_" + std::to_string(p.s_p.selection_strength).substr(0, 3) +
            "_sf_" + std::to_string(p.s_p.selection_freq).substr(0, 4) +
            "_sp_" + std::to_string(p.s_p.selection_duration_prop_to_freq).substr(0,3) +
            "_fA_" + p.e_p.name_func_A +
            "_g_" + std::to_string(p.s_p.n_generations) +
            "_p_" + std::to_string(p.p_p.number_of_inds) +
            "_s" + std::to_string(p.s_p.seed) + ".json";

    std::cout   << "response type: " <<convert_response_type_to_string(p.i_p.net_par.resp_type) << std::endl;
    std::cout << "save name: " << name << std::endl;

    return name;
}

simulation<> create_simple_simulation(int n_gen, int n_inds, bool all_different)
{
    all_params a_p;
    a_p.i_p.net_par.function = sigmoid;
    a_p.i_p.net_par.max_arc = {1,1};
    a_p.i_p.net_par.net_arc = {1,1};

    a_p.e_p.cue_range = {1,1};
    a_p.e_p.env_function_A = sigmoid_env;

    a_p.s_p.n_generations = n_gen;

    auto s  = simulation{a_p};
    if(n_inds == 1)
        return s;

    s.get_pop_non_const() = create_simple_pop(n_inds, all_different);
    return s;
}

double squared_distance_from_best_combination(const fit_and_phen_sens_t &ind, const fit_and_phen_sens_t &best)
{
    return (ind.m_fitness_sens - best.m_fitness_sens) * (ind.m_fitness_sens - best.m_fitness_sens) +
            (ind.m_phenotype_sens - best.m_phenotype_sens) * (ind.m_phenotype_sens - best.m_phenotype_sens);
}

void save_mut_spec_obs(const observer<> &o)
{
    auto filename = create_mut_spec_save_name(o.get_params());
    save_mut_spectrums(o, filename);
}


#ifndef NDEBUG
void test_observer()
{
#define FIX_ISSUE_47
#ifdef FIX_ISSUE_47
    ///An observer stores the sim_param of a simulation
    /// at initialization
    {
        observer o_default;
        //Give sim some non-default params

        env_param e_p{"2", "1"};
        all_params params = {e_p,{},{},{}};

        simulation s{params};
        assert(o_default.get_params() != params);

        observer o_init{obs_param{},params};
        assert(o_init.get_params() == params);
    }
#endif

#define FIX_ISSUE_81
#ifdef FIX_ISSUE_81
    ///The observer stores inputs and optimal values
    {
        obs_param o_p;
        o_p.m_top_proportion = 1;
        observer o(o_p, {});

        //Give sim some non-default inputs and optimal
        simulation s{};
        int n_repeats = 100;

        for(int i = 0; i != n_repeats; ++i){
            sim::tick(s);

            o.store_inputs_and_optimals(s);

            assert(get_nth_gen_inputs(o, s.get_time())  == s.get_stored_inputs());
            assert(get_nth_gen_optimals(o, s.get_time()) == s.get_stored_optimals());
        }

        auto name = "obs_save_test";
        save_json(o, name);
        auto loaded_o = load_json<observer<>>(name);
        assert(o == loaded_o);

        auto o1 = o;
        assert(o1.get_avg_fitness().empty());
        //store some inds in observer (not stored for now)
        o1.store_data_based_on_sensibilities(s);
        assert(o1 != o);

        //store some fitnesses in observer (not stored for now)
        auto o2 = o1;
        assert(o2.get_var_fitness().empty());
        o2.store_data_based_on_sensibilities(s);
        assert(o2 != o1);
    }
#endif

#define FIX_ISSUE_103
#ifdef FIX_ISSUE_103
    ///Observers are compared for
    /// vector of fitness
    /// vector of variances
    /// vector of individuals
    /// vector of input values
    /// vector of optimal values
    /// params
    {
        obs_param o_p;
        o_p.m_top_proportion = 1;
        observer o(o_p, {});
        simulation s{};

        sim::tick(s);

        auto o1 = o;
        assert(o1.get_avg_fitness().empty());
        o1.store_avg_fit(s);
        assert(o1 != o);

        auto o2 = o1;
        assert(o2.get_var_fitness().empty());
        o2.store_var_fit(s);
        assert(o2 != o1);

        auto o3 = o2;
        assert(o3.get_top_inds().empty());
        o3.store_data_based_on_sensibilities(s);
        assert(o3 != o2);

        auto o4 = o3;
        assert(o4.get_inputs_and_optimals().empty());
        o4.store_inputs_and_optimals(s);
        assert(o3 != o2);

        auto o5 = o4;
        assert(o5.get_env_funcs().empty());
        o5.store_env_func(s);
        assert(o5 != o4);
    }
#endif

    ///Observer can be saved and loaded with correct templates
    {
        using net_t = network<mutation_type::NRduplication>;
        using ind_t = individual<net_t>;
        using pop_t = population<ind_t>;
        using sim_t = simulation<pop_t,
        env_change_symmetry_type::asymmetrical,
        env_change_freq_type::regular,
        selection_type::sporadic>;

        observer<sim_t> o;

        save_json(o, "obs_test");
        auto loaded_o = load_default_observer_json("obs_test");

        assert( o.get_avg_fitness() == loaded_o.get_avg_fitness());
        assert( o.get_env_funcs() == loaded_o.get_env_funcs());
        assert( o.get_inputs_and_optimals() == loaded_o.get_inputs_and_optimals());
        assert( o.get_params() == loaded_o.get_params());
        assert( o.get_top_inds().size() == loaded_o.get_top_inds().size());
        for( size_t i = 0 ;  i != o.get_top_inds().size(); i++)
            for( size_t j = 0 ;  i != o.get_top_inds()[i].size(); i++)
            {
                assert( o.get_top_inds()[i][j].m_ind.get_net().get_net_weights() == loaded_o.get_top_inds()[i][j].m_ind.get_net().get_net_weights());
            }
        //assert(typeid (o) == typeid(loaded_o)); ideally this test should pass but cannot think of a way

    }

    ///It is possible to calculate the spectrum of individuals from a given generation
    /// and adds it to the vector of Ind_spectrum data in observer
    {
        //make mutation happen but they do not change anything,
        //It is just required that the ouput of the two operations are the same
        double mut_step = 0.00000000000001;
        int length_of_simulation = 3;
        int n_inds = 1;

        all_params a_p;
        a_p.s_p.n_generations = length_of_simulation;
        a_p.p_p.number_of_inds = n_inds;
        a_p.p_p.mut_step = mut_step;
        a_p.i_p.net_par.net_arc = {1,1};
        a_p.i_p.net_par.max_arc = {1,1};
        simulation s{a_p};

        int record_freq = 1;
        int top_inds_recorded = 1;
        int n_data_points = 1;
        int n_mutations = 1;
        observer o({top_inds_recorded,
                    record_freq,
                    record_freq,
                    n_data_points,
                    n_mutations},
                   s.get_params());

        exec(s, o);
        assert(o.get_top_inds().size() == size_t(length_of_simulation) &&
               o.get_top_spectrums().size() == size_t(length_of_simulation));

        for(int i = 0; i != s.get_n_gen(); i++)
        {
            auto mut_spectrum = o.calculate_mut_spectrums_for_gen(i);
            assert(mut_spectrum == o.get_top_spectrums_gen(i));
        }

        ///It is possible to load an observer and calculate the mutational spetrum of all the recorded individuals
        save_json<observer<>>(o, create_save_name_from_params(o.get_params()));
        auto loaded_o = calculate_mut_spec_from_loaded_observer_data(o.get_params());
        for(int i = 0; i != s.get_n_gen(); i++)
        {
            auto original = o.get_top_spectrums().at(i);
            auto calculated_from_upload = loaded_o.get_top_spectrums().at(i);
            ///the upload will also calculate spectrums of sampled individuals
            /// therefore we just check if the spectrum of the top individual
            /// calculate while the simulaiton was running is also present
            assert(!std::ranges::search(original,calculated_from_upload).empty());
        }
    }

    ///It is possible o exclusively save the mutational spectrums recorded by the observer
    {

        auto s = create_simple_simulation();
        observer o{{}, s.get_params()};
        exec(s,o);
        calculate_mut_spec_from_obs(o);

        save_mut_spec_obs(o);
        auto loaded_spectrums = load_mut_specs<individual<>>(create_mut_spec_save_name(o.get_params()));
        assert(o.get_top_spectrums() == loaded_spectrums);
    }

    ///It is possible to decide how many generations
    /// mutational spectrums will be produced
    {
        int n_gen = 4;
        int n_of_spectrums = n_gen/2;

        auto s = create_simple_simulation(n_gen);
        observer o{{}, s.get_params()};
        exec(s,o);

        calculate_mut_spec_from_obs(o, n_of_spectrums);
        assert(o.get_top_spectrums().size() == n_of_spectrums);
    }

    ///A simulation can be run so that nth (= adaptation_period_proportion) of the time
    /// is an 'adaptation period' where evolution happens in a stable environment
    {
        using sim_t = simulation<population<>,
        env_change_symmetry_type::asymmetrical,
        env_change_freq_type::regular,
        selection_type::constant,
        adaptation_period::on>;

        ///When change is allowed
        /// the enviromental function
        /// will alwys be b
        sim_param s_p;
        s_p.change_freq_A = 1;
        s_p.change_freq_B = 0;
        s_p.n_generations = s_p.adaptation_period_proportion;

        obs_param o_p;
        o_p.m_spectrum_reg_freq = 0;
        o_p.m_top_ind_reg_freq = 0;

        sim_t s{{{},{},{},s_p}};
        observer<sim_t> o{o_p, s.get_params()};

        exec(s,o);

        for(size_t i = 0; i != o.get_env_funcs().size(); i++ )
        {
            if(i < o.get_env_funcs().size() / s_p.adaptation_period_proportion)
            {
                assert(o.get_env_funcs()[i] == 'A');
                continue;
            }
            assert(o.get_env_funcs()[i] == 'B');
        }

    }

    ///If selection is sporadic and happens at regular intervals
    /// Then the individuals will be recorded
    /// following the indicated recording frequency
    /// but with an appropriate delay, so to capture individuals
    /// that have eundergone selection
    {
        using sim_t = simulation<population<>,
        env_change_symmetry_type::asymmetrical,
        env_change_freq_type::regular,
        selection_type::sporadic,
        adaptation_period::off>;

        int selection_freq = 2;

        sim_param s_p;
        s_p.sel_type = selection_type::sporadic;
        s_p.change_freq_type = env_change_freq_type::regular;
        s_p.n_generations = 10;
        s_p.selection_freq = selection_freq;
        s_p.selection_duration = selection_freq / 2;

        obs_param o_p;
        o_p.m_top_ind_reg_freq = selection_freq;


        sim_t s{{{},{},{},{s_p}}};
        observer<sim_t> o(o_p, s.get_params());

        exec(s,o);

        bool is_registered_at_the_end_of_sleection_period = false;

        ///BAD!!!
        for(const auto& top_ind : o.get_top_inds())
        {
            if(top_ind[0].generation % (s_p.selection_freq + o.get_selection_duration() - 1) == 0)
            {
                is_registered_at_the_end_of_sleection_period = true;
            }
        }
        assert(is_registered_at_the_end_of_sleection_period);

    }

    //An observer can store the inputs and the optimal values of a the last generation run in the simulation
    {
        simulation s;
        observer o;
        exec(s,o);
        o.store_inputs_and_optimals(s);
        assert(get_nth_gen_inputs(o,s.get_time()) == s.get_stored_inputs());
        assert(get_nth_gen_optimals(o,s.get_time()) == s.get_stored_optimals());
        assert(o.get_inputs_and_optimals().back().m_gen == s.get_time());
    }

    ///Inputs and optimals are stored only when individuals are stored
    {
        simulation s;
        observer o;
        exec(s,o);
        auto n_recorded_inds = o.get_top_inds().size();
        auto n_recorded_inputs_outputs = o.get_inputs_and_optimals().size();
        assert( n_recorded_inds == n_recorded_inputs_outputs);
    }

    ///Observers can record the avg robustness of the population
    {
        double robust_weight = 10;
        double frail_weight = 0;

        simulation robust_sim;
        robust_sim.changel_all_inds_weights(robust_weight);

        simulation frail_sim;
        assert(pop::all_inds_weights_have_value(frail_sim.get_pop_non_const(),
                                                frail_weight));

        observer o_frail;
        observer o_robust;

        assert(o_frail.get_avg_mutation_sensibility().empty());
        o_frail.store_avg_mut_sensibility(frail_sim);
        assert(o_frail.get_avg_mutation_sensibility().size() == 1);

        assert(o_robust.get_avg_mutation_sensibility().empty());
        o_robust.store_avg_mut_sensibility(robust_sim);
        assert(o_robust.get_avg_mutation_sensibility().size() == 1);

        assert(pairwise_comparison_for_majority(o_robust.get_avg_mutation_sensibility(),
                                                o_frail.get_avg_mutation_sensibility()));

    }

    //    ///The average sensibility to mutations is recorded every generation
    //    {
    //        observer o;
    //        simulation s;
    //        exec(s,o);
    //        assert(o.get_avg_mutation_sensibility().size() == s.get_time());
    //    }

    ///It is possible to calculate the avg robustness of a population in a simulation
    {
        double robust_weight = 10;
        double frail_weight = 0;
        int n_mutations = 10;

        simulation robust_sim;
        robust_sim.changel_all_inds_weights(robust_weight);

        simulation frail_sim;
        assert(pop::all_inds_weights_have_value(frail_sim.get_pop_non_const(),
                                                frail_weight));

        double robust = calc_avg_mutation_sensibility(robust_sim, n_mutations);
        double frail = calc_avg_mutation_sensibility(frail_sim, n_mutations);

        assert(robust > frail);
    }


    ///The avg population mutation sensibility is saved
    {
        sim_param s_p;
        s_p.n_generations = 1;
        simulation s{{env_param{}, ind_param{}, pop_param{}, s_p}};
        observer o({}, s.get_params());
        exec(s,o);
        save_json(o, "test_save");
        auto loaded_o = load_default_observer_json("test_save");
        assert(o.get_avg_mutation_sensibility() == loaded_o.get_avg_mutation_sensibility());
    }

    ///The fitness and phenotype mutation sensibility can be stored in observer
    /// #1 when sensibilities are stored they are stored in 2 separate vectors
    {
        simulation s;
        observer o;

        assert( o.get_fit_phen_mut_sensibility().empty());
        o.store_data_based_on_sensibilities(s);

        assert( !o.get_fit_phen_mut_sensibility().empty());
    }

    ///The fitness and phenotype mutation sensibility can be stored in observer
    /// #2 a value for each individual in the population is stored
    {

        all_params a_p;
        a_p.p_p.number_of_inds = 2;
        simulation s{a_p};
        observer o;

        o.store_data_based_on_sensibilities(s);

        for(const auto& sensibilities : o.get_fit_phen_mut_sensibility())
            assert(sensibilities.size() == s.get_inds().size());

    }

    ///The fitness and phenotype mutation sensibility can be stored in observer
    /// #3 a population of inds with hihger weights should have a hihger mutational robustness
    {
        double robust_weight = 10;
        double frail_weight = 0;

        simulation robust_sim = create_simple_simulation();
        robust_sim.changel_all_inds_weights(robust_weight);

        simulation frail_sim = create_simple_simulation();
        assert(pop::all_inds_weights_have_value(frail_sim.get_pop_non_const(),
                                                frail_weight));

        observer o_frail;
        observer o_robust;

        o_robust.store_data_based_on_sensibilities(robust_sim);
        o_frail.store_data_based_on_sensibilities(frail_sim);

        assert(lhs_has_lower_phen_mutation_sensibility_than_rhs(o_robust, o_frail));
    }

    ///The fitness and phenotype mutation sensibility can be stored in observer
    /// #4 a population of inds on the fitness peak should have a lower mutational robustness
    /// than a population not on the fitness peak if the weights are around the same megnitude
    {
        double optimal_weight = 1;
        double non_optimal_weight = -1;

        simulation optimal_sim = create_simple_simulation();
        optimal_sim.changel_all_inds_weights(optimal_weight);

        simulation non_optimal_sim = create_simple_simulation();
        non_optimal_sim.changel_all_inds_weights(non_optimal_weight);

        optimal_sim.calc_fitness();
        non_optimal_sim.calc_fitness();

        assert(sim::all_inds_have_fitness(1, optimal_sim));
        assert(!sim::all_inds_have_fitness(1, non_optimal_sim));

        observer o_non_optimal;
        observer o_optimal;

        o_optimal.store_data_based_on_sensibilities(optimal_sim);
        o_non_optimal.store_data_based_on_sensibilities(non_optimal_sim);

        assert(lhs_is_more_sensible_to_mutation_effects_on_fitness_than_rhs(o_optimal, o_non_optimal));
    }


    ///The fitness and phenotype mutation sensibility can be stored in observer
    /// #5 A population in which all inds are on the fitness peak
    /// have all individuals with a negative fitness sensibility
    {
        double optimal_weight = 1;
        auto optimal_sim = create_simple_simulation();
        optimal_sim.changel_all_inds_weights(optimal_weight);
        optimal_sim.calc_fitness();
        assert(sim::all_inds_have_fitness(1, optimal_sim));

        observer o;
        o.store_data_based_on_sensibilities(optimal_sim);

        assert(all_fit_mut_sens_are_negative(o));
    }

    ///The fitness and phenotype sensibilities to mutations of the population
    /// are saved every N generations
    /// toghether with the top individuals
    {
        auto s = create_simple_simulation();

        int rec_freq_of_top_inds = 1;
        obs_param o_p(1,rec_freq_of_top_inds);
        observer o(o_p, s.get_params());

        assert(o.get_top_inds().empty());
        assert(o.get_fit_phen_mut_sensibility().empty());

        exec(s,o);

        assert(!o.get_top_inds().empty());
        assert(!o.get_fit_phen_mut_sensibility().empty());

        assert(o.get_first_fit_phen_mut_sens().m_generation == o.get_first_top_ind_of_first_gen().generation);
    }

    ///The fitness and phenotype sensibilities to mutations
    /// can be saved by the observer
    {
        auto s = create_simple_simulation();

        int rec_freq_of_top_inds = 1;
        obs_param o_p(1,rec_freq_of_top_inds);
        observer o(o_p, s.get_params());

        exec(s,o);

        std::string filename{"test"};
        save_json(o, filename);
        auto loaded_o = load_default_observer_json(filename);

        assert(loaded_o.get_fit_phen_mut_sensibility() == o.get_fit_phen_mut_sensibility());
    }

    ///The best individuals also store their sensibilities to mutation
    /// #1 Ind_data contain sensibilities and it is possible to retrieve them
    ///    Ind_data when created also store a non default insantioation of fit_and_phen_sens_t
    {
        auto s = create_simple_simulation();
        observer o;
        sim::calc_fitness_of_pop(s);

        o.store_data_based_on_sensibilities(s);

        auto sensibilities_of_first_top_ind = get_first_top_ind_of_first_record(o).m_sensibilities;
        assert(sensibilities_of_first_top_ind != fit_and_phen_sens_t{});
    }
    ///The best individuals also store their sensibilities to mutation
    /// The sensibilities stored in ind_data correspond to the one actually belonging to the best individual
    {
        auto s = create_simple_simulation();
        observer o;
        s.calc_fitness();

        o.store_data_based_on_sensibilities(s);

        auto best_top_ind_of_first_record = get_first_top_ind_of_first_record(o);
        auto sensibilities_first_record = get_inds_sensibilities_of_first_record(o);
        auto sensibility_of_ind_with_highest_ranking_fitness = find_sensibilities_from_highest_ranking_ind(sensibilities_first_record);

        assert(best_top_ind_of_first_record.m_sensibilities == sensibility_of_ind_with_highest_ranking_fitness);
    }

    ///Sensibilities also store the fitness of the individual
    {
        auto s = create_simple_simulation();
        observer o;
        s.calc_fitness();

        o.store_data_based_on_sensibilities(s);
        auto sensibilities_first_record = get_inds_sensibilities_of_first_record(o);
        auto sensibility_of_ind_with_highest_ranking_fitness = find_sensibilities_from_highest_ranking_ind(sensibilities_first_record);

        assert(sensibility_of_ind_with_highest_ranking_fitness.m_fitness ==
               get_first_top_ind_of_first_record(o).m_ind.get_fitness());
    }
    ///The sensibilities are saved with the rest of the Ind_Data
    {
        fit_and_phen_sens_t non_default_sensibilities(78945,456987);
        Ind_Data<individual<>> i;
        i.m_sensibilities = non_default_sensibilities;

        save_json(i, "test");
        auto loaded_i = load_json<Ind_Data<individual<>>>("test");
        assert(loaded_i.m_sensibilities == i.m_sensibilities);
    }

    ///It is possible to record indivudals based on their sensibilities value
    /// #1 3 individuals are recorded each time
    {
        int n_inds = 3;
        auto s = create_simple_simulation();
        s.get_pop_non_const() = create_simple_pop(n_inds);
        observer o;

        o.store_data_based_on_sensibilities(s);
        o.store_top_mid_low_sens_inds(s);

        assert(o.get_sampled_inds().back().size() == 3);
    }

    /// #2 The individuals sampled are those with the highest, the lowest
    /// and the most median sensibilites
    {
        int n_inds = 3;
        auto s = create_simple_simulation();
        s.get_pop_non_const() = create_simple_pop(n_inds);
        observer o;

        s.sort_and_assign_ranks_by_fitness();
        o.store_data_based_on_sensibilities(s);

        auto top_mid_low_sens_inds = sample_top_mid_low_sens_inds(s, o);
        assert(top_mid_low_sens_inds.size() == 3);

        auto distance_top_ind = distance_from_best_sens_comb(top_mid_low_sens_inds[0],
                o.get_first_fit_phen_mut_sens());
        auto distance_mid_ind = distance_from_best_sens_comb(top_mid_low_sens_inds[1],
                o.get_first_fit_phen_mut_sens());
        auto distance_low_ind = distance_from_best_sens_comb(top_mid_low_sens_inds[2],
                o.get_first_fit_phen_mut_sens());

        assert(distance_top_ind <  distance_mid_ind);
        assert(distance_top_ind <  distance_low_ind);
        assert(distance_mid_ind <  distance_low_ind);
    }

    /// #2.1 it is possible to find the best fitness and phenotype sensibility combination
    /// looking at a vector of fitness and phen sensibilities
    /// the best combination is equal to the
    /// {max recorded fitness sensibility, 0}
    /// since 0 is the best possible phenotype sensibility
    {
        fit_and_phen_sens_t best{1,0};
        fit_and_phen_sens_t worst{-1,0};
        sensibilities_to_mut vec{1,{best,worst}};
        assert(find_best_fit_phen_combination(vec) == best);

        fit_and_phen_sens_t best_but_not_otpimal{1,1};
        sensibilities_to_mut  best_but_not_optimal_vec{1,{best_but_not_otpimal,worst}};
        assert(find_best_fit_phen_combination(best_but_not_optimal_vec) == best);
    }

    ///It is possible to calculate the distance of and individual sensibilites from
    /// the best possible combination of sensibilities
    /// this is the euclidean distance if  phen and fit sensibilities are the cartesian axes
    {
        fit_and_phen_sens_t best{1,0};
        fit_and_phen_sens_t mid{1,1};
        fit_and_phen_sens_t worst{0,1};

        sensibilities_to_mut vec{1,{best, worst}};

        auto best_comb = find_best_fit_phen_combination(vec);

        auto distance_best = squared_distance_from_best_combination(best, best_comb);

        assert(distance_best == 0);

        auto distance_mid = squared_distance_from_best_combination(mid, best_comb);

        assert(distance_mid == 1);

        auto distance_worst = squared_distance_from_best_combination(worst, best_comb);

        assert(distance_worst == 2);
    }
    /// #3 inds are sampled every 10 times the top inds are recorded by default,
    /// as determined by observer params
    {
        obs_param o_p{};
        int n_inds = 3;
        int n_gens = o_p.m_sample_ind_record_freq_to_sens;
        auto s = create_simple_simulation(n_gens);
        s.get_pop_non_const() = create_simple_pop(n_inds);
        observer o(o_p, s.get_params());
        exec(s,o);

        assert(o.get_fit_phen_mut_sensibility().size() ==
               o.get_sampled_inds().size() * o_p.m_sample_ind_record_freq_to_sens);
    }

    ///#4 Sampled individuals are saved
    {
        obs_param o_p{};
        int n_gens = o_p.m_sample_ind_record_freq_to_sens;
        auto s = create_simple_simulation(n_gens);
        s.get_pop_non_const() = create_simple_pop();
        observer o({}, s.get_params());
        exec(s,o);

        save_json(o,"test");
        observer load_o = load_default_observer_json("test");

        assert(o.get_sampled_inds() == load_o.get_sampled_inds());
    }

    ///The population that reproduced (stored in m_new_inds_vec after sim::tick(s))
    /// is sorted by fitness and individuals are assigned a rank
    /// when data about all inds in the population is saved
    {
        simulation s;
        observer o;
        s.get_pop_non_const() = create_simple_pop(5);

        sim::tick(s);
        assert(pop::all_fitnesses_are_not_equal(s.get_new_inds()));

        assert(!pop::is_sorted_by_fitness(s.get_new_inds()));
        assert(!pop::is_sorted_by_rank(s.get_new_inds()));

        assert(is_time_to_record_inds(o,s));
        o.store_ind_data(s);

        assert(pop::is_sorted_by_fitness(s.get_new_inds()));
        assert(pop::is_sorted_by_rank(s.get_new_inds()));
    }

    ///Individuals are assigned an ancestor ID equal to their parents fitness rank
    /// when it is the time to record them
    /// this will be a unique ID that toghether with the the generation
    /// in which the individuals is recorded will make it possible to identify their lineage
    /// #1
    {
        auto s = create_simple_simulation(10);
        assert(s.get_pop_size() ==  1);

        obs_param o_p;
        int frequency_of_recording = 1;
        o_p.m_top_ind_reg_freq = frequency_of_recording;
        observer o(o_p, s.get_params());

        exec(s,o);

        for(int gen = 1; gen < o.get_fit_phen_mut_sensibility().size(); gen++)
        {
            auto ancestor_inds = o.get_fit_phen_mut_sensibility()[gen - 1].m_sensibilities;
            auto inds = o.get_fit_phen_mut_sensibility()[gen].m_sensibilities;

            for(int ind = 0; ind < inds.size(); ind++)
            {
                assert(inds[ind].get_ancestor_ID() == ancestor_inds[ind].get_ID());
            }
        }
    }

    ///It is possible for the observer to record all the individuals' reaction norms in a population
    {
        int n_inds = 10;
        int n_gens = 1;
        auto s = create_simple_simulation(n_gens, n_inds, false);
        assert(sim::all_inds_have_same_net(s));
        observer o({}, s.get_params());

        auto range = o.get_params().e_p.cue_range;
        auto n_points = o.get_params().s_p.m_reac_norm_n_points;

        auto ind = s.get_inds().back();
        auto ind_rn = calculate_reaction_norm(ind.get_net(),range, n_points);

        o.store_all_inds_rn(s);
        auto last_recorded_all_inds_rn = o.get_all_inds_rn().back();

        assert(all_inds_rns_are_equal_to(last_recorded_all_inds_rn.m_reac_norm, ind_rn));
    }

    ///It is possible to record all the reaction norms of all individuals every n generations
    {
        int n_inds = 10;
        int n_gens = 6;
        auto s = create_simple_simulation(n_gens, n_inds, false);

        int record_all_inds_rn_every_n_gens = 2;
        obs_param o_p(1,1,0,1, record_all_inds_rn_every_n_gens);
        observer o(o_p, s.get_params());

        exec(s,o);

        int expected_number_of_records = n_gens / record_all_inds_rn_every_n_gens;
        assert(expected_number_of_records == o.get_all_inds_rn().size());
    }

    ///It is possible to save th reaction norms with the rest of the dataof the observer
    {
        int n_inds = 2;
        int n_gens = 4;
        auto s = create_simple_simulation(n_gens, n_inds, false);

        int record_all_inds_rn_every_n_gens = 2;
        obs_param o_p(1,1,0,1, record_all_inds_rn_every_n_gens);
        observer o(o_p, s.get_params());

        exec(s,o);

        std::string filename{"test_all_inds_rn"};
        save_json(o, filename);
        auto loaded_o =load_json<observer<>>(filename);

        assert(loaded_o.get_all_inds_rn() == o.get_all_inds_rn());
    }

    //it is possible to extract the top and sampled individuals for a generation all toghether
    {
        int sampling_freq = 1;
        obs_param o_p(1,sampling_freq,0,1,0,sampling_freq);
        observer o{o_p};
        simulation s;

        exec(s,o);
        for (int i = 0; i != s.get_n_gen(); i++)
        {
            auto inds = o.get_top_and_sampled_inds_for_generation(i);
        assert(!std::ranges::search(inds, o.get_top_inds_gen(i)).empty() &&
               !std::ranges::search(inds, o.get_sampled_inds_gen(i)).empty());
        }

    }

}
#endif


