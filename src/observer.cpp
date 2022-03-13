#include "observer.h"
#include <fstream>

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
    auto same_env_inputs = lhs.get_input() == rhs.get_input();
    auto same_env_optimum =lhs.get_optimal() == rhs.get_optimal();
    auto same_avg_fitness = lhs.get_avg_fitness() == rhs.get_avg_fitness();
    auto same_fit_var = lhs.get_var_fitness() == rhs.get_var_fitness();
    auto same_top_inds = lhs.get_top_inds() == rhs.get_top_inds();
    auto same_env_func = lhs.get_env_funcs() == rhs.get_env_funcs();

    return same_par &&
            same_env_inputs &&
            same_env_optimum &&
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


///load an observer of correct type based on mutation_type
std::unique_ptr<base_observer> load_observer_json_mutation_type(const all_params &pars)
{
    auto m_t = pars.i_p.m_mutation_type;
    switch (m_t) {
    case mutation_type::NRaddition :
        return load_observer_selection_type<mutation_type::NRaddition>(pars);
        break;
    case mutation_type::NRduplication :
        return load_observer_selection_type<mutation_type::NRaddition>(pars);
        break;
    case mutation_type::activation :
        return load_observer_selection_type<mutation_type::NRaddition>(pars);
        break;
    case mutation_type::addition :
        return load_observer_selection_type<mutation_type::NRaddition>(pars);
        break;
    case mutation_type::duplication :
        return load_observer_selection_type<mutation_type::NRaddition>(pars);
        break;
    case mutation_type::weights :
        return load_observer_selection_type<mutation_type::NRaddition>(pars);
        break;
    case mutation_type::weights_and_activation :
        return load_observer_selection_type<mutation_type::NRaddition>(pars);
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

observer<> calculate_mut_spec_from_observer_data(const all_params& params)
{
    auto o = load_json<observer<>>(create_save_name_from_params(params));
    auto gens = extract_gens(o.get_top_inds());


#ifndef NDEBUG
#pragma omp parallel for ordered
#else
#pragma omp parallel for
#endif
    for (int i = 0 ; i < gens.size(); i++)
    {
        auto spectrum = o.calculate_mut_spectrums_for_gen(gens[i]);
#pragma omp critical
        {
            o.add_spectrum(spectrum);
        }
    }
    return o;
}

std::string create_save_name_from_params(const all_params& p)
{

    return "mut_type_" + convert_mut_type_to_string(p.i_p.m_mutation_type) +
            "_start_arc" + convert_arc_to_string(p.i_p.net_par.net_arc) +
            "_act_r" + std::to_string(p.p_p.mut_rate_activation).substr(0, 8) +
            "_dup_r" + std::to_string(p.p_p.mut_rate_duplication).substr(0, 8) +
            "_ch_A" + std::to_string(p.s_p.change_freq_A).substr(0, 8) +
            "_ch_B" + std::to_string(p.s_p.change_freq_B).substr(0, 8) +
            "_ch_type" + convert_change_symmetry_type_to_string(p.s_p.change_sym_type) +
            "_ch_type" + convert_change_freq_type_to_string(p.s_p.change_freq_type) +
            "_sel_str" + std::to_string(p.s_p.selection_strength).substr(0, 3) +
            "_max_arc" + convert_arc_to_string(p.i_p.net_par.max_arc) +
            "_sel_type" + convert_selection_type_to_string(p.s_p.sel_type) +
            "_" + std::to_string(p.s_p.seed) + ".json";

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

        env_param e_p{env_func_2, env_func_1};
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
        observer o;

        //Give sim some non-default inputs and optimal
        simulation s{};
        int n_repeats = 100;

        for(int i = 0; i != n_repeats; ++i){
            sim::tick(s);

            if(!o.get_input().empty() && !o.get_optimal().empty())
            {
                assert(o.get_input().back() != s.get_input());
                assert(o.get_optimal().back() != s.get_optimal());
            }

            o.store_input(s);
            o.store_optimal(s);

            assert(o.get_input().back() == s.get_input());
            assert(o.get_optimal().back() == s.get_optimal());
        }

        auto name = "obs_save_test";
        save_json(o, name);
        auto loaded_o = load_json<observer<>>(name);
        assert(o == loaded_o);

        auto o1 = o;
        assert(o1.get_avg_fitness().empty());
        //store some inds in observer (not stored for now)
        o1.store_top_n_inds(s,1);
        assert(o1 != o);

        //store some fitnesses in observer (not stored for now)
        auto o2 = o1;
        assert(o2.get_var_fitness().empty());
        o2.store_top_n_inds(s,1);
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
        observer o;

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
        o3.store_top_n_inds(s, 1);
        assert(o3 != o2);

        auto o4 = o3;
        assert(o4.get_input().empty());
        o4.store_input(s);
        assert(o3 != o2);

        auto o5 = o4;
        assert(o5.get_optimal().empty());
        o5.store_optimal(s);
        assert(o5 != o4);

        auto o6 = o5;
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
        assert( o.get_input() == loaded_o.get_input());
        assert( o.get_optimal() == loaded_o.get_optimal());
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
        //make mutation hapen but they do not change anything,
        //It is just required that the ouput of the two operations are the same
        double mut_step = 0.00000000000001;
        int length_of_simulation = 3;
        int n_inds = 1;

        all_params a_p;
        a_p.s_p.n_generations = length_of_simulation;
        a_p.p_p.number_of_inds = n_inds;
        a_p.p_p.mut_step = mut_step;
        simulation s{a_p};

        int record_freq = 1;
        int top_inds_recorded = 1;
        int n_data_points = 4;
        int n_mutations = 10;
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
        auto loaded_o = calculate_mut_spec_from_observer_data(o.get_params());
        for(int i = 0; i != s.get_n_gen(); i++)
        {
            auto original = o.get_top_spectrums().at(i);
            auto calculated_from_upload = loaded_o.get_top_spectrums().at(i + length_of_simulation);
            assert(original == calculated_from_upload);
        }
    }

}
#endif


