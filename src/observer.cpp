#include "observer.h"
#include <fstream>


bool operator==(const all_params& lhs, const all_params& rhs)
{
    return lhs.i_p.net_par.net_arc ==  rhs.i_p.net_par.net_arc &&
            lhs.p_p.mut_rate_weight == rhs.p_p.mut_rate_weight &&
            lhs.p_p.mut_step == rhs.p_p.mut_step &&
            lhs.p_p.number_of_inds == rhs.p_p.number_of_inds &&
            lhs.s_p.change_freq_A == rhs.s_p.change_freq_A &&
            lhs.s_p.change_freq_B == rhs.s_p.change_freq_B &&
            lhs.s_p.seed == rhs.s_p.seed &&
            lhs.s_p.selection_strength == rhs.s_p.selection_strength &&
            are_same_env_functions(lhs.e_p.env_function_A, rhs.e_p.env_function_A)&&
            are_same_env_functions(lhs.e_p.env_function_B, rhs.e_p.env_function_B);
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

#ifndef NDEBUG
void test_observer()
{
#define FIX_ISSUE_47
#ifdef FIX_ISSUE_47
    ///An observer can store the sim_param of a simulation
    {
        observer o;
        //Give sim some non-default params

        env_param e_p{env_func_2, env_func_1};
        all_params params = {e_p,{},{},{}};

        simulation s{params};
        assert(o.get_params() != params);

        o.store_par(s);

        assert(o.get_params() == params);
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


}
#endif