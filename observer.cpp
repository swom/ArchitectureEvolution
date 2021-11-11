#include "observer.h"
#include <fstream>
#include "Stopwatch.hpp"

observer::observer()
{
}


bool operator==(const all_params& lhs, const all_params& rhs)
{
    return lhs.i_p.net_par.net_arc ==  rhs.i_p.net_par.net_arc &&
           lhs.p_p.mut_rate == rhs.p_p.mut_rate &&
           lhs.p_p.mut_step == rhs.p_p.mut_step &&
           lhs.p_p.number_of_inds == rhs.p_p.number_of_inds &&
           lhs.s_p.change_freq == rhs.s_p.change_freq &&
           lhs.s_p.seed == rhs.s_p.seed &&
           lhs.s_p.selection_strength == rhs.s_p.selection_strength &&
           are_same_env_functions(lhs.e_p.env_function_A, rhs.e_p.env_function_A)&&
           are_same_env_functions(lhs.e_p.env_function_B, rhs.e_p.env_function_B);

}

bool operator!=(const all_params& lhs, const all_params& rhs)
{
    return !(lhs == rhs);
}



void observer::store_avg_fit(const simulation& s)
{
    m_avg_fitnesses.push_back(avg_fitness(s));
    m_var_fitnesses.push_back(var_fitness(s));
}

void observer::save_best_n_inds(const simulation &s, int n)
{
    m_top_inds.push_back(get_best_n_inds(s, n));
}

void save_json(const observer& o, const std::string& filename)
{
    std::ofstream  f(filename);
    nlohmann::json json_out;
    json_out = o;
    f << json_out;
}

void exec(simulation& s , observer& o)
{
    stopwatch::Stopwatch sw;
    o.store_par(s);
    for (int i = 0; i < s.get_n_gen(); i++)
    {
        tick (s);
        o.store_avg_fit(s);
        if(i % 1000 == 0)
        {
            o.save_best_n_inds(s,10);
        }
        if(i % 1000 == 0)
        {
            std::cout << "Cycle " << i << ". Elapsed: " << sw.lap<stopwatch::s>() << " seconds." << std::endl;
        }
    }
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
}
#endif
