#ifndef OBSERVER_H
#define OBSERVER_H
#include "simulation.h"
#include "Stopwatch.hpp"


template<class Sim = simulation<>>
class observer
{
private:

    std::vector<double> m_avg_fitnesses;
    std::vector<double> m_var_fitnesses;
    std::vector<char> m_env_functions;
    int m_top_proportion;
    using Pop = typename Sim::pop_t;
    using Ind = typename Pop::ind_t;
    std::vector<std::vector<Ind>> m_top_inds;
    all_params m_params = {};
    std::vector<std::vector<double>> m_input;
    std::vector<double> m_optimal;

public:
    observer(int top_proportion = 10):
        m_top_proportion{top_proportion}
    {
    }


    NLOHMANN_DEFINE_TYPE_INTRUSIVE(observer,
                                   m_avg_fitnesses,
                                   m_var_fitnesses,
                                   m_top_inds,
                                   m_env_functions,
                                   m_params,
                                   m_input,
                                   m_optimal,
                                   m_top_proportion)

    ///returns const ref to m_avg_fitness
    const std::vector<double>& get_avg_fitness() const noexcept{return m_avg_fitnesses;}

    ///returns const ref to vector of env_functions' names
    const std::vector<char>& get_env_funcs() const noexcept {return m_env_functions;}

    ///returns const ref to m_var_fitnesses
    const std::vector<double>& get_var_fitness() const noexcept{return m_var_fitnesses;}

    ///returns const ref to best_ind vector
    const std::vector<std::vector<Ind>>& get_top_inds() const noexcept{return m_top_inds;}

    ///Saves the avg fitness
    void store_avg_fit(const Sim &s)
    {
        m_avg_fitnesses.push_back(sim::avg_fitness(s));
    }

    ///Saves the variance of the fitness
    void store_var_fit(const Sim& s)
    {
        m_var_fitnesses.push_back(sim::var_fitness(s));
    }

    ///Saves the top_proportion nth best individuals in the population
    void store_top_n_inds(const Sim& s)
    {
        m_top_inds.push_back(sim::get_best_n_inds(s, m_top_proportion));
    }

    ///Saves the nth best individuals in the population
    void store_top_n_inds(const Sim& s, int proportion)
    {
        m_top_inds.push_back(sim::get_best_n_inds(s, proportion));
    }

    const all_params& get_params() const noexcept {return m_params;};

    void store_env_func (const Sim& s) noexcept {m_env_functions.push_back(sim::get_name_current_function(s));}

    void store_par (const Sim& s) noexcept {m_params = s.get_params();}

    void store_input(const Sim& s) noexcept {m_input.push_back(s.get_input());}

    void store_optimal(const Sim& s) noexcept {m_optimal.push_back(s.get_optimal());}

    const std::vector<std::vector<double>>& get_input() const noexcept {return m_input;}

    const std::vector<double>& get_optimal() const noexcept {return m_optimal;}


};

template<class Ind>
bool operator==(const observer<Ind>& lhs, const observer<Ind>& rhs);

bool operator==(const all_params& lhs, const all_params& rhs);

bool operator!=(const all_params& lhs, const all_params& rhs);


///Executes a simulation for n generations
template<class Sim>
void exec(Sim& s , observer<Sim>& o)
{
    o.store_par(s);

    namespace sw = stopwatch;
    sw::Stopwatch my_watch;

    for (int i = 0; i < s.get_n_gen(); i++)
    {
        sim::tick (s);

        o.store_avg_fit(s);
        o.store_env_func(s);
        o.store_var_fit(s);
        o.store_input(s);
        o.store_optimal(s);

        if(i % 1000 == 0)
        {
            o.store_top_n_inds(s);
        }
        if(i % 1000 == 0)
        {
            auto lap_ms = my_watch.lap<sw::ms>();
            std::cout << "Cycle " << i << " --Lap time in ms: " << lap_ms << std::endl;;
        }
    }
}

void test_observer();

#endif // OBSERVER_H
