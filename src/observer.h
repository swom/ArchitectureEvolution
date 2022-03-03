#ifndef OBSERVER_H
#define OBSERVER_H
#include "simulation.h"
#include "Stopwatch.hpp"


template<class Sim = simulation<>>
class observer
{
    using Pop = typename Sim::pop_t;
    using Ind = typename Pop::ind_t;
    using Net = typename Ind::net_t;
    using Net_Spect = network_spectrum<Net>;

private:

    std::vector<double> m_avg_fitnesses;
    std::vector<double> m_var_fitnesses;
    std::vector<char> m_env_functions;
    int m_top_proportion;
    std::vector<std::vector<Ind_Data<Ind>>> m_top_inds;
    std::vector<std::vector<Net_Spect>> m_top_spectrums;
    all_params m_params = {};
    std::vector<std::vector<double>> m_input;
    std::vector<double> m_optimal;

public:
    observer(int top_proportion = 1):
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
    const std::vector<std::vector<Ind_Data<Ind>>>& get_top_inds() const noexcept{return m_top_inds;}

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
        m_top_inds.push_back(calculate_reaction_norms(sim::get_best_n_inds(s, m_top_proportion),
                                                      s.get_env_cue_range(),
                                                      100,
                                                      s.get_time()
                                                      )
                             );
    }

    ///Saves the nth best individuals in the population
    void store_top_n_inds(const Sim& s, int proportion)
    {
        m_top_inds.push_back(calculate_reaction_norms(
                                 sim::get_best_n_inds(s, proportion),
                                 s.get_env_cue_range(),
                                 100,
                                 s.get_time()
                                 )
                             );
    }

    ///Stores the network spectrum of the top n best individuals
    void store_network_spectrum_n_best(Sim& s)
    {

        m_top_spectrums.emplace_back(std::vector<Net_Spect>{});
        std::vector<Net_Spect> spectrums = calculate_mut_spectrums(sim::get_best_n_inds(s, m_top_proportion),
                                                              s.get_mut_step(),
                                                              s.get_rng(),
                                                              1000,
                                                              s.get_env_cue_range(),
                                                              100);
         spectrums.swap(m_top_spectrums.back());
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

    while(s.get_time() !=  s.get_n_gen())
    {
        sim::tick(s);

        o.store_env_func(s);
        o.store_var_fit(s);
        o.store_input(s);
        o.store_optimal(s);
        o.store_avg_fit(s);


        if(s.get_time() % 1000 == 0)
        {
            o.store_top_n_inds(s);
        }
        if(s.get_time() % 1 == 0)
        {
            o.store_network_spectrum_n_best(s);
        }
        if(s.get_time() % 1 == 0)
        {
            auto lap_ms = my_watch.lap<sw::ms>();
            std::cout << "Cycle " << s.get_time() << " --Lap time in ms: " << lap_ms << std::endl;
        }

    }
}

template<class Sim>
std::string create_save_name_from_observer_data(const observer<Sim>& o)
{

    return "mut_type_" + convert_mut_type_to_string(o.get_params().i_p.m_mutation_type) +
            "_start_arc" + convert_arc_to_string(o.get_params().i_p.net_par.net_arc) +
            "_act_r" + std::to_string(o.get_params().p_p.mut_rate_activation).substr(0, 8) +
            "_dup_r" + std::to_string(o.get_params().p_p.mut_rate_duplication).substr(0, 8) +
            "_ch_A" + std::to_string(o.get_params().s_p.change_freq_A).substr(0, 8) +
            "_ch_B" + std::to_string(o.get_params().s_p.change_freq_B).substr(0, 8) +
            "_ch_type" + convert_change_type_to_string(o.get_params().s_p.change_type) +
            "_sel_str" + std::to_string(o.get_params().s_p.selection_strength).substr(0, 3) +
            "_max_arc" + convert_arc_to_string(o.get_params().i_p.net_par.max_arc) +
            "_sel_type" + convert_selection_type_to_string(o.get_params().s_p.sel_type) +
            "_" + std::to_string(o.get_params().s_p.seed) + ".json";

}

void test_observer();

#endif // OBSERVER_H