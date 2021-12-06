#ifndef OBSERVER_H
#define OBSERVER_H
#include "simulation.h"

template<mutation_type M = mutation_type::weights>
class observer
{
public:
    observer(int top_proportion = 1);

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
    const std::vector<std::vector<individual>>& get_top_inds() const noexcept{return m_top_inds;}

    ///Saves the avg fitness
    void store_avg_fit(const simulation<M>& s);

    ///Saves the variance of the fitness
    void store_var_fit(const simulation<M>& s);

    ///Saves the top_proportion nth best individuals in the population
    void store_top_n_inds(const simulation<M>& s);

    ///Saves the nth best individuals in the population
    void store_top_n_inds(const simulation<M>& s, int proportion);

    const all_params& get_params() const noexcept {return m_params;};

    void store_env_func (const simulation<M>& s) noexcept {m_env_functions.push_back(get_name_current_function(s));}

    void store_par (const simulation<M>& s) noexcept {m_params = s.get_params();}

    void store_input(const simulation<M>& s) noexcept {m_input.push_back(s.get_input());}

    void store_optimal(const simulation<M>& s) noexcept {m_optimal.push_back(s.get_optimal());}

    const std::vector<std::vector<double>>& get_input() const noexcept {return m_input;}

    const std::vector<double>& get_optimal() const noexcept {return m_optimal;}

private:

    std::vector<double> m_avg_fitnesses;
    std::vector<double> m_var_fitnesses;
    std::vector<char> m_env_functions;
    int m_top_proportion;
    std::vector<std::vector<individual>> m_top_inds;
    all_params m_params = {};
    std::vector<std::vector<double>> m_input;
    std::vector<double> m_optimal;
};

template<mutation_type M>
bool operator==(const observer<M>& lhs, const observer<M>& rhs);

bool operator==(const all_params& lhs, const all_params& rhs);

bool operator!=(const all_params& lhs, const all_params& rhs);


///Executes a simulation for n generations
template<mutation_type M>
void exec(simulation<M>& s , observer<M>& o);

///Saves the enitre GODDDAM SIMULATIONNNN!!!!!!! WHOO NEEDS MEMORRYYYY
template<mutation_type M>
void save_json(const observer<M> &o, const std::string& filename);

///Loads the observer back from json file.
template<mutation_type M>
observer<M> load_observer_json(const std::string& filename);

void test_observer();

#endif // OBSERVER_H
