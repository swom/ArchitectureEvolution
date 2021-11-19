#ifndef OBSERVER_H
#define OBSERVER_H
#include "simulation.h"

class observer
{
public:
    observer();

    NLOHMANN_DEFINE_TYPE_INTRUSIVE(observer,
                                   m_avg_fitnesses,
                                   m_var_fitnesses,
                                   m_top_inds,
                                   m_env_values,
                                   m_params,
                                   m_input,
                                   m_optimal)

    ///Saves the avg fitness and current environment value
    void store_avg_fit(const simulation& s);

    ///Saves the 100 best individuals in the population
    void save_best_n_inds(const simulation& s, int n);

    const all_params& get_params() const noexcept {return m_params;};

    void store_par (const simulation& s) {m_params = s.get_params();}

    void store_input(const simulation& s){m_input.push_back(s.get_input());}

    void store_optimal(const simulation& s){m_optimal.push_back(s.get_optimal());}

    const std::vector<std::vector<double>>& get_input() const noexcept {return m_input;}

    const std::vector<double>& get_optimal() const noexcept {return m_optimal;}

private:

    std::vector<double> m_avg_fitnesses;
    std::vector<double> m_var_fitnesses;
    std::vector<double> m_env_values;
    std::vector<std::vector<individual>> m_top_inds;
    all_params m_params = {};
    std::vector<std::vector<double>> m_input;
    std::vector<double> m_optimal;
};

bool operator==(const observer& lhs, const observer& rhs);

bool operator==(const all_params& lhs, const all_params& rhs);

bool operator!=(const all_params& lhs, const all_params& rhs);


///Executes a simulation for n generations
void exec(simulation& s , observer& o);

///Saves the enitre GODDDAM SIMULATIONNNN!!!!!!! WHOO NEEDS MEMORRYYYY
void save_json(const observer &o, const std::string& filename);

///Loads the observer back from json file.
observer load_observer_json(const std::string& filename);

void test_observer();

#endif // OBSERVER_H
