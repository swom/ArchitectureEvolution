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
                                   m_params)

    ///Saves the avg fitness and current environment value
    void store_avg_fit_and_env(const simulation& s);

    ///Saves the 100 best individuals in the population
    void save_best_n_inds(const simulation& s, int n);

    const all_params& get_params() const noexcept {return m_params;};

    void store_par (const simulation& s) {m_params = s.get_params();}

private:

    std::vector<double> m_avg_fitnesses;
    std::vector<double> m_var_fitnesses;
    std::vector<std::vector<individual>> m_top_inds;
    std::vector<double> m_env_values;
    all_params m_params = {};
};

bool operator==(const all_params& lhs, const all_params& rhs);

bool operator!=(const all_params& lhs, const all_params& rhs);


///Executes a simulation for n generations
void exec(simulation& s , observer& o);

///Saves the enitre GODDDAM SIMULATIONNNN!!!!!!! WHOO NEEDS MEMORRYYYY
void save_json(const observer &o, const std::string& filename);

void test_observer();

#endif // OBSERVER_H
