#ifndef ENVIRONMENT_H
#define ENVIRONMENT_H

#include <vector>
#include <random>
#include "json.hpp"
#include "utilities.h"

static double env_func_A(std::vector<double> input){
  return input[0];
 }


static std::map<std::string, std::function<double(std::vector<double>)>> string_env_function_A_map
{
{"A", env_func_A}
};



struct env_param
{
    NLOHMANN_DEFINE_TYPE_INTRUSIVE(env_param,
                                   targetA,
                                   targetB)
double targetA;
double targetB;
std::function<double(std::vector<double>)> env_function_A{env_func_A};
};


class environment
{
public:
    ///deprecated(sort of)
    environment(double target_valueA, double target_valueB,
                std::function<double(std::vector<double>)> env_functionA = &env_func_A);

    environment(env_param e_p);

    std::uniform_real_distribution<double> get_dist() {return m_cue_distribution;}

    NLOHMANN_DEFINE_TYPE_INTRUSIVE(environment,
                                   m_ref_target_values,
                                   m_current_target_value);
    ///Returns the target value of the environment
    double get_current_target_value() const noexcept {return m_current_target_value;}
    const std::vector<double>& get_ref_target_values() const noexcept {return m_ref_target_values;}

    ///Sets current target value
    void set_current_target_value(double target_value) {m_current_target_value = target_value;}

    ///Returns the current inputs
    const std::vector<double> &get_input() const noexcept {return m_input;}

    ///Returns the current optimal output
    const double &get_optimal() const noexcept {return m_optimal_output;}

    ///Returns the cue distribution of the environment
    const std::uniform_real_distribution<double>&  get_cue_distribtion() const noexcept
        {return m_cue_distribution;}

    ///Updates the n first inputs by drawing random ones from the distribution
    std::vector<double> update_n_inputs(std::mt19937_64 &rng, const size_t n);

    const std::function<double(std::vector<double>)> &get_env_function_A() const {return m_env_function_A;}

    void change_uniform_dist(std::uniform_real_distribution<double> new_dist) {m_cue_distribution = new_dist;}

    double calculate_optimal();

    void update_optimal();


private:

    ///The target value of the environment
    std::vector<double> m_ref_target_values;

    double m_current_target_value;

    /// A distribution to be used for determining cues
    std::uniform_real_distribution<double> m_cue_distribution;


    ///The current inputs that the networks of individuals will recieve
    std::vector<double> m_input;

    ///The optimal output at a given moment; depends on inputs and environmental function
    double m_optimal_output;

    ///Points to The first function linking input to optimal output
    std::function<double(std::vector<double>)> m_env_function_A;

};

///checks if 2 environments are equal
bool operator== (const environment& lhs, const environment& rhs);

void switch_target (environment &e);

void test_environment() noexcept;

///Create a vector of a given number of inputs with value fixed to 1
std::vector<double> create_n_inputs(int n_inputs);

///Create a vector of a given number of inputs from a distribution
std::vector<double> create_n_inputs(std::uniform_real_distribution<double> dist, const int &n_inputs, std::mt19937_64 &rng);

///Create a vector of a given number of inputs from the distribution member of the environment
std::vector<double> create_n_inputs(environment e, const int &n_inputs, std::mt19937_64 &rng);

///Creates a vector of a given number of inputs for a distribution
std::vector<double> create_n_inputs(std::uniform_real_distribution<double> dist,
                                    const int &n_inputs, std::mt19937_64 &rng);


#endif // ENVIRONMENT_H

