#ifndef ENVIRONMENT_H
#define ENVIRONMENT_H

#include <vector>
#include <random>
#include "json.hpp"

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

    const std::function<double(std::vector<double>)> &get_env_function_A() const {return m_env_function_A;}


private:

    ///The target value of the environment
    std::vector<double> m_ref_target_values;

    double m_current_target_value;

    /// A distribution to be used for determining cues
    std::uniform_real_distribution<double> m_cue_distribution;

    ///Points to The first function linking input to optimal output
    std::function<double(std::vector<double>)> m_env_function_A;

    ///The actual function
    //double env_function_A(const std::vector<double> &input);


};

///checks if 2 environments are equal
bool operator== (const environment& lhs, const environment& rhs);

void switch_target (environment &e);

void test_environment() noexcept;

///Create a vector of a given number of inputs with value fixed to 1
std::vector<double> create_n_inputs(int n_inputs);

///Create a vector of a given number of inputs from the distribution member of the environment
std::vector<double> create_n_inputs(environment e, const int &n_inputs, std::mt19937_64 &rng);


#endif // ENVIRONMENT_H

