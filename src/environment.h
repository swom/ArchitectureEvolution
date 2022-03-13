#ifndef ENVIRONMENT_H
#define ENVIRONMENT_H

#include <vector>
#include <random>
#include "range.h"
#include "utilities.h"

static double env_func_1(std::vector<double> input){
  return input[0] * input[0];
}

static double env_func_2(std::vector<double> input){
  return input[0]  * input[0] * input[0];
}




static std::map<std::string, std::function<double(std::vector<double>)>> string_env_function_map
{
{"1", env_func_1},
{"2", env_func_2}
};



struct env_param
{
    env_param(std::function<double(std::vector<double>)> fun_A = env_func_1,
              std::function<double(std::vector<double>)> fun_B = env_func_2,
              std::vector<double> env_cue_range = std::vector<double>{-1,1}) :
        env_function_A{fun_A},
        env_function_B{fun_B},
        cue_range{env_cue_range.front(),env_cue_range.back()}
    {}
std::function<double(std::vector<double>)> env_function_A;
std::function<double(std::vector<double>)> env_function_B;
range cue_range;
};

bool operator==(const env_param& lhs, const env_param& rhs);

class environment
{
public:
    ///deprecated(sort of)
    environment(std::function<double(std::vector<double>)> env_functionA = &env_func_1,
                std::function<double(std::vector<double>)> env_functionB = &env_func_2);

    environment(const env_param &e_p);

    std::uniform_real_distribution<double>& get_dist() noexcept {return m_cue_distribution;}

    ///Returns the cue distribution of the environment
    const std::uniform_real_distribution<double>&  get_cue_distribtion() const noexcept
    {return m_cue_distribution;}

    ///Returns the cue range
    const range& get_cue_range() const noexcept {return m_cue_range;}

    ///Updates the n first inputs by drawing random ones from the distribution
    std::vector<double> update_n_inputs(std::mt19937_64 &rng, const size_t n);

    ///Returns the environmental function A
    const std::function<double(std::vector<double>)> &get_env_function_A() const {return m_env_function_A;}

    ///Returns the environmental function B
    const std::function<double(std::vector<double>)> &get_env_function_B() const {return m_env_function_B;}

    ///Returns the current environmental function
    const std::function<double(std::vector<double>)> &get_current_function() const {return m_current_function;}

    ///Changes the cue distribution to a new given uniform distribution
    void change_uniform_dist(std::uniform_real_distribution<double> new_dist) {m_cue_distribution = new_dist;}

    ///Changes the current environmental function to a new given one.
    void change_env_function(std::function<double(std::vector<double>)> new_func) {m_current_function = new_func;}

    /// Returns the name of the current function
    const char &get_name_current_function() const {return m_name_current_function;}

    ///Switches the name of the current function from A to B or B to A
    void switch_name_current_function();

    ///Returns the environment's rng
    std::mt19937_64 &get_rng() {return m_rng;}



private:

    ///The range of the cue distribution
    range m_cue_range;

    /// A distribution to be used for determining cues
    std::uniform_real_distribution<double> m_cue_distribution;

    ///Points to The first function linking input to optimal output
    std::function<double(std::vector<double>)> m_env_function_A;

    ///Points to The second function linking input to optimal output
    std::function<double(std::vector<double>)> m_env_function_B;

    ///The environmental function that is currently used to determine the optimal output
    std::function<double(std::vector<double>)> m_current_function;

    /// The name of the current function
    char m_name_current_function;

    /// The environment's rng
    std::mt19937_64 m_rng;

};

///checks if 2 environments are equal
bool operator== (const environment& lhs, const environment& rhs);

namespace env {

///Create a vector of a given number of inputs with value fixed to 1
std::vector<double> create_n_inputs(int n_inputs);

///Create a vector of a given number of inputs from a distribution
std::vector<double> create_n_inputs(std::uniform_real_distribution<double> dist, const int &n_inputs, std::mt19937_64 &rng);

///Create a vector of a given number of inputs from the distribution member of the environment
std::vector<double> create_n_inputs(environment& e, const size_t &n_inputs, std::mt19937_64 &rng);


///Calculates the optimal output, given input, using the env function
double calculate_optimal(const environment &e, std::vector<double> input);

///Switches the current environmental function from A to B or B to A
void switch_env_function(environment &e);

}


///homebrew version of json loading and saving to load and save cue distributions for environment
void to_json(nlohmann::json& j, const environment& e);
void from_json(const nlohmann::json& j, environment& e);

void test_environment() noexcept;

#endif // ENVIRONMENT_H

