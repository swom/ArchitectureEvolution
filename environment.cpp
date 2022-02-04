#include "environment.h"
#include "utilities.h"
#include <cassert>
#include <iostream>


environment::environment(std::function<double(std::vector<double>)> env_functionA,
                         std::function<double(std::vector<double>)> env_functionB):
    m_cue_range{-1,1},
    m_env_function_A{env_functionA},
    m_env_function_B{env_functionB},
    m_current_function{env_functionA},
    m_name_current_function{'A'}
{
    m_rng.seed(0);




}

environment::environment(const env_param& e_p):
    m_cue_range{e_p.cue_distrib},
    m_cue_distribution{m_cue_range.m_start, m_cue_range.m_end},
    m_env_function_A{e_p.env_function_A},
    m_env_function_B{e_p.env_function_B},
    m_current_function{e_p.env_function_A},
    m_name_current_function{'A'}
{
    m_rng.seed(0);
}



bool operator== (const environment& lhs, const environment& rhs)
{
  std::uniform_real_distribution<double> lhs_dist = lhs.get_cue_distribtion();
  std::uniform_real_distribution<double> rhs_dist = rhs.get_cue_distribtion();

  bool cue_distrib = are_same_distribution(lhs_dist, rhs_dist);

  bool env_function_A = are_same_env_functions(lhs.get_env_function_A(), rhs.get_env_function_A());

  bool env_function_B = are_same_env_functions(lhs.get_env_function_B(), rhs.get_env_function_B());

  return cue_distrib && env_function_A && env_function_B;
}

void environment::switch_name_current_function()
{
    try
    {
        if(get_name_current_function() == 'A')
            m_name_current_function = 'B';
        else if(get_name_current_function() == 'B')
            m_name_current_function = 'A';
        else throw std::runtime_error("Problem in switching functions: current function has an invalid name");
    }  catch (const std::exception& exc)
    {
        std::cerr << exc.what() << std::endl;
        abort();
    }

}

namespace env {


double calculate_optimal(const environment &e, std::vector<double> input)
{
  return e.get_current_function()(input);
}

std::vector<double> create_n_inputs(int n_inputs)
{
  std::vector<double> input_vector;

  for(int i = 0; i != n_inputs; ++i){
      input_vector.push_back(1234.0);
    }
  return input_vector;
}

std::vector<double> create_n_inputs(environment& e, const int &n_inputs, std::mt19937_64 &rng)
{
  std::vector<double> input_vector(n_inputs);

  for(auto& cue : input_vector){
      cue = e.get_dist()(rng);
    }

  return input_vector;
}

std::vector<double> create_n_inputs(std::uniform_real_distribution<double> dist,
                                    const int &n_inputs, std::mt19937_64 &rng)
{
    std::vector<double> input_vector(n_inputs);

    for(auto& cue : input_vector){
        cue = dist(rng);
    }

    return input_vector;
}


void switch_env_function(environment &e)
{
    try
    {
        if(e.get_name_current_function()=='A'){
            e.change_env_function(e.get_env_function_B());
            e.switch_name_current_function();
          }
        else if (e.get_name_current_function()=='B'){
          e.change_env_function(e.get_env_function_A());
          e.switch_name_current_function();
        }
        else throw std::runtime_error("Error while switching functions: current function has invalid name");
    }
    catch (const std::exception& exc)
    {
        std::cerr << exc.what() << std::endl;
        abort();
    }

}

}


void to_json(nlohmann::json& j, const environment& e)
{
    j = nlohmann::json{{"start", e.get_cue_distribtion().min()}, {"end", e.get_cue_distribtion().max()}};
}
void from_json(const nlohmann::json& j, environment& e)
{
    env_param e_p;
    e_p.cue_distrib = {j.at("start"), j.at("end")};
    e = {e_p};
}

#ifndef NDEBUG
void test_environment() noexcept
{


#define FIX_ISSUE_5
    #ifdef FIX_ISSUE_5
        ///It is possible to create an arbitrary number of inputs #5
        {
            int n_inputs = 3;
            auto inputs = env::create_n_inputs(n_inputs);
            assert(size_t(n_inputs) == inputs.size());
        }
    #endif
#define FIX_ISSUE_9
#ifdef FIX_ISSUE_9
    {
        ///As a first thing, make sure that now environment has a member variable that is a distribution
        ///so that when we construct an environment it already has a built in distribution
        environment e{env_param{}};

        ///Let's create a distribution that we can use to compare the distribution
        std::uniform_real_distribution<double> test_dist(-1,1);

        ///this is a random engine it is the source of randomness that you can plug inside distribution to generate random numbers with certain characteristic
        std::mt19937_64 rng;

        ///We are going to draw numbers from both from the env distribution and the test
        /// and then we are going to store them into vectors
        int repeats = 300000;
        std::vector<double> test_distr_values;
        std::vector<double> env_distr_values;

        for(int i = 0; i != repeats; i++)
        {
            test_distr_values.push_back(test_dist(rng));
            env_distr_values.push_back(e.get_dist()(rng));
        }

        ///We then calculate the mean and standard deviation of
        ///both the numbers drawn from the test and env distribution
        ///and check that they are approximately the same

        auto mean_test = calc_mean(test_distr_values);
        auto stdev_test = calc_stdev(test_distr_values);

        auto mean_env = calc_mean(env_distr_values);
        auto stdev_env = calc_stdev(env_distr_values);

        assert(are_equal_with_more_tolerance(mean_env,mean_test) &&
               are_equal_with_more_tolerance(stdev_env,stdev_test));
    }
#endif

#define FIX_ISSUE_14
#ifdef FIX_ISSUE_14
    {
        environment e{env_param{}};
        int n_inputs = 3;
        std::mt19937_64 rng1;

        auto tester_dist = e.get_dist();
        std::mt19937_64 rng2;

        std::vector<std::vector<double>> env_series(0, std::vector<double>(n_inputs));
        std::vector<std::vector<double>> tester_series(0, std::vector<double>(n_inputs));
        std::vector<std::vector<double>> tester_series1(0, std::vector<double>(n_inputs));

        int repeats = 10000;
        for(int i = 0; i != repeats; i++)
        {
            env_series.push_back(env::create_n_inputs(e, n_inputs, rng1));

            auto tester_cues = std::vector<double>(n_inputs);
            for(auto& cue : tester_cues){ cue = tester_dist(rng2);}
            tester_series.push_back(tester_cues);
        }

        assert(env_series == tester_series);
    }
#endif



#define FIX_ISSUE_52
#ifdef FIX_ISSUE_52

    ///Environment creates new inputs based on its own distribution
    {
        std::mt19937_64 rng;
        auto test_rng = rng;

        environment e{env_param{}};
        auto test_dist = e.get_dist();

        int n_of_inputs_requested = 1;

        std::vector<double> store_test_inputs;
        std::vector<double> store_new_inputs;

        int repeats = 30000;
        for(int i = 0; i != repeats; i++)
        {
            const auto new_inputs = env::create_n_inputs(e, n_of_inputs_requested, rng);
            store_new_inputs.insert(store_new_inputs.end(), new_inputs.begin(), new_inputs.end());

            for(int j = 0; j != n_of_inputs_requested; j++)
            {
                store_test_inputs.push_back(test_dist(test_rng));
            }
        }

        assert(are_from_same_distribution(store_new_inputs,store_test_inputs));

    }
#endif

  #define FIX_ISSUE_11
  #ifdef FIX_ISSUE_11
      {
          environment e{env_param{}};
           std::function<double(std::vector<double>)> env_function = e.get_env_function_A();
        std::vector<double> silly_argument{0.123456,0.98765443};
          env_function(silly_argument);
      }
  #endif



    #define FIX_ISSUE_28
    #ifdef FIX_ISSUE_28
            {
              //Two equal environments returns true
                environment lhs{env_func_1, env_func_2};
                environment rhs = lhs;

                assert(lhs == lhs);

                //Two environments that differ in their cue distribution returns false
                std::uniform_real_distribution<double> new_dist{1.23, 4.56};
                rhs.change_uniform_dist(new_dist);
                assert (!(lhs == rhs));


              //Two environments that differ in their function A or B returns false
                environment rhs2{env_func_2, env_func_2};
                environment rhs3{env_func_1, env_func_1};

                assert (!(lhs == rhs2));
                assert (!(lhs == rhs3));

            }
    #endif


 #define FIX_ISSUE_34
  #ifdef FIX_ISSUE_34
          {
              environment e{env_param{}};
              int n_inputs = 1;
              std::mt19937_64 rng1;
              std::mt19937_64 rng2;
              std::mt19937_64 rng3;

              std::vector<double> env_t0_series;
              std::vector<double> env_t1_series;
              std::vector<double> tester_t1_series;

              int repeats = 10000;
              for(int i = 0; i != repeats; i++)
              {
                  env_t0_series.push_back(env::create_n_inputs(e, n_inputs, rng1)[0]);
              }

              std::uniform_real_distribution<double> new_dist{1.23, 4.56};
              e.change_uniform_dist(new_dist);

              for(int i = 0; i != repeats; i++)
              {
                  env_t1_series.push_back(env::create_n_inputs(e, n_inputs, rng2)[0]);
                  tester_t1_series.push_back(env::create_n_inputs(new_dist, n_inputs, rng3)[0]); //This function working from a distribution is on a newer branch
              }

              auto env_t0_mean = calc_mean(env_t0_series);
              auto env_t1_mean = calc_mean(env_t1_series);
              auto tester_t1_mean = calc_mean(tester_t1_series);

              auto env_t0_stdev = calc_stdev(env_t0_series);
              auto env_t1_stdev = calc_stdev(env_t1_series);
              auto tester_t1_stdev = calc_stdev(tester_t1_series);

                      assert(are_not_equal_with_more_tolerance(env_t1_mean,env_t0_mean) &&
                             are_not_equal_with_tolerance(env_t0_stdev, env_t1_stdev));
                      assert(are_equal_with_more_tolerance(env_t1_mean,tester_t1_mean) &&
                             are_equal_with_tolerance(tester_t1_stdev, env_t1_stdev));

          }
  #endif


#define FIX_ISSUE_57
#ifdef FIX_ISSUE_57

    ///An enviroment can switch between optimum/al functions
    {
      environment e{env_param{}};
      assert(are_same_env_functions(e.get_current_function(), e.get_env_function_A()));
      env::switch_env_function(e); //Changed the name to match what we've been using
      assert(are_same_env_functions(e.get_current_function(), e.get_env_function_B()));
      assert(!are_same_env_functions(e.get_current_function(), e.get_env_function_A()));
      env::switch_env_function(e);
      assert(are_same_env_functions(e.get_current_function(), e.get_env_function_A()));
      assert(!are_same_env_functions(e.get_current_function(), e.get_env_function_B()));
    }
#endif


  #define FIX_ISSUE_58
  #ifdef FIX_ISSUE_58

      ///An environment can store 2 different optimal functions
      {
          environment e{env_param{}};
          auto func_A = e.get_env_function_A();
          auto func_B = e.get_env_function_B();
          assert(!are_same_env_functions(func_A,func_B));
      }
  #endif

#define FIX_ISSUE_77
#ifdef FIX_ISSUE_77

    ///The current function matches with the name of the current function
  {
          //Create an environment with non-default functions
          env_param param{};
          param.env_function_A = env_func_2;
          param.env_function_B = env_func_1;
          environment e{param};

          std::function<double(std::vector<double>)> current_function = e.get_current_function();
          assert(e.get_name_current_function() == 'A' && are_same_env_functions(current_function, env_func_2));

          env::switch_env_function(e);

          current_function = e.get_current_function();
          assert(e.get_name_current_function() == 'B' && are_same_env_functions(current_function, env_func_1));


      }
  #endif

#define FIX_ISSUE_143
#ifdef FIX_ISSUE_143

    ///A cue distribution param in environmental parameters can be used to generate an environment with the given distribution
    {
        env_param param{};
        range distrib{-213,123};
        param.cue_distrib = distrib;
        environment e{param};

        std::uniform_real_distribution<double> test_distrib(distrib.m_start, distrib.m_end);
        assert(are_same_distribution(e.get_cue_distribtion(), test_distrib));

    }
#endif


}
#endif
