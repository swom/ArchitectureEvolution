#include "environment.h"
#include "utilities.h"
#include <cassert>


environment::environment(double target_valueA, double target_valueB, std::function<double(std::vector<double>)> env_functionA):
    m_ref_target_values{target_valueA,target_valueB},
    m_current_target_value {target_valueA},
    m_input(3, 0.5),
    m_env_function_A{env_functionA}
{




}

environment::environment(env_param e_p):
    m_ref_target_values{e_p.targetA,e_p.targetB},
    m_current_target_value {e_p.targetA},
    m_cue_distribution{0., 1.},
    m_input(3, 0.5),
    m_env_function_A{e_p.env_function_A}
{




}



bool operator== (const environment& lhs, const environment& rhs)
{
  bool ref_t_values = lhs.get_ref_target_values() == rhs.get_ref_target_values();
  bool current_t_value = are_equal_with_tolerance(lhs.get_current_target_value(), rhs.get_current_target_value());

  return ref_t_values && current_t_value;
}


std::vector<double> environment::update_n_inputs(std::mt19937_64 &rng, const size_t n)
{
  std::vector<double> new_inputs = create_n_inputs(m_cue_distribution , n, rng);
  m_input = new_inputs;
  return new_inputs;
}

double environment::calculate_optimal()
{
  return get_env_function_A()(get_input());
}

void environment::update_optimal()
{
  m_optimal_output = calculate_optimal();
}


double get_target_valueA(const environment& e)
{
  return e.get_ref_target_values()[0];
}

double get_target_valueB(const environment& e)
{
  return e.get_ref_target_values()[1];
}

void switch_target(environment &e){
  //Check which target value is the current one and switch it over to the other

  if (are_equal_with_tolerance(e.get_current_target_value(),
                               get_target_valueA(e))
      )
    {
      e.set_current_target_value(get_target_valueB(e));
    }
  else
    {
      e.set_current_target_value(get_target_valueA(e));
    }
}

std::vector<double> create_n_inputs(int n_inputs)
{
  std::vector<double> input_vector;

  for(int i = 0; i != n_inputs; ++i){
      input_vector.push_back(1234.0);
    }
  return input_vector;
}

std::vector<double> create_n_inputs(environment e, const int &n_inputs, std::mt19937_64 &rng)
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


#ifndef NDEBUG
void test_environment() noexcept
{

  //an environment has a m_current_target_value member
  {
    double target_valueA = 0.123456;
    double target_valueB = 0.654321;
    environment e{target_valueA, target_valueB};
    assert(e.get_current_target_value() < 0 || e.get_current_target_value() > 0);
  }


  //an env has 2 reference target values;
  {
    double target_valueA = 0.123456;
    double target_valueB = 0.654321;
    environment e{target_valueA, target_valueB};
    assert(e.get_ref_target_values().size() == 2);
  }

  //an env can be initialized with 2 reference target values
  {
    double target_valueA = 0.123456;
    double target_valueB = 0.654321;
    environment e{target_valueA, target_valueB};

    assert(get_target_valueA(e) - target_valueA < 0.0001
           && get_target_valueA(e) - target_valueA > -0.0001);

    assert(get_target_valueB(e) - target_valueB < 0.0001
           && get_target_valueB(e) - target_valueB > -0.0001);
  }



  //Current target value is initialized to the first of the 2 target values
  {
    double targetA = 0.123456;
    double targetB = 0.654321;
    environment e{targetA,targetB};
    assert(are_equal_with_tolerance(e.get_current_target_value(), targetA));
    assert(are_not_equal_with_tolerance(e.get_current_target_value(), targetB));
  }


  //ISSUE_25
  //An environment can switch target values
  {
    double targetA = 0.123456;
    double targetB = 0.654321;
    environment e{targetA,targetB};
    assert(are_equal_with_tolerance(e.get_current_target_value(), targetA));
    switch_target(e);
    assert(are_equal_with_tolerance(e.get_current_target_value(), targetB));
    switch_target(e);
    assert(are_equal_with_tolerance(e.get_current_target_value(), targetA));
  }


  //#define FIX_ISSUE_35

  {
    double targetA = 123456;
    double targetB = 46589;

    env_param e_p{targetA, targetB};
    environment e{e_p};
    assert(are_equal_with_tolerance(get_target_valueA(e), targetA));
    assert(are_equal_with_tolerance(get_target_valueB(e), targetB));

  }

#define FIX_ISSUE_5
    #ifdef FIX_ISSUE_5
        ///It is possible to create an arbitrary number of inputs #5
        {
            int n_inputs = 3;
            auto inputs = create_n_inputs(n_inputs);
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
        std::uniform_real_distribution<double> test_dist(0,1);

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
            env_series.push_back(create_n_inputs(e, n_inputs, rng1));

            auto tester_cues = std::vector<double>(n_inputs);
            for(auto& cue : tester_cues){ cue = tester_dist(rng2);}
            tester_series.push_back(tester_cues);
        }

        assert(env_series == tester_series);
    }
#endif

#define FIX_ISSUE_25
#ifdef FIX_ISSUE_25
    ///Environment can create n new inputs and update them #25
    {
        std::mt19937_64 rng;
        environment e{env_param{}};
        auto env_inp_t1 = e.get_input();

        ///environment shouldn't know how many inputs individuals require
        /// nor we want to construct it with a certain number of inputs
        /// so we will have to specify it
        int n_of_inputs_requested = env_inp_t1.size();

        e.update_n_inputs(rng, n_of_inputs_requested);
        auto env_inp_t2 = e.get_input();

        assert(env_inp_t1 != env_inp_t2);
        assert(env_inp_t1.size() == env_inp_t2.size());

        n_of_inputs_requested++;
        e.update_n_inputs(rng, n_of_inputs_requested);
        auto env_inp_t3 = e.get_input();

        assert(env_inp_t2.size() != env_inp_t3.size());
    }
#endif

#define FIX_ISSUE_26
#ifdef FIX_ISSUE_26
    ///Environment creates new inputs based on its own distribution
    {
        std::mt19937_64 rng;
        auto test_rng = rng;

        environment e{env_param{}};
        auto test_dist = e.get_dist();

        int n_of_inputs_requested = 1;
        std::vector<double> test_inputs;

        std::vector<double> store_env_inputs;
        std::vector<double> store_test_inputs;

        int repeats = 30000;
        for(int i = 0; i != repeats; i++)
        {
            e.update_n_inputs(rng, n_of_inputs_requested);
            auto env_inputs = e.get_input();
            store_env_inputs.insert(store_env_inputs.end(), env_inputs.begin(), env_inputs.end());

            for(int j = 0; j != n_of_inputs_requested; j++)
            {
                store_test_inputs.push_back(test_dist(test_rng));
            }
        }

        assert(are_from_same_distribution(store_env_inputs,store_test_inputs));

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
  
 #define FIX_ISSUE_23
  #ifdef FIX_ISSUE_23
      {
          environment e{env_param{}};

          std::vector<double> input = e.get_input();

      }
  #endif
  
 #define FIX_ISSUE_29
  #ifdef FIX_ISSUE_29
        {            
  
             environment e{env_param{}};

            auto optimal_output = e.get_optimal();

            assert(optimal_output > 0 || optimal_output < 0);

        }
    #endif

    //#define FIX_ISSUE_28
    #ifdef FIX_ISSUE_28
            {
              //Two equal environments returns true
                environment lhs{321, 654};
                environment rhs{123, 456};

                assert(lhs == lhs);

              //Two environments that differ in their reference target values returns false
              //Testing this would also require a function to change it
              //but idk if that's super relevant to have
              //So maybe for now this test can be slightly less rigorous


              //Two environments that differ in their current target values returns false
                rhs.set_current_target_value(lhs.get_current_target_value() + 123);
                //I checked and there is no way to do this via the current constructor. The set function already existed, I didn't create it.

                assert (!(lhs == rhs));
                rhs = lhs;


              //Two environments that differ in their cue distribution returns false
                std::uniform_real_distribution<double> new_dist{1.23, 4.56};
                rhs.change_uniform_dist(new_dist);
                assert (!(lhs == rhs));
                rhs = lhs;

              //Two environments that differ in their input  & optimal output returns false
                std::mt19937_64 rng;
                rhs.update_n_inputs(rng, 3);
                rhs.update_output();
                assert (!(lhs == rhs));
                rhs = lhs;

              //Two environments that differ in their function returns false
                double function(std::vector<double> input);
                //I am not sure this will work as the function has not been defined and cannot be here;
                //if it doesn't a function defined somewhere else might be needed
                environment rhs2{321, 654, function};

                assert (!(lhs == rhs2));

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
                  env_t0_series.push_back(create_n_inputs(e, n_inputs, rng1)[0]);
              }

              std::uniform_real_distribution<double> new_dist{1.23, 4.56};
              e.change_uniform_dist(new_dist);

              for(int i = 0; i != repeats; i++)
              {
                  env_t1_series.push_back(create_n_inputs(e, n_inputs, rng2)[0]);
                  tester_t1_series.push_back(create_n_inputs(new_dist, n_inputs, rng3)[0]); //This function working from a distribution is on a newer branch
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


#define FIX_ISSUE_30
#ifdef FIX_ISSUE_30
    //Environment can use its function to calculate the optimal value and store it
  {
    environment e{env_param{}};
   auto function = e.get_env_function_A();

    std::vector<double> inputs = e.get_input();

    double theoretical_optimal_output = function(inputs);
    double calculated_optimal_output = e.calculate_optimal();

    assert(theoretical_optimal_output == calculated_optimal_output);

    e.update_optimal();
    assert(theoretical_optimal_output == e.get_optimal());

  }
#endif

}
#endif
