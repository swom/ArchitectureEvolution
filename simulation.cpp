#include "simulation.h"

#include <cassert>
#include <vector>
#include <fstream>


simulation::simulation(double targetA, double targetB,
                       int init_pop_size,
                       int seed,
                       double t_change_interval,
                       std::vector<int> net_arch,
                       double sel_str,
                       int number_of_generations):
    m_environment{targetA, targetB},
    m_population{init_pop_size},
    m_n_generations{number_of_generations},
    m_seed{seed},
    m_t_change_env_distr{static_cast<double>(t_change_interval)},
    m_sel_str{sel_str},
    m_change_freq {static_cast<double>(t_change_interval)},
    m_input(net_arch[0], 0.5),
    m_optimal_output{1}
{
    m_rng.seed(m_seed);
    for(auto& ind : m_population.get_inds())
    {
        ind.get_net() = net_arch;
    }
}


simulation::simulation(all_params params):
    m_environment{params.e_p},
    m_population{params.p_p, params.i_p},
    m_n_generations{params.s_p.n_generations},
    m_seed{params.s_p.seed},
    m_t_change_env_distr{static_cast<double>(params.s_p.change_freq)},
    m_sel_str{params.s_p.selection_strength},
    m_change_freq {static_cast<double>(params.s_p.change_freq)},
    m_params {params},
    m_input(params.i_p.net_par.net_arc[0], 0.5),
    m_optimal_output{1}
{
    m_rng.seed(m_seed);
    for(auto& ind : m_population.get_inds())
    {
        ind.get_net() = params.i_p.net_par.net_arc;
    }
}

const std::vector<double>& simulation::get_inds_input() const
{
  assert(all_individuals_have_same_input(*this));
  return get_inds()[0].get_input_values();
}

const std::vector<individual> &simulation::get_inds() const
{
  return simulation::get_pop().get_inds();
}


bool operator ==(const simulation& lhs, const simulation& rhs)
{
    bool pop = lhs.get_pop() == rhs.get_pop();
    bool env = lhs.get_env() == rhs.get_env();
    bool time = lhs.get_time() == rhs.get_time();
    bool sel_str = are_equal_with_tolerance(lhs.get_sel_str(), rhs.get_sel_str());
    bool change_freq = are_equal_with_tolerance(lhs.get_change_freq(), rhs.get_change_freq());

    return pop && env && time && sel_str && change_freq;
}

double avg_fitness(const simulation& s)
{
    return avg_fitness(s.get_pop());
}

void calc_fitness(simulation& s)
{
    s.get_pop() = calc_fitness(s.get_pop(),
                               get_current_env_value(s),/*get_current_cues(), get_optimal env_value_based_on_cues*/
                               s.get_sel_str());
}

void change_all_weights_nth_ind(simulation& s, size_t ind_index, double new_weight)
{
    auto new_net = change_all_weights(get_nth_ind_net(s, ind_index), new_weight);
    change_nth_ind_net(s, ind_index, new_net);
}

void change_current_target_value(simulation& s, double new_target_value)
{
    s.get_env().set_current_target_value(new_target_value);
}

void change_nth_ind_net(simulation& s, size_t ind_index, const network& n)
{
    change_nth_ind_net(s.get_pop(), ind_index, n) ;
}

std::vector<individual> get_best_n_inds(const simulation& s, int n)
{
    return get_best_n_inds(s.get_pop(), n);
}

double get_current_env_value(const simulation&s)
{
    return s.get_env().get_current_target_value();
}

double get_current_env_value(simulation&s)
{
    return s.get_env().get_current_target_value();
}

const std::vector<individual>& get_inds(const simulation&s)
{
    return s.get_pop().get_inds();
}

const individual& get_nth_ind(const simulation& s, size_t ind_index)
{
    return get_nth_ind(s.get_pop(), ind_index);
}

double get_nth_ind_fitness(const simulation& s, const size_t ind_index)
{
    return get_nth_ind_fitness(s.get_pop(), ind_index);
}

const network& get_nth_ind_net(const simulation& s, size_t ind_index)
{
    return get_nth_ind_net(s.get_pop(), ind_index);
}

network& get_nth_ind_net( simulation& s, size_t ind_index)
{
    return get_nth_ind_net(s.get_pop(), ind_index);
}

double find_min_fitness(const simulation&s)
{
    auto inds = s.get_pop().get_inds();

    auto min_ind =
            std::min_element(inds.begin(), inds.end(), [](const individual& lhs, const individual& rhs){
        return lhs.get_fitness() < rhs.get_fitness();});


    return min_ind->get_fitness();
}

bool is_environment_changing (simulation &s) {

    std::bernoulli_distribution distro = s.get_t_change_env_distr();
    return distro (s.get_rng());

}

simulation load_json(
        const std::string& filename
        )
{
    std::ifstream f(filename);
    nlohmann::json json_in;
    simulation s;
    f >> json_in;
    return s = json_in;
}

void reproduce(simulation& s)
{
    reproduce(s.get_pop(), s.get_rng());
}

void tick(simulation &s)
{
    s.increase_time();

    if(is_environment_changing(s)){

        switch_target(s.get_env());
    }

    if(s.get_inds().size()){

      //assign_new_inputs(s);
    }

    select_inds(s);


}

void save_json(const simulation& s, const std::string& filename)
{
    std::ofstream  f(filename);
    nlohmann::json json_out;
    json_out = s;
    f << json_out;
}

void select_inds(simulation& s)
{
    /*pop_perceives_en*/
    calc_fitness(s);
    reproduce(s);
}

double var_fitness(const simulation&s)
{
    return var_fitness(s.get_pop());
}

void update_inputs(population &p, const std::vector<double> &inputs)
{
  for(auto& ind : p.get_inds()){
      ind.assign_input(inputs);
    }
}

bool all_individuals_have_same_input(const simulation &s)
{
  population p = s.get_pop();

  return all_individuals_have_same_input(p);
}

const std::vector<double> &get_nth_individual_input(const simulation &s, const int n)
{

  return get_nth_individual_input(s.get_pop(), n);
}

const std::vector<double> &get_current_input(const simulation &s)
{
 assert(all_individuals_have_same_input(s));
 return get_nth_individual_input(s, 0);
}

double calculate_optimal(const simulation &s)
{
  return(calculate_optimal(s.get_env(), s.get_input()));
}

void simulation::update_inputs(){
  std::vector<double> new_inputs = create_inputs(*this);
  m_input = new_inputs;
}

std::vector<double> create_inputs(simulation &s)
{
  environment &e = s.get_env();
  return(create_n_inputs(e, s.get_inds_input_size(), s.get_rng() ));
}



#ifndef NDEBUG
void test_simulation() noexcept//!OCLINT test may be many
{

    ///A simulation has a member of type population
    ///The population has a vector of individuals of size 1 by default
    {
        simulation s;
        assert(s.get_pop().get_inds().size() == 1u);
    }


#ifdef FIX_ISSUE_27
    ///A simulation can be initialized by
    /// target values for the environment
    /// and the number of individuals in
    /// the populaation
    {
        double target_valueA = 1.23456;
        double target_valueB = 6.12345;
        int init_pop_size = 123456;
        simulation s{target_valueA, target_valueB ,init_pop_size};
        assert(s.get_env().get_target_valueA() == target_valueA
               && s.get_env().get_target_valueB() == target_valueB
               && s.get_pop().get_ind_vec().size() == static_cast<unsigned int>(init_pop_size));
    }
#endif
    ////A simulation should have a random engine, intialized with a seed that is fed to simulation

    {
        double targetA = 1.0;
        double targetB = 0;
        int pop_size = 1;
        int seed = 123456789;
        simulation s{targetA,
                    targetB,
                    pop_size,
                    seed};
        std::mt19937_64 copy_rng(seed);
        assert ( s.get_rng()() == copy_rng());

    }

    ///A simulation should have a t_change_env_distr bernoulli distribution
    /// that determines the interval of environmental change
    /// Default initialization value is 10
    {
        double targetA = 1.0;
        double targetB = 0;
        int pop_size = 1;
        int seed = 123456789;
        double t_change_interval = 0.2;

        simulation s { targetA, targetB,
                    pop_size, seed, t_change_interval};
        std::bernoulli_distribution mockdistrotchange(static_cast<double>(t_change_interval));
        assert (s.get_t_change_env_distr() == mockdistrotchange);

    }



    //A simulation has an absolute counter and you can access it
    {
        simulation s;
        assert((s.get_time() > 0) | (s.get_time() <= 0));
    }



    //Every tick simulation timer increases by one
    {
        simulation s;
        auto init_timer_value = s.get_time();
        tick(s);
        assert(s.get_time() == init_timer_value + 1);

        simulation s2;
        init_timer_value = s2.get_time();
        int repeats = 123;
        for(int i = 0; i != repeats; i++)
        {
            tick(s2);
        }
        assert(s2.get_time() == init_timer_value + repeats);
    }


    ///Simulation can be initialized with network architecture for inds in pop
    {
        std::vector<int> net_arch{1,33,2,1};
        simulation s{0,0,1,0,0, net_arch};

        assert(get_nth_ind_net(s, 0) == network{net_arch});
    }

#define FIX_ISSUE_30
#ifdef FIX_ISSUE_30
    ///individuals in a pop are selected based on how closely they match the current_env_target_value
    {
        int pop_size = 2;
        simulation s{0,0, pop_size};

        //change target value to match output of ind 0 net
        size_t best_ind = 0;
        change_current_target_value(s, response(get_nth_ind(s, best_ind))[0]);
        auto best_net = get_nth_ind_net(s, best_ind);
        auto resp_best = response(get_nth_ind(s, best_ind))[0];

        size_t worst_ind = 1;
        auto worst_net = change_all_weights(get_nth_ind_net(s,worst_ind), 100);
        change_nth_ind_net(s, worst_ind, worst_net);
        auto resp_worst = response(get_nth_ind(s, worst_ind))[0];

        assert(resp_best == get_current_env_value(s));
        assert(resp_worst != get_current_env_value(s));

        select_inds(s);

        //all inds should now have the network that matches the target values
        for(const auto& ind :get_inds(s))
        {
            assert(ind.get_net() == best_net);
        }
    }
#endif

    ///Fitness of individuals is calculated based on how close they are to the current target value
    {
        int pop_size = 2;
        simulation s{0,0,pop_size};

        size_t first_ind = 0;
        size_t second_ind = 1;
        change_all_weights_nth_ind(s, first_ind, 1);
        change_all_weights_nth_ind(s, second_ind, 0.99);


        //change target value to match output of ind 0 net
        change_current_target_value(s, response(get_nth_ind(s, 0))[0]);

        calc_fitness(s);

        ///ind 0 response should match exactly the target value therefore it will have fitness 1 (max)
        auto first_ind_fit =  get_nth_ind_fitness(s,0) ;
        assert(are_equal_with_tolerance( first_ind_fit, 1));

        ///ind 1 response is 0, therefore its fitness would be the lowest in all the population
        auto first_response = response(get_nth_ind(s, 0))[0];
        auto second_response = response(get_nth_ind(s, 1))[0];
        assert(!are_equal_with_tolerance(first_response, second_response));

        auto second_ind_fit =  get_nth_ind_fitness(s,1) ;
        auto min_fit = find_min_fitness(s);

        assert(are_equal_with_tolerance(min_fit,second_ind_fit));

    }

    //#define FIX_ISSUE_34
    {
        simulation s;
        int repeats = 100000;
        int n_switches = 0;
        for(int i = 0; i != repeats; i++)
        {
            if(is_environment_changing(s))
            {
                n_switches++;
            }
        }

        auto expected_repeats = s.get_change_freq() * repeats;
        assert(n_switches - expected_repeats < 20 &&
               n_switches - expected_repeats > -20);
    }


    //#define FIX_ISSUE_38

    {
        //sim_par
        int seed = 10126789;
        double change_freq = 123789;
        double selection_strength = 0.321546;
        int n_generations = 123465;

        sim_param  s_p{seed, change_freq, selection_strength, n_generations};
        all_params params{{}, {}, {}, s_p};
        simulation s{params};

        //test sim
        assert(are_equal_with_tolerance(s.get_change_freq(), change_freq) &&
               are_equal_with_tolerance(s.get_sel_str(), selection_strength) &&
               s.get_seed() == seed &&
               s.get_n_gen() == n_generations);
    }


#define FIX_ISSUE_39
#ifdef FIX_ISSUE_39

    {
        simulation s{0, 0.1, 0};
        int repeats =  100000;
        auto previous_env_value = get_current_env_value(s);

        int number_of_env_change = 0;

        for( int i = 0; i != repeats; i++)
        {
            tick(s);
            if(previous_env_value != get_current_env_value(s))
            {
                previous_env_value = get_current_env_value(s);
                number_of_env_change++;
            }
        }

        auto expected_changes = s.get_change_freq() * repeats;
        assert( number_of_env_change - expected_changes < repeats / 1000 &&
                number_of_env_change - expected_changes > -repeats / 1000);
    }
#endif

#define FIX_ISSUE_40
#ifdef FIX_ISSUE_40
    {
        //create a non-default simulaiton
        simulation s{13, 32, 2, 132, 548, {1,2,3,4,5,6}, 3.14};
        auto name = "sim_save_test";
        save_json(s, name);
        auto loaded_s = load_json(name);
        assert(s == loaded_s);
    }
#endif

#define FIX_ISSUE_4
 #ifdef FIX_ISSUE_4
     {
         population p;
         int n_inputs = 3;
         auto inputs = create_n_inputs(n_inputs);
         update_inputs(p,inputs);
         for(const auto& ind : p.get_inds())
         {
             assert(ind.get_input_values() == inputs);
         }
     }
   #endif

//#define FIX_ISSUE_17
#ifdef FIX_ISSUE_17
    {
        simulation s;
        tick(s);

        int repeats = 5;
        while(repeats != 0)
        {
        auto t1_inputs = get_current_input(s);
        tick(s);
        auto t2_inputs = get_current_input(s);
        assert(t1_inputs != t2_inputs);

        repeats--;
        }

    }
#endif

//#define FIX_ISSUE_18
#ifdef FIX_ISSUE_18
    {
        simulation s;
        assert(all_individuals_have_same_input(s));
        auto input_t1 = get_current_input(s);

        assign_new_inputs(s);
        assert(all_individuals_have_same_input(s));
        auto input_t2 = get_current_input(s);

        assert(input_t1 != input_t2);

    }
#endif

//#define FIX_ISSUE_54
#ifdef FIX_ISSUE_54
    {
        simulation s;
        s.update_inputs();
        s.update_ind_inputs();
        assert(s.get_input() == s.get_inds_input());
    }
#endif

//#define FIX_ISSUE_27
#ifdef FIX_ISSUE_27
    {
        simulation s;
        auto test_e = s.get_env();
        auto n_inputs_requested = s.get_inds_input_size();

        int repeats = 100000;
        std::vector<double> sim_env_values;
        std::vector<double> test_values;

        for(int i = 0; i != repeats; i++)
        {
            const auto env_inputs_t1 = s.get_env_inputs();
            create_inputs(s);
            const auto env_inputs_t2 = s.get_env_inputs();

            assert(env_inputs_t1 != env_inputs_t2);
            assert(env_inputs_t2.size() == s.get_inds_input_size());

            sim_env_values.insert(sim_env_values.end(), env_inputs_t2.begin(), env_inputs_t2.end());

            auto test_inputs = test_e.update_n_inputs(s.get_rng(), n_inputs_requested);
            test_values.insert(test_values.end(), test_inputs.begin(), test_inputs.end());
        }

        assert(are_from_same_distribution(sim_env_values, test_values));
    }
#endif


  
#define FIX_ISSUE_46
//Simulation has a private member that stores the input that are gonna be fed to individuals environment
 #ifdef FIX_ISSUE_46
     {
         simulation s{};

         std::vector<double> input = s.get_input();

     }
 #endif
  
 #define FIX_ISSUE_47
    //Simulation has a private member where it can store the optimal value of its inputs.
   #ifdef FIX_ISSUE_47
         {

              simulation s{};

             auto optimal_output = s.get_optimal();

             assert(optimal_output >= 0 || optimal_output < 0);

         }
     #endif

#define FIX_ISSUE_48
#ifdef FIX_ISSUE_48
    //Simulation can use the environment's function to calculate the optimal value and store it
  {
    simulation s{};
    auto function = s.get_env_function_A();

    std::vector<double> inputs = s.get_input();

    double theoretical_optimal_output = function(inputs);
    double calculated_optimal_output = calculate_optimal(s);

    assert(theoretical_optimal_output == calculated_optimal_output);

    s.update_optimal(calculated_optimal_output);
    assert(theoretical_optimal_output == s.get_optimal());
  }
#endif

#define FIX_ISSUE_49
#ifdef FIX_ISSUE_49
    ///Simulation can ask environment to create n new inputs to update its input with
    {
        simulation s {};
        auto sim_inp_t1 = s.get_input();

        s.update_inputs();
        auto sim_inp_t2 = s.get_input();

        assert(sim_inp_t1 != sim_inp_t2);
        assert(sim_inp_t1.size() == sim_inp_t2.size());
        assert(sim_inp_t2.size() == s.get_inds_input_size());
    }
#endif

  
}
#endif
