#include "simulation.h"

#include <cassert>
#include <vector>

template<class Pop, enum env_change_type E, enum selection_type S>
simulation<Pop, E, S>::simulation(int init_pop_size,
                                  int seed,
                                  double t_change_interval,
                                  std::vector<int> net_arch,
                                  double sel_str,
                                  int number_of_generations):
    m_environment{},
    m_population{init_pop_size},
    m_n_generations{number_of_generations},
    m_seed{seed},
    m_t_change_env_distr_A{static_cast<double>(t_change_interval)},
    m_t_change_env_distr_B{static_cast<double>(t_change_interval)},
    m_sel_str{sel_str},
    m_change_freq_A{static_cast<double>(t_change_interval)},
    m_change_freq_B{static_cast<double>(t_change_interval)},
    m_input(net_arch[0], 1),
    m_optimal_output{1}
{
    m_rng.seed(m_seed);
    for(auto& ind : m_population.get_inds_nonconst())
    {
        ind = individual{net_param{net_arch, linear, net_arch}};
    }
}

namespace sim {

template<class Pop>
bool operator ==(const simulation<Pop> &lhs, const simulation<Pop> &rhs)
{
    bool pop = lhs.get_pop() == rhs.get_pop();
    bool env = lhs.get_env() == rhs.get_env();
    bool time = lhs.get_time() == rhs.get_time();
    bool sel_str = are_equal_with_tolerance(lhs.get_sel_str(), rhs.get_sel_str());
    bool change_freq_A = are_equal_with_tolerance(lhs.get_change_freq_A(), rhs.get_change_freq_A());
    bool change_freq_B = are_equal_with_tolerance(lhs.get_change_freq_B(), rhs.get_change_freq_B());

    return pop && env && time && sel_str && change_freq_A && change_freq_B;
}

template<class Sim>
void change_all_weights_nth_ind(Sim& s, size_t ind_index, double new_weight)
{
    auto new_net = change_all_weights_values_and_activations(get_nth_ind_net(s, ind_index), new_weight);
    change_nth_ind_net(s, ind_index, new_net);
}

template<class Sim>
const typename Sim::pop_t::ind_t & get_nth_ind(const Sim& s, size_t ind_index)
{
    return pop::get_nth_ind(s.get_pop(), ind_index);
}

template<class Sim>
double get_nth_ind_fitness(const Sim& s, const size_t ind_index)
{
    return pop::get_nth_ind_fitness(s.get_pop(), ind_index);
}

template<class Sim>
const typename Sim::pop_t::ind_t::net_t & get_nth_ind_net(const Sim& s, size_t ind_index)
{
    return pop::get_nth_ind_net(s.get_pop(), ind_index);
}

template<class Sim>
double find_min_fitness(const Sim &s)
{
    std::vector<double> inds = pop::extract_fitnesses( s.get_pop().get_inds());

    auto min_inds_fitness = std::min_element(inds.begin(), inds.end());

    return *min_inds_fitness;
}

template<class Sim>
const std::vector<double> &get_nth_individual_input(const Sim &s, const int n)
{

    return pop::get_nth_individual_input(s.get_pop(), n);
}

template<class Sim>
const std::vector<double> &get_current_input(const Sim &s)
{
    return s.get_input();
}

template<class Sim>
std::function<double(std::vector<double>)> get_current_env_function(const Sim &s)
{
    auto e = s.get_env();
    return e.get_current_function();
}

}

double identity_first_element(const std::vector<double> &vector)
{
    return vector[0];
}


#ifndef NDEBUG
void test_simulation() noexcept//!OCLINT test may be many
{

    using namespace sim;
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
        int pop_size = 1;
        int seed = 123456789;
        simulation s{pop_size,
                    seed};
        std::mt19937_64 copy_rng(seed);
        assert ( s.get_rng()() == copy_rng());

    }

    ///A simulation should have a t_change_env_distr bernoulli distribution
    /// that determines the interval of environmental change
    /// Default initialization value is 10
    {
        int pop_size = 1;
        int seed = 123456789;
        double t_change_interval = 0.2;

        simulation s {pop_size, seed, t_change_interval};
        std::bernoulli_distribution mockdistrotchange(static_cast<double>(t_change_interval));
        assert (s.get_t_change_env_distr_A() == mockdistrotchange);

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
        std::vector<int> net_arch{1,1,2,1};
        net_param n_p{net_arch, linear, net_arch};
        simulation s{all_params{{}, {n_p}, {1, 0, 0, 0, 0}, {}}};

        auto sim_1st_net = get_nth_ind_net(s, 0);
        auto expected_net = network{n_p};
        assert(sim_1st_net == expected_net);
    }

#define FIX_ISSUE_68
#ifdef FIX_ISSUE_68
    ///Ex issue #30 test
    ///#define FIX_ISSUE_30
    ///#ifdef FIX_ISSUE_30

    ///individuals in a pop are selected based on how closely they match the current_env_target_value
    {
        std::function<double(std::vector<double>)> identity{identity_first_element};

        ///default func_A is already identity
        env_param identity_env {};
        identity_env.env_function_A = identity;
        assert(are_same_env_functions(identity_env.env_function_A, identity));

        auto potential_identity_net_param = net_param{{1,1}, linear, {1,1}};
        network potential_identity_net {potential_identity_net_param};
        auto potental_identity_ind = ind_param{potential_identity_net_param};
        assert(!net_behaves_like_the_function(potential_identity_net, identity));

        auto identity_net = change_all_weights_values(network{potential_identity_net}, 1);
        assert(net_behaves_like_the_function(identity_net, identity));


        int pop_size = 2;
        auto minimal_pop = pop_param{pop_size, 0, 0, 0, 0};

        auto sim_p = sim_param{};
        sim_p.selection_strength = 2;


        simulation s{all_params{
                identity_env,
                        potental_identity_ind,
                        minimal_pop,
                        sim_p
            }};

        //give simulation a simple input
        //this will be used to calculate optimum
        //and also fed to individuals
        s.update_inputs({1});
        //give inputs to inds
        assign_inputs(s);

        //change target value to match output of ind 0 net
        size_t best_ind = 0;
        sim::change_nth_ind_net(s, best_ind, identity_net);
        auto best_net = get_nth_ind_net(s, best_ind);

        size_t worst_ind = 1;
        change_nth_ind_net(s, worst_ind, change_all_weights_values_and_activations(potential_identity_net,-1));
        auto worst_net = get_nth_ind_net(s,worst_ind);

        assert(net_behaves_like_the_function(best_net, identity_env.env_function_A));
        assert(!net_behaves_like_the_function(worst_net, identity_env.env_function_A));

        s.select_inds();

        //all inds should now have the network that matches the target values
        for(const auto& ind :get_inds(s))
        {
            assert(ind.get_net() == best_net);
        }
    }
#endif


#define FIX_ISSUE_73
#ifdef FIX_ISSUE_73

    ///Fitness of individuals is calculated based on how close they are to the current optimal output
    {

        std::function<double(std::vector<double>)> identity{identity_first_element};

        env_param identity_env_par {};
        identity_env_par.env_function_A = identity;

        auto potential_identity_net_param = net_param{{1,1}, linear, {1,1}};
        network potential_identity_net {potential_identity_net_param};
        auto identity_net = change_all_weights_values_and_activations(potential_identity_net, 1);

        assert(net_behaves_like_the_function(identity_net, identity));

        auto potential_identity_ind_par = ind_param{potential_identity_net_param};


        int pop_size = 2;
        auto minimal_pop = pop_param{pop_size, 0, 0, 0, 0};

        auto sim_p = sim_param{};
        sim_p.selection_strength = 2;


        simulation s{all_params{
                identity_env_par,
                        potential_identity_ind_par,
                        minimal_pop,
                        sim_p
            }};
        //give simulation a simple input
        //this will be used to calculate optimum
        //and also fed to individuals
        s.update_inputs({1});
        //give inputs to inds
        assign_inputs(s);

        size_t first_ind = 0;
        size_t second_ind = 1;
        change_nth_ind_net(s, first_ind, identity_net);

        s.calc_fitness();

        ///ind 0 response should match exactly the optimal output therefore it will have fitness 1 (max)
        auto first_ind_fit =  get_nth_ind_fitness(s, first_ind) ;
        assert(are_equal_with_tolerance( first_ind_fit, 1));

        ///ind 1 response is not the optimal output, therefore its fitness should be the lowest in all the population
        auto first_response = ind::response(get_nth_ind(s, 0), s.get_input())[0];
        auto second_response = ind::response(get_nth_ind(s, 1), s.get_input())[0];
        assert(!are_equal_with_tolerance(first_response, second_response));

        auto second_ind_fit =  get_nth_ind_fitness(s, second_ind) ;
        auto min_fit = find_min_fitness(s);

        assert(!are_equal_with_tolerance(first_ind_fit, second_ind_fit));
        assert(are_equal_with_tolerance(min_fit,second_ind_fit));

    }
#endif


    //#define FIX_ISSUE_34
    {
        simulation s;
        int repeats = 100000;
        int n_switches = 0;
        for(int i = 0; i != repeats; i++)
        {
            if(s.is_environment_changing())
            {
                n_switches++;
            }
        }

        auto expected_repeats = s.get_change_freq_A() * repeats;
        assert(n_switches - expected_repeats < 20 &&
               n_switches - expected_repeats > -20);
    }


    //#define FIX_ISSUE_38

    {
        //sim_par
        int seed = 10126789;
        double change_freq_A = 123789;
        double change_freq_B = 87654321;
        double selection_strength = 0.321546;
        int n_generations = 123465;

        sim_param  s_p{seed,
                    change_freq_A,
                    change_freq_B,
                    selection_strength,
                    n_generations};

        all_params params{{}, {}, {}, s_p};
        simulation s{params};

        //test sim
        assert(are_equal_with_tolerance(s.get_change_freq_A(), change_freq_A) &&
               are_equal_with_tolerance(s.get_change_freq_B(), change_freq_B) &&
               are_equal_with_tolerance(s.get_sel_str(), selection_strength) &&
               s.get_seed() == seed &&
               s.get_n_gen() == n_generations);
    }


    //#define FIX_ISSUE_39
    //#ifdef FIX_ISSUE_39

    {
        simulation s{0};
        environment &e = s.get_env();
        int repeats =  100000;
        auto previous_env_function = e.get_name_current_function();

        int number_of_env_change = 0;

        for( int i = 0; i != repeats; i++)
        {
            tick(s);
            if(previous_env_function != e.get_name_current_function())
            {
                previous_env_function = e.get_name_current_function();
                number_of_env_change++;
            }
        }


        auto expected_changes = s.get_change_freq_A() * repeats;
        assert( number_of_env_change - expected_changes < repeats / 1000 &&
                number_of_env_change - expected_changes > -repeats / 1000);
    }

#define FIX_ISSUE_40
#ifdef FIX_ISSUE_40
    {
        //create a non-default simulaiton
        env_param e_p{};
        e_p.cue_distrib = range{-0.123,0.123};
        net_param n_p;
        n_p.net_arc = {1,2,3,4,5,6};
        n_p.max_arc = {1,2,3,4,5,6};
        ind_param i_p{};
        i_p.net_par = n_p;
        pop_param p_p{};
        p_p.mut_rate_activation = 0.1234;
        sim_param s_p;
        s_p.change_freq_A = 0.12345;

        all_params a_p {e_p, i_p, p_p, s_p};
        simulation s{a_p};
        auto name = "sim_save_test";
        save_json(s, name);
        simulation loaded_s = load_json<simulation<>>(name);
        assert(s == loaded_s);
    }
#endif

#define FIX_ISSUE_4
#ifdef FIX_ISSUE_4
    {
        population p;
        int n_inputs = 3;
        auto inputs = env::create_n_inputs(n_inputs);
        pop::assign_new_inputs_to_inds(p,inputs);
        for(const auto& ind : p.get_inds())
        {
            assert(ind.get_input_values() == inputs);
        }
    }
#endif

#define FIX_ISSUE_17
#ifdef FIX_ISSUE_17
    {
        simulation s;
        tick(s);

        int repeats = 5;
        while(repeats != 0)
        {
            auto t1_inputs = sim::get_current_input(s);
            tick(s);
            auto t2_inputs = sim::get_current_input(s);
            assert(t1_inputs != t2_inputs);

            repeats--;
        }

    }
#endif

#define FIX_ISSUE_18
#ifdef FIX_ISSUE_18
    {
        simulation s;
        assert(all_individuals_have_same_input(s));
        auto input_t1 = get_inds_input(s);

        std::vector<double> new_input;

        assign_new_inputs_to_inds(s, new_input);
        assert(all_individuals_have_same_input(s));
        auto input_t2 = get_inds_input(s);

        assert(input_t1 != input_t2);

    }
#endif

#define FIX_ISSUE_54
    //Simulation passes on its inputs to individuals after updating them
#ifdef FIX_ISSUE_54
    {
        simulation s;
        std::vector<double> new_input{1,2,3};

        s.update_inputs(new_input);
        assign_inputs(s);
        assert(s.get_input() == get_inds_input(s));
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

        s.update_inputs(create_inputs(s));
        auto sim_inp_t2 = s.get_input();

        assert(sim_inp_t1 != sim_inp_t2);
        assert(sim_inp_t1.size() == sim_inp_t2.size());
        assert(sim_inp_t2.size() == get_inds_input_size(s));
    }
#endif

#define FIX_ISSUE_50
    ///Simulation can make environment create new inputs to update its own inputs with, based on the number of inputs individuals have
#ifdef FIX_ISSUE_50
    {
        simulation s;

        int repeats = 100000;
        std::vector<double> sim_values;
        std::vector<double> test_values;

        for(int i = 0; i != repeats; i++)
        {
            const auto sim_inputs_t1 = s.get_input();
            const auto new_inputs = create_inputs(s);

            assert(sim_inputs_t1 != new_inputs);
            assert(new_inputs.size() == get_inds_input_size(s));

            s.update_inputs(new_inputs);
            const auto sim_inputs_t2 = s.get_input();

            assert(sim_inputs_t2 == new_inputs);

            sim_values.insert(sim_values.end(), sim_inputs_t2.begin(), sim_inputs_t2.end());

            environment e = s.get_env();
            auto test_inputs = env::create_n_inputs(e, get_inds_input_size(s), s.get_rng());
            test_values.insert(test_values.end(), test_inputs.begin(), test_inputs.end());
        }

        assert(are_from_same_distribution(sim_values, test_values));
    }
#endif

#define FIX_ISSUE_63
#ifdef FIX_ISSUE_63
    ///A simulation can switch between optimum/al functions
    {
        simulation s;

        std::vector<double> silly_inputs{1.23456, 9.87654};
        s.update_inputs(silly_inputs);

        assert(are_equal_with_tolerance(calculate_optimal(s), env_func_1(silly_inputs)));
        switch_optimal_function(s);
        assert(!are_equal_with_tolerance(calculate_optimal(s), env_func_1(silly_inputs)));
        assert(are_equal_with_tolerance(calculate_optimal(s), env_func_2(silly_inputs)));
        switch_optimal_function(s);
        assert(are_equal_with_tolerance(calculate_optimal(s), env_func_1(silly_inputs)));
        assert(!are_equal_with_tolerance(calculate_optimal(s), env_func_2(silly_inputs)));
    }
#endif

#define FIX_ISSUE_138
#ifdef FIX_ISSUE_138

    ///There should be an input to signal whihc environment function is being used to calculate the optima
    {
        std::vector<int> net_arch{2,2,1};
        all_params params{{},{{net_arch, linear, net_arch}}, {}, {}};
        simulation s{params};

        environment& e = s.get_env();

        assert(e.get_name_current_function() == 'A' && s.get_input().back() == 1);
        perform_environment_change(s);
        assign_inputs(s);
        assert(e.get_name_current_function() == 'B' && s.get_input().back() == -1);
    }
#endif

#define FIX_ISSUE_152
#ifdef FIX_ISSUE_152

    ///Network response depends on the environmental indicator
    {
        net_param n_p{{2,2,1}, linear, {2,2,1}};
        ind_param i_p{n_p};
        all_params params{{},i_p, {1,0,0,0,0}, {}}; //without the constructed pop param, it initializes an empty pop :(
        simulation s{params};

        //Otherwise all weights are 0 and the response is always 0
        change_all_weights_nth_ind(s, 0, 1);

        ///The response should change when the environment changes.

        std::vector<double> responseA = ind::response(get_nth_ind(s, 0), s.get_input());
        perform_environment_change(s);
        assign_inputs(s);
        std::vector<double> responseB = ind::response(get_nth_ind(s, 0), s.get_input());

        assert(responseA != responseB);

    }
#endif

#define FIX_ISSUE_158
#ifdef FIX_ISSUE_158

    ///Simulations have two rngs, one for population and one for environment(always seeded with 0)
    {
        sim_param s_p{1, 0, 0, 0, 0}; //Simulation rng is seeded with 1
        all_params params{{}, {}, {}, s_p};
        simulation s{params};
        std::mt19937_64 rng_before = s.get_rng();

        assert(s.get_rng() != s.get_env_rng());

        s.is_environment_changing();
        assert(s.get_rng() == rng_before);

    }
#endif

    ///#263
    ///Selection can happen sporadically every n generations
    {
        int selection_freq = 5;
        int repeats = 10;

        sim_param s_p{};
        s_p.n_generations = repeats;
        s_p.selection_freq = selection_freq;
        s_p.selection_strength = 10;
        pop_param p_p;
        p_p.number_of_inds = 10000;
        p_p.mut_rate_weight = 0.5;
        p_p.mut_step = 0.1;
        all_params a_p{{},{}, p_p, s_p};

        simulation<population<>,
                env_change_type::symmetrical,
                selection_type::sporadic> s{a_p};

        double stdev_prev_pop;
        double stdev_pop;
        double avg_pop;
        double avg_prev_pop;

        for (int i = 0; i != repeats; i++)
        {
            tick(s);

            stdev_pop = pop::stdev_fitness(s.get_pop());
            avg_pop = pop::avg_fitness(s.get_pop());
            stdev_prev_pop = pop::stdev_fitness(s.get_pop().get_new_inds());
            avg_prev_pop = pop::avg_fitness(s.get_pop().get_new_inds());

            if(s.get_time() % s.get_sel_freq() == 0)
            {
                assert(stdev_pop < stdev_prev_pop);
                assert(avg_prev_pop < avg_pop);

                assert(!are_equal_with_high_tolerance(stdev_pop, stdev_prev_pop));
                assert(!are_equal_with_high_tolerance(avg_prev_pop, avg_pop));
            }
            else if(s.get_time() % s.get_sel_freq() == s.get_sel_freq() - 1)

            {
                assert(are_equal_with_high_tolerance(stdev_pop, stdev_prev_pop));
                assert(are_equal_with_high_tolerance(avg_prev_pop, avg_pop));
            }
        }
    }

#define FIX_ISSUE_249
#ifdef FIX_ISSUE_249
    ///Simulation changes environment/function according to the specific change frequency of each function
    {
        sim_param s_p{};
        s_p.change_freq_A = 0.1;
        s_p.change_freq_B = 0.01;
        all_params a_p{{},{}, {}, s_p};

        simulation<population<>, env_change_type::asymmetrical> s_asym{a_p};
        int repeats = 1000000;

        int n_switches_A = 0;
        for(int i = 0; i != repeats; i++)
        {
            if(s_asym.is_environment_changing())
            {
                n_switches_A++;
            }
        }
        auto expected_repeats_A = s_asym.get_change_freq_A() * repeats;
        assert(n_switches_A - expected_repeats_A < 200 &&
               n_switches_A - expected_repeats_A > -200);

        switch_optimal_function(s_asym);

        int n_switches_B = 0;
        for(int i = 0; i != repeats; i++)
        {
            if(s_asym.is_environment_changing())
            {
                n_switches_B++;
            }
        }
        auto expected_repeats_B = s_asym.get_change_freq_B() * repeats;
        assert(n_switches_B - expected_repeats_B < 200 &&
               n_switches_B - expected_repeats_B > -200);

        ///Test that in symmetrical mode changes from one env to the other
        /// happen with same frequency for both environments

        simulation<population<>, env_change_type::symmetrical> s_sym{a_p};

        n_switches_A = 0;
        for(int i = 0; i != repeats; i++)
        {
            if(s_sym.is_environment_changing())
            {
                n_switches_A++;
            }
        }
        expected_repeats_A = s_sym.get_change_freq_A() * repeats;
        assert(n_switches_A - expected_repeats_A < 200 &&
               n_switches_A - expected_repeats_A > -200);

        switch_optimal_function(s_sym);

        n_switches_B = 0;
        for(int i = 0; i != repeats; i++)
        {
            if(s_sym.is_environment_changing())
            {
                n_switches_B++;
            }
        }
        expected_repeats_B = s_sym.get_change_freq_A() * repeats;
        assert(n_switches_B - expected_repeats_B < 200 &&
               n_switches_B - expected_repeats_B > -200);
        //check they are more or less the same number
        assert(n_switches_A - n_switches_B < 300 &&
               n_switches_A - n_switches_B > -300 );
    }
#endif

#define FIX_ISSUE_261
#ifdef FIX_ISSUE_261
    {
        sim_param s_p{};
        s_p.change_freq_A = 0.1;
        s_p.change_freq_B = 0.01;
        s_p.n_generations = 10000;
        all_params a_p{{},{}, {}, s_p};

        simulation<population<>, env_change_type::regular> regular_sim{a_p};
        int repeats = 1000000;

        auto current_function_name = get_name_current_function(regular_sim);
        for(int i = 0; i != repeats; i++)
        {
            tick(regular_sim);
            if(std::fmod(regular_sim.get_time(), 1.0/s_p.change_freq_A)  == 0)
            {
                assert(current_function_name != get_name_current_function(regular_sim));
                current_function_name = get_name_current_function(regular_sim);
            }

        }
    }
#endif
}
#endif
