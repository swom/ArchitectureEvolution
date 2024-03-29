#include "simulation.h"

#include <cassert>
#include <vector>


bool operator==(const sim_param& lhs, const sim_param& rhs)
{
    bool seeds = lhs.seed == rhs.seed;
    bool change_freq_A = lhs.change_freq_A == rhs.change_freq_A;
    bool change_freq_B = lhs.change_freq_B == rhs.change_freq_B;
    bool selection_strength = lhs.selection_strength == rhs.selection_strength;
    bool n_generations = lhs.n_generations == rhs.n_generations;
    bool selection_freq = lhs.selection_freq == rhs.selection_freq;
    bool change_sym_type = lhs.change_sym_type == rhs.change_sym_type;
    bool change_freq_type = lhs.change_freq_type == rhs.change_freq_type;
    bool sel_type = lhs.sel_type == rhs.sel_type;

    return seeds && change_freq_A && change_freq_B &&
            selection_strength && n_generations &&
            selection_freq && change_sym_type &&
            change_freq_type && sel_type;
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

simulation<> assign_random_IDs_to_inds(simulation<> s, rndutils::xorshift128& rng)
{
    std::uniform_int_distribution dist(-100,100);
    std::for_each(s.get_inds_non_const().begin(), s.get_inds_non_const().end(),
                  [&](auto& ind){ind.set_rank(rng());});
    return s;
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

bool ancestor_ID_is_parent_ID(const simulation<>& s)
{
    return  ancestor_IDs(sim::get_inds(s)) == pop_IDs(sim::get_parents(s));
}

} // namespace sim::

double identity_first_element(const std::vector<double> &vector)
{
    return vector[0];
}

///Creates a simple simualtion with trial evaluation type
simulation<population<>,
env_change_symmetry_type::symmetrical,
env_change_freq_type::stochastic,
selection_type::constant,
adaptation_period::off,
evaluation_type::trial
> create_trial_simulation(const all_params& a_p = all_params{})
{
    return         simulation<population<>,
            env_change_symmetry_type::symmetrical,
            env_change_freq_type::stochastic,
            selection_type::constant,
            adaptation_period::off,
            evaluation_type::trial
            >{a_p};
}

#ifndef NDEBUG
void test_simulation() noexcept//!OCLINT test may be many
{

    using namespace sim;
    ///A simulation has a member of type population
    ///The population has a vector of individuals of size 1 by default
    {
        simulation s;
        assert(s.get_pop_non_const().get_inds().size() == 1u);
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
        auto s = create_trial_simulation();
        auto init_timer_value = s.get_time();
        tick(s);
        assert(s.get_time() == init_timer_value + 1);

        auto s2 = create_trial_simulation();
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
        auto minimal_pop = pop_param{pop_size, 0, 0, 0, 0, 1};

        auto sim_p = sim_param{};
        sim_p.selection_strength = 20;


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

        s.calc_fitness();
        s.reproduce();

        //all inds should now have the network that matches the target values
        for(const auto& ind : sim::get_inds(s))
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
        auto s = create_trial_simulation();
        int repeats = 100000;
        int n_switches = 0;
        for(int i = 0; i != repeats; i++)
        {
            s.increase_time();
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
        double change_freq_A = 0.123789;
        double change_freq_B = 0.87654321;
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
        all_params a_p{{},{}, pop_param(1),{}};
        auto s = create_trial_simulation(a_p);
        environment &e = s.get_env();
        int repeats =  100000;
        auto previous_env_function = e.get_name_current_function();

        int number_of_env_change = 0;

        for( int i = 0; i != repeats; i++)
        {
            s.increase_time();

            if(s.is_environment_changing()){
                perform_environment_change(s);
            }

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
        e_p.cue_range = range{-0.123,0.123};
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
        auto s = create_trial_simulation();
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
        auto s = create_trial_simulation();
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

        s.update_inputs(s.create_inputs());
        auto sim_inp_t2 = s.get_input();

        assert(sim_inp_t1 != sim_inp_t2);
        assert(sim_inp_t1.size() == sim_inp_t2.size());
        assert(sim_inp_t2.size() == get_inds_input_size(s));
    }
#endif

#define FIX_ISSUE_50
    ///Simulation can make environment create new inputs to update its own inputs with,
    /// based on the number of inputs individuals have
#ifdef FIX_ISSUE_50
    {
        std::cout << "   test issue 50" << std::endl;
        simulation s;

        int repeats = 1000000;
        std::vector<double> sim_values;
        std::vector<double> test_values;

        for(int i = 0; i != repeats; i++)
        {
            const auto sim_inputs_t1 = s.get_input();
            const auto new_inputs = s.create_inputs();

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
        std::cout << "   end test issue 50" << std::endl;

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

#define FIX_ISSUE_152
#ifdef FIX_ISSUE_152

    ///Network response depends on the environmental indicator
    {
        net_param n_p{{2,2,1}, linear, {2,2,1}};
        ind_param i_p{n_p};
        all_params params{{},i_p, {1,0,0,0,0}, {}}; //without the constructed pop param, it initializes an empty pop :(
        simulation s{params};

        //Otherwise all weights are 0 and the response is always 0
        s.change_all_weights_nth_ind( 0, 1);

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
        //Simulation rng is seeded with 1
        sim_param s_p{1};
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
        std::cout << "   test issue 263" << std::endl;

        int selection_freq = 5;
        int repeats = 10;

        sim_param s_p{};
        s_p.n_generations = repeats;
        s_p.selection_freq = selection_freq;
        s_p.selection_strength = 10;
        pop_param p_p;
        p_p.number_of_inds = 4000;
        p_p.mut_rate_weight = 0.5;
        p_p.mut_step = 0.1;
        all_params a_p{{},{}, p_p, s_p};

        simulation<population<>,
                env_change_symmetry_type::symmetrical,
                env_change_freq_type::stochastic,
                selection_type::sporadic,//make selection happen sporadically
                adaptation_period::off,
                evaluation_type::trial> s{a_p};


        double avg_pop;
        double avg_prev_pop;

        for (int i = 0; i != repeats; i++)
        {
            tick(s);

            avg_pop = pop::avg_fitness(s.get_pop_non_const());
            avg_prev_pop = pop::avg_fitness(s.get_pop_non_const().get_new_inds());

            if(s.get_time() % s.get_sel_freq() >= 0 &&
                    s.get_time() % s.get_sel_freq() < s.get_selection_duration())
            {
                assert(!are_equal_with_high_tolerance(avg_prev_pop,
                                                      avg_pop)
                       );
                assert(avg_prev_pop < avg_pop);
            }
            else if(s.get_time() % s.get_sel_freq() == s.get_sel_freq() - 1)
            {
                assert(are_equal_with_high_tolerance(avg_prev_pop,
                                                     avg_pop)
                       );
            }
        }
        std::cout << "   end test issue 263" << std::endl;

    }

#define FIX_ISSUE_249
#ifdef FIX_ISSUE_249
    ///Simulation changes environment/function according to the specific change frequency of each function
    {
        sim_param s_p{};
        s_p.change_freq_A = 0.1;
        s_p.change_freq_B = 0.01;
        all_params a_p{{},{}, {}, s_p};

        simulation<population<>, env_change_symmetry_type::asymmetrical> s_asym{a_p};
        int repeats = 1000000;

        int n_switches_A = 0;
        for(int i = 0; i != repeats; i++)
        {
            s_asym.increase_time();
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
            s_asym.increase_time();
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

        simulation<population<>, env_change_symmetry_type::symmetrical> s_sym{a_p};

        n_switches_A = 0;
        for(int i = 0; i != repeats; i++)
        {
            s_sym.increase_time();
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
            s_sym.increase_time();
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
    ///Environment can change regularly instead of stochastically
    {

        int repeats = 1000;
        sim_param s_p{};
        s_p.change_freq_A = 0.1;
        s_p.change_freq_B = 0.2;
        s_p.n_generations = repeats;
        all_params a_p{{},{}, {}, s_p};

        simulation<population<>,
                env_change_symmetry_type::symmetrical,
                env_change_freq_type::regular,//make environment change deterministically
                selection_type::constant,
                adaptation_period::off,
                evaluation_type::trial> regular_sim{a_p};

        auto current_function_name = get_name_current_function(regular_sim);
        for(int i = 0; i !=  s_p.n_generations; i++)
        {
            if(std::fmod(regular_sim.get_time() - 1, 1.0/s_p.change_freq_A) == 0 &&
                    regular_sim.get_time() != 1)
            {
                assert(current_function_name != get_name_current_function(regular_sim));
                current_function_name = get_name_current_function(regular_sim);
            }
            tick(regular_sim);
        }

    }
#endif

#define FIX_ISSUE_275
#ifdef FIX_ISSUE_275
    //The number of trials (output-optimal_output distance)
    //on which fitness can be calculated can be decided in the sim_parameters
    //by defualt the number of trials is 1
    //The final fitness of an individual is the cumulative sum of the rescaled distances
    //of an individual output to the optimal output
    {
        std::cout << "   test issue 275" << std::endl;
        int number_of_trials = 5;
        pop_param p_p1{};
        pop_param p_p2{};

        p_p1.number_of_inds = 2;
        p_p2.number_of_inds = p_p1.number_of_inds;

        assert(p_p1.n_trials == 1);
        p_p2.n_trials = number_of_trials;

        assert(p_p1.n_trials < p_p2.n_trials);

        env_param e_p;
        e_p.cue_range = {1,1};
        all_params a_p1{e_p, ind_param{}, p_p1, sim_param{}};
        all_params a_p2{e_p, ind_param{}, p_p2, sim_param{}};

        auto s1 = create_trial_simulation(a_p1);
        auto s2 = create_trial_simulation(a_p2);

        auto fitnesses_s1 = s1.evaluate_inds();
        auto fitnesses_s2 = s2.evaluate_inds();
        assert(fitnesses_s1 != fitnesses_s2);

        assert(std::equal(fitnesses_s1.begin(),
                          fitnesses_s1.end(),
                          fitnesses_s2.begin(),
                          [p_p1, p_p2](const double& v_s1, const double& v_s2)
        {return v_s1 * p_p2.n_trials / p_p1.n_trials == v_s2;})
               );
        std::cout << "   end test issue 275" << std::endl;

    }

    ///A simulation that is templated with a plastic response_type
    /// will create an extra_input that refers to the current environmental function
    {

        using net_t = network<mutation_type::weights, response_type::plastic>;
        using ind_t = individual<net_t>;
        using pop_t = population<ind_t>;
        using sim_t = simulation<pop_t,
        env_change_symmetry_type::symmetrical,
        env_change_freq_type::regular,
        selection_type::constant,
        adaptation_period::off,
        evaluation_type::trial
        >;

        sim_param s_p;
        s_p.n_generations = 10;
        s_p.change_freq_A = 1;

        pop_param p_p;
        p_p.number_of_inds = 1;

        std::vector<int> arc = {3,1};
        net_param n_p;
        n_p.net_arc = arc;
        n_p.max_arc = arc;

        sim_t s{all_params{{},{}, p_p, s_p}};

        while(s.get_time() != s.get_n_gen())
        {
            tick(s);
            auto created_input = s.create_inputs();
            assert(created_input.back() == s.get_number_for_current_env_function());
            assert(created_input.size() == get_inds_input_size(s) + 1);
        }
    }
#endif

}

#endif
