#include "parser.h"
#include <algorithm>
#include <cassert>
#include <map>



///NOT tested!!!
env_param convert_env_args(const cxxopts::ParseResult& results)
{
    return env_param{
        results["env_func_A"].as<std::string>(),
                results["env_func_B"].as<std::string>(),
                results["cue_range"].as<std::vector<double>>()

    };
    }

    ///NOT tested!!!
    ind_param convert_ind_args(const cxxopts::ParseResult& results)
    {
    return ind_param{
    convert_net_args(results),
    string_to_mut_type_map.find(results["mutation_type"].as<std::string>())->second
};
}

///NOT tested!!!
net_param convert_net_args(const cxxopts::ParseResult& results)
{
    return net_param{
        results["net_arc"].as<std::vector<int>>(),
                string_to_act_func_map.find(results["act_func"].as<std::string>())->second,
                results["max_arc"].as<std::vector<int>>(),
                string_to_response_type_map.find(results["response_type"].as<std::string>())->second,
                results["cue_range"].as<std::vector<double>>(),
                results["n_reac_norm_points"].as<int>()

    };
}

///NOT tested!!!
pop_param convert_pop_args(const cxxopts::ParseResult& results)
{
    return pop_param{
        results["pop_size"].as<int>(),
                results["mut_rate_weight"].as<double>(),
                results["mut_step"].as<double>(),
                results["mut_rate_act"].as<double>(),
                results["mut_rate_dup"].as<double>(),
                results["num_trials"].as<int>(),
    };
}

///NOT tested!!!
sim_param convert_sim_args(const cxxopts::ParseResult& results)
{
    return sim_param{
        results["seed"].as<int>(),
                results["change_freq_A"].as<double>(),
                results["change_freq_B"].as<double>(),
                results["sel_str"].as<double>(),
                results["num_gens"].as<int>(),
                results["sel_freq"].as<int>(),
                results["selection_duration_prop_to_freq"].as<int>(),
                results["n_reac_norm_points"].as<int>(),
                results["adaptation_period_proportion"].as<int>(),
                string_to_env_change_symmetry_type_map.find(results["env_change_sym_type"].as<std::string>())->second,
                string_to_env_change_freq_type_map.find(results["env_change_freq_type"].as<std::string>())->second,
                string_to_sel_type_map.find(results["sel_type"].as<std::string>())->second,
                string_to_adapt_period_map.find(results["adapt_p"].as<std::string>())->second,
                string_to_eval_type_map.find(results["evaluation_type"].as<std::string>())->second
    };
}
///NOT tested!!!
obs_param convert_obs_args(const cxxopts::ParseResult& results)
{
    return obs_param{
        results["top_inds_proportion"].as<int>(),
                results["top_inds_registration_freq"].as<int>(),
                results["top_spec_registration_freq"].as<int>(),
                results["n_mutations"].as<int>(),
                results["all_inds_rn_rec_freq"].as<int>(),
                results["sample_ind_record_freq_to_sens"].as<int>()
    };
}

///NOT tested!!!
cxxopts::Options create_parser(){

    cxxopts::Options options("Switch Simulation",
                             "Insert the parameters for the simualtion and see if you can get a mutational switch to evolve");
    options.allow_unrecognised_options();
    options.add_options()
            ("a,env_func_A",
             "the starting env function A",
             cxxopts::value<std::string>()->default_value("1"))
            ("A,mut_rate_act",
             "the probability with whihc an activation mutation can happen",
             cxxopts::value<double>()->default_value("0.001"))
            ("b,env_func_B",
             "the starting env function B",
             cxxopts::value<std::string>()->default_value("2"))
            ("B,selection_duration_prop_to_freq",
             "the proportion of generations in which selection will happen"
             " over the interval of time between selection events",
             cxxopts::value<int>()->default_value("100"))
            ("C,change_freq_A",
             "the probability with which the target function A will change",
             cxxopts::value<double>()->default_value("0.01"))
            ("c,change_freq_B",
             "the probability with which the target function B will change",
             cxxopts::value<double>()->default_value("0.01"))
            ("D,mut_rate_dup",
             "the probability with whihc a duplication mutation can happen",
             cxxopts::value<double>()->default_value("0.0005"))
            ("d,cue_range",
             "the minimum and maximum of the distribution used to generate environmental cues",
             cxxopts::value<std::vector<double>>()->default_value("-1,1"))
            ("E,sample_ind_record_freq_to_sens",
             "The multiplier used to determine how often to sample individuals"
             " based on sensibilities based on the record frequency of the best individual",
             cxxopts::value<int>()->default_value("10"))
            ("e,env_change_sym_type",
             "type of symmetry of the environmental change that a simulation will undergo",
             cxxopts::value<std::string>()->default_value("symmetrical"))
            ("F,act_func",
             "the string representing the name of the activation function of the net",
             cxxopts::value<std::string>()->default_value("sigmoid"))
            ("f,sel_freq",
             "the number of generations after which selection happens in the sporadic selection scenario",
             cxxopts::value<int>()->default_value("1"))
            ("G,num_gens",
             "number of generations for which the simulation has to run",
             cxxopts::value<int>()->default_value("1000000"))
            ("g,adaptation_period_proportion",
             "the proportion of the total time in the simulation where selection will be constant,"
             " only active if adaptaion_period is on",
             cxxopts::value<int>()->default_value("10"))
            ("H,adapt_p",
             "if set to 'on' half of the runtime of a simulation"
             " will be spent in a stable environmetn and the other "
             "half will have changing environments",
             cxxopts::value<std::string>()->default_value("off")) //'h' is taken for --help
            ("I, all_inds_rn_rec_freq",
               "The number of generations after which reaction norms"
               "of all individuals are saved",
               cxxopts::value<int>()->default_value("0"))
            ("i,top_inds_proportion",
             "the number of the top n indiivduals that will eb stored",
             cxxopts::value<int>()->default_value("1"))
            ("M,mut_step",
             "the variance of the normal distribution from which mutation size is drawn",
             cxxopts::value<double>()->default_value("0.1"))
            ("m,mutation_type",
              "type of mutation that a network will undergo",
              cxxopts::value<std::string>()->default_value("weights"))
            ("N,net_arc",
             "the network architecture",
             cxxopts::value<std::vector<int>>()->default_value("1,2,1"))
            ("n,num_trials",
              "the of trials individuals undergo to calculate their fitness/performance score",
              cxxopts::value<int>()->default_value("1"))
            ("P,pop_size",
            "the numebr of individuals in the simulation",
            cxxopts::value<int>()->default_value("1000"))
            ("p,n_reac_norm_points",
               "the number of inputs on which the reaction norm of an individual will be measured",
               cxxopts::value<int>()->default_value("100"))
            ("Q, evaluation_type",
             "The way individuals performances are evaluated, "
             "either on some random values in a cue range -> trial,"
             "or on the entire range values -> full_rn",
             cxxopts::value<std::string>()->default_value("full_rn"))
            ("q, response_type",
             "the type of response the network will have:"
             "'consitutive' if the network does not receive a signal about the environmental function"
             "'plastic' if the network receives an additional signal representing the environmetnal function",
             cxxopts::value<std::string>()->default_value("constitutive"))
            ("R,top_inds_registration_freq",
              "the number of generations after which top individuals will be selected",
              cxxopts::value<int>()->default_value("1000"))
             ("r,top_spec_registration_freq",
              "the number of generations after which the mutational spectrum of individuals will be recorded",
              cxxopts::value<int>()->default_value("0"))
            ("S,seed",
             "the seed of the rng",
             cxxopts::value<int>()->default_value("0"))
            ("s,sel_type",
             "the type of seelction regime of the simulation can be 'constant' or 'sporadic'",
             cxxopts::value<std::string>()->default_value("constant"))
            ("T,sel_str",
              "the strenght of selection",
              cxxopts::value<double>()->default_value("2")) // 't' is taken for --test
            ("u,n_mutations",
             "the number of mutations each locus will undergo when measuring a mutational spectrum",
             cxxopts::value<int>()->default_value("100"))
            ("W,mut_rate_weight",
              "the probability with whihc a weight mutation can happen",
              cxxopts::value<double>()->default_value("0.01"))
            ("X,max_arc",
             "the maximum size of the network architecture",
             cxxopts::value<std::vector<int>>()->default_value("1,8,1"))
            ("z,env_change_freq_type",
             "type of frequency environmental change that a simulation will undergo",
             cxxopts::value<std::string>()->default_value("stochastic"))

            ("t,test",
             "run all tests")
            ("h, help",
             "explains the stuff")
            ;
    return options;
}


#ifndef NDEBUG

#endif

all_params convert_all_params(const cxxopts::ParseResult &results)
{
    auto env = convert_env_args(results);
    auto ind = convert_ind_args(results);
    auto pop = convert_pop_args(results);
    auto sim = convert_sim_args(results);

    return all_params {
        env, ind, pop, sim
    };

}
