#include "parser.h"
#include <algorithm>
#include <cassert>
#include <map>

void run_simulation_given_arguments(const cxxopts::ParseResult& results)
{
    auto mut_type = convert_ind_args(results).m_mutation_type;

    if(mut_type == mutation_type::weights)
    {

        observer<mutation_type::weights> o;
        auto s = create_simulation<mutation_type::weights>(results);
        exec<mutation_type::weights>(s, o) ;
        save_json(o,
                  convert_arc_to_string(o.get_params().i_p.net_par.net_arc) +
                  "_" + std::to_string(o.get_params().s_p.seed) + ".json");
    }
    else if (mut_type == mutation_type::activation) {

        observer<mutation_type::activation> o;
        auto s = create_simulation<mutation_type::activation>(results);
        exec<mutation_type::activation>(s, o) ;
        save_json(o,
                  convert_arc_to_string(o.get_params().i_p.net_par.net_arc) +
                  "_" + std::to_string(o.get_params().s_p.seed) + ".json");
    }
    else if (mut_type == mutation_type::weights_and_activation) {
        observer<mutation_type::weights_and_activation> o;
        auto s = create_simulation<mutation_type::weights_and_activation>(results);
        exec<mutation_type::weights_and_activation>(s, o) ;
        save_json(o,
                  convert_arc_to_string(o.get_params().i_p.net_par.net_arc) +
                  "_" + std::to_string(o.get_params().s_p.seed) + ".json");
    }
    else
    {
        throw std::runtime_error{"unknown mutation type"};
    }
}

///NOT tested!!!
env_param convert_env_args(const cxxopts::ParseResult& results)
{
    return env_param{
                    string_env_function_map.find(results["env_func_A"].as<std::string>())->second,
                    string_env_function_map.find(results["env_func_B"].as<std::string>())->second,
                    results["cue_distrib"].as<std::vector<double>>()

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
                string_to_act_func_map.find(results["act_func"].as<std::string>())->second
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
    };
}

///NOT tested!!!
sim_param convert_sim_args(const cxxopts::ParseResult& results)
{
    return sim_param{
        results["seed"].as<int>(),
                results["change_freq"].as<double>(),
                results["sel_str"].as<double>(),
                results["num_gens"].as<int>()
    };
}
///NOT tested!!!
cxxopts::Options create_parser(){

    cxxopts::Options options("Switch Simulation",
                             "Insert the parameters for the simualtion and see if you can get a mutational switch to evolve");
    options.allow_unrecognised_options();
    options.add_options()
            ("a,env_func_A", "the starting env function A",cxxopts::value<std::string>()->default_value("1"))
            ("b,env_func_B", "the starting env function B",cxxopts::value<std::string>()->default_value("2"))
            ("N,net_arc", "the network architecture", cxxopts::value<std::vector<int>>()->default_value("1,2,1"))
            ("F,act_func",
             "the string representing the name of the activation function of the net",
             cxxopts::value<std::string>()->default_value("sigmoid"))
            ("W,mut_rate_weight",
             "the probability with whihc a weight mutation can happen",
             cxxopts::value<double>()->default_value("0.01"))
            ("A,mut_rate_act",
             "the probability with whihc an activation mutation can happen",
             cxxopts::value<double>()->default_value("0.001"))
            ("M,mut_step",
             "the variance of the normal distribution from which mutation size is drawn",
             cxxopts::value<double>()->default_value("0.1"))
            ("P,pop_size","the numebr of individuals in the simulation",cxxopts::value<int>()->default_value("1000"))
            ("C,change_freq",
             "the probability with which the target value will change",
             cxxopts::value<double>()->default_value("0.01"))
            ("S,seed",
             "the seed of the rng",
             cxxopts::value<int>()->default_value("0"))
            ("T,sel_str",
             "the strenght of selection",
             cxxopts::value<double>()->default_value("2"))
            ("G,num_gens",
             "number of generations for which the simulation has to run",
             cxxopts::value<int>()->default_value("1000000"))
            ("m,mutation_type",
"type ofg mutation that a network will undergo",
             cxxopts::value<std::string>()->default_value("weights"))
            ("d,cue_distrib",
             "the minimum and maximum of the distribution used to generate environmental cues",
             cxxopts::value<std::vector<double>>()->default_value("-1,1"))
            ("t,test",
             "run all tests")
            ("h, help",
             "explains the stuff")
            ;
    return options;
}


#ifndef NDEBUG

#endif
