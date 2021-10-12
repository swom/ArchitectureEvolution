#include "parser.h"
#include <algorithm>
#include <cassert>
#include <map>

static const double not_parsed = -123456.3142;

std::vector<int> arch_str_to_arch_vec(std::string net_arc)
{
    std::vector<int> net_arc_vec;
    size_t pos = 0;
    std::string delimiter = "-";
    std::string token;
    while ((pos = net_arc.find(delimiter)) != std::string::npos) {
        token = net_arc.substr(0, pos);
        net_arc_vec.push_back(std::stod(token));
        net_arc.erase(0, pos + delimiter.length());
    }
    net_arc_vec.push_back(std::stod(net_arc));

    return net_arc_vec;
}

///NOT tested!!!
env_param convert_env_args(const cxxopts::ParseResult& results)
{
    return env_param{
        results["targetA"].as<double>(),
                results["targetB"].as<double>()
    };
}

///NOT tested!!!
ind_param convert_ind_args(const cxxopts::ParseResult& results)
{
    return ind_param{
        convert_net_args(results),
                0
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
                results["mut_rate"].as<double>(),
                results["mut_step"].as<double>()
    };
}

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
            ("A,targetA", "the value fo env target A", cxxopts::value<double>()->default_value("0.1"))
            ("B,targetB", "the value fo env target B", cxxopts::value<double>()->default_value("0.75"))
            ("N,net_arc", "the network architecture", cxxopts::value<std::vector<int>>()->default_value("1,2,1"))
            ("F,act_func",
             "the string representing the name of the activation function of the net",
             cxxopts::value<std::string>()->default_value("sigmoid"))
            ("R,mut_rate",
             "the probability with whihc a mutation can happen",
             cxxopts::value<double>()->default_value("0.01"))
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
            ("h, help",
             "explains the stuff")
            ;
    return options;
}
std::function<double(double)> parse_act_func(const std::vector<std::string>& args)
{
    auto value = std::find(args.begin(), args.end(),"--mut_rate");
    if(value != args.end())
    {
        return string_to_act_func_map.find(*(value + 1))->second;
    }
    else{
        return linear;
    }
}

double parse_change_freq(const std::vector<std::string>& args)
{
    auto value = std::find(args.begin(), args.end(),"--change_freq");
    if(value != args.end())
    {
        return std::stod(*(value + 1));
    }
    else{
        return not_parsed;
    }
}

env_param parse_env_param(const std::vector<std::string>& args)
{
    env_param e_p{};
    e_p.targetA = parse_targetA(args);
    e_p.targetB = parse_targetB(args);
    return e_p;
}

ind_param parse_ind_param(const std::vector<std::string>& args)
{
    return ind_param{parse_net_param(args), 0};

}

double parse_mut_rate(const std::vector<std::string>& args)
{
    auto value = std::find(args.begin(), args.end(),"--mut_rate");
    if(value != args.end())
    {
        return std::stod(*(value + 1));
    }
    else{
        return not_parsed;
    }
}

double parse_mut_step(const std::vector<std::string>& args)
{
    auto value = std::find(args.begin(), args.end(),"--mut_step");
    if(value != args.end())
    {
        return std::stod(*(value + 1));
    }
    else{
        return not_parsed;
    }
}

std::vector<int> parse_net_arc(const std::vector<std::string>& args)
{
    auto value = std::find(args.begin(), args.end(),"--net_arc");
    std::vector<int> net_arc;
    if(value != args.end())
    {
        return arch_str_to_arch_vec(*(value + 1));
    }
    else{
        return net_arc;
    }
}

net_param parse_net_param(const std::vector<std::string>& args)
{
    net_param n_p;
    n_p.function = parse_act_func(args);
    n_p.net_arc = parse_net_arc(args);
    return n_p;
}

int parse_n_generations(const std::vector<std::string>& args)
{
    auto value = std::find(args.begin(), args.end(),"--n_generations");
    if(value != args.end())
    {
        return std::stod(*(value + 1));
    }
    else{
        return not_parsed;
    }
}

pop_param parse_pop_param(const std::vector<std::string>& args)
{
    pop_param p_p{};
    p_p.mut_rate = parse_mut_rate(args);
    p_p.mut_step = parse_mut_step(args);
    p_p.number_of_inds = parse_pop_size(args);
    return p_p;
}

int parse_pop_size(const std::vector<std::string>& args)
{
    auto value = std::find(args.begin(), args.end(),"--pop_size");
    if(value != args.end())
    {
        return std::stoi(*(value + 1));
    }
    else{
        return not_parsed;
    }
}

int parse_seed(const std::vector<std::string>& args)
{
    auto value = std::find(args.begin(), args.end(),"--seed");
    if(value != args.end())
    {
        return std::stoi(*(value + 1));
    }
    else{
        return not_parsed;
    }
}


int parse_sel_str(const std::vector<std::string>& args)
{
    auto value = std::find(args.begin(), args.end(),"--sel_str");
    if(value != args.end())
    {
        return std::stoi(*(value + 1));
    }
    else{
        return not_parsed;
    }
}

sim_param parse_sim_param(const std::vector<std::string>& args)
{
    sim_param s_p{};
    s_p.change_freq = parse_change_freq(args);
    s_p.seed = parse_seed(args);
    s_p.selection_strength = parse_sel_str(args);
    s_p.n_generations = parse_n_generations(args);
    return s_p;
}

double parse_targetA(const std::vector<std::string>& args)
{
    auto value = std::find(args.begin(), args.end(),"--targetA");
    if(value != args.end())
    {
        return std::stod(*(value + 1));
    }
    else{
        throw "Target A not found";
        return not_parsed;
    }
}

double parse_targetB(const std::vector<std::string>& args)
{
    auto value = std::find(args.begin(), args.end(),"--targetB");
    if(value != args.end())
    {
        return std::stod(*(value + 1));
    }
    else{
        return not_parsed;
    }
}

#ifndef NDEBUG
void test_parser()
{
    ///The change of frequency of environment in the simualtion
    /// can be taken as an argument
    /// from a commnad line
    {

        double expected_change_freq = 0.321456;
        std::string string_expected_change_freq = std::to_string(expected_change_freq);

        //will not take it because -- is missing from change_freq
        const std::vector<std::string>& args_incorrect{
            "random",
            "gibber",
            "change_freq",
            string_expected_change_freq
        };


        assert(parse_change_freq(args_incorrect) != expected_change_freq);

        const std::vector<std::string>& args_correct{
            "random",
            "gibber",
            "--change_freq",
            string_expected_change_freq
        };
        assert(parse_change_freq(args_correct) == expected_change_freq);
    }

    /// The architecture of the network
    /// can be taken as an argument
    /// from a commnad line
    {

        std::vector<int> expected_net_architecture{1,23,45};
        std::string string_expected_architecture = "1-23-45";

        //will not take it because -- is missing from change_freq
        const std::vector<std::string>& args_incorrect{
            "random",
            "gibber",
            "net_arc",
            string_expected_architecture
        };


        assert(parse_net_arc(args_incorrect) != expected_net_architecture);

        const std::vector<std::string>& args_correct{
            "random",
            "gibber",
            "--net_arc",
            string_expected_architecture
        };
        assert(parse_net_arc(args_correct) == expected_net_architecture);
    }
    /// The targetA environment value
    /// can be taken as an argument
    /// from a commnad line
    {

        double expected_target_A{123.45};
        std::string string_expected_targetA = std::to_string(expected_target_A);

        //will not take it because -- is missing from change_freq
        const std::vector<std::string>& args_incorrect{
            "random",
            "gibber",
            "targetA",
            string_expected_targetA
        };


        assert(parse_targetA(args_incorrect) != expected_target_A);

        const std::vector<std::string>& args_correct{
            "random",
            "gibber",
            "--targetA",
            string_expected_targetA
        };
        assert(parse_targetA (args_correct) == expected_target_A);
    }

    /// The targetB environment value
    /// can be taken as an argument
    /// from a commnad line
    {

        double expected_target_B{123.45};
        std::string string_expected_targetB = std::to_string(expected_target_B);

        //will not take it because -- is missing from change_freq
        const std::vector<std::string>& args_incorrect{
            "random",
            "gibber",
            "targetA",
            string_expected_targetB
        };


        assert(parse_targetB(args_incorrect) != expected_target_B);

        const std::vector<std::string>& args_correct{
            "random",
            "gibber",
            "--targetB",
            string_expected_targetB
        };
        assert(parse_targetB(args_correct) == expected_target_B);
    }
    /// The population size
    /// can be taken as an argument
    /// from a commnad line
    {

        double expected_pop_size{12345};
        std::string string_expected_pop_size = std::to_string(expected_pop_size);

        //will not take it because -- is missing from change_freq
        const std::vector<std::string>& args_incorrect{
            "random",
            "gibber",
            "pop_size",
            string_expected_pop_size
        };


        assert(parse_pop_size(args_incorrect) != expected_pop_size);

        const std::vector<std::string>& args_correct{
            "random",
            "gibber",
            "--pop_size",
            string_expected_pop_size
        };
        assert(parse_pop_size(args_correct) == expected_pop_size);
    }
    /// The seed
    /// can be taken as an argument
    /// from a commnad line
    {

        double expected_seed{12345};
        std::string string_expected_seed = std::to_string(expected_seed);

        //will not take it because -- is missing from change_freq
        const std::vector<std::string>& args_incorrect{
            "random",
            "gibber",
            "seed",
            string_expected_seed
        };


        assert(parse_seed(args_incorrect) != expected_seed);

        const std::vector<std::string>& args_correct{
            "random",
            "gibber",
            "--seed",
            string_expected_seed
        };
        assert(parse_seed(args_correct) == expected_seed);
    }
    /// The selection strength
    /// can be taken as an argument
    /// from a commnad line
    {

        double expected_selection_strength{12345};
        std::string string_expected_selection_strength = std::to_string(expected_selection_strength);

        //will not take it because -- is missing from change_freq
        const std::vector<std::string>& args_incorrect{
            "random",
            "gibber",
            "sel_str",
            string_expected_selection_strength
        };


        assert(parse_sel_str(args_incorrect) != expected_selection_strength);

        const std::vector<std::string>& args_correct{
            "random",
            "gibber",
            "--sel_str",
            string_expected_selection_strength
        };
        assert(parse_sel_str(args_correct) == expected_selection_strength);
    }

    /// The mutation rate
    /// can be taken as an argument
    /// from a commnad line
    {
        double mutation_rate{0.001};
        std::string string_expected_mutation_rate = std::to_string(mutation_rate);

        //will not take it because -- is missing from change_freq
        const std::vector<std::string>& args_incorrect{
            "random",
            "gibber",
            "mut_rate",
            string_expected_mutation_rate
        };


        assert(parse_mut_rate(args_incorrect) != mutation_rate);

        const std::vector<std::string>& args_correct{
            "random",
            "gibber",
            "--mut_rate",
            string_expected_mutation_rate
        };
        assert(parse_mut_rate(args_correct) == mutation_rate);
    }

    /// The activation function
    /// can be taken as an argument
    /// from a commnad line
    {
        auto activation_function = sigmoid;
        std::string string_expected_mutation_rate = "sigmoid";

        //will not take it because -- is missing from change_freq
        const std::vector<std::string>& args_incorrect{
            "random",
            "gibber",
            "activation_function",
            string_expected_mutation_rate
        };


        std::vector<double> random_values{0.23,564,123};
        for(const auto& value : random_values)
            assert(parse_act_func(args_incorrect)(value) == activation_function(value));

        const std::vector<std::string>& args_correct{
            "random",
            "gibber",
            "--activation_function",
            string_expected_mutation_rate
        };

        for(const auto& value : random_values)
            assert(parse_act_func(args_correct)(value) == activation_function(value));
    }
}
#endif
