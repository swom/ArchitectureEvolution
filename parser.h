#ifndef PARSER_H
#define PARSER_H
#include <string>
#include <vector>
#include "observer.h"

//Version 2.1 of https://github.com/jarro2783/cxxopts
#include <cxxopts.hpp>

std::vector<int> arch_str_to_arch_vec(std::string net_arc);

///Creates a parser using the cxxopts parser
cxxopts::Options create_parser();

///Functions to convert the parsed arguments into parameter objects
///NOT TESTED!!!
env_param convert_env_args(const cxxopts::ParseResult& results);

ind_param convert_ind_args(const cxxopts::ParseResult& results);

net_param convert_net_args(const cxxopts::ParseResult& results);

pop_param convert_pop_args(const cxxopts::ParseResult& results);

sim_param convert_sim_args(const cxxopts::ParseResult& results);

///runs a simulation given a parsed set of command line arguments
/// that are read as parameters
void run_simulation_given_arguments(const cxxopts::ParseResult& results);

///Given parameters, creates a simulation
template<class Pop, env_change_type E>
simulation<Pop, E> create_simulation(const cxxopts::ParseResult& parameters)
{
  auto env = convert_env_args(parameters);
  auto ind = convert_ind_args(parameters);
  auto pop = convert_pop_args(parameters);
  auto sim = convert_sim_args(parameters);

  all_params params{
      env, ind, pop, sim
  };

  simulation<Pop, E> s{params};
  return s;
}


void test_parser();
#endif // PARSER_H
