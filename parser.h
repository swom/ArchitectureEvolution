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


void test_parser();
#endif // PARSER_H
