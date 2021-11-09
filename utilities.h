#ifndef UTILITIES_H
#define UTILITIES_H
#include<vector>
#include<sstream>
#include <random>
#include <functional>

bool are_equal_with_tolerance(double lhs, double rhs);

bool are_not_equal_with_tolerance(double lhs, double rhs);

bool are_equal_with_more_tolerance(double lhs, double rhs);

bool are_not_equal_with_more_tolerance(double lhs, double rhs);

///Claculates mean of a vector of doubles
double calc_mean(const std::vector<double> &numbers);

///Calculates stdev of vector of doubles
double calc_stdev(const std::vector<double>& numbers);

///Converts a net_arc vec of int to a string
const std::string convert_arc_to_string(const std::vector<int>& v);

///Checks if two vectors of doubles come from the same distribution
bool are_from_same_distribution(const std::vector<double> &lhs, const std::vector<double> &rhs);

///Checks if two distributions are the same (quicker than are_from_same_distribution)
bool are_same_distribution(std::uniform_real_distribution<double> lhs,
                           std::uniform_real_distribution<double> rhs, int n_repeat = 1000);

///Checks if two environmental functions are the same
bool are_same_env_functions(const std::function<double(std::vector<double>)> &lhs,
                            const std::function<double(std::vector<double>)> &rhs);

#endif // UTILITIES_H
