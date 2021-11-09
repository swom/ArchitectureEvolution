#include "utilities.h"
#include<algorithm>
#include<cmath>
#include<numeric>

bool are_equal_with_tolerance(double lhs, double rhs)
{
  return lhs - rhs > -0.0001 && lhs - rhs < 0.0001;
}

bool are_not_equal_with_tolerance(double lhs, double rhs)
{
  return !are_equal_with_tolerance(lhs, rhs);
}

bool are_equal_with_more_tolerance(double lhs, double rhs)
{
  return lhs - rhs > -0.001 && lhs - rhs < 0.001;
}

bool are_not_equal_with_more_tolerance(double lhs, double rhs)
{
  return !are_equal_with_more_tolerance(lhs, rhs);
}

double calc_mean(const std::vector<double>& numbers){
  return std::accumulate(numbers.begin(),
                              numbers.end(), 0.0)/numbers.size();
}

double calc_stdev(const std::vector<double>& numbers)
{
  double accum = 0.0;
  auto mean = calc_mean(numbers);
  std::for_each (std::begin(numbers),
                 std::end(numbers),
                 [&](const double weight) {
      accum += (weight - mean) * (weight - mean);});

  return sqrt(accum / (numbers.size()-1));
}

const std::string convert_arc_to_string(const std::vector<int>& v)
{
std::stringstream ss;
for(size_t i = 0; i < v.size(); ++i)
{
  if(i != 0)
    ss << "-";
  ss << v[i];
}
return ss.str();
}

bool are_from_same_distribution(const std::vector<double> &lhs, const std::vector<double> &rhs)
{
  auto lhs_mean = calc_mean(lhs);
  auto rhs_mean = calc_mean(rhs);

  auto lhs_stdev = calc_stdev(lhs);
  auto rhs_stdev = calc_stdev(rhs);

  if(are_equal_with_more_tolerance(lhs_mean, rhs_mean) &&
         are_equal_with_tolerance(lhs_stdev, rhs_stdev))
    return true;
  else
    return false;
}

bool are_same_distribution(std::uniform_real_distribution<double>lhs,
                           std::uniform_real_distribution<double> rhs, int n_repeat)
{
  std::vector<double> lhs_cues;
  std::vector<double> rhs_cues;
  std::mt19937_64 rng1;
  std::mt19937_64 rng2;

  int repeats = n_repeat;
  for(int i = 0; i != repeats; i++)
  {
      lhs_cues.push_back(lhs(rng1)) ;
      rhs_cues.push_back(rhs(rng2)) ;
  }

  return lhs_cues == rhs_cues;
}

bool are_same_env_functions(const std::function<double(std::vector<double>)> &lhs,
                            const std::function<double(std::vector<double>)> &rhs)
{
    std::vector<double> weird_numbers_1{1234.5678, 98765.432, 7687687686};
    std::vector<double> weird_numbers_2{12345.678, 98.765432, 7687.687686};


  double lhs_optimal_1 = lhs(weird_numbers_1);
  double rhs_optimal_1 = rhs(weird_numbers_1);

  double lhs_optimal_2 = lhs(weird_numbers_2);
  double rhs_optimal_2 = rhs(weird_numbers_2);

  return lhs_optimal_1 == rhs_optimal_1 && lhs_optimal_2 == rhs_optimal_2;

}
