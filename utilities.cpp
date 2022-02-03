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
                           std::uniform_real_distribution<double> rhs,
                           int n_repeat)
{
  std::vector<double> lhs_cues;
  std::vector<double> rhs_cues;
  rndutils::xorshift128 rng1;
  rndutils::xorshift128 rng2;

  int repeats = n_repeat;
  for(int i = 0; i != repeats; i++)
  {
      lhs_cues.push_back(lhs(rng1)) ;
      rhs_cues.push_back(rhs(rng2)) ;
  }

  return lhs_cues == rhs_cues;
}

bool are_same_env_functions(const std::function<double(std::vector<double>)> &lhs,
                            const std::function<double(std::vector<double>)> &rhs, int n_repeats)
{
  std::vector<std::vector<double>> input_series;

  for(int i = 0; i != n_repeats; ++i)
    {
      double j = i;
      std::vector<double> input{j+1, j+2, j+3};
      input_series.push_back(input);
    }

  std::vector<double> lhs_optimal;
  std::vector<double> rhs_optimal;

  int repeats = n_repeats;
  for(int i = 0; i != repeats; i++)
  {
      lhs_optimal.push_back(lhs(input_series[i]));
      rhs_optimal.push_back(rhs(input_series[i]));
  }

  return lhs_optimal == rhs_optimal;

}

