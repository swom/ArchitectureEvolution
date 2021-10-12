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
