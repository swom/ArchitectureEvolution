#include "weight.h"

#include <cassert>
#include <vector>

weight::weight(double weight_init, bool is_active):
  m_weight{weight_init},
  m_is_active{is_active}
{

}

bool operator== (const weight& lhs, const weight& rhs)
{
    return lhs.get_weight() == rhs.get_weight() &&
      lhs.is_active() == rhs.is_active();
}

bool operator!= (const weight& lhs, const weight& rhs)
{
    return !(lhs == rhs);
}

double operator* (double& number, const weight& weight_to_multiply)
{
  double new_weight = weight_to_multiply.get_weight()*number;
  return new_weight;
}

double operator+ (double& number, const weight& weight_to_add)
{
  double new_weight = weight_to_add.get_weight()+number;
  return new_weight;
}

double get_weight(const weight &w)
{
 return w.get_weight();
}

bool is_active(const weight &w)
{
 return w.is_active();
}

std::vector<double> convert_to_double
  (const std::vector<weight> &weights)
{
  std::vector<double> double_vector;

  for(size_t i = 0; i != weights.size(); i++)
    {
      double_vector.push_back(weights[i].get_weight());
    }
  return double_vector;
}

#ifndef NDEBUG
void test_weight() noexcept
{
  #define FIX_ISSUE_80
  #ifdef FIX_ISSUE_80
  {
   weight w{1};
   const double& actual_weight = w.get_weight();
   const bool& var_is_active = w.is_active();
   assert(var_is_active && actual_weight);

   //and then maybe make the free functions
   const double& free_actual_weight = get_weight(w);
   const bool& free_is_active = is_active(w);
   assert(free_is_active  && free_actual_weight );

  }
  #endif

  #define FIX_ISSUE_89
  #ifdef FIX_ISSUE_89
  ///The weight of a weight object can be changed (without changing its activation state)
  {
    weight w{1, true};
    weight w2{2, true};
    weight w3{2, false};
    assert (w != w2);

    w.change_weight(2);

    assert (w == w2);
    assert (w != w3);
  }
  #endif

#define FIX_ISSUE_88
#ifdef FIX_ISSUE_88
 /// We can return a vector of double with weights from a vector of weights
{
int vector_length = 10;
std::vector<weight> weights_vector;


for (int i = 0; i != vector_length; ++i){
    weight w{(double)i};
    weights_vector.push_back(w);
  }

std::vector<double> weights_as_double;

weights_as_double = convert_to_double(weights_vector);

for(int i = 0; i != vector_length; i++)
  {
          assert(weights_vector[i].get_weight() == weights_as_double[i]);
  }
}
#endif

#define FIX_ISSUE_93
#ifdef FIX_ISSUE_93
 /// We can return a vector of double with weights from a vector of weights,
 /// with 0 instead of the weight when the connection is inactive
{
int vector_length = 10;
std::vector<weight> weights_vector;


for (int i = 0; i != vector_length; ++i){
    if(i % 2){
      weight w{(double)i, true};
      weights_vector.push_back(w);
      }
    else{
      weight w{(double)i, false};
      weights_vector.push_back(w);
      }
  }

std::vector<double> weights_as_double;

weights_as_double = convert_to_double_or_zero(weights_vector);

for(int i = 0; i != vector_length; i++)
  {
    if(weights_vector[i].is_active() == true)
          assert(weights_vector[i].get_weight() == weights_as_double[i]);
    else
          assert(weights_as_double[i] == 0);
  }
}
#endif

}
#endif
