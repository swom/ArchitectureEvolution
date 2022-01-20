#ifndef WEIGHT_H
#define WEIGHT_H
#include "json.hpp"
#include <vector>
#include "rndutils.hpp"

class weight
{
public:
  weight(double weight_init = 0, bool is_active = true);

  NLOHMANN_DEFINE_TYPE_INTRUSIVE(weight,
                                 m_weight,
                                 m_is_active);

  ///Returns the weight of a connection
  const double &get_weight() const noexcept {return m_weight;}

  ///Returns a bool indicating whether this connection is active
  const bool &is_active() const noexcept {return m_is_active;}

  void change_weight(double new_weight) {m_weight = new_weight;}

  void change_activation(double new_activation) {m_is_active = new_activation;}

private:
  double m_weight;
  bool m_is_active;
};

bool operator== (const weight& lhs, const weight& rhs);
bool operator!= (const weight& lhs, const weight& rhs);

//I am not sure when this would be needed but some files
//included by json.hpp would throw a fit when i didn't have it
double operator* (double& number, const weight& weight_to_multiply);
double operator+ (double& number, const weight& weight_to_add);


///Free function that returns weight
double get_weight(const weight &w);

///Free function that whether a connection is active
bool is_active(const weight &w);

///Converts a vector of weights to a vector of double with the weights
std::vector<double> convert_to_double(const std::vector<weight> &weights);

///Converts a vector of weights to a vector of double with the weights and weights of inactive connections replaced with 0
std::vector<double> convert_to_double_or_zero(const std::vector<weight> &weights);

void test_weight() noexcept;

#endif // WEIGHT_H
