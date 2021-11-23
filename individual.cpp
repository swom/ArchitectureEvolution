#include "individual.h"

#include <algorithm>
#include <cassert>

individual::individual(std::vector<int> net_arch, mutation_type mut) :
  ///!!!!Attention!!!! input values are for now a fixed amount
  m_input_values(net_arch[0], 1.0)
{
 switch (mut)
 {
 case mutation_type::activation :
     std::make_unique<mutator_network<mutation_type::activation>>(net_arch);
     break;
 default:
     throw std::runtime_error("Unkwon mutation type");
 }

}

individual::individual(ind_param i_p) :
  ///!!!!Attention!!!! input values are for now a fixed amount
  m_input_values(i_p.net_par.net_arc[0], 1.0)
{
 switch (i_p.mutation_type)
 {
 case mutation_type::activation :
     std::make_unique<mutator_network<mutation_type::activation>>(i_p.net_par);
     break;
 default:
     throw std::runtime_error("Unkwon mutation type");
 }

}





bool operator== (const individual& lhs, const individual& rhs)
{
  bool fitness = are_equal_with_tolerance(lhs.get_fitness(), rhs.get_fitness());
  bool network = lhs.get_net() == rhs.get_net();
  bool inputs = lhs.get_input_values() == rhs.get_input_values();

  return fitness && network && inputs;
}

double calc_sqr_distance(const individual& i, double env_value)
{
   return (response(i)[0] - env_value) * (response(i)[0] - env_value);
}

void individual::mutate(double mut_rate, double mut_step, std::mt19937_64& rng)
{
  m_network->mutate(mut_rate, mut_step, rng);
}

std::vector<double> response(const individual& ind)
{
    return response(ind.get_net(),ind.get_input_values(), &sigmoid);
}

#ifndef NDEBUG
void test_individual()
{

  ///An indiivdual is initialized with a network architecture,
  /// by default 1,2,1
  {
    std::vector<int> net_arch{1,2,1};
    individual i{net_arch};
    assert(i.get_net() == network{net_arch});
  }

  ///Individuals have a vector of fixed input values, always equal to 1, for their network
  {
    int n_input = 456;
    std::vector<int> net_arch{n_input};
    individual i{net_arch};
    assert(i.get_input_values().size() == static_cast<size_t>(n_input));
    for(const auto& value : i.get_input_values())
      {
        assert(are_equal_with_tolerance(value, 1.0));
      }
  }

  ///When an individual responds to environment it uses its input values as input
  {
    individual i;
    assert( response(i) == response(i.get_net(),i.get_input_values(), &linear));
  }

//#define FIX_ISSUE_36
  {
    net_param net_par;
    ind_param i_p{net_par};
    individual i{i_p};
    assert(i.get_net() == network{net_par});
  }
}
#endif
