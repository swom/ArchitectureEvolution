#include "individual.h"

#include <algorithm>
#include <cassert>

#ifndef NDEBUG
void test_individual()
{

  ///An indiivdual is initialized with a network architecture,
  /// by default 1,2,1
  {
    std::vector<int> net_arch{1,2,1};
    net_param n_p{net_arch};
    ind_param i_p{};
    individual i{i_p};
    network n{net_arch};
    assert(i.get_net() == network{n_p});
  }

  ///Individuals have a vector of fixed input values, always equal to 1, for their network
  {
    int n_input = 456;
    std::vector<int> net_arch{n_input};
    ind_param i_p{};
    i_p.net_par.net_arc = net_arch;
    i_p.net_par.max_arc = net_arch;
    individual i{i_p};

    assert(i.get_input_values().size() == static_cast<size_t>(n_input));
    for(const auto& value : i.get_input_values())
      {
        assert(are_equal_with_tolerance(value, 1.0));
      }
  }

  ///When an individual responds to environment it uses its input values as input
  {
        individual i{ind_param{}};
        assert( response(i) == output(i.get_net(),i.get_input_values(), &linear));
  }

//#define FIX_ISSUE_36
  {
    net_param net_par;
    ind_param i_p{net_par};
    individual i{i_p};
    assert(i.get_net() == network{net_par});
  }

#define FIX_ISSUE_124
#ifdef FIX_ISSUE_124
  //individaul constructs mutator_network and assigns it to its network pointer
    //by default templated with mutation_type = mutation_type::weights
  {
    net_param net_par;
    ind_param i_p{net_par};
    individual i{i_p};

    network n {net_par};
    network<mutation_type::activation> mutator_net(net_par);


    //Not same template
    assert(!is_same_mutator_network(i.get_net(), mutator_net));
    //Same template
    assert(is_same_mutator_network(i.get_net(), n));
  }
  #endif
  
#define FIX_ISSUE_120
#ifdef FIX_ISSUE_120
///ind_param contains a mutation_type member, with which network can be templated
  {
    ind_param i_p;
    enum mutation_type mut_type = i_p.m_mutation_type;

    assert(mut_type == mutation_type::activation);
  }
#endif

}
#endif
