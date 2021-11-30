#include "individual.h"

#include <algorithm>
#include <cassert>


individual::individual(const ind_param &i_p) :
  ///!!!!Attention!!!! input values are for now a fixed amount
  m_input_values(i_p.net_par.net_arc[0], 1.0)
{
 switch (i_p.mutation_type)
 {
 case mutation_type::activation :
     m_network = std::make_unique<mutator_network<mutation_type::activation>>(i_p.net_par);
     break;
 case mutation_type::weights :
     m_network = std::make_unique<mutator_network<mutation_type::weights>>(i_p.net_par);
     break;
 case mutation_type::weights_and_activation :
     m_network = std::make_unique<mutator_network<mutation_type::weights_and_activation>>(i_p.net_par);
     break;
 default:
     throw std::runtime_error("Unkwon mutation type");
 }

}

individual::individual(const individual& i) noexcept :
    m_fitness{i.get_fitness()},
    m_input_values{i.get_input_values()},
    m_network{std::make_unique<network>(i.get_net())}
{}

bool operator== (const individual& lhs, const individual& rhs)
{
  bool fitness = are_equal_with_tolerance(lhs.get_fitness(), rhs.get_fitness());
  bool network = lhs.get_net() == rhs.get_net();
  bool inputs = lhs.get_input_values() == rhs.get_input_values();

  return fitness && network && inputs;
}

double calc_sqr_distance(const individual& i, double env_value)
{
   auto output = response(i);
   return (output[0] - env_value) * (output[0] - env_value);
}

void individual::change_net(const network& n)
{
    *m_network = n;
}

void individual::mutate(double mut_rate, double mut_step, std::mt19937_64& rng)
{
  m_network->mutate(mut_rate, mut_step, rng);
}

std::vector<double> response(const individual& ind)
{
    return response(ind.get_net(),ind.get_input_values());
}

using json = nlohmann::json;

void to_json(json& j, const individual& ind) {
    j = json{    {"fitness", ind.get_fitness()},
    {"input_values", ind.get_input_values()},
    {"network", ind.get_net()}
};
}

void from_json(const json& j, individual& ind) {
    j.at("fitness").get_to(ind.get_to_fitness());
    j.at("input_values").get_to(ind.get_to_input_values());
    j.at("network").get_to(ind.get_to_net());
}

#ifndef NDEBUG
void test_individual()
{

  ///An indiivdual is initialized with a network architecture,
  /// by default 1,2,1
  {
    std::vector<int> net_arch{1,2,1};
    ind_param i_p{};
    individual i{i_p};
    assert(i.get_net() == network{net_arch});
  }

  ///Individuals have a vector of fixed input values, always equal to 1, for their network
  {
    int n_input = 456;
    std::vector<int> net_arch{n_input};
    ind_param i_p{};
    i_p.net_par.net_arc = net_arch;
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
        assert( response(i) == response(i.get_net(),i.get_input_values(), &linear));
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
  {
    net_param net_par;
    ind_param i_p{net_par};
    individual i{i_p};

    network n {net_par};
    mutator_network<mutation_type::activation> mutator_net(net_par);


    assert(is_same_mutator_network(i.get_net(), mutator_net));
    assert(!is_same_mutator_network(i.get_net(), n));

  #define FIX_ISSUE_125
  #ifdef FIX_ISSUE_125
  ///individual stores a pointer to a base network
    {

      net_param net_par;
      network n{net_par};

      std::unique_ptr<network> n_ptr(new network(net_par));

      ind_param i_p{net_par};
      individual i{i_p};


      assert(i.get_net() == *n_ptr);
      assert(!(i.get_net_ptr() == n_ptr));
    }
#endif

#define FIX_ISSUE_120
#ifdef FIX_ISSUE_120
///ind_param contains a mutation_type member, with which network can be templated
  {
    ind_param i_p;
    enum mutation_type mut_type = i_p.mutation_type;

    assert(mut_type == mutation_type::activation);
  }
#endif

}
#endif
