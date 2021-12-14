#include "node.h"

node::node(std::vector<weight> vector_weights, bool is_active):
    m_active{is_active},
    m_weights{vector_weights}
{

}

bool operator== (const node& lhs, const node& rhs)
{
    return lhs.get_vec_weights() == rhs.get_vec_weights() &&
      lhs.is_active() == rhs.is_active();
}

bool operator!= (const node& lhs, const node& rhs)
{
    return !(lhs == rhs);
}

bool node_is_inactive(node node)
{
    return !node.is_active();
}

#ifndef NDEBUG
void test_node() noexcept
{
#define FIX_ISSUE_197
#ifdef FIX_ISSUE_197
  {
     //non-default vector of weights
     std::vector<weight> weights (3, {1, false});

     node test_node{weights, true};

     assert(test_node.is_active());
     assert(test_node.get_vec_weights() == weights);
  }
#endif

//#define FIX_ISSUE_199
#ifdef FIX_ISSUE_199
    ///biases are in the node class now
  {
     std::vector<weight> weights (1, weight{});
     node test_node{weights, true};

     assert(test_node.get_bias() == 0);
  }
#endif

}
#endif
