#include "node.h"

node::node(std::vector<weight> vector_weights, bool is_active):
    m_active{is_active},
    m_weights{vector_weights}
{

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

}
#endif
