#include "weight.h"

#include <cassert>

weight::weight(double weight_init, bool is_active):
  m_weight{weight_init},
  m_is_active{is_active}
{

}

double get_weight(const weight &w)
{
 return w.get_weight();
}

bool is_active(const weight &w)
{
 return w.is_active();
}

#ifndef NDEBUG
void test_weight() noexcept
{
  #define FIX_ISSUE_80
  #ifdef FIX_ISSUE_80
  {
   weight w;
   const double& actual_weight = w.get_weight();
   const bool& var_is_active = w.is_active();
   assert(var_is_active && actual_weight);

   //and then maybe make the free functions
   const double& free_actual_weight = get_weight(w);
   const bool& free_is_active = is_active(w);
   assert(free_is_active  && free_actual_weight );

  }
  #endif

}
#endif
