#include "environment.h"
#include "utilities.h"
#include <cassert>


environment::environment(double target_valueA, double target_valueB):
  m_ref_target_values{target_valueA,target_valueB},
  m_current_target_value {target_valueA}
{




}

environment::environment(env_param e_p):
  m_ref_target_values{e_p.targetA,e_p.targetB},
  m_current_target_value {e_p.targetA}
{




}



bool operator== (const environment& lhs, const environment& rhs)
{
  bool ref_t_values = lhs.get_ref_target_values() == rhs.get_ref_target_values();
  bool current_t_value = are_equal_with_tolerance(lhs.get_current_target_value(), rhs.get_current_target_value());

  return ref_t_values && current_t_value;
}

double get_target_valueA(const environment& e)
{
  return e.get_ref_target_values()[0];
}

double get_target_valueB(const environment& e)
{
  return e.get_ref_target_values()[1];
}

void switch_target(environment &e){
  //Check which target value is the current one and switch it over to the other

  if (are_equal_with_tolerance(e.get_current_target_value(),
                               get_target_valueA(e))
      )
    {
      e.set_current_target_value(get_target_valueB(e));
    }
  else
    {
      e.set_current_target_value(get_target_valueA(e));
    }
}




#ifndef NDEBUG
void test_environment() noexcept
{

  //an environment has a m_current_target_value member
  {
    double target_valueA = 0.123456;
    double target_valueB = 0.654321;
    environment e{target_valueA, target_valueB};
    assert(e.get_current_target_value() < 0 || e.get_current_target_value() > 0);
  }


  //an env has 2 reference target values;
  {
    double target_valueA = 0.123456;
    double target_valueB = 0.654321;
    environment e{target_valueA, target_valueB};
    assert(e.get_ref_target_values().size() == 2);
  }

  //an env can be initialized with 2 reference target values
  {
    double target_valueA = 0.123456;
    double target_valueB = 0.654321;
    environment e{target_valueA, target_valueB};

    assert(get_target_valueA(e) - target_valueA < 0.0001
           && get_target_valueA(e) - target_valueA > -0.0001);

    assert(get_target_valueB(e) - target_valueB < 0.0001
           && get_target_valueB(e) - target_valueB > -0.0001);
  }



  //Current target value is initialized to the first of the 2 target values
  {
    double targetA = 0.123456;
    double targetB = 0.654321;
    environment e{targetA,targetB};
    assert(are_equal_with_tolerance(e.get_current_target_value(), targetA));
    assert(are_not_equal_with_tolerance(e.get_current_target_value(), targetB));
  }


  //ISSUE_25
  //An environment can switch target values
  {
    double targetA = 0.123456;
    double targetB = 0.654321;
    environment e{targetA,targetB};
    assert(are_equal_with_tolerance(e.get_current_target_value(), targetA));
    switch_target(e);
    assert(are_equal_with_tolerance(e.get_current_target_value(), targetB));
    switch_target(e);
    assert(are_equal_with_tolerance(e.get_current_target_value(), targetA));
  }


  //#define FIX_ISSUE_35

  {
    double targetA = 123456;
    double targetB = 46589;

    env_param e_p{targetA, targetB};
    environment e{e_p};
    assert(are_equal_with_tolerance(get_target_valueA(e), targetA));
    assert(are_equal_with_tolerance(get_target_valueB(e), targetB));

  }


}
#endif
