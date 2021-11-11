#ifndef WEIGHT_H
#define WEIGHT_H


class weight
{
public:
  weight(double weight_init = 1, bool is_active = true);

  ///Returns the weight of a connection
  const double &get_weight() const noexcept {return m_weight;}

  ///Returns a bool indicating whether this connection is active
  const bool &is_active() const noexcept {return m_is_active;}

private:
  double m_weight;
  bool m_is_active;
};



///Free function that returns weight
double get_weight(const weight &w);

///Free function that whether a connection is active
bool is_active(const weight &w);

void test_weight() noexcept;

#endif // WEIGHT_H
