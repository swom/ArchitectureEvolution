#ifndef NODE_H
#define NODE_H
#include <vector>
#include <weight.h>

class node
{
public:


    node(std::vector<weight> vector_weights = {}, bool is_active = false);

    NLOHMANN_DEFINE_TYPE_INTRUSIVE(node,
                                   m_active,
                                   m_weights,
                                   m_bias)

    const bool &is_active() const noexcept {return m_active;}
    const std::vector<weight> &get_vec_weights() const noexcept {return m_weights;}
    const double &get_bias() const noexcept {return m_bias;}
    void change_nth_weight(weight new_weight, size_t index) {m_weights[index] = new_weight;}

    ///Sets the value of the nodes to the given value
    void change_bias(double new_bias) {m_bias = new_bias;}
    void activate() {m_active = true;}

private:
    bool m_active;
    std::vector<weight> m_weights;
    double m_bias;
};

bool operator== (const node& lhs, const node& rhs);

bool operator!= (const node& lhs, const node& rhs);

///returns true if the node is inactive
bool node_is_inactive(node node);

void test_node() noexcept;

#endif // NODE_H
