#ifndef NODE_H
#define NODE_H
#include <vector>
#include "weight.h"

class node
{
public:


    node(std::vector<weight> vector_weights = {}, bool is_active = false);

    NLOHMANN_DEFINE_TYPE_INTRUSIVE(node,
                                   m_active,
                                   m_weights,
                                   m_bias)

    const bool &is_active() const noexcept {return m_active;}

    bool is_inactive() const noexcept {return !m_active;}

    const double &get_bias() const noexcept {return m_bias;}

    ///returns the const reference to the woeight vector
    const std::vector<weight> &get_vec_weights() const noexcept {return m_weights;}

    ///Returns the reference to the weight vector should be used only when
    /// calculating reaction norms with a mutated weight
    std::vector<weight> &get_vec_mutable_weights() noexcept {return m_weights;}

    void change_nth_weight(weight new_weight, size_t index) {m_weights[index] = new_weight;}

    void mutate_bias(const double& mut){m_bias += mut;}

    void reverse_mutate_bias(const double& mut){m_bias -= mut;}

    ///Sets the value of the nodes to the given value
    void change_bias(double new_bias) {m_bias = new_bias;}
    void activate() {m_active = true;}
    void deactivate() {m_active = false;}

private:
    bool m_active;
    std::vector<weight> m_weights;
    double m_bias;
};

bool operator== (const node& lhs, const node& rhs);

bool operator!= (const node& lhs, const node& rhs);

///returns true if the node is active
bool is_active(const node& node);

///returns true if the node is inactive
bool is_inactive(const node& node);

void test_node() noexcept;

#endif // NODE_H
