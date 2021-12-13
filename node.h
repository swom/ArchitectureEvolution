#ifndef NODE_H
#define NODE_H
#include <vector>
#include <weight.h>

class node
{
public:
    node(std::vector<weight> vector_weights, bool is_active = false);

    const bool &is_active() const noexcept {return m_active;}
    const std::vector<weight> &get_vec_weights() const noexcept {return m_weights;}
    void change_nth_weight(weight new_weight, size_t index) {m_weights[index] = new_weight;}

private:
    bool m_active;
    std::vector<weight> m_weights;
};

bool operator== (const node& lhs, const node& rhs);

void test_node() noexcept;

#endif // NODE_H
