#ifndef NODE_H
#define NODE_H
#include <vector>
#include <weight.h>

class node
{
public:
    node(std::vector<weight> vector_weights, bool is_active);

    const bool &is_active() const noexcept {return m_active;}
    const std::vector<weight> &get_vec_weights() const noexcept {return m_weights;}

private:
    bool m_active;
    std::vector<weight> m_weights;
};

void test_node() noexcept;

#endif // NODE_H
