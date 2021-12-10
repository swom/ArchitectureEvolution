#ifndef NODE_H
#define NODE_H
#include <vector>
#include <weight.h>

class node
{
public:
    node();
    bool m_active;
    std::vector<weight> m_weights;
};

#endif // NODE_H
