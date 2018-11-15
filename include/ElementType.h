
#ifndef FDTD_ELEMENTTYPE_H_
#define FDTD_ELEMENTTYPE_H_

#include <string>
#include <map>

enum class ElementType {
    ConstantScalar,
    EdgeE,
    EdgeH,
    NodeScalar,
    NodeVector
};


static std::map<std::string, ElementType> stringToElementTypeMap{
        {"EdgeE", ElementType::EdgeE},
        {"EdgeH", ElementType::EdgeH},
        {"NodeScalar", ElementType::NodeScalar},
        {"NodeVector", ElementType::NodeVector}
};

static std::map<std::string, int> stringDirectionToIntDirectionMap{
        {"x", 0},
        {"y", 1},
        {"z", 2}
};

#endif // FDTD_ELEMENTTYPE_H_
