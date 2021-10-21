#include <iostream>
#include <iterator>
#include "deque.h"

int main() {
    Deque<int> d;
    d.push_back(1);
    auto it = d.begin();
    std::advance(it, 1);
    return 0;
}
