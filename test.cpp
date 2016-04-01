#include "structured_grid.hpp"

#include <iostream>

using value = int;

typedef StructuredGrid<value>::size_type size_type;

// DEBUG -- Print an array
template <typename V>
void print_array(size_type n, size_type m, V v) {
  // row
  for (size_type i=0; i<n; ++i) {
    // col
    for (size_type j=0; j<m-1; ++j) {
      std::cout << v[i*m+j] << ",";
    }
    std::cout << v[i*m+n] << std::endl;
  }
  return;
}

int main() {
  size_type n = 4;
  size_type m = 5;
  std::vector<value>
            v = { 1, 2, 3,
                  4, 5, 6 };

  StructuredGrid<value> grid(n,m,v);
  StructuredGrid<value>::Access val = grid.access();

  // std::cout << val(1,1) << std::endl;

  // DEBUG -- print out the grid
  for (size_type i=0; i<n; ++i) {
    for (size_type j=0; j<m-1; ++j) {
      std::cout << val(i,j) << ",";
    }
    std::cout << val(i,m-1) << std::endl;
  }

  return 0;
}