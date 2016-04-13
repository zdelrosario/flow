#include "structured_grid.hpp"

#include <iostream>
#include <array>

#define len 2
using value = std::array<int,len>;

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
  size_type n = 5;
  size_type m = 6;

  std::vector<value> v = { {1,1}, {2,2}, {3,3}, {4,4},
                           {5,5}, {6,6}, {7,7}, {8,8},
                           {9,9}, {10,10}, {11,11}, {12,12} };

  // value a_zeros = {0,0};
  // value a_ones  = {1,1};
  // std::vector<value> v( (n-2)*(m-2), a_zeros );
  // size_type i = 2, j = 2;
  // v[i*(m-2)+j] = a_ones;

  StructuredGrid<value> grid(n,m,v,len);
  StructuredGrid<value>::Access val = grid.access();

  // value t = val(2,2);
  // std::cout << "(" << t[0] << "," << t[1] << ")" << std::endl;

  // grid.printv();

  for (size_type i=0; i<n; ++i) {
    for (size_type j=0; j<m; ++j) {
      grid.print_state(val(i,j));
      std::cout << " , ";
    }
    std::cout << std::endl;
  }

  return 0;
}