#include "structured_grid.hpp"

#include <iostream>
#include <tuple>

#define len 2
using value = std::tuple<int,int>;

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
  size_type n = 7;
  size_type m = 7;

  std::vector<value> v( (n-2)*(m-2), std::make_tuple(0,0) );

  size_type i = 2, j = 2;
  v[i*(m-2)+j] = std::make_tuple(1,1);

  StructuredGrid<value> grid(n,m,v);
  std::cout << "flag" << std::endl;
  StructuredGrid<value>::Access val = grid.access();

  auto t = val(0,0);
  std::cout << std::get<0>(v[0]) << std::endl;

  return 0;
}