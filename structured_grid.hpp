#ifndef GRID // Include guard
#define GRID

#include <vector>
#include <algorithm>
#include <iostream>

// #include <boost/range/adaptor/strided.hpp>
// #include <boost/range/algorithm/copy.hpp>
// using namespace boost::adaptors;
// using namespace boost::assign;

/** 2-D Structured Grid class
 * TODO -- template on grid data type (for state vector)
 */
template <typename Val>
class StructuredGrid {
  //
  // PUBLIC TYPE DEFINITIONS
  //
public:
  typedef unsigned size_type;
  typedef StructuredGrid grid_type;
  typedef Val value_type;
private:
  //
  // PRIVATE DATA MEMBERS
  //
  size_type n_,m_; // Grid dimensions
  // TODO -- define this statically
  std::vector<value_type> v_; // Grid data
  std::vector<value_type> left_;  // Left boundary, size()==n_
  std::vector<value_type> right_; // Right boundary, size()==n_
  std::vector<value_type> top_;   // Top boundary, size()==m_-2
  std::vector<value_type> bot_;   // Bot boundary, size()==m_-2
  // 
  // PRIVATE HELPER FUNCTIONS
  // 
  /** Return grid value at index pair
   * 
   * @param i Vertical index
   * @param j Horizontal index
   * 
   * @pre 0<=i<=n_-1
   * @pre 0<=j<=m_-1
   * 
   * TODO -- Handle different boundary conditions,
   *         currently assumes Dirichlet
   */
  value_type value(size_type i, size_type j) {
    // Left Boundary
    if (j==0) {
      return left_[i];
    }
    // Right Boundary
    else if (j==m_-1) {
      return right_[i];
    }
    // Top Boundary
    else if (i==0) {
      return top_[j-1];
    }
    // Bot Boundary
    else if (i==n_-1) {
      return bot_[j-1];
    }
    // Interior point
    else {
      return v_[(i-1)*(m_-2)+(j-1)];
    }
  }
public:
  // 
  // PUBLIC MEMBER FUNCTIONS
  // 
  /** Public constructor
   *  Constructs a grid of given dimensions and initial conditions
   * 
   * @param n Height of grid
   * @param m Width of grid
   * @param v Initial conditions of grid
   * 
   * @pre (n-2)*(m-2) == v.size()
   *
   * TODO -- pass boundary conditions
   * TODO -- pass curvilinear mapping
   */
  StructuredGrid(size_type n, size_type m, std::vector<value_type>& v) {
    n_ = n;
    m_ = m;
    v_.resize(v.size());
    std::copy(v.begin(),v.end(),v_.begin());

    // TODO -- allow for non-Dirichlet BC's

    left_.resize(n_);
    right_.resize(n_);
    // TODO -- replace loops with striding iterator
    for (size_type i=1; i<(n_-1); ++i) {
      left_[i]  = v_[(i-1)*(m_-2)+0];
      right_[i] = v_[(i-1)*(m_-2)+m_-3];
    }
    // DEBUG -- repeat top and bottom elements
    left_[0] = left_[1]; left_[n_-1] = left_[n_-2];
    right_[0] = right_[1]; right_[n_-1] = right_[n_-2];

    top_.resize(m_-2);
    std::copy(v_.begin(),v_.begin()+n_-1,top_.begin());
    bot_.resize(m_-2);
    std::copy(v_.end()-(n_-2)-1,v_.end(),bot_.begin());
  }
  // 
  // PUBLIC PROXY OBJECT
  // 
  /** @class StructuredGrid::Access
   *  @brief Proxy object for accessing grid data
   */
  class Access {
    // Allow StructuredGrid to access this object
    friend class StructuredGrid;
    // Pointer back to grid
    StructuredGrid* grid_;
    // Private constructor
    Access(const StructuredGrid* grid)
        : grid_(const_cast<StructuredGrid*>(grid)) {}
    // Public Member functions
  public:
    value_type operator()(size_type i, size_type j) {
      return grid_->value(i,j);
    }
  };
  // Return access object
  Access access() {
    return Access(this);
  }
  /** Access return type
   */
  // class Value {

  // };
  // 
  // DEBUG METHODS
  // 
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
  // Print interior and boundary points
  void printv() {
    // print interior
    print_array(n_-2,m_-2,v_);
    // Print left
    for (auto i: left_)
      std::cout << i << ',';
    std::cout << std::endl;
    // Print right
    for (auto i: right_)
      std::cout << i << ',';
    std::cout << std::endl;
    // Print top
    for (auto i: top_)
      std::cout << i << ',';
    std::cout << std::endl;
    // Print bot
    for (auto i: bot_)
      std::cout << i << ',';
    std::cout << std::endl;
  }

};

#endif // GRID
