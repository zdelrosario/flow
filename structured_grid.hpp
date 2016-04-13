#ifndef GRID // Include guard
#define GRID

#include <vector>
#include <algorithm>
#include <iostream>

/** 2-D Structured Grid class
 * @tparam Val Grid value type
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
  size_type n_,m_,d_; // Grid dimensions
  // TODO -- define these statically
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
  /** Set grid value at index pair
   * 
   * @param i Vertical index
   * @param j Horizontal index
   * @param val Value to assign
   * 
   * @pre 0<=i<=n_-1
   * @pre 0<=j<=m_-1
   */
  value_type set(size_type i, size_type j, value_type val) {
    // Left Boundary
    if (j==0) {
      left_[i] = val;
      return left_[i];
    }
    // Right Boundary
    else if (j==m_-1) {
      right_[i] = val;
      return right_[i];
    }
    // Top Boundary
    else if (i==0) {
      top_[j-1] = val;
      return top_[j-1];
    }
    // Bot Boundary
    else if (i==n_-1) {
      bot_[j-1] = val;
      return bot_[j-1];
    }
    // Interior point
    else {
      v_[(i-1)*(m_-2)+(j-1)] = val;
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
  StructuredGrid(size_type n, size_type m, 
                 std::vector<value_type>& v, size_type d) {
    n_ = n; m_ = m; d_ = d;
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
    bot_.resize(m_-2);
    // Iterate horizontally over 2D grid
    for (size_type i=0; i+2<m_; ++i) {
      top_[i] = *(v_.begin()+i);
      bot_[i] = *(v_.end()-(m_-1)+i);
    }
  }
  // 
  // PUBLIC PROXY OBJECT
  // 
  /** Access return type
   */
  class Value {
    friend class Access;
    friend class StructuredGrid;
    size_type i_,j_;
    StructuredGrid* grid_;
  public:
    Value(size_type i, size_type j, StructuredGrid* grid)
      : i_(i), j_(j), grid_(grid) {}
    value_type operator=(value_type val) {
      return grid_->set(i_,j_,val);
    }
    operator value_type() const { return grid_->value(i_,j_); }
  };
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
    Value operator()(size_type i, size_type j) {
      return Value(i,j,grid_);
    }
  };
  // Return access object
  Access access() {
    return Access(this);
  }
  // 
  // DEBUG METHODS
  // 
  /* Print a 2D array's values */
  template <typename V>
  void print_array(size_type n, size_type m, V v) {
    // row
    for (size_type i=0; i<n; ++i) {
      // col
      for (size_type j=0; j<m-1; ++j) {
        std::cout << v[i*m+j] << " , ";
      }
      std::cout << v[i*m+n] << std::endl;
    }
    return;
  }
  /* Print interior points */
  void printv() {
    // print interior
    print_array(n_-2,m_-2,v_);
  }
  /* Print boundary points */
  void printb() {
    // Print left
    for (auto i: left_)
      std::cout << i << " , ";
    std::cout << std::endl;
    // Print right
    for (auto i: right_)
      std::cout << i << " , ";
    std::cout << std::endl;
    // Print top
    for (auto i: top_)
      std::cout << i << " , ";
    std::cout << std::endl;
    // Print bot
    for (auto i: bot_)
      std::cout << i << " , ";
    std::cout << std::endl;
  }

};

#endif // GRID
