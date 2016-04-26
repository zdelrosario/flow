#ifndef GRID // Include guard
#define GRID

#include <vector>     // data handling
#include <algorithm>  // std::copy
#include <iostream>   // debug
#include <cmath>      // floor
#include <fstream>    // file handling

#include "Util.hpp"   // totally_ordered

/** 2-D Structured Grid class
 * @tparam Val Grid value type
 */
template <typename X, typename Val>
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
  // Grid points
  std::vector<X> x_; // Grid points
  // Grid cell values
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
// std::cout << "grid_value (" << i << "," << j << ")" << std::endl;
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
  /** Vector file writeout
   * @brief Writes a set of vector points 
   *        to a formatted data file
   *
   * @tparam V Array-like container, n-Dimensional
   * 
   * @param outputfile String which defines output filename
   * @param v Vector which defines gridpoints
   * 
   * @post A file with name 'outputfile' is written to the
   *       local directory with the gridpoints
   */
  template <typename V>
  void writeout(std::string outputfile, V& v) {
    std::ofstream f_out(outputfile.c_str());
    // f_out.precision(5);
    // Write out elements
    for (auto it=v.begin(); it!=v.end(); ++it) {
      // f_out << (*it)[0] << "," << (*it)[1] << std::endl;
      for (auto jt=it->begin(); (jt+1)!=it->end(); ++jt)
        f_out << (*jt) << ",";
      f_out << (*it).back() << std::endl;
    }
  }
public:
  // 
  // PUBLIC MEMBER FUNCTIONS
  // 
  /** Public constructor
   *  Constructs a cell grid of given dimensions and initial conditions
   * 
   * @param n Number of total vertical cells
   * @param m Number of total horizontal cells
   * @param v Initial conditions of cells
   * @param x Vector of physical cell corner points
   * 
   * @pre (n-2)*(m-2) == v.size(), number of interior cells
   * @pre (n-1)*(m-1) == x.size(), number of cell corner points
   *
   * TODO -- pass boundary conditions (handle here?)
   * TODO -- pass curvilinear mapping
   */
  StructuredGrid(size_type n, size_type m, 
                 std::vector<value_type>& v,
                 std::vector<X>& x) {
    n_ = n; m_ = m;
    v_.resize(v.size());
    std::copy(v.begin(),v.end(),v_.begin());
    x_.resize(x.size());
    std::copy(x.begin(),x.end(),x_.begin());

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
      bot_[i] = *(v_.end()-(m_-2)+i);
    }
  }
  // 
  // PUBLIC PROXY OBJECTS
  // 
  /** Access return type
   */
  class Value {
    friend class Access;
    friend class StructuredGrid;
    size_type i_,j_;
    StructuredGrid* grid_;
  public:
    /* Invalid Constructor */
    Value() {}
    /* Public Constructor */
    Value(size_type i, size_type j, StructuredGrid* grid)
      : i_(i), j_(j), grid_(grid) {}
    /* Value Assignment Operator */
    Value operator=(Value val) {
      i_ = val.i_;
      j_ = val.j_;
      grid_ = val.grid_;
    }
    /* Forwarding Assignment Operator */
    value_type operator=(value_type val) {
      return grid_->set(i_,j_,val);
    }
    /* Implicit Conversion Operator */
    operator value_type() const { 
      return grid_->value(i_,j_); 
    }
    /* Const Subscript Operator */
    // TODO -- define type as template
    float operator[](size_type ind) const {
      return grid_->value(i_,j_)[ind];
    }
    // DEBUG -- Print value to console
    void print() {
      grid_->print_state(grid_->value(i_,j_));
    }
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
  /** @class StructuredGrid::Cell
   *  @brief Proxy object for a grid cell
   * 
   * @param grid Parent grid
   * @param idx  Cell id number
   * 
   * @pre 0 <= idx <= v_.size()
   */
  class Cell {
    // Allow StructuredGrid to access this object
    friend class StructuredGrid;
    // Pointer back to grid
    StructuredGrid* grid_;
    // Cell index
    size_type i_,j_;
    // Private constructor
    Cell(const StructuredGrid* grid, size_type i, size_type j)
        : grid_(const_cast<StructuredGrid*>(grid)), i_(i), j_(j) {}
  public:
    // Public state vector type
    typedef value_type CellValue;
    // Public Member functions
    size_type iy() {
      return i_;
    }
    size_type jx() {
      return j_;
    }
    Value value() {
      return Value(i_,j_,grid_);
    }
    /* Interior cell index */
    size_type idx() const {
      return (i_-1)*(grid_->m_-2)+j_-1;
    }
    // Returns value of cell at relative index
    Value value(int di, int dj) {
      // Compute indices
      di = di+i_;
      dj = dj+j_;
      // Bound indices
      if (di<0) di=0;
      else if (di>int(grid_->n_-1)) di=int(grid_->n_-1);
      if (dj<0) dj=0;
      else if (dj>int(grid_->m_-1)) dj=int(grid_->m_-1);
      // Access value
      return Value(di,dj,grid_);
    }
    // Returns neighbor cell at relative index
    Cell neighbor(int di, int dj) {
      // Compute indices
      di = di+i_;
      dj = dj+j_;
      // Bound indices
      if (di<0) di=0;
      else if (di>int(grid_->n_-1)) di=int(grid_->n_-1);
      if (dj<0) dj=0;
      else if (dj>int(grid_->m_-1)) dj=int(grid_->m_-1);
      // Access cell
      return grid_->cell(di,dj);
    }
    /** Returns cell corner point
     * @param n Corner point index
     * @pre 1<=n<=4
     */
    X x(size_type n) {
      if (n==1) {
        return grid_->x_[(i_-1)*(grid_->m_-1)+(j_-1)];
      }
      else if (n==2) {
        return grid_->x_[(i_-1)*(grid_->m_-1)+j_];
      }
      else if (n==3) {
        return grid_->x_[i_*(grid_->m_-1)+j_];
      }
      else if (n==4) {
        return grid_->x_[i_*(grid_->m_-1)+j_-1];
      }
      else {
        assert(0);
      }
    }
    /* Returns average cell width */
    float dx() {
      return (this->x(2)[0]+this->x(3)[0]-this->x(1)[0]-this->x(4)[0])/2.0;
    }
    /* Returns average cell height */
    float dy() {
      return (this->x(1)[1]+this->x(2)[1]-this->x(3)[1]-this->x(4)[1])/2.0;
    }
   };

   // Return cell object by interior id number
   Cell cell(size_type idx) {
    return Cell(this,floor(idx/(m_-2))+1,idx%(m_-2)+1);
   }
   /** Return cell by id pair
    * @pre 0 <= i <= n_-1
    * @pre 0 <= j <= m_-1
    */
   Cell cell(size_type i, size_type j) {
    return Cell(this,i,j);
   }
  //
  // CELL ITERATOR
  //
  class cell_iterator : private totally_ordered<cell_iterator> {
  private:
    StructuredGrid* grid_;
    size_type idx_;
  public:
    // Iterator traits (magic)
    typedef Cell                      value_type;
    typedef Cell*                     pointer;
    typedef Cell&                     reference;
    typedef std::ptrdiff_t            difference_type;
    typedef std::forward_iterator_tag iterator_category;
    // Public constructor
    cell_iterator(const StructuredGrid* grid, 
                  size_type idx) :
      grid_(const_cast<StructuredGrid*>(grid)), idx_(idx) {}
    // Implement minimal methods
    Cell operator*() const {
      return grid_->cell(idx_);
    }
    cell_iterator& operator++() { ++idx_; return *this; }
    bool operator==(const cell_iterator& iter) const {
      return idx_ == iter.idx_;
    }
  };
  // Return a cell_iterator
  cell_iterator cell_begin() const {
    return cell_iterator(this,0);
  }
  cell_iterator cell_end() const {
    return cell_iterator(this,v_.size());
  }

  // 
  // FILE HANDLING METHODS
  // 
  /* Writes the current cell values to output */
  void write_values(std::string outputfile) {
    writeout(outputfile, v_);
  }
  /* Writes the current cell values to output */
  void write_grid(std::string outputfile) {
    writeout(outputfile, x_);
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
        print_state(v[i*m+j]);
        std::cout << " , ";
      }
      print_state(v[i*m+n]);
      std::cout << std::endl;
    }
    return;
  }
  /* Print a state vector */
  void print_state(value_type v) {
    std::cout << "(";
    for (size_type i=0; i+1<v.size(); ++i) {
      std::cout << v[i] << ",";
    }
    std::cout << v.back() << ")";
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
      print_state(i);
      std::cout << " , ";
    std::cout << std::endl;
    // Print right
    for (auto i: right_)
      print_state(i);
      std::cout << " , ";
    std::cout << std::endl;
    // Print top
    for (auto i: top_)
      print_state(i);
      std::cout << " , ";
    std::cout << std::endl;
    // Print bot
    for (auto i: bot_)
      print_state(i);
      std::cout << " , ";
    std::cout << std::endl;
  }

};

#endif // GRID
