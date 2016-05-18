#ifndef GRID // Include guard
#define GRID

/** Structured Grid code
 *  TODO:
 *    - Add cell normals
 *    - Add BC's: inflow, outflow, wall, symmetry
 */

#include <vector>     // data handling
#include <valarray>   // std::valarray, std::begin, std::end
#include <array>      // std::array
#include <algorithm>  // std::copy
#include <iostream>   // debug
#include <cmath>      // floor, sqrt
#include <fstream>    // file handling

#include "Util.hpp"   // totally_ordered
#include "gas_dynamics.hpp" // 

#define stages_ 5 // # of RK stages, plus one

/** 2-D Structured Grid class
 * @tparam Val Grid value type
 */
template <typename S, typename X, typename Val, typename Flag>
class StructuredGrid {
  //
  // PUBLIC TYPE DEFINITIONS
  //
public:
  typedef unsigned size_type;
  typedef StructuredGrid grid_type;
  typedef Val value_type;
  typedef Flag flag_type;
  typedef S scalar_type;
private:
  //
  // PRIVATE DATA MEMBERS
  //
  size_type n_,m_; // Grid dimensions
  // Grid points
  std::vector<X> x_; // Grid points
  // Grid cell values
  std::vector<value_type> v_;     // Grid data
  std::array<std::vector<value_type>,stages_> k_; // Space for RK stages
  std::vector<value_type> left_;  // Left boundary, size()==n_
  std::vector<value_type> right_; // Right boundary, size()==n_
  std::vector<value_type> top_;   // Top boundary, size()==m_-2
  std::vector<value_type> bot_;   // Bot boundary, size()==m_-2
  // Boundary condition flags
  std::vector<flag_type> left_b_;  // Left boundary flag, size()==n_
  std::vector<flag_type> right_b_; // Right boundary flag, size()==n_
  std::vector<flag_type> top_b_;   // Top boundary flag, size()==m_-2
  std::vector<flag_type> bot_b_;   // Bot boundary flag, size()==m_-2
  // Boundary cell normals and paralleles
  std::vector<value_type> left_n_;  // Left boundary [normal,parallel]
  std::vector<value_type> right_n_; // Right boundary [normal,parallel]
  std::vector<value_type> top_n_;   // Top boundary [normal,parallel]
  std::vector<value_type> bot_n_;   // Bot boundary [normal,parallel]
  // Freestream properties
  // DEBUG -- fix freestream conditions
  scalar_type rho_inf = 1.1462;
  scalar_type p_inf = 9.725e4;  // same as exit???
  scalar_type u_inf = 68.93;    // 
  scalar_type v_inf = 0.0;      // 
  scalar_type p_e = 9.725e4;    // 
  // 
  // PRIVATE HELPER FUNCTIONS
  // 
  /** Project velocity
   *
   * @param v Input state vector
   * @param n Normal against which to project
   * 
   * @return u_n Velocity projected on normal
   */
  scalar_type proj_vel(const value_type& v, const value_type& n) {
    return uf(v)*n[0] + vf(v)*n[1];
  }
  /** Boundary condition helper function
   * @brief Switches between neumann and dirichlet boundary
   *        condition based on provided flag. Writes to
   *        provided output vector.
   *
   * @param f   Boundary condition vector flag
   * @param v_o Outer value
   * @param v_i Inner value
   * 
   * Boundary conditions are specified by a unsigned vector
   * which specifies the condition per- state vector element.
   *    0 = dirichlet
   *    1 = neumann
   *    2 = mirror
   * The following will apply the BC type to the
   * entire state vector
   *    3 = inflow
   *    4 = outflow
   * 
   * @return State vector obeying boundary conditions
   */
  value_type bc_helper( const flag_type& f, 
                        const value_type& v_o, 
                        const value_type& n,
                        const value_type& v_i ) {
    // Reserve some space
    value_type res; res.resize(f.size());
    
    // DEBUG -- print flag
    for (auto it=f.begin(); it!=f.end(); ++it) {
      std::cout<<*it<<",";
    }
    std::cout<<std::endl;

    // Characteristic (Full-state) BC
    unsigned bc_type = 0;
    for (size_type i=0; i!=f.size(); ++i) {
      if ( (f[i]==3) or (f[i]==4) )
        bc_type = f[i];
    }
    // Apply characteristic BC
    if (bc_type!=0) {
      // Project velocity
      scalar_type U_i,U_o;
      U_i = proj_vel(v_i,n);
      U_o = proj_vel(v_o,n);
      // Compute Riemann invariants
      scalar_type C_i = cf(v_i);
      scalar_type C_o = cf(v_o);
      scalar_type R_i = U_i+2.0/(gam-1)*C_i;
      scalar_type R_o = U_o-2.0/(gam-1)*C_o;
      scalar_type C_s = 0.25*(gam-1)*(R_i-R_o);
      scalar_type U_s = 0.5*(R_i+R_o);

      scalar_type s_inv, U_r, V_r;
      // Inflow
      if (bc_type==3) {
        s_inv = pow(v_o[0],gam)/pf(v_o);
        U_r = uf(v_o) + (U_s-U_o)*n[0];
        V_r = vf(v_o) + (U_s-U_o)*n[1];
      }
      // Outflow
      else if (bc_type==4) {
        s_inv = pow(v_i[0],gam)/pf(v_i);
        U_r = uf(v_i) + (U_s-U_i)*n[0];
        V_r = vf(v_i) + (U_s-U_i)*n[1];
      }
      else {assert(false);}
      // Use invariants to compute ghost state
      res[0] = pow( s_inv*pow(C_s,2)/gam, 1/(gam-1) );
      scalar_type P_r = res[0]*pow(C_s,2)/gam;
      res[1] = res[0]*U_r;
      res[2] = res[0]*V_r;
      res[3] = P_r/(gam-1) + 0.5*res[0]*(pow(U_r,2)+pow(V_r,2));
      return res;
    }

    // Fall through -- per-element BC
    // Choose element based on flag value
    for (size_type i=0; i!=f.size(); ++i) {
      if (f[i] == 1)      // 1 == neumann
        res[i] = v_i[i];
      else if (f[i] == 2) // 2 == mirror
        res[i] = -v_i[i];
      else                // 0 == dirichlet
        res[i] = v_o[i];
    }
    // Return the result
    return res;
  }

  /** Return cell value at index pair
   * 
   * @param i Vertical index
   * @param j Horizontal index
   * 
   * @pre 0<=i<=n_-1, interior points are strict inequality
   * @pre 0<=j<=m_-1
   * 
   * TODO -- Handle different boundary conditions,
   *         currently assumes Dirichlet
   */
  value_type value(size_type i, size_type j) {

    // DEBUG bound assertions
    assert(i<n_);
    assert(j<m_);
    // Left Boundary
    if (j==0) {
      // Top corner
      if (i==0)
        return bc_helper(left_b_[i],left_[i],left_n_[i],
                         value(i+1,j+1));
      // Bottom corner
      else if (i==n_-1)
        return bc_helper(left_b_[i],left_[i],left_n_[i],
                         value(i-1,j+1));
      // Side wall
      else
        return bc_helper(left_b_[i],left_[i],left_n_[i],
                         value(i,j+1));
    }
    // Right Boundary
    else if (j==m_-1) {
      // Top corner
      if (i==0)
        return bc_helper(right_b_[i],right_[i],right_n_[i],
                         value(i+1,j-1));
      // Bottom corner
      else if (i==n_-1)
        return bc_helper(right_b_[i],right_[i],right_n_[i],
                         value(i-1,j-1));
      // Side wall
      else
        return bc_helper(right_b_[i],right_[i],right_n_[i],
                         value(i,j-1));
    }
    // Top Boundary
    else if (i==0) {
      // return top_[j-1];
      return bc_helper(top_b_[j],top_[j],top_n_[j],
                       value(i+1,j));
    }
    // Bot Boundary
    else if (i==n_-1) {
      // return bot_[j-1];
      return bc_helper(bot_b_[j],bot_[j],bot_n_[j],
                       value(i-1,j));
    }
    // Interior point
    else {
      return v_[(i-1)*(m_-2)+(j-1)];
    }
  }
  /** Set cell value at index pair
   * 
   * Note that this function ignores any 
   * sort of boundary conditions.
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
  /* --- RK Stage Helpers --- */
  /** Return RK4 stage cell value at index pair
   * 
   * @param i Vertical index
   * @param j Horizontal index
   * @param ind Stage index
   * 
   * @pre 0<=i<=n_-1, interior points are strict inequality
   * @pre 0<=j<=m_-1
   * @pre 0<=ind<stages_
   */
  value_type stage_value(size_type i, size_type j, size_type ind) {

    // DEBUG bound assertions
    assert(i<n_);
    assert(j<m_);
    // Left Boundary
    if (j==0) {
      // Top corner
      if (i==0)
        return bc_helper(left_b_[i],left_[i],left_n_[i],
                         stage_value(i+1,j+1,ind));
      // Bottom corner
      else if (i==n_-1)
        return bc_helper(left_b_[i],left_[i],left_n_[i],
                         stage_value(i-1,j+1,ind));
      // Side wall
      else
        return bc_helper(left_b_[i],left_[i],left_n_[i],
                         stage_value(i,j+1,ind));
    }
    // Right Boundary
    else if (j==m_-1) {
      // Top corner
      if (i==0)
        return bc_helper(right_b_[i],right_[i],right_n_[i],
                         stage_value(i+1,j-1,ind));
      // Bottom corner
      else if (i==n_-1)
        return bc_helper(right_b_[i],right_[i],right_n_[i],
                         stage_value(i-1,j-1,ind));
      // Side wall
      else
        return bc_helper(right_b_[i],right_[i],right_n_[i],
                         stage_value(i,j-1,ind));
    }
    // Top Boundary
    else if (i==0) {
      // return top_[j-1];
      return bc_helper(top_b_[j],top_[j],top_n_[j],
                       stage_value(i+1,j,ind));
    }
    // Bot Boundary
    else if (i==n_-1) {
      // return bot_[j-1];
      return bc_helper(bot_b_[j],bot_[j],bot_n_[j],
                       stage_value(i-1,j,ind));
    }
    // Interior point
    else {
      return k_[ind][(i-1)*(m_-2)+(j-1)];
    }
  }
  /** Set RK stage cell value at index pair
   * 
   * Note that this function ignores any 
   * sort of boundary conditions.
   * 
   * @param i Vertical index
   * @param j Horizontal index
   * @param ind RK stage index
   * @param val Value to assign
   * 
   * @pre 0<=i<=n_-1
   * @pre 0<=j<=m_-1
   * @pre 0<=ind<stages_
   */
  value_type stage_set(size_type i, size_type j, size_type ind, value_type val) {
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
      k_[ind][(i-1)*(m_-2)+(j-1)] = val;
      return k_[ind][(i-1)*(m_-2)+(j-1)];
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
    for (auto it=std::begin(v); it!=std::end(v); ++it) {
      for (auto jt=std::begin(*it); (jt+1)!=std::end(*it); ++jt)
        f_out << (*jt) << ",";
      f_out << (*it)[it->size()-1] << std::endl;
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
   * Boundary conditions are specified by a unsigned vector
   * which specifies the condition per- state vector element.
   *    0 = dirichlet
   *    1 = neumann
   *    2 = mirror
   * The following will apply the BC type to the
   * entire state vector
   *    3 = inflow
   *    4 = outflow
   * 
   * @param left_b  Left boundary condition
   * @param right_b Right boundary condition
   * @param top_b   Top boundary condition
   * @param bot_b   Bottom boundary condition
   * 
   * @pre (n)*(m) == v.size(), number of total cells
   * @pre (n-1)*(m-1) == x.size(), number of cell corner points
   *
   * TODO -- pass curvilinear mapping
   */
  StructuredGrid(size_type n, size_type m, 
                 std::vector<value_type>& v,
                 std::vector<X>& x,
                 std::vector<flag_type>& left_b,
                 std::vector<flag_type>& right_b,
                 std::vector<flag_type>& top_b,
                 std::vector<flag_type>& bot_b) {
    // Copy grid dimensions
    n_ = n; m_ = m;
    // Copy interior cell values
    v_.resize((n_-2)*(m_-2));
    for (size_type i=0; i<n_-2; ++i) {
      for (size_type j=0; j<m_-2; ++j) {
        v_[i*(m_-2)+j] = v[(i+1)*m_+(j+1)];
      }
    }
    // Resize each RK stage container
    for (size_type ind=0; ind<stages_; ++ind) {
      k_[ind].resize((n_-2)*(m_-2));
    }
    // Copy grid point values
    x_.resize(x.size());
    std::copy(x.begin(),x.end(),x_.begin());
    // Copy boundary condition flags
    left_b_.resize(left_b.size());
    std::copy(left_b.begin(),left_b.end(),left_b_.begin());
    right_b_.resize(right_b.size());
    std::copy(right_b.begin(),right_b.end(),right_b_.begin());
    top_b_.resize(top_b.size());
    std::copy(top_b.begin(),top_b.end(),top_b_.begin());
    bot_b_.resize(bot_b.size());
    std::copy(bot_b.begin(),bot_b.end(),bot_b_.begin());

    // Generate boundary values
    left_.resize(n_);
    right_.resize(n_);
    // TODO -- replace loops with striding iterator
    for (size_type i=0; i<n_; ++i) {
      left_[i]  = v[i*m_];
      right_[i] = v[i*m_+m_-1];
    }

    top_.resize(m_-2);
    bot_.resize(m_-2);
    // Iterate horizontally over top and bottom,
    // skip the left and right edges
    for (size_type i=0; i+2<m_; ++i) {
      top_[i] = *(v.begin()+i+1);
      bot_[i] = *(v.end()-m_+i+1);
    }

    // DEBUG -- fix boundary normals for a box
    left_n_.resize(left_b_.size());
    right_n_.resize(right_b_.size());
    top_n_.resize(top_b_.size());
    bot_n_.resize(bot_b_.size());

    // 
    std::fill(left_n_.begin(),left_n_.end(),  
              value_type({-1.0,+0.0,+0.0,+1.0}));
    std::fill(right_n_.begin(),right_n_.end(),
              value_type({+1.0,+0.0,+0.0,-1.0}));
    std::fill(top_n_.begin(),top_n_.end(),    
              value_type({+0.0,+1.0,+1.0,+0.0}));
    std::fill(bot_n_.begin(),bot_n_.end(),    
              value_type({+0.0,-1.0,-1.0,+0.0}));
  }
  // 
  // STAGE HANDLING FUNCTIONS
  //
  /* Fill RK stage containers
   * @brief Fills all RK stage containers with
   *        a provided state vector value
   * @param val State vector fill value
   */
  void fill_stages(value_type val) {
    for (size_type ind=0; ind<stages_; ++ind) {
      std::fill(k_[ind].begin(),k_[ind].end(),val);
    }
  }

  // 
  // PUBLIC PROXY OBJECTS
  // 
  /** @class StructuredGrid::Proxy
   * @brief Provides read and write access to grid
   *        cell values through proxy object
   */
  class Proxy {
    friend class Access;
    friend class StructuredGrid;
    size_type i_,j_;
    StructuredGrid* grid_;
    short ind_;
  public:
    /* Invalid Constructor */
    Proxy() {}
    /* Public Constructor */
    Proxy(size_type i, size_type j, StructuredGrid* grid, short ind=-1)
      : i_(i), j_(j), grid_(grid), ind_(ind) {}
    /* Value Assignment Operator */
    Proxy operator=(Proxy val) {
      i_ = val.i_;
      j_ = val.j_;
      grid_ = val.grid_;
      ind_ = val.ind_;
    }
    /* Forwarding Assignment Operator */
    value_type operator=(value_type val) {
      if (ind_==-1)
        return grid_->set(i_,j_,val);
      else
        return grid_->stage_set(i_,j_,ind_,val);
    }
    /* Implicit Conversion Operator */
    operator value_type() const { 
      if (ind_==-1)
        return grid_->value(i_,j_); 
      else
        return grid_->stage_value(i_,j_,ind_); 
    }
    /* Const Subscript Operator */
    scalar_type operator[](size_type ind) const {
      if (ind_==-1)
        return grid_->value(i_,j_)[ind];
      else
        return grid_->stage_value(i_,j_,ind_)[ind];
    }
    /* Size of referenced container */
    size_type size() {
      if (ind_==-1)
        return grid_->value(i_,j_).size();
      else
        return grid_->stage_value(i_,j_,ind_).size();
    }
    // DEBUG -- Print value to console
    void print() {
      if (ind_==-1)
        grid_->print_state(grid_->value(i_,j_));
      else
        grid_->print_state(grid_->stage_value(i_,j_,ind_));
    }
  };
  /** @class StructuredGrid::Cell
   *  @brief Proxy object for a gkrid cell
   * 
   * @param grid Parent grid
   * @param idx  Cell id number
   * @param ind  RK stage number, -1 for solution values
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
    short ind_;
    // Private constructor
    Cell(const StructuredGrid* grid, size_type i, size_type j, short ind=-1)
        : grid_(const_cast<StructuredGrid*>(grid)), i_(i), j_(j), ind_(ind) {}
  public:
    // Public types
    typedef value_type CellValue;
    typedef scalar_type CellScalar;
    // Forward assignment operator
    CellValue operator=(CellValue val) {
      return Proxy(i_,j_,grid_,ind_)=val;
    }
    // Public Member functions
    size_type iy() {
      return i_;
    }
    size_type jx() {
      return j_;
    }
    Proxy value() {
      return Proxy(i_,j_,grid_,ind_);
    }
    /* Interior cell index */
    size_type idx() const {
      return (i_-1)*(grid_->m_-2)+j_-1;
    }
    // Returns value of cell at relative index
    Proxy value(int di, int dj) {
      // Compute indices
      di = di+i_;
      dj = dj+j_;
      // Bound indices
      if (di<0) di=0;
      else if (di>int(grid_->n_-1)) di=int(grid_->n_-1);
      if (dj<0) dj=0;
      else if (dj>int(grid_->m_-1)) dj=int(grid_->m_-1);
      // Access value
      return Proxy(di,dj,grid_,ind_);
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
      return grid_->cell(di,dj,ind_);
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
    scalar_type dx() {
      return (this->x(2)[0]+this->x(3)[0]-this->x(1)[0]-this->x(4)[0])/2.0;
    }
    /* Returns average cell height */
    scalar_type dy() {
      return (this->x(1)[1]+this->x(2)[1]-this->x(3)[1]-this->x(4)[1])/2.0;
    }
   };

   // Return cell object by interior id number
   Cell cell(size_type idx, short id=-1) {
    return Cell(this,floor(idx/(m_-2))+1,idx%(m_-2)+1,id);
   }
   /** Return cell by id pair
    * @pre 0 <= i <= n_-1
    * @pre 0 <= j <= m_-1
    */
   Cell cell(size_type i, size_type j, short id=-1) {
    return Cell(this,i,j,id);
   }
  //
  // CELL ITERATOR
  //
  class CellIterator : private totally_ordered<CellIterator> {
  private:
    StructuredGrid* grid_;  // Grid pointer
    size_type idx_;         // Cell index
    short id_;              // RK stage index (-1 for true grid)
  public:
    // Value type from grid
    typedef Val                       Value;
    // Iterator traits (magic)
    typedef Cell                      value_type;
    typedef Cell*                     pointer;
    typedef Cell&                     reference;
    typedef std::ptrdiff_t            difference_type;
    typedef std::forward_iterator_tag iterator_category;
    // Public constructor
    CellIterator(const StructuredGrid* grid, 
                  size_type idx,
                  short id=-1) :
      grid_(const_cast<StructuredGrid*>(grid)), idx_(idx), id_(id) {}
    // Implement minimal methods
    Cell operator*() const {
      return grid_->cell(idx_,id_);
    }
    CellIterator& operator++() { ++idx_; return *this; }
    bool operator==(const CellIterator& iter) const {
      return idx_ == iter.idx_;
    }
  };
  // StructuredGrid methods to return a CellIterator
  CellIterator cell_begin(short id=-1) const {
    return CellIterator(this,0,id);
  }
  CellIterator cell_end(short id=-1) const {
    return CellIterator(this,v_.size(),id);
  }
  /** @class StructuredGrid::Access
   *  @brief Lightweight proxy object for accessing grid data
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
    // Types
    typedef value_type GridValue;
    typedef scalar_type GridScalar;
    // Return cell value
    Proxy operator()(size_type i, size_type j, short ind=-1) {
      return Proxy(i,j,grid_,ind);
    }
    // Return cell iterator
    StructuredGrid::CellIterator cell_begin(short id=-1) const {
      return grid_->cell_begin(id);
    }
    StructuredGrid::CellIterator cell_end(short id=-1) const {
      return grid_->cell_end(id);
    }
  };
  
  // Return access object
  Access access() {
    return Access(this);
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

    std::ofstream f_out(outputfile.c_str());
    f_out << (n_-1) << "," << (m_-1) << std::endl;
    // Write out elements
    for (auto it=x_.begin(); it!=x_.end(); ++it) {
      for (auto jt=it->begin(); (jt+1)!=it->end(); ++jt)
        f_out << (*jt) << ",";
      f_out << (*it).back() << std::endl;
    }
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
    std::cout << v[v.size()-1] << ")";
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
