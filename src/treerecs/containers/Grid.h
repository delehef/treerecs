// Copyright (C) 2018  INRIA
//
// Treerecs is free software: you can redistribute it and/or modify
// it under the terms of the GNU Affero General Public License as
// published by the Free Software Foundation, either version 3 of the
// License, or (at your option) any later version.
//
// Treerecs is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Affero General Public License for more details.
//
// You should have received a copy of the GNU Affero General Public License
// along with this program.  If not, see <https://www.gnu.org/licenses/>.

#ifndef GRID_GRID_H
#define GRID_GRID_H

#include "AbstractGrid.h"

namespace treerecs {

template<typename T>
class Grid: public AbstractGrid<T, std::size_t, std::size_t> {
  /*!
   * \class Grid
   * \brief 2D vector, as a matrix.
   * \details Grid is a template matrix which can store a template type.
   */
private:

protected:

/****************
 * Protected getters
 */
  inline std::size_t _rowKeyToIndex(const std::size_t& key) const { return key; }

  inline std::size_t _colKeyToIndex(const std::size_t& key) const { return key; }

  inline bool _isRowKey(const std::size_t& key) const { return key < this->nrow_; }

  inline bool _isColKey(const std::size_t& key) const { return key < this->ncol_; }

  std::vector<std::size_t> _getRowKeys() const {
    std::vector<std::size_t> ret;
    ret.reserve(this->nrow_);
    std::size_t i = 0;
    while(i < this->nrow_) ret.emplace_back(i++);
    return ret;
  }

  std::vector<std::size_t> _getColKeys() const {
    std::vector<std::size_t> ret;
    ret.reserve(this->ncol_);
    std::size_t i = 0;
    while(i < this->ncol_) ret.emplace_back(i++);
    return ret;
  }

/****************
 * Protected setters
 */
  // Empty definitions because there is no structure to manage indexes.
  inline void _addRowKey(const std::size_t&, const std::size_t&) {}
  inline void _addColKey(const std::size_t&, const std::size_t&) {}
  inline void _removeRowKey(const std::size_t&) {}
  inline void _removeColKey(const std::size_t&) {}
  inline virtual void _eraseAllRowKeys() {}
  inline virtual void _eraseAllColKeys() {}

public:
/****************
 * Constructors
 */
  Grid(): AbstractGrid<T, std::size_t, std::size_t>() {}
  Grid(const std::size_t nrow, const std::size_t ncol) : AbstractGrid<T, std::size_t, std::size_t>(nrow, ncol) {}
  Grid(const std::size_t nrow, const std::size_t ncol, const T& model) : AbstractGrid<T, std::size_t, std::size_t>(nrow, ncol, model) {}
  Grid(const std::vector<std::vector<T>>& model) : AbstractGrid<T, std::size_t, std::size_t>(model) {}
  Grid(const Grid<T>& model) : AbstractGrid<T, std::size_t, std::size_t>(model) {}
  Grid(Grid<T>&& model) : AbstractGrid<T, std::size_t, std::size_t>(model) {}

/****************
 * Destructor
 */
  ~Grid() {};

/****************
 * Setters
 */
  void addRow(const std::vector<T> &row){ AbstractGrid<T, std::size_t, std::size_t>::addRow(row, this->nrow_); }

  void addCol(const std::vector<T> &col){ AbstractGrid<T, std::size_t, std::size_t>::addCol(col, this->ncol_); }

/****************
 * Getters
 */
  void write(std::ostream& os) const {
    /// Write the matrix in an output stream.
    for(std::size_t i = 0; i < this->nrow_ ; i++){
      for(std::size_t j = 0; j < this->ncol_ ; j++) {
        os << this->data_[i][j];
      }
      os << std::endl;
    }
    os << std::endl;
  }

  friend std::ostream& operator<<(std::ostream& os, const Grid<T> grid){
    grid.write(os);
    return os;
  }
};

} // namespace treerecs

#endif //GRID_GRID_H
