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

#ifndef MABADILIKO_TABLE_H
#define MABADILIKO_TABLE_H


#include <map>
#include <iomanip>
#include <treerecs/tools/utils.h>
#include "AbstractGrid.h"

namespace treerecs {
template<typename T, typename rowKey, typename colKey>
class Table : public AbstractGrid<T, rowKey, colKey>{
  /*!
   * \class Table
   * \brief Table is a template container, as a 2 vector but elements can be accessed with any object.
   */

  std::unordered_map<rowKey, std::size_t> rowIndexes_;
  std::unordered_map<colKey, std::size_t> colIndexes_;

private:
/****************
 * Private getters
 */
  /// Returns the unsigned integer in association with the key.
  inline std::size_t _rowKeyToIndex(const rowKey& key) const {
    assert(_isRowKey(key)); return rowIndexes_.at(key);
  }

  /// Returns the unsigned integer in association with the key.
  inline std::size_t _colKeyToIndex(const colKey& key) const {
    assert(_isColKey(key)); return colIndexes_.at(key);
  }

  /// Check if the key, in row, exists.
  inline bool _isRowKey(const rowKey& key) const { return rowIndexes_.find(key) != rowIndexes_.end(); }

  /// Check if the key, in column, exists.
  inline bool _isColKey(const colKey& key) const { return colIndexes_.find(key) != colIndexes_.end(); }

  /// Returns a vector which contains keys of a map.
  template<typename O>
  std::vector<O> _getKeys(const std::unordered_map<O, std::size_t>& map) const{
    std::vector<O> ret;
    ret.reserve(map.size());
    for(auto pair = map.begin(); pair != map.end(); pair++)
      ret.push_back(pair->first);
    return ret;
  }

  /// Returns a vector which contains all row keys.
  std::vector<rowKey> _getRowKeys() const {
    return _getKeys(rowIndexes_);
  }

  /// Returns a vector which contains all col keys.
  std::vector<colKey> _getColKeys() const {
    return _getKeys(colIndexes_);
  }

/****************
 * Private setters
 */
  /// Add a key. Is the key exists, the previous index will be erased (and content lost).
  inline void _addRowKey(const rowKey& key, const std::size_t& index) {
    rowIndexes_[key] = index;
  }

  /// Add a key. Is the key exists, the previous index will be erased (and content lost).
  inline void _addColKey(const colKey& key, const std::size_t& index) {
    colIndexes_[key] = index;
  }

  /// Remove a key.
  inline void _removeRowKey(const rowKey& key) {
    if(not _isRowKey(key)) {
      std::cerr << "The row key is not correct." << std::endl;
      return;
    }
    auto pos = rowIndexes_.at(key);
    for(auto it = rowIndexes_.begin(); it != rowIndexes_.end(); it++){
      if(it->second > pos) it->second--;
    }
    rowIndexes_.erase(key);
  }

  /// Remove a key.
  inline void _removeColKey(const colKey& key) {
    if(not _isColKey(key)) {
      std::cerr << "The column key is not correct." << std::endl;
      return;
    }
    auto pos = colIndexes_.at(key);
    for(auto it = colIndexes_.begin(); it != colIndexes_.end(); it++){
      if(it->second > pos) it->second--;
    }
    colIndexes_.erase(key);
  }

  /// Remove all keys in row.
  virtual inline void _eraseAllRowKeys() { rowIndexes_.clear(); }

  /// Remove all keys in column.
  virtual inline void _eraseAllColKeys() { colIndexes_.clear(); }

  /// Fill a map of keys with new keys.
  template<typename O>
  void _setMap(std::unordered_map<O, std::size_t>& map, const std::vector<O>& keys){
    map.clear();
    for(std::size_t i = 0; i < keys.size(); ++i) map[keys.at(i)] = i;
  }

protected:

public:
/****************
 * Constructors
 */
  Table(): AbstractGrid<T, rowKey, colKey>() {};
  Table(const std::size_t nrow, const std::size_t ncol): AbstractGrid<T, rowKey, colKey>(nrow, ncol) {};
  Table(const Table<T, rowKey, colKey>& model) : AbstractGrid<T, rowKey, colKey>(model),
                                                 rowIndexes_(model.rowIndexes_),
                                                 colIndexes_(model.colIndexes_) {}
  Table(Table<T, rowKey, colKey>&& model) noexcept :
      AbstractGrid<T, rowKey, colKey>(model)
      , rowIndexes_(std::move(model.rowIndexes_))
      , colIndexes_(std::move(model.colIndexes_)){
  }

  Table(const std::size_t nrow, const std::size_t ncol, const T& model) :
                                                 AbstractGrid<T, rowKey, colKey>(nrow, ncol, model)
  {}

/****************
 * Destructor
 */
  ~Table() { rowIndexes_.clear(); colIndexes_.clear(); }


/****************
 * Getters
 */

  std::unordered_map<rowKey, std::size_t> getRowIndexesMap() const {
    /// Returns a map of objects with position in the grid (order), used to get elements in row.
    return rowIndexes_;
  }

  std::unordered_map<colKey, std::size_t> getColIndexesMap() const {
    /// Returns a map of objects with position in the grid (order), used to get elements in column.
    return colIndexes_;
  }

  Table<T, rowKey, colKey> extractRows(const std::vector<rowKey>& keys) const {
    /// Returns a Table which is an extract of the current one.
    Table<T, rowKey, colKey> res;
    for(auto key: keys){
      if(not this->_isRowKey(key))
        std::cerr << __FUNCTION__ << " error: " << key << " is not a correct key of this Table." << std::endl;
      res.addRow(this->getRow(key), key);
    }
    res.setColIndexes(this->getColIndexesMap());
    return res;
  }

  Table<T, rowKey, colKey> extractCols(const std::vector<colKey>& keys) const {
    /// Returns a Table which is an extract of the current one.
    Table<T, rowKey, colKey> res;
    for(auto key: keys){
      if(not this->_isColKey(key))
        std::cerr << __FUNCTION__ << " error: " << key << " is not a correct key of this Table." << std::endl;
      res.addCol(this->getCol(key), key);
    }
    res.setRowIndexes(this->getRowIndexesMap());
    return res;
  }

/****************
 * Setters
 */
 /// Set all row indexes with the key as the key of the row and the value, the row number.
  void setRowIndexes(const std::unordered_map<rowKey, std::size_t>& keys) {
    rowIndexes_ = keys;
  }

  /// Set all col indexes with the key as the key of the column and the value, the column number.
  void setColIndexes(const std::unordered_map<rowKey, std::size_t>& keys) {
    colIndexes_ = keys;
  }

  void setRowIndexes(const std::vector<rowKey>& keys) {
    if(keys.size() != this->nrow_) std::cerr << "Warning Table::setRowIndexes(...): incorrect number of keys." << std::endl;
    _setMap(rowIndexes_, keys);
  }

  void setColIndexes(const std::vector<colKey>& keys) {
    if(keys.size() != this->ncol_) std::cerr << "Warning Table::setColIndexes(...): incorrect number of keys." << std::endl;
    _setMap(colIndexes_, keys);
  }

  Table<T, rowKey, colKey>& operator=(const Table<T, rowKey, colKey>& model) {
    AbstractGrid<T, rowKey, colKey>::operator=(model);
    rowIndexes_ = model.rowIndexes_;
    colIndexes_ = model.colIndexes_;
    return *this;
  }

  Table<T, rowKey, colKey>& operator=(Table<T, rowKey, colKey>&& model) noexcept {
    AbstractGrid<T, rowKey, colKey>::operator=(model);
    rowIndexes_ = std::move(model.rowIndexes_);
    colIndexes_ = std::move(model.colIndexes_);
    return *this;
  }

  friend std::ostream& operator<<(std::ostream& os, const Table<T, rowKey, colKey>& tab) {
     /// Streams all elements of a table.
     tab.write(os, "\t", true, true, true);
     return os;
  }

};

} // namespace treerecs

#endif //MABADILIKO_TABLE_H
