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

#include "ReconcilableRootedTree.h"

#include <treerecs/tools/PhyloTreeToolBox.h>
#include <treerecs/tools/utils.h>

namespace treerecs {

void ReconcilableRootedTree::setSubTreeRootOccurrence_recursive_propagation(const Node &node, const std::size_t &n) {
  if(node != rooted_multifurcating_genetree_->getRoot()) {
    if(rooted_multifurcating_genetree_->getSons(node).size() <= 2)
      setSubTreeRootOccurrence_recursive_propagation(rooted_multifurcating_genetree_->getFather(node), n);
  }

  if(subtree_root_occurrences_.find(node) == subtree_root_occurrences_.end())
    subtree_root_occurrences_[node] = 1;

  subtree_root_occurrences_[node] *= n;
}

ReconcilableRootedTree::ReconcilableRootedTree(
    const std::shared_ptr<bpp::PhyloTree> &genetree, const std::shared_ptr<bpp::PhyloTree> &guidetree
    , const SpeciesGeneMap &speciesGeneMap, const std::size_t &bifurcation_duplication_number
    , const std::size_t &bifurcation_loss_number, const double &duplication_cost, const double &loss_cost
    , const std::unordered_map<Node, CostTable> &costMatrices, const bool rerooted, const double &execution_time
) :
    rooted_multifurcating_genetree_(genetree)
    , guidetree_(guidetree)
    , speciesGeneMap_(speciesGeneMap)
    , bifurcation_duplication_number_(bifurcation_duplication_number)
    , bifurcation_loss_number_(bifurcation_loss_number)
    , duplication_cost_(duplication_cost)
    , loss_cost_(loss_cost)
    , costMatrices_(costMatrices)
    , rerooted_(rerooted)
    , execution_time_(execution_time)
{
  if(costMatrices.size() > 0) {
    for (auto &costmatrix: costMatrices) {
      const std::shared_ptr<bpp::PhyloNode> &polytomy_root = costmatrix.first;
      std::size_t polytomy_root_occurrence = utils::sum(
          costmatrix.second(speciesGeneMap.getAssociatedSpecies(polytomy_root), 1).occurence);

      setSubTreeRootOccurence(polytomy_root, polytomy_root_occurrence);
    }
  } else {
    setSubTreeRootOccurence(genetree->getRoot(), 1);
  }
}

double ReconcilableRootedTree::cost() const {
  double res = 0.0;

  // Add bifurcation costs
  res += ((double)bifurcation_duplication_number_ * duplication_cost_ + (double)bifurcation_loss_number_ * loss_cost_);

  // Add multifurcation costs
  for(auto& mf_matrix: costMatrices_) {
    auto& polytomy_root = mf_matrix.first;
    auto& matrix = mf_matrix.second;
    res += matrix(speciesGeneMap_.getAssociatedSpecies(polytomy_root), 1).value;
  }

  return res;
}

bool ReconcilableRootedTree::has_costtable_at(const Node &node) const {
  return costMatrices_.find(node) != costMatrices_.end();
}

const CostTable &ReconcilableRootedTree::get_costtable_at(const Node &node) const {
  auto costMatrix_it = costMatrices_.find(node);
  if (costMatrix_it == costMatrices_.end())
    std::cerr << "Cost table does not exist for node " << node << "." << std::endl;

  return costMatrix_it->second;
}

double ReconcilableRootedTree::execution_time() const {
  return execution_time_;
}

std::shared_ptr<const bpp::PhyloTree> ReconcilableRootedTree::genetree() const {
  return this->rooted_multifurcating_genetree_;
}

std::shared_ptr<const bpp::PhyloTree> ReconcilableRootedTree::guidetree() const {
  return this->guidetree_;
}

const SpeciesGeneMap &ReconcilableRootedTree::map() const {
  return this->speciesGeneMap_;
}

std::size_t ReconcilableRootedTree::n_multifurcations() const { return costMatrices_.size(); }

void ReconcilableRootedTree::setSubTreeRootOccurence(const Node &node, const std::size_t &n) {

  if(subtree_root_occurrences_.find(node) == subtree_root_occurrences_.end())
    subtree_root_occurrences_[node] = 1;

  subtree_root_occurrences_[node] *= n;

  if(node != rooted_multifurcating_genetree_->getRoot()) {
    setSubTreeRootOccurrence_recursive_propagation(rooted_multifurcating_genetree_->getFather(node), subtree_root_occurrences_.at(node));
  }
}

std::size_t ReconcilableRootedTree::occurrence(void) const {
  if(not rooted_multifurcating_genetree_){
    std::cerr << "Error: there is no original tree in solution." << std::endl;
    print();
  }

  std::size_t res = subtree_root_occurrences_.at(rooted_multifurcating_genetree_->getRoot());

  return res;
}

std::size_t ReconcilableRootedTree::bifurcation_duplication_number() const { return bifurcation_duplication_number_; }

std::size_t ReconcilableRootedTree::bifurcation_loss_number() const { return bifurcation_loss_number_; }

double ReconcilableRootedTree::duplication_cost() const { return duplication_cost_; }

double ReconcilableRootedTree::loss_cost() const { return loss_cost_; }

bool ReconcilableRootedTree::rerooted() const { return rerooted_; }

void ReconcilableRootedTree::print() const {
  std::cout << "Solution: " << std::endl;
  std::cout << "> total_cost = " << cost() << std::endl;
  std::cout << "> multifurcations = " << costMatrices_.size() << std::endl;
  std::cout << "> bifurcation duplications = " << bifurcation_duplication_number_ << std::endl;
  std::cout << "> bifurcation losses = " << bifurcation_loss_number_ << std::endl;
  //std::cout << "> Matrices:" << std::endl;
  //for(auto matrix: costMatrices_){
  //  std::cout << "At polytomy: " << matrix.first << ":" << matrix.second.get(speciesGeneMap_.getAssociatedSpecies(matrix.first), 1) << std::endl;
  //}
}

} // namespace treerecs
