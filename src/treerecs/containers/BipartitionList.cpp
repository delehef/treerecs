//
// File: BipartitionList.cpp
// Created by: Nicolas Galtier and Julien Dutheil
// Created on: Tue Apr 13 15:09 2007
// Modified by Nicolas Comte

/*
   Copyright or © or Copr. Bio++ Development Team, (November 16, 2004)

   This software is a computer program whose purpose is to provide classes
   for phylogenetic data analysis.

   This software is governed by the CeCILL  license under French law and
   abiding by the rules of distribution of free software.  You can  use,
   modify and/ or redistribute the software under the terms of the CeCILL
   license as circulated by CEA, CNRS and INRIA at the following URL
   "http://www.cecill.info".

   As a counterpart to the access to the source code and  rights to copy,
   modify and redistribute granted by the license, users are provided only
   with a limited warranty  and the software's author,  the holder of the
   economic rights,  and the successive licensors  have only  limited
   liability.

   In this respect, the user's attention is drawn to the risks associated
   with loading,  using,  modifying and/or developing or reproducing the
   software by the user in light of its specific status of free software,
   that may mean  that it is complicated to manipulate,  and  that  also
   therefore means  that it is reserved for developers  and  experienced
   professionals having in-depth computer knowledge. Users are therefore
   encouraged to load and test the software's suitability as regards their
   requirements in conditions enabling the security of their systems and/or
   data to be ensured and,  more generally, to use and operate it in the
   same conditions as regards security.

   The fact that you are presently reading this means that you have had
   knowledge of the CeCILL license and that you accept its terms.
 */

#include "BipartitionList.h"
#include <treerecs/tools/BipartitionTools.h>
#include "stddef.h"

#include <Bpp/Exceptions.h>
#include <Bpp/Text/TextTools.h>
#include <Bpp/Io/FileTools.h>

#include <iostream>
#include <climits> // defines CHAR_BIT
#include <treerecs/tools/PhyloTreeToolBox.h>

using namespace bpp;

namespace treerecs {

/****************************************************************/
/* utilitary classes required for sorting elements/bipartitions */
/****************************************************************/

class StringAndInt
{
public:
  int ind;
  std::string str;

public:
  StringAndInt() : ind(0),
                   str() {}
};

bool operator<(StringAndInt sai1, StringAndInt sai2)
{
  if (sai1.str < sai2.str)
    return true;
  return false;
}

/******************************************************************************/

class IntAndInt
{
public:
  size_t ind;
  int val;
};

bool operator<(IntAndInt iai1, IntAndInt iai2)
{
  if (iai1.val < iai2.val)
    return true;
  return false;
}


/******************************************************************************/

BipartitionList::BipartitionList(const bpp::PhyloTree& tr, bool sorted, std::vector<int>* index) :
    bitBipartitionList_(),
    elements_(),
    sorted_(sorted)
{
  std::size_t nbbip;

  elements_ = tr.getAllLeavesNames();

  if (tr.isRooted() and tr.getSons(tr.getRoot()).size() <= 2) {
    nbbip = tr.getNumberOfNodes() - 2;
  }
  else {
    nbbip = tr.getNumberOfNodes() - 1;
  }

  if (sorted)
    std::sort(elements_.begin(), elements_.end());

  std::size_t lword  = static_cast<std::size_t>(BipartitionTools::LWORD);
  std::size_t nbword = (elements_.size() + lword - 1) / lword;
  std::size_t nbint  = nbword * lword / (CHAR_BIT * sizeof(int));

  // std::cout << "nbbip = " << nbbip << std::endl;
  // std::cout << "lword = " << lword << std::endl;
  // std::cout << "nbword = " << nbword << std::endl;
  // std::cout << "nbint = " << nbword << std::endl;

  for (std::size_t i = 0; i < nbbip; i++)
  {
    bitBipartitionList_.push_back(new int[nbint]);
    for (std::size_t j = 0; j < nbint; j++)
    {
      bitBipartitionList_[i][j] = 0;
    }
  }

  std::size_t cpt = 0;
  std::vector<std::string> underlyingNames;
  buildBitBipartitions(tr, tr.getRoot(), bitBipartitionList_, elements_, &cpt, index);
}

/******************************************************************************/

BipartitionList::BipartitionList(
    const std::vector<std::string>& elements,
    const std::vector<int*>& bitBipL) :
    bitBipartitionList_(),
    elements_(elements),
    sorted_()
{
  size_t lword  = static_cast<size_t>(BipartitionTools::LWORD);
  size_t nbword = (elements.size() + lword - 1) / lword;
  size_t nbint  = nbword * lword / (CHAR_BIT * sizeof(int));

  for (size_t i = 0; i < bitBipL.size(); i++)
  {
    bitBipartitionList_.push_back(new int[nbint]);
    for (size_t j = 0; j < nbint; j++)
    {
      bitBipartitionList_[i][j] = bitBipL[i][j];
    }
  }

  std::vector<std::string> cpelements_ = elements;
  std::sort(cpelements_.begin(), cpelements_.end());
  if (cpelements_ == elements)
    sorted_ = true;
  else
    sorted_ = false;
}

/******************************************************************************/

BipartitionList::BipartitionList(const BipartitionList& bipL) :
    bitBipartitionList_(),
    elements_(bipL.elements_),
    sorted_(bipL.sorted_)
{
  size_t lword  = static_cast<size_t>(BipartitionTools::LWORD);
  size_t nbword = (bipL.getNumberOfElements() + lword - 1) / lword;
  size_t nbint  = nbword * lword / (CHAR_BIT * sizeof(int));

  bitBipartitionList_.resize(bipL.getNumberOfBipartitions());
  std::vector<int*> bitBipL = bipL.getBitBipartitionList();
  for (size_t i = 0; i < bipL.getNumberOfBipartitions(); i++)
  {
    bitBipartitionList_[i] = new int[nbint];
    for (size_t j = 0; j < nbint; j++)
    {
      bitBipartitionList_[i][j] = bitBipL[i][j];
    }
  }
}

/******************************************************************************/

BipartitionList& BipartitionList::operator=(const BipartitionList& bipL)
{
  size_t lword  = static_cast<size_t>(BipartitionTools::LWORD);
  size_t nbword = (bipL.getNumberOfElements() + lword - 1) / lword;
  size_t nbint  = nbword * lword / (CHAR_BIT * sizeof(int));

  for (size_t i = 0; i < bitBipartitionList_.size(); i++)
  {
    delete[] bitBipartitionList_[i];
  }
  bitBipartitionList_.resize(bipL.getNumberOfBipartitions());
  std::vector<int*> bitBipL = bipL.getBitBipartitionList();
  for (size_t i = 0; i < bipL.getNumberOfBipartitions(); i++)
  {
    bitBipartitionList_[i] = new int[nbint];
    for (size_t j = 0; j < nbint; j++)
    {
      bitBipartitionList_[i][j] = bitBipL[i][j];
    }
  }

  elements_ = bipL.elements_;
  sorted_   = bipL.sorted_;
  return *this;
}

/******************************************************************************/

BipartitionList::~BipartitionList()
{
  for (size_t i = 0; i < bitBipartitionList_.size(); i++)
  {
    delete[] bitBipartitionList_[i];
  }
}

/******************************************************************************/

std::map<std::string, bool> BipartitionList::getBipartition(size_t i) const
noexcept(false)
{
  std::map<std::string, bool> bip;

  if (i >= bitBipartitionList_.size())
    throw Exception("Bipartition index exceeds BipartitionList size");

  for (size_t j = 0; j < elements_.size(); j++)
  {
    if (BipartitionTools::testBit(bitBipartitionList_[i], static_cast<int>(j)))
      bip[elements_[j]] = true;
    else
      bip[elements_[j]] = false;
  }
  return bip;
}

/******************************************************************************/

int* BipartitionList::getBitBipartition(size_t i) noexcept(false)
{
  if (i >= bitBipartitionList_.size())
    throw Exception("Bipartition index exceeds BipartitionList size");

  return bitBipartitionList_[i];
}

/******************************************************************************/

bool BipartitionList::haveSameElementsThan(std::map<std::string, bool>& bipart) const
{
  std::vector<std::string> elements = elements_;
  std::vector<std::string> keys;

  std::map<std::string, bool>::iterator it;

  for (it = bipart.begin(); it != bipart.end(); it++)
  {
    keys.push_back(it->first);
  }

  std::sort(elements.begin(), elements.end());
  std::sort(keys.begin(), keys.end());

  if (elements == keys)
    return true;
  return false;
}

/******************************************************************************/

void BipartitionList::addBipartition(std::map<std::string, bool>& bipart, bool checkElements)
noexcept(false)
{
  if (checkElements && !BipartitionList::haveSameElementsThan(bipart))
    throw Exception("Distinct bipartition element sets");

  size_t lword  = static_cast<size_t>(BipartitionTools::LWORD);
  size_t nbword = (elements_.size() + lword - 1) / lword;
  size_t nbint  = nbword * lword / (CHAR_BIT * sizeof(int));
  bitBipartitionList_.push_back(new int[nbint]);
  size_t ind    = bitBipartitionList_.size() - 1;
  for (size_t j = 0; j < nbint; j++)
  {
    bitBipartitionList_[ind][j] = 0;
  }

  for (size_t i = 0; i < elements_.size(); i++)
  {
    if (bipart[elements_[i]] == true)
      BipartitionTools::bit1(bitBipartitionList_[ind], static_cast<int>(i));
    else
      BipartitionTools::bit0(bitBipartitionList_[ind], static_cast<int>(i));
  }
}

/******************************************************************************/

void BipartitionList::deleteBipartition(size_t i) noexcept(false)
{
  if (i >= bitBipartitionList_.size())
    throw Exception("Bipartition index exceeds BipartitionList size");

  delete[] bitBipartitionList_[i];
  bitBipartitionList_.erase(bitBipartitionList_.begin() + static_cast<ptrdiff_t>(i));
}

/******************************************************************************/

bool BipartitionList::containsBipartition(std::map<std::string, bool>& bipart, bool checkElements) const
noexcept(false)
{
  size_t i, j;
  bool dac, padac;

  if (checkElements && !BipartitionList::haveSameElementsThan(bipart))
    throw Exception("Distinct bipartition element sets");

  for (i = 0; i < bitBipartitionList_.size(); i++)
  {
    dac = padac = false;
    for (j = 0; j < elements_.size(); j++)
    {
      if (BipartitionTools::testBit(bitBipartitionList_[i], static_cast<int>(j)))
      {
        if (bipart[elements_[j]])
          dac = true;
        else
          padac = true;
      }
      else
      {
        if (bipart[elements_[j]])
          padac = true;
        else
          dac = true;
      }
      if (dac && padac)
        break;
    }
    if (j == elements_.size())
      return true;
  }
  return false;
}

/******************************************************************************/

bool BipartitionList::areIdentical(size_t k1, size_t k2) const noexcept(false)
{
  bool dac, padac;

  if (k1 >= bitBipartitionList_.size())
    throw Exception("Bipartition index exceeds BipartitionList size");
  if (k2 >= bitBipartitionList_.size())
    throw Exception("Bipartition index exceeds BipartitionList size");

  dac = padac = false;
  for (size_t j = 0; j < elements_.size(); j++)
  {
    if (BipartitionTools::testBit(bitBipartitionList_[k1], static_cast<int>(j)))
    {
      if (BipartitionTools::testBit(bitBipartitionList_[k2], static_cast<int>(j)))
        dac = true;
      else
        padac = true;
    }
    else
    {
      if (BipartitionTools::testBit(bitBipartitionList_[k2], static_cast<int>(j)))
        padac = true;
      else
        dac = true;
    }
    if (dac && padac)
      return false;
  }
  return true;
}

/******************************************************************************/

bool BipartitionList::areCompatible(size_t k1, size_t k2) const noexcept(false)
{
  bool uu, uz, zu, zz;

  if (k1 >= bitBipartitionList_.size())
    throw Exception("Bipartition index exceeds BipartitionList size");
  if (k2 >= bitBipartitionList_.size())
    throw Exception("Bipartition index exceeds BipartitionList size");

  uu = uz = zu = zz = false;

  for (size_t j = 0; j < elements_.size(); j++)
  {
    if (BipartitionTools::testBit(bitBipartitionList_[k1], static_cast<int>(j)))
    {
      if (BipartitionTools::testBit(bitBipartitionList_[k2], static_cast<int>(j)))
        uu = true;
      else
        uz = true;
    }
    else
    {
      if (BipartitionTools::testBit(bitBipartitionList_[k2], static_cast<int>(j)))
        zu = true;
      else
        zz = true;
    }
    if (uu && uz && zu && zz)
      return false;
  }

  return true;
}

/******************************************************************************/

bool BipartitionList::areAllCompatible() const
{
  for (size_t i = 0; i < bitBipartitionList_.size(); i++)
  {
    for (size_t j = i + 1; j < bitBipartitionList_.size(); j++)
    {
      if (!BipartitionList::areCompatible(i, j))
        return false;
    }
  }
  return true;
}

/******************************************************************************/

bool BipartitionList::areAllCompatibleWith(std::map<std::string, bool>& bipart, bool checkElements) const
noexcept(false)
{
  if (checkElements && !haveSameElementsThan(bipart))
    throw Exception("Distinct bipartition element sets");
  size_t nbBip = bitBipartitionList_.size();
  const_cast<BipartitionList*>(this)->addBipartition(bipart, false);

  for (size_t i = 0; i < nbBip; i++)
  {
    if (!areCompatible(i, nbBip))
    {
      const_cast<BipartitionList*>(this)->deleteBipartition(nbBip);
      return false;
    }
  }
  const_cast<BipartitionList*>(this)->deleteBipartition(nbBip);
  return true;
}

/******************************************************************************/

void BipartitionList::sortElements()
{
  std::vector<StringAndInt> relements_;
  StringAndInt sai;
  size_t nbbip;

  for (size_t i = 0; i < elements_.size(); i++)
  {
    sai.str = elements_[i];
    sai.ind = static_cast<int>(i);
    relements_.push_back(sai);
  }

  std::sort(relements_.begin(), relements_.end());

  for (size_t i = 0; i < elements_.size(); i++)
  {
    elements_[i] = relements_[i].str;
  }

  nbbip = bitBipartitionList_.size();
  bitBipartitionList_.resize(2 * nbbip);
  size_t lword  = static_cast<size_t>(BipartitionTools::LWORD);
  size_t nbword = (elements_.size() + lword - 1) / lword;
  size_t nbint  = nbword * lword / (CHAR_BIT * sizeof(int));

  for (size_t j = nbbip; j < 2 * nbbip; j++)
  {
    bitBipartitionList_[j] = new int[nbint];
    for (size_t k = 0; k < nbint; k++)
    {
      bitBipartitionList_[j][k] = 0;
    }
    for (size_t i = 0; i < elements_.size(); i++)
    {
      if (BipartitionTools::testBit(bitBipartitionList_[j - nbbip], relements_[i].ind))
        BipartitionTools::bit1(bitBipartitionList_[j], static_cast<int>(i));
      else
        BipartitionTools::bit0(bitBipartitionList_[j], static_cast<int>(i));
    }
  }

  for (size_t j = 0; j < nbbip; j++)
  {
    delete[] bitBipartitionList_[j];
  }

  bitBipartitionList_.erase(bitBipartitionList_.begin(), bitBipartitionList_.begin() + static_cast<ptrdiff_t>(nbbip));
  sorted_ = true;
}

/******************************************************************************/

size_t BipartitionList::getPartitionSize(std::size_t k) const noexcept(false)
{
  std::size_t size = 0;
  if (k >= bitBipartitionList_.size())
    throw Exception("Bipartition index exceeds BipartitionList size");

  for (std::size_t i = 0; i < elements_.size(); i++)
  {
    if (BipartitionTools::testBit(bitBipartitionList_[k], static_cast<int>(i)))
      size++;
  }

  if (size <= elements_.size() / 2)
    return size;
  else
    return elements_.size() - size;
}

/******************************************************************************/

void BipartitionList::removeTrivialBipartitions()
{
  size_t size = bitBipartitionList_.size();
  for (size_t i = size; i > 0; i--)
  {
    if (BipartitionList::getPartitionSize(i - 1) < 2)
      BipartitionList::deleteBipartition(i - 1);
  }
}

/******************************************************************************/

void BipartitionList::addTrivialBipartitions(bool checkExisting)
{
  std::map<std::string, bool> bip;

  for (size_t i = 0; i < elements_.size(); i++)
  {
    bip[elements_[i]] = false;
  }
  for (size_t i = 0; i < elements_.size(); i++)
  {
    bip[elements_[i]] = true;
    if (checkExisting && BipartitionList::containsBipartition(bip, false))
      continue;
    BipartitionList::addBipartition(bip, false);
    bip[elements_[i]] = false;
  }
}

/******************************************************************************/

void BipartitionList::sortByPartitionSize()
{
  std::vector<int*> sortedBitBipL;
  std::vector<IntAndInt> iaiVec;
  IntAndInt iai;

  for (size_t i = 0; i < bitBipartitionList_.size(); i++)
  {
    iai.ind = i;
    iai.val = static_cast<int>(BipartitionList::getPartitionSize(i));
    iaiVec.push_back(iai);
  }

  std::sort(iaiVec.begin(), iaiVec.end());

  for (size_t i = 0; i < bitBipartitionList_.size(); i++)
  {
    sortedBitBipL.push_back(bitBipartitionList_[iaiVec[i].ind]);
  }

  bitBipartitionList_ = sortedBitBipL;
}

/******************************************************************************/

void BipartitionList::flip(std::size_t k) noexcept(false)
{
  if (k >= bitBipartitionList_.size())
    throw Exception("Bipartition index exceeds BipartitionList size");
  size_t lword = static_cast<size_t>(BipartitionTools::LWORD);
  size_t nbword = (elements_.size() + lword - 1) / lword;
  size_t nbint = nbword * lword / (CHAR_BIT * sizeof(int));
  int* flipbip = new int[nbint];
  for (size_t i = 0; i < nbint; i++)
  {
    flipbip[i] = 0;
  }
  BipartitionTools::bitNot(flipbip, bitBipartitionList_[k], nbint);
  delete[] bitBipartitionList_[k];
  bitBipartitionList_[k] = flipbip;
}

/******************************************************************************/

void BipartitionList::removeRedundantBipartitions()
{
  bool deletion = true;

  while (deletion)
  {
    deletion = false;
    for (size_t i = 0; i < bitBipartitionList_.size(); i++)
    {
      for (size_t j = i + 1; j < bitBipartitionList_.size(); j++)
      {
        if (BipartitionList::areIdentical(i, j))
        {
          BipartitionList::deleteBipartition(j);
          deletion = true;
          break;
        }
      }
      if (deletion)
        break;
    }
  }
}

/******************************************************************************/

std::shared_ptr<bpp::PhyloTree> BipartitionList::toTree(const bool rooted) const
noexcept(false)
{

  BipartitionList* sortedBipL;
  std::vector<int*> sortedBitBipL;
  int* bip;
  std::vector<std::shared_ptr<bpp::PhyloNode>> vecNd, sonNd;
  std::vector<bool> alive;
  std::size_t lword, nbword, nbint, ii;
  bool verbose = false;
  std::size_t ancestor_index = 0;

  std::shared_ptr<bpp::PhyloTree> resulting_tree(new bpp::PhyloTree(false));

  // check, copy and prepare bipartition list


  if (!BipartitionList::areAllCompatible())
    throw bpp::Exception("Trying to build a tree from incompatible bipartitions");

  sortedBipL = dynamic_cast<BipartitionList*>(clone());
  for (std::size_t i = 0; i < sortedBipL->getNumberOfBipartitions(); i++)
  {
    auto bipartitionSize_i = sortedBipL->getPartitionSize(i);
    auto numberOfElements = sortedBipL->getNumberOfElements();
    if (bipartitionSize_i > numberOfElements / 2)
      sortedBipL->flip(i);
  }
  sortedBipL->sortByPartitionSize();
  sortedBipL->removeRedundantBipartitions();

  if(rooted) {
    auto last_bipartition = sortedBipL->getBipartition(sortedBipL->getNumberOfBipartitions() - 1);
    sortedBipL->addBipartition(last_bipartition);
    sortedBipL->flip(sortedBipL->getNumberOfBipartitions() - 1);
  }

  sortedBitBipL = sortedBipL->getBitBipartitionList();

  for (std::size_t i = 0; i < sortedBipL->getNumberOfBipartitions(); i++)
  {
    alive.push_back(true);
  }
  vecNd.resize(sortedBipL->getNumberOfBipartitions() + 1);
  lword  = static_cast<std::size_t>(BipartitionTools::LWORD);
  nbword = (elements_.size() + lword - 1) / lword;
  nbint  = nbword * lword / (CHAR_BIT * sizeof(int));
  bip    = new int[1]; bip[0] = 0;

  if(verbose) std::cout << "SortedBipL: " << std::endl << sortedBipL->toMatrix() << std::endl;

  //main loop: create one node per bipartition
  for (std::size_t i = 0; i < sortedBipL->getNumberOfBipartitions(); i++)
  {
    if(verbose) std::cout << "> Bipartition: " << i << std::endl;
    if(verbose) std::cout << "  state of vdNd = " << vecNd << std::endl;
    if(verbose) std::cout << "  " << sortedBipL->getBipartition(i) << std::endl;
    if(verbose) std::cout << "  Bipartition size = " << sortedBipL->getPartitionSize(i) << std::endl;
    if (sortedBipL->getPartitionSize(i) == 1)
    { // terminal, leaf
      for (std::size_t j = 0; j < sortedBipL->getNumberOfElements(); j++)
      {
        if (BipartitionTools::testBit(sortedBitBipL[i], static_cast<int>(j)))
        {
          vecNd[i] = std::shared_ptr<bpp::PhyloNode>(new bpp::PhyloNode(elements_[j]));
          resulting_tree->createNode(vecNd[i]);
          if(verbose) std::cout << "  > Add leaf: " << vecNd[i] << std::endl;
          break;
        }
      }
    }
    else
    { // internal node
      sonNd.clear();
      for (std::size_t j = 0; j < i; j++) // We are looking for bifurcations which are not connected to the final graph
      {
        if (alive[j])
        {
          for (ii = 0; ii < nbint; ii++)
          {
            BipartitionTools::bitOr(bip, sortedBitBipL[j] + ii, sortedBitBipL[i] + ii, 1);
            if (bip[0] != sortedBitBipL[i][ii]) {
              break;
            }
          }
          if (ii == nbint)
          {
            sonNd.push_back(vecNd[j]);
            alive[j] = false;
          }
        }
      }
      if(verbose) std::cout << "  > Add ancestral node which connects" << sonNd << std::endl;
      vecNd[i] = std::shared_ptr<bpp::PhyloNode>(new bpp::PhyloNode("ancestor " + std::to_string(ancestor_index++)));
      resulting_tree->createNode(vecNd[i]);

      for(auto& son: sonNd) {
        resulting_tree->link(vecNd[i], son);//, std::shared_ptr<bpp::PhyloBranch>(new bpp::PhyloBranch()));
      }
    }
    if(verbose) std::cout << "  state of vedNd = " << vecNd << std::endl;
  }

  sonNd.clear();

  if(verbose) std::cout << "Rooting..." << std::endl;
  auto rootNd = std::shared_ptr<bpp::PhyloNode>(new bpp::PhyloNode("root"));
  resulting_tree->createNode(rootNd);
  for (std::size_t i = 0; i < sortedBipL->getNumberOfBipartitions(); i++)
  {
    if (alive[i] and vecNd[i] != nullptr) {
      resulting_tree->link(rootNd, vecNd[i]);
    }
  }

  // Because the tree is undirected, we need to create one which is directed.
  auto result = PhyloTreeToolBox::cloneTree(*resulting_tree, rootNd, false, false);
  if(verbose) std::cout << "done" << std::endl;

  //construct tree and return
  delete sortedBipL;
  return result;
}

/******************************************************************************/
std::vector<std::string> BipartitionList::buildBitBipartitions(const bpp::PhyloTree& tree, const std::shared_ptr<bpp::PhyloNode>& nd, std::vector<int*>& bitbip, const std::vector<std::string>& elements, std::size_t* cpt, std::vector<int>* index) const
{
  std::vector<std::string> underelements_, retelements_;


  // Fill underelements with nodes under the current.
  if (tree.getNumberOfSons(nd) == 0)
    underelements_.push_back(nd->getName());

  for (std::size_t i = 0; i < tree.getNumberOfSons(nd); i++)
  {
    retelements_ = BipartitionList::buildBitBipartitions(tree, tree.getSons(nd)[i], bitbip, elements, cpt, index);
    for (std::size_t j = 0; j < retelements_.size(); j++)
    {
      underelements_.push_back(retelements_[j]);
    }
  }

  if (!tree.hasFather(nd))
    return underelements_;  // root node

  auto nd_father = tree.getFather(nd);

  if (!tree.hasFather(nd_father)) // If the node is the direct right son of the root, stop
  {
    std::size_t nbrootson = tree.getNumberOfSons(nd_father);
    if (nbrootson == 2 && nd == tree.getSons(nd_father)[1])
      return underelements_;  // son 2 of root node when root node has 2 sons
  }

  bool ones;
  if (underelements_.size() <= elements.size() / 2)
    ones = true;
  else
    ones = false;

  for (std::size_t i = 0; i < elements.size(); i++)
  {
    if (ones)
      BipartitionTools::bit0(bitbip[*cpt], static_cast<int>(i));
    else
      BipartitionTools::bit1(bitbip[*cpt], static_cast<int>(i));
  }

  for (std::size_t i = 0; i < underelements_.size(); i++)
  {
    std::size_t taxa_ind = 0;
    while (underelements_[i] != elements[taxa_ind])
      taxa_ind++;
    if (ones)
      BipartitionTools::bit1(bitbip[*cpt], static_cast<int>(taxa_ind));
    else
      BipartitionTools::bit0(bitbip[*cpt], static_cast<int>(taxa_ind));
  }

  (*cpt)++;

  if (index)
    index->push_back(tree.getNodeIndex(nd));

  return underelements_;
}

/******************************************************************************/
Table<int, std::string, std::size_t> BipartitionList::toMatrix() const
{
  std::vector< std::map<std::string, bool> > bipl;
  for (std::size_t i = 0; i < getNumberOfBipartitions(); i++)
  {
    bipl.push_back(getBipartition(i));
  }

  std::vector<std::string> el = getElementNames();

  //RowMatrix<int> mat(el.size(), getNumberOfBipartitions());

  Table<int, std::string, std::size_t> mat(el.size(), getNumberOfBipartitions());
  mat.setRowIndexes(el);
  std::vector<std::size_t> col_keys(getNumberOfBipartitions());
  std::iota(col_keys.begin(), col_keys.end(), 1);
  mat.setColIndexes(col_keys);

  for (std::size_t j = 0; j < el.size(); j++)
  {
    auto element = el[j];
    for (std::size_t i = 0; i < getNumberOfBipartitions(); i++)
    {
      mat(element, i + 1) = bipl[i][el[j]];
    }
  }
  return mat;
}
/******************************************************************************/

} // namespace treerecs
