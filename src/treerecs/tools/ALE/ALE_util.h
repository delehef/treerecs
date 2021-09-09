//all code by Szollosi GJ et al.; ssolo@elte.hu; CC BY-SA 3.0;
//modified by N.Comte
#pragma once
#include "ALE.h"

namespace treerecs {

inline bool fexists(const std::string& filename)
{
  std::ifstream ifile(filename);
  return ifile.good();
}

approx_posterior * observe_ALE_from_nexus(std::string fname,int burnin=100, int every=1,int until=-1);
approx_posterior * observe_ALE_from_file(std::string fname,int burnin=100, int every=1,int until=-1); // NO del-loc
approx_posterior * observe_ALE_from_file(std::vector<std::string> fnames,int burnin=100, int every=1,int until=-1); // NO del-loc
approx_posterior * observe_ALE_from_strings(std::vector<std::string> trees); // NO del-loc
approx_posterior * observe_ALE_from_string(std::string tree); // NO del-loc
approx_posterior * load_ALE_from_file(std::string fname); // NO del-loc
std::string save_ALE_to_file(std::string fname);

std::string canonical_branch_lengths(std::string Sstring);
void canonical_branch_lengths( tree_type * S );

} // namespace treerecs
