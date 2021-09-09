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
// Ale includes

#include <memory>
#include <vector>
#include <string>

// Treerecs containers includes
#include <treerecs/containers/BipartitionList.h>

// Treerecs tools includes
#include <treerecs/tools/ALE/ALE.h>
#include <treerecs/tools/ALE/exODT.h>
#include <treerecs/tools/utils.h>
#include <treerecs/tools/IO/RefreshablePrinter.h>
#include <treerecs/tools/Timer.h>
#include <treerecs/tools/IO/IO.h>
#include <treerecs/tools/PhyloTreeToolBox.h>
#include <treerecs/tools/BipartitionTools.h>
#include <treerecs/tools/SpeciesGeneMapper.h>

using namespace treerecs;

/// @brief Test if a file exists according to a given filename.
/// @param filename filename.
inline bool fexists(const std::string& filename) {
  std::ifstream ifile(filename);
  return ifile.good();
}

/// @brief Extract tree in parenthesis format from a file.
/// @param fname filename.
std::string readTreeFromFile(const std::string& fname) {
  if (!IO::exists(fname)) {
    std::cerr << "Error, file "<< fname
              << " does not seem accessible." << std::endl;
    exit(1);
  }
  std::ifstream file_stream(fname.c_str());
  std::string tree_str = "";
  if (file_stream.is_open()) {  //  ########## read trees ############
    while (not file_stream.eof()) {  // read all file to check a tree.
      std::string line;
      std::getline(file_stream, line);
      if (line.find("(") != line.npos )
      {
        tree_str = line;
        break;
      }
    }
  }
  return tree_str;
}

int main(int argc, char ** argv) {
  std::cout << "ALEestimate using ALE v" << ALE_VERSION << std::endl;

  if (argc < 2) {
    std::cout << "usage:\n ./Ale species_tree_file gene_tree_file "
        "[separators=gene_name_separator [position=(after, before)]] "
        "[smap=leaves_map_file] O_R=OriginationAtRoot "
        "delta=DuplicationRate tau=TransferRate lambda=LossRate "
        "beta=weight_of_sequence_evidence outputFiles=n" << std::endl;
    return 1;
  }

  // Species tree business
  std::string species_tree_file = argv[1];
  utils::trim_str(species_tree_file);
  // Reading the species tree from within the file
  std::string species_tree_str = readTreeFromFile(species_tree_file);
  std::cout << "\n\tRead species tree from: " << argv[1] << std::endl;

  // Gene tree business
  std::string gene_tree_file = argv[2];
  utils::trim_str(gene_tree_file);

  std::string head = gene_tree_file;
  std::string ale_name = head + ".ale";

  // Reading the gene tree from within the file
  std::string gene_tree_str = readTreeFromFile(gene_tree_file);

  // Create an approx posterior instance using gene_tree.
  std::shared_ptr<approx_posterior> ale(new approx_posterior(gene_tree_str));

  std::vector<std::string> gene_tree_strs;
  gene_tree_strs.reserve(1);
  // Silly: we need to produce a vector with a single element...
  gene_tree_strs.push_back(gene_tree_str);
  ale->observation(gene_tree_strs);
  std::cout << "\n\tObserved " << ale->observations
            << " gene tree(s) from: " <<  argv[2] << std::endl;

  // We initialise a coarse grained reconciliation model for calculating the sum
  exODT_model* model = new exODT_model();

  // Getting the other options
  scalar_type samples = 100;
  scalar_type O_R = 1, beta = 1;
  scalar_type delta = 0.01, tau = 0.01, lambda = 0.1;
  std::string fractionMissingFile = "";
  bool outputyn = false;  // output files option.

  MappingMethod mapping_method = MappingMethod::TREES;
  std::string smap_file_name = "";

  for (int i = 3; i < argc; i++) {
    std::string next_field=argv[i];
    std::vector <std::string> tokens;
    // boost::split(tokens,next_field,boost::is_any_of("="),
    // boost::token_compress_on);
    tokens = utils::splitString(next_field, "=");
    if (tokens[0] == "sample") {
      samples = atoi(tokens[1].c_str());
    } else if ( tokens[0] == "separators" ) {
      model->set_model_parameter("gene_name_separators", tokens[1]);
      mapping_method = MappingMethod::GENE_NAMES;
    } else if ( tokens[0] == "position" ) {
      if ( not (utils::string_comp(tokens[1], "after", false)
              or utils::string_comp(tokens[1], "before", false)) ) {
        std::cout << "Please define a correct species name position:"
            " after or before (the separator)." << std::endl;
      }
      model->set_model_parameter("gene_name_position", tokens[1]);
      mapping_method = MappingMethod::GENE_NAMES;
    } else if ( tokens[0] == "smap" ) {
      smap_file_name = tokens[1];
      mapping_method = MappingMethod::SMAP;
    } else if ( tokens[0] == "delta" ) {
      delta = atof(tokens[1].c_str());
      std::cout << "\n\tDelta fixed to " << delta << std::endl;
    } else if ( tokens[0] == "tau" ) {
      tau = atof(tokens[1].c_str());
      std::cout << "\n\tTau fixed to " << tau << std::endl;
    } else if ( tokens[0] == "lambda" ) {
      lambda = atof(tokens[1].c_str());
      std::cout << "Lambda fixed to " << lambda << std::endl;
    } else if ( tokens[0] == "O_R" ) {
      O_R = atof(tokens[1].c_str());
      std::cout << "\n\tO_R set to " << O_R << std::endl;
    } else if ( tokens[0] == "beta" ) {
      beta = atof(tokens[1].c_str());
      std::cout << "\n\tBeta set to " << beta << std::endl;
    } else if ( tokens[0] == "fraction_missing" ) {
      fractionMissingFile = tokens[1];
      std::cout << "\n\tFile containing fractions of missing genes set to "
                << fractionMissingFile << std::endl;
    } else if ( tokens[0] == "outputFiles" ) {
      if (tokens[1] == "y"
          || tokens[1] == "yes"
          || tokens[1] == "Y"
          || tokens[1] == "YES") {
        outputyn = true;
      }
    }
  }

  GeneMap<std::string, std::string> map;

  if ( mapping_method == MappingMethod::GENE_NAMES ) {
    // before tells where is the species name.
    bool before = true;  // default value
    if ( model->string_parameter.find("gene_name_position" )
       != model->string_parameter.end())
      before = utils::string_comp(
          "before", model->string_parameter.at("gene_name_position"), false);

    map = SpeciesGeneMapper::mapWithGenesNames(
        gene_tree_str, species_tree_str,
        model->string_parameter.at("gene_name_separators"), before);

  } else if ( mapping_method == MappingMethod::TREES ) {
    map = SpeciesGeneMapper::mapWithTrees(gene_tree_str, species_tree_str);

  } else if ( mapping_method == MappingMethod::SMAP ) {
    map = SpeciesGeneMapper::mapFromFile(gene_tree_str, species_tree_str,
                                         smap_file_name);
  }

  // Constructing the ALE_undated object and computing the logLk.
  model->setMap(map);
  model->set_model_parameter("BOOTSTRAP_LABELS", "yes");
  model->construct_undated(species_tree_str, fractionMissingFile);

  // Set model parameters.
  model->set_model_parameter("seq_beta", beta);
  model->set_model_parameter("O_R", O_R);
  // a set of inital rates
  model->set_model_parameter("delta", delta);
  model->set_model_parameter("tau", tau);
  model->set_model_parameter("lambda", lambda);

  // calculate_EGb() must always be called after changing rates to calculate E-s
  // and G-s cf. http://arxiv.org/abs/1211.4606
  model->calculate_undatedEs();
  double loglk = log(model->pun(ale, true));

  std::cout << "\n\tReconciliation model likelihood computed, logLk: "
            << loglk << std::endl;

  // Output
  if ( outputyn ) {
    std::cout << "\n\tSampling reconciled gene trees.." << std::endl;
    std::vector<std::string> sample_strings;
    std::vector <std::shared_ptr<tree_type>> sample_trees;

    Timer<int> sampling_progression(static_cast<int>(samples));

    for ( int i = 0 ; i < samples ; i++ ) {
      // Sample the undated tree.
      std::string sample_tree = model->sample_undated();
      sample_strings.push_back(sample_tree);

      if ( ale->last_leafset_id > 3 ) {
        // Generate the tree.
        std::shared_ptr<tree_type> G(IO::newickToPhyloTree(sample_tree, false));

        // Set a name for each leaf.
        std::vector<Node> leaves = G->getAllLeaves();
        for ( auto it = leaves.begin() ; it != leaves.end() ; it++ ) {
          std::string name = (*it)->getName();
          std::vector<std::string> tokens;
          // boost::split(tokens,name,boost::is_any_of(".@"),
          // boost::token_compress_on);
          tokens = utils::splitString(name, ".@");
          (*it)->setName(tokens[0]);
          tokens.clear();
        }
        leaves.clear();

        // Add the tree into the sample.
        sample_trees.push_back(G);
      }

      sampling_progression.next();
      utils::progressionBar(std::cout, sampling_progression,
                            true, "Sampling trees");  // print progression bar.
    }

    RefreshablePrinter::clean(std::cout);
    std::cout << std::endl;  // endl of progression bar output.

    std::vector<std::string> tokens = utils::splitString(gene_tree_file, "/");

    ale_name = tokens[tokens.size()-1];

    std::string outname = ale_name+".uml_rec";
    std::ofstream fout(outname.c_str());
    fout << "#ALEevaluate using ALE v" << ALE_VERSION
         << "; CC BY-SA 3.0;" << std::endl << std::endl;
    fout << "S:\t" << model->string_parameter["S_with_ranks"] << std::endl;
    fout << std::endl;
    fout << "Gene tree from:\t" << gene_tree_file << std::endl;
    fout << ">logl: " << loglk << std::endl;
    fout << "rate of\t Duplications\tTransfers\tLosses" << std::endl;
    fout << "\t" << delta << "\t" << tau << "\t" << lambda << std::endl;
    fout << std::endl;
    fout << samples << " reconciled G-s:\n" << std::endl;
    for ( int i = 0 ; i < samples ; i++ )
      fout << sample_strings[i] << std::endl;

    fout << "# of\t Duplications\tTransfers\tLosses\tSpeciations" << std::endl;
    fout << "Total \t"<< model->MLRec_events["D"]/samples
         << "\t" << model->MLRec_events["T"]/samples
         << "\t" << model->MLRec_events["L"]/samples
         << "\t" << model->MLRec_events["S"]/samples << std::endl;
    fout << std::endl;
    fout << "# of\t Duplications\tTransfers\tLosses\tOriginations\tcopies"
         << std::endl;
    fout << model->counts_string_undated(samples);
    fout.close();

    // Generate a consensus tree according to sampling.
    std::cout << "Results in: " << outname << std::endl;
    if (ale->last_leafset_id > 3) {
      std::cout << "Calculating consensus tree." << std::endl;
      auto con_tree = PhyloTreeToolBox::thresholdConsensus(sample_trees, 0.5);

      std::string con_name = ale_name+".ucons_tree";

      std::ofstream con_out(con_name.c_str());
      con_out <<  "#ALEsample using ALE v" << ALE_VERSION
              << " by Szollosi GJ et al.; ssolo@elte.hu; CC BY-SA 3.0;"
              << std::endl;
      PhyloTreeToolBox::computeBootstrapValues(*con_tree, sample_trees);
      std::string con_tree_sup = IO::PhyloTreeToNewick(*con_tree);
      con_out << con_tree_sup << std::endl;
      con_out.close();
      std::cout << std::endl << "Consensus tree in " << con_name<< std::endl;
    }

    std::string t_name = ale_name + ".uTs";
    std::ofstream tout(t_name.c_str());
    tout << "#from\tto\tfreq.\n";

    for ( int e = 0 ; e < model->last_branch ; e++ )
      for ( int f = 0 ; f < model->last_branch ; f++ )
        if  ( model->T_to_from[e][f] > 0 )
        {
          if ( e < model->last_leaf )
            tout << "\t" << model->node_name[model->id_nodes[e]];
          else
            tout << "\t" << e;
          if (f < model->last_leaf)
            tout << "\t" << model->node_name[model->id_nodes[f]];
          else
            tout << "\t" << f;
          tout << "\t" << model->T_to_from[e][f]/samples << std::endl;
        }

    tout.close();
    std::cout << "Transfers in: " << t_name << std::endl;
  }

  delete model;

  return 0;
}
