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

#ifndef TREERECS_STATISTICS_H
#define TREERECS_STATISTICS_H

#include <treerecs/containers/ReconciledRootedTree.h>
#include <treerecs/containers/Table.h>

namespace treerecs {

/*!
 * \brief Statistics provides functions to analyse data and Treerecs results.
 */
class Statistics {
private:

protected:

public:

/****************
 * Getters
 */
  ///
  /// \brief Computes the average score of results by threshold.
  /// \details The first column of the resulting Table is the "total_cost" according to the number of duplications and
  ///          losses events in gene tree.
  ///          Formula is :
  ///            total_score = (number_of_duplications * duplication_cost) + (number_of_losses * loss_cost).
  ///
  ///          The second column "ale_loglk" is the ALE loglikelihood evaluation.
  ///          Each row of the Table has the corresponding threshold as name.
  /// \return A table with a first column ("total_cost") the score obtained by the number of duplications and losses, and the second one ("ale_loglk") which contains ALE log likelihood.
  static Table<double, double, std::string> average_scores_by_thresholds(
      const std::unordered_map<std::shared_ptr<bpp::PhyloTree>
          , std::map<double,std::vector<ReconciledRootedTree>>>& input);

  ///
  /// \brief Create class frequencies of values in a given container according to the number of classes to constitute.
  ///   The key of the map is the maximal value of the class (if the minimal value of the container is 0, its indicated
  ///   also the step between two classes). The value associated to the key is the frequency. Sum of frequencies is
  ///   equal to 1.0.
  /// \return A map of frequencies by values classes.
  template<typename Iterator>
  static std::map<typename std::iterator_traits<Iterator>::value_type, double>
  frequencies(
      const Iterator& begin
      , const Iterator& end
      , const long &nclasses
  ) {

    std::map<typename std::iterator_traits<Iterator>::value_type, double> distribution;

    long nvalues = std::distance(begin, end);

    if(nvalues == 0) return distribution;

    // First of all: sort the vector
    std::list<typename std::iterator_traits<Iterator>::value_type> values(
        begin, end);
    values.sort();
    typename std::iterator_traits<Iterator>::value_type min = values.front();
    typename std::iterator_traits<Iterator>::value_type max = values.back();

    auto n_resulting_classes = static_cast<decltype(min)>((nclasses < nvalues) ?
                                                          nclasses : nvalues);

    double range = static_cast<double>(max-min);

    double step = range / n_resulting_classes;


    double class_i = min + step;
    distribution[class_i] = 0.0;
    for (auto value: values) {
      while (not utils::double_equal_or_inferior(value, class_i)) {
        class_i += step;
        distribution[class_i] = 0.0;

      }
      distribution[class_i] += 1.0;
    }

    for (auto &pair: distribution) {
      pair.second /= static_cast<double>(nvalues);
    }

    return distribution;
  }

  ///
  /// \brief Barplot in ASCII of frequencies. Should be used with Statistics::frequencies.
  /// \brief Plot has this visual:
  ///
  ///     [0.00 - 0.16] :  **********           (12.12%)
  ///     [0.16 - 0.29] :  *********            (10.60%)
  ///     [0.29 - 0.43] :  ******************** (22.72%)
  ///     [0.43 - 0.56] :  ********             ( 9.09%)
  ///     [0.56 - 0.70] :  ****                 ( 4.54%)
  ///     [0.70 - 0.83] :  ****************     (18.18%)
  ///     [0.83 - 0.97] :  ******************** (22.72%)
  ///
  template<class T>
  static void plotFrequencies(
      const std::map<T, double> freq /// Frequencies.
      , std::ostream &os /// Output stream where the barplot is going to be printed.
  ) {
    // Prints parameters
    int barlength = MAXIMAL_PLOT_WIDTH;
    std::size_t maxLabelSize = 0;
    std::size_t maxFreqPrintSize = 0;

    // During freq exploration, we will looking for the maximum value.
    double maxFreq = 0.0;

    // Collect keys of the freq map
    std::vector<T> freq_keys;
    freq_keys.resize(freq.size());

    // Collect their values print, as intervals.
    std::vector<std::string> value_interval_prints;
    value_interval_prints.reserve(freq.size());

    std::size_t i = 0;
    for(auto& pair: freq){
      std::string value_interval_print =  "[";
      std::string born0 = (i == 0 ? std::to_string((T)0.0) : std::to_string(freq_keys.back()));
      born0 = utils::trunc_string_number(born0);
      value_interval_print+= born0 + " - ";
      std::string born1 = std::to_string(pair.first);
      born1 = utils::trunc_string_number(born1);
      value_interval_print+= born1 + "]";

      std::size_t label_size = utils::getStreamObjectSize(value_interval_print, os);

      if (maxLabelSize < label_size)
        maxLabelSize = label_size;

      if (pair.second > maxFreq)
        maxFreq = pair.second;

      value_interval_prints.push_back(value_interval_print);
      freq_keys.push_back(pair.first);
      i++;
    }

    maxFreqPrintSize = utils::trunc_string_number(std::to_string(maxFreq * 100.0)).size();
    i = 0;
    for (auto &pair: freq) {
      //pair.first is the class
      //pair.second is the frequency
      // A print shows like that:
      // class_label : ******      ( frequency_print%)

      auto& class_label = value_interval_prints.at(i);
      auto& frequency = pair.second;

      auto relative_frequency = frequency/ maxFreq;

      auto bar_length_used = std::size_t(relative_frequency * (double)barlength);
      auto frequency_print = utils::trunc_string_number(std::to_string(frequency * 100.0));

      os << class_label
         << std::string(maxLabelSize - utils::getStreamObjectSize(class_label, os) + 1, ' ')
         << ":  " << std::string(bar_length_used, '*')
         << std::string(barlength - bar_length_used, ' ')
         << " (" << std::string(maxFreqPrintSize - frequency_print.size(), ' ') << frequency_print << "%)"
         << std::endl;
      i++;
    }
  }
};

} // namespace treerecs

#endif //TREERECS_STATISTICS_H
