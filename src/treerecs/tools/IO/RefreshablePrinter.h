//
// Created by comten on 07/09/18.
//

#ifndef TREERECS_REFRESHABLEPRINTER_H
#define TREERECS_REFRESHABLEPRINTER_H

#include "time.h"

#include <iostream>
#include <cmath>
#include <chrono>
#include <cassert>

namespace treerecs {

/*!
 * @class RefreshablePrinter
 * @brief allows clean prints by saving the size of the last element printed.
 */
class RefreshablePrinter {
public:
  /// Clean the current terminal line by erasing each character by a space.
  static std::ostream& clean(std::ostream& os);

  /// Stream the message in an std::ostream.
  static void print(std::ostream& os, std::string message, bool force = false);

  /// Check if refresh is possible according to the refresh frequency.
  static bool printable();

  /// Change refreshing time between two prints.
  static unsigned int refresh_time(const unsigned int value);

protected:
  /// Indicates the size of the last element printed by a RefreshablePrinter.
  static std::size_t previous_message_size_used_;
  /// Time of the last print.
  static std::chrono::time_point<std::chrono::system_clock> last_print_time_;
  /// Refreshing time in ms.
  static unsigned int refresh_time_;
};

} // namespace treerecs

#endif //TREERECS_REFRESHABLEPRINTER_H
