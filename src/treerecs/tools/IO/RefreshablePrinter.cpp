//
// Created by comten on 07/09/18.
//

#include "RefreshablePrinter.h"

namespace treerecs {

// RefreshablePrinter definitions
unsigned int RefreshablePrinter::refresh_time_ = 25;

std::chrono::time_point<std::chrono::system_clock> RefreshablePrinter::last_print_time_ = std::chrono::system_clock::now();

std::size_t RefreshablePrinter::previous_message_size_used_ = 0;

void RefreshablePrinter::print(std::ostream &os, std::string message, bool force) {
  if(printable() or force) {
    // Manage overwritings : the new progression bar printed at t+1 may be
    // shorter than the previous one at t.
    // So, artefacts from previous prints could be noisy. We add blanks above to
    // erase previous messages.
    auto message_size = message.size();

    // Erase previous last characters with spaces " " for a clean print
    if (previous_message_size_used_ > message_size) {
      message += std::string(previous_message_size_used_ - message_size, ' ');
      assert(message.size() == previous_message_size_used_);
    }

    message += std::string("\r");
    os << message;
    os.flush();

    previous_message_size_used_ = message_size;

    last_print_time_ = std::chrono::system_clock::now();
  }
}

bool RefreshablePrinter::printable() {
  std::size_t ms_since_last_print =
      static_cast<unsigned int>(
          std::chrono::duration_cast<std::chrono::milliseconds>(
      std::chrono::system_clock::now() - last_print_time_).count());
  return RefreshablePrinter::refresh_time_ < ms_since_last_print;
}

std::ostream &RefreshablePrinter::clean(std::ostream &os) {
  print(os, "", true);
  return os;
}

unsigned int RefreshablePrinter::refresh_time(const unsigned int refresh_time) {
  unsigned int old_refresh_time = refresh_time_;
  refresh_time_ = refresh_time;
  return old_refresh_time;
}

} // namespace treerecs
