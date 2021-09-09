//
// Created by comten on 07/09/18.
//

#ifndef TREERECS_TIMER_H
#define TREERECS_TIMER_H

#include "time.h"

#include <iostream>
#include <cmath>
#include <chrono>
#include <cassert>

namespace treerecs {

template <typename T>
class Timer {
private:
  /// Step of the iteration.
  T step_;
  /// Current value of the iteration.
  T current_;
  /// Ending value of the iteration.
  T max_;
  /// Starting value.
  T min_;
  /// Starting clocks.
  clock_t starting_clocks_;
  /// Current clocks.
  clock_t current_clocks_;
  /// Starting time.
  std::chrono::time_point<std::chrono::system_clock> starting_time_;
  /// Current time.
  std::chrono::time_point<std::chrono::system_clock> current_time_;

protected:

public:
  /****************
   * Constructors
   */
  explicit Timer(
      const T& max /// Ending value of the iterator.
      , const T& min = 0 /// Starting value of the iterator (default = 0).
      , const T& step = 1 /// Increasing (or decreasing) value of the iterator (default = 1).
  ) : step_(step)
      , current_(min)
      , max_(max)
      , min_(0)
      , starting_clocks_(clock())
      , current_clocks_(clock())
      , starting_time_(std::chrono::system_clock::now())
      , current_time_(std::chrono::system_clock::now())
  {}

  /****************
   * Getters
   */
  /// Current state of the timer
  T current() const;

  /// Maximum state of the timer
  T max() const;

  /// Minimum state of the timer
  T min() const;

  /// Returns the total number of steps recorded.
  std::size_t total_number_of_steps() const;

  /// Returns the average clocks spent for one step (in seconds).
  double average_clocks_per_step() const;

  /// Returns the average time spent for one step (in seconds).
  double average_time_per_step() const;

  /// Returns the total time spent until the last update (in seconds).
  double time_past_since_last_update() const;

  /// Returns the total clocks spent until the last update (in seconds).
  double clocks_past_since_last_update() const;

  /// Returns the total time spent since the beginning (in seconds).
  double time_past() const;

  /// Returns the total clocks spent since the beginning (in seconds).
  double clocks_past() const;

  /// Returns remaining steps until the end.
  T remaining_steps() const;

  /// Returns remaining estimated time until the end (in seconds).
  double remaining_time() const;

  /// Returns remaining estimated clocks until the end (in seconds).
  double remaining_clocks() const;

  /// Returns current progression (from 0 to 1).
  double progress() const;

  /// Indicates if the loop has terminated.
  bool ended(void) const;

  /****************
   * Setters
   */

  /// Change maximal value
  void setMaximum(const T& max);

  /// Update the current state in the progression.
  void update(const T &i);

  /// Increases the progression to one step.
  void next(void);

  /// Increases the progression to n steps.
  void next(const T& steps);
};

// Progression definitions
template <typename T>
T Timer<T>::current() const {
  return current_;
}

template <typename T>
T Timer<T>::min() const {
  return min_;
}

template <typename T>
T Timer<T>::max() const {
  return max_;
}

template <typename T>
std::size_t Timer<T>::total_number_of_steps() const {
  return (std::size_t) floor((double)(max_ - min_)/(double)step_);
}

template <typename T>
double Timer<T>::average_clocks_per_step() const {
  return clocks_past_since_last_update() / ((double)(current_) - min_);
}

template <typename T>
double Timer<T>::average_time_per_step() const {
  return time_past_since_last_update() / ((double)(current_) - min_);
}

template <typename T>
double Timer<T>::time_past_since_last_update() const {
  return 0.001 * (double)std::chrono::duration_cast<std::chrono::milliseconds>(current_time_ - starting_time_).count();
}

template <typename T>
double Timer<T>::clocks_past_since_last_update() const {
  return (double) (current_clocks_ - starting_clocks_)/CLOCKS_PER_SEC;
}

template <typename T>
double Timer<T>::time_past() const {
  return 0.001 *
         (double) std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now() - starting_time_).count();
}

template <typename T>
double Timer<T>::clocks_past() const {
  return (double) (clock() - starting_clocks_)/CLOCKS_PER_SEC;
}

template <typename T>
T Timer<T>::remaining_steps() const {
  return (max_ - current_)/step_;
}

template <typename T>
double Timer<T>::remaining_time() const {
  //return (double) time_past_since_last_update() * (1.0 - progress())/progress();
  return (double) remaining_steps() * average_time_per_step();
}

template <typename T>
double Timer<T>::remaining_clocks() const {
  //return (double) clocks_past_since_last_update() * (1.0 - progress())/progress();
  return (double) remaining_steps() * average_clocks_per_step();
}

template <typename T>
double Timer<T>::progress() const {
  double p = (double) (current_ - min_) / (double) (max_ - min_) ;
  return p;
}

template <typename T>
bool Timer<T>::ended(void) const { return current_ >= max_; }

template <typename T>
void Timer<T>::setMaximum(const T &max) { max_ = max; }

template <typename T>
void Timer<T>::update(const T &i) {
  this->current_ = i;
  this->current_clocks_ = clock();
  this->current_time_ = std::chrono::system_clock::now();
}

template <typename T>
void Timer<T>::next(void) {
  this->update(this->current_ + this->step_);
  this->current_clocks_ = clock();
  this->current_time_ = std::chrono::system_clock::now();
}

template <typename T>
void Timer<T>::next(const T &steps) {
  this->update(this->current_ + steps);
  this->current_clocks_ = clock();
  this->current_time_ = std::chrono::system_clock::now();
}

} // namespace treerecs

#endif //TREERECS_TIMER_H
