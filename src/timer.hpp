#ifndef TIMER_HPP
#define TIMER_HPP

#include <chrono>
#include <iostream>
#include <string>
#include <unordered_map>

namespace naiad
{

class Timer_element
{
  public:

    Timer_element() {}
    Timer_element(const std::string & name_) : name{name_} {}

    void start();
    void stop();

    using Clock = std::chrono::steady_clock;
    using Point = std::chrono::time_point<Clock>;
    using Duration = std::chrono::duration<Clock::rep, Clock::period>;

    std::string name;
    double elapsed{0.0};
    Point time_start{Clock::now()};
    Duration time_elapsed{0};
    bool running{false};
};

class Timer
{
  public:

    Timer();
    void start(const std::string & name);
    void stop(const std::string & name);
    void summary(std::ostream & os);

  private:

    std::unordered_map<std::string, Timer_element> timers;
};

} // namespace naiad

#endif
