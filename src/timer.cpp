#include "timer.hpp"

#include <algorithm>
#include <vector>

namespace naiad
{

void Timer_element::start()
{
  time_start = Clock::now();
  running = true;
}

void Timer_element::stop()
{
  if (running)
  {
    time_elapsed += Clock::now() - time_start;
    elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(time_elapsed).count() * 1e-3;
  }
  running = false;
}

Timer::Timer() { start("total"); }

void Timer::start(const std::string & name)
{
  const auto find{timers.find(name)};
  if (find == timers.end())
  {
    const auto tt{timers.insert({name, Timer_element{name}})};
    if (tt.second)
      tt.first->second.start();
    else
    {
      std::cerr << "failed to insert timer: " << name << std::endl;
      std::abort();
    }
  }
  else
  {
    find->second.start();
  }
}

void Timer::stop(const std::string & name)
{
  const auto find{timers.find(name)};
  if (find != timers.end())
    find->second.stop();
  else
  {
    std::cerr << "failed to find and stop timer: " << name << std::endl;
    std::abort();
  }
}

void Timer::summary(std::ostream & os)
{
  for (auto & t : timers)
    t.second.stop();

  // I make a copy so that I can sort them in the output
  std::vector<Timer_element> tsort;
  tsort.reserve(timers.size());
  for (const auto & t : timers)
    tsort.emplace_back(t.second);

  constexpr auto sort_runtime{[](const Timer_element & t1, const Timer_element & t2)
                              { return t1.elapsed > t2.elapsed; }};
  std::sort(tsort.begin(), tsort.end(), sort_runtime);

  const double total{timers.at("total").elapsed};

  std::size_t max_str_len{0};
  for (const auto & t : tsort)
    max_str_len = std::max(max_str_len, t.name.length());

  os << "=== TIMING SUMMARY ===" << std::endl;
  // TODO maybe put this in some sort of order?
  // TODO appropriate string padding
  os << "  Name " << std::string(max_str_len - 3, ' ') << " Elapsed [s]    Fraction [%] " << std::endl;
  os << " " << std::string(max_str_len + 2, '-') << "  -------------  --------------" << std::endl;
  for (const auto & t : tsort)
    os << "  " << std::string(max_str_len - t.name.length(), ' ') << t.name << "      "
       << ((t.elapsed > 1e3) ? std::format("{:8.2e}", t.elapsed) : std::format("{:8.3f}", t.elapsed)) << "        "
       << std::format("{:6.2f}%", t.elapsed / total * 1e2) << std::endl;
  os << std::endl;
}

} // namespace naiad
