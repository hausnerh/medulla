#include <iostream>
#include <functional>

#ifndef YELL_TRY_DIE_H
#define YELL_TRY_DIE_H

void yell(const std::string& message)
{
  std::cerr << message << std::endl;
}

[[noreturn]] void die(const std::string& message)
{
  throw std::runtime_error(message);
}

template <typename Func, typename... Args>
  auto try_call(const std::string& description, Func&& function, Args&&... arguments) ->
    decltype(std::invoke(std::forward<Func>(function), std::forward<Args>(arguments)...))
    {
      try
      {
        return std::invoke(std::forward<Func>(function), std::forward<Args>(arguments)...);
      }
      catch (const std::exception& e)
      {
        std::string msg = "cannot " + description + ": " + e.what();
        die(msg);
      }
      catch (...)
      {
        std::string msg = "cannot " + description + ": Unknown Error";
        die(msg);
      }
    }

#endif // YELL_TRY_DIE_H
