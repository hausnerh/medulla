/**
 * @file yell_try_die.h
 * @brief Process controll to make it a little easier to do error/process control
 * @author hhausner@fnal.gov
 */

#ifndef YELL_TRY_DIE_H
#define YELL_TRY_DIE_H

#include <iostream>
#include <functional>

/**
 * @brief Yell and error message
 */
void yell(const std::string& message)
{
  std::cerr << message << std::endl;
}

/**
 * @brief thow a runtime error with a descriptive failure
 */
[[noreturn]] void die(const std::string& message)
{
  throw std::runtime_error(message);
}

/**
 * @brief Attempt to call a function with some arguments, die if it fails.
 */
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