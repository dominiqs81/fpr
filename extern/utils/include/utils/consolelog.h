/**
 * @file consolelog.h
 * @brief Console logging functions
 *
 * @author Domenico Salvagnin <dominiqs at gmail dot com>
 * 2019
 */

#ifndef CONSOLELOG_H
#define CONSOLELOG_H

#include <fmt/format.h>
#include <fmt/color.h>
#include <unistd.h>

namespace dominiqs {

template<typename ...Args>
void consoleLog(Args&&... args)
{
	fmt::print(std::forward<Args>(args)...);
	fmt::print("\n");
}


template<typename ...Args>
void consoleInfo(Args&&... args)
{
	if (isatty(STDOUT_FILENO))  fmt::print(fmt::fg(fmt::color::green), std::forward<Args>(args)...);
	else                        fmt::print(std::forward<Args>(args)...);
	fmt::print("\n");
}


template<typename ...Args>
void consoleWarn(Args&&... args)
{
	if (isatty(STDOUT_FILENO))  fmt::print(fmt::fg(fmt::color::yellow), std::forward<Args>(args)...);
	else                        fmt::print(std::forward<Args>(args)...);
	fmt::print("\n");
}


template<typename ...Args>
void consoleError(Args&&... args)
{
	if (isatty(STDOUT_FILENO))  fmt::print(fmt::fg(fmt::color::red), std::forward<Args>(args)...);
	else                        fmt::print(std::forward<Args>(args)...);
	fmt::print("\n");
}


#ifdef DEBUG_LOG

/* The debug level can be any non-negative integer.
 *
 * The current level is expected to be allocated by the caller (somewhere,
 * either globally or locally) and be named DEBUG_LEVEL
 *
 * The reason the current level is not just a define is so that we can
 * change the value at runtime (e.g., within a debugger).
 */
#define consoleDebug(level, ...)  if (level <= DEBUG_LEVEL) { consoleLog(__VA_ARGS__); }

#else

#define consoleDebug(level, ...)  do {} while (false)

#endif //<DEBUG_LOG

} // namespace dominiqs

#endif /* CONSOLELOG_H */
