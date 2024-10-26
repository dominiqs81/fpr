#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "doctest.h"

#define DEBUG_LOG
#include <utils/consolelog.h>
#include <algorithm>

using namespace dominiqs;

TEST_CASE("check it compiles")
{
	consoleLog("The answer is {}", 42);
	consoleError("The answer is not {}", 40.0);
}

TEST_CASE("test debug log")
{
	int DEBUG_LEVEL = 2;
	consoleDebug(2, "The answer is {}", 42);
	consoleDebug(3, "The answer is not {}", 40.0);
}
