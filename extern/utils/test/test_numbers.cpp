#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "doctest.h"

#include <utils/numbers.h>

using namespace dominiqs;

const double eps = 1e-4;

TEST_CASE("gcd")
{
	CHECK( gcd<int>(2, 3) == 1 );
	CHECK( gcd<int>(2, 4) == 2 );
	CHECK( gcd<int>(-2, 3) == 1 );
	CHECK( gcd<int>(-2, 4) == 2 );
}
