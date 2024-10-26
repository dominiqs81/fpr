#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "doctest.h"

#include <utils/floats.h>

using namespace dominiqs;

const double eps = 1e-4;

TEST_CASE("float comparisons")
{
	for (double x: {-10000.0, -1.5, 0.0, +1.5, +10000.0}) {
		double y = x + 2.0 * eps;
		CHECK( !equal(x, y, eps) );
		CHECK( different(x, y, eps) );
		CHECK( lessThan(x, y, eps) );
		CHECK( lessEqualThan(x, y, eps) );
		CHECK( !greaterThan(x, y, eps) );
		CHECK( !greaterEqualThan(x, y, eps) );

		double z = x + 0.5 * eps;
		CHECK( equal(x, z, eps) );
		CHECK( !different(x, z, eps) );
		CHECK( !lessThan(x, z, eps) );
		CHECK( lessEqualThan(x, z, eps) );
		CHECK( !greaterThan(x, z, eps) );
		CHECK( greaterEqualThan(x, z, eps) );
	}
}


TEST_CASE("float sign")
{
	CHECK( isNull(eps / 2, eps) );
	CHECK( isNull(-eps / 2, eps) );
	CHECK( !isNotNull(eps / 2, eps) );
	CHECK( !isNotNull(-eps / 2, eps) );
	CHECK( !isPositive(eps / 2, eps) );
	CHECK( !isPositive(-eps / 2, eps) );
	CHECK( !isNegative(eps / 2, eps) );
	CHECK( !isNegative(-eps / 2, eps) );
	CHECK( !isNull(2 * eps, eps) );
	CHECK( !isNull(-2 * eps, eps) );
	CHECK( isNotNull(2 * eps, eps) );
	CHECK( isNotNull(-2 * eps, eps) );
	CHECK( isPositive(2 * eps, eps) );
	CHECK( !isPositive(-2 * eps, eps) );
	CHECK( !isNegative(2 * eps, eps) );
	CHECK( isNegative(-2 * eps, eps) );
}


TEST_CASE("float round")
{
	double x = 10.5;
	CHECK( floorEps(x, eps) == 10.0 );
	CHECK( ceilEps(x, eps) == 11.0 );
	CHECK( fractionalPart(x, eps) == 0.5 );
	CHECK( integralityViolation(x, eps) == 0.5 );
	CHECK( !isInteger(x, eps) );

	double half = eps / 2;
	x = 10.0 - half;
	CHECK( floorEps(x, eps) == 10.0 );
	CHECK( ceilEps(x, eps) == 10.0 );
	CHECK( fractionalPart(x, eps) <= half );
	CHECK( integralityViolation(x, eps) <= half );
	CHECK( isInteger(x, eps) );

	x = -10.0 + half;
	CHECK( floorEps(x, eps) == -10.0 );
	CHECK( ceilEps(x, eps) == -10.0 );
	CHECK( fractionalPart(x, eps) <= half );
	CHECK( integralityViolation(x, eps) <= half );
	CHECK( isInteger(x, eps) );
}


TEST_CASE("float to int")
{
	CHECK( double2int<int>(0.0) == 0 );
	CHECK( double2int<int>(eps / 2) == 0 );
	CHECK( double2int<int>(-eps / 2) == 0 );
	CHECK( double2int<int>(10.1) == 10 );
	CHECK( double2int<int>(10.7) == 11 );
	CHECK( double2int<int>(-10.1) == -10 );
	CHECK( double2int<int>(-10.7) == -11 );
}
