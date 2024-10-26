#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "doctest.h"

#include <utils/fileconfig.h>

using namespace dominiqs;

TEST_CASE("load cfg file")
{
	FileConfig cfg;
	cfg.load("data/simple.cfg");

	CHECK( cfg.get<std::string>("name",  "") == "SimpleTest" );
	CHECK( cfg.get<double>("numerical.tolerance", 0.0) == 1e-6 );
	CHECK( cfg.get<int>("res.limit", -1) == 10000 );

	SUBCASE("load with merge") {
		cfg.load("data/overload.cfg");
		CHECK( cfg.get<std::string>("name",  "") == "OverloadedTest" );
		CHECK( cfg.get<double>("numerical.tolerance", 0.0) == 1e-3 );
		CHECK( cfg.get<int>("res.limit", -1) == 10000 );
		CHECK( cfg.get<int>("res.otherLimit", -1) == 20000 );
	}
}
