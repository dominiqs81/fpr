#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "doctest.h"

#include <utils/index_partition.h>

using namespace dominiqs;


TEST_CASE("index partition empty")
{
	int n = 5; //< # indices
	int ns = 2; //< # sets
	IndexPartition part(n, ns);
	CHECK( part.getNumIndices() == n );
	CHECK( part.getNumSets() == ns );
	for (int i = 0; i < n; i++) {
		CHECK( part.isInNoSet(i) );
	}
	for (int s = 0; s < ns; s++) {
		CHECK( part.isEmpty(s) );
		CHECK( part.begin(s) == part.end(s) );
	}
}


TEST_CASE("index partition basics")
{
	int n = 5; //< # indices
	int ns = 2; //< # sets
	IndexPartition part(n, ns);
	part.add(0, 0);
	part.add(1, 1);
	part.add(2, 0);
	part.add(3, 1);
	part.add(4, 0);
	CHECK( part.begin(0) != part.end(0) );
	CHECK( part.begin(1) != part.end(1) );
	std::vector<int> set0;
	std::copy(part.begin(0), part.end(0), std::back_inserter(set0));
	CHECK( set0.size() == 3 );
	CHECK( set0[0] == 4 );
	CHECK( set0[1] == 2 );
	CHECK( set0[2] == 0 );
	std::vector<int> set1;
	std::copy(part.begin(1), part.end(1), std::back_inserter(set1));
	CHECK( set1.size() == 2 );
	CHECK( set1[0] == 3 );
	CHECK( set1[1] == 1 );
	// remove
	part.remove(2);
	set0.clear();
	std::copy(part.begin(0), part.end(0), std::back_inserter(set0));
	CHECK( set0.size() == 2 );
	CHECK( set0[0] == 4 );
	CHECK( set0[2] == 0 );
	// clear
	part.clear();
	for (int i = 0; i < n; i++) {
		CHECK( part.isInNoSet(i) );
	}
	for (int s = 0; s < ns; s++) {
		CHECK( part.isEmpty(s) );
		CHECK( part.begin(s) == part.end(s) );
	}
}
