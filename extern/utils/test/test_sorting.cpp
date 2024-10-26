#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "doctest.h"

#include <utils/sorting.h>
#include <vector>

using namespace dominiqs;

TEST_CASE("shellsort")
{
	std::vector<int> vi = {4, -3, 2, -9};
	std::vector<int> sorted_vi = {-9, -3, 2, 4};

	SUBCASE("int as vector") {
		shellSort(vi);
		CHECK( std::equal(vi.begin(), vi.end(), sorted_vi.begin()) );
	}

	SUBCASE("int as pointer") {
		shellSort(vi.data(), vi.size());
		CHECK( std::equal(vi.begin(), vi.end(), sorted_vi.begin()) );
	}

	std::vector<double> vd = {4.3, -3.5, 2.7, -9.0};
	std::vector<double> sorted_vd = {-9.0, -3.5, 2.7, 4.3};

	SUBCASE("double as vector") {
		shellSort(vd);
		CHECK( std::equal(vd.begin(), vd.end(), sorted_vd.begin()) );
	}

	SUBCASE("double as pointer") {
		shellSort(vd.data(), vd.size());
		CHECK( std::equal(vd.begin(), vd.end(), sorted_vd.begin()) );
	}
}


TEST_CASE("permshellsort")
{
	std::vector<int> vi = {4, -3, 2, -9};
	std::vector<int> perm(vi.size());
	std::vector<int> sorted_vi = {-9, -3, 2, 4};

	permShellSort(vi.data(), perm.data(), vi.size(), std::less<int>() );
	permute(vi.data(), perm.data(), vi.size() );
	CHECK( std::equal(vi.begin(), vi.end(), sorted_vi.begin()) );

	std::vector<double> vd = {4.3, -3.5, 2.7, -9.0};
	std::vector<double> sorted_vd = {-9.0, -3.5, 2.7, 4.3};

	permShellSort(vd.data(), perm.data(), vd.size(), std::less<int>() );
	permute(vd.data(), perm.data(), vd.size() );
	CHECK( std::equal(vd.begin(), vd.end(), sorted_vd.begin()) );
}


TEST_CASE("insertionsort")
{
	std::vector<std::pair<int,int>> v = {{0,4}, {0,-3}, {0,2}, {0,-9}};
	std::vector<std::pair<int,int>> sorted_v = v;

	insertionSort(v.data(), v.size(), [](const auto& x, const auto& y){
		return (x.first < y.first);
	});
	/* insertion sort is stable! */
	CHECK( std::equal(v.begin(), v.end(), sorted_v.begin()) );
}
