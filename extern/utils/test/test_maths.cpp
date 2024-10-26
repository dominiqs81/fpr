#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "doctest.h"

#include <utils/maths.h>
#include <algorithm>

using namespace dominiqs;

TEST_CASE("SparseVector basic operations")
{
	SparseVector sv;

	CHECK( sv.empty() );
	CHECK( sv.size() == 0 );
	CHECK( sv.idx() == nullptr );
	CHECK( sv.coef() == nullptr );

	SUBCASE("push-pop") {
		sv.push(2, 1.0);
		CHECK( !sv.empty() );
		CHECK( sv.size() == 1 );
		CHECK( sv.idx() != nullptr );
		CHECK( sv.coef() != nullptr );
		CHECK( sv.idx()[0] == 2 );
		CHECK( sv.coef()[0] == 1.0 );

		sv.pop();
		CHECK( sv.empty() );
		CHECK( sv.size() == 0 );
		CHECK( sv.idx() == nullptr );
		CHECK( sv.coef() == nullptr );
	}
}


TEST_CASE("SparseVector constructors and operators")
{
	SparseVector sv = {{1, 3.0}, {4, -2.0}};

	CHECK( sv.size() == 2 );
	CHECK( sv.idx() != nullptr );
	CHECK( sv.coef() != nullptr );

	const int idx[] = {1, 4};
	const double coef[] = {3.0, -2.0};
	CHECK( std::equal(sv.idx(), sv.idx() + sv.size(), idx) );
	CHECK( std::equal(sv.coef(), sv.coef() + sv.size(), coef) );

	SparseVector sv2;
	sv2.push(1, 3.0);
	sv2.push(4, -2.0);
	CHECK( sv == sv2 );
	CHECK( sv.view() == sv2.view() );
	sv2.push(5, 1.0);
	CHECK( sv != sv2 );
	CHECK( sv.view() != sv2.view() );

	SparseVector sv3 = as_span({0.0, 3.0, 0.0, 0.0, -2.0});
	CHECK( sv == sv3 );
	CHECK( sv.view() == sv3.view() );
}


TEST_CASE("SparseVector advanced")
{
	const double dense[] = {1.0, 0.0, 0.0, -3.0, 5.4};
	SparseVector sv;

	// gather
	sv.gather(dense);
	CHECK( !sv.empty() );
	CHECK( sv.size() == 3 );
	CHECK( sv.idx() != nullptr );
	CHECK( sv.coef() != nullptr );
	const int idx[] = {0, 3, 4};
	const double coef[] = {1.0, -3.0, 5.4};
	CHECK( std::equal(sv.idx(), sv.idx() + sv.size(), idx) );
	CHECK( std::equal(sv.coef(), sv.coef() + sv.size(), coef) );

	// scatter
	double dense2[] = {0.0, 0.0, 0.0, 0.0, 0.0};
	sv.scatter(dense2);
	CHECK( std::equal(dense, dense+5, dense2) );

	// unscatter
	sv.unscatter(dense2);
	CHECK( std::all_of(dense2, dense2+5, [](double x){ return x == 0.0;}) );

	// scale
	sv.scale(0.5);
	const double coef2[] = {0.5, -1.5, 2.7};
	CHECK( std::equal(sv.idx(), sv.idx() + sv.size(), idx) );
	CHECK( std::equal(sv.coef(), sv.coef() + sv.size(), coef2) );

	// negate
	sv.negate();
	const double coef3[] = {-0.5, 1.5, -2.7};
	CHECK( std::equal(sv.idx(), sv.idx() + sv.size(), idx) );
	CHECK( std::equal(sv.coef(), sv.coef() + sv.size(), coef3) );
}


TEST_CASE("SparseMatrix")
{
	SparseMatrix sm;
	CHECK( sm.k == 0 );
	CHECK( sm.U == 0 );
	CHECK( sm.nnz == 0 );

	SUBCASE("transpose empty") {
		SparseMatrix transposed = sm.transpose();
		CHECK( sm == transposed );
	}
	SUBCASE("transpose: no rows but cols") {
		sm.U = 3;
		CHECK( sm.k == 0 );
		CHECK( sm.U == 3 );
		CHECK( sm.nnz == 0 );
		SparseMatrix transposed = sm.transpose();
		CHECK( transposed.k == 3 );
		CHECK( transposed.U == 0 );
		CHECK( transposed.nnz == 0 );
		for (int i = 0; i < 3; i++) {
			CHECK( transposed[i].size() == 0 );
		}
	}
	SUBCASE("transpose: rows but no cols") {
		for (int i = 0; i < 3; i++)  sm.add(SparseVector{});
		CHECK( sm.k == 3 );
		CHECK( sm.U == 0 );
		CHECK( sm.nnz == 0 );
		SparseMatrix transposed = sm.transpose();
		CHECK( transposed.k == 0 );
		CHECK( transposed.U == 3 );
		CHECK( transposed.nnz == 0 );
	}
	SUBCASE("3x3 identity") {
		sm.U = 3;
		for (int i = 0; i < 3; i++) {
			sm.add(SparseVector{{i, 1.0}});
		}
		CHECK( sm.k == 3 );
		CHECK( sm.U == 3 );
		CHECK( sm.nnz == 3 );
		for (int i = 0; i < 3; i++) {
			CHECK( sm[i] == SparseVector{{i, 1.0}} );
		}

		SparseMatrix transposed = sm.transpose();
		CHECK( sm == transposed );
	}
	SUBCASE("3x3 identity again") {
		SparseMatrix sm2 = {
			{1.0, 0.0, 0.0},
			{0.0, 1.0, 0.0},
			{0.0, 0.0, 1.0}
		};
		CHECK( sm2.k == 3 );
		CHECK( sm2.U == 3 );
		CHECK( sm2.nnz == 3 );
		for (int i = 0; i < 3; i++) {
			CHECK( sm2[i] == SparseVector{{i, 1.0}} );
		}

		SparseMatrix transposed = sm2.transpose();
		CHECK( sm2 == transposed );
	}
	SUBCASE("2x4") {
		SparseMatrix sm = {
			{1.0, 2.0, 0.0, 0.0},
			{0.0, 1.0, 1.0, 3.0}
		};
		CHECK( sm.k == 2 );
		CHECK( sm.U == 4 );
		CHECK( sm.nnz == 5 );
	
		SparseMatrix transposed = sm.transpose();
		CHECK( sm.k == transposed.U );
		CHECK( sm.U == transposed.k );
		CHECK( sm.nnz == transposed.nnz );
		SparseMatrix expected = {
			{1.0, 0.0},
			{2.0, 1.0},
			{0.0, 1.0},
			{0.0, 3.0},
		};
		CHECK( expected == transposed );
	}
}
