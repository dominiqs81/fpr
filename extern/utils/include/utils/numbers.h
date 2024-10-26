/**
 * @file numbers.h
 * @brief Numerical utilities Header
 *
 * @author Domenico Salvagnin dominiqs@gmail.com
 * 2008
 */

#ifndef NUMBERS_H
#define NUMBERS_H

#include <algorithm>
#include <type_traits>
#include <iostream>

#include <cassert>

namespace dominiqs {

/**
 * returns the smallest number eps > 0, for which 1.0 + eps > 1.0 on the current machine
 */

template<typename FloatType=double>
FloatType getMachineEps()
{
	static_assert(std::is_floating_point_v<FloatType>);
	FloatType one = 1.0;
	FloatType eps = 1.0;
	FloatType res;
	do {
		res = eps;
		eps /= 2.0;
	}
	while ((one + eps) > one);
	return res;
}

/**
 * calculate the greatest common divisor of two integral types n and m
 */

template <typename IntType> IntType gcd(IntType n, IntType m)
{
	static_assert(std::is_integral_v<IntType>);
	IntType zero(0);
	if (n < zero) n = -n;
	if (m < zero) m = -m;
	// check for zeros
	if (m == zero) return n;
	if (n == zero) return m;
	assert( n > zero );
	assert( m > zero );
	// binary Stein algorithm
	IntType const one(1);
	// the greatest divisors of m, n that are powers of 2
	IntType mOdds(zero);
	IntType nOdds(zero);
	IntType const mask(one); //< binary mask to check for % 2
	// n,m stripped of trailing zeros
	for (; !(m & mask) ; m >>= 1, mOdds++);
	for (; !(n & mask) ; n >>= 1, nOdds++);
	for (; n != m ;) {
		// n == m is gcd
		if (n < m) for (m -= n, m >>= 1; !(m & mask) ; m >>= 1); //< strips trailing zeros from m
		else       for (n -= m, n >>= 1; !(n & mask) ; n >>= 1);
	}
	return m << (std::min)( mOdds, nOdds );
}

/**
 * calculates the least common multiple of two integral types n and m
 */

template <typename IntType> IntType lcm(IntType n, IntType m)
{
	static_assert(std::is_integral_v<IntType>);
	return n / gcd<IntType>(n, m) * m;
}

/**
 * Calculate the binomial coefficient (n choose m)
 * @return the binomial coefficient
 */

template<typename IntType>
double choose(IntType m, IntType n)
{
	static_assert(std::is_integral_v<IntType>);
	if (n > m) return 0;
	if (n == m) return 1;
	double r = 1;
	for (IntType d = 1; d <= n; d++) {
		r *= m--;
		r /= d;
	}
	return r;
}

} // namespace dominiqs

#endif /* NUMBERS_H */
