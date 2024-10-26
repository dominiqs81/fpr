/**
 * @file maths.h
 * @brief Mathematical utilities Header
 *
 * @author Domenico Salvagnin dominiqs@gmail.com
 */

#ifndef MATHS_H
#define MATHS_H

#include <cmath>
#include <vector>
#include <span>
#include <string>
#include <memory>
#include <functional>

#include <cassert>
#include "floats.h"


namespace dominiqs {

/* Small helpers to be able to pass initializer_lists
 * to functions (including constructors) that
 * expect an std::span, as this does not provide
 * the necessary constructor
 */
template <class Pointer>
constexpr auto as_span_impl(Pointer* p, std::size_t s) { return std::span<Pointer>(p, s); }

template <class Container>
auto as_span(Container& c) { return as_span_impl(std::data(c), std::size(c)); }

template <class Container>
auto as_span(const Container& c) { return as_span_impl(std::data(c), std::size(c)); }

template <class T, int N>
constexpr auto as_span(T (&arr)[N]) noexcept { return as_span_impl(arr, N); }

template <class T, int N>
constexpr auto as_span(const T (&arr)[N]) noexcept { return as_span_impl(arr, N); }


/** Sparse Vector Iterator */
template<typename IdxType, typename CoefType>
class SparseVectorIteratorType
{
public:
	using iterator_category = std::random_access_iterator_tag;
	using difference_type   = std::ptrdiff_t;
	using value_type        = std::pair<IdxType,CoefType>;
	using const_reference   = std::pair<const IdxType&, const CoefType&>;
	// Constructor
	SparseVectorIteratorType(const IdxType* idx, const CoefType* coef) : idx_ptr(idx), coef_ptr(coef) {}
	// Derefence operator
	const_reference operator*() const { return const_reference(*idx_ptr, *coef_ptr); }
	// Prefix increment
	SparseVectorIteratorType& operator++() { idx_ptr++; coef_ptr++; return *this; }  
	// Postfix increment
	SparseVectorIteratorType operator++(int) { SparseVectorIteratorType tmp = *this; ++(*this); return tmp; }
	// Operators
	bool operator==(const SparseVectorIteratorType& other) const { return (idx_ptr == other.idx_ptr) && (coef_ptr == other.coef_ptr); }
	bool operator!=(const SparseVectorIteratorType& other) const { return !(*this == other); }
private:
	const IdxType* idx_ptr;
	const CoefType* coef_ptr;
};


/** Sparse Vector View
 *
 * Non-owning view on a sparse vector
 */
template<typename IdxType, typename CoefType>
class SparseVectorViewType
{
public:
	SparseVectorViewType(const IdxType* idx, const CoefType* coef, std::size_t n) : indices(idx), coefs(coef), count(n) {}
	using iterator = SparseVectorIteratorType<IdxType,CoefType>;
	iterator begin() const { return iterator(indices, coefs); }
	iterator end() const { return iterator(indices+count, coefs+count); }
	const IdxType* idx() const { return indices; }
	const CoefType* coef() const { return coefs; }
	std::size_t size() const { return count; }
	// Operators
	bool operator==(const SparseVectorViewType& other) const
	{
		if (size() != other.size())  return false;
		if (!std::equal(idx(), idx() + size(), other.idx()))  return false;
		if (!std::equal(coef(), coef() + size(), other.coef()))  return false;
		return true;
	}
	bool operator!=(const SparseVectorViewType& other) const
	{
		return !(*this == other);
	}
private:
	const IdxType* indices;
	const CoefType* coefs;
	const std::size_t count;
};

/**
 * Sparse Vector
 * Stores a vector in sparse form
 */

template<typename IdxType, typename CoefType>
class SparseVectorType
{
public:
	using size_type = std::size_t;
	using view_type = SparseVectorViewType<IdxType,CoefType>;
	using iterator = typename view_type::iterator;
	// Constructors
	SparseVectorType() = default;
	// Initialize a sparse vector from a pair of C arrays
	SparseVectorType(std::span<const IdxType> idx, std::span<const CoefType> coef) : indices(idx.begin(), idx.end()), coefs(coef.begin(), coef.end())
	{
		assert( indices.size() == coefs.size() );
	}
	// Initialize a sparse vector from a dense vector given as a span
	SparseVectorType(std::span<const CoefType> dense, double eps = defaultEPS)
	{
		gather(dense, eps);
	}
	// Initialize a sparse vector from a list of (i,v) pairs
	SparseVectorType(std::span<const std::pair<IdxType,CoefType>> sp)
	{
		for (const auto& [i,v]: sp)  push(i,v);
	}
	// Initialize a sparse vector from a list of (i,v) pairs (as initializer lists)
	SparseVectorType(std::initializer_list<const std::pair<IdxType,CoefType>> lst)
	{
		for (const auto& [i,v]: lst)  push(i,v);
	}
	// Operators
	bool operator==(const SparseVectorType& other) const
	{
		return (indices == other.indices) && (coefs == other.coefs);
	}
	bool operator!=(const SparseVectorType& other) const
	{
		return !(*this == other);
	}
	operator view_type() const { return view_type(idx(), coef(), size()); }
	// size
	inline size_type size() const { return indices.size(); }
	inline size_type capacity() const { return indices.capacity(); }
	void resize(size_type newSize) { indices.resize(newSize); coefs.resize(newSize); }
	inline bool empty() const { return indices.empty(); }
	void reserve(size_type n) { indices.reserve(n); coefs.reserve(n); }
	inline void clear() { indices.clear(); coefs.clear(); }
	// push/pop
	inline void push(IdxType i, CoefType v)
	{
		indices.emplace_back(i);
		coefs.emplace_back(v);
	}
	inline void pop()
	{
		indices.pop_back();
		coefs.pop_back();
	}
	// copy from pair of C arrays
	void copy(std::span<const IdxType> idx, std::span<const CoefType> coef)
	{
		assert( idx.size() == coef.size() );
		resize(idx.size());
		std::copy(idx.begin(), idx.end(), indices.begin());
		std::copy(coef.begin(), coef.end(), coefs.begin());
	}
	// get data
	//@{
	view_type view() const { return view_type(idx(), coef(), indices.size()); }
	IdxType* idx() { return indices.empty() ? nullptr : &(indices[0]); }
	const IdxType* idx() const { return indices.empty() ? nullptr : &(indices[0]); }
	CoefType* coef() { return coefs.empty() ? nullptr : &(coefs[0]); }
	const CoefType* coef() const { return coefs.empty() ? nullptr : &(coefs[0]); }
	//@}
	// iterators
	// @{
	iterator begin() const { return view().begin(); }
	iterator end() const { return view().end(); }
	// @}
	/** conversions to/from dense vectors */
	/** Shrink a dense vector into a sparse one */
	void gather(std::span<const CoefType> in, double eps = defaultEPS)
	{
		clear();
		for (IdxType i = 0; i < (IdxType)in.size(); ++i) if (isNotNull(in[i], eps)) push(i, in[i]);
	}
	/** Expand a sparse vector into a dense one */
	void scatter(std::span<CoefType> out, bool reset = false)
	{
		if (reset) std::fill(out.begin(), out.end(), CoefType(0));
		size_type cnt = size();
		for (size_type i = 0; i < cnt; ++i) out[indices[i]] = coefs[i];
	}
	/** Zero out the coefficients of a dense vector corresponding to the support of a sparse vector */
	void unscatter(std::span<CoefType> out)
	{
		size_type cnt = size();
		for (size_type i = 0; i < cnt; ++i) out[indices[i]] = CoefType(0);
	}
	/** Scale sparse vector by a given multiplier */
	void scale(double lambda = 1.0)
	{
		size_type cnt = size();
		for (size_type i = 0; i < cnt; ++i) coefs[i] *= lambda;
	}
	/** Invert signs (e.g., multiply by -1) */
	void negate()
	{
		size_type cnt = size();
		for (size_type i = 0; i < cnt; ++i) coefs[i] = -coefs[i];
	}
private:
	std::vector<IdxType> indices;
	std::vector<CoefType> coefs;
};

using SparseVector = SparseVectorType<int,double>;


/**
 * Linear Constraint
 */

class Constraint
{
public:
	/** @name Data */
	//@{
	std::string name; //< constraint name
	SparseVector row; /**< coefficient row \f$ a_i \f$ */
	double rhs; /**< right hand side \f$ b_i \f$ */
	double range = 0.0; /**< range value for ranged row: linear expression in [rhs-range,rhs] */
	char sense; /**< constraint sense */
	//@}
	Constraint* clone() const { return new Constraint(*this); }
	/**
	 * Check if this constraint is satisfied by assignment x, with tolerance eps
	 */
	bool satisfiedBy(const double* x, double eps = defaultEPS) const;
	/**
	 * Compute the violation of constraint by assignment x
	 * A positive value means a constraint violation,
	 * while a negative one means constraint is slack
	 */
	double violation(const double* x) const;
	/**
	 * Check if a constraint is slack w.r.t. assignment x, with tolerance eps
	 * @return true if constraint is slack, false otherwise
	 */
	bool isSlack(const double* x, double eps = defaultEPS) const;
};

inline bool Constraint::satisfiedBy(const double* x, double eps) const
{
	return (!isPositive(violation(x), eps));
}

inline bool Constraint::isSlack(const double* x, double eps) const
{
	return isNegative(violation(x), eps);
}

typedef std::shared_ptr<Constraint> ConstraintPtr;

/**
 * Performe the operation: v <- v + lambda w
 * where both v and w are dense vectors and lambda is a scalar
 */

void accumulate(double* v, const double* w, int n, double lambda = 1.0);

/**
 * Performe the operation: v <- v + lambda w
 * where v is a dense vector, w is sparse and lambda is a scalar
 */

void accumulate(double* v, const int* wIdx, const double* wCoef, int n, double lambda = 1.0);

inline void accumulate(double* v, const SparseVector::view_type& w, double lambda = 1.0)
{
	accumulate(v, w.idx(), w.coef(), w.size(), lambda);
}

/**
 * Scale a dense vector v multiplying it by lambda: v <- lambda v
 */

void scale(double* v, int n, double lambda = 1.0);

/**
 * Dot Product between dense vectors (manual loop unrolling)
 */

double dotProduct(const double* a1, const double* a2, int n);

/**
 * Dot Product between a sparse and a dense vector
 */

double dotProduct(const int* idx1, const double* a1, int n, const double* a2);

inline double dotProduct(const SparseVector::view_type& a1, const double* a2)
{
	return dotProduct(a1.idx(), a1.coef(), a1.size(), a2);
}

/**
 * Checks if two dense vectors have disjoint support
 */

bool disjoint(const double* a1, const double* a2, int n);

/**
 * Euclidian norm on a dense vector
 */

double euclidianNorm(const double* v, int n);

/**
 * Euclidian distance between two dense vectors
 */

double euclidianDistance(const double* a1, const double* a2, int n);

/**
 * Lexicographically compare two array of doubles
 * @param s1 first array
 * @param s2 second array
 * @param n arrays size
 * @return comparison result {-1, 0, 1} if {<, =, >} respectively
 */

int lexComp(const double* s1, const double* s2, int n);


/**
 * Sparse Matrix
 *
 * A sparse matrix is basically a list of k sparse vectors,
 * each of which has index set [0,U). The index range is used
 * for transpose operations.
 *
 * The data structure is quite standard: we store all sparse
 * vectors one after the other in two arrays (ind and val)
 * and keep track of where each vector starts (beg) and
 * how many items it has (cnt). We could do with beg only
 * assuming there are no holes, but having an explicit
 * cnt is more flexible (and allows for holes, for example).
 */
class SparseMatrix
{
public:
	int k = 0; //< number of sparse vectors
	int U = 0; //< index range for vectors
	int nnz = 0; //< total number of nonzeros
	std::vector<int> beg; //< beginning of each vector in data arrays
	std::vector<int> cnt; //< length of each vector
	std::vector<int> ind; //< indices
	std::vector<double> val; //< values
	// constructors
	SparseMatrix() = default;
	SparseMatrix(std::initializer_list<std::initializer_list<double>> lst)
	{
		for (const auto& row: lst)
		{
			U = std::max(U, (int)(row.size()));
			add(SparseVector(as_span(row)));
		}
	}
	using view_type = SparseVectorViewType<int,double>;
	view_type operator[](int i) const
	{
		return view_type(ind.data() + beg[i], val.data() + beg[i], cnt[i]);
	}
	/* Append a new sparse vector to the matrix */
	void add(SparseMatrix::view_type row);
	/* Construct the transposed matrix */
	SparseMatrix transpose() const;
	// Operators
	bool operator==(const SparseMatrix& other) const
	{
		if ((k != other.k) || (U != other.U) || (nnz != other.nnz))  return false;
		for (int i = 0; i < U; i++)
		{
			if ((*this)[i] != other[i])  return false;
		}
		return true;
	}
	bool operator!=(const SparseMatrix& other) const
	{
		return !(*this == other);
	}
};


/**
 * Unary predicate that incrementally compute the variance of a list of numbers
 * The results can be obtained with result()
 * The flag fromSample decides if the variance was from the entire population
 * or a subset of observation, i.e. if fromSample is true the sum of squares is
 * divided by n - 1, otherwise by n
 */

class IncrementalVariance
{
public:
	IncrementalVariance() : cnt(0), mean(0.0), sumsq(0.0) {}
	void operator() (double x)
	{
		cnt++;
		double delta = x - mean;
		mean += delta / cnt;
		sumsq += delta * (x - mean);
	}
	int count() const { return cnt; }
	double result(bool fromSample = false) const { return (cnt > 1) ? (sumsq / (cnt - int(fromSample))) : 0.0; }
protected:
	int cnt;
	double mean;
	double sumsq;
};

/**
 * Return the arithmetic mean of a list of numbers
 */

template<typename ForwardIterator, typename T = double>
T mean(ForwardIterator first, ForwardIterator last)
{
	T sum = 0.0;
	int cnt = 0;
	std::for_each(first, last, [&](T value) { sum += value; ++cnt; });
	return (cnt > 0) ? (sum / cnt) : sum;
}

/**
 * Return the geometric mean of a list of positive numbers
 */

template<typename ForwardIterator, typename T = double>
T geomMean(ForwardIterator first, ForwardIterator last)
{
	T sum = 0.0;
	int cnt = 0;
	std::for_each(first, last, [&](T value) { sum += std::log(value); ++cnt; });
	return (cnt > 0) ? std::exp(sum / cnt) : sum;
}

/**
 * Return the variance of a list of numbers
 * The flag fromSample decides if the variance was from the entire population
 * or a subset of observation, i.e. if fromSample is true the sum of squares is
 * divided by n - 1, otherwise by n
 */

template<typename ForwardIterator, typename T = double>
T variance(ForwardIterator first, ForwardIterator last, bool fromSample = false)
{
	T m = mean(first, last);
	int cnt = 0;
	T sum = 0.0;
	std::for_each(first, last, [&](T value) { sum += std::pow(value - m, 2.0); ++cnt; });
	return (cnt > 1) ? (sum / (cnt - int(fromSample))) : sum;
}

/**
 * Return the standard deviation of a list of numbers
 * The flag fromSample decides if the variance was from the entire population
 * or a subset of observation, i.e. if fromSample is true the sum of squares is
 * divided by n - 1, otherwise by n
 */

template<typename ForwardIterator, typename T = double>
T stdev(ForwardIterator first, ForwardIterator last, bool fromSample = false)
{
	return std::sqrt(variance(first, last, fromSample));
}

/**
 * Compute the hash of a set of objects incrementally
 */
template <class T>
inline void hash_combine(std::size_t& seed, const T& v)
{
	std::hash<T> hasher;
	seed ^= hasher(v) + 0x9e3779b9 + (seed<<6) + (seed>>2);
}

/**
 * Compute the hash of the elements in the range [first,last)
 */
template <class ForwardIterator>
inline void hash_range(std::size_t& seed, ForwardIterator first, ForwardIterator last)
{
	while (first != last) hash_combine(seed, *first++);
}


} // namespace dominiqs

#endif /* MATHS_H */
