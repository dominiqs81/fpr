/**
 * @file vectorlist.h
 * @brief Vector List data structure
 *
 * @author Domenico Salvagnin dominiqs@gmail.com
 */

#ifndef VECTORLIST_H
#define VECTORLIST_H

#include <span>
#include <vector>
#include <cassert>

/**
 * Vector List
 *
 * Stores a list of sparse vectors of integer values in the range [0, maxVal)
 */
class VectorList
{
public:
	VectorList() : beg{0} {}
	inline int nVectors() const { return numVec; }
	inline int maxValue() const { return maxVal; }
	inline size_t nNonzeros() const { return beg[numVec]; }
	/* Get a view on the i-th vector */
	using view_type = std::span<const int>;
	inline view_type operator[](int i) const
	{
		assert((i >= 0) && (i < numVec));
		return view_type(&data[beg[i]], beg[i+1]-beg[i]);
	}
	inline void setMaxValue(int mv) { maxVal = mv; }
	/* Add a vector (passed as a span) to the list */
	inline void add(std::span<const int> vec)
	{
		add(vec.begin(), vec.end());
	}
	/* Add a vector (passed as a initializer_list) to the list */
	template<typename T>
	inline void add(std::initializer_list<T> vec)
	{
		add(vec.begin(), vec.end());
	}
	/* Construct a transposed vector list */
	VectorList transpose() const;
	/* Comparison operators */
	bool operator==(const VectorList& other) const
	{
		if (nVectors() != other.nVectors())  return false;
		if (maxValue() != other.maxValue())  return false;
		if (nNonzeros() != other.nNonzeros())  return false;
		if (!std::equal(beg.begin(), beg.end(), other.beg.begin()))  return false;
		if (!std::equal(data.begin(), data.end(), other.data.begin()))  return false;
		return true;
	}
	bool operator!=(const VectorList& other) const
	{
		return !(*this == other);
	}
private:
	int numVec = 0; //< number of vectors
	int maxVal = 1; //< upper bound on integer values
	std::vector<size_t> beg; //< size k+1
	std::vector<int> data; //< size #nnz
	// helper function to handle addition
	template<typename ForwardIterator>
	inline void add(ForwardIterator first, ForwardIterator last)
	{
		data.insert(data.end(), first, last);
		beg.push_back(data.size());
		numVec++;
	}
};

#endif /* VECTORLIST_H */
