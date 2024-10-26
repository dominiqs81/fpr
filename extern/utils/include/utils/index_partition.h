/**
 * @file index_partition.h
 * @brief Indexed partition data structure
 *
 * @author Domenico Salvagnin dominiqs@gmail.com
 */

#ifndef INDEX_PARTITION_H
#define INDEX_PARTITION_H

#include <cassert>
#include <numeric>
#include <iterator>
#include <cstddef>


namespace dominiqs
{

/**
 * @brief Efficient data structure for keeping linked sets of indices
 *
 * The data structure stores numIndices indices from 0 to numIndices - 1 
 * in up to numSets sets of indices 0 to numSets - 1
 * A given index can be in at most one set at a given time.
 */

class IndexPartition
{
public:
	/* Creates an index partition with indices up to @param n and sets up to @param ns */
	IndexPartition(int n = 0, int ns = 0) : numIndices(n), numSets(ns), flink(numIndices + numSets), blink(numIndices + numSets)
	{
		std::iota(flink.begin(), flink.end(), 0);
		std::iota(blink.begin(), blink.end(), 0);
	}
	/* Resizes the partition to indices up to @param n and to sets up to @param ns (previous content is reset!) */
	inline void resize(int n, int ns)
	{
		numIndices = n;
		numSets = ns;
		flink.resize(numIndices + numSets);
		std::iota(flink.begin(), flink.end(), 0);
		blink.resize(numIndices + numSets);
		std::iota(blink.begin(), blink.end(), 0);
	}
	/* Reset the content of the set, without resizing */
	inline void clear()
	{
		std::iota(flink.begin(), flink.end(), 0);
		std::iota(blink.begin(), blink.end(), 0);
	}
	/* add index @param index to set @param set (does NOT check whether the item was already in other set!!!) */
	inline void add(int index, int set)
	{
		assert( index >= 0 && index < numIndices );
		assert( set >= 0 && set < numSets );
		assert( isInNoSet(index) );
		int k = flink[numIndices + set];
		flink[numIndices + set] = index;
		flink[index] = k;
		blink[k] = index;
		blink[index] = numIndices + set;
	}
	/* remove index @param index from its current set (safe even if the index was already in no set) */
	inline void remove(int index)
	{
		assert( index >= 0 && index < numIndices );
		flink[blink[index]] = flink[index];
		blink[flink[index]] = blink[index];
		flink[index] = index;
		blink[index] = index;
	}
	inline bool isEmpty(int set) const
	{
		return (flink[numIndices + set] == (numIndices + set));
	}
	inline int top(int set) const
	{
		assert( !isEmpty(set) );
		return flink[numIndices + set];
	}
	inline bool isInNoSet(int index) const
	{
		return (flink[index] == index);
	}
	inline int getNumIndices() const { return numIndices; }
	inline int getNumSets() const { return numSets; }
	/* Iterator */
	struct const_iterator
	{
	public:
		using iterator_category = std::bidirectional_iterator_tag;
		using difference_type   = std::ptrdiff_t;
		using value_type        = const int;
		using pointer           = value_type*;
		using reference         = value_type&;
		// Constructors
		const_iterator() {}
		explicit const_iterator(const IndexPartition* p, int pos) : part{p}, index{pos} {}
		// Dereference
		reference operator*() const { return index; }
		pointer operator->() const { return &index; }
		// Increment
		const_iterator& operator++() { index = part->flink[index]; return *this; }
		const_iterator operator++(int) { const_iterator tmp = *this; ++(*this); return tmp; }
		// Decrement
		const_iterator& operator--() { index = part->blink[index]; return *this; }
		const_iterator operator--(int) { const_iterator tmp = *this; --(*this); return tmp; }
		// Comparison
		friend bool operator==(const const_iterator& a, const const_iterator& b) { return (a.part == b.part) && (a.index == b.index); };
		friend bool operator!=(const const_iterator& a, const const_iterator& b) { return !(a == b); };
	private:
		const IndexPartition* part = nullptr;
		int index = 0;
	};
	const_iterator begin(int set) const { return const_iterator(this, flink[numIndices + set]); }
	const_iterator end(int set) const { return const_iterator(this, numIndices + set); }
protected:
	// data
	int numIndices;
	int numSets;
	std::vector<int> flink;
	std::vector<int> blink;
};


} // namespace

#endif /** INDEX_PARTITION_H */
