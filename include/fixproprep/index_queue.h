#ifndef INDEX_QUEUE_H
#define INDEX_QUEUE_H

#include <vector>
#include <cassert>


/** Circular queue to store up to k indices of type Integer in the range [0,n) */
template<typename Integer>
class IndexQueue
{
public:
	IndexQueue(Integer _k, Integer _n) : k{_k}, n{_n}, indices(k), ismember(n, false)
	{
		assert( k > 0 );
		assert( n > 0 );
	}
	bool has(Integer x) const
	{
		assert( 0 <= x && x < n );
		return ismember[x];
	}
	bool empty() const { return (count == 0); }
	bool full() const { return (count == k); }
	size_t size() const
	{
		return (size_t)(count);
	}
	void clear()
	{
		if (count) {
			for (int itr = 0; itr < count; itr++)  ismember[indices[(first + itr) % k]] = false;
			count = 0;
			first = 0;
			last = 0;
		}
		else {
			assert( first == last );
		}
	}
	Integer operator[](Integer itr) const
	{
		assert( (0 <= itr) && (itr < count) );
		return indices[(first + itr) % k];
	}
	void push(Integer x) {
		assert( 0 <= x && x < n );
		// do nothing if the index is already in the queue
		if (has(x))  return;
		// if full pop the oldest element (this is basically overwriting old data)
		if (full())  pop();
		// add x
		assert( !full() );
		ismember[x] = true;
		indices[last] = x;
		increment(last);
		count++;
		assert( (0 < count) && (count <= k) );
		assert( (0 <= first) && (first < k) );
		assert( (0 <= last) && (last < k) );
	}
	Integer pop() {
		assert( !empty() );
		Integer ret = indices[first];
		ismember[ret] = false;
		increment(first);
		count--;
		assert( (0 <= count) && (count < k) ); 
		assert( (0 <= first) && (first < k) );
		assert( (0 <= last) && (last < k) );
		return ret;
	}
private:
	Integer k;
	Integer first = 0;
	Integer last = 0;
	Integer count = 0;
	Integer n;
	std::vector<Integer> indices;
	std::vector<bool> ismember;
	void increment(Integer& itr) const {
		itr++;
		if (itr == k)  itr = 0;
	} 
};

#endif /* INDEX_QUEUE_H */
