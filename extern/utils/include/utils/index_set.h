/**
 * @file index_set.h
 * @brief Indexed set data structure
 *
 * @author Domenico Salvagnin dominiqs@gmail.com
 */

#ifndef INDEX_SET_H
#define INDEX_SET_H

#include <vector>
#include <span>


namespace dominiqs {

/** @brief Unordered container for integer indices in the range [0,size)
 *
 * - Insertion, removal and membership are O(1)
 * - Indices are stored contiguosly in memory (fast traversal)
 * - Memory consumption is proportial to the size of the range,
 *   not to the actual number of elements in the container.
 */
template<typename Integer>
class IndexSet
{
public:
	IndexSet(Integer n) : pos(n, -1) {}
	Integer maxSize() const { return (Integer)pos.size(); }
	Integer size() const { return (Integer)items.size(); }
	bool empty() const { return (size() == 0); }
	bool has(Integer x) const { return (pos[x] != -1); }
	Integer operator[](Integer index) const { return items[index]; }
	void add(Integer x)
	{
		if (!has(x)) {
			items.push_back(x);
			pos[x] = (Integer)items.size()-1;
		}
	}
	void remove(Integer x)
	{
		if (has(x)) {
			// swap with item in last position
			Integer last = items.back();
			items.pop_back();
			pos[last] = pos[x];
			items[pos[last]] = last;
			pos[x] = -1;
		}
	}
	Integer pop()
	{
		Integer last = items.back();
		items.pop_back();
		pos[last] = -1;
		return last;
	}
	void clear() {
		while (!empty())  pop();
	}
	std::span<const Integer> data() const { return items; }
private:
	std::vector<Integer> items;
	std::vector<Integer> pos;
};


} // namespace

#endif /* INDEX_SET_H */
