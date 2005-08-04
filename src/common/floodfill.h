#ifndef __floodfill_h
#define __floodfill_h

#include <astro/math/vector.h>

#include <ext/hash_set>
#include <queue>

namespace __gnu_cxx
{
	template<> struct hash<I2> { size_t operator()(const I2 &v) const { return (v.x << 16) + v.y; } };
};

template <typename T>
void neighbors(std::queue<T> &q, const T &v)
{
	ASSERT(0); // this function must be overloaded for every class that you flood with
}

template <>
void neighbors(std::queue<I2> &q, const I2 &v)
{
	q.push(I2(v.x+1, v.y  ));
	q.push(I2(v.x-1, v.y  ));
	q.push(I2(v.x,   v.y+1));
	q.push(I2(v.x,   v.y-1));
}

template <typename T>
class FloodFill
{
public:
	__gnu_cxx::hash_set<T, __gnu_cxx::hash<T> > volume;
	std::queue<T> q;

public:
	FloodFill() { }

	virtual void progress() {};

	virtual bool test(const T &v) = 0;
	virtual void action(const T &v) = 0;

	void reset() { volume.clear(); q = std::queue<T>(); }

	inline void FloodFill::flood(const T &v_initial)
	{
		// flood search volume
		q.push(v_initial);

		// sanity check
		ASSERT(test(v_initial));

		while(q.size()) {
			const T v = q.front(); q.pop();

			if(volume.count(v)) { continue; }
			progress();

			if(!test(v)) { continue; }
			action(v);

			volume.insert(v);

			neighbors(q, v);
		}
	}
};

#endif
