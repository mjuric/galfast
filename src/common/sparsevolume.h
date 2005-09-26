#include "binarystream.h"

#include <ext/hash_map>
#include <map>
#include <astro/math/vector.h>
#include <vector>
#include "xcat.h"

//
// The result of the binning will be stored in a map, since the volume is sparsly sampled
//

struct sparse_volume
{
	typedef peyton::math::S3 K;
	typedef float V;

	struct hash_S3 { size_t operator()(const K &v) const { return (v.x << 22) + (v.y << 11) + v.z; } };
	struct less_S3 { bool operator()(const K &v1, const K &v2) const { return abs(v1) < abs(v2); } };
	typedef __gnu_cxx::hash_map<K, V, hash_S3> volume_map;
//	typedef std::map<K, V, less_S3> volume_map;

	volume_map volume;				///< map between quantized 3D space and volume covered at given pixel
	float dx;						///< the physical linear dimension of pixels stored in vm [pc]

	double r0, r1;					///< dimensions of the mapped volume, in pixel units
};

inline obinarystream &operator <<(obinarystream &out, const sparse_volume::volume_map &vm)
{
	out << vm.size();
	FOREACH(vm) { out << (*i).first << (*i).second; }
	return out;
}

inline void read_sv_element(ibinarystream &in, sparse_volume::K &k, sparse_volume::V &v)
{
	in >> k >> v;
}

inline ibinarystream &operator >>(ibinarystream &in, sparse_volume::volume_map &vm)
{
	size_t size; sparse_volume::K k; sparse_volume::V v;

	vm.clear();

	in >> size;
	FOR(0, size) { read_sv_element(in, k, v); vm[k] = v; }
	return in;
}

inline obinarystream &operator <<(obinarystream &out, const sparse_volume &vm)
{
	out << vm.dx << vm.r0 << vm.r1 << vm.volume;
	return out;
}

inline ibinarystream &operator >>(ibinarystream &in, sparse_volume &vm)
{
	in >> vm.dx >> vm.r0 >> vm.r1 >> vm.volume;
	return in;
}

inline void set_or_add(sparse_volume &vm, const peyton::math::I3 &p, double v)
{
	sparse_volume::volume_map::iterator it;
	if((it = vm.volume.find(p)) == vm.volume.end())
	{
		vm.volume[p] = v;
	} else {
		(*it).second += v;
	}
}

class sparse_volume_streamer : public sparse_volume
{
protected:
	ibinarystream &in;
	header h;
	size_t size;
public:
	sparse_volume_streamer(ibinarystream &in_) : in(in_)
	{
		in >> h >> dx >> r0 >> r1 >> size;
	}
	bool next(sparse_volume::K &k, sparse_volume::V &v)
	{
		if(size == 0) return false;
		read_sv_element(in, k, v);
		size--;
		return true;
	}
};

///////

struct sparse_volume_info
{
	typedef peyton::math::S3 K;
	typedef std::vector<int> V;

	struct hash_S3 { size_t operator()(const K &v) const { return (v.x << 22) + (v.y << 11) + v.z; } };
	struct less_S3 { bool operator()(const K &v1, const K &v2) const { return abs(v1) < abs(v2); } };
	typedef __gnu_cxx::hash_map<K, V, hash_S3> volume_map;

	volume_map runs;				///< map between quantized 3D space and runs covering given pixel
};

inline obinarystream &operator <<(obinarystream &out, const sparse_volume_info::volume_map &vm)
{
	out << vm.size();
	FOREACH(vm) { out << (*i).first << (*i).second; }
	return out;
}

inline void read_sv_element(ibinarystream &in, sparse_volume_info::K &k, sparse_volume_info::V &v)
{
	in >> k >> v;
}

inline ibinarystream &operator >>(ibinarystream &in, sparse_volume_info::volume_map &vm)
{
	size_t size; sparse_volume_info::K k; sparse_volume_info::V v;

	vm.clear();

	in >> size;
	FOR(0, size) { read_sv_element(in, k, v); vm[k] = v; }
	return in;
}

inline obinarystream &operator <<(obinarystream &out, const sparse_volume_info &vm)
{
	out << vm.runs;
	return out;
}

inline ibinarystream &operator >>(ibinarystream &in, sparse_volume_info &vm)
{
	in >> vm.runs;
	return in;
}

