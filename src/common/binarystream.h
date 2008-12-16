#ifndef __binarystream_XXX_h
#define __binarystream_XXX_h

#include <iosfwd>
#include <string>
#include <cstring>
#include <ctime>
#include <astro/util.h>

//
// simple wrapper class to handle binary IO
//

class obinarystream {
public:
	void write(const char *v, size_t size);
public:
	std::ostream &f;
	std::string filename;
public:
	obinarystream(std::ostream &f_, const std::string &filename_ = "") : f(f_), filename(filename_) {}
};

template<typename T>
inline obinarystream &operator <<(obinarystream &out, const T &v)
{
	out.write(reinterpret_cast<const char *>(&v), sizeof(T));
	return out;
}

// String specializations

template<>
inline obinarystream &operator <<(obinarystream &out, const std::string &v)
{
	out << v.size();
	out.write(v.c_str(), v.size());
	return out;
}

inline obinarystream &operator <<(obinarystream &out, const char *&v)
{
	size_t len = strlen(v);
	out << len;
	out.write(v, len);
	return out;
}

///////////////////////////////////////////////////////////////////

class ibinarystream {
public:
	void read(char *v, size_t size);
public:
	std::istream &f;
	std::string filename;
public:
	ibinarystream(std::istream &f_, const std::string &filename_ = "") : f(f_), filename(filename_) {}
};

template<typename T>
inline ibinarystream &operator >>(ibinarystream &in, T &v)
{
	in.read(reinterpret_cast<char *>(&v), sizeof(T));
	return in;
}

// String specializations

template<>
inline ibinarystream &operator >>(ibinarystream &in, std::string &v)
{
	int len;
	in >> len;

	char buf[len+1];
	buf[len] = 0;
	in.read(buf, len);

	v = buf;

	return in;
}

#define BOSTREAM(T...) obinarystream &operator <<(obinarystream &out, T)
#define BISTREAM(T...) ibinarystream &operator >>(ibinarystream &in, T)

/////////////////////////////////////////////////////////////////

#include <map>
#include <astro/exceptions.h>

SIMPLE_EXCEPTION(EBinaryIO);

class header
{
protected:
	static const int magic;
	static const int current_version;
	int version;

	friend obinarystream &operator <<(obinarystream &out, const header &h);
	friend ibinarystream &operator >>(ibinarystream &in, header &h) throw (EBinaryIO);
	friend OSTREAM(const header &h);
public:
	std::string description;
	time_t datetime;

	typedef std::map<std::string, std::string> data_map;
	data_map data;
public:
	std::string &operator[](const std::string &s) { return data[s]; }

	header(std::string description);
	header();
};

obinarystream &operator <<(obinarystream &out, const header::data_map &data);
ibinarystream &operator >>(ibinarystream &in, header::data_map &data);

obinarystream &operator <<(obinarystream &out, const header &h);
ibinarystream &operator >>(ibinarystream &in, header &h) throw (EBinaryIO);

OSTREAM(const header &h);

#endif
