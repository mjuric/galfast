#ifndef __analysis_h
#define __analysis_h

#include <astro/io/gzstream/fstream.h>
#include "textstream.h"
#include <set>
#include <valarray>

////////////////////////////
//
//   Data analysis macros
//

#define output_or_die(var, filename) \
	std::ofstream var((std::string(filename)).c_str()); \
	if(!var.good()) { die("Could not open output file [" + std::string(filename) + "]"); }

#define input_or_die(var, filename) \
	std::ifstream var((std::string(filename)).c_str()); \
	if(!var.good()) { die("Could not open input file [" + std::string(filename) + "]"); }

#define text_input_or_die(var, filename) \
	input_or_die(var##_stream, filename); \
	itextstream var(var##_stream);

#define text_output_or_die(var, filename) \
	output_or_die(var##_stream, filename); \
	otextstream var(var##_stream);

#define binary_input_or_die(var, filename) \
	input_or_die(var##_stream, filename); \
	ibinarystream var(var##_stream, filename);

#define binary_output_or_die(var, filename) \
	output_or_die(var##_stream, filename); \
	obinarystream var(var##_stream, filename);

//// gzip streams

#define gz_output_or_die(var, filename) \
	peyton::io::gzstream::ofstream var((std::string(filename)+".gz").c_str()); \
	if(!var.good()) { die("Could not open output file [" + std::string(filename) + ".gz]"); }

#define gz_input_or_die(var, filename) \
	peyton::io::gzstream::ifstream var((std::string(filename)+".gz").c_str()); \
	if(!var.good()) { die("Could not open input file [" + std::string(filename) + ".gz]"); }

#define gz_text_input_or_die(var, filename) \
	gz_input_or_die(var##_stream, filename); \
	itextstream var(var##_stream);

#define gz_text_output_or_die(var, filename) \
	gz_output_or_die(var##_stream, filename); \
	otextstream var(var##_stream);

#define gz_binary_input_or_die(var, filename) \
	gz_input_or_die(var##_stream, filename); \
	ibinarystream var(var##_stream);

#define gz_binary_output_or_die(var, filename) \
	gz_output_or_die(var##_stream, filename); \
	obinarystream var(var##_stream);

#define TN(i) typename T##i
#define TP(i) T##i &v##i, int p##i
#define TB(i) in.bind(v##i, p##i)

template<TN(1), TN(2), TN(3), TN(4), TN(5), TN(6), TN(7), TN(8), TN(9), TN(10)> void bind(itextstream &in, TP(1), TP(2), TP(3), TP(4), TP(5), TP(6), TP(7), TP(8), TP(9), TP(10)) { TB(1); TB(2); TB(3); TB(4); TB(5); TB(6); TB(7); TB(8); TB(9); TB(10); }
template<TN(1), TN(2), TN(3), TN(4), TN(5), TN(6), TN(7), TN(8), TN(9)> void bind(itextstream &in, TP(1), TP(2), TP(3), TP(4), TP(5), TP(6), TP(7), TP(8), TP(9)) { TB(1); TB(2); TB(3); TB(4); TB(5); TB(6); TB(7); TB(8); TB(9); }
template<TN(1), TN(2), TN(3), TN(4), TN(5), TN(6), TN(7), TN(8)> void bind(itextstream &in, TP(1), TP(2), TP(3), TP(4), TP(5), TP(6), TP(7), TP(8)) { TB(1); TB(2); TB(3); TB(4); TB(5); TB(6); TB(7); TB(8); }
template<TN(1), TN(2), TN(3), TN(4), TN(5), TN(6), TN(7)> void bind(itextstream &in, TP(1), TP(2), TP(3), TP(4), TP(5), TP(6), TP(7)) { TB(1); TB(2); TB(3); TB(4); TB(5); TB(6); TB(7); }
template<TN(1), TN(2), TN(3), TN(4), TN(5), TN(6)> void bind(itextstream &in, TP(1), TP(2), TP(3), TP(4), TP(5), TP(6)) { TB(1); TB(2); TB(3); TB(4); TB(5); TB(6); }
template<TN(1), TN(2), TN(3), TN(4), TN(5)> void bind(itextstream &in, TP(1), TP(2), TP(3), TP(4), TP(5)) { TB(1); TB(2); TB(3); TB(4); TB(5); }
template<TN(1), TN(2), TN(3), TN(4)> void bind(itextstream &in, TP(1), TP(2), TP(3), TP(4)) { TB(1); TB(2); TB(3); TB(4); }
template<TN(1), TN(2), TN(3)> void bind(itextstream &in, TP(1), TP(2), TP(3)) { TB(1); TB(2); TB(3); }
template<TN(1), TN(2)> void bind(itextstream &in, TP(1), TP(2)) { TB(1); TB(2); }
template<TN(1)> void bind(itextstream &in, TP(1)) { TB(1); }


//
// The template machinery which follows allows load() to gobble up almost anything.
// load() accesses the variable passed through an adaptor class, which provides a
// push_back function. If a class already has a push_back function (eg., a vector).
// the default adaptor just call it. If the class has no push_back function, the
// user must specialize the adapt_c template for the class and write the appropriate
// push_back function.
// For STL classes with no push_back (eg. set), the specialized adaptors have been
// written (you can look to them as an example).
//

//
// Default load adaptor - just call the inner class push_back
//
template<class A>
class adapt_c {
protected:
	A &a;
	typedef typename A::value_type value_type;
public:
	adapt_c(A &a_) : a(a_) {}

	void push_back(const typename A::value_type &s) { a.push_back(s); }
	void clear() { a.clear(); }
};

//
// Set template adaptor - call set<>::insert
//
template<class A>
class adapt_c< std::set<A> > {
protected:
	std::set<A> &a;
	typedef A value_type;
public:
	adapt_c(std::set<A> &a_) : a(a_) {}

	void push_back(const value_type &s) { a.insert(s); }
	void clear() { a.clear(); }
};

//
// valarray class adaptor - it assumes that _the array has enough space
// to contain the inserted values_
//
template<class A>
class adapt_c< std::valarray<A> > {
protected:
	std::valarray<A> &a;
	typedef A value_type;
	int at;
public:
	adapt_c(std::valarray<A> &a_) : a(a_), at(0) {}

	void push_back(const value_type &s) { a[at]=s; at++; }
	void clear() { at = 0; }
};

#undef TB
#define TB(i) \
	adapt_c<T##i> vr##i(v##i); \
	typename T##i::value_type s##i; \
	if(p##i >= 0) { in.bind(s##i, p##i); }
#define TS(i) vr##i.push_back(s##i);
#define TU(i) in.unbind(s##i);

class ticker {
public:
	int tickk;
	int step;

	typedef int value_type;
public:
	ticker(int step_) : tickk(-1) { open("", step_); }
	ticker(const std::string &title, int step_) : tickk(-1) { open(title, step_); }
	~ticker() { close(); }

	void open(const std::string &title, int step_)
	{
		close();
		
		step = step_;
		tickk = 0;

		if(title.size())
		{
			std::cerr << title << "...\n";
		}
	}
	
	void close()
	{
		if(tickk > 0)
		{
			std::cerr << " [" << tickk << "].\n";
		}
		tickk = 0;
	}

	int count() const { return tickk; }

	void tick()
	{
		tickk++;
		if(tickk % step == 0) { std::cerr << "#"; }
		if(tickk % (step*50) == 0) { std::cerr << " [" << tickk << "]\n"; }
	}

	void push_back(const int &s) { tick(); }
};

class spinner
{
public:
	const static int NCHARS = 4;
	static char chars[NCHARS+1];

public:
	int at;
	bool first;
		
	spinner() : at(0), first(true) {}
	void tick();
	void start();
	void stop();
};

//
// This is the type of function that the heavily macro-ed functions below generate:
//
#if 0
template<typename T1>
void load(itextstream &in, T1 &vr1, int p1)
{
	adapt_c<T1> v1(vr1);
	typename T1::value_type s1;
	if(p1 >= 0) { in.bind(s1, p1); }

	while(in.next()) {
		v1.push_back(s1);
	}
	in.unbind(s1);
}
#endif

template<TN(1), TN(2), TN(3), TN(4), TN(5), TN(6), TN(7), TN(8), TN(9), TN(10)>
int load(itextstream &in, TP(1), TP(2), TP(3), TP(4), TP(5), TP(6), TP(7), TP(8), TP(9), TP(10))
{
	TB(1); TB(2); TB(3); TB(4); TB(5); TB(6); TB(7); TB(8); TB(9); TB(10);
	int cnt = 0;
	while(in.next()) {
		TS(1); TS(2); TS(3); TS(4); TS(5); TS(6); TS(7); TS(8); TS(9); TS(10);
		cnt++;
	}
	TU(1); TU(2); TU(3); TU(4); TU(5); TU(6); TU(7); TU(8); TU(9); TU(10);
	return cnt;
}

template<TN(1), TN(2), TN(3), TN(4), TN(5), TN(6), TN(7), TN(8), TN(9)>
int load(itextstream &in, TP(1), TP(2), TP(3), TP(4), TP(5), TP(6), TP(7), TP(8), TP(9))
{
	TB(1); TB(2); TB(3); TB(4); TB(5); TB(6); TB(7); TB(8); TB(9);
	int cnt = 0;
	while(in.next()) {
		TS(1); TS(2); TS(3); TS(4); TS(5); TS(6); TS(7); TS(8); TS(9);
		cnt++;
	}
	TU(1); TU(2); TU(3); TU(4); TU(5); TU(6); TU(7); TU(8); TU(9);
	return cnt;
}

template<TN(1), TN(2), TN(3), TN(4), TN(5), TN(6), TN(7), TN(8)>
int load(itextstream &in, TP(1), TP(2), TP(3), TP(4), TP(5), TP(6), TP(7), TP(8))
{
	TB(1); TB(2); TB(3); TB(4); TB(5); TB(6); TB(7); TB(8);
	int cnt = 0;
	while(in.next()) {
		TS(1); TS(2); TS(3); TS(4); TS(5); TS(6); TS(7); TS(8);
		cnt++;
	}
	TU(1); TU(2); TU(3); TU(4); TU(5); TU(6); TU(7); TU(8);
	return cnt;
}

template<TN(1), TN(2), TN(3), TN(4), TN(5), TN(6), TN(7)>
int load(itextstream &in, TP(1), TP(2), TP(3), TP(4), TP(5), TP(6), TP(7))
{
	TB(1); TB(2); TB(3); TB(4); TB(5); TB(6); TB(7);
	int cnt = 0;
	while(in.next()) {
		TS(1); TS(2); TS(3); TS(4); TS(5); TS(6); TS(7);
		cnt++;
	}
	TU(1); TU(2); TU(3); TU(4); TU(5); TU(6); TU(7);
	return cnt;
}

template<TN(1), TN(2), TN(3), TN(4), TN(5), TN(6)>
int load(itextstream &in, TP(1), TP(2), TP(3), TP(4), TP(5), TP(6))
{
	TB(1); TB(2); TB(3); TB(4); TB(5); TB(6);
	int cnt = 0;
	while(in.next()) {
		TS(1); TS(2); TS(3); TS(4); TS(5); TS(6);
		cnt++;
	}
	TU(1); TU(2); TU(3); TU(4); TU(5); TU(6);
	return cnt;
}

template<TN(1), TN(2), TN(3), TN(4), TN(5)>
int load(itextstream &in, TP(1), TP(2), TP(3), TP(4), TP(5))
{
	TB(1); TB(2); TB(3); TB(4); TB(5);
	int cnt = 0;
	while(in.next()) {
		TS(1); TS(2); TS(3); TS(4); TS(5);
		cnt++;
	}
	TU(1); TU(2); TU(3); TU(4); TU(5);
	return cnt;
}

template<TN(1), TN(2), TN(3), TN(4)>
int load(itextstream &in, TP(1), TP(2), TP(3), TP(4))
{
	TB(1); TB(2); TB(3); TB(4);
	int cnt = 0;
	while(in.next()) {
		TS(1); TS(2); TS(3); TS(4);
		cnt++;
	}
	TU(1); TU(2); TU(3); TU(4);
	return cnt;
}

template<TN(1), TN(2), TN(3)>
int load(itextstream &in, TP(1), TP(2), TP(3))
{
	TB(1); TB(2); TB(3);
	int cnt = 0;
	while(in.next()) {
		TS(1); TS(2); TS(3);
		cnt++;
	}
	TU(1); TU(2); TU(3);
	return cnt;
}

template<TN(1), TN(2)>
int load(itextstream &in, TP(1), TP(2))
{
	TB(1); TB(2);
	int cnt = 0;
	while(in.next()) {
		TS(1); TS(2);
		cnt++;
	}
	TU(1); TU(2);
	return cnt;
}

template<TN(1)>
int load(itextstream &in, TP(1))
{
	TB(1);
	int cnt = 0;
	while(in.next()) {
		TS(1);
		cnt++;
	}
	TU(1);
	return cnt;
}

/*
template<TN(1), TN(2), TN(3), TN(4), TN(5), TN(6), TN(7), TN(8), TN(9)> void bind(itextstream &in, TP(1), TP(2), TP(3), TP(4), TP(5), TP(6), TP(7), TP(8), TP(9)) { TB(1); TB(2); TB(3); TB(4); TB(5); TB(6); TB(7); TB(8); TB(9); }
template<TN(1), TN(2), TN(3), TN(4), TN(5), TN(6), TN(7), TN(8)> void bind(itextstream &in, TP(1), TP(2), TP(3), TP(4), TP(5), TP(6), TP(7), TP(8)) { TB(1); TB(2); TB(3); TB(4); TB(5); TB(6); TB(7); TB(8); }
template<TN(1), TN(2), TN(3), TN(4), TN(5), TN(6), TN(7)> void bind(itextstream &in, TP(1), TP(2), TP(3), TP(4), TP(5), TP(6), TP(7)) { TB(1); TB(2); TB(3); TB(4); TB(5); TB(6); TB(7); }
template<TN(1), TN(2), TN(3), TN(4), TN(5), TN(6)> void bind(itextstream &in, TP(1), TP(2), TP(3), TP(4), TP(5), TP(6)) { TB(1); TB(2); TB(3); TB(4); TB(5); TB(6); }
template<TN(1), TN(2), TN(3), TN(4), TN(5)> void bind(itextstream &in, TP(1), TP(2), TP(3), TP(4), TP(5)) { TB(1); TB(2); TB(3); TB(4); TB(5); }
template<TN(1), TN(2), TN(3), TN(4)> void bind(itextstream &in, TP(1), TP(2), TP(3), TP(4)) { TB(1); TB(2); TB(3); TB(4); }
template<TN(1), TN(2), TN(3)> void bind(itextstream &in, TP(1), TP(2), TP(3)) { TB(1); TB(2); TB(3); }
template<TN(1), TN(2)> void bind(itextstream &in, TP(1), TP(2)) { TB(1); TB(2); }
template<TN(1)> void bind(itextstream &in, TP(1)) { TB(1); }
*/

#undef TN
#undef TP
#undef TB

inline bool die(const std::string &s)
{
	std::cout << s << "\n";
	exit(-1);
	return false;
}

// initialize a scalar value to zero upon creation
// useful in arrays and hashtables
template <typename T>
struct zero_init
{
	T val;
	zero_init() { val = T(0); }
	operator T& () { return val; }
	operator T () const { return val; }
	T& operator =(const T& a) { val = a.val; return *this; }
/*	T& operator ++(int) { return val++; }
	T& operator ++() { return ++val; }*/
};

#endif
