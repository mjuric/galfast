#ifndef __textstream_h
#define __textstream_h

#include <map>
#include <iostream>
#include <string>
#include <sstream>
#include <astro/util.h>

class textstream {
public:
	enum { UNDEF, INT, DOUBLE, FLOAT, STRING, SHORT };
public:
	class field {
	public:
		int type;
		int pos;
		union {
			int *v_int;
			short *v_short;
			double *v_double;
			float *v_float;
			std::string *v_string;
			void *v_void;
		} v;
	public:
		field() : type(textstream::UNDEF), pos(-1) { v.v_void = NULL; }

		field(int p, int &var) : type(textstream::INT), pos(p) { v.v_int = &var; }
		field(int p, short &var) : type(textstream::SHORT), pos(p) { v.v_short = &var; }
		field(int p, double &var) : type(textstream::DOUBLE), pos(p) { v.v_double =&var; }
		field(int p, float &var) : type(textstream::FLOAT), pos(p) { v.v_float = &var; }
		field(int p, std::string &var) : type(textstream::STRING), pos(p) { v.v_string =&var; }
	};
};




class itextstream : public textstream {
protected:
	std::istream &f;
	typedef std::map< int, field > field_map;
	field_map fields;
	int max_field;
	bool ret_comments;
public:
	int nread;
	std::string line;
	bool is_comment;
public:
	
public:
	itextstream(std::istream &f_) : f(f_), nread(0), max_field(0), ret_comments(false) {}

	void bind(std::string &s, int p) { fields[p] = field(p, s); max_field = std::max(max_field, p); }
	void bind(double &s, int p) { fields[p] = field(p, s); max_field = std::max(max_field, p); }
	void bind(float &s, int p) { fields[p] = field(p, s); max_field = std::max(max_field, p); }
	void bind(int &s, int p) { fields[p] = field(p, s); max_field = std::max(max_field, p); }
	void bind(short &s, int p) { fields[p] = field(p, s); max_field = std::max(max_field, p); }

	template<typename T>
	bool unbind(T &s)
	{
		void *ptr = (void *)&s;
		FOREACH(fields) {
			if((*i).second.v.v_void != ptr) continue;

			fields.erase(i);
			
			max_field = 0;
			FOREACH(fields) { max_field = std::max(max_field, (*i).first); }

			return true;
		}
		return false;
	}

	bool returncomments(bool rc) { ret_comments = rc; }
	
	itextstream &skip(int n = 1)
	{
		while(n > 0)
		{
			getline(f, line);
			n--;

			if(f.eof()) return *this;
			if(ret_comments && line[0] == '#') { n++; }
		}
	}

	bool isallspace(const std::string &s)
	{
		for(int i = 0; i != s.size(); i++)
		{
			if(!isspace(s[i])) { return false; }
		}
		return true;
	}

	itextstream &next()
	{
		nread = -1;
		do {
			getline(f, line);
//			std::cerr << "LINE: [" << line << "] " << (bool)*this << "\n";
//		} while(!ret_comments && (line[0] == '#' && !f.eof()));
		} while(!f.eof() && ((line[0] == '#' && !ret_comments) || isallspace(line)));

		if(f.eof()) { return *this; }
		if(line[0] == '#') { nread = 0; is_comment = true; return *this; }
		else { is_comment = false; }

		// parse line
		std::istringstream str(line.c_str());
		std::string token;
		int cnt = -1;
		nread = 0;
		while(str >> token) {
			cnt++;
			if(cnt > max_field) break;
			if(fields.find(cnt) == fields.end()) continue;

			switch(fields[cnt].type) {
				case INT: *(fields[cnt].v.v_int) = atoi(token.c_str()); break;
				case SHORT: *(fields[cnt].v.v_short) = atoi(token.c_str()); break;
				case DOUBLE: *(fields[cnt].v.v_double) = atof(token.c_str()); break;
				case FLOAT: *(fields[cnt].v.v_float) = atof(token.c_str()); break;
				case STRING: *(fields[cnt].v.v_string) = token; break;
			}
			nread++;
		}
		return *this;
	}

/*	itextstream &next2()
	{
		nread = -1;
		do {
			getline(f, line);
		} while(!ret_comments && (line[0] == '#' && !f.eof()));

		if(f.eof()) return *this;
		if(line[0] == '#') { nread = 0; is_comment = true; return *this; }
		else { is_comment = false; }

		// parse line
		std::istringstream str(line.c_str());
		std::string token;
		int cnt = -1;
		nread = 0;
		while(str >> token) {
			cnt++;
			if(cnt > max_field) break;
			if(fields.find(cnt) == fields.end()) continue;

			switch(fields[cnt].type) {
				case INT: *(fields[cnt].v.v_int) = atoi(token.c_str()); break;
				case SHORT: *(fields[cnt].v.v_short) = atoi(token.c_str()); break;
				case DOUBLE: *(fields[cnt].v.v_double) = atof(token.c_str()); break;
				case FLOAT: *(fields[cnt].v.v_float) = atof(token.c_str()); break;
				case STRING: *(fields[cnt].v.v_string) = token; break;
			}
			nread++;
		}
		return *this;
	}*/

	operator bool() { return nread != -1; }
	bool iscomment() { return is_comment; }
};





class otextstream : textstream {
public:
	std::ostream &f;
	bool start;
public:
	otextstream(std::ostream &f_) : f(f_), start(true) {}
};

template<typename T>
otextstream & operator <<(otextstream &o, const T &v)
{
	if(!o.start) { o.f << " "; }
	else { o.start = false; }
	o.f << v;
	return o;
}

class nl {
public:
};

inline otextstream & operator <<(otextstream &o, const nl &v)
{
	o.f << "\n";
	o.start = true;
	return o;
}

#endif
