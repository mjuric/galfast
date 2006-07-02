#include "binarystream.h"
#include <iostream>
#include <ctime>
#include <astro/util.h>
#include <iostream>
#include "analysis.h"

using namespace std;

// this should be put into its own file
char spinner::chars[spinner::NCHARS+1] = "-/|\\";
void spinner::tick()
{
	if(!first) {
		char buf[3] = {8, chars[at], 0};
		std::cerr << buf;
	} else {
		first = false;
		std::cerr << chars[at];
	}

	at++;
	if(at == NCHARS) { at = 0; }
}
void spinner::stop()
{
	if(first) return;
	std::cerr << (char)8;
	first = true;
	at = 0;
	std::cerr << " ";
}
void spinner::start()
{
	stop();
}


// random magic number
const int header::magic = 0x37592664;
const int header::current_version = 1;

void ibinarystream::read(char *v, size_t n)
{
	f.read(v, n);
}

void obinarystream::write(const char *v, size_t n)
{
	f.write(v, n);
}

header::header(std::string description_)
: description(description_), datetime(time(NULL)), version(current_version)
{
}

header::header()
: description("Unititialized header"), datetime(0), version(-1)
{
}

obinarystream &operator <<(obinarystream &out, const header::data_map &data)
{
	out << data.size();
	FOREACH(data) { out << (*i).first << (*i).second; }
	return out;
}

ibinarystream &operator >>(ibinarystream &in, header::data_map &data)
{
	size_t size; string k, v;

	data.clear();

	in >> size;
	FOR(0, size) { in >> k >> v; data[k] = v; }
	return in;
}

obinarystream &operator <<(obinarystream &out, const header &h)
{
	out << header::magic << h.version << h.description << h.datetime << h.data;
	return out;
}

ibinarystream &operator >>(ibinarystream &in, header &h) throw (EBinaryIO)
{
	int magic = 0;
	in >> magic;
	if(magic != header::magic) { THROW(EBinaryIO, "This file does not start with a standard binary header. Perhaps the file has no header information, is compressed or corrupted?"); }

	in >> h.version >> h.description >> h.datetime >> h.data;

	return in;
}

OSTREAM(const header &h)
{
	cout << h.description << "\n\n";

	cout << "Header keywords:" << "\n";
	FOREACH(h.data) { cout << "    " << (*i).first << " = " << (*i).second << "\n"; }
	cout << "\n";
	
	cout << "File saved on " << ctime(&h.datetime) << "\n";
	cout << "Internal header version: " << h.version << "\n";
	cout << "This code can read headers up to version: " << header::current_version << "\n";
}
