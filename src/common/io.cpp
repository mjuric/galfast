/***************************************************************************
 *   Copyright (C) 2004 by Mario Juric                                     *
 *   mjuric@astro.Princeton.EDU                                            *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 *   This program is distributed in the hope that it will be useful,       *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
 *   GNU General Public License for more details.                          *
 *                                                                         *
 *   You should have received a copy of the GNU General Public License     *
 *   along with this program; if not, write to the                         *
 *   Free Software Foundation, Inc.,                                       *
 *   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.             *
 ***************************************************************************/


#include <boost/iostreams/filter/bzip2.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/iostreams/device/file.hpp>

#include <astro/system/log.h>

#include "io.h"

flex_output::~flex_output()
{
	if(stream != &std::cout)
	{
		delete stream;
	}
	delete sbout;
}

std::ostream *flex_output::open(const std::string &fn)
{
	if(!fn.size()) { return NULL; }
	
	using namespace boost::iostreams;

	stream = NULL; sbout = NULL;

	if(fn.size() > 4 && fn.rfind(".bz2") == fn.size()-4)
	{
		sbout = new filtering_streambuf<output>(bzip2_compressor() | file_sink(fn));
		stream = new std::ostream(sbout);
		MLOG(verb2) << "Outputing bzip2 compressed data to " << fn << ".";
	}
	else if(fn.size() > 3 && fn.rfind(".gz") == fn.size()-3)
	{
		sbout = new filtering_streambuf<output>(gzip_compressor() | file_sink(fn));
		stream = new std::ostream(sbout);
		MLOG(verb2) << "Outputing gzip compressed data to " << fn << ".";
	}
	else if(fn == "-")
	{
		MLOG(verb2) << "Outputing plain text to standard output.";
		stream = &std::cout;
	}
	else
	{
		MLOG(verb2) << "Outputing plain text to " << fn << ".";
		stream = new std::ofstream(fn.c_str());
	}

	return stream;
}




flex_input::~flex_input()
{
	if(stream != &std::cin)
	{
		delete stream;
	}
	delete sbin;
}

std::istream *flex_input::open(const std::string &fn)
{
	if(!fn.size()) { return NULL; }
	
	using namespace boost::iostreams;

	stream = NULL; sbin = NULL;

	if(fn.size() > 4 && fn.rfind(".bz2") == fn.size()-4)
	{
		sbin = new filtering_streambuf<input>(bzip2_decompressor() | file_source(fn));
		stream = new std::istream(sbin);
		MLOG(verb2) << "Reading bzip2 compressed data from " << fn << ".";
	}
	else if(fn.size() > 3 && fn.rfind(".gz") == fn.size()-3)
	{
		sbin = new filtering_streambuf<input>(gzip_decompressor() | file_source(fn));
		stream = new std::istream(sbin);
		MLOG(verb2) << "Reading gzip compressed data from " << fn << ".";
	}
	else if(fn == "-")
	{
		MLOG(verb2) << "Reading plain text from standard input.";
		stream = &std::cin;
	}
	else
	{
		MLOG(verb2) << "Reading plain text from " << fn << ".";
		stream = new std::ifstream(fn.c_str());
	}

	return stream;
}

