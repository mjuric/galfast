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

#include <iostream>
#include <boost/iostreams/filtering_stream.hpp>

class flex_output
{
protected:
	std::ostream *stream;
	boost::iostreams::filtering_streambuf<boost::iostreams::output> *sbout;

public:
	flex_output(const std::string &fn = "") { open(fn); }
	~flex_output();

	std::ostream *open(const std::string &fn);
	std::ostream &out() { return *this->stream; }
};

class flex_input
{
protected:
	std::istream *stream;
	boost::iostreams::filtering_streambuf<boost::iostreams::input> *sbin;

public:
	flex_input(const std::string &fn = "") { open(fn); }
	~flex_input();

	std::istream *open(const std::string &fn);
	std::istream &in() { return *this->stream; }
};
