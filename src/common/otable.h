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

#ifndef otable_h__
#define otable_h__

#include "column.h"

#include <astro/exceptions.h>
#include <astro/util.h>

#include <boost/shared_ptr.hpp>
#include <vector>

/**
	fmtout -- Formatter for textual serialization of otable.

	This is basically a type-safe(r) version of printf, that falls back
	to iostreams for non-builtin types. Probably an overkill given we
	usually only allow ints/doubles/floats/strings in otable anyway.
*/
class fmtout
{
protected:
	static const size_t BUFMAX = 20000;
	char buf[BUFMAX+1];	// line buffer
	size_t pos;
public:
	fmtout() : pos(0) {}
	const char *c_str() const { return buf; }

	size_t prep_buf()
	{
		if(pos == BUFMAX)
		{
			// This should really never happen ...
			buf[BUFMAX] = 0;
			THROW(peyton::exceptions::EAny, "Line buffer exhausted");
		}

		// Spaces between fields
		if(pos != 0) { buf[pos] = ' '; pos++; }
		return BUFMAX - pos;
	}

	template<typename T>
	int printf_aux(char *dest, size_t len, const char *fmt, const T &v)
	{
		// Default action: explode, because this has to be overloaded for
		// every printf-legal type
		THROW(peyton::exceptions::EAny, "Internal error");
	}

	int printf_aux(char *dest, size_t maxlen, const char *fmt, const float &v) 	{ pos += snprintf(dest, maxlen, fmt, v); }
	int printf_aux(char *dest, size_t maxlen, const char *fmt, const int &v) 	{ pos += snprintf(dest, maxlen, fmt, v); }
	int printf_aux(char *dest, size_t maxlen, const char *fmt, const double &v) 	{ pos += snprintf(dest, maxlen, fmt, v); }

	template<typename T>
	int printf(const std::string &fmt, const T &v)
	{
		size_t len = prep_buf();

		if(!fmt.size())
		{
			// No format specification -- revert to iostreams
			std::ostringstream ss;
			ss << v;
			std::string out = ss.str();
			strncpy(buf+pos, out.c_str(), std::min(len, out.size()));
			pos += out.size();
		}
		else
		{
			// sprintf format specified, use it
			printf_aux(buf+pos, len, fmt.c_str(), v);
		}
	}
};

/**
	column_type_traits -- Type traits of column types. These have to
	be specialized for each allowed column type (see the notes below)
*/
struct column_type_traits
{
	const std::string typeName;
	const size_t elementSize;

	virtual void  serialize(fmtout &out, const std::string &format, const void *val) const = 0;
	virtual void  unserialize(void *val, std::istream &in) const = 0;
	virtual void* constructor(void *p) const = 0;
	virtual void  destructor(void *val) const = 0;
	virtual char  fits_tform() const = 0;

	static const column_type_traits *get(const std::string &datatype);
	template<typename T> static const column_type_traits *get() { ASSERT(0); }
protected:
	static std::map<std::string, boost::shared_ptr<column_type_traits> > defined_types;
	column_type_traits(const std::string &name, const size_t size) : typeName(name), elementSize(size) {}
};
/**
	These are C type->traits mappings. Specialize them for each datatype supported by column_type_traits::get
*/
template<> const column_type_traits *column_type_traits::get<float>();
template<> const column_type_traits *column_type_traits::get<int>();
template<> const column_type_traits *column_type_traits::get<double>();
template<> const column_type_traits *column_type_traits::get<char>();

/**
*/
class otable
{
protected:
	class kv
	{
	public:
		std::string what;

		std::map<std::string, std::string> m_properties;	// arbitrary (key,value) pairs that the user/modules can set and inspect
	protected:
		friend class otable;

		virtual void set_property(const std::string &key, const std::string &value) = 0;
		virtual const std::string &get_property(const std::string &key, bool *exists = NULL) const = 0;

		virtual void serialize_def(std::ostream &out) const = 0;

		kv(const std::string &what_) : what(what_) {}
	};

	class columnclass : public kv
	{
	protected:
		friend class otable;

		otable &parent;
		std::string className;			// primary name of this class (e.g., "photometry", "color", "astrometry", ...)

		std::string formatString;		// default format of fields of this class
		const column_type_traits *typeProxy;	// default datatype of field of this class

		columnclass(otable &parent_);

		virtual void set_property(const std::string &key, const std::string &value);
		virtual void serialize_def(std::ostream &out) const
		{
			out << className << "{";

			// aliases
			FOREACH(parent.cclasses)
			{
				if(i->second.get() != this) { continue; }
				if(i->first == className) { continue; }
				out << "alias=" << i->first << ";";
			}

			// keywords
			if(!formatString.empty()) { out << "fmt=" << formatString << ";"; }

			// properties
			FOREACH(m_properties) { out << i->first << "=" << i->second << ";"; }

			out << "}";
		}

		virtual const std::string &get_property(const std::string &key, bool *exists = NULL) const
		{
			std::map<std::string, std::string>::const_iterator it = m_properties.find(key);
			if(it == m_properties.end())
			{
				if(exists != NULL) { *exists = false; }
				static const std::string empty;
				return empty;
			}
			if(exists != NULL) { *exists = true; }
			return it->second;
		}
	};

	std::string galfast_version;
public:
	struct columndef : public kv
	{
	protected:
		friend class otable;
		friend struct save_column_default;

		std::string columnName;				// primary name of the column

		const columnclass *columnClass;			// class of this column (note: it's never NULL)
		const column_type_traits *typeProxy;		// a proxy for type's serialization/construction/properties (note: must be accessed through type())
		std::string formatString;			// io::formatter format string of the column (note: must be accessed through getFormatString())

		bool m_hidden;
		struct {
			std::map<std::string, int> str2idx;
			std::map<int, std::string> idx2str;
		} fieldNames;

		otable &parent;					// parent table of this column

	protected:
		boost::shared_ptr<columndef> clone(const std::string &newColumnName, const std::map<int, std::string> *newFieldNames = NULL) const
		{
			boost::shared_ptr<columndef> c(new columndef(parent));
			c->columnName = newColumnName;

			c->ptr.reshape(ptr);

			c->columnClass = columnClass;
			c->formatString = formatString;
			c->typeProxy = typeProxy;

			c->m_hidden = m_hidden;

			c->setFieldNames(newFieldNames ? *newFieldNames : fieldNames.idx2str);

			c->m_properties = m_properties;

			return c;
		}

	public:
		virtual const std::string &get_property(const std::string &key, bool *exists = NULL) const
		{
			// return the requested property, if exists, delegate to class otherwise
			std::map<std::string, std::string>::const_iterator it = m_properties.find(key);
			if(it == m_properties.end())
			{
				return columnClass->get_property(key, exists);
			}
			if(exists != NULL) { *exists = true; }
			return it->second;
		}

		const std::string &getPrimaryName() const { return columnName; }

		const std::string &getFormatString() const
		{
			if(!formatString.empty()) { return formatString; }
			return columnClass->formatString;
		}

		bool hidden() const
		{
			return this->m_hidden;
		}
		void set_hidden(bool h)
		{
			this->m_hidden = h;
		}

		const column_type_traits *type() const
		{
			return typeProxy ? typeProxy : columnClass->typeProxy;
		}

		/*
			Note: The data is stored in 'structure of arrays format'. E.g., if
			this column is a vector of double[3], and the length of the table is 4,
			the memory layout is:
				11 21 31 41 xx
				12 22 32 42 xx
				13 23 33 43 xx
			where xx is some possible padding to ensure proper word-alignment
			of adjacent rows.

			The memory location of element j in row i, use:
				Aij = data + pitch*i + elementSize*j
		*/
 		column<char>	ptr;

		friend struct cmp_in;
		friend struct cmp_out;

	protected:
		columndef(otable &parent_);

		void alloc(const size_t len);
		void dealloc();

		virtual void set_property(const std::string &key, const std::string &value);
		template<typename T> column<T> &dataptr()
		{
			ASSERT(column_type_traits::get<T>() == type())
			{
				std::cerr << "Attempting to access a " << type()->typeName << " column as " << column_type_traits::get<T>()->typeName << "\n";
			}
			return (column<T> &)ptr;
		}
		template<typename T> const column<T> &dataptr() const
		{
			ASSERT(column_type_traits::get<T>() == type())
			{
				std::cerr << "Attempting to access a " << type()->typeName << " column as " << column_type_traits::get<T>()->typeName << "\n";
			}
			return (const column<T> &)ptr;
		}
	public:
		void *rawdataptr(int &elementSize, int &nfields, size_t &pitch)
		{
			nfields = ptr.width();
			pitch = ptr.pitch();
			elementSize = ptr.elementSize();

			column<char>::host_t h = ptr;
			return &h(0);
		}
		~columndef();
		columndef &add_alias(const std::string &name)		// add an additional name (alias) to this column
		{
			set_property("alias", name);
		}
		size_t getAliases(std::set<std::string> &result) const;	// get the list of all aliases of this column. Returns the number of aliases.
		void serialize(fmtout &line, const size_t row) const;	// write out the element at row row
		void unserialize(std::istream &in, const size_t row);	// read in an element into row row
		virtual void serialize_def(std::ostream &out) const;	// write out the definition of this column
		size_t capacity() const { return ptr.nrows(); }
		const std::string &className() const { static std::string empty; return columnClass ? columnClass->className : empty; }
		size_t setFieldNames(const std::map<int, std::string> &names);
		size_t getFieldNames(std::set<std::string> &names) const;
		size_t getFieldNames(std::map<int, std::string> &names) const;
		int getFieldIndex(const std::string &name) const
		{
			if(!fieldNames.str2idx.count(name)) { return -1; }
			return fieldNames.str2idx.at(name);
		}
	};

protected:
	friend struct cmp_in;
	friend struct cmp_out;
	friend struct record_loaded_columns;
	friend struct save_column_default;

	std::map<std::string, boost::shared_ptr<columnclass> > cclasses;
	std::map<std::string, boost::shared_ptr<columndef> > columns;
	size_t nrows_capacity;	// maximum number of rows in the table
	size_t nrows;	// rows actually in the table
	std::vector<std::string> colInput,	// columns to unserialize from file on the next unserialize_body call
				colOutput;	// columns to serialize to file, if they're in use
	std::vector<const columndef*>
				outColumns;	// columns to serialize (filled out by serialize_head, used by serialize_body)
public:
	size_t size() const { return nrows; }
	void set_size(size_t newsize)
	{
		nrows = newsize;
		if(nrows > capacity())
		{
			std::ostringstream ss;
			ss << "Attempted to set more rows than capacity allows. capacity=" << capacity() << ", nrows=" << nrows;
			THROW(peyton::exceptions::EAny, ss.str());
		}
	}

	size_t capacity() const { return nrows_capacity; }
	void clear() { nrows = 0; }
	size_t add_row()
	{
		size_t tmp = nrows++;
		if(nrows > capacity()) { THROW(peyton::exceptions::EAny, "Maximum number of rows in otable reached"); }
		return tmp;
	}

public:
	struct parse_callback { virtual bool operator()(kv *kvobj) = 0; };
	
protected:
	void init(const std::string &galfast_version);
	kv *parse(const std::string &defs, parse_callback *cback = NULL);

	void getColumnsForOutput(std::vector<const columndef*> &cols) const;	// aux helper
	void getColumnsForInput(std::vector<columndef*> &cols);			// aux helper

public:
	int getSortedColumnsForOutput(std::vector<const columndef*> &outColumns) const;

	// column lookup by name
	columndef &getColumn(const std::string &name);
	const columndef &getColumn(const std::string &name) const;
	void set_column_property(const std::string &name, const std::string &key, const std::string &value)
	{
		getColumn(name).set_property(key, value);
	}

	template<typename T> column<T>       &col(const std::string &name)       { return getColumn(name).dataptr<T>(); }
	template<typename T> const column<T> &col(const std::string &name) const { return getColumn(name).dataptr<T>(); }
	//bool have_column(const std::string &name) const { return columns.count(name); }
	bool using_column(const std::string &name) const;
	columndef &use_column(const std::string &coldef, bool setOutput = true);
	columndef &use_column_by_cloning(const std::string &newColumnName, const std::string &existingColumnName, std::map<int, std::string> *newFieldNames = NULL, bool setOutput = true);
	size_t get_used_columns(std::set<std::string> &cols) const;		// returns the list of columns in use
	size_t get_used_columns_by_class(std::set<std::string> &cols, const std::string &className) const;
	void alias_column(const std::string &column, const std::string &alias)
	{
		getColumn(column).add_alias(alias);
	}
	void drop_column(const std::string &name)
	{
		// remove this column from the list of columns.
		// NOTE: For now, please don't use this for nothing other than
		// dropping aliases.
		columns.erase(name);
	}
//	bool del_column(const std::string &name);

	otable(const size_t len, const std::string &galfast_version_)
	{
		nrows_capacity = len;
		nrows = 0;

		init(galfast_version_);
	}

	struct mask_functor { virtual bool shouldOutput(int row) const = 0; };
	struct default_mask_functor : public mask_functor { virtual bool shouldOutput(int row) const { return true; } };

	// serialization/unserialization routines
	std::ostream& serialize_header(std::ostream &out);
	std::istream& unserialize_header(std::istream &in, std::set<std::string> *columns = NULL);
	size_t serialize_body(std::ostream& out, size_t from = 0, size_t to = -1, const mask_functor &mask = default_mask_functor()) const;
	std::istream& unserialize_body(std::istream& in);
	size_t set_output(const std::string &colname, bool output);
	size_t set_output_all(bool output = true);
};

#endif // otable_h__
