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

#include "config.h"

#include "otable.h"

#include <astro/useall.h>

///////////////////////////////////////////////
otable::kv *otable::parse(const std::string &defs, otable::parse_callback *cback)
{
	// parse stringDef and create a corresponding tag instance
	// the tag definition format is: ((class|column)) name [n] { field1=value1; field2=value2; ...}
	std::istringstream ss(defs);
	char c;
	kv *kvobj = NULL;
	std::string what, name;
	int ncolumn = 0;
	while(ss >> c)
	{
		// find out the metatype of object being parsed (class or column).
		ss.unget();
		if(c == '(') { ss >> what; }
		else { what = "(column)"; }

		// load the name of the object being parsed
		name.clear();
		ss >> c;	// eats any whitespace to next character
		while(ss && (isalnum(c) || c == '_' || c == ':')) { name += c; ss.get(c); }
		ss.unget();
		if(name.empty())
		{
			THROW(EAny, "Error reading column name at character " + str((size_t)ss.tellg()) + " of the line.");
		}

		// get or instantiate this object
		if(what == "(class)")
		{
			if(!cclasses.count(name))
			{
				cclasses[name].reset(new columnclass(*this));
			}
			kvobj = cclasses[name].get();
		}
		else if(what == "(column)")
		{
			if(!columns.count(name))
			{
				columns[name].reset(new columndef(*this));
			}
			kvobj = columns[name].get();

/*			if(setInputOn)
			{
				columns[name]->input = ncolumn;
				ncolumn++;
			}*/
		}
		else
		{
			THROW(EAny, "Expected 'class' or 'column', got " + what);
		}
		kvobj->set_property("__name__", name);

		ss >> c;	// look for '[' or '{'
		// special (traditional) syntax for arrays
		if(c == '[' && what == "(column)")
		{
			size_t width;
			ss >> width;
			kvobj->set_property("n", str(width));
//			std::cerr << what << " " << name << ": [] =" << width << "\n";

			ss >> c;	// look for ']'
			if(c != ']') { THROW(EAny, "Expected ']', got " + str(c)); }

			ss >> c;	// look for '{'
		}

		if(c != '{')
		{
			if(what == "(column)") { ss.unget(); if(cback) (*cback)(kvobj); continue; }	// bare column defininiton (without details)
			THROW(EAny, "Expected '{', got " + str(c));
		}

		ss >> c;
		while(c != '}')
		{
			ss.unget();

			std::string key;
			ss >> c; while(ss && (isalnum(c) || c == '_')) { key += c; ss.get(c); }
			if(!ss) { THROW(EAny, "End of file while reading key name"); }
			if(c != '=') {
				THROW(EAny, "Expected '=', got " + str(c));
			}

			std::string value;
			ss >> c; while(ss && (c != ';' && c != '}')) { value += c; ss.get(c); }
			if(!ss) { THROW(EAny, "End of file while reading field name"); }
			if(c == '}') { ss.unget(); }

//			std::cerr << what << " " << name << ": " << key << "=" << value << "\n";
			kvobj->set_property(key, value);

			ss >> c;
		}

		if(cback) (*cback)(kvobj);
	}
	return kvobj;
}

template<typename T>
struct default_column_type_traits : public column_type_traits
{
	const char m_tform_code;

	virtual void  serialize(fmtout &out, const std::string &format, const void *val) const { const T *v = reinterpret_cast<const T*>(val); out.printf(format, *v); }
	virtual void  unserialize(void *val, std::istream &in) const { T *v = reinterpret_cast<T*>(val); in >> *v; }
	virtual void* constructor(void *p) const { return new (p) T(); }
	virtual void  destructor(void *val) const { reinterpret_cast<T*>(val)->~T(); }
	virtual char  fits_tform() const { return m_tform_code; }

	default_column_type_traits(const std::string &name, const char tform_code) : column_type_traits(name, sizeof(T)), m_tform_code(tform_code) {}
};

// These are C type->traits mappings. Specialize them for each datatype supported by column_type_traits::get
// and declare them in model.h (or else it won't work!!)
template<> const column_type_traits *column_type_traits::get<float>()  { return column_type_traits::get("float"); }
template<> const column_type_traits *column_type_traits::get<int>()    { return column_type_traits::get("int"); }
template<> const column_type_traits *column_type_traits::get<double>() { return column_type_traits::get("double"); }
template<> const column_type_traits *column_type_traits::get<char>()   { return column_type_traits::get("char"); }

std::map<std::string, boost::shared_ptr<column_type_traits> > column_type_traits::defined_types;
const column_type_traits *column_type_traits::get(const std::string &datatype)
{
	static bool initialized = false;
	if(!initialized)
	{
		#define ADDTYPE(strT, fitsT, T) defined_types[strT].reset(new default_column_type_traits<T>(strT, fitsT));
		ADDTYPE("int",    'J', int);
		ADDTYPE("double", 'D', double);
		ADDTYPE("char",   'A', char);
		ADDTYPE("float",  'E', float);
		#undef CREATETYPE
		initialized = true;
	}

	if(!defined_types.count(datatype)) { THROW(EAny, "Unknown tag data type '" + datatype + "'"); }
	return defined_types[datatype].get();
}

otable::columnclass::columnclass(otable &parent_)
	: parent(parent_), kv("(class)")
{
	typeProxy = column_type_traits::get("float");	// default column type
	ASSERT(typeProxy);
}

otable::columndef::~columndef()
{
	dealloc();
}

otable::columndef::columndef(otable &parent_)
	: parent(parent_), kv("(column)")
{
	// defaults
	columnClass = parent.cclasses["default"].get();
	typeProxy = NULL;	// default to class type
	m_hidden = false;		// default to outputing the column
}

void otable::columndef::alloc(const size_t nrows)
{
	// do nothing if the length is OK
	if(nrows == ptr.nrows()) { return; }

	dealloc();

	const column_type_traits *tt = type();
	size_t elementSize = tt->elementSize;
	ptr.resize(nrows, ptr.width(), elementSize);

	// call constructors
	char *base = ptr.get();
	size_t pitch = ptr.pitch();
	for(size_t i = 0; i != ptr.width(); i++)
	{
		for(size_t j = 0; j != ptr.nrows(); j++)
		{
			void *Aij = base + pitch*i + elementSize*j;
			tt->constructor(Aij);
		}
	}
}

void otable::columndef::dealloc()
{
	if(!ptr.size()) { return; }

	// call destructors
	const column_type_traits *tt = type();
	size_t elementSize = ptr.elementSize();
	size_t pitch = ptr.pitch();
	char *base = ptr.get();
	for(size_t i = 0; i != ptr.width(); i++)
	{
		for(size_t j = 0; j != ptr.nrows(); j++)
		{
			void *Aij = base + pitch*i + elementSize*j;
			tt->destructor(Aij);
		}
	}

	// deallocate the memory without erasing the metadata
	ptr.resize(0);
}

void otable::columnclass::set_property(const std::string &key, const std::string &value)
{
	if(key == "__name__")
	{
		if(className.empty()) { className = value; }
		return;
	}

	if(key == "fmt") { formatString = value; return; }

	if(key == "type")
	{
		typeProxy = column_type_traits::get(value);
		ASSERT(typeProxy);
		return;
	}

	// default: store the property in m_properties map
	m_properties[key] = value;
}

void otable::columndef::set_property(const std::string &key, const std::string &value)
{
	if(key == "fmt") { formatString = value; return; }

	if(key == "__name__")
	{
		ASSERT(columnName.empty() || columnName == value);
		columnName = value;

		return;
	}

	if(key == "alias")
	{
		if(!parent.columns.count(value))
		{
			ASSERT(!columnName.empty());
			parent.columns[value] = parent.columns[columnName];
		}
		ASSERT(parent.columns[value] == parent.columns[columnName]);
		return;
	}

	if(key == "class")
	{
		ASSERT(parent.cclasses.count(value));
		dealloc();
		columnClass = parent.cclasses[value].get();
		return;
	}

	if(key == "n") // vector width
	{
		dealloc();
		int width = atoi(value.c_str());
		ASSERT(width > 1);
		ptr.resize(ptr.nrows(), width);
		return;
	}

	if(key == "type")
	{
		dealloc();
		typeProxy = column_type_traits::get(value);
		return;
	}

	if(key == "hidden")
	{
		m_hidden = value == "true" || atoi(value.c_str()) != 0;
		return;
	}

	if(key == "fieldNames")
	{
		// value = "idx:fieldname,idx:fieldname,..."
		size_t at = 0;
		size_t len;
		do
		{
			len = value.find(",", at);
			if(len != std::string::npos) { len -= at; }
			std::string pair = value.substr(at, len);
			at += len+1;

			int semi = pair.find(':');
			ASSERT(semi != std::string::npos);

			std::string sidx = pair.substr(0, semi);
			std::string name = pair.substr(semi+1);
			int idx = atoi(sidx.c_str());

			fieldNames.str2idx[name] = idx;
			fieldNames.idx2str[idx] = name;
		} while(len != std::string::npos);
		return;
	}

	// default: store the property in m_properties map
	m_properties[key] = value;
}

void otable::columndef::serialize(fmtout &line, const size_t row) const
{
	const column_type_traits *tt = type();
// 	char *at = (char*)data + tt->elementSize*row;
	const char *at = ((column<char> &)ptr).get() + tt->elementSize*row;
	const std::string &fmt = getFormatString();

	FOR(0, ptr.width())
	{
		tt->serialize(line, fmt, at);
		at += ptr.pitch();
	}
}
void otable::columndef::unserialize(std::istream &in, const size_t row)
{
	const column_type_traits *tt = type();
//	char *at = (char*)data + tt->elementSize*row;
	char *at = ptr.get() + tt->elementSize*row;
	FOR(0, ptr.width())
	{
		tt->unserialize(at, in);
		at += ptr.pitch();
	}
}

size_t otable::columndef::setFieldNames(const std::map<int, std::string> &names)
{
	fieldNames.str2idx.clear();
	fieldNames.idx2str.clear();

	FOREACH(names)
	{
		assert(!fieldNames.idx2str.count(i->first));
		assert(!fieldNames.str2idx.count(i->second));

		fieldNames.idx2str[i->first] = i->second;
		fieldNames.str2idx[i->second] = i->first;
	}
}

size_t otable::columndef::getFieldNames(std::map<int, std::string> &names) const
{
	names = fieldNames.idx2str;
}

size_t otable::columndef::getFieldNames(std::set<std::string> &names) const
{
	FOREACH(fieldNames.str2idx)
	{
		names.insert(i->first);
	}
}

size_t otable::columndef::getAliases(std::set<std::string> &result) const
{
	// aliases
	result.clear();
	FOREACH(parent.columns)
	{
		if(i->second.get() != this) { continue; }
		if(i->first == columnName) { continue; }
		result.insert(i->first);
	}

	return result.size();
}

void otable::columndef::serialize_def(std::ostream &out) const
{
	out << columnName;
	if(ptr.width() > 1) { out << "[" << ptr.width() << "]"; }

	const otable::columndef *dflt = NULL;
	if(parent.columns.count("default::" + columnName))
	{
		dflt = parent.columns.at("default::" + columnName).get();
	}

	std::stringstream ss;
	// keywords that have changed from their defaults
	#define DFLT(var) (dflt && dflt->var == var)
	if(typeProxy                           && !DFLT(typeProxy))    { ss << "type=" << typeProxy->typeName << ";"; }
	if(columnClass->className != "default" && !DFLT(columnClass))  { ss << "class=" << columnClass->className << ";"; }
	if(!formatString.empty()               && !DFLT(formatString)) { ss << "fmt=" << formatString << ";"; }
	if(dflt                                && !DFLT(m_hidden))     { ss << "hidden=" << m_hidden << ";"; }
	#undef DFLT

	// fieldNames
	if(fieldNames.idx2str.size())
	{
		ss << "fieldNames=";
		bool first = true;
		FOREACH(fieldNames.idx2str)
		{
			if(!first) { ss << ","; }
			ss << i->first << ":" << i->second;
			first = false;
		}
		ss << ";";
	}

	// aliases
	FOREACH(parent.columns)
	{
		if(i->second.get() != this) { continue; }
		if(i->first == columnName) { continue; }
		ss << "alias=" << i->first << ";";
	}

	// properties
	FOREACH(m_properties)
	{
		const std::string &key = i->first, &val = i->second;

		// ignore if this property is there by default
		if(dflt && dflt->get_property(key) == val) { continue; }

		ss << i->first << "=" << i->second << ";";
	}

	// output details only if they differ from defaults
	std::string details = ss.str();
	if(!details.empty())
	{
		out << "{" << details << "}";
	}
}

void otable::getColumnsForOutput(std::vector<const columndef*> &outColumns) const
{
	outColumns.clear();
	FOREACH(colOutput)
	{
		if(!columns.count(*i)) continue;		// don't display columns that have not been created
		if(columns.at(*i)->hidden()) continue;		// don't display hidden columns

		outColumns.push_back(columns.at(*i).get());
	}
}


struct pred_col_less
{
	std::map<std::string, int> &field;

	pred_col_less(std::map<std::string, int> &f) : field(f) { }

	int get_rank(const otable::columndef *a) const
	{
		// find alias with the smallest rank
		std::map<std::string, int>::iterator it = field.find(a->getPrimaryName());
		int rank = it != field.end() ? it->second : 1000;

		std::set<std::string> aliases;
		a->getAliases(aliases);
		FOREACH(aliases)
		{
			it = field.find(*i);
			if(it == field.end()) { continue; }
			rank = std::min(rank, it->second);
		}
		return rank;
	}

	bool operator()(const otable::columndef *a, const otable::columndef *b) const
	{
		int ranka = get_rank(a), rankb = get_rank(b);
		bool altb = ranka == rankb ? strcasecmp(a->getPrimaryName().c_str(), b->getPrimaryName().c_str()) < 0 : ranka < rankb;
/*		std::cout << "rank(" << a->getPrimaryName() << ")=" << ranka << "  ";
		std::cout << "rank(" << b->getPrimaryName() << ")=" << rankb << "  ";
		std::cout << (altb ? "a < b" : "a > b") << "\n";*/
		return altb;
	}
};

int otable::getSortedColumnsForOutput(std::vector<const columndef*> &outColumns) const
{
	getColumnsForOutput(outColumns);

	// sort the output columns into a sane
	std::map<std::string, int> field;
	int ord = 1;
	field["lb"] = ord++;
	field["radec"] = ord++;
	field["XYZ"] = ord++;
	field["DM"] = ord++;
	field["absmag"] = ord++;
	field["comp"] = ord++;
	field["FeH"] = ord++;
	field["vcyl"] = ord++;
	field["pmlb"] = ord++;
	field["pmradec"] = ord++;
	field["Am"] = ord++;
	std::sort(outColumns.begin(), outColumns.end(), pred_col_less(field));

	return outColumns.size();
}
	
std::ostream& otable::serialize_header(std::ostream &out)
{
	getSortedColumnsForOutput(outColumns);

	FOREACH(outColumns)
	{
		(*i)->serialize_def(out);
		out << " ";
	}
	return out;
}

// serialization/unserialization routines
size_t otable::serialize_body(std::ostream& out, size_t from, size_t to, const mask_functor &mask) const
{
	ASSERT(from >= 0);
	if(to > size()) { to = size(); }

	size_t cnt = 0;
	FORj(row, from, to)
	{
		if(!mask.shouldOutput(row)) { continue; }
		cnt++;

		fmtout line;
		FOREACH(outColumns)
		{
			(*i)->serialize(line, row);
		}
		out << line.c_str() << "\n";
	}
	return cnt;
};

void otable::getColumnsForInput(std::vector<columndef*> &inColumns)
{
	FOREACH(colInput)
	{
		ASSERT(columns.count(*i));
		inColumns.push_back(&getColumn(*i));
	}
}

struct record_loaded_columns : public otable::parse_callback
{
	std::vector<otable::columndef *> columns;

	virtual bool operator()(otable::kv *kvobj)
	{
//		std::cout << "META: " << kvobj->what << "\n";
		if(kvobj->what != "(column)") return true;
		columns.push_back(static_cast<otable::columndef*>(kvobj));

		return true;
	}
};

std::istream& otable::unserialize_header(std::istream &in, std::set<std::string> *columns)
{
	// gobble-up an optional the comment sign
	char c;
	in >> c;
	if(c != '#') { in.unget(); }

	// gobble up the rest of the line, where we expect the definitions to sit
	std::string line;
	getline(in, line);

	record_loaded_columns lc;
	parse(line, &lc);
	FOREACH(lc.columns)
	{
		colInput.push_back((*i)->columnName);
		if(columns) { columns->insert((*i)->columnName); }
	}

	colOutput = colInput;	// Output all unserialized columns by default

	return in;
}

std::istream& otable::unserialize_body(std::istream& in)
{
	std::vector<columndef*> inColumns;
	getColumnsForInput(inColumns);

	nrows = 0;
	FORj(row, 0, capacity())	// read at most length rows
	{
		bool first = true;
		FOREACH(inColumns)
		{
			(*i)->unserialize(in, row);
			if(!in) { break; }
			first = false;
		}
		if(!in) {
			if(!first) { THROW(EAny, "Incomplete last line."); }
			break;
		}
		nrows++;
	}

	if(!in && !in.eof()) { THROW(EAny, "Error after reading " + str(nrows) + " rows."); }

	return in;
};

size_t otable::set_output(const std::string &colname, bool output)
{
	std::vector<std::string>::iterator it = find(colOutput.begin(), colOutput.end(), colname);
	if(it != colOutput.end())
	{
		if(output) { return it - colOutput.begin(); }
		colOutput.erase(it);
		return -1;
	}

	if(output)
	{
		colOutput.push_back(colname);
		return colOutput.size()-1;
	}

	return -1;
}

size_t otable::set_output_all(bool output)
{
	// output all columns that are in use
	colOutput.clear();
	FOREACH(columns)
	{
		if(i->first != i->second->columnName) { continue; }
		if(i->second->capacity() == 0) { continue; }
		set_output(i->first, output);
	}
}

otable::columndef &otable::use_column(const std::string &coldef, bool setOutput)
{
	columndef *col = (columndef *)parse("(column) " + coldef);
	ASSERT(col != NULL);

	col->alloc(capacity()); // Ensure the column is allocated

	if(setOutput) { set_output(col->columnName, true); }

	return *col;
}

otable::columndef &otable::use_column_by_cloning(const std::string &newColumnName, const std::string &existingColumnName, std::map<int, std::string> *newFieldNames, bool setOutput)
{
	ASSERT(columns.count(newColumnName) == 0);

	columndef &exCol = getColumn(existingColumnName);
	boost::shared_ptr<columndef> col = columns[newColumnName] = exCol.clone(newColumnName, newFieldNames);

	col->alloc(capacity()); // Ensure the column is allocated

	if(setOutput) { set_output(col->columnName, true); }

	return *col;
}

void otable::init()
{
	// definition of built-in classes and column defaults
	parse(
	"(class) default      {fmt=% 7.3f;}" 			// NOTE: This class must come be defined before any columns are ever instantiated
	"(class) magnitude    {fmt=% 7.3f;}" 			// -12.345
	"(class) color        {fmt=% 6.3f;}"			// -12.345
	"(class) astrometry   {fmt=% 13.8f; type=double;}"	// -123.12345678
	"(class) position     {fmt=% 10.2f;}"			// -123456.78
	"(class) propermotion {fmt=% 8.2f;}"			// -1234.12
	"(class) velocity     {fmt=% 7.1f;}"			// -1234.1
	"(class) flags        {fmt=% 4d; type=int;}"		// 1234

	// definition of built-in fields
	"(column) comp          {type=int; fmt=%3d;}"
	"(column) radec[2]      {class=astrometry;}"
	"(column) lb[2]         {class=astrometry;}"
	"(column) XYZ [3]       {class=position;}"
	"(column) DM            {class=magnitude;}"
	"(column) FeH           {fmt=% 6.3f;}"
	"(column) vcyl[3]       {class=velocity;}"
	"(column) pmlb[3]       {class=propermotion;}"
	"(column) pmradec[3]    {class=propermotion;}"
	"(column) star_name[40] {type=char;}"
	"(column) hidden	{type=int;hidden=true;}"
	"(column) projIdx       {type=int;hidden=true;}"
	"(column) projXY[2]     {type=float;hidden=true;}"
	"(column) Am            {type=float;fmt=%5.2f;}"
	);

	// store these column definitions as defaults
	typeof(columns) tmp(columns);
	FOREACH(tmp)
	{
		if(i->second->columnName != i->first) { continue; } // skip aliases

		std::string newName = "default::" + i->first;
		columns[newName] = i->second->clone(newName);
	}
}

bool otable::using_column(const std::string &name) const
{
	if(!columns.count(name)) return false;
	return columns.at(name)->capacity() != 0;
}

size_t otable::get_used_columns(std::set<std::string> &cols) const
{
	cols.clear();
	FOREACH(columns)
	{
		if(i->second->capacity() == 0) { continue; }
		cols.insert(i->first);
	}
	return cols.size();
}

size_t otable::get_used_columns_by_class(std::set<std::string> &cols, const std::string &className) const
{
	cols.clear();
	FOREACH(columns)
	{
		if(i->second->capacity() == 0) { continue; }
		if(i->second->className() != className) { continue; }
		cols.insert(i->first);
	}
	return cols.size();
}

otable::columndef &otable::getColumn(const std::string &name)
{
	// Auto-create if needed
	if(!columns.count(name))
	{
		// autocreate
		use_column(name);
	}
	ASSERT(columns.count(name))
	{
		std::cerr << "Column " << name << " doesn't exist and couldn't be autocreated?!";
	}
 	columndef &col = *columns[name].get();

	// Auto-create column data
	if(col.capacity() != capacity())
	{
		col.alloc(capacity());
		ASSERT(col.capacity() == capacity());
	}

	return col;
}

