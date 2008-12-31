#ifndef __gslcc_exceptions_h
#define __gslcc_exceptions_h

#include <astro/exceptions.h>

namespace peyton {
namespace exceptions {

	class EGSL : public EAny 
	{
	public:
		int status;
	public:
		EGSL(int status, const std::string &fun, const std::string &f, int l, const std::string &st);
	};

	#define GSLTHROW1(ex, status) if((status)) { throw ex((status), __PRETTY_FUNCTION__, __FILE__, __LINE__, peyton::system::stacktrace()); }
	#define GSLTHROW(status) if((status)) { throw EXC((status), __PRETTY_FUNCTION__, __FILE__, __LINE__, peyton::system::stacktrace()); }

	#define GSL_EXCEPTION(ex) \
		class ex : public EGSL { \
		public: \
			ex(int status, const std::string &fun, const std::string &f, int l, const std::string &st) \
				: EGSL(status, fun, f, l, st) {} \
		}

	inline EGSL::EGSL(int status, const std::string &fun, const std::string &f, int l, const std::string &st)
		: EAny(gsl_strerror(status), fun, f, l, st)
	{}

}
}

#endif
