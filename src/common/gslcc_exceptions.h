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
		EGSL(int status, std::string fun, std::string f, int l);
	};

	#define GSLTHROW1(ex, status) if((status)) { throw ex((status), __PRETTY_FUNCTION__, __FILE__, __LINE__); }
	#define GSLTHROW(status) if((status)) { throw EXC((status), __PRETTY_FUNCTION__, __FILE__, __LINE__); }

	#define GSL_EXCEPTION(ex) \
		class ex : public EGSL { \
		public: \
			ex(int status, std::string fun, std::string f, int l) \
				: EGSL(status, fun, f, l) {} \
		}

	inline EGSL::EGSL(int status, std::string fun, std::string f, int l)
		: EAny(gsl_strerror(status), fun, f, l)
	{}

}
}

#endif
