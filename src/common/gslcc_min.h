#ifndef __gslcc_min_h
#define __gslcc_min_h

#include <gsl/gsl_min.h>
#include <gsl/gsl_errno.h>

#include "gslcc_exceptions.h"

namespace peyton {
namespace exceptions {

	GSL_EXCEPTION(EGSLMinimizer);

}
}


namespace gsl {

	double mmizer_function(double y, void *vdata);

	class mmizer : public gsl_function
	{
	protected:
		gsl_min_fminimizer *s;
		double errabs, errrel;
		int max_iter;
		double xmin, xmax;
		int status;

		typedef peyton::exceptions::EGSLMinimizer EXC;

	public:
		mmizer(double xmin_, double xmax_, double errabs_ = 0, double errrel_ = 1e-5)
			: xmin(xmin_), xmax(xmax_), max_iter(100), errabs(errabs_), errrel(errrel_), status(GSL_SUCCESS)
		{
			s = gsl_min_fminimizer_alloc(gsl_min_fminimizer_brent);

			function = &mmizer_function;
			params = this;
		}

		// function to minimize - overload to specify it
		virtual double fn_to_minimize(double x) = 0;

		double evaluate(double xguess)
		{
//			if(!(xmin < xguess && xguess < xmax)) {
//				std::cerr << "Error: " << xmin << " " << xmax << " " << xguess << "\n";
//			}
			status = gsl_min_fminimizer_set(s, this, xguess, xmin, xmax);
//			if(status) { std::cerr << xguess << "\n"; }
			GSLTHROW(status);

			int iter = 0;
			do
			{
				iter++;
				status = gsl_min_fminimizer_iterate(s);
				GSLTHROW(status);

				double m = gsl_min_fminimizer_x_minimum(s);
				double a = gsl_min_fminimizer_x_lower(s);
				double b = gsl_min_fminimizer_x_upper(s);

				status = gsl_min_test_interval(a, b, errabs, errrel);

				if(status == GSL_SUCCESS) break;
			}
			while (status == GSL_CONTINUE && iter < max_iter);
			GSLTHROW(status);

			return gsl_min_fminimizer_x_minimum(s);
		}

		bool good() { return status == GSL_SUCCESS; }

		virtual ~mmizer()
		{
			gsl_min_fminimizer_free(s);
		}
	};

	inline double mmizer_function(double y, void *vdata)
	{
		mmizer *data = (mmizer *)vdata;
		return data->fn_to_minimize(y);
	}

}

#endif
