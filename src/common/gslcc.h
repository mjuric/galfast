#ifndef __gslcc_h
#define __gslcc_h

#include <gsl/gsl_multifit.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>

#include <iostream>

#include <astro/util.h>

namespace gsl
{
	class matrix
	{
	protected:
		gsl_matrix *m;
	public:
		matrix(int r, int c)
			: m(NULL)
		{
			m = gsl_matrix_alloc(r, c);
		}

		matrix operator=(const matrix& w) { gsl_matrix_memcpy(m, w.m); return *this; }

		operator gsl_matrix*()       { return m; }
		operator const gsl_matrix*() const { return m; }

		~matrix()
		{
			if(m != NULL) { gsl_matrix_free(m); }
		}

		double &operator()(int r, int c)       { return m->data[r*m->tda + c]; }
		double  operator()(int r, int c) const { return m->data[r*m->tda + c]; }

		size_t cols() const { return m->size2; }
		size_t rows() const { return m->size1; }
		size_t size() const { return rows()*cols(); }
	private:
		matrix() : m(NULL) { }
		friend class matrix_wrap;
	};
	inline OSTREAM(const matrix &m)
	{
		FORj(i, 0, m.rows()) {
			FORj(j, 0, m.cols()) {
				if(j != 0) std::cout << ", ";
				std::cout << m(i, j);
			}
			std::cout << "\n";
		}
		return out;
	}

	class matrix_wrap : public matrix
	{
	public:
		matrix_wrap(gsl_matrix *m_) { m = m_; }
		~matrix_wrap() { m = NULL; }
	};

	class vector
	{
	protected:
		gsl_vector *v;
	public:
		vector operator=(const vector& w) { gsl_vector_memcpy(v, w.v); return *this; }

		operator gsl_vector*()       { return v; }
		operator const gsl_vector*() const { return v; }

		double &operator()(int i)       { return v->data[i*v->stride]; }
		double  operator()(int i) const { return v->data[i*v->stride]; }

		size_t size() const { return v->size; }

		vector(int n)
		{
			v = gsl_vector_alloc(n);
		}

		~vector()
		{
			if(v != NULL) { gsl_vector_free(v); }
		}
	private:
		vector() : v(NULL) {}
		friend class vector_wrap;
	};

	class vector_wrap : public vector
	{
	public:
		vector_wrap(gsl_vector *v_) { v = v_; }
		~vector_wrap() { v = NULL; }
	};




	class multifit
	{
	protected:
		gsl_multifit_linear_workspace *wspc;
		int n, p;

	public:
		multifit(size_t n_, size_t p_) : n(n_), p(p_) { wspc = gsl_multifit_linear_alloc(n, p); }
		~multifit() { gsl_multifit_linear_free(wspc); }

		int linear(const gsl::matrix &X, const gsl::vector &y, gsl::vector &c, gsl::matrix &cov, double &chisq)
		{
			return gsl_multifit_linear(X, y, c, cov, &chisq, wspc);
		}
			
		int wlinear(const gsl::matrix &X, const gsl::vector &w, const gsl::vector &y, gsl::vector &c, gsl::matrix &cov, double &chisq)
		{
			return gsl_multifit_wlinear(X, w, y, c, cov, &chisq, wspc);
		}
	};

};


#endif
