#ifndef _paralax_h__
#define _paralax_h__

#include <vector>
#include <cmath>

#include <gsl/gsl_poly.h>

#include <astro/io/binarystream.h>

extern double Rg;	// galactic center distance

namespace stardist
{
	inline double logD(float m, double M) { return 1. + (m - M)/5.; }		///< log10 distance, based on absolute and apparent i magnitudes
	inline double D(float m, double M) { return std::pow(10.0, logD(m, M)); }	///< distance, based on absolute and apparent i magnitudes
	inline float m(double D, double M) { return M + 5*log10(D/10); }
};

class plx_gri_locus_ng
{
protected:
	std::vector<double> Mrc;
public:
	void setParalaxCoefficients(const std::vector<double> &Mcoef);

	double Mr(const float ri) const
	{
		return gsl_poly_eval(&Mrc[0], Mrc.size(), ri);
	}

	friend inline BOSTREAM2(const plx_gri_locus_ng &plx) { return out << plx.Mrc; }
	friend inline BISTREAM2(plx_gri_locus_ng &plx) { return in >> plx.Mrc; }
};

#endif // _paralax_h__
