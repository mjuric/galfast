
//--- libpeyton includes. Add any libpeyton includes _before_ including
//--- astro/useall.h

#include "config.h"

#include <astro/system/options.h>
#include <astro/math.h>
#include <astro/io/format.h>
#include <astro/useall.h>
#include <astro/system/fs.h>

#include "textstream.h"
#include "analysis.h"
#include "paralax.h"

#include <iostream>
#include <fstream>
#include <iomanip>

using namespace std;

void ml_magnitudes(float &gml, float &rml, float &iml,
	const float g, const float r, const float i,
	const float gErr, const float rErr, const float iErr,
	const float ml_ri, const float ml_gr)
{
	float wg = 1/sqr(gErr);
	float wr = 1/sqr(rErr);
	float wi = 1/sqr(iErr);
	
	rml = (wr*r + wg*(g-ml_gr) + wi*(i+ml_ri)) / (wr + wi + wg);
	iml = rml - ml_ri;
	gml = rml + ml_gr;
}

int main(int argc, char **argv)
{
try
{
	VERSION_DATETIME(version, "$Id: locusfit.cpp,v 1.4 2006/07/13 23:27:30 mjuric Exp $");

	Options opts(
		argv[0],
		"This program has not been described",
		version,
		Authorship::majuric
	);

	//# add any arguments your program needs. eg:
	opts.argument("g_column", "The column index of g magnitude (1 based)");
	opts.argument("r_column", "The column index of r magnitude (1 based)");
	opts.argument("i_column", "The column index of i magnitude (1 based)");
	opts.argument("gErr_column", "The column index of g magnitude error (1 based)");
	opts.argument("rErr_column", "The column index of r magnitude error (1 based)");
	opts.argument("iErr_column", "The column index of i magnitude error (1 based)");
	opts.argument("Ar_column", "The column index of r band extinction (1 based, set to 0 if the extinction is already taken into account.)");

	// add any options your program might need. eg:
	// opts.option("meshFactor", "meshFactor", 0, "--", Option::required, "4", "Resolution decrease between radial steps");
	opts.option("prior").value("ri_prior.txt").param_optional().desc("Applies a prior along the locus. The parameter is the file containing the P(r-i) prior");

	parse_options(opts, argc, argv);

	/////// Start your application code here
	gsl_set_error_handler_off ();

//	plx_gri_locus as;

	// load density prior
	string priorfile = opts["prior"];
	if(priorfile.size())
	{
		if(!file_exists(priorfile))
		{
			cerr << priorfile << " does not exist. Aborting.\n";
			exit(-1);
		}
		
		paralax.mlri.setprior(priorfile);
		cout << "# Using prior from " << priorfile << "\n";

#if 0
		plx_gri_locus as0;
		{
			double y0 = 0.07;
			double x0 = as.gr(y0) + 0.05;
			double s = .05;

			float lnL;
			double ri1 = as0.mlri(y0, x0, s, s, s, &lnL);
			double prob1 = -as0.mlri.fn_to_minimize(ri1);
			double ri2 = as.mlri(y0, x0, s, s, s);
			double prob2 = -as.mlri.fn_to_minimize(ri2);
			double norm = prob1 - prob2;
			cerr << "# ml_ri(noprior) = " << ri1 << ", ml_ri(prior) = " << ri2 << "\n";

			for(double ri = y0 - 0.2; ri < y0 + 0.2; ri+=0.001)
			{
				double prob1 = -as0.mlri.fn_to_minimize(ri);
				double prob2 = -as.mlri.fn_to_minimize(ri);
				prob2 += norm;
				cout << ri << " " << prob1 << " " << prob2 << "\n";
			}
		}
		return 0;
#endif
	} else {
		cout << "# Using uniform prior along the locus\n";
	}

	int g_column = (int)opts["g_column"] - 1;
	int r_column = (int)opts["r_column"] - 1;
	int i_column = (int)opts["i_column"] - 1;
	int gErr_column = (int)opts["gErr_column"] - 1;
	int rErr_column = (int)opts["rErr_column"] - 1;
	int iErr_column = (int)opts["iErr_column"] - 1;
	int Ar_column = (int)opts["Ar_column"] - 1;

	itextstream in(cin);
	in.returncomments(true);

	float g, r, i, gErr, rErr, iErr, Ar;
	bind(in, g, g_column, r, r_column, i, i_column, gErr, gErr_column, rErr, rErr_column, iErr, iErr_column);
	if(Ar_column != -1) { bind(in, Ar, Ar_column); }
	else { Ar = 0; }

	cout << setprecision(6);

	float ml_ri, ml_gr, lnL;
	while(in.next())
	{
		if(in.iscomment()) { cout << in.line << "\n"; continue; }

		g -= Ar*3.793/2.751;
		r -= Ar*1;
		i -= Ar*2.086/2.751;
		
		//cerr << g << " " << r << " " << i << " " << gErr << " " << rErr << " " << iErr << " " << r - i << " " << g - r << "\n";

		try {
			ml_ri = paralax.mlri(r - i, g - r, gErr, rErr, iErr, &lnL);
			ml_gr = paralax.gr(ml_ri);
			ml_magnitudes(g, r, i, g, r, i, gErr, rErr, iErr, ml_ri, ml_gr);

			g += Ar*3.793/2.751;
			r += Ar*1;
			i += Ar*2.086/2.751;
		}
		catch(EGSLMinimizer &e)
		{
			// the star is unfittable
			g = r = i = ml_ri = ml_gr = -1;
		}

		cout << in.line << ' ' << g << ' ' << r << ' ' << i << ' ' << ml_ri << ' ' << ml_gr << ' ' << lnL << '\n';
	}
}
catch(EAny &e)
{
	e.print();
}
}
