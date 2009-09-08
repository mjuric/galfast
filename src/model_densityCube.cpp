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

#include "analysis.h"
#include <fstream>

#include "model_densityCube.h"
#include "gdc.h"
#include "skyconfig_impl.h"

#include <astro/system/config.h>

void densityCube::prerun(host_state_t &hstate, bool draw)
{
	// bind the luminosity function texture to texture reference
	densityCubeLF.bind(hstate.lf);
	densityCubeTex.bind(hstate.den);
}

void densityCube::postrun(host_state_t &hstate, bool draw)
{
	// unbind LF texture reference
	densityCubeLF.unbind();
	densityCubeTex.unbind();
}

int get_delta(double &delta, double &xmin, double &xmax, const std::vector<double> &x)
{
	std::vector<double> tmp(x);
	std::sort(tmp.begin(), tmp.end());
	tmp.erase(std::unique(tmp.begin(), tmp.end()), tmp.end());
	
	if(tmp.size() == 0)
	{
		delta = 1.;
	}
	else
	{
		std::vector<double> del;
		del.reserve(tmp.size());
		for(int i = 1; i != tmp.size(); i++)
		{
			del[i-1] = tmp[i] - tmp[i-1];
		}
		std::sort(del.begin(), del.end());
		delta = del[0];
	}
	xmin = tmp.front();
	xmax = tmp.back();

	double fn = ((xmax - xmin) / delta) + 1;
	int n = (int)rint(fn);

#if 0
	// sanity check (works only on grids where there are no gaps)
	int ns = tmp.size();
	assert(n == ns);
#endif

	return n;
}

cuxTexture<float, 3> load_3D_texture_from_text(const std::string &denfn, float offs[3], float scale[3])
{
	cuxTexture<float, 3> den;

	text_input_or_die(in, denfn);
	std::vector<double> x, y, z, val;
	load(in, x, 0, y, 1, z, 2, val, 3);
	for(int i = 0; i != x.size(); i++)
	{
		x[i] = x[i]*scale[0] + offs[0];
		y[i] = y[i]*scale[1] + offs[1];
		z[i] = z[i]*scale[2] + offs[2];
	}

	// auto-detect dx,dy,dz and the extent and add a zero-density border
	// around the cube
	double dx, xmin, xmax;
	double dy, ymin, ymax;
	double dz, zmin, zmax;
	int nx = get_delta(dx, xmin, xmax, x)+2;
	int ny = get_delta(dy, ymin, ymax, y)+2;
	int nz = get_delta(dz, zmin, zmax, z)+2;
	den.coords[0] = texcoord_from_range(0+1, nx-1-1, xmin, xmax);
	den.coords[1] = texcoord_from_range(0+1, ny-1-1, ymin, ymax);
	den.coords[2] = texcoord_from_range(0+1, nz-1-1, zmin, zmax);

	den.tex = cuxSmartPtr<float>(nx, ny, nz);
	FOREACH(den.tex) { *i = 0.f; }
	FORj(at, 0, x.size())
	{
		double xx = x[at], yy = y[at], zz = z[at];
		int i = 1 + (int)rint((xx - xmin) / dx);
		int j = 1 + (int)rint((yy - ymin) / dy);
		int k = 1 + (int)rint((zz - zmin) / dz);
		double v = val[at];

		den.tex(i, j, k) = v;
	}

#if 0
	FORj(k, 0, den.tex.depth())
	{
		FORj(j, 0, den.tex.height())
		{
			FORj(i, 0, den.tex.width())
			{
				std::cerr << den.tex(i, j, k) << " ";
			}
			std::cerr << "\n";
		}
		std::cerr << "\n";
	}
#endif

	return den;
}


cuxTexture<float, 3> load_3D_texture_from_gdc(const std::string &denfn)
{
	FILE* densityFile = fopen(denfn.c_str(), "r");
	if(densityFile==NULL)
	{
		THROW(EAny, "Error: can not open density file " + denfn);
	}

	gdc::GdcFileRead gfr;
	gfr.in = densityFile;
	int status = gfr.ReadGdcHeader();
	if(status != 0)
	{
		fclose(densityFile);
		THROW(EAny, "Error reading GDC file header. Bytes read: " + str(gfr.bytesRead));
	}

	// Prepare the texture ranges
	int 	nx = gfr.par.xBinN,
		ny = gfr.par.yBinN,
		nz = gfr.par.zBinN;
	float	dx = gfr.par.kpcBin*1000;
	float	xmax = nx * dx / 2.,
		ymax = ny * dx / 2.,
		zmax = nz * dx / 2.,
		xmin = -xmax,
		ymin = -ymax,
		zmin = -zmax;

	printf("GDC header data:\n");
	printf("    xSize: %d   ySize: %d   zSize: %d\n", gfr.par.xBinN, gfr.par.yBinN, gfr.par.zBinN);
	printf("    binKpc: %.2f    maxError %.4f\n",gfr.par.kpcBin, gfr.par.stepPerc/2 ); 
	printf("    bounds x (kpc): %.2f %.2f\n", xmin/1000., xmax/1000. );
	printf("    bounds y (kpc): %.2f %.2f\n", ymin/1000., ymax/1000. );
	printf("    bounds z (kpc): %.2f %.2f\n", zmin/1000., zmax/1000. );

	cuxTexture<float, 3> den(cuxSmartPtr<float>(nx+2, ny+2, nz+2));
	FOREACH(den.tex) { *i = 0.f; }

#if 1
	MLOG(verb1) << "Loading density (GDC file)...";
	for(int iz=0; iz < nz; iz++)
	{
		#if VERBOSE_OUTPUT
			printf("iz: %d  ",iz);
			fflush(stdout);
			gfr.ResetStatistics();
		#endif
		for (int ix=0; ix < nx; ix++)
		{
			for (int iy=0; iy < ny; iy++)
			{
				float v = gfr.ReadGdcData();
				den.tex(ix+1, iy+1, iz+1) = v;
			}
		}
		#if VERBOSE_OUTPUT
			gfr.PrintStatistics();
		#endif
	}
	fclose(densityFile);

	MLOG(verb1) << "Total density in density file: " << std::setprecision(4) << gfr.totalDensity ;

	//divide by volume
	float dV = cube(dx);
	FOREACH(den.tex) { *i /= dV; }
#endif

	// texture coordinates (pixel 1 corresponds to xmin, pixel nx to xmax, because of padding (see above))
	den.coords[0] = texcoord_from_range(1, nx, xmin, xmax);
	den.coords[1] = texcoord_from_range(1, ny, ymin, ymax);
	den.coords[2] = texcoord_from_range(1, nz, zmin, zmax);

	return den;
#if 0
	if (cfg.count("densityFactor")!=0 && cfg.count("totalDensity")!=0) {
		printf("Warning: in fileDensityModel: Both densityFactor and totalDensity given, totalDensity will take precedence\n");	
		}
	
	float densityFactor	=1.0f;
	if (cfg.count("totalDensity")!=0) {
		float totalDensityRequested=cfg.get("totalDensity");
		densityFactor=totalDensityRequested/gfr.totalDensity;
		}
	else if (cfg.count("densityFactor")!=0){
		densityFactor=cfg.get("densityFactor");
		}	
	MLOG(verb1) << "Density factor: " << std::setprecision(5) << densityFactor;	

	//divide by volume
	densityFactor/=parsecBin*parsecBin*parsecBin;
		
	for (int ix=0;ix<xBinN;ix++)
		for (int iy=0;iy<yBinN;iy++)
			for (int iz=0;iz<zBinN;iz++)
				GetBinDenHost(ix,iy,iz)*=densityFactor;
#endif
}

inline std::string extension(const std::string &fn)
{
	size_t dot = fn.rfind('.');
	if(dot == std::string::npos) { return ""; }
	return fn.substr(dot+1);
}

void densityCube::load(host_state_t &hstate, const peyton::system::Config &cfg)
{
	// density distribution parameters
	cfg.get(f,  	"f",		1.f);
	cfg.get(comp,	"comp",		0);

	// 3D density cube
	float offs[3], scale[3];
	cfg.get(offs[0], "xoffs", 0.f);
	cfg.get(offs[1], "yoffs", 0.f);
	cfg.get(offs[2], "zoffs", 0.f);
	cfg.get(scale[0], "xscale", 1.f);
	cfg.get(scale[1], "yscale", 1.f);
	cfg.get(scale[2], "zscale", 1.f);

	std::string denfn = cfg.get("density");
	std::string ext = tolower(extension(denfn));
	if(ext == "txt")
	{
		hstate.den = load_3D_texture_from_text(denfn, offs, scale);
	}
	else if(ext == "gdc")
	{
		hstate.den = load_3D_texture_from_gdc(denfn);
	}
	else
	{
		THROW(EAny, "Unknown extension '" + ext + "' for density model file");
	}

	// apply scaling to textures
//		float xi = (x - tcx.x) * tcx.y + 0.5f;
	FOR(0, 3)
	{
		hstate.den.coords[i].x  = (hstate.den.coords[i].x + offs[i]) * scale[i];
		hstate.den.coords[i].y /= scale[i];
	}

	// luminosity function
	if(cfg.count("lumfunc"))
	{
		hstate.lf = load_and_resample_1D_texture(cfg["lumfunc"].c_str());
	}
	else
	{
		hstate.lf = load_constant_texture(1.f);
	}
}

extern "C" skyConfigInterface *create_model_densitycube()
{
	return new skyConfig<densityCube>();
}
