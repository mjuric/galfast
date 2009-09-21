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

#include "skygen.h"
#include "analysis.h"
#include "sph_polygon.h"
#include "spline.h"
#include "pipeline.h"

#include <fstream>
#include <iomanip>

#include <dlfcn.h>

#include <astro/io/format.h>
#include <astro/system/config.h>
#include <astro/useall.h>

// Clips out stars not within the requested observation area
// Currently only used as support to os_skygen
class partitioned_skymap;
class os_clipper : public osink
{
protected:
	// hemisphere -- pixelized north/south sky (aux class).
	struct hemisphere
	{
		peyton::math::lambert 	proj;	// projection
		partitioned_skymap 	*sky;	// skymap of the hemisphere polygon

		hemisphere() : sky(NULL) {}
		~hemisphere() { delete sky; }
	private:
		// disallow copy, copy constructable
		hemisphere &operator=(const hemisphere&);
		hemisphere(const hemisphere&);
	};

	hemisphere hemispheres[2];	// pixelized northern and southern sky

public:
	struct pixel
	{
		double l, b; // in Radians
		double X, Y; // lambert coordinates
		int projIdx; // index of the projection for (X,Y)<->(l,b) transform (hemispheres[projIdx].proj)
		float pixelArea, coveredArea;	// area of the nominal pixel, area covered by the footprint within the pixel
		
		pixel(double l_, double b_, double X_, double Y_, int projIdx_, float pixelArea_, float coveredArea_)
			: X(X_), Y(Y_), l(l_), b(b_), projIdx(projIdx_), pixelArea(pixelArea_), coveredArea(coveredArea_) {}
	};

public:
	virtual size_t process(otable &in, size_t begin, size_t end, rng_t &rng);
	virtual bool construct(const peyton::system::Config &cfg, otable &t, opipeline &pipe); // NOTE: overriden as it's abstract, but should NEVER be called directly. Use construct_from_hemispheres() instead.

	virtual const std::string &name() const { static std::string s("clipper"); return s; }
	virtual int priority() { return PRIORITY_INSTRUMENT; } // ensure this module is placed near the end of the pipeline

	// constructs the clipper object from north/south hemispheres in projection proj
	void construct_from_hemispheres(float dx, const peyton::math::lambert &nproj, const std::pair<gpc_polygon, gpc_polygon> &sky);

	int getPixelCenters(std::vector<os_clipper::pixel> &pix) const;			// returns the centers of all pixels
	int getProjections(std::vector<std::pair<double, double> > &ppoles) const;	// returns the poles of all used projections

	os_clipper() : osink()
	{
		req.insert("projIdx");
		req.insert("projXY");
		prov.insert("hidden");
	}
};
extern "C" opipeline_stage *create_module_clipper() { return new os_clipper(); }	// Factory; called by opipeline_stage::create()

// GPU generator input
class os_skygen : public osource
{
protected:
	std::vector<boost::shared_ptr<skyConfigInterface> > kernels;
	size_t nstarLimit;	// maximum number of stars to generate
	float nstars;		// the mean number of stars to generate (if nstars=0, the number will be determined by the model)
	cuxTexture<float, 3>	ext_north, ext_south;	// north/south extinction maps

public:
	virtual bool construct(const peyton::system::Config &cfg, otable &t, opipeline &pipe);
	virtual size_t run(otable &t, rng_t &rng);
	virtual const std::string &name() const { static std::string s("skygen"); return s; }
	virtual const std::string &type() const { static std::string s("input"); return s; }

protected:
	const os_clipper &load_footprints(const std::string &footprints, float dx, opipeline &pipe);
	skyConfigInterface *create_kernel_for_model(const std::string &model);
	int load_models(otable &t, skygenConfig &sc, const std::string &model_cfg_list, const os_clipper &clipper);
	void load_pdf(float &dx, skygenConfig &sc, otable &t, const std::string &cfgfn);
	void load_extinction_maps(const std::string &econf);
};
extern "C" opipeline_stage *create_module_skygen() { return new os_skygen(); }	// Factory; called by opipeline_stage::create()

/***********************************************************************/

cuxTexture<float, 1> load_constant_texture(float val, float X0, float X1)
{
	cuxTexture<float, 1> tex(2);
	tex(0U) = val;
	tex(1U) = val;
	tex.coords[0].x = X0;
	tex.coords[0].y = 1./(X1 - X0);
	return tex;
}

cuxTexture<float, 1> construct_1D_texture_by_resampling(double *X, double *Y, int ndata, int nsamp)
{
	spline tx;
	tx.construct(X, Y, ndata);

	// resample to texture
	cuxTexture<float, 1> tex(nsamp);
	float X0 = X[0], X1 = X[ndata-1], dX = (X1 - X0) / (nsamp-1);
	for(int i=0; i != nsamp; i++)
	{
		tex(i) = tx(X0 + i*dX);
	}

	// construct 
	tex.coords[0].x = X0;
	tex.coords[0].y = 1./dX;
	return tex;
}

cuxTexture<float, 1> load_and_resample_1D_texture(const char *fn, int nsamp)
{
	// load the points from the file, and construct
	// a spline to resample from
	text_input_or_die(txin, fn);
	std::vector<double> X, Y;
	::load(txin, X, 0, Y, 1);

	return construct_1D_texture_by_resampling(&X[0], &Y[0], X.size(), nsamp);
}

//
// Return texcoord that will map x to imgx and y to imgy
//
float2 texcoord_from_range(float imgx, float imgy, float x, float y)
{
	float2 tc;

	tc.x = (-imgy * x + imgx * y)/(imgx - imgy);
	tc.y = (imgx - imgy) / (x - y);

	return tc;
}


int add_lonlat_rect(sph_polygon &poly, Radians l0, Radians b0, Radians l1, Radians b1, Radians dl, gpc_op op = GPC_UNION);
Radians clip_zone_of_avoidance(std::list<sph_polygon> &poly, int npoly, const Config &cfg)
{
	Radians bmin;
	cfg.get(bmin,	"bmin", 0.);
	bmin = rad(bmin);

	if(bmin)
	{
		MLOG(verb1) << "Zone of avoidance: clipping |b| < " << deg(bmin) << " deg";

		// add a "clip" contour to last npoly polygons
		REVEACH(poly)
		{
			if(!npoly--) { break; }

			const double dl = rad(.1);	// TODO: Make this user settable?
			add_lonlat_rect(*i, 0, -bmin, 0, bmin, dl, GPC_DIFF);
		}
	}

	return bmin;
}

int load_footprint_rect(std::list<sph_polygon> &foot, const peyton::system::Config &cfg);
int load_footprint_beam(std::list<sph_polygon> &foot, const peyton::system::Config &cfg);

int load_footprint(std::list<sph_polygon> &foot, const Config &cfg)
{
	if(!cfg.count("type")) { THROW(EAny, "The footprint configuration file must contain the 'type' keyword."); }

	int k = -1;
	     if(cfg["type"] == "beam") { k = load_footprint_beam(foot, cfg); }
	else if(cfg["type"] == "rect") { k = load_footprint_rect(foot, cfg); }

	if(k >= 0)
	{
		clip_zone_of_avoidance(foot, k, cfg);
		return k;
	}

	THROW(EAny, "The requested footprint type " + cfg["type"] + " is unknown.");
}

// Calculate the 'skymap', for fast point-in-polygon lookups of the observed
// sky footprint.
partitioned_skymap *make_skymap_slow(Radians dx, gpc_polygon sky)
{
	using boost::shared_ptr;
	using std::cout;

	partitioned_skymap *skymap = new partitioned_skymap;
	skymap->dx = dx;

 	poly_bounding_box(skymap->x0, skymap->x1, skymap->y0, skymap->y1, sky);
	std::cerr << "Creating fast skymap lookup (this may take a while) ... ";
	//std::cerr << skymap->x0 << " " <<  skymap->x1 << " " <<  skymap->y0 << " " <<  skymap->y1 << "\n";

	int X = 0;
	for(double x = skymap->x0; x < skymap->x1; x += skymap->dx, X++) // loop over all x values in the bounding rectangle
	{
		int Y = 0;
		double xa = x, xb = x+skymap->dx;
		for(double y = skymap->y0; y < skymap->y1; y += skymap->dx, Y++) // loop over all y values for a given x
		{
			double ya = y, yb = y+skymap->dx;
			gpc_polygon r = poly_rect(xa, xb, ya, yb);

			gpc_polygon poly;
			gpc_polygon_clip(GPC_INT, &sky, &r, &poly);
			if(poly.num_contours == 0) continue; // if there are no observations in this direction

			// store the polygon into a fast lookup map
			partitioned_skymap::pixel_t &pix = skymap->skymap[std::make_pair(X, Y)];
			pix.poly = poly;
			pix.coveredArea = polygon_area(poly);
			pix.pixelArea = sqr(skymap->dx);
		}
	}

	std::cerr << " done.\n";
	return skymap;
}

void make_skymap_piece(gpc_polygon sky, partitioned_skymap *skymap, int X0, int X1, int Y0, int Y1)
{
	// test for intersection
	double xa = skymap->x0 + X0*skymap->dx;
	double xb = skymap->x0 + X1*skymap->dx;
	double ya = skymap->y0 + Y0*skymap->dx;
	double yb = skymap->y0 + Y1*skymap->dx;
	gpc_polygon r = poly_rect(xa, xb, ya, yb);

	gpc_polygon poly;
	gpc_polygon_clip(GPC_INT, &sky, &r, &poly);
	if(poly.num_contours == 0)
	{
		return; // if there are no observations in this region
	}

	int DX = X1-X0, DY = Y1-Y0;
	if(DX == 1 && DY == 1) // leaf
	{
		partitioned_skymap::pixel_t &pix = skymap->skymap[std::make_pair(X0, Y0)];
		pix.poly = poly;
		pix.coveredArea = polygon_area(poly);
		pix.pixelArea = sqr(skymap->dx);
		return;
	}

	// subdivide further
	make_skymap_piece(poly, skymap, X0 +    0, X0 + DX/2, Y0 +    0, Y0 + DY/2);
	make_skymap_piece(poly, skymap, X0 + DX/2, X0 + DX,   Y0 +    0, Y0 + DY/2);
	make_skymap_piece(poly, skymap, X0 + DX/2, X0 + DX,   Y0 + DY/2, Y0 +   DY);
	make_skymap_piece(poly, skymap, X0 +    0, X0 + DX/2, Y0 + DY/2, Y0 +   DY);
	gpc_free_polygon(&poly);
}

// Calculate the 'skymap', for fast point-in-polygon lookups of the observed
// sky footprint.
partitioned_skymap *make_skymap(Radians dx, gpc_polygon sky)
{
	using boost::shared_ptr;
	using std::cout;

	partitioned_skymap *skymap = new partitioned_skymap;
	skymap->dx = dx;

 	poly_bounding_box(skymap->x0, skymap->x1, skymap->y0, skymap->y1, sky);

	int NX = (int)ceil((skymap->x1 - skymap->x0) / skymap->dx);
	int NY = (int)ceil((skymap->y1 - skymap->y0) / skymap->dx);
	// round up to nearest power-of-two
	NX = 1 << (int)ceil(log2(NX));
	NY = 1 << (int)ceil(log2(NY));
	if(NX == 1) NX++;
	if(NY == 1) NY++;
	NX = NY = std::max(NX, NY);	// make things really simple...
	assert(skymap->x0 + NX*skymap->dx >= skymap->x1);
	assert(skymap->y0 + NY*skymap->dx >= skymap->y1);

	// hierarchically subdivide
	make_skymap_piece(sky, skymap,    0, NX/2,    0, NY/2);
	make_skymap_piece(sky, skymap, NX/2, NX,      0, NY/2);
	make_skymap_piece(sky, skymap, NX/2, NX,   NY/2,   NY);
	make_skymap_piece(sky, skymap,    0, NX/2, NY/2,   NY);

	return skymap;
}

////////////////////////////////////////////////////////////////////////////

std::pair<gpc_polygon, gpc_polygon> project_to_hemispheres(const std::list<sph_polygon> &foot, const peyton::math::lambert &proj, Radians dx);

const os_clipper &os_skygen::load_footprints(const std::string &footprints, float dx, opipeline &pipe)
{
	// Load footprints. The result is the list of spherical polygons on the sky;
	std::list<sph_polygon> foot;

	std::vector<std::string> footstr;
	split(footstr, footprints);
	FOREACH(footstr)
	{
		Config cfg(*i);									// load footprint config
		load_footprint(foot, cfg);
	}

	// Project the footprint onto north/south hemispheres
	peyton::math::lambert proj(rad(90), rad(90));
	std::pair<gpc_polygon, gpc_polygon> sky = project_to_hemispheres(foot, proj, dx);

	// setup clipper for the footprint
	boost::shared_ptr<opipeline_stage> clipper_s(opipeline_stage::create("clipper"));	// clipper for this footprint
	os_clipper &clipper = *static_cast<os_clipper*>(clipper_s.get());
	clipper.construct_from_hemispheres(dx, proj, sky);
	pipe.add(clipper_s);

	gpc_free_polygon(&sky.first);
	gpc_free_polygon(&sky.second);

	return clipper;
}

typedef skyConfigInterface *(*modelFactory_t)();

skyConfigInterface *os_skygen::create_kernel_for_model(const std::string &model)
{
	void *me = dlopen(NULL, RTLD_LAZY);
	if(me == NULL)
	{
		const char *err = dlerror();
		THROW(EAny, err);
	}

	std::string factory_name = "create_model_";
	FOREACH(model)
	{
		if(!isalnum(*i)) { continue; }
		factory_name += tolower(*i);
	}

	modelFactory_t factory = (modelFactory_t)dlsym(me, factory_name.c_str());
	if(factory)
	{
		skyConfigInterface *ret = factory();
		if(ret) { return ret; }
	}

	THROW(EAny, "Unknow density model '" + model + "'");
	return NULL;
}

int os_skygen::load_models(otable &t, skygenConfig &sc, const std::string &model_cfg_list, const os_clipper &clipper)
{
	// projections the pixels are bound to
	std::vector<std::pair<double, double> > lb0;	// projection poles
	int nproj = clipper.getProjections(lb0);
	assert(nproj == 2);
	for(int i=0; i != lb0.size(); i++)
	{
		sc.proj[i].init(lb0[i].first, lb0[i].second);
	}

	// prepare the skypixels to be processed by subsequently loaded models
	std::vector<os_clipper::pixel> pix;
	sc.npixels = clipper.getPixelCenters(pix);
	std::vector<skypixel> skypixels(sc.npixels);
	FOR(0, sc.npixels)
	{
		os_clipper::pixel &p = pix[i];
		float coveredFraction = p.coveredArea / p.pixelArea;
		skypixels[i] = skypixel(p.l, p.b, p.X, p.Y, p.projIdx, sqrt(p.pixelArea), coveredFraction);
	}

	// Load density model configurations, instantiate
	// and configure the correct skyConfig kernels
	std::vector<std::string> models;
	split(models, model_cfg_list);
	FOREACH(models)
	{
		Config cfg(*i);		// model configuration file
		boost::shared_ptr<skyConfigInterface> kernel(create_kernel_for_model(cfg.get("model")));

		// setup the sky pixels to be processed by this model
		kernel->init(t, cfg, sc, &skypixels[0]);

		kernels.push_back(kernel);
	}
}

void os_skygen::load_pdf(float &dx, skygenConfig &sc, otable &t, const std::string &cfgfn)
{
	// load PDF config
	Config cfg(cfgfn);

	// sky binning scale
	cfg.get(dx,	"dx", 1.f);
	dx = rad(dx);

	// Density binning parameters
	sc.M0 = cfg.get_any_of("M0", "ri0");
	sc.M1 = cfg.get_any_of("M1", "ri1");
	sc.dM = cfg.get_any_of("dM", "dri");
	sc.m0 = cfg.get_any_of("m0", "r0");
	sc.m1 = cfg.get_any_of("m1", "r1");
	sc.dm = cfg.get("dm");
	assert(sc.dM == sc.dm);

	sc.nM = (int)round((sc.M1 - sc.M0) / sc.dM);
	sc.nm = (int)round((sc.m1 - sc.m0) / sc.dm);
	sc.m0 += 0.5*sc.dm;
	sc.M1 -= 0.5*sc.dM;

	// prepare the table for output
	t.use_column("lb");
	t.use_column("projIdx");
	t.use_column("projXY");
	t.use_column("Am");
	t.use_column("XYZ");
	t.use_column("comp");

	// prepare the PDF-related table columns for output
	std::string 	absmagName, bandName;
	bandName = cfg.get("band");
	absmagName = "abs" + bandName;

	t.use_column(absmagName); t.alias_column(absmagName, "absmag");
				  t.alias_column(absmagName, "M1");
				  t.set_column_property("absmag", "band", bandName);
	t.use_column("DM");
}

#if HAVE_LIBCFITSIO
#include "fitsio2.h"

void FITS_ERRCHECK(const std::string &h, int s)
{
	if((s) != 0)
	{
		char errmsg[160];
		std::ostringstream ss;
		ss << h << "\n";
		ss << "CFITSIO: ";
		bool first = true;
		while(fits_read_errmsg(errmsg))
		{
			if(!first) { ss << "         "; }
			ss << errmsg << "\n";
			first = false;
		}
		THROW(EAny, ss.str());
	}
}

float read_fits_key(fitsfile *fptr, const std::string &key, const std::string &fn, int *status_out = NULL)
{
	float ret;
	int status = 0;
	fits_read_key(fptr, TFLOAT, (char*)key.c_str(), &ret, NULL, &status);
	if(status_out)
		*status_out = status;
	else
		FITS_ERRCHECK("Error reading " + key + " key.", status);
	return ret;
}

float2 texcoord_from_wcs(fitsfile *fptr, int n, const std::string &fn, int *status_out = NULL)
{
	// read the mapping from pixels to coordinates from FITS header
	// This is stored in old-style CRPIXn CRVALn CDELTn keys (see
	// http://www.aanda.org/index.php?option=article&access=bibcode&bibcode=2002A%2526A...395.1061GPDF
	// for details)
	float imgx, x, dx, dummy;

	imgx = read_fits_key(fptr, "CRPIX" + str(n), fn, status_out) - 1; // -1 because of the FITS vs. C indexing convention
	if(status_out && *status_out) { float2 tc; return tc; }

	   x = read_fits_key(fptr, "CRVAL" + str(n), fn);
	  dx = read_fits_key(fptr, "CDELT" + str(n), fn);

	int status = 0;
	fits_read_key(fptr, TFLOAT, (char*)("CROTA" + str(n)).c_str(), &dummy, NULL, &status);
	if(status && status != KEY_NO_EXIST)
	{
		THROW(EAny, "Axis rotations not supported (CROTA" + str(n) + " keyword present in " + fn + ")");
	}

	float2 tc = {x - imgx * dx, 1.f/dx};
	return tc;
}

/*
	Generating a test extinction cube in IDL:
	=========================================

	arr = findgen(20, 30, 40)
	dim = size(arr)

	mkhdr, hdr, arr

	sxaddpar, hdr, 'CRPIX1',  1, 'Pixel coord'
	sxaddpar, hdr, 'CRVAL1', -2, 'Lambert coord'
	sxaddpar, hdr, 'CDELT1', (4.0/(dim[1]-1)), 'Coord increment'

	sxaddpar, hdr, 'CRPIX2',  1, 'Pixel coord'
	sxaddpar, hdr, 'CRVAL2', -2, 'Lambert coord'
	sxaddpar, hdr, 'CDELT2', (4.0/(dim[2]-1)), 'Coord increment'

	sxaddpar, hdr, 'CRPIX3',  1, 'Pixel coord'
	sxaddpar, hdr, 'CRVAL3',  0, 'Distance modulus coord'
	sxaddpar, hdr, 'CDELT3', (30.0/(dim[3]-1)), 'Coord increment'

	writefits, 'myfile.fits', arr, hdr
*/

//
// (l0,phi1) will be loaded from (LAMBDA0,PHI1) keywords. If any of the two are not present,
//           both shall be set to -100 on return from the subroutine.
//
cuxTexture<float, 3> load_extinction_map(const std::string &fn, Radians &l0, Radians &phi1)
{
	fitsfile *fptr;
	int status = 0;
	fits_open_file(&fptr, fn.c_str(), READONLY, &status);
	FITS_ERRCHECK("Error opening extinction map.", status);

	int bitpix, naxis;
	long naxes[3];
	fits_get_img_param(fptr, 3, &bitpix, &naxis, naxes, &status);
	FITS_ERRCHECK("Error reading extinction map header.", status);
	if(naxis != 3) { THROW(EAny, "Supplied extinction map file is " + str(naxis) + ", instead of 3-dimensional."); }
	if(bitpix != FLOAT_IMG) { THROW(EAny, "Supplied extinction map file is not a data cube of 32-bit floats"); }

	cuxSmartPtr<float> img(naxes[0], naxes[1], naxes[2], -1, 1); // allocate a 1-byte aligned floating point 3D cube
	uint32_t nelem = img.size();
	long fpixel[3] = { 1, 1, 1 };
	fits_read_pix(fptr, TFLOAT, fpixel, nelem, NULL, &img(0, 0, 0), NULL, &status);
	FITS_ERRCHECK("Error reading extinction map.", status);

	// setup texture coordinates
	float2 tc[3];
	tc[0] = texcoord_from_wcs(fptr, 1, fn, &status);
	if(status == 0)
	{
		tc[1] = texcoord_from_wcs(fptr, 2, fn);
		tc[2] = texcoord_from_wcs(fptr, 3, fn);
	}
	else
	{
		MLOG(verb1) << "WARNING: Did not find coordinate system keys in " << fn << ", will be using the defaults.";
		tc[0] = texcoord_from_range(0, img.extent(0)-1, -2,  2);	// X range
		tc[1] = texcoord_from_range(0, img.extent(1)-1, -2,  2);	// Y range
		tc[2] = texcoord_from_range(0, img.extent(2)-1, 10, 15);	// DM range
	}

	// load lambert projection pole (if specified)
	l0   = rad(read_fits_key(fptr, "LAMBDA0", fn, &status));
	phi1 = rad(read_fits_key(fptr, "PHI1",    fn, &status));
	if(status) { l0 = phi1 = -100; status = 0; }

	fits_close_file(fptr, &status);
	FITS_ERRCHECK("Error closing extinction map file.", status);

	return cuxTexture<float, 3>(img, tc);
}
#endif

cuxSmartPtr<float4> resample_extinction_texture(cuxTexture<float, 3> &tex, float2 crange[3], int npix[3], ::lambert *proj);
void resample_and_output_texture(const std::string &outfn, cuxTexture<float, 3> &tex, float2 crange[3], int npix[3], ::lambert *proj)
{
	cuxSmartPtr<float4> res = resample_extinction_texture(tex, crange, npix, proj);

	std::cerr << "Writing output to " << outfn << "\n";
	std::ofstream out(outfn.c_str());
	FOREACH(res)
	{
		float4 v = *i;
		out << io::format("%12.8f %12.8f %12.8f   %8.5f\n") << v.x << v.y << v.z << v.w;
	}
}

void resample_texture(const std::string &outfn, const std::string &texfn, float2 crange[3], int npix[3], bool deproject, Radians l0req, Radians b0req)
{
	Radians l0, b0;
	cuxTexture<float, 3> tex = load_extinction_map(texfn, l0, b0);

	// compute input texture ranges
	// autodetect crange and npix if not given
	float2 irange[3];
	for(int i = 0; i != 3; i++)
	{
		irange[i].x = tex.coords[i].x;
		irange[i].y = tex.coords[i].x + (tex.extent(i)-1) / tex.coords[i].y;

		if(npix[i] == 0) { npix[i] = tex.extent(i); }

		// convert to radians, if deprojecting
		if(deproject && i != 2)
		{
			crange[i].x = rad(crange[i].x);
			crange[i].y = rad(crange[i].y);
		}

		if(crange[i].x == crange[i].y)
		{
			if(!deproject || i == 2)
			{
				crange[i] = irange[i];
			}
			else
			{
				// deprojecting. i=0 is longitude, i=1 is latitude
				switch(i)
				{
				case 0:
					crange[0] = make_float2(0, ctn::twopi);
					break;
				case 1:
					crange[1] = make_float2(-ctn::halfpi, ctn::halfpi);
					break;
				}
			}
		}
	}

	::lambert proj;

	MLOG(verb1) << " Input texture : x = [" << irange[0].x << ", " << irange[0].y << "] (" << tex.extent(0) << " pixels)\n";
	MLOG(verb1) << "               : y = [" << irange[1].x << ", " << irange[1].y << "] (" << tex.extent(1) << " pixels)\n";
	MLOG(verb1) << "               : z = [" << irange[2].x << ", " << irange[2].y << "] (" << tex.extent(2) << " pixels).\n";
	MLOG(verb1) << "     Deproject : " << (deproject ? "yes" : "no");

	if(deproject)
	{
		std::string pole;
		if(l0 == -100.)
		{
			if(l0req == -100.)
			{
				l0 = rad(90);
				b0 = rad(90);
				pole = "(default)";
			}
			else
			{
				l0 = l0req;
				b0 = b0req;
				pole = "(manually set)";
			}
		}
		else
		{
			pole = "(read from file header)";
		}

		proj.init(l0, b0);

		MLOG(verb1) << "     Proj. pole : l0=" << deg(proj.l0) << ", b0=" << deg(proj.b0) << " " << pole;
		MLOG(verb1) << "Output texture : l = [" << deg(crange[0].x) << ", " << deg(crange[0].y) << "] (" << npix[0] << " pixels)\n";
		MLOG(verb1) << "               : b = [" << deg(crange[1].x) << ", " << deg(crange[1].y) << "] (" << npix[1] << " pixels)\n";
	}
	else
	{
		MLOG(verb1) << "Output texture : x = [" << crange[0].x << ", " << crange[0].y << "] (" << npix[0] << " pixels)\n";
		MLOG(verb1) << "               : y = [" << crange[1].x << ", " << crange[1].y << "] (" << npix[1] << " pixels)\n";
	}

	MLOG(verb1) << "               : z = [" << crange[2].x << ", " << crange[2].y << "] (" << npix[2] << " pixels).\n";

	resample_and_output_texture(outfn, tex, crange, npix, deproject ? &proj : NULL);

	MLOG(verb2) << "Resampled.";
}

bool load_extinction_maps(cuxTexture<float, 3> &ext_north, cuxTexture<float, 3> &ext_south, const std::string &econf)
{
	if(!econf.size())	// no extinction
	{
		cuxSmartPtr<float> tex(2, 2, 2);
		FORj(i, 0, 2) FORj(j, 0, 2) FORj(k, 0, 2)
			tex(i, j, k) = 0.f;
		float2 tc = make_float2(-2, 1./4.);

		ext_south = ext_north = cuxTexture<float, 3>(tex, tc);

		MLOG(verb2) << "Extinction maps not given, assuming no extinction.";

		return false;
	}

	// load 3D FITS file
#if HAVE_LIBCFITSIO
	// econf format: <ext_north.fits> <ext_south.fits> [scaling_factor]
	std::string northfn, southfn;
	std::istringstream ss(econf);
	if(!(ss >> northfn >> southfn))
	{
		THROW(EAny, "'<ext_north.fits> <ext_south.fits> [scaling_factor]' format expected for extinction map specification (keyword 'econf')");
	}
	float scale = 1.f;
	ss >> scale;

	Radians l0, b0;
	ext_north = load_extinction_map(northfn, l0, b0);
	if(
		(l0 != -100. && fabs(l0 - ctn::halfpi) > 1e-5) ||
		(b0 != -100. && fabs(b0 - ctn::halfpi) > 1e-5)
	)
	{
		MLOG(verb1) << "WARNING: Expecting l=90, b=90 as pole of northern hemisphere projection. Got " << deg(l0) << " " << deg(b0) << " instead.";
	}
	ext_south = load_extinction_map(southfn, l0, b0);
	if(
		(l0 != -100. && fabs(l0 + ctn::halfpi) > 1e-5) ||
		(b0 != -100. && fabs(b0 + ctn::halfpi) > 1e-5)
	)
	{
		MLOG(verb1) << "WARNING: Expecting l=-90, b=-90 as pole of southern hemisphere projection. Got " << deg(l0) << " " << deg(b0) << " instead.";
	}

	FOREACH(ext_north) { *i *= scale; }
	FOREACH(ext_south) { *i *= scale; }

	MLOG(verb1) << "Extinction maps: " << northfn << " (north), " << southfn << " (south).";
	MLOG(verb2) << "Extinction north: " << northfn << " [ X x Y x DM = " << ext_north.width() << " x " << ext_north.height() << " x " << ext_north.depth() << "]\n";
	MLOG(verb2) << "Extinction south: " << southfn << " [ X x Y x DM = " << ext_south.width() << " x " << ext_south.height() << " x " << ext_south.depth() << "]\n";
	MLOG(verb2) << "Ext. scale factor: " << scale << "\n";

#if 0	// debug -- resample the north sky into a text file
	float2 crange[3];
	int npix[3];
	crange[0] = make_float2(-2, 2);
	crange[1] = make_float2(-2, 2);
	crange[2] = make_float2(0, 30);
	npix[0] = 20;
	npix[1] = 30;
	npix[2] = 40;
	resample_and_output_texture("northAr.txt", ext_north, tc_ext_north, crange, npix, NULL);
	abort();
#endif

#else // HAVE_LIBCFITSIO
	THROW(EAny, "Cannot load extinction maps from FITS files; recompile with FITS I/O support.");
#endif // HAVE_LIBCFITSIO

	return true;
}

void os_skygen::load_extinction_maps(const std::string &econf)
{
	::load_extinction_maps(ext_north, ext_south, econf);
}

bool os_skygen::construct(const Config &cfg, otable &t, opipeline &pipe)
{
	// maximum number of stars skygen is allowed to generate (0 for unlimited)
	cfg.get(nstarLimit, "nstarlimit", (size_t)100*1000*1000);
	// mean number of stars skygen should generate (0 to leave it to the model to determine this)
	cfg.get(nstars, "nstars", 0.f);

	// load extinction volume maps
	load_extinction_maps(cfg["extinction"]);

	// load PDF/sky pixelization config and prepare the table for output
	skygenConfig sc;
	float dx;
	load_pdf(dx, sc, t, cfg.get("pdf"));

	// load footprints and construct the clipper
	const os_clipper &clipper = load_footprints(cfg.get("foot"), dx, pipe);

	// load models
	load_models(t, sc, cfg.get("model"), clipper);

	return true;
}

DECLARE_TEXTURE(ext_north, float, 3, cudaReadModeElementType);
DECLARE_TEXTURE(ext_south, float, 3, cudaReadModeElementType);

size_t os_skygen::run(otable &in, rng_t &rng)
{
	cuxTextureBinder tb_north(::ext_north, ext_north);
	cuxTextureBinder tb_south(::ext_south, ext_south);

	double nstarsExpected = 0;
	FOREACH(kernels)
	{
		(*i)->initRNG(rng);
		nstarsExpected += (*i)->integrateCounts();
	}

	//
	// Adjust normalization, if nstars is given
	//
	if(nstars != 0)
	{
		double norm = nstars / nstarsExpected;
		FOREACH(kernels) { (*i)->setDensityNorm(norm); }
		nstarsExpected = nstars;

		MLOG(verb1) << "Stars expected: " << nstarsExpected << " (forced)\n";
	}
	else
	{
		MLOG(verb1) << "Stars expected: " << nstarsExpected << "\n";
	}

	//
	// Stop if the expected number of stars is more than the limit
	//
	if(nstarsExpected > nstarLimit && nstarLimit != 0)
	{
		THROW(EAny, "The expected number of generated stars (" + str(nstarsExpected) + ") "
			"exceeds the configuration limit (" + str(nstarLimit) + "). Increase "
			"nstarlimit or set to zero for unlimited."
			);
	}

	size_t starsGenerated = 0;
	FOREACH(kernels)
	{
		starsGenerated += (*i)->run(in, nextlink);
	}
	return starsGenerated;
}

////////////////////////////////////////////////////////////////////////////
//
//	os_clipper -- clip the output to observed sky footprint
//
////////////////////////////////////////////////////////////////////////////

bool os_clipper::construct(const Config &cfg, otable &t, opipeline &pipe)
{
	THROW(EAny, "Module 'clipper' must never be instantiated directly.");

	return true;
}

// Construct a pixelization of the sky given footprint polygons projected
// to north and south Galactic hempsiphere
void os_clipper::construct_from_hemispheres(float dx, const peyton::math::lambert &proj, const std::pair<gpc_polygon, gpc_polygon> &sky)
{
	// set the north hemisphere projection to input map projection
	hemispheres[0].proj = proj;
	hemispheres[1].proj = peyton::math::lambert(modulo(proj.l0 + ctn::pi, ctn::twopi), -proj.phi1);

	// construct north/south skymaps
	hemispheres[0].sky = make_skymap(dx, sky.first);
	hemispheres[1].sky = make_skymap(dx, sky.second);

	DLOG(verb1) << "Sky pixels in the north: " << hemispheres[0].sky->skymap.size();
	DLOG(verb1) << "Sky pixels in the south: " << hemispheres[1].sky->skymap.size();
}

// returns the poles of all used projections
int os_clipper::getProjections(std::vector<std::pair<double, double> > &ppoles) const
{
	ppoles.clear();
	for(int i = 0; i != 2; i++)
	{
		const peyton::math::lambert &proj = hemispheres[i].proj;
		ppoles.push_back(std::make_pair(proj.l0, proj.phi1));
	}
	return ppoles.size();
}

// returns the centers of all sky pixels (pencil beams into which
// the sky has been pixelized by construct_from_hemispheres())
int os_clipper::getPixelCenters(std::vector<os_clipper::pixel> &pix) const
{
	pix.clear();

	FORj(projIdx, 0, 2)
	{
		partitioned_skymap *skymap = hemispheres[projIdx].sky;
		FOREACH(skymap->skymap)
		{
			Radians x, y, l, b;
			x = skymap->x0 + skymap->dx*(i->first.first  + 0.5);
			y = skymap->y0 + skymap->dx*(i->first.second + 0.5);
			hemispheres[projIdx].proj.deproject(l, b, x, y);
			pix.push_back(pixel(l, b, x, y, projIdx, i->second.pixelArea, i->second.coveredArea));
		}
	}

	return pix.size();
}

// ::process() override -- set hidden=1 for every row that is outside
// the exact input sky footprint
size_t os_clipper::process(otable &in, size_t begin, size_t end, rng_t &rng)
{
	// fetch prerequisites
	cint_t::host_t pIdx     = in.col<int>("projIdx");
	cfloat_t::host_t projXY = in.col<float>("projXY");
	cint_t::host_t	hidden  = in.col<int>("hidden");

	// debugging statistics
	int nstars[2] = { 0, 0 };
	for(size_t row=begin; row < end; row++)
	{
		// clip everything outside the footprint polygon
		int projIdx = pIdx(row);
		nstars[projIdx]++;

		Radians x, y;
 		x = projXY(row, 0);
 		y = projXY(row, 1);

		// immediately reject if in the southern hemisphere (for this projection)
		if(sqr(x) + sqr(y) > 2.)
		{
			hidden(row) = 1;
			continue;
		}

		partitioned_skymap *skymap = hemispheres[projIdx].sky;
		std::pair<int, int> XY;
		XY.first = (int)((x - skymap->x0) / skymap->dx);
		XY.second = (int)((y - skymap->y0) / skymap->dx);

		typeof(skymap->skymap.begin()) it = skymap->skymap.find(XY);
		if(it == skymap->skymap.end()) { continue; }

		// check that the star is inside survey footprint, reject if it's not
		gpc_polygon &poly = it->second.poly;
		gpc_vertex vtmp = { x, y };

		bool infoot = in_polygon(vtmp, poly);
		hidden(row) = !infoot;
	}

	DLOG(verb1) << "nstars north: " << nstars[0];
	DLOG(verb1) << "nstars south: " << nstars[1];

	return nextlink->process(in, begin, end, rng);
}
