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

#include "gpc_cpp.h"

#include <boost/shared_ptr.hpp>
#include <boost/lambda/lambda.hpp>
#include <sstream>

#include "simulate.h"
#include "projections.h"
#include "model.h"
#include "paralax.h"
#include "analysis.h"
//#include "dm.h"
#include "io.h"
#include "gpu.h"
#include "simulate_base.h"

#include <vector>
#include <map>
#include <string>

#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>

#include <astro/io/format.h>
#include <astro/math.h>
#include <astro/util.h>
#include <astro/system/log.h>
#include <astro/useall.h>

using namespace boost::lambda;
namespace ct = column_types;
#include "observe.h"

int split_fvec(std::vector<float>& arr, const std::string &text)
{
	std::vector<float>::value_type tmp;
	std::stringstream ss(text);

	arr.clear();
	while(ss >> tmp) { arr.push_back(tmp); }
	return arr.size();
}

std::vector<float> split_fvec(const std::string &text)
{
	std::vector<float> ret;
	split_fvec(ret, text);
	return ret;
}

void fvecToFarray(const std::vector<float>& src, farray5& dst)
{
	for (int i=0;i<5;i++) {
		dst[i]=src[i];
		}
}

void farray_to_iarray(farray5& fa, iarray5& ia)
{
	for (int i=0;i<5;i++)
		ia[i]=int(fa[i]*100);
}

void farray_to_i8array(farray5& fa, i8array5& ia)
{
	for (int i=0;i<5;i++)
		ia[i]=char(fa[i]/10);
}

extern os_kinTMIII_data os_kinTMIII_par;

DECLARE_KERNEL(
	os_kinTMIII_kernel(otable_ks ks, gpu_rng_t rng, ct::cint::gpu_t comp, ct::cfloat::gpu_t XYZ, ct::cfloat::gpu_t vcyl))
size_t os_kinTMIII::process(otable &in, size_t begin, size_t end, rng_t &rng)
{
	// ASSUMPTIONS:
	//	- Bahcall-Soneira component tags exist in input
	//	- galactocentric XYZ coordinates exist in input

	// fetch prerequisites
	using namespace column_types;
	cint   &comp  = in.col<int>("comp");
	cfloat &XYZ   = in.col<float>("XYZ");
	cfloat &vcyl   = in.col<float>("vcyl");

// 	os_kinTMIII_data_int par_int;
// 	
// 	farray_to_iarray(vR, par_int.vR);
//     farray_to_iarray(vPhi1, par_int.vPhi1);
//     farray_to_iarray(vZ, par_int.vZ);
//     farray_to_iarray(sigmaRR, par_int.sigmaRR);
//     farray_to_iarray(sigmaRPhi, par_int.sigmaRPhi);
//     farray_to_iarray(sigmaRZ, par_int.sigmaRZ);
//     farray_to_iarray(sigmaPhiPhi1, par_int.sigmaPhiPhi1);
//     farray_to_iarray(sigmaPhiPhi2, par_int.sigmaPhiPhi2);
//     farray_to_iarray(sigmaZPhi, par_int.sigmaZPhi);
//     farray_to_iarray(sigmaZZ, par_int.sigmaZZ);
// 
// 	farray_to_i8array(HvR, par_int.HvR);
// 	farray_to_i8array(HvPhi, par_int.HvPhi);
// 	farray_to_i8array(HvZ, par_int.HvZ);
// 	farray_to_i8array(HsigmaRR, par_int.HsigmaRR);
// 	farray_to_i8array(HsigmaRPhi, par_int.HsigmaRPhi);
// 	farray_to_i8array(HsigmaRZ, par_int.HsigmaRZ);
// 	farray_to_i8array(HsigmaPhiPhi, par_int.HsigmaPhiPhi);
// 	farray_to_i8array(HsigmaZPhi, par_int.HsigmaZPhi);
// 	farray_to_i8array(HsigmaZZ, par_int.HsigmaZZ);
// 
// 	par_int.fk=fk;
// 	par_int.DeltavPhi=DeltavPhi;

#if HAVE_CUDA
	{
		activeDevice dev(gpuExecutionEnabled("os_kinTMIII_kernel")? 0 : -1);
		if(gpuGetActiveDevice() >= 0)
		{
			cudaError err = cudaMemcpyToSymbol("os_kinTMIII_par", (os_kinTMIII_data*)this, sizeof(os_kinTMIII_data));
			if(err != cudaSuccess) { MLOG(verb1) << "CUDA Error: " << cudaGetErrorString(err); abort(); }
		}
	}
#endif
	os_kinTMIII_par = *this;

	CALL_KERNEL(os_kinTMIII_kernel, otable_ks(begin, end), rng, comp, XYZ, vcyl);
	return nextlink->process(in, begin, end, rng);
}


template<typename T> inline OSTREAM(const std::vector<T> &v) { FOREACH(v) { out << *i << " "; }; return out; }

bool os_kinTMIII::construct(const Config &cfg, otable &t, opipeline &pipe)
{
	cfg.get(fk           , "fk"           , 3.0f);
	cfg.get(DeltavPhi    , "DeltavPhi"    , 34.0f);
	fk = fk / (1. + fk);	// renormalize to probability of drawing from the first gaussian

	typedef std::vector<float> fvec;
		fvec 	vPhi1_fvec, vPhi2_fvec, vR_fvec, vZ_fvec,
			sigmaPhiPhi1_fvec, sigmaPhiPhi2_fvec, sigmaRR_fvec, sigmaZZ_fvec, sigmaRPhi_fvec, sigmaZPhi_fvec, sigmaRZ_fvec,
			HvPhi_fvec, HvR_fvec, HvZ_fvec,
			HsigmaPhiPhi_fvec, HsigmaRR_fvec, HsigmaZZ_fvec, HsigmaRPhi_fvec, HsigmaZPhi_fvec, HsigmaRZ_fvec;

	cfg.get(vR_fvec           , "vR"           , split_fvec("0 0 0 0 0"));	
	cfg.get(vPhi1_fvec        , "vPhi"         , split_fvec("-194 19.2 1.25 0 0"));	
	cfg.get(vZ_fvec           , "vZ"           , split_fvec("0 0 0 0 0"));	
	cfg.get(sigmaRR_fvec      , "sigmaRR"      , split_fvec("40 5 1.5 0 0"));	
	cfg.get(sigmaRPhi_fvec    , "sigmaRPhi"    , split_fvec("0 0 0 0 0"));	
	cfg.get(sigmaRZ_fvec      , "sigmaRZ"      , split_fvec("0 0 0 0 0"));	
	cfg.get(sigmaPhiPhi1_fvec , "sigmaPhiPhi1" , split_fvec("12 1.8 2 0 11"));		// dynamically changeable
	cfg.get(sigmaPhiPhi2_fvec , "sigmaPhiPhi2" , split_fvec("34 1.2 2 0 0"));
	cfg.get(sigmaZPhi_fvec    , "sigmaZPhi"    , split_fvec("0 0 0 0 0"));	
	cfg.get(sigmaZZ_fvec      , "sigmaZZ"      , split_fvec("25 4 1.5 0 0"));	

	fvecToFarray(vR_fvec, vR);
    fvecToFarray(vPhi1_fvec, vPhi1);
    fvecToFarray(vZ_fvec, vZ);
    fvecToFarray(sigmaRR_fvec, sigmaRR);
    fvecToFarray(sigmaRPhi_fvec, sigmaRPhi);
    fvecToFarray(sigmaRZ_fvec, sigmaRZ);
    fvecToFarray(sigmaPhiPhi1_fvec, sigmaPhiPhi1);
    fvecToFarray(sigmaPhiPhi2_fvec, sigmaPhiPhi2);
    fvecToFarray(sigmaZPhi_fvec, sigmaZPhi);
    fvecToFarray(sigmaZZ_fvec, sigmaZZ);	

	cfg.get(HvR_fvec          , "HvR"          , split_fvec("0 0 0 0 0"));	
	cfg.get(HvPhi_fvec        , "HvPhi"        , split_fvec("0 0 0 0 0"));	
	cfg.get(HvZ_fvec          , "HvZ"          , split_fvec("0 0 0 0 0"));	
	cfg.get(HsigmaRR_fvec     , "HsigmaRR"     , split_fvec("135 0 0 0 0"));	
	cfg.get(HsigmaRPhi_fvec   , "HsigmaRPhi"   , split_fvec("0 0 0 0 0"));	
	cfg.get(HsigmaRZ_fvec     , "HsigmaRZ"     , split_fvec("0 0 0 0 0"));	
	cfg.get(HsigmaPhiPhi_fvec , "HsigmaPhiPhi" , split_fvec("85 0 0 0 0"));	
	cfg.get(HsigmaZPhi_fvec   , "HsigmaZPhi"   , split_fvec("0 0 0 0 0"));	
	cfg.get(HsigmaZZ_fvec     , "HsigmaZZ"     , split_fvec("85 0 0 0 0"));	

	fvecToFarray(HvR_fvec, HvR);
	fvecToFarray(HvPhi_fvec, HvPhi);
	fvecToFarray(HvZ_fvec, HvZ);
	fvecToFarray(HsigmaRR_fvec, HsigmaRR);
	fvecToFarray(HsigmaRPhi_fvec, HsigmaRPhi);
	fvecToFarray(HsigmaRZ_fvec, HsigmaRZ);
	fvecToFarray(HsigmaPhiPhi_fvec, HsigmaPhiPhi);
	fvecToFarray(HsigmaZPhi_fvec, HsigmaZPhi);
	fvecToFarray(HsigmaZZ_fvec, HsigmaZZ);

 	vPhi2 = vPhi1;
 	vPhi2[0] += DeltavPhi;

	// some info
	MLOG(verb2) << "Disk gaussian normalizations: " << fk << " : " << (1-fk);
	MLOG(verb2) << "Second disk gaussian offset:  " << DeltavPhi;

	MLOG(verb2) << "vR coefficients:              " << vR_fvec;
	MLOG(verb2) << "vZ coefficients:              " << vZ_fvec;
	MLOG(verb2) << "sigmaRR coefficients:         " << sigmaRR_fvec;
	MLOG(verb2) << "sigmaRPhi coefficients:       " << sigmaRPhi_fvec;
	MLOG(verb2) << "sigmaRZ coefficients:         " << sigmaRZ_fvec;
	MLOG(verb2) << "sigmaPhiPhi1 coefficients:    " << sigmaPhiPhi1_fvec;
	MLOG(verb2) << "sigmaPhiPhi2 coefficients:    " << sigmaPhiPhi2_fvec;
	MLOG(verb2) << "sigmaZPhi coefficients:       " << sigmaZPhi_fvec;
	MLOG(verb2) << "sigmaZZ coefficients:         " << sigmaZZ_fvec;

	MLOG(verb2) << "HvR coefficients:             " << HvR_fvec;
	MLOG(verb2) << "HvZ coefficients:             " << HvZ_fvec;
	MLOG(verb2) << "HsigmaRR coefficients:        " << HsigmaRR_fvec;
	MLOG(verb2) << "HsigmaRPhi coefficients:      " << HsigmaRPhi_fvec;
	MLOG(verb2) << "HsigmaRZ coefficients:        " << HsigmaRZ_fvec;
	MLOG(verb2) << "HsigmaPhiPhi coefficients:    " << HsigmaPhiPhi_fvec;
	MLOG(verb2) << "HsigmaZPhi coefficients:      " << HsigmaZPhi_fvec;
	MLOG(verb2) << "HsigmaZZ coefficients:        " << HsigmaZZ_fvec;

	return true;
}
