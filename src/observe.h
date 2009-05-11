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

// add Fe/H information
class os_FeH : public osink, os_FeH_data
{
public:
	virtual size_t process(otable &in, size_t begin, size_t end, rng_t &rng);
	virtual bool init(const Config &cfg, otable &t);
	virtual const std::string &name() const { static std::string s("FeH"); return s; }

	os_FeH() : osink()
	{
		prov.insert("FeH");
		req.insert("comp");
		req.insert("XYZ");
	}
};

// add Fe/H information
class os_fixedFeH : public osink
{
	protected:
		float fixedFeH;

	public:
		virtual size_t process(otable &in, size_t begin, size_t end, rng_t &rng);
		virtual bool init(const Config &cfg, otable &t);
		virtual const std::string &name() const { static std::string s("fixedFeH"); return s; }

		os_fixedFeH() : osink(), fixedFeH(0)
		{
			prov.insert("FeH");
		}
};

// convert velocities to proper motions
class os_vel2pm : public osink , public os_vel2pm_data
{
public:
	virtual size_t process(otable &in, size_t begin, size_t end, rng_t &rng);
	virtual bool init(const Config &cfg, otable &t);
	virtual const std::string &name() const { static std::string s("vel2pm"); return s; }

	os_vel2pm() : osink()
	{
		coordsys=GAL;
		req.insert("lb");
		req.insert("XYZ");
		req.insert("vcyl");
	}
};
	

// add kinematic information
class os_kinTMIII : public osink, os_kinTMIII_data, os_kinTMIII_data_groupedA
{	
	public:
		void add_dispersion(float v[3], float Rsquared, float Z, farray5 *ellip[6], gpu_rng_t &rng);
		void compute_means(float v[3], float Rsquared, float Z, farray5 *means[3]);

		void get_disk_kinematics(float v[3], float Rsquared, float Z, gpu_rng_t &rng, bool &firstGaussian);
		void get_halo_kinematics(float v[3], float Rsquared, float Z, gpu_rng_t &rng);

	public:
		virtual size_t process(otable &in, size_t begin, size_t end, rng_t &rng);
		virtual bool init(const Config &cfg, otable &t);
		virtual const std::string &name() const { static std::string s("kinTMIII"); return s; }

		os_kinTMIII() : osink()
		{
			prov.insert("vcyl");
			req.insert("comp");
			req.insert("XYZ");
		}
};



