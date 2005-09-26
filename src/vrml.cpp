	void vrml_model()
	{
		sdss::RunGeometry geom;

		set<int> runs;
		loadRuns(runs);

		output_or_die (out, "/u/mjuric/public_html/model_bright.wrl")
		out << vrml::header() << "\n";

		vrml::Transform root;
		root.rotation = V4(1, 0, 0, -ctn::pi/2);
		root.scale = V3(-1, -1, 1);

		// lights
		vrml::PointLight pl;
		pl.attenuation = V3(0, 0, 1);
		pl.clone(&root)->setLocation(0, 0, 1);
		pl.clone(&root)->setLocation(0, 0, -1);
		pl.clone(&root)->setLocation(0, 1, 0);
		pl.clone(&root)->setLocation(0, -1, 0);
		pl.clone(&root)->setLocation(1, 0, 0);
		pl.clone(&root)->setLocation(-1, 0, 0);

		// run geometries
		FOREACH(runs) {
			geomDB.getGeometry(*i, geom);
//			if(geom.run != 745) { continue; }

			vrml::IndexedFaceSet *ifs = new vrml::IndexedFaceSet(&root);
			vrml_model_of_run(geom, *ifs);
		}

		//
		// coordinate axis
		//
		V3 colors[] = { V3(1, 0, 0), V3(0, 1, 0), V3(0, 0, 1) };
		V4 rotations[] = { V4(0, 0, 1, -1.570795), V4(1, 0, 0, 1.570795), V4(0, 0, 1, 0) };
		FOR(0, 3) {
			root.addChild(
				vrml::transform()->setRotation(rotations[i])->addChild(
					vrml::cylinder(.02, 5)->setEmissiveColor(colors[i]) )->addChild(
					vrml::transform()->setTranslation(0, 2.5, 0)->addChild(
						vrml::cone(.04, .2)->setEmissiveColor(colors[i])
					)
				)
			);
		}

		//
		// draw galaxy
		//
		root.addChild(
			vrml::transform()->setRotation(1, 0, 0, ctn::pi/2)->addChild(
				vrml::transform()->setTranslation(8, 0, 0)->addChild(
					vrml::cylinder(15, .02)->setEmissiveColor(.3, .3, .5)->setDiffuseColor(.3, .3, .5)->setTransparency(.6) )->addChild(
					vrml::transform()->setScale(1, .46, 1)->addChild(
					 	vrml::sphere(1.5)->setEmissiveColor(1, .5, .5)->setDiffuseColor(1, 1, 0)->setTransparency(0)
					)
				)
			)
		);

		cout << "Outputing...\n";
		out << root << "\n";
	}


	void vrml_galaxy_plane(vrml::IndexedFaceSet &ifs)
	{
		int nSplits = 120;
		double step = rad(360.) / nSplits;

		vector<V3> poly; V3 v;
		for(int i = 0; i != nSplits; i++) {
			v.spherical(3., ctn::pi/2, i*step);
			poly.push_back(v);
		}
		ifs.add(poly);
		ifs.solid = false;
	}

	void vrml_model_of_run(sdss::RunGeometry &geom, vrml::IndexedFaceSet &ifs)
	{
		cout << "Model: " << geom.run << "\n";
		sdss::Mask mask(geom.muStart, geom.muEnd, geom.nu);

		int nStep = (int)(geom.length() / rad(3.));
		Radians step = geom.length() / nStep;

		// front side
//		double d0 = 86.69/1000, d1 = 1397.2/1000;		// m dwarfs, r-i=1-1.1
		double d0 = 1432./1000, d1 = 16621.8/1000;	// bright stars, g-r=0.4-0.5
		Radians l, b;
		FORj(col, 0, 6) {
			Radians nu0 = mask.lo(col), nu1 = mask.hi(col); //+ rad(10);
			Radians nu[] = { nu0, nu0, nu1, nu1 };
			Radians d[] = { d1, d1, d0, d0 };

			vector<V3> poly0(4), poly1(4), polyt0, polyt1, polyb0, polyb1;
			FOR(0, nStep)
			{
				Radians mu0 = geom.muStart + i*step;
				Radians mu1 = mu0 + step;
				
				Radians mu[] = { mu0, mu1, mu1, mu0 };
				Radians l[4], b[4];

				FORj(k, 0, 4) { // front and back faces of the wedge
					coordinates::gcsgal(geom.node, geom.inc, mu[k], nu[k], l[k], b[k]);
					poly0[3-k].spherical(d0, ctn::pi/2 - b[k], l[k]);
					poly1[k].spherical(d1, ctn::pi/2 - b[k], l[k]);
				}
				
				ifs.add(poly0);
				ifs.add(poly1);

				// top and bottom building
				polyt1.push_back(poly1[3]);
				polyt0.push_back(poly0[1]);
				polyb1.push_back(poly1[1]);
				polyb0.push_back(poly0[3]);

				if(i == 0) { // starting side
					vector<V3> poly(4);
					poly[0] = poly0[3]; poly[1] = poly1[0]; poly[2] = poly1[3]; poly[3] = poly0[0];
					ifs.add(poly);
				} else if(i == nStep-2) { // ending side
					vector<V3> poly(4);
					poly[0] = poly0[1]; poly[1] = poly1[2]; poly[2] = poly1[1]; poly[3] = poly0[2];
					ifs.add(poly);
				}
				
			}
			
			// construct top and bottom (respecting the CCW vertex ordering)
			polyt1.insert(polyt1.end(), polyt0.rbegin(), polyt0.rend());
			ifs.add(polyt1);

			polyb0.insert(polyb0.end(), polyb1.rbegin(), polyb1.rend());
			ifs.add(polyb0);
		}
	}
