	void test()
	{
		int n = 100;
		gsl::matrix X(n, 3), cov(3, 3);
		gsl::vector y(n), p(3);

		double a = 1.3, b = 2.33, c = -1.11;
		int k = 0;
		FORj(i, 0, 10) FORj(j, 0, 10)
		{
			y(k) = a + b*i + c*j + 1*(float(rand())/RAND_MAX - 1.);
			X(k, 0) = 1;
			X(k, 1) = i;
			X(k, 2) = j;
			k++;
		}
		
		double chisq;
		gsl::multifit (n, 3) . linear(X, y, p, cov, chisq);

		cout << "Fit : " << p(0) << " " << p(1) << " " << p(2) << "\n";
		cout << "Covariance : ";
		gsl_matrix_fprintf(stdout, cov, "%f");
		cout << "Chisq : " << chisq << "\n";
	}

	void test2()
	{
		float zero = 0.;
		float inf = 1./zero;
		cout << "inf : " << inf << "\n";
		cout << "1/inf : " << 1./inf << "\n";
		cout << "inf+1 : " << 1. + inf << "\n";
	}

	void test1()
	{
		text_input_or_die  (ts, "mdwarfs.dat");

		int size = 4000000;
		valarray<int> run(size);
		valarray<float> u(size), g(size), r(size), i(size), z(size), Ar(size), l(size), b(size);

		cout << "Loading...\n";

		load(ts, run,0, u,10, g,11, r,12, i,13, z,14, Ar,9, l,7, b,8, tickk,-1);

		cout << "Loaded " << u.size() << " records\n";
	}
	
	void test3()
	{
		//
		// Used while testing the system::Config class
		//
		
		peyton::system::Config c("x.conf", "test = gaga\nt2 = gugu");
		std::string bla = c["bla"];
		cout << "bla = " << bla << "\n";
		cout << util::expand_text("This is a [$xx] test.", c) << "\n";
		
		exit(0);
	}
#if 0	
	void test4(int argc, char **argv)
	{
		//
		// Used to test peyton::system::Options
		//
		
		Options options;
		//           key      longname   short value  argument          default              description
		options.option("remote", "noremote",  0,  "0",   Option::none,       Option::nodefault,   "Do not start remote servers");
		options.option("remote", "remote",   'r', "1",   Option::none,       Option::nodefault,   "Start remote servers (default)");
		options.option("test",   "test",     't', "1",   Option::none,       "0",   "Just do testing, no actual running");
		options.option("usage",  "usage",    'u', "all", Option::optional, "",    "Do not start remote servers");

		try {
			options.parse(argc, argv);
		} catch(EOptions &e) {
			//e.print();
			cerr << "\n" << options.usage();
			exit(0);
		}

		FOREACH(options)
		{
			cout << (*i).first << " -> '" << (*i).second << "'\n";
		}
		exit(0);
	}
#endif

#if 0
			magick::LogTransfer <typeof(counts)> trans(counts);
			magick::SineFilter red;

//			FOREACH(counts) { *i = i.x + i.y; }

			magick::write(counts, countsfile + ".png", &trans, &red);
//			magick::write(counts, countsfile + ".png", (magick::LogTransfer <typeof(counts)> *)NULL, &red);
#endif
