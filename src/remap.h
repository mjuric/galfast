#pragma once

typedef std::list<V2> Vertices;

class Poly {
public:
	Vertices vs;
	double xmin, xmax;
public:
	void push_back(V2 v);

	Poly cutflip(int x);
	double area() const;

	inline Poly cutMyself();
	void deintegerize();
	Poly unflip();
};

inline std::ostream &operator <<(std::ostream &o, const Poly &p) { FOREACHj(v, p.vs) o << *v << " "; o << '[' << p.area() << ']'; return o; }

template<typename T, typename MapSpecification>
T mapPixel(V2 pixel[4], MapSpecification &ms)
{
	Poly base;
	FOR(0,4) { base.push_back(pixel[i]); }

//	cout << "input base : " << base << "\n";
	ASSERT(base.area() >= 0);

	T mass = 0;
	base.deintegerize();
	for(int i = (int)base.xmin + 1; !base.vs.empty(); i++) {
		Poly left = base.cutflip(i);

//		cout << "left : " << left << "\n";
		ASSERT(left.area() >= 0);

		left.deintegerize();
		for(int j = (int)left.xmin + 1; !left.vs.empty(); j++) {
			Poly bottom = left.cutflip(j);

			ASSERT(bottom.area() >= 0);

//			cout << "   ---- j = " << j << "\n";
//			cout << "   bottom : " << bottom.unflip().unflip() << "\n";
//			cout << "   base   : " << left << "\n";
			mass = mass + ms(i, j) * bottom.area();
		}
	}
	return mass;
}

//
//	Public functions
//

template<typename T>
struct SafeMapSpecification {
protected:
	const ImageBase<T> &img;
public:
	SafeMapSpecification(const ImageBase<T> &i) : img(i) {}

	// return the value of the requested pixel
	T operator()(int i, int j) const {
		if(i < 0 || j < 0 || i >= img.width() || j >= img.height()) { return 0; }
		return img(i, j);
	}

	// map pixel from destination image to original image
	// return value are pixel coordinates
	virtual V2 invMap(const V2 p) = 0;

	// volume transformation rho(new) = volumeFactor * rho(old)
	virtual double volumeFactor(const V2 p) = 0;
};


template<typename T, typename MapSpecification>
ImageBase<T> &mapTransform(
	ImageBase<T> &out,
	MapSpecification &ms,				// map specification
	V2 minn, V2 maxx,				// output rectangle
	double dx					// output scale factor
	)
{
	V2 size = round(maxx - minn) / dx;		// calculate image size
	maxx = minn + size*dx;				// adjust max values to rounded up versions

	out.resize(size);
	out.setOrigin(minn);
	out.setScale(dx);

	DEBUG(verbose, "map [ size = " << size << ", origin = " << minn << ", scale = " << dx << " ]");

	double offset[4][2] = { {0, 0}, {dx, 0}, {dx, dx}, {0, dx} };

	V2 r[4];
	for(int j = 0; j != int(size.y); j++) {
		const double y = minn.y + j * dx;

		for(int i = 0; i != int(size.x); i++) {
			const double x = minn.x + i * dx;

			FORj(k, 0, 4) { // map a destination pixel back to the originating image
				V2 p(x + offset[k][0], y + offset[k][1]);
				r[k] = ms.invMap(p);
			}

			out(i, j) = mapPixel<T, MapSpecification>(r, ms);			 // map the density into destination pixel
			out(i, j) = out(i, j) / ms.volumeFactor(V2(x + dx/2, y + dx/2)); // weighting of the space
		}
	}

	return out;
}
