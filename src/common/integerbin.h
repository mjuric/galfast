#ifndef __integer_bin_h
#define __integer_bin_h

#include <valarray>
#include <astro/image/indexers.h>
#include <astro/util.h>
#include "ximage.h"

/*!
	\brief 	Bins an image of size (x, y) with bin size being (xbin, ybin) pixels, 
			and keeps the center of the old image in the center of the new image,
			while truncating the edges.
			Assuming that physical coordinates (0, 0) are at the bottom left corner of
			the central pixel of the original image, they stay there in the binned image.
*/
template<typename T>
void integer_bin(std::valarray<T> &outa, std::valarray<T> &imga, int &x, int &y, const int xbin, const int ybin)
{
	const int xc = x / 2;
	const int yc = y / 2;

	// calculate binned image size
	// note that if the halfsize is not an exact multiple of
	// bin factors, there'll be some truncation of the large image
	// at the edges, which must be taken into account in the binning loop
	// Also not that this is NOT equivalent to x / xbin, because of integer arithmetic
	const int xb = (xc / xbin) * 2;
	const int yb = (yc / ybin) * 2;
	const int xo = xc % xbin;
	const int yo = yc % ybin;

//	cerr << io::format(" ", true) << xb << yb << xo << yo << "\n";

	outa.resize(xb*yb);
	outa = T(0);

	peyton::image::ind2<T> img(imga, x, y);
	peyton::image::ind2<T> out(outa, xb, yb);

	FORj(i, 0, xb*xbin) FORj(j, 0, yb*ybin)
	{
//		if(i/xbin - xb/2  == 0 && j/ybin -yb / 2  ==0)
//		{
//			cout << xo << " " << yo << " " << xo + i << " " << yo + j << "\n";
//		}

		out(i/xbin,j/ybin) += img(xo + i, yo + j);
	}
	
	x = xb;
	y = yb;

//	exit(-1);
}

inline void integer_bin(XImage &img, const int xbin, const int ybin = -1)
{
	int x = img.x(), y = img.y();
	std::valarray<float> tmp;

	integer_bin(tmp, img.array, x, y, xbin, ybin != -1 ? ybin : xbin);

	img.resize(x, y);
	img.array = tmp;
}

#endif
