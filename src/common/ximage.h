#ifndef __ximage_h
#define __ximage_h

#include <astro/image/indexers.h>

class XImage : public peyton::image::ind2<float>
{
public:
	std::valarray<float> array;

	XImage(int x = 0, int y = 0, float initial = 0.)
		: array(initial, x*y), peyton::image::ind2<float>(array, x, y)
	{}

	bool is_compatible_with(const XImage &img) const
	{
		return img.array.size() == array.size()
			&& img.x() == x()
			&& img.y() == y();
	}
	
	void be_compatible_with(const XImage &img)
	{
		array.resize(img.array.size());
		x() = img.x();
		y() = img.y();
	}


	bool ensure_compat(const XImage &img)
	{
		if(is_compatible_with(img)) return false;
		be_compatible_with(img);
		return true;
	}
	
	void resize(int x_, int y_, float initial = 0.)
	{
		x() = x_;
		y() = y_;
		array.resize(x_*y_);
		array = initial;
	}
	
	XImage &operator=(const XImage &img)
	{
		ensure_compat(img);
		array = img.array;
		return *this;
	}
};

#endif
