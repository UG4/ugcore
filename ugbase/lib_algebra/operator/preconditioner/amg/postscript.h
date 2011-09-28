/**
 * \file postscript.h
 *
 * \author Martin Rupp
 *
 * \date 01.08.2010
 *
 * Goethe-Center for Scientific Computing 2010.
 */


#ifndef __H__UG__LIB_DISC__POSTSCRIPT_H__
#define __H__UG__LIB_DISC__POSTSCRIPT_H__

#include <fstream>
#include <string>

namespace ug{

class postscript
{
public:
	postscript()
	{
		bounds_left = 1e18;
		bounds_right = -1e18;
		bounds_top = 1e18;
		bounds_bottom = -1e18;
		mindist=1e18;
	}

	~postscript()
	{
		double scale = 1;
		scale = std::max( 400/(bounds_right-bounds_left), 400/(bounds_top-bounds_bottom) );
		file <<		"%%BoundingBox: " << bounds_left*scale*1.05 << " " << bounds_top*scale*1.05 << " " << bounds_right*scale*1.05 << " " << bounds_bottom*scale*1.05 << "\n"
					"%%Pages: 1\n"
					"%%DocumentsFonts: Monaco\n"
					"%%Copyright 2010 G-CSC - All Rights Reserved Worldwide\n"
					"%%EndComments\n\n";

		file <<		"1 setlinejoin\n"
					"1 setlinecap\n"
					"/Monaco findfont 10 scalefont setfont\n\n";

		file << 	"/M {moveto} def\n"
					"/S {lineto stroke} def\n"
					"/L {lineto} def\n"
					"/C {closepath fill} def\n"
					"/N {newpath} def\n"
					"/R {setrgbcolor} def\n"
					"/W {setlinewidth} def\n\n";

		file <<		"%%Endprolog\n"
					"%\n"
					"%%Page: 1 1\n"
					"%\n\n";

		mindist = sqrt(mindist);
		file << 	mindist * 0.01 << " W\n"
					"/Monaco findfont " << mindist * 0.2 << " scalefont setfont\n0 0 0 R\n";


		file << 	scale << " " << scale << " scale\n";
		file << out.str();
		file << "showpage\n\n%%Trailer";
	}

	bool create(const char *filename)
	{
		file.open((std::string(filename) + ".ps").c_str(), std::ios::out);
		file << 	"%!PS-Adobe-2.0 EPSF-1.2\n"
					"%%Title: " << filename << "\n"
					"%%Creator: ug postscript output\n"
					"%%CreationDate:\n";
		return true;
	}

	void setcolor(double r, double g, double b)
	{
		out << r << " " << g << " " << b << " R\n";
	}
	void move_to(double x, double y)
	{
		extend_bounds(x, y);
		out << "N " << x << " " << y << " M\n";
		last_movetox = x;
		last_movetoy = y;
	}

	void line_to(double x, double y)
	{
		extend_bounds(x, y);
		out << x << " " << y << " S\n";
		double d = ((x-last_movetox)*(x-last_movetox) + (y-last_movetoy)*(y-last_movetoy));
		if(d < mindist) mindist = d;
	}

	void extend_bounds(double x, double y)
	{
		if(x < bounds_left) bounds_left = x;
		if(x > bounds_right) bounds_right = x;
		if(y < bounds_top) bounds_top = y;
		if(y > bounds_bottom) bounds_bottom = y;
	}

	void line(double x1, double y1, double x2, double y2)
	{
		move_to(x1, y1);
		line_to(x2, y2);
	}

	void set_line_width(double width)
	{
		out << width << " W\n";
	}

	void print_text(const char *text)
	{
		out << "(" << text << ") show\n";
	}

	void print_text(const std::string &text)
	{
		out << "(" << text << ") show\n";
	}



private:
	double mindist;
	double last_movetox;
	double last_movetoy;
	double bounds_left, bounds_right, bounds_top, bounds_bottom;
	std::ostringstream out;

	std::fstream file;
};

}
#endif // __H__UG__LIB_DISC__POSTSCRIPT_H__
