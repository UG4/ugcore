
#ifndef GNUPLOT_OUTPUT_H_
#define GNUPLOT_OUTPUT_H_

template<typename Vector_type, typename postype>
void WriteVectorGnuplot(std::string filename, const Vector_type &v,
		postype *positions, int dimensions, const Vector_type *compareVec=NULL)
{
	size_t N = A.num_rows();
	std::fstream f((filename+".sh").c_str(), std::ios::out);
	f << "#!/bin/bash\n"
			"cat > gnuplotTemporaryFile <<EOF\n"
			"set dgrid3d N,N\n"
			"set style data lines\n"
			"set pm3d\n"
			"splot \"-\" pal\n";
	if(dimensions == 1)
		for(size_t i=0; i < N; i++)
			f << positions[i][0] << " " << v[i] << "\n";
	else if(dimensions == 2)
		for(size_t i=0; i < N; i++)
			f << positions[i][0] << " " << positions[i][1] << " " << v[i] << "\n";
	else
		for(size_t i=0; i < N; i++)
			f << positions[i][0] << " " << positions[i][1] << " " << positions[i][2] << " " << v[i] << "\n";
	f <<	"e\n"
			"EOF\n"
			"cat gnuplotTemporaryFile | gnuplot -persist\n"
			"rm gnuplotTemporaryFile\n";
}


#endif /* GNUPLOT_OUTPUT_H_ */
