/*
 * vtkoutput.cpp
 *
 *  Created on: 25.04.2012
 *      Author: andreasvogel
 */

#include "vtkoutput.h"

namespace ug{

void VTKOutput::
write_empty_grid_piece(VTKFileWriter& File)
{
//	write that no elements are in the grid
	int n = 0;
	File.write("    <Piece NumberOfPoints=\"0\" NumberOfCells=\"0\">\n");
	File.write("      <Points>\n");
	File.write("        <DataArray type=\"Float32\" NumberOfComponents=\"3\" format=\"binary\">\n");
	File.write_base64(n);
	File.write("\n        </DataArray>\n");
	File.write("      </Points>\n");
	File.write("      <Cells>\n");
	File.write("        <DataArray type=\"Int32\" Name=\"connectivity\" format=\"binary\">\n");
	File.write_base64(n);
	File.write("\n        </DataArray>\n");
	File.write("        <DataArray type=\"Int32\" Name=\"offsets\" format=\"binary\">\n");
	File.write_base64(n);
	File.write("\n        </DataArray>\n");
	File.write("        <DataArray type=\"Int8\" Name=\"types\" format=\"binary\">\n");
	File.write_base64(n);
	File.write("\n        </DataArray>\n");
	File.write("      </Cells>\n");
	File.write("    </Piece>\n");
}

////////////////////////////////////////////////////////////////////////////////
// FileNames
////////////////////////////////////////////////////////////////////////////////


void VTKOutput::
vtu_filename(std::string& nameOut, std::string nameIn, int rank,
             int si, int maxSi, int step)
{
//	copy name
	nameOut = nameIn.substr(0, nameIn.find_first_of('.'));

#ifdef UG_PARALLEL
// 	process index
	if(pcl::GetNumProcesses() > 1)
		AppendCounterToString(nameOut, "_p", rank, pcl::GetNumProcesses() - 1);
#endif

// 	subset index
	if(si >= 0)
		AppendCounterToString(nameOut, "_s", si, maxSi);

// 	time index
	if(step >= 0)
		AppendCounterToString(nameOut, "_t", (int) step);

// 	add file extension
	nameOut.append(".vtu");
}


void VTKOutput::
pvtu_filename(std::string& nameOut, std::string nameIn,
              int si, int maxSi, int step)
{
//	copy name
	nameOut = nameIn.substr(0, nameIn.find_first_of('.'));

// 	subset index
	if(si >= 0)
		AppendCounterToString(nameOut, "_s", si, maxSi);

// 	time index
	if(step >= 0)
		AppendCounterToString(nameOut, "_t", (int) step);

// 	add file extension
	nameOut.append(".pvtu");
}


void VTKOutput::
pvd_filename(std::string& nameOut, std::string nameIn)
{
//	copy name
	nameOut = nameIn.substr(0, nameIn.find_first_of('.'));

// 	add file extension
	nameOut.append(".pvd");
}


void VTKOutput::
pvd_time_filename(std::string& nameOut, std::string nameIn, int step)
{
//	copy name
	nameOut = nameIn.substr(0, nameIn.find_first_of('.'));

// 	time index
	if(step >= 0)
		AppendCounterToString(nameOut, "_t", (int) step);

// 	add file extension
	nameOut.append(".pvd");
}

void VTKOutput::
write_subset_pvd(int numSubset, const std::string& filename, int step, number time)
{
//	file pointer
	FILE* file;

//	to store name of file
	std::string name;

//	get rank, outproc bool and number of processes
	bool isOutputProc = GetLogAssistant().is_output_process();
	int rank = 0;
	int numProcs = 1;

#ifdef UG_PARALLEL
	rank = pcl::GetProcRank();
	numProcs = pcl::GetNumProcesses();
#endif

//	only output proc writes this file
	if (isOutputProc)
	{
	//	get file name
		if(step >= 0) pvd_time_filename(name, filename, step);
		else pvd_filename(name, filename);

	//	open file
		file = fopen(name.c_str(), "w");
		if (file == NULL)
			UG_THROW("VTKOutput: Cannot print to file.");

	// 	Write beginning of file
		fprintf(file, "<?xml version=\"1.0\"?>\n");
		fprintf(file, "<VTKFile type=\"Collection\" version=\"0.1\">\n");
		fprintf(file, "  <Collection>\n");

	// 	Include files for all subsets
		for(int si = 0; si < numSubset; ++si)
		{
			vtu_filename(name, filename, rank, si, numSubset-1, step);
			if(numProcs > 1) pvtu_filename(name, filename, si, numSubset-1, step);

			name = FilenameWithoutPath(name);
			fprintf(file, "  <DataSet timestep=\"%g\" part=\"%d\" file=\"%s\"/>\n",
			        		time, si, name.c_str());
		}

	//	write closing tag
		fprintf(file, "  </Collection>\n");
		fprintf(file, "</VTKFile>\n");
		fclose(file);
	}

	if (isOutputProc && numProcs > 1)
	{
		std::string procName(filename);
		procName.append("_processwise");

	//	get file name
		if(step >= 0) pvd_time_filename(name, filename, step);
		else pvd_filename(name, filename);

	//	open File
		file = fopen(name.c_str(), "w");
		if (file == NULL)
			UG_THROW("VTKOutput: Cannot print to file.");

	// 	Write to file
		fprintf(file, "<?xml version=\"1.0\"?>\n");
		fprintf(file, "<VTKFile type=\"Collection\" version=\"0.1\">\n");
		fprintf(file, "  <Collection>\n");

	// 	include files from all procs
		for (int r = 0; r < numProcs; r++)
			for(int si = 0; si < numSubset; ++si)
			{
				vtu_filename(name, filename, rank, si, numSubset-1, step);
				if(numProcs > 1) pvtu_filename(name, filename, si, numSubset-1, step);

				name = FilenameWithoutPath(name);
				fprintf(file, "  <DataSet timestep=\"%g\" part=\"%d\" file=\"%s\"/>\n",
				        	time, r, name.c_str());
			}

	//	end file
		fprintf(file, "  </Collection>\n");
		fprintf(file, "</VTKFile>\n");
		fclose(file);
	}
}


void VTKOutput::
select_nodal_scalar(const char* fctName, const char* name)
{
	std::vector<std::string> tokens;
	std::string fctString(fctName);
	TokenizeString(fctString, tokens, ',');
	if(tokens.size() != 1)
		UG_THROW("VTK:select_nodal_scalar: In order to select"
				" a nodal scalar for output to vtk,"
				" exactly one function components must be chosen.");

//	skip already selected
	for(size_t i = 0; i < tokens.size(); ++i)
	{
		RemoveWhitespaceFromString(tokens[i]);
		for(size_t j = 0; j < m_vSymbFct.size(); ++j)
		{
		//	skip function selected twice, if and only if same name given
			if(m_vSymbFct[j].first.size() == tokens[i].size() &&
				m_vSymbFct[j].first == tokens[i])
			{
				if(m_vSymbFct[j].second == name &&
					m_vSymbFct[j].second.size() == strlen(name)){
					return;
				}
				else{
					UG_THROW("VTK:select_nodal_scalar: Selecting component "
							<< tokens[i] << " again, but with different name " <<
							name << " instead of already scheduled "
							<< m_vSymbFct[j].second);
				}
			}

		//	check if name is not in use
			if(m_vSymbFct[j].second == name &&
				m_vSymbFct[j].second.size() == strlen(name)){
				UG_THROW("VTK:select_nodal_scalar: Selecting component "
						<< tokens[i] << ", but with already used name " <<
						name << ". This is not allowed, use different name.");
			}
		}
	}

	m_vSymbFct.push_back(std::pair<std::string, std::string>(fctName, name));
}


void VTKOutput::
select_nodal_vector(const char* fctNames, const char* name)
{
	std::vector<std::string> tokens;
	std::string fctString(fctNames);
	TokenizeString(fctString, tokens, ',');
	if(tokens.size() > 3)
		UG_THROW("VTK:select_nodal_vector: In order to select"
				" a nodal vector for output to vtk,"
				" maximal #dim function components must be chosen.");

	for(size_t j = 0; j < m_vSymbFct.size(); ++j)
	{
	//	check if name is not in use
		if(m_vSymbFct[j].second == name &&
			m_vSymbFct[j].second.size() == strlen(name))
			UG_THROW("VTK:select_nodal_vector: Using name " << name <<
			       " that is already used by other data is not allowed.");
	}

	m_vSymbFct.push_back(std::pair<std::string, std::string>(fctNames, name));
}


} // end namespace ug
