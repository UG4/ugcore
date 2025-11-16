/*
 * Copyright (c) 2015:  G-CSC, Goethe University Frankfurt
 * Author: Markus Breit
 * 
 * This file is part of UG4.
 * 
 * UG4 is free software: you can redistribute it and/or modify it under the
 * terms of the GNU Lesser General Public License version 3 (as published by the
 * Free Software Foundation) with the following additional attribution
 * requirements (according to LGPL/GPL v3 §7):
 * 
 * (1) The following notice must be displayed in the Appropriate Legal Notices
 * of covered and combined works: "Based on UG4 (www.ug4.org/license)".
 * 
 * (2) The following notice must be displayed at a prominent place in the
 * terminal output of covered works: "Based on UG4 (www.ug4.org/license)".
 * 
 * (3) The following bibliography is recommended for citation and must be
 * preserved in all covered files:
 * "Reiter, S., Vogel, A., Heppner, I., Rupp, M., and Wittum, G. A massively
 *   parallel geometric multigrid solver on hierarchically distributed grids.
 *   Computing and visualization in science 16, 4 (2013), 151-164"
 * "Vogel, A., Reiter, S., Rupp, M., Nägel, A., and Wittum, G. UG4 -- a novel
 *   flexible software system for simulating pde based models on high performance
 *   computers. Computing and visualization in science 16, 4 (2013), 165-179"
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU Lesser General Public License for more details.
 */

#include <string>
#include <fstream>
#include <sstream>
#include <vector>
#include <queue>
#include <map>
#include "file_io_swc.h"
#include "common/error.h"
#include "common/util/string_util.h"
#include "lib_grid/global_attachments.h"
#include "lib_grid/algorithms/subset_util.h" // EraseEmptySubsets
#include "lib_grid/algorithms/element_side_util.h" // GetOpposingSide
#include <boost/lexical_cast.hpp>

using namespace std;

namespace ug{


bool FileReaderSWC::load_file(const char* fileName)
{
    m_vPts.clear();

    std::ifstream inFile(fileName);
    if (!inFile)
    {
    	UG_LOGN("SWC input file '" << fileName << "' could not be opened for reading.");
    	return false;
    }

    // read line by line
    std::string line;
    size_t lineCnt = 0;
    size_t curInd = 0;
    std::map<int, size_t> indexMap;
    while (std::getline(inFile, line))
    {
        ++lineCnt;

        // trim whitespace
        line = TrimString(line);

        // ignore anything from possible '#' onwards
        size_t nChar = line.size();
        for (size_t i = 0; i < nChar; ++i)
        {
            if (line.at(i) == '#')
            {
                line = line.substr(0, i);
                break;
            }
        }

        // empty lines can be ignored
        if (line.empty()) continue;

        // split the line into tokens
        std::istringstream buf(line);
        std::istream_iterator<std::string> beg(buf), end;
        std::vector<std::string> strs(beg, end);

        // assert number of tokens is correct
       if (strs.size() != 7)
       {
    	   UG_LOGN("Error reading SWC file '" << fileName << "': Line " << lineCnt
    		   << " does not contain exactly 7 values.");
    	   return false;
       }


        // collect data for SWC point
        m_vPts.resize(m_vPts.size() + 1);
        SWCPoint& pt = m_vPts.back();

        // get index from file and map to our index
        indexMap[boost::lexical_cast<int>(strs[0])] = curInd;

        // type
        int type = boost::lexical_cast<int>(strs[1]);
        switch (type)
        {
            case 0: pt.type = swc_types::SWC_UNDF; break;
            case 1: pt.type = swc_types::SWC_SOMA; break;
            case 2: pt.type = swc_types::SWC_AXON; break;
            case 3: pt.type = swc_types::SWC_DEND; break;
            case 4: pt.type = swc_types::SWC_APIC; break;
            case 5: pt.type = swc_types::SWC_FORK; break;
            case 6: pt.type = swc_types::SWC_END; break;
            default: pt.type = swc_types::SWC_CUSTOM;
        }

        // coordinates
        pt.coords.x() = boost::lexical_cast<number>(strs[2]);
        pt.coords.y() = boost::lexical_cast<number>(strs[3]);
        pt.coords.z() = boost::lexical_cast<number>(strs[4]);

        // radius
        pt.radius = boost::lexical_cast<number>(strs[5]);

        // connections
        int conn = boost::lexical_cast<int>(strs[6]);
        if (conn >= 0)
        {
            std::map<int, size_t>::const_iterator it = indexMap.find(conn);
            if (it == indexMap.end())
            {
            	UG_LOGN("Error reading SWC file '" << fileName << "': Line " << lineCnt
            		<< " refers to unknown parent index " << conn << ".");
            	return false;
            }

            size_t parentID = indexMap[conn];
            pt.conns.push_back(parentID);
            m_vPts[parentID].conns.push_back(curInd);
        }

        // increase current point index
        ++curInd;
    }

    return true;
}



bool FileReaderSWC::create_grid(Grid& g, ISubsetHandler* pSH, number scale_length)
{
	if (!g.has_vertex_attachment(aPosition))
		g.attach_to_vertices(aPosition);
	Grid::VertexAttachmentAccessor<APosition> aaPos(g, aPosition);

	// handle diameter attachment
	if (!GlobalAttachments::is_declared("diameter"))
	{
		try {GlobalAttachments::declare_attachment<ANumber>("diameter");}
		catch (ug::UGError& err)
		{
			UG_LOGN(err.get_msg());
			return false;
		}
	}

	ANumber aDiam = GlobalAttachments::attachment<ANumber>("diameter");
	if (!g.has_vertex_attachment(aDiam))
		g.attach_to_vertices(aDiam, true);

	Grid::AttachmentAccessor<Vertex, ANumber> aaDiam(g, aDiam);

	// create grid
	const size_t nP = m_vPts.size();
	std::vector<Vertex*> vrts(nP, nullptr);
	for (size_t i = 0; i < nP; ++i)
	{
		const SWCPoint& pt = m_vPts[i];

		// create vertex and save
		Vertex* v = vrts[i] = *g.create<RegularVertex>();
		VecScale(aaPos[v], pt.coords, scale_length);
		pSH->assign_subset(v, pt.type - 1);
		aaDiam[v] = 2 * pt.radius * scale_length;

		// create edge connections to already created vertices
		for (size_t j = 0; j < pt.conns.size(); ++j)
		{
			if (pt.conns[j] < i)
			{
				Edge* e = *g.create<RegularEdge>(EdgeDescriptor(vrts[pt.conns[j]], v));
				pSH->assign_subset(e, pt.type - 1);
			}
		}
	}

	// final subset management
	AssignSubsetColors(*pSH);
	pSH->set_subset_name("soma", 0);
	pSH->set_subset_name("axon", 1);
	pSH->set_subset_name("dend", 2);
	pSH->set_subset_name("apic", 3);
	pSH->set_subset_name("fork", 4);
	pSH->set_subset_name("end", 5);
	pSH->set_subset_name("custom", 6);
	EraseEmptySubsets(*pSH);

	return true;
}



const std::vector<swc_types::SWCPoint>& FileReaderSWC::swc_points() const
{
	return m_vPts;
}



std::vector<swc_types::SWCPoint>& FileReaderSWC::swc_points()
{
	return m_vPts;
}



bool LoadGridFromSWC(Grid& g, ISubsetHandler* pSH, const char* fileName, AVector3& aPos)
{
	// now read the file
	FileReaderSWC reader;
	return reader.load_file(fileName) && reader.create_grid(g, pSH, 1.0);
}



bool ExportGridToSWC(Grid& g, ISubsetHandler* pSH, const char* fileName, AVector3& aPos)
{
	// get access to positions
	if (!g.has_vertex_attachment(aPosition))
	{
		UG_LOGN("Position attachment not attached to grid.");
		return false;
	}
	Grid::VertexAttachmentAccessor<APosition> aaPos(g, aPosition);

	// get access to diameter attachment
	if (!GlobalAttachments::is_declared("diameter"))
	{
		try {GlobalAttachments::declare_attachment<ANumber>("diameter");}
		catch (ug::UGError& err)
		{
			UG_LOGN(err.get_msg());
			return false;
		}
	}
	ANumber aDiam = GlobalAttachments::attachment<ANumber>("diameter");
	if (!g.has_vertex_attachment(aDiam))
	{
		UG_LOGN("WARNING: No diameter attachment attached to grid. "
			"Will use 1.0 as default diameter.");
		g.attach_to_vertices_dv(aDiam, 1.0);
	}
	Grid::AttachmentAccessor<Vertex, ANumber> aaDiam(g, aDiam);

	// analyze subset names to find out corresponding swc-types
	size_t nss = pSH->num_subsets();
	std::vector<size_t> vType(nss);
	bool soma_subset_present = false;
	for (size_t i = 0; i < nss; ++i)
	{
		std::string name(pSH->get_subset_name(i));
		std::transform(name.begin(), name.end(), name.begin(), ::toupper);
		if (name.find("SOMA") != std::string::npos)
		{
			soma_subset_present = true;
			vType[i] = 1;
		}
		else if (name.find("AXON") != std::string::npos)
			vType[i] = 2;
		else if (name.find("APIC") != std::string::npos)
			vType[i] = 4;
		else if (name.find("DEND") != std::string::npos)
			vType[i] = 3;
    else if (name.find("END") != std::string::npos)
      vType[i] = 6;
    else if (name.find("FORK") != std::string::npos)
      vType[i] = 5;
    else if (name.find("CUSTOM") != std::string::npos)
      vType[i] = 7;
		else vType[i] = 0;
	}

	if (!soma_subset_present)
		UG_LOGN("Warning: No somatic subset could be identified.")

	if (g.begin<Vertex>() == g.end<Vertex>())
	{
		UG_LOGN("Warning: No vertices contained in grid.")
		return true;
	}

	// find soma vertex (if identifiable)
	Vertex* start = *g.begin<Vertex>();
	if (soma_subset_present)
	{
		g.begin_marking();
		std::queue<Vertex*> q; // corresponds to breadth-first
		q.push(start);
		while (!q.empty())
		{
			Vertex* v = q.front();
			if (vType[pSH->get_subset_index(v)] == 1) break;
			g.mark(v);
			q.pop();

			// push neighboring elems to queue
			Grid::traits<Edge>::secure_container edges;
			g.associated_elements(edges, v);

			size_t sz = edges.size();
			for (size_t e = 0; e < sz; ++e)
			{
				Vertex* otherEnd = GetOpposingSide(g, edges[e], v);
				if (!g.is_marked(otherEnd))
					q.push(otherEnd);
			}
		}
		g.end_marking();

		if (q.empty())
			UG_LOGN("Warning: No soma vertex could be found in the requested neuron.")
		else
			start = q.front();
	}

	// write the neuron to file
	std::ofstream outFile(fileName, std::ios::out);
	UG_COND_THROW(!outFile.is_open(), "Could not open output file '" << fileName << "'.");

	outFile << "# This file has been generated by UG4." << std::endl;

	std::stack<std::pair<Vertex*, int> > stack; // corresponds to depth-first
	stack.push(std::make_pair(start, -1));

	g.begin_marking();
	int ind = 0;   // by convention, swc starts with index 1
	bool all_types_identified = true;
	while (!stack.empty())
	{
		// get all infos regarding vertex
		std::pair<Vertex*, int>& info = stack.top();
		Vertex* v = info.first;
		int conn = info.second;
		stack.pop();

		// mark curr vrt
		g.mark(v);

		size_t type = vType[pSH->get_subset_index(v)];
		if (!type) all_types_identified = false;

		const vector3& coord = aaPos[v];

		number radius = 0.5*aaDiam[v];

		// write line to file
		outFile << ++ind << " " << type << " "
			<< coord[0] << " " << coord[1] << " " << coord[2] << " "
			<< radius << " " << conn << std::endl;

		// push neighboring elems to queue
		Grid::traits<Edge>::secure_container edges;
		g.associated_elements(edges, v);

		size_t sz = edges.size();
		for (size_t e = 0; e < sz; ++e)
		{
			Vertex* otherEnd = GetOpposingSide(g, edges[e], v);
			if (!g.is_marked(otherEnd))
				stack.push(std::make_pair(otherEnd, ind));
		}
	}
	g.end_marking();

	if (!all_types_identified)
		UG_LOGN("WARNING: Some vertex type(s) - soma, dendrite, axon, etc. -\n"
			"could not be identified using the subset names.\n"
			<< "To ensure correct types in the resulting swc file, the ugx subset names\n"
			"need to contain one of the strings \"SOMA\", \"AXON\", \"DEND\", \"APIC\",\n"
			"upper/lower case can be ignored.");

	outFile.close();

	return true;
}



}//	end of namespace
