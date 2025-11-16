/*
 * Copyright (c) 2009-2015:  G-CSC, Goethe University Frankfurt
 * Author: Sebastian Reiter
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

#ifndef __LIBMESH_LOADER_OBJ__
#define __LIBMESH_LOADER_OBJ__

#include <vector>
#include <string>
#include "../../math/ugmath.h"

namespace ug
{

/// \addtogroup ugbase_common_util
/// \{

class LoaderObj
{
	public:
		class Object
		{
			public:
				Object() : m_iMaterialIndex(-1)	{}

				std::string	m_strName;
				std::string	m_strMaterialName;
				int			m_iMaterialIndex;
				std::vector<int>        m_vQuadList;
				std::vector<int>		m_vTriangleList;
				std::vector<int>        m_vQuadListTex;
				std::vector<int>		m_vTriangleListTex;
				std::vector<int>		m_vEdgeList;
		};

		class Material
		{
			public:
				std::string 	m_strName;
				std::string 	m_strTextureDiffuse;
				ug::vector4		m_vDiffuse;
				float			m_fAlpha;
		};

		using ObjectVector = std::vector<Object*>;
		using ObjectIterator = ObjectVector::iterator;
		using MaterialVector = std::vector<Material>;

		~LoaderObj();
		bool load_file(const char* strFilename, bool convertQuadsToTris = true);
		void clear();

		inline ObjectIterator objects_begin()	{return m_vObjects.begin();}
		inline ObjectIterator objects_end()		{return m_vObjects.end();}
		inline int num_objects() const	{return m_vObjects.size();}
		inline const Object* get_object(int index) const {return m_vObjects[index];};

		inline int num_materials() const {return m_vMaterials.size();}
		inline const Material& get_material(int index) const {return m_vMaterials[index];}

		inline std::vector<ug::vector3>::const_iterator points_begin()	{return m_vPoints.begin();}
		inline std::vector<ug::vector3>::const_iterator points_end()	{return m_vPoints.end();}
		inline const ug::vector3* point(int index) const				{return &m_vPoints[index];}
		inline int num_points()										{return m_vPoints.size();}

		inline std::vector<ug::vector2>::const_iterator uv_begin()	{return m_vTexCoords.begin();}
		inline std::vector<ug::vector2>::const_iterator uv_end()	{return m_vTexCoords.end();}
		inline const ug::vector2* uv(int index) const				{return &m_vTexCoords[index];}
		inline int num_uvs()										{return m_vTexCoords.size();}

	protected:
		int get_material_index_by_name(const char* name) const;

		ObjectVector	m_vObjects;
		MaterialVector	m_vMaterials;

		std::vector<ug::vector3>	m_vPoints;
		std::vector<ug::vector2>	m_vTexCoords;
};

// end group ugbase_common_util
/// \}

}

#endif
