/*
 * discretizationscheme_impl.h
 *
 *  Created on: 28.10.2009
 *      Author: andreasvogel
 */

#ifndef __H__LIBDISCRETIZATION__DISCRETIZATIONSCHEME_IMPL__
#define __H__LIBDISCRETIZATION__DISCRETIZATIONSCHEME_IMPL__

namespace ug{

//////////////////////
////// Dirichlet
//////////////////////

template <int d>
bool DirichletValues<d>::add_dirichlet_nodes(NumericalSolution<d>& u, int nr_func, DirichletBNDCond<d>* dirichbnd, SubsetHandler& sh, uint subsetIndex)
{
	VertexBaseIterator iterBegin, iterEnd, iter;
	iterBegin = sh.begin<VertexBase>(subsetIndex);
	iterEnd = sh.end<VertexBase>(subsetIndex);

	Grid* grid = sh.get_assigned_grid();
	Grid::VertexAttachmentAccessor<Attachment<MathVector<d> > > aaPos(*grid, *(u.get_domain()->get_position_attachment()));

	int index;
	number val;
	MathVector<d> corner;

	/* loop over all Vertices */
	int nvertices = 0, ndirnodes = 0;
	for(iter = iterBegin; iter != iterEnd; iter++)
	{
		VertexBase *vert = *iter;
		corner = aaPos[vert];
		nvertices++;

		if(IsBoundaryVertex2D(*grid, vert))
		{
			index = (int) u.get_pattern()->get_index(vert, nr_func);
			dirichbnd->BNDValueFunction(corner, val);
			m_vector_values.push_back(val);
			m_vector_indices.push_back(index);
			m_matrixrow_indices.push_back(index);
			ndirnodes++;
		}
	}
	//std::cout << nvertices << " Vertices checked: Nr. Dirichlet nodes found:" << ndirnodes << "(now total: "<< m_vector_indices.size() <<")"<< std::endl;

	return true;
}

template <int d>
bool DirichletValues<d>::set_values(Vector& vec)
{
	double* valueArray = new double[m_vector_indices.size()];
	int* indexArray = new int[m_vector_indices.size()];

	for(uint i=0; i<m_vector_indices.size(); i++)
	{
		valueArray[i] = m_vector_values[i];
		indexArray[i] = m_vector_indices[i];
	}

	if(vec.set_values(m_vector_indices.size(), indexArray, valueArray) != true)
		return false;
	delete valueArray;
	delete indexArray;

	std::cout << m_vector_indices.size() << " Boundary nodes set to Dirichlet value" << std::endl;
	return true;
}

template <int d>
bool DirichletValues<d>::set_zero_values(Vector& vec)
{
	double* valueArray = new double[m_vector_indices.size()];
	int* indexArray = new int[m_vector_indices.size()];

	for(uint i=0; i<m_vector_indices.size(); i++)
	{
		valueArray[i] = 0.0;
		indexArray[i] = m_vector_indices[i];
	}

	if(vec.set_values(m_vector_indices.size(), indexArray, valueArray) != true)
		return false;
	delete valueArray;
	delete indexArray;

	std::cout << m_vector_indices.size() << " Boundary nodes set to zero value" << std::endl;
	return true;
}


template <int d>
bool DirichletValues<d>::set_rows(Matrix& mat)
{
	int* indexArray = new int[m_matrixrow_indices.size()];

	for(uint i=0; i<m_matrixrow_indices.size(); i++)
	{
		indexArray[i] = m_matrixrow_indices[i];
	}

	if(mat.set_dirichletrows(m_matrixrow_indices.size(), indexArray) != true)
		return false;
	delete indexArray;

	std::cout << m_matrixrow_indices.size() << " Matrix rows set to Identity row" << std::endl;
	return true;
}


//////////////////////
////// FE1
//////////////////////

template <typename TElem, int d>
FE1Discretization<TElem, d>::FE1Discretization()
{
	InitializeIntegrationPoints();
}

template <typename TElem, int d>
void FE1Discretization<TElem, d>::prepareElement(TElem* elem, typename ug::NumericalSolution<d>& u, int nr_func, SubsetHandler& sh, int SubsetIndex)
{
	DoFPattern* pattern = u.get_pattern();
	m_TrialSpace = &(u.template get_TrialSpace<TElem>(nr_func));

	Grid* grid = sh.get_assigned_grid();
	Grid::VertexAttachmentAccessor<Attachment<MathVector<d> > > m_aaPos(*grid, *(u.get_domain()->get_position_attachment()));

	u.get_local_DoFValues(elem, nr_func, m_loc_u);

	// get corner coordinates and global DoF numbers
	for(int i=0; i< m_nsh; i++)
	{
		VertexBase* vert = elem->vertex(i);
		m_corners[i] = m_aaPos[vert];
		m_rows[i] = (int) (pattern->get_index(vert, nr_func));
	}

	// get global ip's ...
	for(int ip=0; ip<m_nip; ip++)
	{
		m_RefElem.template mapLocalToGlobal<d>(m_corners, m_LocalIP[ip], m_GlobalIP[ip]);
		m_RefElem.Trafo(m_corners, m_LocalIP[ip], m_IPTrafo[ip], m_det[ip]);

		// and Shape and ShapeGrad of Solution u at those points
		for(int i = 0; i< m_nsh; i++)
		{
			// TODO: This does not work anymore
			//m_TrialSpace->evaluateShapeGrad(i, m_LocalIP[ip], m_LocShapeGrad[ip][i]);
			//m_TrialSpace->evaluateShape(i, m_LocalIP[ip], m_Shape[ip][i]);
			MatVecMult(m_GlobShapeGrad[ip][i], m_IPTrafo[ip], m_LocShapeGrad[ip][i]);
		}
	}
}

template <typename TElem, int d>
void FE1Discretization<TElem, d>::add_op_to_local_jacobian(TimeOperator<d>* op, TElem* elem, NumericalSolution<d>& u)
{
	number JacobianArray[m_nsh];

	// compute shape fkt, grads and Trafos at ip
	for(int ip=0; ip<m_nip; ip++)
	{
		// evalutate operator
		op->compute_jacobian_at_ip(m_GlobalIP[ip], m_Shape[ip], m_GlobShapeGrad[ip], JacobianArray, m_nsh);

		// add entry to LocalJacobian
		for(int i=0; i< m_nsh; i++)
			for(int j=0; j< m_nsh; j++)
			{
				m_LocalJacobian[i][j] +=  m_Shape[ip][i] * JacobianArray[j] * m_det[ip] *m_RefElem.size();
			}
	}
}

template <typename TElem, int d>
void FE1Discretization<TElem, d>::add_op_to_local_jacobian(ScalarDifferentialOperator<d>* op, TElem* elem, NumericalSolution<d>& u)
{
	number JacobianArray[m_nsh];

	// compute shape fkt, grads and Trafos at ip
	for(int ip=0; ip<m_nip; ip++)
	{
		// evalutate operator
		op->compute_jacobian_at_ip(m_GlobalIP[ip], m_Shape[ip], m_GlobShapeGrad[ip], JacobianArray, m_nsh);

		// add entry to LocalJacobian
		for(int i=0; i< m_nsh; i++)
			for(int j=0; j< m_nsh; j++)
			{
				m_LocalJacobian[i][j] +=  m_Shape[ip][i] * JacobianArray[j] * m_det[ip] * m_RefElem.size();
			}
	}
}

template <typename TElem, int d>
void FE1Discretization<TElem, d>::add_op_to_local_jacobian(DivergenzDifferentialOperator<d>* op, TElem* elem, NumericalSolution<d>& u)
{
	MathVector<d> JacobianArray[m_nsh];


	// compute shape fkt, grads and Trafos at ip
	for(int ip=0; ip<m_nip; ip++)
	{
		// evalutate operator
		op->compute_jacobian_at_ip(m_GlobalIP[ip], m_Shape[ip], m_GlobShapeGrad[ip], JacobianArray, m_nsh);

		// add entry to LocalJacobian
		for(int i=0; i< m_nsh; i++)
			for(int j=0; j< m_nsh; j++)
			{
				m_LocalJacobian[i][j] +=  VecDot(m_GlobShapeGrad[ip][i], JacobianArray[j]) * m_det[ip] * m_RefElem.size();
			}
	}
}

template <typename TElem, int d>
void FE1Discretization<TElem, d>::add_rhs_to_local_defect(RHS<d>* rhs, TElem* elem)
{
	double defect;

	// compute shape fkt, grads and Trafos at ip
	for(int ip=0; ip<m_nip; ip++)
	{
		// evalutate operator
		rhs->compute_defect_at_ip(m_GlobalIP[ip], defect);

		// add entry to LocalJacobian
		for(int i=0; i< m_nsh; i++)
		{
			m_LocalDefect[i] -=  m_Shape[ip][i] * defect * m_det[ip] * m_RefElem.size();
		}
	}
}

template <typename TElem, int d>
void FE1Discretization<TElem, d>::add_op_to_local_defect(TimeOperator<d>* op, TElem* elem, NumericalSolution<d>& u)
{
	MathVector<d> uGrad;
	number uShape;

	number DefectArray;
	MathVector<d> help;

	// compute shape fkt, grads and Trafos at ip
	for(int ip=0; ip<m_nip; ip++)
	{
		uShape = 0.;
		VecSet(uGrad, 0.0);
		for(int i = 0; i< m_nsh; i++)
		{
			uShape += m_loc_u[i] * m_Shape[ip][i];
			VecScale(help, m_GlobShapeGrad[ip][i], m_loc_u[i]);
			uGrad +=  help;
		}

		// evalutate operator
		op->compute_defect_at_ip(m_GlobalIP[ip], uShape, uGrad, DefectArray);

		// add entry to LocalJacobian
		for(int i=0; i< m_nsh; i++)
		{
			m_LocalDefect[i] +=  m_Shape[ip][i] * DefectArray * m_det[ip] * m_RefElem.size();
		}
	}
}

template <typename TElem, int d>
void FE1Discretization<TElem, d>::add_op_to_local_defect(ScalarDifferentialOperator<d>* op, TElem* elem, NumericalSolution<d>& u)
{
	MathVector<d> uGrad;
	number uShape;

	number DefectArray;
	MathVector<d> help;

	// compute shape fkt, grads and Trafos at ip
	for(int ip=0; ip<m_nip; ip++)
	{
		uShape = 0.;
		VecSet(uGrad, 0.0);
		for(int i = 0; i< m_nsh; i++)
		{
			uShape += m_loc_u[i] * m_Shape[ip][i];
			VecScale(help, m_GlobShapeGrad[ip][i], m_loc_u[i]);
			uGrad +=  help;
		}

		// evalutate operator
		op->compute_defect_at_ip(m_GlobalIP[ip], uShape, uGrad, DefectArray);

		// add entry to LocalJacobian
		for(int i=0; i< m_nsh; i++)
		{
			m_LocalDefect[i] +=  m_Shape[ip][i] * DefectArray * m_det[ip] * m_RefElem.size();
		}
	}
}

template <typename TElem, int d>
void FE1Discretization<TElem, d>::add_op_to_local_defect(DivergenzDifferentialOperator<d>* op, TElem* elem, NumericalSolution<d>& u)
{
	MathVector<d> uGrad;
	number uShape;

	MathVector<d> DefectArray;
	MathVector<d> help;

	// compute shape fkt, grads and Trafos at ip
	for(int ip=0; ip<m_nip; ip++)
	{
		uShape = 0.;
		VecSet(uGrad, 0.0);
		for(int i = 0; i< m_nsh; i++)
		{
			uShape += m_loc_u[i] * m_Shape[ip][i];
			VecScale(help, m_GlobShapeGrad[ip][i], m_loc_u[i]);
			uGrad +=  help;
		}

		// evalutate operator
		op->compute_defect_at_ip(m_GlobalIP[ip], uShape, uGrad, DefectArray);

		// add entry to LocalJacobian
		for(int i=0; i< m_nsh; i++)
		{
			m_LocalDefect[i] +=  VecDot(m_GlobShapeGrad[ip][i], DefectArray) * m_det[ip] * m_RefElem.size();
		}
	}
}

template <typename TElem, int d>
void FE1Discretization<TElem, d>::reset_local_jacobian()
{
	MatSet(m_LocalJacobian, 0.0);
}

template <typename TElem, int d>
void FE1Discretization<TElem, d>::reset_local_defect()
{
	VecSet(m_LocalDefect, 0.0);
}

template <typename TElem, int d>
bool FE1Discretization<TElem, d>::send_local_jacobian_to_global_jacobian(TElem* tr, Matrix& mat)
{
	int nrows = m_nsh;
	int ncols[m_nsh];
	for(int i=0; i< m_nsh; ++i)
	{
		ncols[i] = m_nsh;
	}
	int cols[m_nsh*m_nsh];
	double values[m_nsh*m_nsh];

	for(int i=0; i<m_nsh; i++)
		for(int j=0; j<m_nsh; j++)
		{
			cols[m_nsh*i + j] = m_rows[j];
			values[m_nsh*i + j] = m_LocalJacobian[i][j];
		}

	if(mat.add_values(nrows, ncols, m_rows, cols, values) == true) return true;

	return false;
}


template <typename TElem, int d>
bool FE1Discretization<TElem, d>::send_local_jacobian_to_global_jacobian(TElem* tr, Matrix& mat, number scale)
{
	int nrows = m_nsh;
	int ncols[m_nsh];
	for(int i=0; i< m_nsh; ++i)
	{
		ncols[i] = m_nsh;
	}
	int cols[m_nsh*m_nsh];
	double values[m_nsh*m_nsh];

	for(int i=0; i<m_nsh; i++)
		for(int j=0; j<m_nsh; j++)
		{
			cols[m_nsh*i + j] = m_rows[j];
			values[m_nsh*i + j] = scale*m_LocalJacobian[i][j];
		}

	if(mat.add_values(nrows, ncols, m_rows, cols, values) == true) return true;

	return false;
}

template <typename TElem, int d>
bool FE1Discretization<TElem, d>::send_local_defect_to_global_defect(TElem* tr, Vector& vec)
{
	int nentries = m_nsh;
	double values[m_nsh];

	for(int i=0; i< m_nsh; i++)
	{
		values[i] = m_LocalDefect[i];
	}

	if(vec.add_values(nentries, m_rows, values) == true) return true;

	return false;
}

template <typename TElem, int d>
bool FE1Discretization<TElem, d>::send_local_defect_to_global_defect(TElem* tr, Vector& vec, number scale)
{
	int nentries = m_nsh;
	double values[m_nsh];

	for(int i=0; i< m_nsh; i++)
	{
		values[i] = scale * m_LocalDefect[i];
	}

	if(vec.add_values(nentries, m_rows, values) == true) return true;

	return false;
}

template <typename TElem, int d>
FE1Discretization<TElem, d>::~FE1Discretization()
{
}

template <typename TElem, int d>
void FE1Discretization<TElem, d>::InitializeIntegrationPoints()
{
	// define integration points
	//m_nip = 1;

	if(reference_element_traits<TElem>::NumberCorners == 3)
	{
		m_LocalIP = new vector2[1];
		m_LocalIP[0].x = 1./3.;
		m_LocalIP[0].y = 1./3.;
	}
	else if(reference_element_traits<TElem>::NumberCorners == 4)
	{
		m_LocalIP = new vector2[1];
		m_LocalIP[0].x = .5;
		m_LocalIP[0].y = .5;
	}
	else assert(0 && "ERROR in InitializeIntegrationPoints: ReferenceElement not detected");

	//m_GlobalIP = new vector3[1];
}

//////////////////////
////// FVE1lump
//////////////////////

template <typename TElem, int d>
FVE1lumpDiscretization<TElem, d>::FVE1lumpDiscretization()
{
	InitializeIntegrationPoints();
}

template <typename TElem, int d>
void FVE1lumpDiscretization<TElem, d>::prepareElement(TElem* elem, NumericalSolution<d>& u, int nr_func, SubsetHandler& sh, int SubsetIndex)
{
	DoFPattern* pattern = u.get_pattern();
	m_TrialSpace = &(u.template get_TrialSpace<TElem>(nr_func));

	Grid* grid = sh.get_assigned_grid();
	Grid::VertexAttachmentAccessor<Attachment<MathVector<d> > > m_aaPos(*grid, *(u.get_domain()->get_position_attachment()));

	// get corner coordinates and global DoF numbers
	for(int i=0; i< m_nsh; i++)
	{
		VertexBase* vert = elem->vertex(i);
		m_corners[i] = m_aaPos[vert];
		m_rows[i] = (int) (pattern->get_index(vert, nr_func));
	}

	// get global ip's ...
	for(int ip=0; ip<m_nip; ip++)
	{
		//TODO: This does not work anymore
		//m_RefElem.template mapLocalToGlobal<d>(m_corners, m_LocalIP[ip], m_GlobalIP[ip]);
		//m_RefElem.Trafo(m_corners, m_LocalIP[ip], m_IPTrafo[ip], m_det[ip]);

		// and Shape and ShapeGrad of Solution u at those points
		for(int i = 0; i< m_nsh; i++)
		{
			m_TrialSpace->evaluateShapeGrad(i, m_LocalIP[ip], m_LocShapeGrad[ip][i]);
			m_TrialSpace->evaluateShape(i, m_LocalIP[ip], m_Shape[ip][i]);
			MatVecMult(m_GlobShapeGrad[ip][i], m_IPTrafo[ip], m_LocShapeGrad[ip][i]);
		}
	}
}

template <typename TElem, int d>
void FVE1lumpDiscretization<TElem, d>::add_op_to_local_jacobian(TimeOperator<d>* op, TElem* elem, NumericalSolution<d>& u)
{
	number JacobianArray[m_nsh];

	// compute shape fkt, grads and Trafos at ip
	for(int ip=0; ip<m_nip; ip++)
	{
		// evalutate operator
		op->compute_jacobian_at_ip(m_GlobalIP[ip], m_Shape[ip], m_GlobShapeGrad[ip], JacobianArray, m_nsh);

		// add entry to LocalJacobian
		for(int i=0; i< m_nsh; i++)
			for(int j=0; j< m_nsh; j++)
			{
				m_LocalJacobian[i][j] +=  m_Shape[ip][i] * JacobianArray[j] * m_det[ip] *m_RefElem.size();
			}
	}
}

template <typename TElem, int d>
void FVE1lumpDiscretization<TElem, d>::add_op_to_local_jacobian(ScalarDifferentialOperator<d>* op, TElem* elem, NumericalSolution<d>& u)
{
	number JacobianArray[m_nsh];

	// compute shape fkt, grads and Trafos at ip
	for(int ip=0; ip<m_nip; ip++)
	{
		// evalutate operator
		op->compute_jacobian_at_ip(m_GlobalIP[ip], m_Shape[ip], m_GlobShapeGrad[ip], JacobianArray, m_nsh);

		// add entry to LocalJacobian
		for(int i=0; i< m_nsh; i++)
			for(int j=0; j< m_nsh; j++)
			{
				m_LocalJacobian[i][j] +=  m_Shape[ip][i] * JacobianArray[j] * m_det[ip] * m_RefElem.size();
			}
	}
}

template <typename TElem, int d>
void FVE1lumpDiscretization<TElem, d>::add_op_to_local_jacobian(DivergenzDifferentialOperator<d>* op, TElem* elem, NumericalSolution<d>& u)
{
	MathVector<d> JacobianArray[m_nsh];


	// compute shape fkt, grads and Trafos at ip
	for(int ip=0; ip<m_nip; ip++)
	{
		// evalutate operator
		op->compute_jacobian_at_ip(m_GlobalIP[ip], m_Shape[ip], m_GlobShapeGrad[ip], JacobianArray, m_nsh);

		// add entry to LocalJacobian
		for(int i=0; i< m_nsh; i++)
			for(int j=0; j< m_nsh; j++)
			{
				m_LocalJacobian[i][j] +=  VecDot(m_GlobShapeGrad[ip][i], JacobianArray[j]) * m_det[ip] * m_RefElem.size();
			}
	}
}

template <typename TElem, int d>
void FVE1lumpDiscretization<TElem, d>::add_rhs_to_local_defect(RHS<d>* rhs, TElem* elem)
{
	double defect;

	// compute shape fkt, grads and Trafos at ip
	for(int ip=0; ip<m_nip; ip++)
	{
		// evalutate operator
		rhs->compute_defect_at_ip(m_GlobalIP[ip], defect);

		// add entry to LocalJacobian
		for(int i=0; i< m_nsh; i++)
		{
			m_LocalDefect[i] -=  m_Shape[ip][i] * defect * m_det[ip] * m_RefElem.size();
		}
	}
}

template <typename TElem, int d>
void FVE1lumpDiscretization<TElem, d>::add_op_to_local_defect(TimeOperator<d>* op, TElem* elem, NumericalSolution<d>& u)
{
	MathVector<d> uGrad;
	number uShape;

	number DefectArray;
	MathVector<d> help;

	// compute shape fkt, grads and Trafos at ip
	for(int ip=0; ip<m_nip; ip++)
	{
		uShape = 0.;
		VecSet(uGrad, 0.0);
		for(int i = 0; i< m_nsh; i++)
		{
			uShape += m_loc_u[i] * m_Shape[ip][i];
			VecScale(help, m_GlobShapeGrad[ip][i], m_loc_u[i]);
			uGrad +=  help;
		}

		// evalutate operator
		op->compute_defect_at_ip(m_GlobalIP[ip], uShape, uGrad, DefectArray);

		// add entry to LocalJacobian
		for(int i=0; i< m_nsh; i++)
		{
			m_LocalDefect[i] +=  m_Shape[ip][i] * DefectArray * m_det[ip] * m_RefElem.size();
		}
	}
}

template <typename TElem, int d>
void FVE1lumpDiscretization<TElem, d>::add_op_to_local_defect(ScalarDifferentialOperator<d>* op, TElem* elem, NumericalSolution<d>& u)
{
	MathVector<d> uGrad;
	number uShape;

	number DefectArray;
	MathVector<d> help;

	// compute shape fkt, grads and Trafos at ip
	for(int ip=0; ip<m_nip; ip++)
	{
		uShape = 0.;
		VecSet(uGrad, 0.0);
		for(int i = 0; i< m_nsh; i++)
		{
			uShape += m_loc_u[i] * m_Shape[ip][i];
			VecScale(help, m_GlobShapeGrad[ip][i], m_loc_u[i]);
			uGrad +=  help;
		}

		// evalutate operator
		op->compute_defect_at_ip(m_GlobalIP[ip], uShape, uGrad, DefectArray);

		// add entry to LocalJacobian
		for(int i=0; i< m_nsh; i++)
		{
			m_LocalDefect[i] +=  m_Shape[ip][i] * DefectArray * m_det[ip] * m_RefElem.size();
		}
	}
}

template <typename TElem, int d>
void FVE1lumpDiscretization<TElem, d>::add_op_to_local_defect(DivergenzDifferentialOperator<d>* op, TElem* elem, NumericalSolution<d>& u)
{
	MathVector<d> uGrad;
	number uShape;

	MathVector<d> DefectArray;
	MathVector<d> help;

	// compute shape fkt, grads and Trafos at ip
	for(int ip=0; ip<m_nip; ip++)
	{
		uShape = 0.;
		VecSet(uGrad, 0.0);
		for(int i = 0; i< m_nsh; i++)
		{
			uShape += m_loc_u[i] * m_Shape[ip][i];
			VecScale(help, m_GlobShapeGrad[ip][i], m_loc_u[i]);
			uGrad +=  help;
		}

		// evalutate operator
		op->compute_defect_at_ip(m_GlobalIP[ip], uShape, uGrad, DefectArray);

		// add entry to LocalJacobian
		for(int i=0; i< m_nsh; i++)
		{
			m_LocalDefect[i] +=  VecDot(m_GlobShapeGrad[ip][i], DefectArray) * m_det[ip] * m_RefElem.size();
		}
	}
}

template <typename TElem, int d>
void FVE1lumpDiscretization<TElem, d>::reset_local_jacobian()
{
	MatSet(m_LocalJacobian, 0.0);
}

template <typename TElem, int d>
void FVE1lumpDiscretization<TElem, d>::reset_local_defect()
{
	VecSet(m_LocalDefect, 0.0);
}

template <typename TElem, int d>
bool FVE1lumpDiscretization<TElem, d>::send_local_jacobian_to_global_jacobian(TElem* tr, Matrix& mat)
{
	int nrows = m_nsh;
	int ncols[m_nsh];
	for(int i=0; i< m_nsh; ++i)
	{
		ncols[i] = m_nsh;
	}
	int cols[m_nsh*m_nsh];
	double values[m_nsh*m_nsh];

	for(int i=0; i<m_nsh; i++)
		for(int j=0; j<m_nsh; j++)
		{
			cols[m_nsh*i + j] = m_rows[j];
			values[m_nsh*i + j] = m_LocalJacobian[i][j];
		}

	if(mat.add_values(nrows, ncols, m_rows, cols, values) == true) return true;

	return false;
}


template <typename TElem, int d>
bool FVE1lumpDiscretization<TElem, d>::send_local_jacobian_to_global_jacobian(TElem* tr, Matrix& mat, number scale)
{
	int nrows = m_nsh;
	int ncols[m_nsh];
	for(int i=0; i< m_nsh; ++i)
	{
		ncols[i] = m_nsh;
	}
	int cols[m_nsh*m_nsh];
	double values[m_nsh*m_nsh];

	for(int i=0; i<m_nsh; i++)
		for(int j=0; j<m_nsh; j++)
		{
			cols[m_nsh*i + j] = m_rows[j];
			values[m_nsh*i + j] = scale*m_LocalJacobian[i][j];
		}

	if(mat.add_values(nrows, ncols, m_rows, cols, values) == true) return true;

	return false;
}

template <typename TElem, int d>
bool FVE1lumpDiscretization<TElem, d>::send_local_defect_to_global_defect(TElem* tr, Vector& vec)
{
	int nentries = m_nsh;
	double values[m_nsh];

	for(int i=0; i< m_nsh; i++)
	{
		values[i] = m_LocalDefect[i];
	}

	if(vec.add_values(nentries, m_rows, values) == true) return true;

	return false;
}

template <typename TElem, int d>
bool FVE1lumpDiscretization<TElem, d>::send_local_defect_to_global_defect(TElem* tr, Vector& vec, number scale)
{
	int nentries = m_nsh;
	double values[m_nsh];

	for(int i=0; i< m_nsh; i++)
	{
		values[i] = scale * m_LocalDefect[i];
	}

	if(vec.add_values(nentries, m_rows, values) == true) return true;

	return false;
}

template <typename TElem, int d>
FVE1lumpDiscretization<TElem, d>::~FVE1lumpDiscretization()
{
}

template <typename TElem, int d>
void FVE1lumpDiscretization<TElem, d>::InitializeIntegrationPoints()
{
	if(reference_element_traits<TElem>::NumberCorners == 3)
	{
		m_LocalIP = new vector2[1];
		m_LocalIP[0].x = 1./3.;
		m_LocalIP[0].y = 1./3.;
	}
	else if(reference_element_traits<TElem>::NumberCorners == 4)
	{
		m_LocalIP = new vector2[1];
		m_LocalIP[0].x = .5;
		m_LocalIP[0].y = .5;
	}
	else assert(0 && "ERROR in InitializeIntegrationPoints: ReferenceElement not detected");

}


}


#endif /* __H__LIBDISCRETIZATION__DISCRETIZATIONSCHEME_IMPL__ */
