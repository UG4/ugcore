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
////// FE1
//////////////////////

template <typename TElem, typename TPosition>
FE1Discretization<TElem, TPosition>::FE1Discretization()
{
	InitializeIntegrationPoints();
}

template <typename TElem, typename TPosition>
void FE1Discretization<TElem, TPosition>::prepareElement(TElem* elem, NumericalSolution& u, int nr_func, SubsetHandler& sh, int SubsetIndex)
{
	NumericalSolutionPattern* pattern = u.get_pattern();
	m_TrialSpace = &u.TrialSpace<TElem>(nr_func, SubsetIndex);

	Grid* grid = sh.get_assigned_grid();
	Grid::VertexAttachmentAccessor<TPosition> m_aaPos(*grid, aPosition);

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
		m_RefElem.mapLocalToGlobal(m_corners, m_LocalIP[ip], m_GlobalIP[ip]);
		m_RefElem.Trafo(m_corners, m_LocalIP[ip], m_IPTrafo[ip], m_det[ip]);

		// and Shape and ShapeGrad of Solution u at those points
		for(int i = 0; i< m_nsh; i++)
		{
			m_TrialSpace->evaluateShapeGrad(i, m_LocalIP[ip], m_LocShapeGrad[ip][i]);
			m_TrialSpace->evaluateShape(i, m_LocalIP[ip], m_Shape[ip][i]);
			MatVecMult(m_GlobShapeGrad[ip][i], m_IPTrafo[ip], m_LocShapeGrad[ip][i]);
		}
	}
}

template <typename TElem, typename TPosition>
void FE1Discretization<TElem, TPosition>::add_op_to_local_jacobian(TimeOperator* op, TElem* elem, NumericalSolution& u)
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

template <typename TElem, typename TPosition>
void FE1Discretization<TElem, TPosition>::add_op_to_local_jacobian(ScalarDifferentialOperator* op, TElem* elem, NumericalSolution& u)
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

template <typename TElem, typename TPosition>
void FE1Discretization<TElem, TPosition>::add_op_to_local_jacobian(DivergenzDifferentialOperator* op, TElem* elem, NumericalSolution& u)
{
	MathVector<2> JacobianArray[m_nsh];


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

template <typename TElem, typename TPosition>
void FE1Discretization<TElem, TPosition>::add_rhs_to_local_defect(RHS* rhs, TElem* elem)
{
	double defect;

	// compute shape fkt, grads and Trafos at ip
	for(int ip=0; ip<m_nip; ip++)
	{
		// evalutate operator
		rhs->compute_defect_at_ip(m_LocalIP[ip], m_GlobalIP[ip], defect);

		// add entry to LocalJacobian
		for(int i=0; i< m_nsh; i++)
		{
			m_LocalDefect[i] -=  m_Shape[ip][i] * defect * m_det[ip] * m_RefElem.size();
		}
	}
}

template <typename TElem, typename TPosition>
void FE1Discretization<TElem, TPosition>::add_op_to_local_defect(TimeOperator* op, TElem* elem, NumericalSolution& u)
{
	MathVector<2> uGrad;
	number uShape;

	number DefectArray;
	MathVector<2> help;

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

template <typename TElem, typename TPosition>
void FE1Discretization<TElem, TPosition>::add_op_to_local_defect(ScalarDifferentialOperator* op, TElem* elem, NumericalSolution& u)
{
	MathVector<2> uGrad;
	number uShape;

	number DefectArray;
	MathVector<2> help;

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

template <typename TElem, typename TPosition>
void FE1Discretization<TElem, TPosition>::add_op_to_local_defect(DivergenzDifferentialOperator* op, TElem* elem, NumericalSolution& u)
{
	MathVector<2> uGrad;
	number uShape;

	MathVector<2> DefectArray;
	MathVector<2> help;

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

template <typename TElem, typename TPosition>
void FE1Discretization<TElem, TPosition>::reset_local_jacobian()
{
	MatSet(m_LocalJacobian, 0.0);
}

template <typename TElem, typename TPosition>
void FE1Discretization<TElem, TPosition>::reset_local_defect()
{
	VecSet(m_LocalDefect, 0.0);
}

template <typename TElem, typename TPosition>
bool FE1Discretization<TElem, TPosition>::send_local_jacobian_to_global_jacobian(TElem* tr, Matrix& mat)
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


template <typename TElem, typename TPosition>
bool FE1Discretization<TElem, TPosition>::send_local_jacobian_to_global_jacobian(TElem* tr, Matrix& mat, number scale)
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

template <typename TElem, typename TPosition>
bool FE1Discretization<TElem, TPosition>::send_local_defect_to_global_defect(TElem* tr, Vector& vec)
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

template <typename TElem, typename TPosition>
bool FE1Discretization<TElem, TPosition>::send_local_defect_to_global_defect(TElem* tr, Vector& vec, number scale)
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

template <typename TElem, typename TPosition>
FE1Discretization<TElem, TPosition>::~FE1Discretization()
{
}

template <typename TElem, typename TPosition>
void FE1Discretization<TElem, TPosition>::InitializeIntegrationPoints()
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

template <typename TElem, typename TPosition>
FVE1lumpDiscretization<TElem, TPosition>::FVE1lumpDiscretization()
{
	InitializeIntegrationPoints();
}

template <typename TElem, typename TPosition>
void FVE1lumpDiscretization<TElem, TPosition>::prepareElement(TElem* elem, NumericalSolution& u, int nr_func, SubsetHandler& sh, int SubsetIndex)
{
	NumericalSolutionPattern* pattern = u.get_pattern();
	m_TrialSpace = &u.TrialSpace<TElem>(nr_func, SubsetIndex);

	Grid* grid = sh.get_assigned_grid();
	Grid::VertexAttachmentAccessor<TPosition> m_aaPos(*grid, aPosition);

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
		m_RefElem.mapLocalToGlobal(m_corners, m_LocalIP[ip], m_GlobalIP[ip]);
		m_RefElem.Trafo(m_corners, m_LocalIP[ip], m_IPTrafo[ip], m_det[ip]);

		// and Shape and ShapeGrad of Solution u at those points
		for(int i = 0; i< m_nsh; i++)
		{
			m_TrialSpace->evaluateShapeGrad(i, m_LocalIP[ip], m_LocShapeGrad[ip][i]);
			m_TrialSpace->evaluateShape(i, m_LocalIP[ip], m_Shape[ip][i]);
			MatVecMult(m_GlobShapeGrad[ip][i], m_IPTrafo[ip], m_LocShapeGrad[ip][i]);
		}
	}
}

template <typename TElem, typename TPosition>
void FVE1lumpDiscretization<TElem, TPosition>::add_op_to_local_jacobian(TimeOperator* op, TElem* elem, NumericalSolution& u)
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

template <typename TElem, typename TPosition>
void FVE1lumpDiscretization<TElem, TPosition>::add_op_to_local_jacobian(ScalarDifferentialOperator* op, TElem* elem, NumericalSolution& u)
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

template <typename TElem, typename TPosition>
void FVE1lumpDiscretization<TElem, TPosition>::add_op_to_local_jacobian(DivergenzDifferentialOperator* op, TElem* elem, NumericalSolution& u)
{
	MathVector<2> JacobianArray[m_nsh];


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

template <typename TElem, typename TPosition>
void FVE1lumpDiscretization<TElem, TPosition>::add_rhs_to_local_defect(RHS* rhs, TElem* elem)
{
	double defect;

	// compute shape fkt, grads and Trafos at ip
	for(int ip=0; ip<m_nip; ip++)
	{
		// evalutate operator
		rhs->compute_defect_at_ip(m_LocalIP[ip], m_GlobalIP[ip], defect);

		// add entry to LocalJacobian
		for(int i=0; i< m_nsh; i++)
		{
			m_LocalDefect[i] -=  m_Shape[ip][i] * defect * m_det[ip] * m_RefElem.size();
		}
	}
}

template <typename TElem, typename TPosition>
void FVE1lumpDiscretization<TElem, TPosition>::add_op_to_local_defect(TimeOperator* op, TElem* elem, NumericalSolution& u)
{
	MathVector<2> uGrad;
	number uShape;

	number DefectArray;
	MathVector<2> help;

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

template <typename TElem, typename TPosition>
void FVE1lumpDiscretization<TElem, TPosition>::add_op_to_local_defect(ScalarDifferentialOperator* op, TElem* elem, NumericalSolution& u)
{
	MathVector<2> uGrad;
	number uShape;

	number DefectArray;
	MathVector<2> help;

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

template <typename TElem, typename TPosition>
void FVE1lumpDiscretization<TElem, TPosition>::add_op_to_local_defect(DivergenzDifferentialOperator* op, TElem* elem, NumericalSolution& u)
{
	MathVector<2> uGrad;
	number uShape;

	MathVector<2> DefectArray;
	MathVector<2> help;

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

template <typename TElem, typename TPosition>
void FVE1lumpDiscretization<TElem, TPosition>::reset_local_jacobian()
{
	MatSet(m_LocalJacobian, 0.0);
}

template <typename TElem, typename TPosition>
void FVE1lumpDiscretization<TElem, TPosition>::reset_local_defect()
{
	VecSet(m_LocalDefect, 0.0);
}

template <typename TElem, typename TPosition>
bool FVE1lumpDiscretization<TElem, TPosition>::send_local_jacobian_to_global_jacobian(TElem* tr, Matrix& mat)
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


template <typename TElem, typename TPosition>
bool FVE1lumpDiscretization<TElem, TPosition>::send_local_jacobian_to_global_jacobian(TElem* tr, Matrix& mat, number scale)
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

template <typename TElem, typename TPosition>
bool FVE1lumpDiscretization<TElem, TPosition>::send_local_defect_to_global_defect(TElem* tr, Vector& vec)
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

template <typename TElem, typename TPosition>
bool FVE1lumpDiscretization<TElem, TPosition>::send_local_defect_to_global_defect(TElem* tr, Vector& vec, number scale)
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

template <typename TElem, typename TPosition>
FVE1lumpDiscretization<TElem, TPosition>::~FVE1lumpDiscretization()
{
}

template <typename TElem, typename TPosition>
void FVE1lumpDiscretization<TElem, TPosition>::InitializeIntegrationPoints()
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
