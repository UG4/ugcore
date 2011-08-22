/*
 * level_set_util.h
 *
 *  Created on: 01.07.2011
 *      Author: Christian Wehner  christian.wehner@gcsc.uni-frankfurt.de
 */

#ifndef LEVEL_SET_UTIL_IMPL_H_
#define LEVEL_SET_UTIL_IMPL_H_

#include "level_set.h"
#include "lib_discretization/local_finite_element/local_shape_function_set.h"
#include "lib_discretization/spatial_discretization/disc_util/finite_volume_geometry.h"
#include "lib_discretization/reference_element/reference_element.h"

namespace ug{

template<typename TGridFunction>
bool FV1LevelSetDisc<TGridFunction>::analytic_velocity(MathVector<dim> & v,number t,MathVector<dim> x)
{
	switch(dim)
	{
	    case 1: 
		   v[0] = 1;
		   return 0;
		case 2: 
		   v[0] = -x[1];
		   v[1] = x[0];
		  // v[0] = 1;
		  // v[1] = 0;
//		   v[0] = 1;
	//	   v[1] = 1;
		   return 0;
		case 3: 
		   v[0] = -x[1];
		   v[1] =  x[0];
		   v[2] =  0;
		   return 0;
	};		
	return true;
};

template<typename TGridFunction>
number FV1LevelSetDisc<TGridFunction>::analytic_source(number t,MathVector<dim> x)
{
	switch(dim)
	{
	    case 1: 
		   return 1;
		case 2: 
		   return 0;
		case 3: 
		   return 1;
	};		
	return 1;
};

template<typename TGridFunction>
number FV1LevelSetDisc<TGridFunction>::analytic_solution(number t,MathVector<dim> x)
{
	MathVector<dim> xnew;
    switch(dim)
	{
        case 1:  
			//number v=1;
			//xnew=x-v*t;
			//return sin(xnew);
        	return 0;
        case 2:
        	// return 0;
        	//return x[1];//-t;
          //  number z;
      //      z=min(min(abs(x[0]-1),abs(x[0]+1)),min(abs(x[1]-1),abs(x[1]+1)));
            //z=min(min(abs(x[0]),abs(1-x[0])),min(abs(x[1]),abs(1-x[1])));
        //	return z;
		    xnew[0] = x[0]*cos(t)+x[1]*sin(t);
		 	xnew[1] = -x[0]*sin(t)+x[1]*cos(t);
		 	return sqrt((xnew[0]-0.5)*(xnew[0]-0.5) + xnew[1]*xnew[1]) - 0.5;
            xnew[0] = x[0] - t;
            xnew[1] = x[1] - t;
			return x[0]-t;
			return sqrt(xnew[0]*xnew[0] + xnew[1]*xnew[1]) - 0.5;
           // return 2*xnew[0]-xnew[1];
        case 3:  
			return 0;
	}

    return -1;
};

template<typename TGridFunction>
bool FV1LevelSetDisc<TGridFunction>::fill_v_vec(TGridFunction& vel,int component)
{
    //	get domain of grid function
	domain_type& domain = vel.get_domain();

//	get grid type of domain
	typedef typename domain_type::grid_type grid_type;

//	get grid of domain
	grid_type& grid = domain.get_grid();

	typedef typename domain_type::position_accessor_type position_accessor_type;

	if (component>dim){
		throw(UGError("fill_v_vec: component > dim."));
	}

	for (int si=0;si<vel.num_subsets();++si)
	{
		for(VertexBaseConstIterator iter = vel.template begin<VertexBase>(si);
									   iter != grid.template end<VertexBase>(si); ++iter)
		{
		//	get vertex
			VertexBase* vrt = *iter;
			MathVector<dim> coord;
			MathVector<dim> vnode;
			position_accessor_type aaPos = domain.get_position_accessor();
			coord = aaPos[vrt];

			typedef typename TGridFunction::multi_index_vector_type index_type;

		//	get vector holding all indices on the vertex
			multi_index_vector_type ind;

		//	read indices on vertex
			typedef typename TGridFunction::dof_distribution_type dof_distribution_type;
			const dof_distribution_type& dd = vel.get_dof_distribution();

			dd.inner_multi_indices(vrt, 0, ind);
			analytic_velocity(vnode,m_time,coord);
			BlockRef(vel[ind[0][0]],ind[0][1]) = vnode[component];
	     }
	};
	return true;
};

template<typename TGridFunction>
bool compute_normal_velocity(TGridFunction& vx,TGridFunction& vy,const TGridFunction& u){
	return true;
}

template<typename TGridFunction>
bool FV1LevelSetDisc<TGridFunction>::limit_grad(TGridFunction& uOld, aaGrad& aaGrad)
{
	//	get domain of grid function
	domain_type& domain = uOld.get_domain();

	//	get grid type of domain
	typedef typename domain_type::grid_type grid_type;

	//	get grid of domain
	grid_type& grid = domain.get_grid();

	typedef typename domain_type::position_accessor_type position_accessor_type;

	position_accessor_type aaPos = domain.get_position_accessor();

	multi_index_vector_type ind;
	typedef typename TGridFunction::dof_distribution_type dof_distribution_type;
	const dof_distribution_type& dd = uOld.get_dof_distribution();

	//	create Attachment for scv-volume size
	ANumber aMax;
	ANumber aMin;

	//	attach to grid
	grid.attach_to_vertices(aMin);
	grid.attach_to_vertices(aMax);

	//	get attachment accessor to access values
	Grid::VertexAttachmentAccessor<ANumber> aaMin(grid, aMin);
	Grid::VertexAttachmentAccessor<ANumber> aaMax(grid, aMax);
	for (int si=0;si<uOld.num_subsets();++si){
		for(VertexBaseConstIterator iter = uOld.template begin<VertexBase>(si) ;iter !=uOld.template end<VertexBase>(si); ++iter)
				{
		    	    VertexBase* vrt = *iter;
		    	    MathVector<dim> coord;
		    	    coord = aaPos[vrt];
		    	    //	read indices on vertex
		    	    typedef typename TGridFunction::multi_index_vector_type index_type;
		      	    //	get vector holding all indices on the vertex
				    dd.inner_multi_indices(vrt, 0, ind);
				    aaMax[vrt] = BlockRef(uOld[ind[0][0]],ind[0][1]);
				    aaMin[vrt] = BlockRef(uOld[ind[0][0]],ind[0][1]);
			    }
	}
	for (int si=0;si<uOld.num_subsets();++si){
		//UG_LOG("si " << si << "\n");
		for(EdgeBaseConstIterator iter = uOld.template begin<EdgeBase>(si) ;iter !=uOld.template end<EdgeBase>(si); ++iter)
		{
			EdgeBase* edge = *iter;
			VertexBase* vi=edge->vertex(0);
			VertexBase* vj=edge->vertex(1);
			dd.inner_multi_indices(vi, 0, ind);
			number ui = BlockRef(uOld[ind[0][0]],ind[0][1]);
			dd.inner_multi_indices(vj, 0, ind);
			number uj = BlockRef(uOld[ind[0][0]],ind[0][1]);
			//UG_LOG("edge " << aaPos[vi] << "-" << aaPos[vj] << " [" << ui << " " << uj << "]\n");
			if (uj<aaMin[vi]){
				aaMin[vi]=uj;
			}
			if (uj>aaMax[vi]){
				aaMax[vi]=uj;
			}
			if (ui<aaMin[vj]){
				aaMin[vj]=ui;
			}
			if (ui>aaMax[vj]){
				aaMax[vj]=ui;
			}
		}

	};
	for (int si=0;si<2;++si){
		for(VertexBaseConstIterator iter = uOld.template begin<VertexBase>(si) ;iter !=uOld.template end<VertexBase>(si); ++iter)
		{
		    VertexBase* vrt = *iter;
			MathVector<dim> coord;
			coord = aaPos[vrt];
			// read indices on vertex
			typedef typename TGridFunction::multi_index_vector_type index_type;
			//	get vector holding all indices on the vertex
			dd.inner_multi_indices(vrt, 0, ind);
			//UG_LOG(coord << " min=" << aaMin[vrt] << " max=" << aaMax[vrt] << "\n");
		}
	}
	for (int si=0;si<1;++si){
	    for(EdgeBaseConstIterator iter = uOld.template begin<EdgeBase>(si) ;iter !=uOld.template end<EdgeBase>(si); ++iter)
   	    {
         EdgeBase* edge = *iter;
         VertexBase* vi=edge->vertex(0);
         VertexBase* vj=edge->vertex(1);
         MathVector<dim> coordi,coordj,coordij,distVec,gradi,gradj;
         gradi = aaGrad[vi];
         gradj = aaGrad[vj];
		 coordi = aaPos[vi];
		 coordj = aaPos[vj];
		 dd.inner_multi_indices(vi, 0, ind);
		 number ui = BlockRef(uOld[ind[0][0]],ind[0][1]);
		 dd.inner_multi_indices(vj, 0, ind);
		 number uj = BlockRef(uOld[ind[0][0]],ind[0][1]);
		 VecScaleAdd(coordij,0.5,coordi,0.5,coordj);
		 VecSubtract(distVec, coordij,coordi);
		 number uij = ui + distVec*gradi;
		 number alpha = 1;
		 if (uij>ui){
		     if (uij>aaMax[vi]) alpha=(aaMax[vi]-ui)/(distVec*gradi);
		     if (alpha<1){
		    	 //UG_LOG("edge " << coordi << " " << coordj << "\n");
		         //UG_LOG(coordi << " u " << ui << " uij "  << uij << " max " << aaMax[vi] << " alpha " << alpha << "\n");
			     aaGrad[vi]*=alpha;
		     };
		 }else{
		     if (uij<aaMin[vi]) alpha=(aaMin[vi]-ui)/(distVec*gradi);
		     if (alpha<1){
		    	 //UG_LOG("edge " << coordi << " " << coordj << "\n");
		    	 //UG_LOG(coordi << " u " << ui << " uij "  << uij << " min " << aaMax[vi] << " alpha " << alpha << "\n");
		    	 aaGrad[vi]*=alpha;
		     };
		 };
		 VecSubtract(distVec, coordij,coordj);
		 uij = uj + distVec*gradj;
		 alpha = 1;
		 if (uij>uj){
		    if (uij>aaMax[vj]) alpha=(aaMax[vj]-uj)/(distVec*gradj);
		    if (alpha<1){
		    	//UG_LOG("-- edge " << coordi << " " << coordj << "\n");
		    	//UG_LOG(coordj << " u " << uj << " uij "  << uij << " max " << aaMax[vj] << " alpha " << alpha << "\n");
		 	    aaGrad[vj]*=alpha;
		    };
		 }else{
		 	if (uij<aaMin[vj]) alpha=(aaMin[vj]-uj)/(distVec*gradj);
		 	if (alpha<1){
		 		//UG_LOG("-- edge " << coordi << " " << coordj << "\n");
		 		//UG_LOG(coordj << " u " << uj << " uij "  << uij << " min " << aaMax[vj] << " alpha " << alpha << "\n");
		 	    aaGrad[vj]*=alpha;
		 	};
		 };
         // UG_LOG(" coord vertex 0 " << aaPos[v0] << " coord vertex 1 " << aaPos[v1] << "\n");
	}
	//	get element iterator type
	typedef typename domain_traits<dim>::const_iterator ElemIterator;
	//	get element type
	typedef typename domain_traits<dim>::geometric_base_object ElemType;
	//	get iterators
	ElemIterator iter = uOld.template begin<ElemType>(si);
	ElemIterator iterEnd = uOld.template end<ElemType>(si);
	//	loop elements of dimension
	for(  ;iter !=iterEnd; ++iter)
	{
	    //	get Elem
		ElemType* elem = *iter;

		//	get vertices of the Elem
		std::vector<VertexBase*> vVrt;
		CollectVertices(vVrt, grid, elem);

		//	get position accessor
		typedef typename domain_type::position_accessor_type position_accessor_type;
		const position_accessor_type& aaPos = domain.get_position_accessor();

		//	resize corners
		std::vector<MathVector<dim> > vCornerCoord;
		//	compute center of mass
		MathVector<dim> center;
		std::vector<MathVector<dim> > grad;
		center=0;
		size_t noc=vVrt.size();
		number u[noc];
    	for(size_t i = 0; i < noc; ++i){
			vCornerCoord.push_back( aaPos[vVrt[i]] );
			grad.push_back( aaGrad[vVrt[i]] );
			VecAppend(center,vCornerCoord[i]);
			dd.inner_multi_indices(vVrt[i], 0, ind);
			u[i]=BlockRef(uOld[ind[0][0]],ind[0][1]);
		};
		center/=noc;
		for (size_t i=0;i<noc;++i){
			number alpha=1;
			MathVector<dim> distVec;
			number uCenter;
    		VecSubtract(distVec,center,vCornerCoord[i]);
			uCenter = u[i] + distVec*grad[i];
			if (uCenter>u[i]){
		        if (uCenter>aaMax[vVrt[i]]) alpha=(aaMax[vVrt[i]]-u[i])/(distVec*grad[i]);
		        if (alpha<1){
		        	// UG_LOG("* " << vCornerCoord[i] << " uCenter " << uCenter << " ui " << u[i] << " max " << aaMax[vVrt[i]] << " alpha " << alpha << "\n");
			        aaGrad[vVrt[i]]*=alpha;
		        };
			}else{
		        if (uCenter<aaMin[vVrt[i]]) alpha=(aaMin[vVrt[i]]-u[i])/(distVec*grad[i]);
		        if (alpha<1){
		        	// UG_LOG("*#* " << vCornerCoord[i] << " uCenter " << uCenter << " ui " << u[i] << " min " << aaMin[vVrt[i]] << " alpha " << alpha << "\n");
			        aaGrad[vVrt[i]]*=alpha;
		        };
		    };
		};
	};
	};
	//	detach from grid
	grid.detach_from_vertices(aMin);
	grid.detach_from_vertices(aMax);
	return true;
};


/// interpolates a function on an element
template<typename TGridFunction>
template <typename TElem>
bool FV1LevelSetDisc<TGridFunction>::assemble_element(TElem& elem, DimFV1Geometry<dim>& geo, grid_type& grid,TGridFunction& uNew,const TGridFunction& uOld,aaGrad& aaGradient, aaSCV& aaVolume )
{
	// 	Type of multi index vector
	typedef typename dof_distribution_type::multi_index_vector_type multi_index_vector_type;

	//	get domain
	domain_type& domain = uNew.get_domain();

	//	get grid of domain
	// typename domain_type::grid_type& grid = domain.get_grid();

	//	create Multiindex
	multi_index_vector_type multInd;

	//	get element iterator type
	typedef typename domain_traits<dim>::const_iterator ElemIterator;

	//	get element type
	typedef typename domain_traits<dim>::geometric_base_object ElemType;

	//	hard code function (fct=0)
	//\todo: generalize
	size_t fct=0;

	//	id of shape functions used
	LFEID id = uOld.local_finite_element_id(fct);

	//	get vertices of the Elem
	std::vector<VertexBase*> vVrt;
	CollectVertices(vVrt, grid, elem);

	//	get position accessor
	typedef typename domain_type::position_accessor_type position_accessor_type;
	const position_accessor_type& aaPos = domain.get_position_accessor();

    //	resize corners
	std::vector<MathVector<dim> > coCoord;

	//	extract corner coordinates
	for(size_t i = 0; i < vVrt.size(); ++i)
	   coCoord.push_back( aaPos[vVrt[i]] );

	// update fv geometry
	geo.update(elem, &(coCoord[0]), &domain.get_subset_handler());

    //UG_LOG("geo.num_scvf() " << geo.num_scvf() << "\n");
    //UG_LOG("geo.num_scv() " << geo.num_scv() << "\n");
		
//	read indices on vertex
	typedef typename TGridFunction::dof_distribution_type dof_distribution_type;
	const dof_distribution_type& dd = uOld.get_dof_distribution();

//  fill node value vector
	std::vector<number> uValue(geo.num_scv());
	for (size_t i=0;i < geo.num_scv();i++){
		// if (dd.template inner_multi_indices<VertexBase>(vVrt[i], 0, multInd) != 1) return false;
		dd.inner_multi_indices(vVrt[i], 0, multInd);
		uValue[i]=BlockRef(uOld[multInd[0][0]],multInd[0][1]);
		//UG_LOG(coCoord[i] << "corner "<< i << " value " << uValue[i] << "\n");
	}

	//UG_LOG("uValues " << uValue[0] << " " << uValue[1] << " " << uValue[2] << "\n");

//  fill grad vector
    MathVector<dim> grad[geo.num_scv()];
	for (size_t i=0;i < geo.num_scv();i++){
		// if (dd.template inner_multi_indices<VertexBase>(vVrt[i], 0, multInd) != 1) return false;
		dd.inner_multi_indices(vVrt[i], 0, multInd);
		grad[i]=aaGradient[vVrt[i]];
		//UG_LOG("corner " << i << " gradient " << grad[i] << "\n");
	}
	
	//  fill corner velocity vector
	MathVector<dim> coVelocity[geo.num_scv()];
	for (size_t i=0;i<geo.num_scv();++i){
		coVelocity[i]=0;
	}
	if (m_gamma!=0){
        switch(m_velocity_type){
            case HardcodedData:
        	    for (size_t i=0;i < geo.num_scv();i++) analytic_velocity(coVelocity[i],m_time,coCoord[i]);
        	    break;
            case FunctorData:
        	    for (size_t i=0;i < geo.num_scv();i++){
        		    m_vel_x_fct(coVelocity[i][0],coCoord[i],m_time);
        	        if (dim>=2) m_vel_y_fct(coVelocity[i][1],coCoord[i],m_time);
        	        if (dim>=3) m_vel_z_fct(coVelocity[i][2],coCoord[i],m_time);
        	    }
        	    break;
            case VectorData:
        	    for (size_t i=0;i < geo.num_scv();i++){
        		    dd.inner_multi_indices(vVrt[i], 0, multInd);
        		    coVelocity[i][0]=BlockRef((*m_vel_x_vec)[multInd[0][0]],multInd[0][1]);
        		    if (dim>=2) coVelocity[i][1]=BlockRef((*m_vel_y_vec)[multInd[0][0]],multInd[0][1]);
        		    if (dim>=3) coVelocity[i][2]=BlockRef((*m_vel_z_vec)[multInd[0][0]],multInd[0][1]);
        	    }
        	    break;
            case ConstantData:
        	    for (size_t i=0;i < geo.num_scv();i++){
        	        coVelocity[i][0] = m_constantv_x;
        		    if (dim>=2) coVelocity[i][1] = m_constantv_y;
        		    if (dim>=3) coVelocity[i][2] = m_constantv_z;
        	    };
        }
        if (m_gamma!=1) for (size_t i=0;i < geo.num_scv();i++) coVelocity[i]*=m_gamma;
	}
    if (m_delta!=0){
    	for (size_t i=0;i < geo.num_scv();i++){
    	    number vnorm = VecLength(grad[i]);
    	    if (vnorm>1e-15) for (int j=0;j<dim;j++) coVelocity[i][j] += m_delta/vnorm*grad[i][j];
    	    //UG_LOG("corner " << i << " gradient" << grad[i] << " velocity " << coVelocity[i] << "\n");
    	};
    };
//  fill ipVel velocity vector
	MathVector<dim> ipVelocity[geo.num_scvf()];
	//UG_LOG("num scvf " << geo.num_scvf() << "\n");
	if ((m_interpolate_v_in_ip==true)||(m_delta!=0)){
		for (size_t ip=0;ip < geo.num_scvf();ip++){
			ipVelocity[ip] = 0;
			const typename DimFV1Geometry<dim>::SCVF& scvf = geo.scvf(ip);
			for (size_t co=0;co < geo.num_scv();co++){
			   // UG_LOG("co " << co << "\n");
				for (int j=0;j<dim;j++){
				    ipVelocity[ip][j] += scvf.shape(co)*coVelocity[co][j];
				    //UG_LOG("ip " << ip << " shape " << co << "[" << j << "]=" << scvf.shape(co) << " co velocity " << j << "=" << coVelocity[co][j]<<"\n");
				};
			};
		}
	} else
	{
		switch(m_velocity_type){
			case HardcodedData:
				for (size_t ip=0;ip < geo.num_scvf();ip++){
					const typename DimFV1Geometry<dim>::SCVF& scvf = geo.scvf(ip);
					MathVector<dim>	ipCoord = scvf.global_ip();
					analytic_velocity(ipVelocity[ip],m_time,ipCoord);
				}
				break;
			case FunctorData:
				for (size_t ip=0;ip < geo.num_scvf();ip++){
					const typename DimFV1Geometry<dim>::SCVF& scvf = geo.scvf(ip);
					MathVector<dim>	ipCoord = scvf.global_ip();
				    m_vel_x_fct(ipVelocity[ip][0],ipCoord,m_time);
			        if (dim>=2) m_vel_y_fct(ipVelocity[ip][1],ipCoord,m_time);
				    if (dim>=3) m_vel_z_fct(ipVelocity[ip][2],ipCoord,m_time);
				};
				break;
			case VectorData:
				UG_LOG("Exact velocity in ip can not be used when using vector valued velocity field. Set m_interpolate_v_in_ip=false.\n");
				break;
			case ConstantData:
				for (size_t ip=0;ip < geo.num_scvf();ip++){
				    ipVelocity[ip][0] = m_constantv_x;
				    if (dim>=2) ipVelocity[ip][1] = m_constantv_y;
				    if (dim>=3) ipVelocity[ip][2] = m_constantv_z;
				};
				break;
		}
	};
    //  fill source vector
	std::vector<number> coSource(geo.num_scv());
    switch (m_source_type){
        case HardcodedData:
          	for (size_t i=0;i<geo.num_scv();i++) coSource[i] = analytic_source(m_time,coCoord[i]);
       	    break;
        case FunctorData:
        	for (size_t i=0;i<geo.num_scv();i++) m_source_fct(coSource[i],coCoord[i],m_time);
        	break;
        case VectorData:
        	for (size_t i=0;i<geo.num_scv();i++){
        		dd.inner_multi_indices(vVrt[i], 0, multInd);
	    	    coSource[i] = BlockRef((*m_source_vec)[multInd[0][0]],multInd[0][1]);
        	}
        	break;
        case ConstantData:
        	for (size_t i=0;i<geo.num_scv();i++) coSource[i] = m_source_constant;
        	break;
    }
    // get finite volume geometry
	// FV1Geometry<TElem,dim> geo;

    //	typename TGridFunction::vector_type& v_vec = *dynamic_cast<typename TGridFunction::vector_type*>(&u);

	size_t base;
	number flux;
	MathVector<dim> distVec;
	// compute fluxes
	for(size_t ip = 0; ip < geo.num_scvf(); ++ip)
	{
	    MathVector<dim> bNode;
	    // 	get current SCVF
	    const typename DimFV1Geometry<dim>::SCVF& scvf = geo.scvf(ip);
	    MathVector<dim>	ipCoord = scvf.global_ip();
	    size_t from = scvf.from();
	    size_t to   = scvf.to();
	    if (scvf.normal()*ipVelocity[ip]>0){
		    base = from;
		} else {
		    base = to;
		};
	    //UG_LOG("from " << coCoord[from] << " to " << coCoord[to] << "\n");
		VecSubtract(distVec, ipCoord,coCoord[base]);
		//UG_LOG("ip vel " << ipVelocity[ip] << " normal " << scvf.normal() << " v*n " << ipVelocity[ip]*scvf.normal() << "\n");
		//UG_LOG("distVec*grad[base] " << distVec*grad[base] << " coSource[base] " << coSource[base] << " grad[base]*coVelocity[base] " << grad[base]*coVelocity[base] << "\n");
		//UG_LOG("u^{n+0.5}_{ip_i} = " << uValue[base] + distVec*grad[base] + 0.5*m_dt*(coSource[base] - grad[base]*coVelocity[base]) << "\n");
		// flux = v * n * u_{ip(i)}^{n+0.5}
		flux = m_dt*(ipVelocity[ip]*scvf.normal())*( uValue[base] + (distVec*grad[base]) + 0.5*m_dt*(coSource[base] - (grad[base]*coVelocity[base])) );
		//UG_LOG(ip << " flux=" << flux << "\n");
		dd.inner_multi_indices(vVrt[from], 0, multInd);
        BlockRef(uNew[multInd[0][0]],multInd[0][1])-=flux/aaVolume[ vVrt[from] ];
		if (m_divFree==false){
		    BlockRef(uNew[multInd[0][0]],multInd[0][1])+=m_dt*(ipVelocity[ip]*scvf.normal())*(uValue[from] + 0.5*m_dt*(coSource[from] - (grad[from]*coVelocity[from])))/aaVolume[ vVrt[from] ];
		};
		//UG_LOG("source flux from " << m_dt*(ipVelocity[ip]*scvf.normal())*(uValue[from] + 0.5*m_dt*(coSource[from] - (grad[from]*coVelocity[from])))/aaVolume[ vVrt[from] ] << "\n");
		//UG_LOG("v*n=" << (ipVelocity[ip]*scvf.normal()) << " uExtra=" <<  (uValue[from] + 0.5*m_dt*(coSource[from] - (grad[from]*coVelocity[from]))) << " vol=" << aaVolume[ vVrt[from] ] << "\n");
        dd.inner_multi_indices(vVrt[to], 0, multInd);
        BlockRef(uNew[multInd[0][0]],multInd[0][1])+=flux/aaVolume[ vVrt[to] ];
		if (m_divFree==false){
		    BlockRef(uNew[multInd[0][0]],multInd[0][1])-=m_dt*(ipVelocity[ip]*scvf.normal())*(uValue[to] + 0.5*m_dt*(coSource[to] - (grad[to]*coVelocity[to])))/aaVolume[ vVrt[to] ];
		};
		//UG_LOG("source flux to " << m_dt*(ipVelocity[ip]*scvf.normal())*(uValue[to] + 0.5*m_dt*(coSource[to] - (grad[to]*coVelocity[to])))/aaVolume[ vVrt[to] ] << "\n");
		//UG_LOG("v*n=" << (ipVelocity[ip]*scvf.normal()) << " uExtra=" <<  (uValue[to] + 0.5*m_dt*(coSource[to] - (grad[to]*coVelocity[to]))) << " vol=" << aaVolume[ vVrt[to] ] << "\n");
        number localCFL = std::max(m_dt*abs(ipVelocity[ip]*scvf.normal())/aaVolume[ vVrt[from] ],m_dt*abs(ipVelocity[ip]*scvf.normal())/aaVolume[ vVrt[to] ] );
        //UG_LOG("localCFL " << localCFL << "\n");
        if (localCFL>m_maxCFL){
            m_maxCFL = localCFL;
        };
        //UG_LOG("global CFL " << m_maxCFL << "\n");
         //		number dist = VecTwoNorm(distVec);
	};
	// boundary
	//	evaluate finite volume geometry

		//UG_LOG("nr of bdry subsets: " << geo.num_boundary_subsets() << " nr of bdry faces: " << geo.num_bf() << "\n");
		if (geo.num_bf()>0){
		  	for (int si=0;si<uNew.num_subsets();++si)
		    {
		   		//UG_LOG("si=" << si << "  num_bf(si)=" << geo.num_bf(si) << "\n");
		   		for(size_t i = 0; i < geo.num_bf(si); ++i)
		   		{
		   	    // 	get current BF
		   			const typename DimFV1Geometry<dim>::BF& bf = geo.bf(si, i);
			    	const size_t nodeID = bf.node_id();
				    MathVector<dim> bipVelocity;
				    number bipU=0;
				    number bipSource=0;
				    MathVector<dim> bipGrad;
				    bipVelocity=0;
				    bipGrad=0;
				    const MathVector<dim>* globalGradVec = bf.global_grad_vector();
				    bipVelocity=0;
				    for (size_t co=0;co<geo.num_scv();co++){
					    //UG_LOG("corner " << co << " num_sh " << bf.num_sh() << "\n");
				        for (int j=0;j<dim;j++){
				            bipVelocity[j] += bf.shape(co)*coVelocity[co][j];
		                    bipGrad[j] += bf.shape(co) * globalGradVec[co][j];
				    	};
				        bipU += bf.shape(co) *uValue[co];
				        bipSource += bf.shape(co) * coSource[co];
				    }
    				//flux = m_dt*(bipVelocity*bf.normal())*( bipU + 0.5*m_dt*(bipSource - (bipGrad*bipVelocity)) );
    				flux = m_dt*(bipVelocity*bf.normal())*uValue[nodeID];// first order approximation
    				dd.inner_multi_indices(vVrt[nodeID], 0, multInd);
    				BlockRef(uNew[multInd[0][0]],multInd[0][1])-=flux/aaVolume[ vVrt[nodeID] ];
    				//UG_LOG("coord=" << bf.global_ip( )<< "vel=" << bipVelocity << " n=" << bf.normal() << " flux=" << flux << " source flux=" << m_dt*(bipVelocity*bf.normal())*(uValue[nodeID] + 0.5*m_dt*(coSource[nodeID] - (grad[nodeID]*coVelocity[nodeID])))/aaVolume[ vVrt[nodeID] ] << "\n");
    				//UG_LOG("v*n=" << (bipVelocity*bf.normal()) << " uExtra=" << (uValue[nodeID] + 0.5*m_dt*(coSource[nodeID] - (grad[nodeID]*coVelocity[nodeID]))) << " volume=" << aaVolume[ vVrt[nodeID] ] << "\n");
    				if (m_divFree==false){
    				    BlockRef(uNew[multInd[0][0]],multInd[0][1])+=m_dt*(bipVelocity*bf.normal())*(uValue[nodeID] + 0.5*m_dt*(coSource[nodeID] - (grad[nodeID]*coVelocity[nodeID])))/aaVolume[ vVrt[nodeID] ];
    				    //BlockRef(uNew[multInd[0][0]],multInd[0][1])+=m_dt*(bipVelocity*bf.normal())*(uValue[nodeID] )/aaVolume[ vVrt[nodeID] ];// first order approximation
    				};
		    	};
		     };
		 };

	// give out corner values for debug
//	for (size_t i=0;i < geo.num_scv();i++){
			// if (dd.template inner_multi_indices<VertexBase>(vVrt[i], 0, multInd) != 1) return false;
//			dd.inner_multi_indices(vVrt[i], 0, multInd);
//			uValue[i]=BlockRef(uNew[multInd[0][0]],multInd[0][1]);
			//UG_LOG(coCoord[i] << "NEW corner "<< i << " value " << uValue[i] << "\n");
//	}

	//UG_LOG("NEW uValues " << uValue[0] << " " << uValue[1] << " " << uValue[2] << "\n");
	// end of "for debug"

//	we're done
	return true;
}

/*************************

COMPUTE VOLUME OF CONTROL VOLUMES

**************************/
template <typename TGridFunction>
bool FV1LevelSetDisc<TGridFunction>::
calculate_vertex_vol(TGridFunction& u,aaSCV& aaScvVolume)
{
	//	get domain
		domain_type& domain = u.get_domain();

	//	get grid of domain
		typename domain_type::grid_type& grid = domain.get_grid();

	//	create a FV Geometry for the dimension
		DimFV1Geometry<dim> geo;

	//	element iterator type
		typedef typename domain_traits<dim>::const_iterator ElemIterator;

	//	element type
		typedef typename domain_traits<dim>::geometric_base_object ElemType;

	//	get position accessor
		typedef typename domain_type::position_accessor_type position_accessor_type;
		const position_accessor_type& aaPos = domain.get_position_accessor();

	//	sum up all contributions of the sub control volumes to one vertex in an attachment
		for(int si = 0; si < u.num_subsets(); ++si)
		{
		//	get iterators
			ElemIterator iter = u.template begin<ElemType>(si);
			ElemIterator iterEnd = u.template end<ElemType>(si);

		//	loop elements of dimension
			for(  ;iter !=iterEnd; ++iter)
			{
			//	get Elem
				ElemType* elem = *iter;
			//	get vertices of the element
				std::vector<VertexBase*> vVrt;
				CollectVertices(vVrt, grid, elem);
			//	resize corners
				std::vector<MathVector<dim> > vCornerCoord;

			//	extract corner coordinates
				for(size_t i = 0; i < vVrt.size(); ++i)
			    	vCornerCoord.push_back( aaPos[vVrt[i]] );

			//	evaluate finite volume geometry
				geo.update(elem, &(vCornerCoord[0]), &domain.get_subset_handler());


			//	loop corners
				for (size_t i=0;i < geo.num_scv();i++)
				{
					//	get scv for sh
					const typename DimFV1Geometry<dim>::SCV& scv = geo.scv(i);

					aaScvVolume[vVrt[i]] += scv.volume();
				};
			}
		}

	return true ;

}

/*************************

COMPUTE VERTEX GRAD AND VOLUME OF CONTROL VOLUME 

**************************/

template<typename TGridFunction>
bool FV1LevelSetDisc<TGridFunction>::
calculate_vertex_grad_vol(TGridFunction& u, aaGrad& aaGradient,aaSCV& aaScvVolume)
{
	//	domain type
	//	typedef typename TGridFunction::domain_type domain_type;

	/// dof distribution type
	//	typedef typename TGridFunction::dof_distribution_type dof_distribution_type;

	// 	dimension of world (and reference element)
	//	static const int dim = domain_type::dim;

	// 	Type of multi index vector
		typedef typename dof_distribution_type::multi_index_vector_type multi_index_vector_type;

	//	get domain
		domain_type& domain = u.get_domain();

	//	get grid of domain
		typename domain_type::grid_type& grid = domain.get_grid();

	//	create Multiindex
		multi_index_vector_type multInd;

	//	create a FV Geometry for the dimension
		DimFV1Geometry<dim> geo;

	//	get element iterator type
		typedef typename domain_traits<dim>::const_iterator ElemIterator;

	//	get element type
		typedef typename domain_traits<dim>::geometric_base_object ElemType;

	//	hard code function (fct=0)
	//\todo: generalize
		size_t fct=0;

		// initialize attachment value
		SetAttachmentValues(aaScvVolume, grid.vertices_begin(), grid.vertices_end(), 0);
		SetAttachmentValues(aaGradient, grid.vertices_begin(), grid.vertices_end(), 0);

	//	sum up all contributions of the sub control volumes to one vertex in an attachment
		for(int si = 0; si < u.num_subsets(); ++si)
		{
		//	get iterators
			ElemIterator iter = u.template begin<ElemType>(si);
			ElemIterator iterEnd = u.template end<ElemType>(si);

		//	loop elements of dimension
			for(  ;iter !=iterEnd; ++iter)
			{
			//	get Elem
				ElemType* elem = *iter;

			//	get vertices of the Elem
				std::vector<VertexBase*> vVrt;
				CollectVertices(vVrt, grid, elem);

			//	get position accessor
				typedef typename domain_type::position_accessor_type position_accessor_type;
				const position_accessor_type& aaPos = domain.get_position_accessor();

			//	resize corners
				std::vector<MathVector<dim> > vCornerCoord;

			//	extract corner coordinates
				for(size_t i = 0; i < vVrt.size(); ++i)
					vCornerCoord.push_back( aaPos[vVrt[i]] );

			//	evaluate finite volume geometry
				geo.update(elem, &(vCornerCoord[0]), &domain.get_subset_handler());

				//UG_LOG("Num Verts loaded: "<<vVrt.size()<<"\n");
				//UG_LOG("Num SCV computed: "<<geo.num_scv()<<"\n");m_vPP.
				UG_ASSERT(vVrt.size() == geo.num_scv(), "Must match.");
				for(size_t i = 0; i < vVrt.size(); ++i)
					if(vVrt[i] == NULL)
						UG_LOG("Vertex "<<i<<" is NULL.\n");

				std::vector<number> uValue(geo.num_scv());

				typedef typename TGridFunction::multi_index_vector_type index_type;

			//	read indices on vertex
				typedef typename TGridFunction::dof_distribution_type dof_distribution_type;
				const dof_distribution_type& dd = u.get_dof_distribution();

				for (size_t i=0;i < geo.num_scv();i++)
				{
				//	get indices of function fct on vertex
					dd.inner_multi_indices(vVrt[i], fct, multInd);

				//	read value of index from vector
					uValue[i]=BlockRef(u[multInd[0][0]],multInd[0][1]);

				//	debug log
					//UG_LOG("corner " << i << " " << uValue[i] << "\n");
				}

			//	storage for global gradient
				MathVector<dim> globalGrad;

			//	loop corners
				for (size_t i=0;i < geo.num_scv();i++)
				{
				//	get scv for sh
					const typename DimFV1Geometry<dim>::SCV& scv = geo.scv(i);

				//	debug log
					//UG_LOG("gradient for corner " << i << "\n");

				//	reset global gradient
					globalGrad = 0.0;

				//	sum up gradients of shape functions in corner
					for(size_t sh = 0 ; sh < geo.num_scv(); ++sh)
					{
						//UG_LOG("local grad " << sh << " : " << scv.local_grad(sh) << "\n");
						//UG_LOG("unscaled global grad " << sh << " = " << scv.global_grad(sh) << "\n");
						//UG_LOG("uvalue(" << sh << ") =" << uValue[sh] << "\n");

						VecScaleAppend(globalGrad, uValue[sh], scv.global_grad(sh));
					}

				//	volume of scv
					number vol = scv.volume();

					//UG_LOG("*** global grad " << i << ": " << globalGrad << "\n");

				//	scale gradient by volume
					globalGrad *= vol;

				//	add both values to attachements
					aaGradient[vVrt[i]] += globalGrad;
					aaScvVolume[vVrt[i]] += vol;
				};
			}
		}

	int count=0;
	typedef typename domain_type::position_accessor_type position_accessor_type;
	position_accessor_type aaPos = u.get_domain().get_position_accessor();
	for (int si=0;si < u.num_subsets();++si){
		//UG_LOG("si " << si << "\n");
	    for(VertexBaseConstIterator iter = u.template begin<VertexBase>(si);
						   iter != grid.template end<VertexBase>(si); ++iter)
	    {
	    //	get vertex
		    VertexBase* vrt = *iter;
		    if (aaScvVolume[vrt]!=0){
		        (aaGradient[vrt]) /= aaScvVolume[vrt];
		    };
		    MathVector<dim> coord = aaPos[vrt];
		    MathVector<dim> exact;
		    exact[0] = 0;
		    exact[1] = 0;
		    //exact[0] = cos(coord[0]);// 6*coord[0];
		    //exact[1] = -4*sin(coord[1]);//-4*coord[1];
		    //number gError = sqrt( (exact[0]-aaGradient[vrt][0])*(exact[0]-aaGradient[vrt][0]) + (exact[1]-aaGradient[vrt][1])*(exact[1]-aaGradient[vrt][1]) );
	        //UG_LOG(count << "[ " << coord[0] << "," << coord[1] << " ] vol= " << aaScvVolume[vrt] << " " << "grad= ["
	        //		<< aaGradient[vrt][0] << "," << aaGradient[vrt][1] << "] exact grad = ["
	        //		<< exact[0] << "," << exact[1] << "] error: " << gError <<  "\n");
	        count++;
	    }
	}
	//UG_LOG("#*#*#*#\n");
	return true ;
}

template<typename TGridFunction>
bool FV1LevelSetDisc<TGridFunction>::assign_dirichlet(TGridFunction& numsol){
	//	get domain of grid function
		domain_type& domain = numsol.get_domain();

	//	get grid type of domain
		typedef typename domain_type::grid_type grid_type;

	//	get grid of domain
		grid_type& grid = domain.get_grid();

		typedef typename domain_type::position_accessor_type position_accessor_type;

		UG_LOG("dirichlet\n");

		UG_LOG("nr dir ss" << m_dirichlet_sg.num_subsets() << "\n");

		for(size_t i = 0; i < m_dirichlet_sg.num_subsets(); ++i)
		{
	        const int si = m_dirichlet_sg[i];
	        UG_LOG("Dirichlet boundary is: "<<si<< "\n");
			for(VertexBaseConstIterator iter = numsol.template begin<VertexBase>(si);
										   iter != grid.template end<VertexBase>(si); ++iter)
			{
			//	get vertex
				VertexBase* vrt = *iter;
				number exactVal;
				MathVector<dim> coord;
				position_accessor_type aaPos = domain.get_position_accessor();
				coord = aaPos[vrt];

				typedef typename TGridFunction::multi_index_vector_type index_type;

			//	get vector holding all indices on the vertex
				multi_index_vector_type ind;

			//	read indices on vertex
				typedef typename TGridFunction::dof_distribution_type dof_distribution_type;
				const dof_distribution_type& dd = numsol.get_dof_distribution();

				const size_t numInd = dd.inner_multi_indices(vrt, 0, ind);

			//	check indices
				if(numInd != 1) {UG_LOG("ERROR: Wrong number of indices!"); return false;}

				switch (m_dirichlet_data_type){
				    case (HardcodedData):
		                exactVal= analytic_solution(m_time,coord);
				        break;
				    case (FunctorData):
				    	m_solution_fct(exactVal,coord,m_time);
				        break;
				    case (ConstantData):
				    	exactVal = m_dirichlet_constant;
				        break;
				    case (VectorData):
				    	break;
				}
				BlockRef(numsol[ind[0][0]],ind[0][1]) = exactVal;
				//UG_LOG("coord " << coord << "\n");

				//if ((coord[0]==-1)||(coord[0]==1)||(coord[1]==-1)||(coord[1]==1)){
			    	//BlockRef(numsol[ind[0][0]],ind[0][1]) = exactVal;
				//};
				//if ((coord[0]==0)||(coord[0]==1)||(coord[1]==0)||(coord[1]==1)){
				//};
		     }
		};
	return true;
}

template<typename TGridFunction>
bool FV1LevelSetDisc<TGridFunction>::compute_error(TGridFunction& numsol)
{
//	get domain of grid function
	domain_type& domain = numsol.get_domain();

//	get grid type of domain
	typedef typename domain_type::grid_type grid_type;

//	get grid of domain
	grid_type& grid = domain.get_grid();

	typedef typename domain_type::position_accessor_type position_accessor_type;

	//	create Attachment for scv-volume size
	ANumber aScvVolume;

	//	attach to grid
	grid.attach_to_vertices(aScvVolume);

	//	get attachment accessor to access values
	Grid::VertexAttachmentAccessor<ANumber> aaScvVolume(grid, aScvVolume);

	// initialize attachment value
	SetAttachmentValues(aaScvVolume, grid.vertices_begin(), grid.vertices_end(), 0);

    number l1Error=0;
	number l2Error=0;
	number maxErr=0;

	bool bRes = true;
	// calculate scv size
	if (calculate_vertex_vol(numsol,aaScvVolume)==false){UG_LOG("ERROR: gradient computation failed in compute_error function!"); };

	//UG_LOG("----------------------------\n");

	if(!bRes) {UG_LOG("Error while calculating CV Volume.\n"); return false;}
	for (int si=0;si<numsol.num_subsets();++si){
		for(VertexBaseConstIterator iter = numsol.template begin<VertexBase>(si);
									   iter != grid.template end<VertexBase>(si); ++iter)
		{
		//	get vertex
			VertexBase* vrt = *iter;
			number exactVal;
			MathVector<dim> coord;
			position_accessor_type aaPos = domain.get_position_accessor();
			coord = aaPos[vrt];
		
			typedef typename TGridFunction::multi_index_vector_type index_type;

		//	get vector holding all indices on the vertex
			multi_index_vector_type ind;

		//	read indices on vertex
			typedef typename TGridFunction::dof_distribution_type dof_distribution_type;
			const dof_distribution_type& dd = numsol.get_dof_distribution();

			const size_t numInd = dd.inner_multi_indices(vrt, 0, ind);

		//	check indices
			if(numInd != 1) {UG_LOG("ERROR: Wrong number of indices!"); return false;}

			switch (m_dirichlet_data_type){
			    case (HardcodedData):
				     exactVal= analytic_solution(m_time,coord);
				     break;
				case (FunctorData):
				   	m_solution_fct(exactVal,coord,m_time);
				    break;
				case (ConstantData):
				   	exactVal = m_dirichlet_constant;
				    break;
				case (VectorData):
					break;
			}
		
			number differ = abs(BlockRef(numsol[ind[0][0]],ind[0][1])-exactVal);
		
			l1Error += aaScvVolume[vrt] * differ;
			l2Error += aaScvVolume[vrt] * differ*differ;

			if (m_print==true){
			    if (differ>0)
			    UG_LOG("coord=" << coord << " value=" << BlockRef(numsol[ind[0][0]],ind[0][1]) << " exact=" << exactVal << " error=" << differ << "\n");
			};
		
			if (differ > maxErr) maxErr = differ;

	     }
	};
	l2Error = sqrt(l2Error);
	UG_LOG("timestep " << m_timestep_nr << " time " << m_time << "\n");
	UG_LOG("l1 error: " << l1Error << "\n");
	UG_LOG("l2 error: " << l2Error << "\n");
	UG_LOG("maximum error: " << maxErr << "\n");
	return true;	
};

template<typename TGridFunction>
bool FV1LevelSetDisc<TGridFunction>::init_function(TGridFunction& u)
{
//	get domain of grid function
	domain_type& domain = u.get_domain();

//	get grid type of domain
	typedef typename domain_type::grid_type grid_type;

//	get grid of domain
	grid_type& grid = domain.get_grid();

	typedef typename domain_type::position_accessor_type position_accessor_type;

	//	read indices on vertex
	position_accessor_type aaPos = domain.get_position_accessor();
	typedef typename TGridFunction::multi_index_vector_type index_type;
	typedef typename TGridFunction::dof_distribution_type dof_distribution_type;
	const dof_distribution_type& dd = u.get_dof_distribution();


	for (int si=0;si<u.num_subsets();++si){
		for(VertexBaseConstIterator iter = u.template begin<VertexBase>(si);
									   iter != grid.template end<VertexBase>(si); ++iter)
		{
		//	get vertex
			VertexBase* vrt = *iter;
			MathVector<dim> coord;
			coord = aaPos[vrt];
		//	get vector holding all indices on the vertex
			multi_index_vector_type ind;
			dd.inner_multi_indices(vrt, 0, ind);
			BlockRef(u[ind[0][0]],ind[0][1]) = analytic_solution(m_time,coord);
	     }
	};
	return true;
};

template<typename TGridFunction>
bool FV1LevelSetDisc<TGridFunction>::runtimetest(TGridFunction& uNew)
{
	//	get domain of grid function
	domain_type& domain = uNew.get_domain();

	//	get grid type of domain
	typedef typename domain_type::grid_type grid_type;

	//	get grid of domain
	grid_type& grid = domain.get_grid();

	//	element iterator type
	typedef typename domain_traits<dim>::const_iterator ElemIterator;

    //	element type
	typedef typename domain_traits<dim>::geometric_base_object ElemType;
	
	//	create a FV Geometry for the dimension
	DimFV1Geometry<dim> geo;

	/*	bool testfvgeom=false;

	if (testfvgeom==true){

	for (int si=0;si<uNew.num_subsets();++si){
	    //	get iterators
	    ElemIterator iter = uNew.template begin<ElemType>(si);
	    ElemIterator iterEnd = uNew.template end<ElemType>(si);
		//	loop elements of dimension
		for(  ;iter !=iterEnd; ++iter)
		{
		//	get Elem
			ElemType* elem = *iter;
		//	get vertices of the Elem
			std::vector<VertexBase*> vVrt;
			CollectVertices(vVrt, grid, elem);
		//	get position accessor
			typedef typename domain_type::position_accessor_type position_accessor_type;
			const position_accessor_type& aaPos = domain.get_position_accessor();
		//	resize corners
			std::vector<MathVector<dim> > vCornerCoord;
		//	extract corner coordinates
			for(size_t i = 0; i < vVrt.size(); ++i)
				vCornerCoord.push_back( aaPos[vVrt[i]] );

		//	evaluate finite volume geometry
			geo.update(elem, &(vCornerCoord[0]), &domain.get_subset_handler());
		};
	};

	}

	for (int si=0;si<uNew.num_subsets();++si){
		//	get iterators
		ElemIterator iter = uNew.template begin<ElemType>(si);
		ElemIterator iterEnd = uNew.template end<ElemType>(si);
		//	loop elements of dimension
		for(  ;iter !=iterEnd; ++iter)
		{
			//	get Elem
			ElemType* elem = *iter;
			//	get vertices of the Elem
			std::vector<VertexBase*> vVrt;
			CollectVertices(vVrt, grid, elem);
			//	get position accessor
			typedef typename domain_type::position_accessor_type position_accessor_type;
			const position_accessor_type& aaPos = domain.get_position_accessor();
			//	resize corners
			std::vector<MathVector<dim> > vCornerCoord;
			//	extract corner coordinates
			for(size_t i = 0; i < vVrt.size(); ++i)
				vCornerCoord.push_back( aaPos[vVrt[i]] );
		};
	};
	
	*/
	for (int si=0;si<uNew.num_subsets();++si){
		//	get iterators
		ElemIterator iter = uNew.template begin<ElemType>(si);
		ElemIterator iterEnd = uNew.template end<ElemType>(si);
		//	loop elements of dimension
		for(  ;iter !=iterEnd; ++iter)
		{
			//	get Elem
			ElemType* elem = *iter;
			//	get vertices of the Elem
			std::vector<VertexBase*> vVrt;
			CollectVertices(vVrt, grid, elem);
		};
	};
/*
	MathVector<dim> coord;
	typedef typename domain_type::position_accessor_type position_accessor_type;
	position_accessor_type aaPos = domain.get_position_accessor();

	for (int si=0;si<uNew.num_subsets();++si){
		VertexBaseConstIterator iter = uNew.template begin<VertexBase>(si);
	    VertexBaseConstIterator iterEnd = uNew.template end<VertexBase>(si);
		for(;iter != iterEnd; ++iter)
		{
		//	get vertex
			VertexBase* vrt = *iter;
			coord = aaPos[vrt];
		}
	};*/
	return true;
}

template<typename TGridFunction>
bool FV1LevelSetDisc<TGridFunction>::advect_lsf(TGridFunction& uNew,TGridFunction& uOld)
{
	//	get domain of grid function
	domain_type& domain = uNew.get_domain();

	//	get grid type of domain
	typedef typename domain_type::grid_type grid_type;

	//	get grid of domain
	grid_type& grid = domain.get_grid();

	typedef typename domain_type::position_accessor_type position_accessor_type;

	//	element iterator type
	typedef typename domain_traits<dim>::const_iterator ElemIterator;

    //	element type
	typedef typename domain_traits<dim>::geometric_base_object ElemType;


	//	get element iterator type
	typedef typename domain_traits<dim>::const_iterator ElemIterator;

	//	create Attachment for scv-volume size
	ANumber aScvVolume;

	//	typedef of gradient attachment
	typedef Attachment<MathVector<dim> > AGradient;

	//	create Attachment for gradient
	AGradient aGradient;

	//	attach to grid
	grid.attach_to_vertices(aScvVolume);
	grid.attach_to_vertices(aGradient);

	//	get attachment accessor to access values
	Grid::VertexAttachmentAccessor<ANumber> aaScvVolume(grid, aScvVolume);
	Grid::VertexAttachmentAccessor<AGradient> aaGradient(grid, aGradient);

	// initialize attachment value
	SetAttachmentValues(aaScvVolume, grid.vertices_begin(), grid.vertices_end(), 0);
	SetAttachmentValues(aaGradient, grid.vertices_begin(), grid.vertices_end(), 0);


	// const dof_distribution_type& dd = uNew.get_dof_distribution();
	m_maxCFL=0;
	
	//UG_LOG("***************************************************\n");
    //UG_LOG("***************************************************\n");
	//UG_LOG("***************************************************\n");

	/*for (int si=0;si<uNew.num_subsets();++si){
		for(VertexBaseConstIterator iter = uNew.template begin<VertexBase>(si);
								   iter != grid.template end<VertexBase>(si); ++iter)
	    {
		    VertexBase* vrt = *iter;
		    typedef typename TGridFunction::multi_index_vector_type index_type;

	   //	get vector holding all indices on the vertex
		    index_type ind;

    	//	read indices on vertex
		    if (dd.inner_multi_indices(vrt, 0, ind)!=1){UG_LOG("ERROR: Wrong number of indices!"); return false;};

		    BlockRef(uNew[ind[0][0]],ind[0][1]) = BlockRef(uOld[ind[0][0]],ind[0][1]);
		    //UG_LOG("uNew: " << BlockRef(uNew[ind[0][0]],ind[0][1]) << "uOld: " << BlockRef(uOld[ind[0][0]],ind[0][1]) << "\n");
	    }
	};*/
	VecAssign(uNew,uOld);
    MathVector<dim> coord;
    position_accessor_type aaPos = domain.get_position_accessor();

	typedef typename TGridFunction::multi_index_vector_type index_type;
    //	get vector holding all indices on the vertex
	multi_index_vector_type ind;
    //	read indices on vertex
	typedef typename TGridFunction::dof_distribution_type dof_distribution_type;
	const dof_distribution_type& dd = uNew.get_dof_distribution();

	for (size_t step=0;step<m_nrOfSteps;step++)
	{
	    // calculate scv size and gradient
	    if (calculate_vertex_grad_vol(uNew,aaGradient, aaScvVolume)==false){UG_LOG("ERROR: gradient computation failed!"); };
	    // SetAttachmentValues(aaGradient, grid.vertices_begin(), grid.vertices_end(), 0); for debug set gradient to 0
	    if (m_limiter==true){
	    	limit_grad(uNew,aaGradient);
	    };

	    for (int si=0;si<uNew.num_subsets();++si){
	        if (m_dirichlet_sg.num_subsets()!=0) if (m_dirichlet_sg.contains(si)==true) continue;
	        if (m_neumann_sg.num_subsets()!=0) if (m_neumann_sg.contains(si)==true) continue;
	        if (m_inactive_sg.num_subsets()!=0) if (m_inactive_sg.contains(si)==true) continue;
		    //UG_LOG("***************************************************\n");
	        //UG_LOG("***********************" << si << "**************************\n");
		    //UG_LOG("***************************************************\n");
	    	// test edge iterator
	//    	for(EdgeBaseConstIterator iter = uNew.template begin<EdgeBase>(si) ;iter !=uNew.template end<EdgeBase>(si); ++iter)
	  //  	{
    //            EdgeBase* edge = *iter;
   //             VertexBase* v0=edge->vertex(0);
     //           VertexBase* v1=edge->vertex(1);
            // UG_LOG(" coord vertex 0 " << aaPos[v0] << " coord vertex 1 " << aaPos[v1] << "\n");
	    //	}

		    //	get iterators
		    ElemIterator iter = uNew.template begin<ElemType>(si);
		    ElemIterator iterEnd = uNew.template end<ElemType>(si);

		    DimFV1Geometry<dim> geo;

		    // flag given neumann bnd subsets at the geometry, such that the
		    //	geometry produces boundaryfaces (BF) for all sides of the
		    //	element, that is in one of the subsets
		    for(size_t i = 0; i < m_neumann_sg.num_subsets(); ++i)
		    {
		        const int bndSi = m_neumann_sg[i];
		        UG_LOG("Neumann boundary is: "<<bndSi<< "\n");
		        geo.add_boundary_subset(bndSi);
		    }

		    //	loop elements of dimension
		    for(  ;iter !=iterEnd; ++iter)
		    {
		        //	get Elem
			    //UG_LOG("*** ELEM ***\n");
			    ElemType* elem = *iter;
			    //UG_LOG("element \n");
			    // uNew = uOld
			    assemble_element(elem, geo, uNew. get_domain().get_grid(),uNew,uOld,aaGradient,aaScvVolume);
		    };
			if (m_source==true){
			    for(VertexBaseConstIterator iter = uNew.template begin<VertexBase>(si);
									   iter != grid.template end<VertexBase>(si); ++iter)
     		    {
        		    //	get vertex
			    	//UG_LOG("*** VERTEX ***\n");
    	    		VertexBase* vrt = *iter;
    	    		coord = aaPos[vrt];
	    			dd.inner_multi_indices(vrt, 0, ind);
	    			//UG_LOG("uNew " << BlockRef(uNew[ind[0][0]],ind[0][1]) << " ");
		    		BlockRef(uNew[ind[0][0]],ind[0][1]) += m_dt*analytic_source(m_time+0.5*m_dt,coord);
		    		//UG_LOG("uNew " << BlockRef(uNew[ind[0][0]],ind[0][1]) << "\n");
		    		//UG_LOG("coord=" << aaPos[vrt] << " value=" << BlockRef(uNew[ind[0][0]],ind[0][1]) << " error=" << abs(BlockRef(uNew[ind[0][0]],ind[0][1])-analytic_solution(m_time+m_dt,aaPos[vrt])) << "\n");
	            }
		   };
	    }
	    m_time += m_dt;
	    UG_LOG("m_dt " << m_dt << " m_time " << m_time << "\n");
	    m_timestep_nr++;
	    assign_dirichlet(uNew);

        //	now process dirichlet values
	    for(size_t i = 0; i < m_vPP.size(); ++i)
	    {
	    //	UG_LOG(i << "\n");
	    //	m_vPP[i]->adjust_solution(uNew, uNew.get_dof_distribution(),m_time);
	    }
	    //if(m_print){
	    	UG_LOG("step: " << m_timestep_nr << "\n");
	        UG_LOG("time: " << m_time << "\n");
            UG_LOG("max CFL nr " << m_maxCFL << "\n");
	    //};
	    if (m_nrOfSteps>1){
	    	// Attention, uOld is overwritten
	    	VecAssign(uOld,uNew);
	    };
	};
    //	detach from grid
	grid.detach_from_vertices(aScvVolume);
	grid.detach_from_vertices(aGradient);

//	done
	return true;
}

} // end namespace ug

#endif /* LEVEL_SET_UTIL_IMPL_H_ */
