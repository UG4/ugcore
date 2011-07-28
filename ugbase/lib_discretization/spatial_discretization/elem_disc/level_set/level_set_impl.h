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
		//   v[0] = 1;
		//   v[1] = 1;
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
		   return 0;
		case 2: 
		   return 0;
		case 3: 
		   return 0;
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
		    xnew[0] = x[0]*cos(t)+x[1]*sin(t);
			xnew[1] = -x[0]*sin(t)+x[1]*cos(t);
        //    xnew[0] = x[0] - t;
        //    xnew[1] = x[1] - t;
		//	return x[0]-t;
		//	return sqrt(xnew[0]*xnew[0] + xnew[1]*xnew[1]) - 0.5;
            return 2*xnew[0]-xnew[1];
        case 3:  
			return 0;
	}

    return -1;
};

/// interpolates a function on an element
template<typename TGridFunction>
template <typename TElem>
bool FV1LevelSetDisc<TGridFunction>::assemble_element(TElem& elem,grid_type& grid,TGridFunction& uNew,const TGridFunction& uOld,aaGrad& aaGradient, aaSCV& aaVolume )
{
	// 	Type of multi index vector
	typedef typename dof_distribution_type::multi_index_vector_type multi_index_vector_type;

	//	get domain
	domain_type& domain = uNew.get_domain();

	//	get grid of domain
	// typename domain_type::grid_type& grid = domain.get_grid();

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

	//	evaluate finite volume geometry
	geo.update(elem, domain.get_subset_handler(), &(coCoord[0]) );
    //UG_LOG("geo.num_scvf() " << geo.num_scvf() << "\n");
    //UG_LOG("geo.num_scv() " << geo.num_scv() << "\n");
		
    //  fill corner velocity vector
    MathVector<dim> coVelocity[geo.num_scv()];
	for (size_t i=0;i < geo.num_scv();i++){
	     if (m_analyticalVelocity){
		     analytic_velocity(coVelocity[i],m_time,coCoord[i]);
		 };
	};
//  fill ipVel velocity vector
	MathVector<dim> ipVelocity[geo.num_scvf()];
	for (size_t i=0;i < geo.num_scvf();i++){
	     if (m_analyticalVelocity){
		     analytic_velocity(ipVelocity[i],m_time,coCoord[i]);
		 };
	};

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
	
//  fill source vector
	std::vector<number> coSource(geo.num_scv());
    if (m_source){
	    if (m_analyticalSource){
	    	for (size_t i=0;i<geo.num_scv();i++){
	    	    coSource[i] = analytic_source(m_time,coCoord[i]);
	    	}
		};
	}else{
	    for (size_t i=0;i<geo.num_scv();i++){
		    coSource[i] = 0;
		};
	};
	
    // get finite volume geometry
	// FV1Geometry<TElem,dim> geo;

    //	typename TGridFunction::vector_type& v_vec = *dynamic_cast<typename TGridFunction::vector_type*>(&u);

	size_t base;
	number flux;
	MathVector<dim> distVec;
		
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
	    //UG_LOG("from " << from << " to " << to << " base " << base << "\n");
		VecSubtract(distVec, ipCoord,coCoord[base]);
		//UG_LOG("ip vel " << ipVelocity[ip] << " normal " << scvf.normal() << " v*n " << ipVelocity[ip]*scvf.normal() << "\n");
		//UG_LOG("distVec*grad[base] " << distVec*grad[base] << " coSource[base] " << coSource[base] << " grad[base]*coVelocity[base] " << grad[base]*coVelocity[base] << "\n");
		//UG_LOG("u^{n+0.5}_{ip_i} = " << uValue[base] + distVec*grad[base] + 0.5*m_dt*(coSource[base] - grad[base]*coVelocity[base]) << "\n");
		// flux = v * n * u_{ip(i)}^{n+0.5}
		flux = m_dt*(ipVelocity[ip]*scvf.normal())*( uValue[base] + distVec*grad[base] + 0.5*m_dt*(coSource[base] - grad[base]*coVelocity[base]) );
		//UG_LOG(ip << " flux=" << flux << "\n");
		dd.inner_multi_indices(vVrt[from], 0, multInd);
        BlockRef(uNew[multInd[0][0]],multInd[0][1])-=flux/aaVolume[ vVrt[from] ];
        dd.inner_multi_indices(vVrt[to], 0, multInd);
        BlockRef(uNew[multInd[0][0]],multInd[0][1])+=flux/aaVolume[ vVrt[to] ];
        number localCFL = std::max(m_dt*abs(ipVelocity[ip]*scvf.normal())/aaVolume[ vVrt[from] ],m_dt*abs(ipVelocity[ip]*scvf.normal())/aaVolume[ vVrt[to] ] );
        //UG_LOG("localCFL " << localCFL << "\n");
        if (localCFL>m_maxCFL){
            m_maxCFL = localCFL;
        };
        //UG_LOG("global CFL " << m_maxCFL << "\n");
         //		number dist = VecTwoNorm(distVec);
	};


	// for debug: new corner values:
	for (size_t i=0;i < geo.num_scv();i++){
			// if (dd.template inner_multi_indices<VertexBase>(vVrt[i], 0, multInd) != 1) return false;
			dd.inner_multi_indices(vVrt[i], 0, multInd);
			uValue[i]=BlockRef(uNew[multInd[0][0]],multInd[0][1]);
			//UG_LOG(coCoord[i] << "NEW corner "<< i << " value " << uValue[i] << "\n");
	}

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
				geo.update(elem, domain.get_subset_handler(), &(vCornerCoord[0]) );


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
				geo.update(elem, domain.get_subset_handler(), &(vCornerCoord[0]) );

				//UG_LOG("Num Verts loaded: "<<vVrt.size()<<"\n");
				//UG_LOG("Num SCV computed: "<<geo.num_scv()<<"\n");
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

				exactVal= analytic_solution(m_time,coord);

				if ((coord[0]==-1)||(coord[0]==1)||(coord[1]==-1)||(coord[1]==1)){
			    	BlockRef(numsol[ind[0][0]],ind[0][1]) = exactVal;
				};
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
		
			exactVal= analytic_solution(m_time,coord);
		
			number differ = abs(BlockRef(numsol[ind[0][0]],ind[0][1])-exactVal);
		
			l1Error += aaScvVolume[vrt] * differ;
			l2Error += aaScvVolume[vrt] * differ*differ;
		
			if (differ > maxErr) maxErr = differ;

	     }
	};
	l2Error = sqrt(l2Error);
	UG_LOG("l1 error: " << l1Error << "\n");
	UG_LOG("l2 error: " << l2Error << "\n");
	UG_LOG("maximum error: " << maxErr << "\n");
	return true;	
};


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
	VecScaleAssign(uNew,1,uOld);

	for (size_t step=0;step<m_nrOfSteps;step++)
	{
	    // calculate scv size and gradient
	    if (calculate_vertex_grad_vol(uNew,aaGradient, aaScvVolume)==false){UG_LOG("ERROR: gradient computation failed!"); };
        // SetAttachmentValues(aaGradient, grid.vertices_begin(), grid.vertices_end(), 0);// for debug set gradient to 0

	    for (int si=0;si<uNew.num_subsets();++si){
		    //UG_LOG("***************************************************\n");
	        //UG_LOG("***********************" << si << "**************************\n");
		    //UG_LOG("***************************************************\n");

		    //	get iterators
		    ElemIterator iter = uNew.template begin<ElemType>(si);
		    ElemIterator iterEnd = uNew.template end<ElemType>(si);
		    //	loop elements of dimension
		    for(  ;iter !=iterEnd; ++iter)
		    {
		        //	get Elem
			    //UG_LOG("*** ELEM ***\n");
			    ElemType* elem = *iter;
			    //UG_LOG("element \n");
			    // uNew = uOld
			    assemble_element(elem, uNew.get_domain().get_grid(),uNew,uOld,aaGradient,aaScvVolume);
		    };
	    }
	    m_time += m_dt;
	    assign_dirichlet(uNew);

        //	now process dirichlet values
	    for(size_t i = 0; i < m_vPP.size(); ++i)
	    {
	    //	m_vPP[i]->adjust_solution(uNew, uNew.get_dof_distribution,m_time);
	    }
	    if(m_print){
	        UG_LOG("time: " << m_time << "\n");
            UG_LOG("max CFL nr " << m_maxCFL << "\n");
	    };
	    if (m_nrOfSteps>1){
	    	VecScaleAssign(uOld,1,uNew);
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
