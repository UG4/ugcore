/*
 * ElementalUserData.h
 *
 *  Created on: 10.03.2025
 *      Author: jwang
 */

#ifndef _H__UG__LIB_DISC__SPATIAL_DISC__HOMORASTER_USERDATA_H_
#define _H__UG__LIB_DISC__SPATIAL_DISC__HOMORASTER_USERDATA_H_
 
 
#include <vector>
#include <string>
 //#include <utility> // std::pair and std::make_pair
#include <math.h> //pow
 //#include <numbers> // C++ 20
 
#include "common/util/raster.h"
 //#include "lib_disc/domain_util.h"
 //#include "lib_disc/domain.h"
 //#include "lib_disc/function_spaces/grid_function.h" //Grid_function
 //#include "lib_disc/spatial_disc/user_data/linker/linker.h"
#include "lib_disc/spatial_disc/user_data/std_glob_pos_data.h"
 
#include "lib_grid/algorithms/space_partitioning/lg_ntree.h"
#include "lib_grid/algorithms/bounding_box_util.h"
 
 namespace ug {
 /*
  * This file is part of UG4 to storage the value of each element in domain
  *
  * input:
  * domain, data, use_relative_position(true: pos=pos-box_center of elem)
  *
  *
  * output:
  * data: attached on the found element.
  *
  *
  *
  */
 
 template<typename TDomain>
 class HomoRasterUserData:public StdGlobPosData<HomoRasterUserData<TDomain>, number, TDomain::dim >
 {
     //const number PI=3.141592653589793238463;
 public:
    typedef number TData;
    static const int dim = TDomain::dim;
	typedef CplUserData<TData, dim> base_type;

     typedef typename grid_dim_traits<dim>::element_type element_type;
     //typedef typename NTreeGridData<world_dim>::position_attachment_t  position_attachment_t;
     typedef lg_ntree<dim, dim, element_type> element_tree_type;
     typedef Raster<number, dim> TRaster;
     typedef SmartPtr<TRaster> input_type;

private:
	SmartPtr<TRaster>   m_spRaster;
    //SmartPtr<TDomain> m_spDom;
	int m_interpOrder; //0: node 1: linear 2: average
    element_tree_type m_elemTree;
    bool m_withTree;

public:
    HomoRasterUserData(SmartPtr<TRaster> myRaster)
    :m_withTree(false){
        m_spRaster=myRaster;
    }

    /*HomoRasterUserData(SmartPtr<TDomain > spDom, const int level=0)
    {
        set_domain(spDom, level);
    }*/

   /* virtual ~HomoRasterUserData()
	{}*/

public:

    void set_order(int order) {m_interpOrder = order;}

    void set_domain(SmartPtr<TDomain> spDom, const int level=0){
        m_elemTree.set_grid(*spDom->grid(), spDom->position_attachment());
        m_elemTree.create_tree(spDom->grid()->template begin<element_type>(level), spDom->grid()->template end<element_type>(level));
        m_withTree=true;
    }

    public:
	inline void evaluate (number& D,
			 		      const MathVector<dim>& globIP,
			 		      number time, int si) const
	{
			if(!evaluate(D, globIP))
                  UG_THROW("It could not be founded an element containing the specified point: " << globIP);
           
	}

     inline bool evaluate(number& value, const MathVector<dim>& x) const
     {

        if(m_withTree){
        element_type* elem=NULL;
        /*if(!m_elemTree.bounding_box(0).contains_point(pos)){
            return false;
        }*/
         if(!FindContainingElement(elem, m_elemTree, x)){
            return false;
         }

         // m_interpOrder=0, how many raster_nodes in elem -> average(node_value)
         // m_interpOrder=1 -> at least one raster_node in elem :
                                // refine raster_box (max. 2-3)-> 
                                //1. all raster_nodes in the elem
                                // 2. for dim==2, at least 3 raster/for dim==3, at least 5 raster_nodes
         // m_interpOrder=2, at least one raster_node in elem -> average(box_average)

        }else{
          get_number(value, x);
        }
     }

     private:
     inline void get_number(number& y, const MathVector<dim>& x) const
	{
		typename TRaster::Coordinate xcoord(x);	
		y = m_spRaster->interpolate(xcoord, m_interpOrder);
		return;
	}


 };
} 

#endif /* _H__UG__LIB_DISC__SPATIAL_DISC__HOMORASTER_USERDATA_H_*/