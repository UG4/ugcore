/*
 * interface_handler_local_impl.h
 *
 *  Created on: 19.01.2015
 *      Author: suze
 */

#ifndef INTERFACE_HANDLER_LOCAL_DIFFUSION_DATA_H_
#define INTERFACE_HANDLER_LOCAL_DIFFUSION_DATA_H_



namespace ug{


///////////////////////////////////////////////////////////////
/// methods for class 'InterfaceHandlerLocalDiffusion'
///////////////////////////////////////////////////////////////


///////////////////////////////////////////////////////////////
/// hard coded boundary conditions
///////////////////////////////////////////////////////////////

template<int dim>
inline double get_jump_value_ex3(const MathVector<dim> position, const MathVector<dim> center, const int orientation)
{

	return 0.0;

	if ( orientation == 1)
		return 0.0;

	double absValue = position[0]*position[0] + position[1]* position[1];
	double returnValue = exp(-absValue);

	return returnValue;
}





template<int dim>
inline double get_jump_value_const(const MathVector<dim> position, const MathVector<dim> center, const int orientation)
{
	return 0.0;

	double absValue = position[0]*position[0] + position[1]* position[1];
	double factor = 8*(2*absValue - position[0] - position[1]);

	return factor*exp(-absValue);

	double sum = position[0] + position[1];

	double returnValue = log(absValue) - sin(sum);

	return returnValue;
}




template<int dim>
inline double get_jump_grad_value_kappa_Frei(const MathVector<dim> position, const MathVector<dim> center, const int orientation)
{
	return 0.0;

	if ( orientation == 1)
		return -1.0;
	else
		return -0.1;

}

template<int dim>
inline double get_jump_grad_value_kappa_Frei_inverse(const MathVector<dim> position, const MathVector<dim> center, const int orientation)
{
	if ( orientation == 1)
		return -1.0;
	else
		return -0.1;

}




template<int dim>
inline double get_source_kappa_konform(const MathVector<dim> position, const MathVector<dim> center, const int orientation)
{
	double center_x = 0.0;
	double center_y = -0.0244;

	if ( orientation == 1)
	{
		return  16.0*16.0;
	}
	else
	{
		if ( orientation != -1)
			UG_THROW("wrong orientation!\n");

		double dist_x = position[0] - center_x;
		double dist_y = position[1] - center_y;
	 	double dist = sqrt(dist_x*dist_x+dist_y*dist_y);

	 	return 200*16*dist*dist;

	}

}



template<int dim>
inline double get_source_kappa_Frei(const MathVector<dim> position, const MathVector<dim> center, const int orientation)
{
	double center_x = 0.0; //0.1;
	double center_y = 0.0; //0.2;
	double kappa_1 = 0.1;
	double kappa_2 = 1.0;

	if ( orientation == 1)
	{
		double dist_x = position[0] - center_x;
		double dist_y = position[1] - center_y;
	 	double distSq = sqrt(dist_x*dist_x+dist_y*dist_y);

		return 2.0*16*kappa_1*kappa_2*distSq;
	}
	else
	{
		if ( orientation != -1)
			UG_THROW("wrong orientation!\n");

	 	return 4.0*kappa_1*kappa_2;

	}

}

    
template<int dim>
inline double get_source_kappa_Frei_inverse(const MathVector<dim> position, const MathVector<dim> center, const int orientation)
{
	double center_x = 0.1;
	double center_y = 0.2;

	if ( orientation == 1)
	{
		return  1.0*1.0*4.0;
	}
	else
	{
		if ( orientation != -1)
			UG_THROW("wrong orientation!\n");

		double dist_x = position[0] - center_x;
		double dist_y = position[1] - center_y;
	 	double dist = sqrt(dist_x*dist_x+dist_y*dist_y);

	 	return 2.0*0.1*0.1*16*dist*dist;

	}

}






} // end namespace ug



#endif /* INTERFACE_HANDLER_LOCAL_DIFFUSION_DATA_H_ */
