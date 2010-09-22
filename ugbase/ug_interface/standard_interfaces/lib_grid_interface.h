//	created by Sebastian Reiter
//	s.b.reiter@googlemail.com
//	y10 m09 d22

#ifndef __H__UG__INTERFACE__LIB_GRID_INTERFACE__
#define __H__UG__INTERFACE__LIB_GRID_INTERFACE__

#include "../ug_interface.h"
#include "lib_grid/lib_grid.h"

namespace ug{
namespace interface
{

///	Represents a ug::Grid.
class GridObject : public ObjectBase<GridObject, IObject>
{
	public:
		static const char* static_type_name()			{return "Grid";}

		GridObject() : m_pGrid(NULL)	{}	
		~GridObject()	{if(m_pGrid) delete m_pGrid;}
		
		virtual IObject* clone()
		{
			GridObject* go = new GridObject;
			go->m_pGrid = new Grid;
			return go;
		}
		
		virtual Grid& get_grid()		{return *m_pGrid;}

		
	private:
		Grid*	m_pGrid;
};

///	Represents a ug::MultiGrid
class MultiGridObject : public ObjectBase<MultiGridObject, GridObject>
{
	public:
		static const char* static_type_name()			{return "MultiGrid";}
		
		MultiGridObject() : m_pMultiGrid(NULL)	{}
		~MultiGridObject()	{if(m_pMultiGrid) delete m_pMultiGrid;}
		
		virtual IObject* clone()
		{
			MultiGridObject* mgo = new MultiGridObject;
			mgo->m_pMultiGrid = new MultiGrid;
			return mgo;
		}
		
		virtual Grid& get_grid()				{return *m_pMultiGrid;}
		virtual MultiGrid& get_multi_grid()		{return *m_pMultiGrid;}
		
	private:
		MultiGrid* m_pMultiGrid;
};

class ISubsetHandlerObject : public AbstractObjectBase<ISubsetHandlerObject, IObject>
{
	public:
		static const char* static_type_name()	{return "ISubsetHandler";}
		
		ISubsetHandlerObject()
		{
			typedef ISubsetHandlerObject THIS;
			
			{//	num_subsets
				MethodDesc& md = add_method("num_subsets", &THIS::num_subsets);
				md.params_out().add_int("num");
			}
		}
		
		virtual ISubsetHandler& get_subset_handler_interface() = 0;
		
	protected:		
		void num_subsets(ParameterList& in, ParameterList& out)
		{
			out.set_int(0, (int)get_subset_handler_interface().num_subsets());
		}
};

///	Represents a ug::SubsetHandler
class SubsetHandlerObject : public ObjectBase<SubsetHandlerObject, ISubsetHandlerObject>
{
	public:
		static const char* static_type_name()	{return "SubsetHandler";}
		
		SubsetHandlerObject()
		{
			typedef SubsetHandlerObject THIS;
			
			{//	set_grid
				MethodDesc& md = add_method("set_grid", &THIS::set_grid);
				md.params_in().add_object("grid", "Grid");
				md.params_out().add_object("this", "SubsetHandler");
			}
		}
		
		virtual ISubsetHandler& get_subset_handler_interface()	{return m_subsetHandler;}
		virtual SubsetHandler& get_subset_handler()				{return m_subsetHandler;}
				
	protected:
		void set_grid(ParameterList& in, ParameterList& out)
		{
			GridObject* go = static_cast<GridObject*>(in.to_object(0));
			
			if(go)
				m_subsetHandler.assign_grid(go->get_grid());
			out.set_object(0, this);
		}

	protected:
		SubsetHandler m_subsetHandler;
};

///	Represents a ug::SubsetHandler
class MGSubsetHandlerObject : public ObjectBase<MGSubsetHandlerObject, ISubsetHandlerObject>
{
	public:
		static const char* static_type_name()	{return "MGSubsetHandler";}

		MGSubsetHandlerObject()
		{
			typedef MGSubsetHandlerObject THIS;
			
			{//	set_grid
				MethodDesc& md = add_method("set_grid", &THIS::set_grid);
				md.params_in().add_object("grid", "MultiGrid");
				md.params_out().add_object("this", "MGSubsetHandler");
			}
		}
		
		virtual ISubsetHandler& get_subset_handler_interface()	{return m_subsetHandler;}
		virtual MGSubsetHandler& get_subset_handler()			{return m_subsetHandler;}
		
	protected:
		void set_grid(ParameterList& in, ParameterList& out)
		{
			MultiGridObject* go = static_cast<MultiGridObject*>(in.to_object(0));
			
			if(go)
				m_subsetHandler.assign_grid(go->get_multi_grid());
			out.set_object(0, this);
		}

	protected:
		MGSubsetHandler m_subsetHandler;
};

///	Represents ug::LoadGridFromFile 
class LoadGridFunc : public IGlobalFunction
{
	public:
		LoadGridFunc()
		{
			m_method.name() = "load_grid";
			m_method.tooltip() = "loads a grid from file.";
			m_method.params_in().add_object("grid", "Grid");
			m_method.params_in().add_object("subset_handler", "ISubsetHandler");
			m_method.params_in().add_string("filename");
			m_method.params_out().add_int("success");
		}
		
	protected:
		virtual void execute(ParameterList& in, ParameterList& out)
		{
			GridObject* go = static_cast<GridObject*>(in.to_object(0));
			ISubsetHandlerObject* sho = static_cast<ISubsetHandlerObject*>(in.to_object(1));
			const char* filename = in.to_string(2);
			
			ISubsetHandler& sh = sho->get_subset_handler_interface();
			Grid& g = go->get_grid();
			
			if(sh.get_assigned_grid() != &g){
				UG_LOG("ERROR in load_grid: No grid assigned to specified subset-handler.");
			//	failure. return 0
				out.set_int(0, 0);
				return;
			}
				
		//	load grid with subset handler
			if(LoadGridFromFile(g, filename, sh)){
			//	success. return 1
				out.set_int(0, 1);
				return;
			}
			
		//	failure. return 0
			out.set_int(0, 0);
			return;
		}
};

///	Represents ug::SaveGridToFile
class SaveGridFunc : public IGlobalFunction
{
	public:
		SaveGridFunc()
		{
			m_method.name() = "save_grid";
			m_method.tooltip() = "saves a grid to a file.";
			m_method.params_in().add_object("grid", "Grid");
			m_method.params_in().add_object("subset_handler", "SubsetHandler");
			m_method.params_in().add_string("filename");
			m_method.params_out().add_int("success");
		}
		
	protected:
		virtual void execute(ParameterList& in, ParameterList& out)
		{
			GridObject* go = static_cast<GridObject*>(in.to_object(0));
			SubsetHandlerObject* sho = static_cast<SubsetHandlerObject*>(in.to_object(1));
			const char* filename = in.to_string(2);
			
			SubsetHandler& sh = sho->get_subset_handler();
			Grid& g = go->get_grid();
			
			if(sh.get_assigned_grid() != &g){
				UG_LOG("ERROR in load_grid: No grid assigned to specified subset-handler.");
			//	failure. return 0
				out.set_int(0, 0);
				return;
			}
				
		//	load grid with subset handler
			if(SaveGridToFile(g, filename, sh)){
			//	success. return 1
				out.set_int(0, 1);
				return;
			}
				
		//	failure. return 0
			out.set_int(0, 0);
			return;
		}
};

}//	end of namespace 
}//	end of namespace 

#endif
