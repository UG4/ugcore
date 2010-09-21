//	created by Sebastian Reiter
//	s.b.reiter@googlemail.com
//	y10 m09 d20

#include "../ug_interface.h"
#include "../registry.h"
#include "lib_grid/lib_grid.h"

using namespace std;

namespace ug{
namespace interface
{

class GridObject : public IObject
{
	public:
		virtual IObject* clone()
		{
			GridObject* go = new GridObject;
			go->m_pGrid = new Grid;
			return go;
		}
		
		virtual const char* type_name()			{return "Grid";}
		virtual const char* supported_types()	{return ":Grid:";}

		virtual Grid& get_grid()		{return *m_pGrid;}
	
	protected:
		Grid*	m_pGrid;
};

class MultiGridObject : public GridObject
{
	public:
		virtual IObject* clone()
		{
			MultiGridObject* mgo = new MultiGridObject;
			mgo->m_pMultiGrid = new MultiGrid;
			mgo->m_pGrid = mgo->m_pMultiGrid;
			return mgo;
		}
		
		virtual const char* type_name()			{return "MultiGrid";}
		virtual const char* supported_types()	{return ":Grid:MultiGrid:";}
		
		virtual MultiGrid& get_multigrid()		{return *m_pMultiGrid;}
		
	protected:
		MultiGrid* m_pMultiGrid;
};

class SubsetHandlerObject : public IObject
{
	public:
		SubsetHandlerObject()
		{
		//	we need a set_grid method
			MethodDesc& md0 = add_method("set_grid");
			md0.params_in().add_object("grid", ":Grid:");
			md0.params_out().add_object("this", ":SubsetHandler:");
			
		//	clone
			MethodDesc& md1 = add_method("clone");
			md1.params_out().add_object("clone", ":SubsetHandler:");
		}
		
		virtual IObject* clone()				{return new SubsetHandlerObject;}
		
		virtual const char* type_name()			{return "SubsetHandler";}
		virtual const char* supported_types()	{return ":SubsetHandler:";}
		
		virtual SubsetHandler& get_subset_handler()		{return m_subsetHandler;}
		
	protected:
		virtual void execute(size_t methodIndex,
						 ParameterList& in, ParameterList& out)
		{
			switch(methodIndex){
				case 0:{ //	set_grid
					GridObject* go = static_cast<GridObject*>(in.to_object(0));
					if(go)
						m_subsetHandler.assign_grid(go->get_grid());
					out.set_object(0, this);
				}break;
				
				case 1:{ //	clone
					out.set_object(0, clone());
				}break;
			}
		}
		
	protected:
		SubsetHandler m_subsetHandler;
};

class LoadGridFunc : public IGlobalFunction
{
	public:
		LoadGridFunc()
		{
			m_method.name() = "load_grid";
			m_method.tooltip() = "loads a grid from file.";
			m_method.params_in().add_object("grid", ":Grid:");
			m_method.params_in().add_object("subset_handler", ":SubsetHandler:");
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
			
			if(sh.get_assigned_grid() != &g)
				sh.assign_grid(g);
				
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

class SaveGridFunc : public IGlobalFunction
{
	public:
		SaveGridFunc()
		{
			m_method.name() = "save_grid";
			m_method.tooltip() = "saves a grid to a file.";
			m_method.params_in().add_object("grid", ":Grid:");
			m_method.params_in().add_object("subset_handler", ":SubsetHandler:");
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
			
			if(sh.get_assigned_grid() != &g)
				sh.assign_grid(g);
				
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

void RegisterLibGridInterface(Registry& reg)
{
	reg.register_object<GridObject>();
	reg.register_object<MultiGridObject>();
	reg.register_object<SubsetHandlerObject>();
	
	reg.register_global_function<LoadGridFunc>();
	reg.register_global_function<SaveGridFunc>();
}

}//	end of namespace 
}//	end of namespace 
