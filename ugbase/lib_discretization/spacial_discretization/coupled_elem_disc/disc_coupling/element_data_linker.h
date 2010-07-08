
#ifndef __H__LIB_DISCRETIZATION__ELEMENT_DATA_LINKER__
#define __H__LIB_DISCRETIZATION__ELEMENT_DATA_LINKER__

#include <typeinfo>

#include "element_data_items.h"
#include "element_data_import_export.h"

namespace ug{

// dummy type
struct EmptyType {};

// typelist
template
<
  typename T1=EmptyType,
  typename T2=EmptyType,
  typename T3=EmptyType,
  typename T4=EmptyType,
  typename T5=EmptyType
> struct typelist;

// implementation of typelist
template
<
  typename T1,
  typename T2,
  typename T3,
  typename T4,
  typename T5
>
struct typelist
{
  typedef T1 head;
  typedef typelist< T2, T3, T4, T5 > tail;
  enum
  {
    length = tail::length+1
  };
};

// empty typelist specialization
template<>
struct typelist< EmptyType, EmptyType, EmptyType, EmptyType >
{
  enum
  {
    length = 0
  };
};

// provider for link function type
template <typename TDataType, typename TPositionType,  typename T0, typename T1, typename T2, typename T3, typename T4, int N = typelist<T0,T1,T2,T3,T4>::length>
struct LinkFunctionType
{
	typedef EmptyType value;
};

template <>
template <typename TDataType, typename TPos,  typename T0, typename T1, typename T2, typename T3, typename T4>
struct LinkFunctionType<TDataType, TPos, T0,T1,T2,T3,T4, 5>
{
	typedef void (*LinkValueFunction)(std::vector<TDataType>&,
										const DataImport<T0, TPos>&, const DataImport<T1, TPos>&, const DataImport<T2, TPos>&, const DataImport<T3, TPos>&, const DataImport<T4, TPos>&);

	typedef void (*LinkDerivativeFunction)(std::vector<std::vector<TDataType> >&, std::size_t glob_sys,
										const DataImport<T0, TPos>&, const DataImport<T1, TPos>&, const DataImport<T2, TPos>&, const DataImport<T3, TPos>&, const DataImport<T4, TPos>&);
};

template <>
template <typename TDataType, typename TPos,  typename T0, typename T1, typename T2, typename T3, typename T4>
struct LinkFunctionType<TDataType, TPos,  T0,T1,T2,T3,T4, 4>
{
	typedef void (*LinkValueFunction)(std::vector<TDataType>&,
										const DataImport<T0, TPos>&, const DataImport<T1, TPos>&, const DataImport<T2, TPos>&, const DataImport<T3, TPos>&);

	typedef void (*LinkDerivativeFunction)(std::vector<std::vector<TDataType> >&, std::size_t glob_sys,
										const DataImport<T0, TPos>&, const DataImport<T1, TPos>&, const DataImport<T2, TPos>&, const DataImport<T3, TPos>&);
};

template <>
template <typename TDataType, typename TPos,  typename T0, typename T1, typename T2, typename T3, typename T4>
struct LinkFunctionType<TDataType, TPos,  T0,T1,T2,T3,T4, 3>
{
	typedef void (*LinkValueFunction)(std::vector<TDataType>&,
										const DataImport<T0, TPos>&, const DataImport<T1, TPos>&, const DataImport<T2, TPos>&);

	typedef void (*LinkDerivativeFunction)(std::vector<std::vector<TDataType> >&, std::size_t glob_sys,
										const DataImport<T0, TPos>&, const DataImport<T1, TPos>&, const DataImport<T2, TPos>&);
};

template <>
template <typename TDataType, typename TPos,  typename T0, typename T1, typename T2, typename T3, typename T4>
struct LinkFunctionType<TDataType, TPos,  T0,T1,T2,T3,T4, 2>
{
	typedef void (*LinkValueFunction)(std::vector<TDataType>&,
										const DataImport<T0, TPos>&, const DataImport<T1, TPos>&);

	typedef void (*LinkDerivativeFunction)(std::vector<std::vector<TDataType> >&, std::size_t glob_sys,
										const DataImport<T0, TPos>&, const DataImport<T1, TPos>&);
};

template <>
template <typename TDataType, typename TPos,  typename T0, typename T1, typename T2, typename T3, typename T4>
struct LinkFunctionType<TDataType, TPos,  T0,T1,T2,T3,T4, 1>
{
	typedef void (*LinkValueFunction)(std::vector<TDataType>&,
										const DataImport<T0, TPos>&);

	typedef void (*LinkDerivativeFunction)(std::vector<std::vector<TDataType> >&, std::size_t glob_sys,
										const DataImport<T0, TPos>&);
};

template <>
template <typename TDataType, typename TPos,  typename T0, typename T1, typename T2, typename T3, typename T4>
struct LinkFunctionType<TDataType, TPos,  T0,T1,T2,T3,T4, 0>
{
	typedef void (*LinkValueFunction)(std::vector<TDataType>&);

	typedef void (*LinkDerivativeFunction)(std::vector<std::vector<TDataType> >&, std::size_t glob_sys);
};

// predeclaration
template <	typename TDataType, typename TPositionType,
			typename T0 = EmptyType, typename T1 = EmptyType, typename T2 = EmptyType, typename T3 = EmptyType, typename T4 = EmptyType>
class DataLinker;

template <typename TDataType, typename TPositionType,  typename T0, typename T1, typename T2, typename T3, typename T4, int N = typelist<T0,T1,T2,T3,T4>::length>
struct InvokeClass{};

template <>
template <typename TDataType, typename TPositionType,  typename T0, typename T1, typename T2, typename T3, typename T4>
struct InvokeClass<TDataType, TPositionType,  T0,T1,T2,T3,T4, 0>{public:
	void compute_values(DataLinker<TDataType, TPositionType,  T0,T1,T2,T3,T4>* me)
	{me->m_linkValueFunction(me->m_values);	}
	void compute_derivatives(DataLinker<TDataType, TPositionType,  T0,T1,T2,T3,T4>* me, std::size_t sys)
	{}
};

template <>
template <typename TDataType, typename TPositionType,  typename T0, typename T1, typename T2, typename T3, typename T4>
struct InvokeClass<TDataType, TPositionType,  T0,T1,T2,T3,T4, 1>{public:
	void compute_values(DataLinker<TDataType, TPositionType,  T0,T1,T2,T3,T4>* me)
	{me->m_linkValueFunction(me->m_values, me->m_import0);	}
	void compute_derivatives(DataLinker<TDataType, TPositionType,  T0,T1,T2,T3,T4>* me, std::size_t sys)
	{me->m_linkDerivativeFunction((me->m_derivatives)[sys], sys, me->m_import0);	}
};

template <>
template <typename TDataType, typename TPositionType,  typename T0, typename T1, typename T2, typename T3, typename T4>
struct InvokeClass<TDataType, TPositionType,  T0,T1,T2,T3,T4, 2>{public:
	void compute_values(DataLinker<TDataType, TPositionType,  T0,T1,T2,T3,T4>* me)
	{me->m_linkValueFunction(me->m_values, me->m_import0, me->m_import1);	}
	void compute_derivatives(DataLinker<TDataType, TPositionType,  T0,T1,T2,T3,T4>* me, std::size_t sys)
	{me->m_linkDerivativeFunction((me->m_derivatives)[sys], sys, me->m_import0, me->m_import1);};
};

template <>
template <typename TDataType, typename TPositionType,  typename T0, typename T1, typename T2, typename T3, typename T4>
struct InvokeClass<TDataType, TPositionType,  T0,T1,T2,T3,T4, 3>{public:
	void compute_values(DataLinker<TDataType, TPositionType,  T0,T1,T2,T3,T4>* me)
	{me->m_linkValueFunction(me->m_values, me->m_import0, me->m_import1, me->m_import2);	}
	void compute_derivatives(DataLinker<TDataType, TPositionType,  T0,T1,T2,T3,T4>* me, std::size_t sys)
	{me->m_linkDerivativeFunction((me->m_derivatives)[sys], sys, me->m_import0, me->m_import1, me->m_import2);};
};

template <>
template <typename TDataType, typename TPositionType,  typename T0, typename T1, typename T2, typename T3, typename T4>
struct InvokeClass<TDataType, TPositionType,  T0,T1,T2,T3,T4, 4>{public:
	void compute_values(DataLinker<TDataType, TPositionType,  T0,T1,T2,T3,T4>* me)
	{me->m_linkValueFunction(me->m_values, me->m_import0, me->m_import1, me->m_import2, me->m_import3);	}
	void compute_derivatives(DataLinker<TDataType, TPositionType,  T0,T1,T2,T3,T4>* me, std::size_t sys)
	{me->m_linkDerivativeFunction((me->m_derivatives)[sys], sys, me->m_import0, me->m_import1, me->m_import2, me->m_import3);};
};

template <>
template <typename TDataType, typename TPositionType,  typename T0, typename T1, typename T2, typename T3, typename T4>
struct InvokeClass<TDataType, TPositionType,  T0,T1,T2,T3,T4, 5>{public:
	void compute_values(DataLinker<TDataType, TPositionType,  T0,T1,T2,T3,T4>* me)
	{me->m_linkValueFunction(me->m_values, me->m_import0, me->m_import1, me->m_import2, me->m_import3, me->m_import4);	}
	void compute_derivatives(DataLinker<TDataType, TPositionType,  T0,T1,T2,T3,T4>* me, std::size_t sys)
	{me->m_linkDerivativeFunction((me->m_derivatives)[sys], sys, me->m_import0, me->m_import1, me->m_import2, me->m_import3, me->m_import4); };
};



template <	typename TDataType, typename TPositionType,
			typename T0 = EmptyType, typename T1 = EmptyType, typename T2 = EmptyType, typename T3 = EmptyType, typename T4 = EmptyType>
class DataLinkerPossibility : public DataPossibilityItem {

	public:
		typedef TDataType out_data_type;
		typedef TPositionType position_type;

	protected:
		typedef typelist<T0, T1, T2, T3, T4> InTypeList;

		typedef typename LinkFunctionType<TDataType, TPositionType,  T0, T1, T2, T3, T4, InTypeList::length>::LinkValueFunction LinkValueFunction;
		typedef typename LinkFunctionType<TDataType, TPositionType,  T0, T1, T2, T3, T4, InTypeList::length>::LinkDerivativeFunction LinkDerivativeFunction;

	public:
		DataLinkerPossibility(std::string name, LinkValueFunction valueFunc, LinkDerivativeFunction derivFunc) :
			DataPossibilityItem(name, InTypeList::length, &typeid(TDataType), &typeid(TPositionType)),
			m_valueFunction(valueFunc), m_derivFunction(derivFunc)
		{};

		DataExportItem* create_data_export()
		{
			DataLinker<TDataType, TPositionType,  T0, T1, T2, T3, T4>* linker =
				new DataLinker<TDataType, TPositionType,  T0, T1, T2, T3, T4>(this->name(), this, m_valueFunction, m_derivFunction);

			m_createdDataExports.push_back(dynamic_cast<DataExportItem*>(linker));
			return dynamic_cast<DataExportItem*>(linker);
		}

		bool link(DataPossibilityItem* posItem, std::size_t slot)
		{
			if(!(slot < InTypeList::length))
			{
				UG_LOG("DataLinkerPossibility::clear_slot(): Slot "<< slot << " does not exist.\n");
				return false;
			}

			// no linking if already linked
			if(is_linked(slot) == true)
			{
				UG_LOG("DataLinkerPossibility::link(): Slot "<< slot << " already linked. Cannot link another possibility to this slot.\n");
				return false;
			}

			// no linking if posItem is not instantiable (currently, maybe relaxe this condition in the future)
			if(posItem->is_instantiable() != true)
			{
				UG_LOG("DataLinkerPossibility::link(): Possibility Item that should be linked is not instantiable. Cannot link.\n");
				return false;
			}

			if(posItem->position_type() != dynamic_cast<DataPossibilityItem*>(this)->position_type())
			{
				UG_LOG("DataLinkerPossibility::link(): Position type does not match. Cannot link.\n");
				return false;
			}

			if(posItem->data_type() != this->data_type(slot))
			{
				UG_LOG("DataLinkerPossibility::link(): Data Type of slot " << slot << " does not match data type of Possibility Item. Cannot link.\n");
				return false;
			}

			this->m_slotPosItems[slot] = posItem;

			LinkedPossibility linkedPos; linkedPos.item = this; linkedPos.slot = slot;
			posItem->add_linked_possibility(linkedPos);

			return true;
		}

		const LinkValueFunction get_eval_function() const {return m_valueFunction;};

	protected:
		const std::type_info* data_type(std::size_t slot)
		{
			switch(slot)
			{
			case 0: return &typeid(T0);
			case 1: return &typeid(T1);
			case 2: return &typeid(T2);
			case 3: return &typeid(T3);
			case 4: return &typeid(T4);
			default: UG_ASSERT(slot < InTypeList::length, "Cannot determine data type. Aborting.\n");
			}
			return NULL;
		}

	protected:
		LinkValueFunction m_valueFunction;
		LinkDerivativeFunction m_derivFunction;
};


template <	typename TDataType, typename TPositionType,
			typename T0, typename T1, typename T2, typename T3, typename T4>
class DataLinker : public DataExport<TDataType, TPositionType> {
	friend class DataImport<TDataType,TPositionType>;
	friend class DataContainer;

	public:
		typedef TDataType out_data_type;
		typedef TPositionType position_type;

	protected:
		typedef typelist<T0, T1, T2, T3, T4> InTypeList;

		typedef typename LinkFunctionType<TDataType, TPositionType,  T0, T1, T2, T3, T4, InTypeList::length>::LinkValueFunction LinkValueFunction;
		typedef typename LinkFunctionType<TDataType, TPositionType,  T0, T1, T2, T3, T4, InTypeList::length>::LinkDerivativeFunction LinkDerivativeFunction;

	protected:

		friend class InvokeClass<TDataType, TPositionType,  T0, T1, T2, T3, T4, InTypeList::length>;
		InvokeClass<TDataType, TPositionType,  T0, T1, T2, T3, T4, InTypeList::length> m_Invoke;

	public:
		DataLinker(std::string name, DataPossibilityItem* possibility, LinkValueFunction func, LinkDerivativeFunction derivFunc) :
			DataExport<TDataType, TPositionType>(name, possibility),
			m_linkValueFunction(func), m_linkDerivativeFunction(derivFunc),
			m_import0(name + ": Slot 0"),
			m_import1(name + ": Slot 1"),
			m_import2(name + ": Slot 2"),
			m_import3(name + ": Slot 3"),
			m_import4(name + ": Slot 4")
		{
			UG_DLOG(LIB_DISC_LINKER, 2, "DataLinker::DataLinker: Creating Data Linker.\n");
			m_impItem[0] = dynamic_cast<DataImportItem*>(&m_import0);
			m_impItem[1] = dynamic_cast<DataImportItem*>(&m_import1);
			m_impItem[2] = dynamic_cast<DataImportItem*>(&m_import2);
			m_impItem[3] = dynamic_cast<DataImportItem*>(&m_import3);
			m_impItem[4] = dynamic_cast<DataImportItem*>(&m_import4);
		}

		virtual void compute(bool computeDerivatives)
		{
			UG_DLOG(LIB_DISC_LINKER, 2, "DataLinker::compute: Compute values.\n");
			m_Invoke.compute_values(this);

			if(computeDerivatives)
			{
				UG_DLOG(LIB_DISC_LINKER, 2, "DataLinker::compute: Compute Derivatives w.r.t " << this->num_sys() << " system(s).\n");
				for(std::size_t sys = 0; sys < this->num_sys(); ++sys)
				{
					UG_DLOG(LIB_DISC_LINKER, 2, "DataLinker::compute: Compute derivatives w.r.t system " << this->sys(sys) << ".\n");
					m_Invoke.compute_derivatives(this, sys);
				}
			}
		}

		// ALL functions derived publically from export

		bool reset_derivative_array()
		{
			if(!is_linked()) return true;
			UG_DLOG(LIB_DISC_LINKER, 2, "DataLinker::reset_derivative_array: Reset derivative array.\n");

			// reset
			this->m_num_sys = 0;
			this->m_sys.clear();
			this->m_num_sh.clear();

			UG_DLOG(LIB_DISC_LINKER, 2, "DataLinker::reset_derivative_array: Checking " << num_slots() << " own imports.\n");
			// get all systems this linker depends on and remember num_sh of each sys
			for(std::size_t i = 0; i < num_slots();  ++i)
			{
				std::size_t num_sys = m_impItem[i]->num_sys();
				UG_DLOG(LIB_DISC_LINKER, 2, "DataLinker::reset_derivative_array: Checking import " << i << ", num_sys = " << num_sys << ".\n");
				for(std::size_t s = 0; s < num_sys; ++s)
				{
					std::size_t t = m_impItem[i]->sys(s);
					std::size_t num_sh = m_impItem[i]->num_sh(s);
					UG_DLOG(LIB_DISC_LINKER, 2, "DataLinker::reset_derivative_array: Adding system " << t << " with " << num_sh << " unknowns (if not already added).\n");

					std::vector<std::size_t>::iterator it, it_num;

					if(this->m_sys.empty())
					{
						UG_DLOG(LIB_DISC_LINKER, 2, "DataLinker::reset_derivative_array: Start Dependency list.\n");
						this->m_num_sys++;
						this->m_sys.push_back(t);
						this->m_num_sh.push_back(num_sh);
					}
					else
					{
						for(it = this->m_sys.begin(), it_num = this->m_num_sh.begin(); ; ++it, ++it_num)
						{
							// insert before greater elements
							if(*it > t)
							{
								UG_DLOG(LIB_DISC_LINKER, 2, "DataLinker::reset_derivative_array: Insert into Dependency list.\n");
								this->m_num_sys++;
								this->m_sys.insert(it, t);
								this->m_num_sh.insert(it_num, num_sh);
								break;
							}

							// if end reached
							if(it == this->m_sys.end())
							{
								UG_DLOG(LIB_DISC_LINKER, 2, "DataLinker::reset_derivative_array: Insert at end of Dependency list.\n");
								this->m_num_sys++;
								this->m_sys.insert(it, t);
								this->m_num_sh.insert(it_num, num_sh);
								break;
							}

							// break if already in list
							if((*it == t)) break;
						}
					}
				}
			}

			m_max_sh = 0;
			this->m_values.resize(this->num_ip());
			this->m_derivatives.resize(this->num_sys());
			for(std::size_t s = 0; s < this->num_sys(); ++s)
			{
				m_max_sh = (this->num_sh(s) > m_max_sh) ? this->num_sh(s) : m_max_sh;
				this->m_derivatives[s].resize(this->num_ip());
				for(std::size_t ip = 0; ip < this->m_derivatives[s].size(); ++ip)
				{
					this->m_derivatives[s][ip].resize(this->num_sh(s));
				}
			}

			return true;
		}

		// number of data exports linked by this linker
		virtual std::size_t num_slots() const {	return InTypeList::length;}

		virtual std::string slot_name(std::size_t slot) const{ return m_impItem[slot]->name();}

		// add a Data Export number i
		virtual bool link(DataExportItem* Export, std::size_t slot)
		{
			UG_ASSERT(slot < num_slots(), "Access to slot 'slot', but linker has not so many slots.");

			if(m_impItem[slot]->link_data_export(Export) != true)
			{
				UG_LOG("DataLinker::link_data_export: Error while linking.\n");
				return false;
			}

			if(reset_derivative_array() != true)
			{
				UG_LOG("DataLinker::link_data_export: Error while resetting derivative array.\n");
				return false;
			}

			return true;
		}

		// remove Data Export number i
		virtual bool clear_slot(std::size_t slot)
		{
			UG_ASSERT(slot < num_slots(), "Access to slot 'slot', but linker has not so many slots.");
			return m_impItem[slot]->clear_data_export();
		}

		// return if an export is set at slot i
		virtual bool is_linked(std::size_t slot) const
		{
			UG_ASSERT(slot < num_slots(), "Access to slot 'slot', but linker has not so many slots.");
			return m_impItem[slot]->is_linked();
		}

		// return if all exports are set
		virtual bool is_linked() const
		{
			for(std::size_t i = 0; i < num_slots(); ++i)
			{
				if(is_linked(i) == false) return false;
			}
			return true;
		}

		// help function: is called to set values in linker imports as well
		virtual bool set_linker_positions(const std::vector<position_type>& positions)
		{
			for(std::size_t i = 0; i < num_slots(); ++i)
			{
				dynamic_cast<DataImportPosition<TPositionType>*>(m_impItem[i])->set_positions(positions);
			}
			return true;
		};

		// get registered export of slot i
		virtual const DataExportItem* get_data_export(std::size_t slot) const
		{
			UG_ASSERT(slot < num_slots(), "Access to slot i, but linker has not so many slots.");

			return m_impItem[slot]->get_data_export();
		}

	protected:
		LinkValueFunction m_linkValueFunction;
		LinkDerivativeFunction m_linkDerivativeFunction;

		std::size_t m_max_sh;

		DataImportItem* m_impItem[5];

		// Data Export, to which import slot is linked
		DataImport<T0, TPositionType> m_import0;
		DataImport<T1, TPositionType> m_import1;
		DataImport<T2, TPositionType> m_import2;
		DataImport<T3, TPositionType> m_import3;
		DataImport<T4, TPositionType> m_import4;
};

// print current values of import to ostream
std::ostream& operator<<(std::ostream& out, const EmptyType& data);

}

#endif /*__H__LIB_DISCRETIZATION__ELEMENT_DATA_LINKER__*/
