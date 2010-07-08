/*
 * elementdata.h
 *
 *  Created on: 30.06.2009
 *      Author: andreasvogel
 */

#ifndef __H__LIB_DISCRETIZATION__ELEMENTDATA__
#define __H__LIB_DISCRETIZATION__ELEMENTDATA__

#include <typeinfo>
#include <string>
#include <vector>

namespace ug{

// predeclaration
template <typename TDataType, typename TPositionType> class DataImport;
template <typename TDataType, typename TPositionType> class DataExport;


class DataExportItem {
	friend class DataContainer;

	public:
		DataExportItem(std::string name);
		std::string name();
		const std::type_info* DataType() const;
		const std::type_info* PositionType() const;

		virtual void compute() = 0;
		virtual ~DataExportItem(){};

	protected:
		std::string _name;
		const std::type_info* _DataType;
		const std::type_info* _PositionType;
};

class DataImportItem {
	friend class DataContainer;

	public:
		DataImportItem(std::string name);
		std::string name();
		const std::type_info* DataType() const;
		const std::type_info* PositionType() const;

		virtual ~DataImportItem(){};

	protected:
		virtual bool registerDataExport(DataExportItem* Export) = 0;

	protected:
		std::string _name;
		const std::type_info* _DataType;
		const std::type_info* _PositionType;
};


template <typename TDataType, typename TPositionType>
class DataExport : DataExportItem{
	friend class DataImport<TDataType,TPositionType>;
	friend class DataContainer;

	public:
		typedef TDataType data_type;
		typedef TPositionType position_type;

	protected:
		struct DataImportInfo
		{
			DataImport<TDataType, TPositionType>* import;
			int n;
			TDataType* valueList;
			TPositionType* positionList;
		};
		typedef void (*EvalFunction)(TDataType*, TPositionType*, int);

	public:
		DataExport(std::string name, EvalFunction func);
		void compute();
		~DataExport();

	protected:
		DataExport(std::string name) : DataExportItem(name){};

		bool registerDataImport(DataImport<TDataType, TPositionType>* Import);
		bool unregisterDataImport(DataImport<TDataType, TPositionType>* Import);
		bool updateDataImport(DataImport<TDataType, TPositionType>* Import);

		void writeDataToImport();

	protected:
		typedef std::vector< DataImportInfo > ImportList;

	protected:
		ImportList _ImportList;
		EvalFunction _evalFunction;
};

struct EmptyType
{

};

template
<
  typename T1=EmptyType,
  typename T2=EmptyType,
  typename T3=EmptyType,
  typename T4=EmptyType,
  typename T5=EmptyType
> struct typelist;

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

template<>
struct typelist< EmptyType, EmptyType, EmptyType, EmptyType >
{
  enum
  {
    length = 0
  };
};


template <	typename TDataType, typename TPositionType,
			typename TInDataType1 = EmptyType, typename TInDataType2 = EmptyType, typename TInDataType3 = EmptyType, typename TInDataType4 = EmptyType, typename TInDataType5 = EmptyType>
class DataLinker : public DataExport<TDataType, TPositionType> {
	friend class DataImport<TDataType,TPositionType>;
	friend class DataContainer;

	public:
		typedef TDataType out_data_type;
		typedef TPositionType position_type;
		typedef DataLinker<TDataType,TPositionType,TInDataType1,TInDataType2,TInDataType3,TInDataType4, TInDataType5> my_type;

	protected:
		typedef typelist<TInDataType1,TInDataType2,TInDataType3,TInDataType4, TInDataType5> InTypeList;

		typedef void (*LinkFunction0)(TDataType*, TPositionType*, int);
		typedef void (*LinkFunction1)(TDataType*, TPositionType*, TInDataType1*, int);
		typedef void (*LinkFunction2)(TDataType*, TPositionType*, TInDataType1*, TInDataType2*, int);
		typedef void (*LinkFunction3)(TDataType*, TPositionType*, TInDataType1*, TInDataType2*, TInDataType3*, int);
		typedef void (*LinkFunction4)(TDataType*, TPositionType*, TInDataType1*, TInDataType2*, TInDataType3*, TInDataType4*, int);
		typedef void (*LinkFunction5)(TDataType*, TPositionType*, TInDataType1*, TInDataType2*, TInDataType3*, TInDataType4*, TInDataType5*, int);

		template<int N, typename T0,typename T1,typename T2,typename T3,typename T4,typename T5> struct SelectFctType { typedef T0 Type; };
		template<typename T0,typename T1,typename T2,typename T3,typename T4,typename T5> struct SelectFctType<1,T0,T1,T2,T3,T4,T5> { typedef T1 Type; };
		template<typename T0,typename T1,typename T2,typename T3,typename T4,typename T5> struct SelectFctType<2,T0,T1,T2,T3,T4,T5> { typedef T2 Type; };
		template<typename T0,typename T1,typename T2,typename T3,typename T4,typename T5> struct SelectFctType<3,T0,T1,T2,T3,T4,T5> { typedef T3 Type; };
		template<typename T0,typename T1,typename T2,typename T3,typename T4,typename T5> struct SelectFctType<4,T0,T1,T2,T3,T4,T5> { typedef T4 Type; };
		template<typename T0,typename T1,typename T2,typename T3,typename T4,typename T5> struct SelectFctType<5,T0,T1,T2,T3,T4,T5> { typedef T5 Type; };

		typedef typename SelectFctType< InTypeList::length, LinkFunction0, LinkFunction1, LinkFunction2, LinkFunction3, LinkFunction4, LinkFunction5 >::Type LinkFunctionType;

		template<int N, typename T1,typename T2,typename T3,typename T4,typename T5> class LinkFunctionWrapper {
		public:
			inline void invoke(my_type* me, LinkFunctionType fct, TDataType* out, TPositionType* pos, int n){};
		};

		template<typename T1,typename T2,typename T3,typename T4,typename T5> class LinkFunctionWrapper<0,T1,T2,T3,T4,T5> {
		public:
			inline void invoke(my_type* me, LinkFunctionType fct, TDataType* out, TPositionType* pos, int n)
			{
				(*fct)(out, pos, n);
			};
		};
		template<typename T1,typename T2,typename T3,typename T4,typename T5> class LinkFunctionWrapper<1,T1,T2,T3,T4,T5> {
		public:
			inline void invoke(my_type* me, LinkFunctionType fct, TDataType* out, TPositionType* pos, int n)
			{
				(*fct)(out, pos, me->_valueList1, n);
			};
		};
		template<typename T1,typename T2,typename T3,typename T4,typename T5> class LinkFunctionWrapper<2,T1,T2,T3,T4,T5> {
		public:
			inline void invoke(my_type* me, LinkFunctionType fct, TDataType* out, TPositionType* pos, int n)
			{
				(*fct)(out, pos, me->_valueList1, me->_valueList2, n);
			};
		};
		template<typename T1,typename T2,typename T3,typename T4,typename T5> class LinkFunctionWrapper<3,T1,T2,T3,T4,T5> {
		public:
			inline void invoke(my_type* me, LinkFunctionType fct, TDataType* out, TPositionType* pos, int n)
			{
				(*fct)(out, pos, me->_valueList1, me->_valueList2, me->_valueList3, n);
			};
		};
		template<typename T1,typename T2,typename T3,typename T4,typename T5> class LinkFunctionWrapper<4,T1,T2,T3,T4,T5> {
		public:
			inline void invoke(my_type* me, LinkFunctionType fct, TDataType* out, TPositionType* pos, int n)
			{
				(*fct)(out, pos, me->_valueList1, me->_valueList2, me->_valueList3, me->_valueList4, n);
			};
		};
		template<typename T1,typename T2,typename T3,typename T4,typename T5> class LinkFunctionWrapper<5,T1,T2,T3,T4,T5> {
		public:
			inline void invoke(my_type* me, LinkFunctionType fct, TDataType* out, TPositionType* pos, int n)
			{
				(*fct)(out, pos, me->_valueList1, me->_valueList2, me->_valueList3, me->_valueList4, me->_valueList5, n);
			};
		};

	public:
		DataLinker(std::string name, LinkFunctionType func) : DataExport<TDataType, TPositionType>(name)
		{
			_LinkFunction = func;
		}

		void compute()
		{
			//compute values
			for(int i = 0; i < this->_ImportList.size(); ++i)
			{
				_LinkFunctionWrapper.invoke(this, _LinkFunction, this->_ImportList[i].valueList, this->_ImportList[i].positionList, this->_ImportList[i].n);
			}

			// write Data to connected Imports
			this->writeDataToImport();
		}


	protected:
		LinkFunctionType _LinkFunction;
		LinkFunctionWrapper<InTypeList::length, TInDataType1, TInDataType2, TInDataType3, TInDataType4, TInDataType5> _LinkFunctionWrapper;

		int _n;
		TInDataType1* _valueList1;
		TInDataType2* _valueList2;
		TInDataType3* _valueList3;
		TInDataType4* _valueList4;
		TInDataType5* _valueList5;
		TPositionType* _positionList;
		DataExport<TInDataType1*, TPositionType>* _Export1;
		DataExport<TInDataType2*, TPositionType>* _Export2;
		DataExport<TInDataType2*, TPositionType>* _Export3;
		DataExport<TInDataType2*, TPositionType>* _Export4;
		DataExport<TInDataType2*, TPositionType>* _Export5;


};


template <typename TDataType, typename TPositionType>
class DataImport :public DataImportItem {
	friend class DataExport<TDataType,TPositionType>;
	friend class DataContainer;

	public:
		typedef TDataType data_type;
		typedef TPositionType position_type;

	public:
		DataImport(std::string name);

		inline TDataType& operator[](size_t i);
		inline const TDataType& operator[](size_t i) const;

		void setPositions(TPositionType* posList, int n);

		bool registerDataExport(DataExport<TDataType, TPositionType>* Export);
		bool unregisterDataExport();

		~DataImport();

	protected:
		bool registerDataExport(DataExportItem *Export);

	protected:
		int _n;
		TDataType* _valueList;
		TPositionType* _positionList;
		DataExport<TDataType, TPositionType>* _Export;
};


class DataContainer {
	public:
		template <typename TDataType, typename TPositionType>
		bool registerItem(DataImport<TDataType, TPositionType>& ImportItem);

		template <typename TDataType, typename TPositionType>
		bool registerItem(DataExport<TDataType, TPositionType>& ExportItem);

		template <typename TDataType, typename TPositionType>
		bool registerItem(DataLinker<TDataType, TPositionType>& LinkerItem);

		template <typename TDataType, typename TPositionType>
		bool link(DataImport<TDataType, TPositionType>& ImportItem, DataExport<TDataType, TPositionType>& ExportItem);

		bool linkByNumber(int ImportNr, int ExportNr);

		bool compute();

		bool printInfo();

	protected:
		typedef std::vector<DataImportItem*> ImportItemList;
		typedef std::vector<DataExportItem*> ExportItemList;
		typedef std::vector<DataExportItem*> LinkerItemList;

	protected:
		ImportItemList _importItemList;

		ExportItemList _exportItemList;
		ExportItemList _neddedExportItemList;

		LinkerItemList _linkerItemList;
};

};

#include "elementdata_impl.h"

#endif /* __H__LIB_DISCRETIZATION__ELEMENTDATA__ */
