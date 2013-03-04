#include "bridge/bridge.h"
#include "bridge/util.h"
#include "bridge/suffix_tag.h"
#include "type_converter.h"
#include "common/common.h"
#include "lib_disc/spatial_disc/user_data/const_user_data.h"
#include "bindings/lua/lua_util.h"
#include "bindings_vrl.h"
#include "type_converter.h"
#include "threading.h"
#include <iostream>
#include <sstream>
#include <boost/function.hpp>
#include <jni.h>
#include "lib_disc/spatial_disc/user_data/user_data.h"
#include "lib_disc/spatial_disc/user_data/std/std_pos_data.h"
#include "lib_disc/spatial_disc/user_data/std/std_linker_data.h"
#include "lib_disc/spatial_disc/user_data/data_linker.h"
#include "lib_disc/spatial_disc/user_data/data_linker_traits.h"

namespace ug {
namespace vrl {

template <std::size_t dim>
inline jdoubleArray ConvertParametersToJava(JNIEnv *env, const MathVector<dim>& x,
                                            const number time)
{
	jdoubleArray array = NULL;
	array = env->NewDoubleArray(dim+1);

	jdouble vTmp[dim + 1];
	for (size_t i = 0; i < dim; i++) vTmp[i] = x[i];
	vTmp[dim] = time;

	env->SetDoubleArrayRegion(array, 0, dim+1, vTmp);
	checkException(env, "CopyParametersToJava: error");

	return array;
}

template <typename TData>
struct vrl_traits;

template <>
struct vrl_traits<number>
{
	typedef jdouble jType;
	static const unsigned int retArrayDim = 0;
	static const int size = 1;
	static std::string name() {return "Number";}
	static std::string errMsg()
	{
		std::stringstream ss;
		ss << "<font color=\"red\">User Data (Number): Invocation error: </font>\n";
		return ss.str();
	}
	static std::string callSignature() {return "([DI)D";}
	static std::string signature() {return "D";}

	static void push(jdouble* vStorage, const number x)
	{
		vStorage[0] = x;
	}

	static void toJava(JNIEnv *env, jdouble& res, const number& x)
	{
		res = x;
	}

	static void toC(JNIEnv *env, number& res, jdouble& from)
	{
		res = from;
	}

	static void call(JNIEnv *env, number& res,
	                 jobject obj, jmethodID method,
	                 jdoubleArray params, jint si)
	{
		res = env->CallDoubleMethod(obj, method, params, si);
		if (checkException(env, errMsg())) {}
	}
};

template <std::size_t dim>
struct vrl_traits<ug::MathVector<dim> >
{
	typedef jdoubleArray jType;
	static const unsigned int retArrayDim = 1;
	static const int size = dim;
	static std::string name() {return "Vector";}
	static std::string errMsg()
	{
		std::stringstream ss;
		ss << "<font color=\"red\">User Data (Vector): Invocation error: </font>\n";
		return ss.str();
	}
	static std::string callSignature() {return "([DI)[D";}
	static std::string signature() {return "[D";}

	static void push(jdouble* vStorage, const MathVector<dim>& x)
	{
		for(size_t i = 0; i < dim; ++i)
			vStorage[i] = x[i];
	}

	static void toJava(JNIEnv *env, jdoubleArray& res, const MathVector<dim>& x)
	{
		jdouble vTmp[size];
		for (jint i = 0; i < dim; i++) vTmp[i] = x[i];
		env->SetDoubleArrayRegion(res, 0, size, vTmp);
	}

	static void toC(JNIEnv *env, MathVector<dim>& res, jdoubleArray& array)
	{
		jint arrayLength = env->GetArrayLength(array);
		if (arrayLength != size)
			UG_THROW(RED_BEGIN <<
			         "User Data (Vector): Wrong return value in Code. "
			         "Required: Vector of size "<<size<<" (passed size: " <<
			         arrayLength << ")" << COLOR_END << std::endl);

		jdouble vTmp[size];
		env->GetDoubleArrayRegion(array, 0, size, vTmp);
		for (jint i = 0; i < (jint)dim; ++i) res[i] = vTmp[i];
	}

	static void call(JNIEnv *env, MathVector<dim>& res,
	                 jobject obj, jmethodID method,
	                 jdoubleArray params, jint si)
	{
		jdoubleArray tmp =
				(jdoubleArray) env->CallObjectMethod(obj, method, params, si);

		if(checkException(env, errMsg()))
		{
			toC(env, res, tmp);
		}
	}
};

template <std::size_t dim>
struct vrl_traits<ug::MathMatrix<dim,dim> >
{
	typedef jobjectArray jType;
	static const unsigned int retArrayDim = 2;
	static const int size = dim*dim;
	static std::string name() {return "Matrix";}
	static std::string errMsg()
	{
		std::stringstream ss;
		ss << "<font color=\"red\">User Data (Matrix): Invocation error: </font>\n";
		return ss.str();
	}
	static std::string callSignature() {return "([DI)[[D";}
	static std::string signature() {return "[[D";}

	static void push(jdouble* vStorage, const MathMatrix<dim, dim>& D)
	{
		for(size_t i = 0; i < dim; ++i){
			for(size_t j = 0; j < dim; ++j){
				vStorage[i*dim+j] = D[i][j];
			}
		}
	}

	static void toJava(JNIEnv *env, jobjectArray& res, const MathMatrix<dim,dim>& x)
	{
		UG_THROW("Not implemented.");
	}

	static void toC(JNIEnv *env, MathMatrix<dim,dim>& mat, jobjectArray& array)
	{
		const int rowSize = env->GetArrayLength(array);
		if(rowSize != dim)
			UG_THROW(RED_BEGIN <<
			         "User Data (Matrix): Wrong return value in Code. "
			         "Required: Matrix of size "<<dim<<"x"<<dim<<
			         " (passed row size: "<<rowSize<< ")" << COLOR_END << std::endl);

		jdouble rowEntrys[dim];

		for(int i=0; i < (int)dim; ++i){
			jdoubleArray row = (jdoubleArray)env->GetObjectArrayElement(array, i);

			const int colSize = env->GetArrayLength(row);
			if(colSize != dim)
				UG_THROW(RED_BEGIN
				         "User Data (Matrix): Wrong return value in Code. "
				         "Required: Matrix of size "<<dim<<"x"<<dim<<
				         " (passed column size: "<<colSize<< ")" << COLOR_END << std::endl);

			env->GetDoubleArrayRegion(row, 0, dim, rowEntrys);
			for(int j=0; j < (int)dim; ++j)
				mat[i][j]= rowEntrys[j];

			// alternative is to force garbageCollector to "pin" or copy internally
//			jdouble* rowEntrys = env->GetDoubleArrayElements(row, 0);
//			env->ReleaseDoubleArrayElements(row, rowEntrys, 0);
		}
	}

	static void call(JNIEnv *env, MathMatrix<dim,dim>& res,
	                 jobject obj, jmethodID method,
	                 jdoubleArray params, jint si)
	{
		jobjectArray tmp =
				(jobjectArray) env->CallObjectMethod(obj, method, params, si);

		if(checkException(env, errMsg()))
		{
			toC(env, res, tmp);
		}
	}

};


////////////////////////////////////////////////////////////////////////////////
//	VRLUserFunction
////////////////////////////////////////////////////////////////////////////////

template <typename TData, int dim, typename TDataIn>
class VRLUserLinker
	: public StdDataLinker<VRLUserLinker<TData, dim, TDataIn>, TData, dim>
{
	public:
	///	type of base class
		typedef DataLinker<TData, dim> base_type;

	protected:
		static jobject compileUserDataString(JNIEnv *env, const char* s)
		{
			jclass cls = env->FindClass("edu/gcsc/vrl/ug/UserDataCompiler");
			if (env->ExceptionCheck()) {env->ExceptionDescribe();}

			jmethodID runMethod = env->GetStaticMethodID(
					cls, "compile",
					"(Ljava/lang/String;I)Ljava/lang/Object;");
			if (env->ExceptionCheck()) {env->ExceptionDescribe();}

			return env->CallStaticObjectMethod(cls, runMethod, stringC2J(env, s),
			                                   retArrayDim);
		}

		static jclass getUserDataClass(JNIEnv *env)
		{
			jclass result = env->FindClass("edu/gcsc/vrl/ug/UserData");
			if (env->ExceptionCheck()) {env->ExceptionDescribe();}

			return result;
		}

		static jmethodID getUserDataRunMethod(JNIEnv *env, jclass cls)
		{
			std::string signature = vrl_traits<TData>::callSignature();
			std::stringstream mName; mName << "run" << retArrayDim;
			jmethodID result = env->GetMethodID(cls, mName.str().c_str(), signature.c_str());

			if (!checkException(env))
			{
				UG_LOG("[VRL-Bindings] "<<name()<<" Error:"
						<< " cannot find userdata method."
						<< " Please check your implementation!" << std::endl);
			}
			return result;
		}

	public:
		static const unsigned int retArrayDim = vrl_traits<TData>::retArrayDim;

		VRLUserLinker() : initialized(false) {}

		static std::string name()
		{
			std::stringstream ss;
			ss << "VRLUserLinker"<<vrl_traits<TData>::name() << vrl_traits<TDataIn>::name() << dim << "d";
			return ss.str();
		}

		static std::string group_name()
		{
			std::stringstream ss;
			ss << "VRLUserLinker"<<vrl_traits<TData>::name()<<vrl_traits<TDataIn>::name();
			return ss.str();
		}

		size_t num_args() const {return m_numArgs;}
		const char* vrl_value_callback() {return m_ValueCode.c_str();}
		const char* vrl_deriv_callback(size_t arg) {return m_vDerivCode[arg].c_str();}

		void set_vrl_value_callback(const char* expression, size_t numArgs)
		{
			m_numArgs = numArgs;
			vUserDataDeriv.resize(numArgs);
			m_vDerivCode.resize(numArgs);
			set_num_input(numArgs);

			JNIEnv* env = threading::getEnv(getJavaVM());

			releaseGlobalRefs();

			m_ValueCode = expression;
			userDataValue = compileUserDataString(env, expression);
			userDataClass = getUserDataClass(env);
			runMethod = getUserDataRunMethod(env, userDataClass);

			checkException(env, name().append(": Cannot setup evaluation class or method."));

			// create thread-safe references
			// (GC won't deallocate them until manual deletion request)
			userDataValue = env->NewGlobalRef(userDataValue);
			userDataClass = (jclass) env->NewGlobalRef((jobject) userDataClass);
			checkException(env, name().append(": Global Reference Error."));
		}

	///	sets the vrl function used to compute the derivative
		void set_vrl_deriv_callback(size_t arg, const char* expression)
		{
			JNIEnv* env = threading::getEnv(getJavaVM());

			m_vDerivCode[arg] = expression;
			vUserDataDeriv[arg] = compileUserDataString(env, expression);
			checkException(env, name().append(": Cannot setup evaluation class or method."));

			// create thread-safe references
			// (GC won't deallocate them until manual deletion request)
			vUserDataDeriv[arg] = env->NewGlobalRef(vUserDataDeriv[arg]);
			checkException(env, name().append(": Global Reference Error."));
		}

	///	set input i
		void set_input(size_t i, SmartPtr<UserData<TDataIn, dim> > data)
		{
			UG_ASSERT(i < m_vpUserData.size(), "Input not needed");
			UG_ASSERT(i < m_vpDependData.size(), "Input not needed");

		//	check input number
			if(i >= this->num_input())
				UG_THROW("DataLinker::set_input: Only " << this->num_input()
								<< " inputs can be set. Use 'set_num_input' to increase"
								" the number of needed inputs.");

		//	remember userdata
			m_vpUserData[i] = data;

		//	cast to dependent data
			m_vpDependData[i] = data.template cast_dynamic<DependentUserData<TDataIn, dim> >();

		//	forward to base class
			base_type::set_input(i, data);
		}

	///	computes the value
		virtual void compute(LocalVector* u, GeometricObject* elem, bool bDeriv = false)
		{
		//	vector of data for all inputs
			std::vector<TDataIn> vDataIn(this->num_input());

			const number t = this->time();
			const int si = this->subset();

			for(size_t s = 0; s < this->num_series(); ++s)
				for(size_t ip = 0; ip < this->num_ip(s); ++ip)
				{
				//	gather all input data for this ip
					for(size_t c = 0; c < vDataIn.size(); ++c)
						vDataIn[c] = m_vpUserData[c]->value(this->series_id(c,s), ip);

				//	evaluate data at ip
					eval_value(this->value(s,ip), vDataIn, this->ip(s, ip), t, si);
				}

		//	check if derivative is required
			if(!bDeriv || this->zero_derivative()) return;

		//	clear all derivative values
			this->clear_derivative_values();

		//	loop all inputs
			for(size_t c = 0; c < vDataIn.size(); ++c)
			{
			//	check if input has derivative
				if(this->zero_derivative(c)) continue;

			//	loop ips
				for(size_t s = 0; s < this->num_series(); ++s)
					for(size_t ip = 0; ip < this->num_ip(s); ++ip)
					{
					//	gather all input data for this ip
						vDataIn[c] = m_vpUserData[c]->value(this->series_id(c,s), ip);

					//	data of derivative w.r.t. one component at ip-values
						TData derivVal;

					//	evaluate data at ip
						eval_deriv(derivVal, vDataIn, this->ip(s, ip), t, si, c);

					//	loop functions
						for(size_t fct = 0; fct < this->input_num_fct(c); ++fct)
						{
						//	get common fct id for this function
							const size_t commonFct = this->input_common_fct(c, fct);

						//	loop dofs
							for(size_t dof = 0; dof < this->num_sh(fct); ++dof)
							{
								linker_traits<TData, TDataIn>::
								mult_add(this->deriv(s, ip, commonFct, dof),
								         derivVal,
								         m_vpDependData[c]->deriv(this->series_id(c,s), ip, fct, dof));
							}
						}
					}
			}
		}


		inline void evaluate (TData& value,
							  const MathVector<dim>& globIP,
							  number time, int si) const
		{
		//	vector of data for all inputs
			std::vector<TDataIn> vDataIn(this->num_input());

		//	gather all input data for this ip
			for(size_t c = 0; c < vDataIn.size(); ++c)
				(*m_vpUserData[c])(vDataIn[c], globIP, time, si);

		//	evaluate data at ip
			eval_value(value, vDataIn, globIP, time, si);
		}


		template <int refDim>
		inline void evaluate (TData& value,
							  const MathVector<dim>& globIP,
							  number time, int si,
							  LocalVector& u,
							  GeometricObject* elem,
							  const MathVector<dim> vCornerCoords[],
							  const MathVector<refDim>& locIP) const
		{
		//	vector of data for all inputs
			std::vector<TDataIn> vDataIn(this->num_input());

		//	gather all input data for this ip
			for(size_t c = 0; c < vDataIn.size(); ++c)
				(*m_vpUserData[c])(vDataIn[c], globIP, time, si, u, elem, vCornerCoords, locIP);

		//	evaluate data at ip
			eval_value(value, vDataIn, globIP, time, si);
		}


		template <int refDim>
		inline void evaluate(TData vValue[],
							 const MathVector<dim> vGlobIP[],
							 number time, int si,
							 LocalVector& u,
							 GeometricObject* elem,
							 const MathVector<dim> vCornerCoords[],
							 const MathVector<refDim> vLocIP[],
							 const size_t nip,
							 const MathMatrix<refDim, dim>* vJT = NULL) const
		{
		//	vector of data for all inputs
			std::vector<TDataIn> vDataIn(this->num_input());

		//	gather all input data for this ip
			for(size_t ip = 0; ip < nip; ++ip)
			{
				for(size_t c = 0; c < vDataIn.size(); ++c)
					(*m_vpUserData[c])(vDataIn[c], vGlobIP[ip], time, si, u, elem, vCornerCoords, vLocIP[ip]);

			//	evaluate data at ip
				eval_value(vValue[ip], vDataIn, vGlobIP[ip], time, si);
			}
		}

	public:
	///	evaluates the data at a given point and time
		void eval_value(TData& value, const std::vector<TDataIn>& dataIn,
						const MathVector<dim>& x, number time, int si) const
		{
			JNIEnv* env = threading::getEnv(getJavaVM());

			const int dataSize = vrl_traits<TDataIn>::size * dataIn.size();
			int argSize = 	dataSize // data
							+ vrl_traits<MathVector<dim> >::size // x - vector
							+ vrl_traits<number>::size; // time

			jdoubleArray params = env->NewDoubleArray(argSize);

			// copy data
			jdouble* vTmp = new jdouble[dataSize];
			for (size_t i = 0; i < dataIn.size(); i++) {
				vrl_traits<TDataIn>::push(vTmp, dataIn[i]);
				env->SetDoubleArrayRegion(params, i*vrl_traits<TDataIn>::size,
				                          	  	 (i+1)*vrl_traits<TDataIn>::size, vTmp);
			}
			delete[] vTmp;
			checkException(env, "CopyParametersToJava: error");

			// copy vector and time
			jdouble vTmp2[dim + 1];
			for (size_t i = 0; i < dim; i++) vTmp2[i] = x[i];
			vTmp2[dim] = time;
			env->SetDoubleArrayRegion(params, dataSize, dataSize+dim+1, vTmp2);
			checkException(env, "CopyParametersToJava: error");

			jint jsi = si;

			if (runMethod != NULL)
				vrl_traits<TData>::call(env, value, userDataValue, runMethod,
				                        params, jsi);
		}

	///	evaluates the data at a given point and time
		void eval_deriv(TData& value, const std::vector<TDataIn>& dataIn,
						const MathVector<dim>& x, number time, int si, size_t arg) const
		{
			JNIEnv* env = threading::getEnv(getJavaVM());

			const int dataSize = vrl_traits<TDataIn>::size * dataIn.size();
			int argSize = 	dataSize // data
							+ vrl_traits<MathVector<dim> >::size // x - vector
							+ vrl_traits<number>::size; // time

			jdoubleArray params = env->NewDoubleArray(argSize);

			// copy data
			jdouble* vTmp = new jdouble[dataSize];
			for (size_t i = 0; i < dataIn.size(); i++) {
				vrl_traits<TDataIn>::push(vTmp, dataIn[i]);
				env->SetDoubleArrayRegion(params, i*vrl_traits<TDataIn>::size,
				                          	  	 (i+1)*vrl_traits<TDataIn>::size, vTmp);
			}
			delete[] vTmp;
			checkException(env, "CopyParametersToJava: error");

			// copy vector and time
			jdouble vTmp2[dim + 1];
			for (size_t i = 0; i < dim; i++) vTmp2[i] = x[i];
			vTmp2[dim] = time;
			env->SetDoubleArrayRegion(params, dataSize, dataSize+dim+1, vTmp2);
			checkException(env, "CopyParametersToJava: error");

			jint jsi = si;

			if (runMethod != NULL)
				vrl_traits<TData>::call(env, value, vUserDataDeriv[arg], runMethod,
				                        params, jsi);
		}

	public:
		void releaseGlobalRefs()
		{
			// deleting thread-safe global references
			if (initialized) {
				JNIEnv* localEnv = threading::getEnv(getJavaVM());
				localEnv->DeleteGlobalRef(userDataValue);
				localEnv->DeleteGlobalRef((jobject) userDataClass);
			}
		}

		~VRLUserLinker()
		{
			releaseGlobalRefs();
		}

	protected:
	///	set number of needed inputs
		void set_num_input(size_t num)
		{
		//	resize arrays
			m_vpUserData.resize(num, NULL);
			m_vpDependData.resize(num, NULL);

		//	forward size to base class
			base_type::set_num_input(num);
		}

	private:
		size_t m_numArgs;
		bool initialized;
		jobject userDataValue;
		std::vector<jobject> vUserDataDeriv;
		jclass userDataClass;
		jmethodID runMethod;

		std::string m_ValueCode;
		std::vector<std::string> m_vDerivCode;

	protected:
	///	data input
		std::vector<SmartPtr<UserData<TDataIn, dim> > > m_vpUserData;

	///	data input casted to dependend data
		std::vector<SmartPtr<DependentUserData<TDataIn, dim> > > m_vpDependData;
};


////////////////////////////////////////////////////////////////////////////////
// VRLUserData
////////////////////////////////////////////////////////////////////////////////

template <typename TData, int dim>
class VRLUserData
	: public StdPositionData<VRLUserData<TData, dim>, TData, dim>
{
	protected:
		static jobject compileUserDataString(JNIEnv *env, const char* s)
		{
			jclass cls = env->FindClass("edu/gcsc/vrl/ug/UserDataCompiler");
			if (env->ExceptionCheck()) {env->ExceptionDescribe();}

			jmethodID runMethod = env->GetStaticMethodID(
					cls, "compile",
					"(Ljava/lang/String;I)Ljava/lang/Object;");
			if (env->ExceptionCheck()) {env->ExceptionDescribe();}

			return env->CallStaticObjectMethod(cls, runMethod, stringC2J(env, s),
			                                   retArrayDim);
		}

		static jclass getUserDataClass(JNIEnv *env)
		{
			jclass result = env->FindClass("edu/gcsc/vrl/ug/UserData");
			if (env->ExceptionCheck()) {env->ExceptionDescribe();}

			return result;
		}

		static jmethodID getUserDataRunMethod(JNIEnv *env, jclass cls)
		{
			std::string signature = vrl_traits<TData>::callSignature();
			std::stringstream mName; mName << "run" << retArrayDim;
			jmethodID result = env->GetMethodID(cls, mName.str().c_str(), signature.c_str());

			if (!checkException(env))
			{
				UG_LOG("[VRL-Bindings] "<<name()<<" Error:"
						<< " cannot find userdata method."
						<< " Please check your implementation!" << std::endl);
			}
			return result;
		}

	public:
		static const unsigned int retArrayDim = vrl_traits<TData>::retArrayDim;

		VRLUserData() : initialized(false) {}

		static std::string params()
		{
			// DO NOT USE underscore for param names (i.e. NO "_x" !!!)
			std::stringstream params;
			if(dim >= 1) params <<  "\"" << "x" << "\"";
			if(dim >= 2) params << ",\"" << "y" << "\"";
			if(dim >= 3) params << ",\"" << "z" << "\"";
						 params << ",\"" << "t" << "\"";
						 params << ",\"" << "si" << "\"";
			return params.str();
		}

		static std::string name()
		{
			std::stringstream ss;
			ss << "VRLUser"<<vrl_traits<TData>::name() << dim << "d";
			return ss.str();
		}

		static std::string group_name()
		{
			std::stringstream ss;
			ss << "VRLUser"<<vrl_traits<TData>::name();
			return ss.str();
		}

		void set_vrl_callback(const char* expression)
		{
			JNIEnv* env = threading::getEnv(getJavaVM());

			releaseGlobalRefs();

			userDataObject = compileUserDataString(env, expression);
			userDataClass = getUserDataClass(env);
			runMethod = getUserDataRunMethod(env, userDataClass);

			checkException(env, name().append(": Cannot setup evaluation class or method."));

			// create thread-safe references
			// (GC won't deallocate them until manual deletion request)
			userDataObject = env->NewGlobalRef(userDataObject);
			userDataClass = (jclass) env->NewGlobalRef((jobject) userDataClass);
			checkException(env, name().append(": Global Reference Error."));
			initialized = true;
		}

		///	evaluates the data at a given point and time
		inline void evaluate(TData& value, const MathVector<dim>& x, number time, int si) const
		{
			JNIEnv* env = threading::getEnv(getJavaVM());

			// convert parameters
			jdoubleArray params = ConvertParametersToJava(env, x, time);
			jint jsi = si;

			if (runMethod != NULL)
				vrl_traits<TData>::call(env, value, userDataObject, runMethod,
				                        params, jsi);
		}

		void releaseGlobalRefs()
		{
			// deleting thread-safe global references
			if (initialized) {
				JNIEnv* localEnv = threading::getEnv(getJavaVM());
				localEnv->DeleteGlobalRef(userDataObject);
				localEnv->DeleteGlobalRef((jobject) userDataClass);
			}
		}

		~VRLUserData()
		{
			releaseGlobalRefs();
		}

	private:
		bool initialized;
		jobject userDataObject;
		jclass userDataClass;
		jmethodID runMethod;
};

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

template <int dim>
class VRLCondUserNumber : public UserData<number, dim, bool>
{
protected:
	jdouble condData2Double(JNIEnv *env, jobject obj) const
	{
		// todo: cache this
		jclass cls = env->FindClass("edu/gcsc/vrl/ug/Cond");
		if (env->ExceptionCheck()) {env->ExceptionDescribe();}

		jmethodID method = env->GetMethodID(cls, "getValue", "()D");
		if (env->ExceptionCheck()) {env->ExceptionDescribe();}

		return env->CallDoubleMethod(obj, method);
	}

	jdouble condData2Boolean(JNIEnv *env, jobject obj) const
	{
		// todo: cache this
		jclass cls = env->FindClass("edu/gcsc/vrl/ug/Cond");
		if (env->ExceptionCheck()) {env->ExceptionDescribe();}

		jmethodID method = env->GetMethodID(cls, "getCondBool", "()Z");
		if (env->ExceptionCheck()) {env->ExceptionDescribe();}

		return env->CallBooleanMethod(obj, method);
	}

	jobject compileCondUserDataString(JNIEnv *env, const char* s) const
	{
		jclass cls = env->FindClass("edu/gcsc/vrl/ug/CondUserDataCompiler");
		if (env->ExceptionCheck()) {env->ExceptionDescribe();}

		jmethodID runMethod = env->GetStaticMethodID(
				cls, "compile",
				"(Ljava/lang/String;)Ljava/lang/Object;");
		if (env->ExceptionCheck()) {env->ExceptionDescribe();}

		return env->CallStaticObjectMethod(cls, runMethod, stringC2J(env, s));
	}

	jclass getCondUserDataClass(JNIEnv *env) const
	{
		jclass result = env->FindClass("edu/gcsc/vrl/ug/CondUserData");
		if (env->ExceptionCheck()) {env->ExceptionDescribe();}

		return result;
	}

	jmethodID getCondUserDataRunMethod(JNIEnv *env, jclass cls) const
	{
		std::string signature = "([DI)Ledu/gcsc/vrl/ug/Cond;";
		std::string mName = "run";

		jmethodID result = env->GetMethodID(cls, mName.c_str(), signature.c_str());
		if (!checkException(env))
		{
			UG_LOG("[VRL-Bindings] Error:"
					<< " cannot find cond-userdata method."
					<< " Please check your implementation!" << std::endl);
		}

		return result;
	}

public:
	static std::string params()
	{
		// DO NOT USE underscore for param names (i.e. NO "_x" !!!)
		std::stringstream params;
		if(dim >= 1) params <<  "\"" << "x" << "\"";
		if(dim >= 2) params << ",\"" << "y" << "\"";
		if(dim >= 3) params << ",\"" << "z" << "\"";
					 params << ",\"" << "t" << "\"";
					 params << ",\"" << "si" << "\"";
		return params.str();
	}

	static std::string name()
	{
		std::stringstream ss;
		ss << "VRLCondUser"<<vrl_traits<number>::name() << dim << "d";
		return ss.str();
	}

	static std::string group_name()
	{
		std::stringstream ss;
		ss << "VRLCondUser"<<vrl_traits<number>::name();
		return ss.str();
	}

	VRLCondUserNumber() : initialized(false)
	{
		std::stringstream stream;
		stream << "<font color=\"red\">VRLCondUserNumber"
				<< dim << "D: invokation error:</font>";
		invocationErrorMsg = stream.str();
	}

	void set_vrl_callback(const char* expression)
	{
		JNIEnv* env = threading::getEnv(getJavaVM());

		releaseGlobalRefs();

		userDataObject = compileCondUserDataString(env, expression);
		userDataClass = getCondUserDataClass(env);
		runMethod = getCondUserDataRunMethod(env, userDataClass);

		// create thread-safe references 
		// (GC won't deallocate them until manual deletion request)
		userDataObject = env->NewGlobalRef(userDataObject);
		userDataClass = (jclass) env->NewGlobalRef((jobject) userDataClass);

		initialized = true;
	}

	///	evaluates the data at a given point and time
	bool operator() (number& c, const MathVector<dim>& x, number time, int si) const
	{
		JNIEnv* env = threading::getEnv(getJavaVM());

		// convert parameters
		jdoubleArray params = ConvertParametersToJava(env, x, time);
		jint jsi = si;

		bool result = false;
		if (runMethod != NULL)
		{
			jobject bndResult = (jdoubleArray) env->CallObjectMethod(
					userDataObject,
					runMethod,
					params, jsi);

			if (checkException(env, invocationErrorMsg)) {
				result = condData2Boolean(env, bndResult);
				c = condData2Double(env, bndResult);
			}
		}

		return result;
	}

	virtual void compute(LocalVector* u, GeometricObject* elem, bool computeDeriv = false)
	{
		// \todo: should remember flag
		for (size_t s = 0; s < this->num_series(); ++s)
		{
			for (size_t i = 0; i < this->num_ip(s); ++i)
			{
				this->operator()(this->value(s, i), this->ip(s, i),
				                 this->time(), this->subset());
			}
		}
	}

	void releaseGlobalRefs()
	{
		// deleting thread-safe global references
		if (initialized) {
			JNIEnv* localEnv = threading::getEnv(getJavaVM());
			localEnv->DeleteGlobalRef(userDataObject);
			localEnv->DeleteGlobalRef((jobject) userDataClass);
		}
	}

	~VRLCondUserNumber() {
		releaseGlobalRefs();
	}

private:
	std::string invocationErrorMsg;
	bool initialized;
	jobject userDataObject;
	jclass userDataClass;
	jmethodID runMethod;
};

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

template <typename TData>
class PrintUserData2d {
public:
	void set(SmartPtr<UserData<TData, 2> > user) {m_spNumber = user;}

	static std::string dataname()
	{
		std::stringstream ss; ss << "User"<<vrl_traits<TData>::name()<<"2d";
		return ss.str();
	}

	static std::string name()
	{
		std::stringstream ss; ss << "PrintUser"<<vrl_traits<TData>::name()<<"2d";
		return ss.str();
	}

	std::string print(number x, number y, number time, int si)
	{
		MathVector < 2 > v(x, y);
		TData ret;

		if (m_spNumber.valid()) (*m_spNumber)(ret, v, time, si);
		else UG_THROW(name()<<": Data not set.");

		std::stringstream ss;
		ss << ret << std::endl;
		return ss.str();
	}

private:
	SmartPtr<UserData<TData, 2> > m_spNumber;
};

template <typename TData>
class PrintCondUserData2d {
public:
	static std::string dataname()
	{
		std::stringstream ss; ss << "CondUser"<<vrl_traits<TData>::name()<<"2d";
		return ss.str();
	}

	static std::string name()
	{
		std::stringstream ss; ss << "PrintCondUser"<<vrl_traits<TData>::name()<<"2d";
		return ss.str();
	}

	void set(SmartPtr<UserData<number, 2, bool> > user) {m_spData = user;}

	std::string print(number x, number y, number time, int si)
	{
		MathVector < 2 > v(x, y);
		number ret;
		bool bndResult = false;

		if (m_spData.valid()) bndResult = (*m_spData)(ret, v, time, si);
		else {
			UG_THROW(name()<<": Data not set.");
			ret = -1;
		}

		std::stringstream stream;
		stream << "[";
		if(bndResult) stream << "true";
		else stream << "false";
		stream << ", " << ret << "]";
		return stream.str();
	}

private:
	SmartPtr<UserData<number, 2, bool> > m_spData;
};

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

template <typename TData, int dim>
void RegisterUserDataType(ug::bridge::Registry& reg, const std::string& grp)
{
	std::string tag = ug::bridge::GetDimensionTag<dim>();

	//	VRLUserType
	{
		typedef VRLUserData<TData, dim> T;
		typedef UserData<TData, dim> TBase;
		std::stringstream options;
		options << "Input:|user-data|dim=" << T::retArrayDim << ";"
				<< "params=["<<T::params()<<"];";
		reg.add_class_<T, TBase>(T::name(), grp)
			.add_constructor()
			.add_method("data", &T::set_vrl_callback, "", options.str().c_str())
			.set_construct_as_smart_pointer(true);
		reg.add_class_to_group(T::name(), T::group_name(), tag);
	}

	//	VRLUserLinkerTypeNumber
	{
		typedef VRLUserLinker<TData, dim, number> T;
		typedef UserData<TData, dim> TBase;
		std::stringstream options;
		reg.add_class_<T, TBase>(T::name(), grp)
			.add_constructor()
			.add_method("num_args", &T::num_args)
			.add_method("data", &T::set_vrl_value_callback)
			.add_method("vrl_value_callback", &T::vrl_value_callback)
			.add_method("deriv", &T::set_vrl_deriv_callback)
			.add_method("vrl_deriv_callback", &T::vrl_deriv_callback)
			.add_method("set_input", static_cast<void (T::*)(size_t, SmartPtr<UserData<number, dim> >)>(&T::set_input))
			.set_construct_as_smart_pointer(true);
		reg.add_class_to_group(T::name(), T::group_name(), tag);
	}

	// PrintUserType2d
	if(dim == 2)
	{
		typedef PrintUserData2d<TData> T;
		reg.add_class_<T > (T::name(), grp)
				.add_constructor()
				.add_method("set", &T::set, "", T::dataname())
				.add_method("print", &T::print, "Result", "x#y#t#si");
	}
}

template <int dim>
void RegisterUserData(ug::bridge::Registry& reg, const char* parentGroup)
{
	// 	get group
	std::string grp = std::string(parentGroup);
	std::string tag = ug::bridge::GetDimensionTag<dim>();

	RegisterUserDataType<number, dim>(reg, grp);
	RegisterUserDataType<MathVector<dim>, dim>(reg, grp);
	RegisterUserDataType<MathMatrix<dim,dim>, dim>(reg, grp);

	//	VRLCondUserNumber
	{
		typedef VRLCondUserNumber<dim> T;
		typedef UserData<number, dim, bool> TBase;
		std::stringstream options;
		options << "Input:|cond-user-data|params=["<<T::params()<<"];";
		reg.add_class_<T, TBase>(T::name(), grp)
			.add_constructor()
			.add_method("data", &T::set_vrl_callback, "", options.str().c_str())
			.set_construct_as_smart_pointer(true);
		reg.add_class_to_group(T::name(), T::group_name(), tag);
	}

	if(dim == 2)
	{
		typedef PrintCondUserData2d<number> T3;
		reg.add_class_<T3 > (T3::name(), grp)
				.add_constructor()
				.add_method("set", &T3::set, "", T3::dataname())
				.add_method("print", &T3::print, "Result", "x#y#t#si");
	}
}

void RegisterUserData(ug::bridge::Registry& reg, const char* parentGroup)
{
#ifdef UG_DIM_1
	ug::vrl::RegisterUserData < 1 > (reg, parentGroup);
#endif
#ifdef UG_DIM_2
	ug::vrl::RegisterUserData < 2 > (reg, parentGroup);
#endif
#ifdef UG_DIM_3
	ug::vrl::RegisterUserData < 3 > (reg, parentGroup);
#endif
}

} // vrl::
} // ug::
