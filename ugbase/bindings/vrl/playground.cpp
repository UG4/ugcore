/*
 * Copyright (c) 2011-2012:  Steinbeis Forschungszentrum (STZ Ölbronn)
 * Author: Michael Hoffer
 * 
 * This file is part of UG4.
 * 
 * UG4 is free software: you can redistribute it and/or modify it under the
 * terms of the GNU Lesser General Public License version 3 (as published by the
 * Free Software Foundation) with the following additional attribution
 * requirements (according to LGPL/GPL v3 §7):
 * 
 * (1) The following notice must be displayed in the Appropriate Legal Notices
 * of covered and combined works: "Based on UG4 (www.ug4.org/license)".
 * 
 * (2) The following notice must be displayed at a prominent place in the
 * terminal output of covered works: "Based on UG4 (www.ug4.org/license)".
 * 
 * (3) The following bibliography is recommended for citation and must be
 * preserved in all covered files:
 * "Reiter, S., Vogel, A., Heppner, I., Rupp, M., and Wittum, G. A massively
 *   parallel geometric multigrid solver on hierarchically distributed grids.
 *   Computing and visualization in science 16, 4 (2013), 151-164"
 * "Vogel, A., Reiter, S., Rupp, M., Nägel, A., and Wittum, G. UG4 -- a novel
 *   flexible software system for simulating pde based models on high performance
 *   computers. Computing and visualization in science 16, 4 (2013), 165-179"
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU Lesser General Public License for more details.
 */

#include "playground.h"
#include "bindings_vrl.h"



namespace ug {
namespace vrl {


TestClass::TestClass() {
	UG_LOG("Constructor TestClass() called.\n");
}

TestClass::TestClass(std::string name) {
	UG_LOG("Constructor TestClass(std::string name) called.\n");
}

int TestClass::add(int a, int b) {
	return a + b;
}

std::string TestClass::getString() {
	UG_LOG("Test123\n");
	return "Test123";
}

SmartPtr<TestClass> TestClass::smartTestImpl() {
	return SmartPtr<TestClass > (new TestClass());
}

ConstSmartPtr<TestClass> TestClass::constSmartTestImpl() {
	return ConstSmartPtr<TestClass > (new TestClass());
}

int TestClass::print_name() {
	UG_LOG("Name is Test\n");
	return 1;
}

int TestClass::print() {
	UG_LOG("Test::print()\n");
	return 0;
}

int TestClass::print2() {
	UG_LOG("Test::print2()\n");
	return 1;
}

int TestClass::print2() const {
	UG_LOG("Test::print2() const\n");
	return 1;
}

TestClass::~TestClass() {
	UG_LOG("~TestClass" << std::endl);
}

int SmartTestFunction(SmartPtr<TestClass> test) {
	UG_LOG("SmartTestFunc: ");

	test->print2();

	return test.refcount();
}

int ConstSmartTestFunction(ConstSmartPtr<TestClass> test) {
	UG_LOG("ConstSmartTestFunc: ");

	test->print2();

	return test.refcount();
}

void TestSmartPtr2ConstPtr(const TestClass* t) {
	t->print2();
}

SmartPtrCls::SmartPtrCls() {
	UG_LOG("Constructor SmartPtrCls() called." << std::endl);
	_name = "noname";
	_data = nullptr;
}

SmartPtrCls::SmartPtrCls(std::string name) {
	_name = name;
	UG_LOG("Constructor SmartPtrCls(std::string name) called." << std::endl);
	_data = nullptr;
}

void SmartPtrCls::print_name() {
	UG_LOG("Name is " << _name << std::endl);
}

void SmartPtrCls::create_data(int size) {
	if (_data ==nullptr) {
		_data = new byte[size];
	} else {
		UG_LOG(_name << ": Data already created! " << std::endl);
	}
}

SmartPtrCls::~SmartPtrCls() {
	UG_LOG("~SmartPtrCls:" << _name << std::endl);

	if (_data !=nullptr ){
		delete[] _data;
		_data = nullptr;
	}
}





//cpoliwoda start playground

	//for debug / understand how to call a c++ function from java
	//the 1st method
	void ChrisPoliTest(){
		//std::cout   << "trunk/ugbase/bindings/vrl/bindings_vrl.cpp : ChrisPoliTest() "
		std::cout   << "trunk/ugbase/bridge/misc_bridges/test_bridge.cpp : ChrisPoliTest() "
					<< std::endl
					<< "THANK YOU  :-)"
					<< std::endl;
	}

	//the 2nd method
	void ChrisPoliTest(int i){
		//std::cout   << "trunk/ugbase/bindings/vrl/bindings_vrl.cpp : ChrisPoliTest() "
		std::cout   << "trunk/ugbase/bridge/misc_bridges/test_bridge.cpp : ChrisPoliTest() "
					<< std::endl
					<< "THANK YOU "<< i<<"-times  :-)"
					<< std::endl;
	}



	//the 3th method
	void ChrisPoliTest(bool array[]){
		// need to be the same size as on java side the size of the array
		// set manually because do not know how to get size automatically
		int arraySize = 2;

		std::cout << "cpp true = " << true << std::endl;
		std::cout << "cpp false = " << false << std::endl;

		for (int i = 0; i < arraySize; ++i) {
			std::cout   << "array[ "<< i <<" ] = "<< array[i] << std::endl;
		}

	}

	//the 4th method
	void ChrisPoliTest(std::vector<bool> vec){

		std::cout << "cpp : trunk/ugbase/bridge/misc_bridges/test_bridge.cpp ."<<
				" ChrisPoliTest(std::vector<bool> vec)"  << std::endl;


		std::cout << "cpp true = " << true << std::endl;
		std::cout << "cpp false = " << false << std::endl;

		for (size_t i = 0; i < vec.size(); ++i) {
			std::cout   << "vec[ "<< i <<" ] = "<< vec[i] << std::endl;
		}

	}

	//the 5th method
	// methods with same name and parameters can NOT be differed by their
	// return type therefore we have to change the functionname too
	// after changing the return type !!!
	std::vector<bool> ChrisPoliTestReturn(std::vector<bool> vec){

		std::cout << "cpp : trunk/ugbase/bridge/misc_bridges/test_bridge.cpp ."<<
				" std::vector<bool> ChrisPoliTestReturn(std::vector<bool> vec)"
				<< std::endl;

		/*std::cout << std::boolalpha()<< std::endl;
		std::cout << "cpp true = " << true << std::endl;
		std::cout << "cpp false = " << false << std::endl;
		std::cout << std::noboolalpha()<< std::endl;
		*/
		std::cout << "cpp true = " << true << std::endl;
		std::cout << "cpp false = " << false << std::endl;

		for (size_t i = 0; i < vec.size(); ++i) {
			std::cout   << "vec[ "<< i <<" ] = "<< vec[i] << std::endl;
		}

		std::cout << "cpp : ChrisPoliTestReturn( ) BEFORE return" << std::endl;

		return vec;
	}


	// this methode is for testing data transfer / communication via file
	// in vrl-ug remote modus
	// the void method with const paths to input & output file
	void ChrisPoliTestCreateFile(){

		std::cout << "cpp : trunk/ugbase/bridge/misc_bridges/test_bridge.cpp ."<<
				" void ChrisPoliTestCreateFile( )"
				<< std::endl;

		std::ifstream infile;
		std::ofstream outfile;

		std::string inputPath = "/Users/christianpoliwoda/Desktop/output/input.txt";
		std::string outputPath = "/Users/christianpoliwoda/Desktop/output/output.txt";

		std::string tmpReadLine;

		    infile.open(inputPath.c_str());
		    outfile.open(outputPath.c_str());

		    if (infile.is_open())
		        while (!infile.eof()) {// To get you all the lines.
		            {
		                getline(infile, tmpReadLine); // Saves the line in STRING.
		                std::cout << tmpReadLine << std::endl; // Prints our STRING. and begin new line for next string
		                outfile << "COPY " << tmpReadLine << std::endl;
		            }
		        }

		    infile.close();
		    outfile.close();
	}

	// this methode is for testing data transfer / communication via file
	// in vrl-ug remote modus
	std::string ChrisPoliTestCreateFile(std::string path){

		std::cout << "cpp : trunk/ugbase/bridge/misc_bridges/test_bridge.cpp ."<<
				" std::string ChrisPoliTestCreateFile(std::string path)"
				<< std::endl;



		std::ifstream infile;
		std::ofstream outfile;

		std::string outputPath = "/Users/christianpoliwoda/Desktop/output/example-output.txt";

		std::string tmpReadLine;

		    infile.open(path.c_str());
		    outfile.open(outputPath.c_str());

		    if (infile.is_open())
		        while (!infile.eof()) {// To get you all the lines.
		            {
		                getline(infile, tmpReadLine); // Saves the line in STRING.
		                std::cout << tmpReadLine << std::endl; // Prints our STRING. and begin new line for next string
		                outfile << "COPY " << tmpReadLine << std::endl;
		            }
		        }

		    infile.close();
		    outfile.close();

		return outputPath;
	}


	// this methode is for testing data transfer / communication via file
	// in vrl-ug remote modus
	std::string ChrisPoliTestCreateFile(std::string inPath, std::string outPath){

		std::cout << "cpp : trunk/ugbase/bindings/vrl/playground.cpp ."<<
				" std::string ChrisPoliTestCreateFile(std::string path)"
				<< std::endl;

		std::ifstream infile;
		std::ofstream outfile;

		std::cout << "cpp: inPath "<< inPath << std::endl;
		std::cout << "cpp: outPath "<< outPath << std::endl;

		std::string tmpReadLine;

		    infile.open(inPath.c_str());
		    outfile.open(outPath.c_str());

		    if (infile.is_open())
		        while (!infile.eof()) {// To get you all the lines.
		            {
		                getline(infile, tmpReadLine); // Saves the line in STRING.
		                std::cout << tmpReadLine << std::endl; // Prints our STRING. and begin new line for next string
		                outfile << "COPY " << tmpReadLine << std::endl;
		            }
		        }

		    infile.close();
		    outfile.close();

		return outPath;
	}


	//cpoliwoda end playground




void registerPlayground(ug::bridge::Registry& reg) {
	reg.add_class_<TestClass > ("TestClass", "ug4/testing")
			.add_constructor()
			.add_constructor<void(*)(std::string)>()
			.add_method("add", &TestClass::add, "result",
			"a|default|min=-3;max=5;value=-12#b|default|min=-1;max=1;value=23")
			.add_method("getString", &TestClass::getString)
			//			.add_method("performTest", &TestClass::performTest)
			.add_method("print", &TestClass::print)
			.add_method("smartTestImpl", &TestClass::smartTestImpl)
			.add_method("constSmartTestImpl", &TestClass::constSmartTestImpl);

	reg.add_function("SmartTestFunction", &SmartTestFunction, "ug4/testing");
	reg.add_function("ConstSmartTestFunction", &ConstSmartTestFunction, "ug4/testing");
	reg.add_function("TestSmartPtr2ConstPtr", &TestSmartPtr2ConstPtr, "ug4/testing");

	reg.add_class_<SmartPtrCls > ("SmartPtrCls", "ug4/testing")
			.add_constructor()
			.add_constructor<void(*)(std::string)>()
			.add_method("print_name", &SmartPtrCls::print_name)
			.add_method("create_data", &SmartPtrCls::create_data)
			.set_construct_as_smart_pointer(true);





	        //
			//  register cpoliwoda start playground
			//

			//reg.add_function("ChrisPoliTest",
			//		&ChrisPoliTest,
			//		grp);//before 2nd method with same name but with other params

			//reg.add_function("ChrisPoliTest",
			//		static_cast<void (*)(void)>(&ChrisPoliTest));
			//		//,grp);//the 1st method changed to this after 2nd method added

			/*
			 reg.add_function("ChrisPoliTest",
					static_cast<void (*)(int)>(&ChrisPoliTest),
					grp,"i");//the 2nd method

			reg.add_function("ChrisPoliTest",
							static_cast<void (*)(bool[])>(&ChrisPoliTest),
							grp,"bool-array");//the 3th method
			*/

			/*
			reg.add_function("ChrisPoliTest",
							static_cast<void (*)(std::vector<bool>)>(&ChrisPoliTest),
							grp,"bool-vec");//the 4th method
			*/


			reg.add_function("ChrisPoliTestReturn",
							( std::vector<bool> (*)(std::vector<bool>) )(
							&ChrisPoliTestReturn),"bool-vec-return");//the 5th method
							//,grp,"bool-vec-return");//the 5th method




			/*reg.add_function("ChrisPoliTestCreateFile",
							static_cast<void (*)(void)>
							(&ChrisPoliTest),grp);//the void method with const paths to input & output file
			*/
			//reg.add_function("ChrisPoliTestCreateFile",
			//				(void (*)(void))
			//				(&ChrisPoliTestCreateFile));//,grp);//the void method with const paths to input & output file



			/*reg.add_function("ChrisPoliTestCreateFile",
							static_cast<std::string (*)(std::string)>
							(&ChrisPoliTestCreateFile),grp,"path");//an other method
			*/
			//reg.add_function("ChrisPoliTestCreateFile",
			//				(std::string (*)(std::string))
			//				(&ChrisPoliTestCreateFile),"path");//an other method
			//				//,grp,"path");//an other method


			reg.add_function("ChrisPoliTestCreateFile",
										(std::string (*)(std::string, std::string))
										(&ChrisPoliTestCreateFile),"outpath","inPath#outPath");//an other method
										//,grp,"path");//an other method


			//
	        //  register cpoliwoda end playground
			//




}

} // end vrl::
}// end ug::


