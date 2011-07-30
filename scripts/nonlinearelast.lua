------------------------------------------------------------------------------
--
--   Lua - Script to perform the Elasticity-Problem for finite deformations
--
--   Author: Raphael Prohl
--
------------------------------------------------------------------------------

ug_load_script("ug_util.lua")

-- always 3d
dim = 3

-- choose dimension and algebra
InitUG(dim, CPUAlgebraSelector());
	
-----------------------------------------------------------------
--  Setup FE Nonlinear Element Discretization
-----------------------------------------------------------------

--creates a new instance of FE1 nonlinear Elasticity disc
elemDisc = FE1NonlinearElasticity3d()
--creates a new instance of Elasticity Tensor User Data
elast = ElasticityTensorUserData3d()
--method-call
elemDisc:set_elasticity_tensor(elast)

--initialize MathVector and MathTensor
--x = MathVector() --muss ich diese Klasse erst noch registrieren?
--Tensortest = MathTensor()
--Tensortest:set(0.0)

--method-call
elast:test_evaluate(0.0,0.0,0.0,0.0)

--MathTensor:operator<<(out,Tensortest)