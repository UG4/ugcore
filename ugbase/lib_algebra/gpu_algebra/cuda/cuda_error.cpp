#include <cuda_runtime.h>
namespace ug{
const char* CUDAError(cudaError_t err)
{
	switch (err)
	{
		case cudaSuccess:
			return "CUDA ERROR"
			"\nThe API call returned with no errors. In the case of query calls, this"
			"\ncan also mean that the operation being queried is complete (see"
			"\n::cudaEventQuery() and ::cudaStreamQuery()).";

		case cudaErrorMissingConfiguration: return
			"CUDA ERROR"
			"\nThe device function being invoked (usually via ::cudaLaunch()) was not"
			"\npreviously configured via the ::cudaConfigureCall() function.";

		case cudaErrorMemoryAllocation: return

			"CUDA ERROR"
			"\nThe API call failed because it was unable to allocate enough memory to"
			"\nperform the requested operation.";

		case cudaErrorInitializationError: return
			"CUDA ERROR"
			"\nThe API call failed because the CUDA driver and runtime could not be"
			"\ninitialized.";

		case cudaErrorLaunchFailure: return
			"CUDA ERROR"
			"\nAn exception occurred on the device while executing a kernel. Common"
			"\ncauses include dereferencing an invalid device pointer and accessing"
			"\nout of bounds shared memory. The device cannot be used until"
			"\n::cudaThreadExit() is called. All existing device memory allocations"
			"\nare invalid and must be reconstructed if the program is to continue"
			"\nusing CUDA.";

		case cudaErrorPriorLaunchFailure: return

			"CUDA ERROR"
			"\nThis indicated that a previous kernel launch failed. This was previously"
			"\nused for device emulation of kernel launches."
			"\ndeprecated"
			"\nThis error return is deprecated as of CUDA 3.1. Device emulation mode was"
			"\nremoved with the CUDA 3.1 release.";

		case cudaErrorLaunchTimeout: return


			"CUDA ERROR"
			"\nThis indicates that the device kernel took too long to execute. This can"
			"\nonly occur if timeouts are enabled - see the device property"
			"\n\ref ::cudaDeviceProp::kernelExecTimeoutEnabled \"kernelExecTimeoutEnabled\""
			"\nfor more information. The device cannot be used until ::cudaThreadExit()"
			"\nis called. All existing device memory allocations are invalid and must be"
			"\nreconstructed if the program is to continue using CUDA.";

		case cudaErrorLaunchOutOfResources: return

			"CUDA ERROR"
			"\nThis indicates that a launch did not occur because it did not have"
			"\nappropriate resources. Although this error is similar to"
			"\n::cudaErrorInvalidConfiguration, this error usually indicates that the"
			"\nuser has attempted to pass too many arguments to the device kernel, or the"
			"\nkernel launch specifies too many threads for the kernel's register count.";

		case cudaErrorInvalidDeviceFunction: return


			"CUDA ERROR"
			"\nThe requested device function does not exist or is not compiled for the"
			"\nproper device architecture.";

		case cudaErrorInvalidConfiguration: return

			"CUDA ERROR"
			"\nThis indicates that a kernel launch is requesting resources that can"
			"\nnever be satisfied by the current device. Requesting more shared memory"
			"\nper block than the device supports will trigger this error, as will"
			"\nrequesting too many threads or blocks. See ::cudaDeviceProp for more"
			"\ndevice limitations.";

		case cudaErrorInvalidDevice: return


			"CUDA ERROR"
			"\nThis indicates that the device ordinal supplied by the user does not"
			"\ncorrespond to a valid CUDA device.";

		case cudaErrorInvalidValue: return

			"CUDA ERROR"
			"\nThis indicates that one or more of the parameters passed to the API call"
			"\nis not within an acceptable range of values.";

		case cudaErrorInvalidPitchValue: return

			"CUDA ERROR"
			"\nThis indicates that one or more of the pitch-related parameters passed"
			"\nto the API call is not within the acceptable range for pitch.";

		case cudaErrorInvalidSymbol: return

			"CUDA ERROR"
			"\nThis indicates that the symbol name/identifier passed to the API call"
			"\nis not a valid name or identifier.";

		case cudaErrorMapBufferObjectFailed: return

			"CUDA ERROR"
			"\nThis indicates that the buffer object could not be mapped.";

		case cudaErrorUnmapBufferObjectFailed: return

			"CUDA ERROR"
			"\nThis indicates that the buffer object could not be unmapped.";

		case cudaErrorInvalidHostPointer: return

			"CUDA ERROR"
			"\nThis indicates that at least one host pointer passed to the API call is"
			"\nnot a valid host pointer.";

		case cudaErrorInvalidDevicePointer: return


			"CUDA ERROR"
			"\nThis indicates that at least one device pointer passed to the API call is"
			"\nnot a valid device pointer.";

		case cudaErrorInvalidTexture: return


			"CUDA ERROR"
			"\nThis indicates that the texture passed to the API call is not a valid"
			"\ntexture.";

		case cudaErrorInvalidTextureBinding: return

			"CUDA ERROR"
			"\nThis indicates that the texture binding is not valid. This occurs if you"
			"\ncall ::cudaGetTextureAlignmentOffset() with an unbound texture.";

		case cudaErrorInvalidChannelDescriptor: return

			"CUDA ERROR"
			"\nThis indicates that the channel descriptor passed to the API call is not"
			"\nvalid. This occurs if the format is not one of the formats specified by"
			"\n::cudaChannelFormatKind, or if one of the dimensions is invalid.";

		case cudaErrorInvalidMemcpyDirection: return

			"CUDA ERROR"
			"\nThis indicates that the direction of the memcpy passed to the API call is"
			"\nnot one of the types specified by ::cudaMemcpyKind.";

		case cudaErrorAddressOfConstant: return

			"CUDA ERROR"
			"\nThis indicated that the user has taken the address of a constant variable,"
			"\nwhich was forbidden up until the CUDA 3.1 release."
			"\ndeprecated"
			"\nThis error return is deprecated as of CUDA 3.1. Variables in constant"
			"\nmemory may now have their address taken by the runtime via"
			"\n::cudaGetSymbolAddress().";

		case cudaErrorTextureFetchFailed: return

			"CUDA ERROR"
			"\nThis indicated that a texture fetch was not able to be performed."
			"\nThis was previously used for device emulation of texture operations."
			"\ndeprecated"
			"\nThis error return is deprecated as of CUDA 3.1. Device emulation mode was"
			"\nremoved with the CUDA 3.1 release.";

		case cudaErrorTextureNotBound: return

			"CUDA ERROR"
			"\nThis indicated that a texture was not bound for access."
			"\nThis was previously used for device emulation of texture operations."
			"\ndeprecated"
			"\nThis error return is deprecated as of CUDA 3.1. Device emulation mode was"
			"\nremoved with the CUDA 3.1 release.";

		case cudaErrorSynchronizationError: return

			"CUDA ERROR"
			"\nThis indicated that a synchronization operation had failed."
			"\nThis was previously used for some device emulation functions."
			"\ndeprecated"
			"\nThis error return is deprecated as of CUDA 3.1. Device emulation mode was"
			"\nremoved with the CUDA 3.1 release.";

		case cudaErrorInvalidFilterSetting: return

			"CUDA ERROR"
			"\nThis indicates that a non-float texture was being accessed with linear"
			"\nfiltering. This is not supported by CUDA.";

		case cudaErrorInvalidNormSetting: return

			"CUDA ERROR"
			"\nThis indicates that an attempt was made to read a non-float texture as a"
			"\nnormalized float. This is not supported by CUDA.";

		case cudaErrorMixedDeviceExecution: return

			"CUDA ERROR"
			"\nMixing of device and device emulation code was not allowed."
			"\ndeprecated"
			"\nThis error return is deprecated as of CUDA 3.1. Device emulation mode was"
			"\nremoved with the CUDA 3.1 release.";

		case cudaErrorCudartUnloading: return

			"CUDA ERROR"
			"\nThis indicates that a CUDA Runtime API call cannot be executed because"
			"\nit is being called during process shut down, at a point in time after"
			"\nCUDA driver has been unloaded.";

		case cudaErrorUnknown: return


			"CUDA ERROR"
			"\nThis indicates that an unknown internal error has occurred.";

		case cudaErrorNotYetImplemented: return

			"CUDA ERROR"
			"\nThis indicates that the API call is not yet implemented. Production"
			"\nreleases of CUDA will never return this error."
			"\ndeprecated"
			"\nThis error return is deprecated as of CUDA 4.1.";

		case cudaErrorMemoryValueTooLarge: return

			"CUDA ERROR"
			"\nThis indicated that an emulated device pointer exceeded the 32-bit address"
			"\nrange."
			"\ndeprecated"
			"\nThis error return is deprecated as of CUDA 3.1. Device emulation mode was"
			"\nremoved with the CUDA 3.1 release.";

		case cudaErrorInvalidResourceHandle: return

			"CUDA ERROR"
			"\nThis indicates that a resource handle passed to the API call was not"
			"\nvalid. Resource handles are opaque types like ::cudaStream_t and"
			"\n::cudaEvent_t.";

		case cudaErrorNotReady: return

			"CUDA ERROR"
			"\nThis indicates that asynchronous operations issued previously have not"
			"\ncompleted yet. This result is not actually an error, but must be indicated"
			"\ndifferently than ::cudaSuccess (which indicates completion). Calls that"
			"\nmay return this value include ::cudaEventQuery() and ::cudaStreamQuery().";

		case cudaErrorInsufficientDriver: return

			"CUDA ERROR"
			"\nThis indicates that the installed NVIDIA CUDA driver is older than the"
			"\nCUDA runtime library. This is not a supported configuration. Users should"
			"\ninstall an updated NVIDIA display driver to allow the application to run.";

		case cudaErrorSetOnActiveProcess: return

			"CUDA ERROR"
			"\nThis indicates that the user has called ::cudaSetValidDevices(),"
			"\n::cudaSetDeviceFlags(), ::cudaD3D9SetDirect3DDevice(),"
			"\n::cudaD3D10SetDirect3DDevice, ::cudaD3D11SetDirect3DDevice(), or"
			"\n::cudaVDPAUSetVDPAUDevice() after initializing the CUDA runtime by"
			"\ncalling non-device management operations (allocating memory and"
			"\nlaunching kernels are examples of non-device management operations)."
			"\nThis error can also be returned if using runtime/driver"
			"\ninteroperability and there is an existing ::CUcontext active on the"
			"\nhost thread.";

		case cudaErrorInvalidSurface: return

			"CUDA ERROR"
			"\nThis indicates that the surface passed to the API call is not a valid"
			"\nsurface.";

		case cudaErrorNoDevice: return

			"CUDA ERROR"
			"\nThis indicates that no CUDA-capable devices were detected by the installed"
			"\nCUDA driver.";

		case cudaErrorECCUncorrectable: return

			"CUDA ERROR"
			"\nThis indicates that an uncorrectable ECC error was detected during"
			"\nexecution.";

		case cudaErrorSharedObjectSymbolNotFound: return

			"CUDA ERROR"
			"\nThis indicates that a link to a shared object failed to resolve.";

		case cudaErrorSharedObjectInitFailed: return


			"CUDA ERROR"
			"\nThis indicates that initialization of a shared object failed.";

		case cudaErrorUnsupportedLimit: return

			"CUDA ERROR"
			"\nThis indicates that the ::cudaLimit passed to the API call is not"
			"\nsupported by the active device.";

		case cudaErrorDuplicateVariableName: return

			"CUDA ERROR"
			"\nThis indicates that multiple global or constant variables (across separate"
			"\nCUDA source files in the application) share the same string name.";

		case cudaErrorDuplicateTextureName: return

			"CUDA ERROR"
			"\nThis indicates that multiple textures (across separate CUDA source"
			"\nfiles in the application) share the same string name.";

		case cudaErrorDuplicateSurfaceName: return

			"CUDA ERROR"
			"\nThis indicates that multiple surfaces (across separate CUDA source"
			"\nfiles in the application) share the same string name.";

		case cudaErrorDevicesUnavailable: return

			"CUDA ERROR"
			"\nThis indicates that all CUDA devices are busy or unavailable at the current"
			"\ntime. Devices are often busy/unavailable due to use of"
			"\n::cudaComputeModeExclusive, ::cudaComputeModeProhibited or when long"
			"\nrunning CUDA kernels have filled up the GPU and are blocking new work"
			"\nfrom starting. They can also be unavailable due to memory constraints"
			"\non a device that already has active CUDA work being performed.";

		case cudaErrorInvalidKernelImage: return

			"CUDA ERROR"
			"\nThis indicates that the device kernel image is invalid.";

		case cudaErrorNoKernelImageForDevice: return

			"CUDA ERROR"
			"\nThis indicates that there is no kernel image available that is suitable"
			"\nfor the device. This can occur when a user specifies code generation"
			"\noptions for a particular CUDA source file that do not include the"
			"\ncorresponding device configuration.";

		case cudaErrorIncompatibleDriverContext: return

			"CUDA ERROR"
			"\nThis indicates that the current context is not compatible with this"
			"\nthe CUDA Runtime. This can only occur if you are using CUDA"
			"\nRuntime/Driver interoperability and have created an existing Driver"
			"\ncontext using the driver API. The Driver context may be incompatible"
			"\neither because the Driver context was created using an older version "
			"\nof the API, because the Runtime API call expects a primary driver "
			"\ncontext and the Driver context is not primary, or because the Driver "
			"\ncontext has been destroyed. Please see ref CUDART_DRIVER \"Interactions "
			"\nwith the CUDA Driver API\" for more information.";

		case cudaErrorPeerAccessAlreadyEnabled: return

			"CUDA ERROR"
			"\nThis error indicates that a call to ::cudaDeviceEnablePeerAccess() is"
			"\ntrying to re-enable peer addressing on from a context which has already"
			"\nhad peer addressing enabled.";

		case cudaErrorPeerAccessNotEnabled: return


			"CUDA ERROR"
			"\nThis error indicates that ::cudaDeviceDisablePeerAccess() is trying to "
			"\ndisable peer addressing which has not been enabled yet via "
			"\n::cudaDeviceEnablePeerAccess().";

		case cudaErrorDeviceAlreadyInUse: return


			"CUDA ERROR"
			"\nThis indicates that a call tried to access an exclusive-thread device that "
			"\nis already in use by a different thread.";

		case cudaErrorProfilerDisabled: return

			"CUDA ERROR"
			"\nThis indicates profiler is not initialized for this run. This can"
			"\nhappen when the application is running with external profiling tools"
			"\nlike visual profiler.";

		case cudaErrorProfilerNotInitialized: return

			"CUDA ERROR"
			"\ndeprecated"
			"\nThis error return is deprecated as of CUDA 5.0. It is no longer an error"
			"\nto attempt to enable/disable the profiling via ::cudaProfilerStart or"
			"\n::cudaProfilerStop without initialization.";

		case cudaErrorProfilerAlreadyStarted: return


			"CUDA ERROR"
			"\ndeprecated"
			"\nThis error return is deprecated as of CUDA 5.0. It is no longer an error"
			"\nto call cudaProfilerStart() when profiling is already enabled.";

		case cudaErrorProfilerAlreadyStopped: return


			"CUDA ERROR"
			"\ndeprecated"
			"\nThis error return is deprecated as of CUDA 5.0. It is no longer an error"
			"\nto call cudaProfilerStop() when profiling is already disabled.";

		case cudaErrorAssert: return

			"CUDA ERROR"
			"\nAn assert triggered in device code during kernel execution. The device"
			"\ncannot be used again until ::cudaThreadExit() is called. All existing "
			"\nallocations are invalid and must be reconstructed if the program is to"
			"\ncontinue using CUDA. ";

		case cudaErrorTooManyPeers: return

			"CUDA ERROR"
			"\nThis error indicates that the hardware resources required to enable"
			"\npeer access have been exhausted for one or more of the devices "
			"\npassed to ::cudaEnablePeerAccess().";

		case cudaErrorHostMemoryAlreadyRegistered: return

			"CUDA ERROR"
			"\nThis error indicates that the memory range passed to ::cudaHostRegister()"
			"\nhas already been registered.";

		case cudaErrorHostMemoryNotRegistered: return

			"CUDA ERROR"
			"\nThis error indicates that the pointer passed to ::cudaHostUnregister()"
			"\ndoes not correspond to any currently registered memory region.";

		case cudaErrorOperatingSystem: return

			"CUDA ERROR"
			"\nThis error indicates that an OS call failed.";

		case cudaErrorPeerAccessUnsupported: return

			"CUDA ERROR"
			"\nThis error indicates that P2P access is not supported across the given"
			"\ndevices.";

		case cudaErrorLaunchMaxDepthExceeded: return

			"CUDA ERROR"
			"\nThis error indicates that a device runtime grid launch did not occur "
			"\nbecause the depth of the child grid would exceed the maximum supported"
			"\nnumber of nested grid launches. ";

		case cudaErrorLaunchFileScopedTex: return


			"CUDA ERROR"
			"\nThis error indicates that a grid launch did not occur because the kernel "
			"\nuses file-scoped textures which are unsupported by the device runtime. "
			"\nKernels launched via the device runtime only support textures created with "
			"\nthe Texture Object API's.";

		case cudaErrorLaunchFileScopedSurf: return

			"CUDA ERROR"
			"\nThis error indicates that a grid launch did not occur because the kernel "
			"\nuses file-scoped surfaces which are unsupported by the device runtime."
			"\nKernels launched via the device runtime only support surfaces created with"
			"\nthe Surface Object API's.";

		case cudaErrorSyncDepthExceeded: return

			"CUDA ERROR"
			"\nThis error indicates that a call to ::cudaDeviceSynchronize made from"
			"\nthe device runtime failed because the call was made at grid depth greater"
			"\nthan than either the default (2 levels of grids) or user specified device "
			"\nlimit ::cudaLimitDevRuntimeSyncDepth. To be able to synchronize on "
			"\nlaunched grids at a greater depth successfully, the maximum nested "
			"\ndepth at which ::cudaDeviceSynchronize will be called must be specified "
			"\nwith the ::cudaLimitDevRuntimeSyncDepth limit to the ::cudaDeviceSetLimit"
			"\napi before the host-side launch of a kernel using the device runtime. "
			"\nKeep in mind that additional levels of sync depth require the runtime "
			"\nto reserve large amounts of device memory that cannot be used for "
			"\nuser allocations.";

		case cudaErrorLaunchPendingCountExceeded: return

			"CUDA ERROR"
			"\nThis error indicates that a device runtime grid launch failed because"
			"\nthe launch would exceed the limit ::cudaLimitDevRuntimePendingLaunchCount."
			"\nFor this launch to proceed successfully, ::cudaDeviceSetLimit must be"
			"\ncalled to set the ::cudaLimitDevRuntimePendingLaunchCount to be higher "
			"\nthan the upper bound of outstanding launches that can be issued to the"
			"\ndevice runtime. Keep in mind that raising the limit of pending device"
			"\nruntime launches will require the runtime to reserve device memory that"
			"\ncannot be used for user allocations.";

		case cudaErrorNotPermitted: return
			"CUDA ERROR"
			"\nThis error indicates the attempted operation is not permitted.";

		case cudaErrorNotSupported: return
			"CUDA ERROR"
			"\nThis error indicates the attempted operation is not supported"
			"\non the current system or device.";

		case cudaErrorStartupFailure: return

			"CUDA ERROR"
			"\nThis indicates an internal startup failure in the CUDA runtime.";

		case cudaErrorApiFailureBase: return

			"CUDA ERROR"
			"\nAny unhandled CUDA driver error is added to this value and returned via"
			"\nthe runtime. Production releases of CUDA should not return such errors."
			"\ndeprecated"
			"\nThis error return is deprecated as of CUDA 4.1.";
		default:
			return "CUDA ERROR: Unknown error ";
	}

};

}
