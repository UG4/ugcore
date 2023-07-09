# Configuration for RPMs

set(CPACK_GENERATOR "Bundle") # For MacOS

#  The name of the generated bundle. This appears in the macOS Finder as the bundle name. Required.
set(CPACK_BUNDLE_NAME "UG4")

#Path to an macOS Property List (.plist) file that will be used for the generated bundle. This assumes that the caller has generated or specified their own Info.plist file. Required.
#CPACK_BUNDLE_PLIST

# Path to an macOS icon file that will be used as the icon for the generated bundle. This is the icon that appears in the macOS Finder for the bundle, and in the macOS dock when the bundle is opened. Required.
#CPACK_BUNDLE_ICON

# Path to a startup script. This is a path to an executable or script that will be run whenever an end-user double-clicks the generated bundle in the macOS Finder. Optional.
##CPACK_BUNDLE_STARTUP_COMMAND

# New in version 3.2.
# The name of your Apple supplied code signing certificate for the application. The name usually takes the form Developer ID Application: [Name] or 3rd Party Mac Developer Application: [Name]. If this variable is not set the application will not be signed.
#CPACK_BUNDLE_APPLE_CERT_APP

# The name of the Property List (.plist) file that contains your Apple entitlements for sandboxing your application. This file is required for submission to the macOS App Store.
#CPACK_BUNDLE_APPLE_ENTITLEMENTS

   

# A list of additional files that you wish to be signed. You do not need to list the main application folder, or the main executable. You should list any frameworks and plugins that are included in your app bundle.
#CPACK_BUNDLE_APPLE_CODESIGN_FILES

 
# Additional parameter that will passed to codesign. Default value: --deep -f   
#CPACK_BUNDLE_APPLE_CODESIGN_PARAMETER

    
# Path to the codesign(1) command used to sign applications with an Apple cert. 
#CPACK_COMMAND_CODESIGN

   
