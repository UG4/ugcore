#!/bin/bash

# call e.g. like this: 'generate_amalgamated_docu.sh doxylog.log 1'
# Please make sure that the DocuGen plugin is installed.
doxylog="${1:-"doxylog.log"}"
with_regdocu="${2:-0}"


# remove old html and tags
function cleanup_old_docu {
	echo "Step 1/8: Removing old docu if existent"
	if [ "$with_regdocu" -eq "1" ]; then
		rm -rf ug4/html ug4/plugins/html ug4/apps/html ug4/regdocu/html *.tags &> "${doxylog}1"
	else
		rm -rf ug4/html ug4/plugins/html ug4/apps/html *.tags &> "${doxylog}1"
	fi
}

# generate tags for ugbase
function generate_ugbase {
	echo "Step 2/8: Generating tags for ugbase"
	doxygen - < doxy_config_ugbase_tags.txt &> "${doxylog}2"
}

# generate tags and html for plugins
function generate_plugins {
	echo "Step 3/8: Generating tags and html for plugins"
	doxygen - < doxy_config_plugins.txt &> "${doxylog}3"
}

function generate_regdocu {
	if [ "$with_regdocu" -eq "1" ]; then
		echo "Step 4/8: Generating tags and html for Registry"
		rm -rf regdocu
		mkdir regdocu
		../../bin/ugshell -call GenerateScriptReferenceDocu\(\"regdocu\", true, false, true, false\) > "${doxylog}4"
		doxygen - < doxy_config_regdocu.txt &>> "${doxylog}4"
	else
		echo "Step 4/8: Skipping tags and html for Registry"
		echo "  WARNING This is not recommended as some links in the final docu will be broken!"
	fi
}

# generate html for apps
function generate_apps {
	echo "Step 5/8: Generating tags and html for apps"
	doxygen - < doxy_config_apps.txt &> "${doxylog}5"
}

# generate html for ugbase
function generate_ug4 {
	echo "Step 6/8: Generating html for ugbase"
	doxygen - < doxy_config_amalgamated.txt &> "${doxylog}6"
}

# prepare amalgamation of html
function prepare_amalgamate {
	echo "Step 7/8: Prepare amalgamation"
	if [ "$with_regdocu" -eq 1 ]; then
		mkdir -p ug4/html/apps ug4/html/plugins ug4/html/regdocu &> "${doxylog}7"
	else
		mkdir -p ug4/html/apps ug4/html/plugins &> "${doxylog}7"
	fi
}

# move all html into a single directory
function amalgamate {
	echo "Step 8/8: Amalgamate html of ubgase, plugins and apps"
	if [ "$with_regdocu" -eq "1" ]; then
		mv -fu apps/html/* ug4/html/apps/.\
			&& mv -fu plugins/html/* ug4/html/plugins/. \
			&& mv -fu regdocu/html/* ug4/html/regdocu/. &> "${doxylog}8"
	else
		mv -fu apps/html/* ug4/html/apps/.\
			&& mv -fu plugins/html/* ug4/html/plugins/. &> "${doxylog}8"
	fi
}

# stitch all steps into a single command list
# (if any fails, all fail)
cleanup_old_docu \
&& generate_ugbase \
&& generate_plugins \
&& generate_regdocu \
&& generate_apps \
&& generate_ug4 \
&& prepare_amalgamate \
&& amalgamate

#&& generate_apps \Automated performance modeling of the UG4 simulation framework

