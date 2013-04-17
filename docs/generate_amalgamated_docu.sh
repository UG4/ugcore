#!/bin/bash

doxylog="$1"

# remove old html and tags
function cleanup_old_docu {
	echo "Step 1/7: Removing old docu if existend"
	rm -rf ug4/html ug4/plugins/html ug4/apps/html *.tags &> $doxylog
}

# generate tags for ugbase
function generate_ugbase {
	echo "Step 2/7: Generating tags for ugbase"
	doxygen - < doxy_config_ugbase_tags.txt 1>> $doxylog 2>> $doxylog
}

# generate tags and html for plugins
function generate_plugins {
	echo "Step 3/7: Generating tags and html for plugins"
	doxygen - < doxy_config_plugins.txt 1>> $doxylog 2>> $doxylog
}

# generate html for apps
function generate_apps {
	echo "Step 4/7: Generating tags and html for apps"
	doxygen - < doxy_config_apps.txt 1>> $doxylog 2>> $doxylog
}

# generate html for ugbase
function generate_ug4 {
	echo "Step 5/7: Generating html for ugbase"
	doxygen - < doxy_config_ug4_mathjax.txt 1>> $doxylog 2>> $doxylog
}

# prepare amalgamation of html
function prepare_amalgamate {
	echo "Step 6/7: Prepare amalgamation"
	mkdir -p ug4/html/apps ug4/html/plugins 1>> $doxylog 2>> $doxylog
}

# move all html into a single directory
function amalgamate {
	echo "Step 7/7: Amalgamate html of ubgase, plugins and apps"
	mv -fu apps/html/* ug4/html/apps/. && mv -fu plugins/html/* ug4/html/plugins/. 1>> $doxylog 2>> $doxylog
}

# stitch all steps into a single command list
# (if any fails, all fail)
cleanup_old_docu \
&& generate_ugbase \
&& generate_plugins \
&& generate_apps \
&& generate_ug4 \
&& prepare_amalgamate \
&& amalgamate
