#!/bin/bash

# remove old html and tags
function cleanup_old_docu {
	rm -rf ug4/html ug4/plugins/html ug4/apps/html *.tags
}

# generate tags for ugbase
function generate_ugbase {
	doxygen - < doxy_config_ugbase_tags.txt
}

# generate tags and html for plugins
function generate_plugins {
	doxygen - < doxy_config_plugins.txt
}

# generate html for apps
function generate_apps {
	doxygen - < doxy_config_apps.txt
}

# generate html for ugbase
function generate_ug4 {
	doxygen - < doxy_config_ug4_mathjax.txt
}

# prepare amalgamation of html
function prepare_amalgamate {
	mkdir -p ug4/html/apps ug4/html/plugins
}

# move all html into a single directory
function amalgamate {
	mv -fu apps/html/* ug4/html/apps/. && mv -fu plugins/html/* ug4/html/plugins/.
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
