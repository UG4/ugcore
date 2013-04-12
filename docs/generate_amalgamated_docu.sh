#!/bin/bash

doxygen - < doxy_config_ugbase_tags.txt && doxygen - < doxy_config_plugins.txt && doxygen - < doxy_config_apps.txt && doxygen - < doxy_config_ug4.txt && mkdir ug4/html/apps ug4/html/plugins && mv apps/html/* ug4/html/apps/. && mv plugins/html/* ug4/html/plugins/.
