cmake_minimum_required(VERSION 2.6)

####
# libng Project
####

project(LIBNG)

# build libng library
add_library(ng
	include/ng.h
	include/ng_info.h
	include/ng_node.h
	include/ng_element.h
	src/ng.c
	src/ng_error.h
	src/ng_error.c
	src/ng_info.c
	src/ng_parser.h
	src/ng_parser.c
	src/ng_writer.c
	
	# tokstream
	src/tokstream/tokstream.h
	src/tokstream/tokstream.c
)
