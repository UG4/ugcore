cmake_minimum_required(VERSION 2.6)

####
# liblgm Project
####

project(LIBLGM)

# build liblgm library
add_library(lgm
	include/lgm.h
	include/lgm_info.h
	include/lgm_line.h
	include/lgm_surface.h
	src/lgm.c
	src/lgm_error.h
	src/lgm_error.c
	src/lgm_info.c
	src/lgm_parser.h
	src/lgm_parser.c
	src/lgm_writer.c
	
	# tokstream
	src/tokstream/tokstream.h
	src/tokstream/tokstream.c
)
