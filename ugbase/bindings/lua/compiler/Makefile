# note: this is the makefile for parser.y -> parser.cpp and lexer.l -> lexer.cpp
all: parser lexer

parser: parser.y
	bison -p ug4_lua_YY_ -d -y parser.y -o parser.cpp

lexer: lexer.l
	flex -P ug4_lua_YY_ -o lexer.cpp lexer.l

