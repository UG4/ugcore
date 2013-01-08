CC      = gcc
CFLAGS  =
LDFLAGS =

OBJ = lexer.o parser.o main.o pclass.o

prog: $(OBJ)
	g++ -g $(CFLAGS) -o prog $(OBJ) $(LDFLAGS)

parser.cpp: parser.y
	bison -d -y parser.y -o parser.cpp

lexer.cpp: lexer.l
	flex -o lexer.cpp lexer.l

%.o: %.c
	$(CC) -g $(CFLAGS) -c $<

%.o: %.cpp
	g++ -g -c $<





