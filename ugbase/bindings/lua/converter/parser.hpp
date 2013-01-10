/* A Bison parser, made by GNU Bison 2.3.  */

/* Skeleton interface for Bison's Yacc-like parsers in C

   Copyright (C) 1984, 1989, 1990, 2000, 2001, 2002, 2003, 2004, 2005, 2006
   Free Software Foundation, Inc.

   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2, or (at your option)
   any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program; if not, write to the Free Software
   Foundation, Inc., 51 Franklin Street, Fifth Floor,
   Boston, MA 02110-1301, USA.  */

/* As a special exception, you may create a larger work that contains
   part or all of the Bison parser skeleton and distribute that work
   under terms of your choice, so long as that work isn't itself a
   parser generator using the skeleton or a modified version thereof
   as a parser skeleton.  Alternatively, if you modify or redistribute
   the parser skeleton itself, you may (at your option) remove this
   special exception, which will cause the skeleton and the resulting
   Bison output files to be licensed under the GNU General Public
   License without this special exception.

   This special exception was added by the Free Software Foundation in
   version 2.2 of Bison.  */

/* Tokens.  */
#ifndef YYTOKENTYPE
# define YYTOKENTYPE
   /* Put the tokens into the symbol table, so that GDB and other debuggers
      know about them.  */
   enum yytokentype {
     YY_INTEGER = 258,
     VARIABLE = 259,
     WHILE = 260,
     IF = 261,
     THEN = 262,
     END = 263,
     ELSE = 264,
     ELSEIF = 265,
     LOCAL = 266,
     FUNCTION = 267,
     RETURN = 268,
     MATH_SIN = 269,
     MATH_COS = 270,
     MATH_EXP = 271,
     MATH_ABS = 272,
     MATH_LOG = 273,
     MATH_LOG10 = 274,
     MATH_SQRT = 275,
     MATH_FLOOR = 276,
     MATH_CEIL = 277,
     MATH_POW = 278,
     MATH_MAX = 279,
     MATH_MIN = 280,
     MATH_PI = 281,
     OR = 282,
     AND = 283,
     NE = 284,
     EQ = 285,
     LE = 286,
     GE = 287,
     UMINUS = 288
   };
#endif
/* Tokens.  */
#define YY_INTEGER 258
#define VARIABLE 259
#define WHILE 260
#define IF 261
#define THEN 262
#define END 263
#define ELSE 264
#define ELSEIF 265
#define LOCAL 266
#define FUNCTION 267
#define RETURN 268
#define MATH_SIN 269
#define MATH_COS 270
#define MATH_EXP 271
#define MATH_ABS 272
#define MATH_LOG 273
#define MATH_LOG10 274
#define MATH_SQRT 275
#define MATH_FLOOR 276
#define MATH_CEIL 277
#define MATH_POW 278
#define MATH_MAX 279
#define MATH_MIN 280
#define MATH_PI 281
#define OR 282
#define AND 283
#define NE 284
#define EQ 285
#define LE 286
#define GE 287
#define UMINUS 288




#if ! defined YYSTYPE && ! defined YYSTYPE_IS_DECLARED
typedef union YYSTYPE
#line 21 "parser.y"
{
    double iValue;                 /* integer value */
    int sIndex;                /* symbol table index */
    nodeType *nPtr;             /* node pointer */
	const char* pFunction;
}
/* Line 1529 of yacc.c.  */
#line 122 "parser.hpp"
	YYSTYPE;
# define yystype YYSTYPE /* obsolescent; will be withdrawn */
# define YYSTYPE_IS_DECLARED 1
# define YYSTYPE_IS_TRIVIAL 1
#endif

extern YYSTYPE yylval;

