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
     LUAPARSER_VARIABLE = 259,
     LUAPARSER_WHILE = 260,
     LUAPARSER_IF = 261,
     LUAPARSER_THEN = 262,
     LUAPARSER_END = 263,
     LUAPARSER_ELSE = 264,
     LUAPARSER_ELSEIF = 265,
     LUAPARSER_LOCAL = 266,
     LUAPARSER_FUNCTION = 267,
     LUAPARSER_RETURN = 268,
     LUAPARSER_MATH_SIN = 269,
     LUAPARSER_MATH_COS = 270,
     LUAPARSER_MATH_EXP = 271,
     LUAPARSER_MATH_ABS = 272,
     LUAPARSER_MATH_LOG = 273,
     LUAPARSER_MATH_LOG10 = 274,
     LUAPARSER_MATH_SQRT = 275,
     LUAPARSER_MATH_FLOOR = 276,
     LUAPARSER_MATH_CEIL = 277,
     LUAPARSER_MATH_POW = 278,
     LUAPARSER_MATH_MAX = 279,
     LUAPARSER_MATH_MIN = 280,
     LUAPARSER_MATH_PI = 281,
     LUAPARSER_DO = 282,
     LUAPARSER_FOR = 283,
     LUAPARSER_BREAK = 284,
     LUAPARSER_OR = 285,
     LUAPARSER_AND = 286,
     LUAPARSER_NE = 287,
     LUAPARSER_EQ = 288,
     LUAPARSER_LE = 289,
     LUAPARSER_GE = 290,
     LUAPARSER_UMINUS = 291
   };
#endif
/* Tokens.  */
#define YY_INTEGER 258
#define LUAPARSER_VARIABLE 259
#define LUAPARSER_WHILE 260
#define LUAPARSER_IF 261
#define LUAPARSER_THEN 262
#define LUAPARSER_END 263
#define LUAPARSER_ELSE 264
#define LUAPARSER_ELSEIF 265
#define LUAPARSER_LOCAL 266
#define LUAPARSER_FUNCTION 267
#define LUAPARSER_RETURN 268
#define LUAPARSER_MATH_SIN 269
#define LUAPARSER_MATH_COS 270
#define LUAPARSER_MATH_EXP 271
#define LUAPARSER_MATH_ABS 272
#define LUAPARSER_MATH_LOG 273
#define LUAPARSER_MATH_LOG10 274
#define LUAPARSER_MATH_SQRT 275
#define LUAPARSER_MATH_FLOOR 276
#define LUAPARSER_MATH_CEIL 277
#define LUAPARSER_MATH_POW 278
#define LUAPARSER_MATH_MAX 279
#define LUAPARSER_MATH_MIN 280
#define LUAPARSER_MATH_PI 281
#define LUAPARSER_DO 282
#define LUAPARSER_FOR 283
#define LUAPARSER_BREAK 284
#define LUAPARSER_OR 285
#define LUAPARSER_AND 286
#define LUAPARSER_NE 287
#define LUAPARSER_EQ 288
#define LUAPARSER_LE 289
#define LUAPARSER_GE 290
#define LUAPARSER_UMINUS 291




#if ! defined YYSTYPE && ! defined YYSTYPE_IS_DECLARED
typedef union YYSTYPE
#line 32 "parser.y"
{
    double iValue;                 /* integer value */
    int sIndex;                /* symbol table index */
    nodeType *nPtr;             /* node pointer */
	const char* pFunction;
}
/* Line 1529 of yacc.c.  */
#line 128 "parser.hpp"
	YYSTYPE;
# define yystype YYSTYPE /* obsolescent; will be withdrawn */
# define YYSTYPE_IS_DECLARED 1
# define YYSTYPE_IS_TRIVIAL 1
#endif

extern YYSTYPE yylval;

