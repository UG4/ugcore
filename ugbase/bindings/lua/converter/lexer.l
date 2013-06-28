/* 
 * \file	lexer.l
 * \author	Martin Rupp
 *
 * Created on 20. November 2012, 10:16
 * 
 * use with
 * 		flex -o lexer.cpp lexer.l
 * to generate lexer.cpp
 */

%{
#include <stdlib.h>
#include "lua_parser_class.h"
#include "parser.hpp"
using namespace ug;
void yyerror(const char *);
using namespace ug;
extern LUAParserClass *globalP;
%}

w              [ \t\v\a]+
o              [ \t\v\a]*
name           [_a-zA-Z][_a-zA-Z0-9]*
n              [0-9]+
exp            [Ee][+-]?{n}
number         ({n}|{n}[.]{n}){exp}?

%%

--[^\n]*    ;       /* lua comment*/

0           {
                yylval.iValue = atof(yytext);
                return YY_INTEGER;
            }

{number}    {
                yylval.iValue = atof(yytext);
                return YY_INTEGER;                
            }

[-()<>=+*/,] {
                return *yytext;
             }

">="            return LUAPARSER_GE;
"<="            return LUAPARSER_LE;
"=="            return LUAPARSER_EQ;
"~="            return LUAPARSER_NE;
"end"			return LUAPARSER_END;
"local"			return LUAPARSER_LOCAL;
"function"		return LUAPARSER_FUNCTION;
"return"		return LUAPARSER_RETURN;

"if"            return LUAPARSER_IF;
"then"          return LUAPARSER_THEN;
"else"          return LUAPARSER_ELSE;
"elseif"        return LUAPARSER_ELSEIF;

"for"			return LUAPARSER_FOR;
"do"			return LUAPARSER_DO;
"break"			return LUAPARSER_BREAK;


"and"			return LUAPARSER_AND;
"or"			return LUAPARSER_OR;

"math.cos"		return LUAPARSER_MATH_COS;
"math.sin"		return LUAPARSER_MATH_SIN;
"math.exp"		return LUAPARSER_MATH_EXP;
"math.abs"		return LUAPARSER_MATH_ABS;
"math.log"		return LUAPARSER_MATH_LOG;
"math.log10"    return LUAPARSER_MATH_LOG10;
"math.sqrt"     return LUAPARSER_MATH_SQRT;
"math.floor"    return LUAPARSER_MATH_FLOOR;
"math.ceil"     return LUAPARSER_MATH_CEIL;

"math.pow"      return LUAPARSER_MATH_POW;
"math.max"      return LUAPARSER_MATH_MAX;
"math.min"      return LUAPARSER_MATH_MIN;

"math.pi"       return LUAPARSER_MATH_PI;


{name}       { 
                yylval.sIndex = globalP->get_id_for_name(yytext);
                return LUAPARSER_VARIABLE;
            }

[ \t\n]+        ;       /* ignore whitespace */

.               {
                    char buf[255];
                    sprintf(buf, "Unknown character '%s' at line %d", yytext, yyget_lineno());
                    yyerror(buf);
                }

%%

int yywrap(void) {
    return 1;
// remove some warnings
if(0) { unput(0); yyinput(); }
}
