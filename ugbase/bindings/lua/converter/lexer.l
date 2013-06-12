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

">="            return GE;
"<="            return LE;
"=="            return EQ;
"~="            return NE;
"end"			return END;
"local"			return LOCAL;
"function"		return FUNCTION;
"return"		return RETURN;

"if"            return IF;
"then"          return THEN;
"else"          return ELSE;
"elseif"        return ELSEIF;

"for"			return TK_FOR;
"do"			return TK_DO;
"break"			return TK_BREAK;


"and"			return AND;
"or"			return OR;

"math.cos"		return MATH_COS;
"math.sin"		return MATH_SIN;
"math.exp"		return MATH_EXP;
"math.abs"		return MATH_ABS;
"math.log"		return MATH_LOG;
"math.log10"    return MATH_LOG10;
"math.sqrt"     return MATH_SQRT;
"math.floor"    return MATH_FLOOR;
"math.ceil"     return MATH_CEIL;

"math.pow"      return MATH_POW;
"math.max"      return MATH_MAX;
"math.min"      return MATH_MIN;

"math.pi"       return MATH_PI;


{name}       { 
                yylval.sIndex = globalP->get_id_for_name(yytext);
                return VARIABLE;
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
