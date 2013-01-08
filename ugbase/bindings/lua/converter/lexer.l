%{
#include <stdlib.h>
#include "pclass.h"
#include "parser.hpp"
using namespace ug;
void yyerror(char *);
using namespace ug;
extern pclass *globalP;
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
"if"            return IF;
"then"          return THEN;
"end"			return END;
"local"			return LOCAL;
"function"		return FUNCTION;
"return"		return RETURN;

"and"			return AND;
"or"			return OR;

"math.cos"		return MATHCOS;
"math.sin"		return MATHSIN;
"math.exp"		return MATHEXP;


{name}       { 
                yylval.sIndex = globalP->get_id_for_name(yytext);
                return VARIABLE;
            }

[ \t\n]+        ;       /* ignore whitespace */

.               yyerror("Unknown character");

%%

int yywrap(void) {
    return 1;
// remove some warnings
if(0) { unput(0); yyinput(); }
}
