/*
 * \file	parser.y
 * \author	Martin Rupp
 *
 * Created on 20. November 2012, 10:16
 * 
 * use with
 * 		bison -d -y parser.y -o parser.cpp
 * to generate parser.cpp
 */

%{
#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>

#include "parser_node.h"
#include "lua_parser_class.h"

#define YYERROR_VERBOSE
/* prototypes */
using namespace ug;
int yylex(void);
typedef struct yy_buffer_state *YY_BUFFER_STATE;
YY_BUFFER_STATE yy_scan_string (const char *yy_str  );

LUAParserClass *globalP;

void yyerror(const char *s);
%}

%union {
    double iValue;                 /* integer value */
    int sIndex;                /* symbol table index */
    nodeType *nPtr;             /* node pointer */
	const char* pFunction;
};



%token <iValue> YY_INTEGER
%token <sIndex> LUAPARSER_VARIABLE
%token LUAPARSER_WHILE LUAPARSER_IF LUAPARSER_THEN LUAPARSER_END LUAPARSER_ELSE LUAPARSER_ELSEIF LUAPARSER_LOCAL LUAPARSER_FUNCTION LUAPARSER_RETURN
%token LUAPARSER_MATH_SIN LUAPARSER_MATH_COS LUAPARSER_MATH_EXP LUAPARSER_MATH_ABS LUAPARSER_MATH_LOG LUAPARSER_MATH_LOG10 LUAPARSER_MATH_SQRT
%token LUAPARSER_MATH_FLOOR LUAPARSER_MATH_CEIL LUAPARSER_MATH_POW 
%token LUAPARSER_MATH_MAX LUAPARSER_MATH_MIN 
%token LUAPARSER_MATH_PI LUAPARSER_DO LUAPARSER_FOR LUAPARSER_BREAK
%nonassoc LUAPARSER_ELSE

%left LUAPARSER_AND LUAPARSER_OR
%left LUAPARSER_GE LUAPARSER_LE LUAPARSER_EQ LUAPARSER_NE '>' '<'
%left '+' '-'
%left '*' '/'

%nonassoc LUAPARSER_UMINUS

%type <nPtr> stmt expr stmt_list funcstat arguments exprlist arguments2 elseblock forStep

%%

program:
        func                { }
        ;

func:
          funcstat			{ globalP->add($1); }
        ;


stmt:
          expr												{ $$ = $1; }
		| LUAPARSER_RETURN exprlist							{ $$ = globalP->opr1('R', $2);}
		| LUAPARSER_LOCAL LUAPARSER_VARIABLE '=' expr		{ $$ = globalP->opr2('=', globalP->id($2), $4); globalP->set_local($2); }
		| LUAPARSER_LOCAL LUAPARSER_VARIABLE				{ $$ = NULL; globalP->set_local($2); }        
		| LUAPARSER_VARIABLE '=' expr						{ $$ = globalP->opr2('=', globalP->id($1), $3); globalP->assert_local($1); }
        | LUAPARSER_IF expr LUAPARSER_THEN stmt_list elseblock LUAPARSER_END			
        													{ $$ = globalP->opr3(LUAPARSER_IF, $2, $4, $5); }        
		| LUAPARSER_FOR LUAPARSER_VARIABLE '=' expr ',' expr forStep LUAPARSER_DO stmt_list LUAPARSER_END 
															{ $$ = globalP->forOp(globalP->id($2), $4, $6, $7, $9); }
		| LUAPARSER_BREAK									{ $$ = globalP->opr0(LUAPARSER_BREAK); }
        ;

forStep:
		',' expr											{ $$ = $2; }
		|													{ $$ = globalP->con(1); }
		;

elseblock:
        LUAPARSER_ELSEIF expr LUAPARSER_THEN stmt_list elseblock 
        													{ $$ = globalP->opr3(LUAPARSER_ELSEIF, $2, $4, $5); }
        | LUAPARSER_ELSE stmt_list                   		{ $$ = globalP->opr1(LUAPARSER_ELSE, $2); }
        |					                                { $$ = NULL; }
        ;

funcstat:
		LUAPARSER_FUNCTION LUAPARSER_VARIABLE '(' arguments ')' stmt_list LUAPARSER_END
													         { globalP->set_name($2); globalP->set_arguments($4); $$ = $6; }
		;


arguments:
		LUAPARSER_VARIABLE ',' arguments					{ $$ = globalP->opr2(',', globalP->id($1), $3); }
		| LUAPARSER_VARIABLE								{ $$ = globalP->id($1); }
		;

arguments2:
		expr',' arguments2			{ $$ = globalP->opr2(',', $1, $3); }
		| expr						{ $$ = $1; }
		;


exprlist:
		expr ',' exprlist			{ $$ = globalP->opr2(',', $1, $3); }
		| expr						{ $$ = $1; }
		;

stmt_list:
          stmt                  { $$ = $1; }
        |  stmt stmt_list { $$ = globalP->opr2(';', $1, $2); }
        ;

expr:
          YY_INTEGER               { $$ = globalP->con($1); }          
        | LUAPARSER_VARIABLE '(' arguments2 ')' 	{ $$ = globalP->function(globalP->id($1), $3); }
        | LUAPARSER_VARIABLE             			{ $$ = globalP->id($1); }
        | '-' expr %prec LUAPARSER_UMINUS { $$ = globalP->opr1(LUAPARSER_UMINUS, $2); }
        | expr '+' expr         { $$ = globalP->opr2('+', $1, $3); }
        | expr '-' expr         { $$ = globalP->opr2('-', $1, $3); }
        | expr '*' expr         { $$ = globalP->opr2('*', $1, $3); }
        | expr '/' expr         { $$ = globalP->opr2('/', $1, $3); }
        | expr '<' expr         { $$ = globalP->opr2('<', $1, $3); }
        | expr '>' expr         { $$ = globalP->opr2('>', $1, $3); }
        | expr LUAPARSER_GE expr          { $$ = globalP->opr2(LUAPARSER_GE, $1, $3); }
        | expr LUAPARSER_LE expr          { $$ = globalP->opr2(LUAPARSER_LE, $1, $3); }
        | expr LUAPARSER_NE expr          { $$ = globalP->opr2(LUAPARSER_NE, $1, $3); }
        | expr LUAPARSER_EQ expr          { $$ = globalP->opr2(LUAPARSER_EQ, $1, $3); }
		| expr LUAPARSER_AND expr          { $$ = globalP->opr2(LUAPARSER_AND, $1, $3); }
		| expr LUAPARSER_OR expr          { $$ = globalP->opr2(LUAPARSER_OR, $1, $3); }
        | '(' expr ')'          { $$ = $2; }
		| LUAPARSER_MATH_COS '(' expr ')' { $$ = globalP->opr1(LUAPARSER_MATH_COS, $3); }
		| LUAPARSER_MATH_SIN '(' expr ')' { $$ = globalP->opr1(LUAPARSER_MATH_SIN, $3); }
		| LUAPARSER_MATH_EXP '(' expr ')' { $$ = globalP->opr1(LUAPARSER_MATH_EXP, $3); }
        | LUAPARSER_MATH_ABS '(' expr ')' { $$ = globalP->opr1(LUAPARSER_MATH_ABS, $3); }
        | LUAPARSER_MATH_LOG '(' expr ')' { $$ = globalP->opr1(LUAPARSER_MATH_LOG, $3); }
        | LUAPARSER_MATH_LOG10 '(' expr ')' { $$ = globalP->opr1(LUAPARSER_MATH_LOG10, $3); }
        | LUAPARSER_MATH_SQRT '(' expr ')' { $$ = globalP->opr1(LUAPARSER_MATH_SQRT, $3); }
        | LUAPARSER_MATH_FLOOR '(' expr ')' { $$ = globalP->opr1(LUAPARSER_MATH_FLOOR, $3); }
        | LUAPARSER_MATH_CEIL '(' expr ')' { $$ = globalP->opr1(LUAPARSER_MATH_CEIL, $3); }
        | LUAPARSER_MATH_POW '(' expr ',' expr ')' { $$ = globalP->opr2(LUAPARSER_MATH_CEIL, $3, $5); }
        | LUAPARSER_MATH_MIN '(' expr ',' expr ')' { $$ = globalP->opr2(LUAPARSER_MATH_MIN, $3, $5); }
        | LUAPARSER_MATH_MAX '(' expr ',' expr ')' { $$ = globalP->opr2(LUAPARSER_MATH_MAX, $3, $5); }
        | LUAPARSER_MATH_PI                  { $$ = globalP->opr0(LUAPARSER_MATH_PI); }
		
        ;

%%


void yyerror(const char *s)
{
	globalP->err << "error: " << s << "\n";
}

void yaccparse(const char*command, ug::LUAParserClass *p)
{	
	globalP = p;
	if(command!=NULL)
		yy_scan_string(command);
    yyparse();
}