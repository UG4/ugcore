%{
#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>

#include "parser_node.h"
#include "pclass.h"

#define YYERROR_VERBOSE
/* prototypes */
using namespace ug;
int yylex(void);
typedef struct yy_buffer_state *YY_BUFFER_STATE;
YY_BUFFER_STATE yy_scan_string (const char *yy_str  );

pclass *globalP;

void yyerror(const char *s);
%}

%union {
    double iValue;                 /* integer value */
    int sIndex;                /* symbol table index */
    nodeType *nPtr;             /* node pointer */
	const char* pFunction;
};



%token <iValue> YY_INTEGER
%token <sIndex> VARIABLE
%token WHILE IF THEN END ELSE ELSEIF LOCAL FUNCTION RETURN
%token MATH_SIN MATH_COS MATH_EXP MATH_ABS MATH_LOG MATH_LOG10 MATH_SQRT
%token MATH_FLOOR MATH_CEIL MATH_POW 
%token MATH_MAX MATH_MIN 
%token MATH_PI TK_DO TK_FOR TK_BREAK
%nonassoc ELSE

%left AND OR
%left GE LE EQ NE '>' '<'
%left '+' '-'
%left '*' '/'

%nonassoc UMINUS

%type <nPtr> stmt expr stmt_list funcstat arguments exprlist arguments2 elseblock forStep

%%

program:
        func                { }
        ;

func:
          funcstat			{ globalP->add($1); }
        ;


stmt:
          expr							{ $$ = $1; }
		| RETURN exprlist				{ $$ = globalP->opr1('R', $2);}
		| LOCAL VARIABLE '=' expr		{ $$ = globalP->opr2('=', globalP->id($2), $4); globalP->set_local($2); }
		| LOCAL VARIABLE				{ $$ = NULL; globalP->set_local($2); }        
		| VARIABLE '=' expr				{ $$ = globalP->opr2('=', globalP->id($1), $3); globalP->assert_local($1); }
        | IF expr THEN stmt_list elseblock END	{ $$ = globalP->opr3(IF, $2, $4, $5); }        
		| TK_FOR VARIABLE '=' expr ',' expr forStep TK_DO stmt_list END { $$ = globalP->forOp(globalP->id($2), $4, $6, $7, $9); }
		| TK_BREAK						{ $$ = globalP->opr0(TK_BREAK); }
        ;

forStep:
		',' expr						{ $$ = $2; }
		|								{ $$ = globalP->con(1); }
		;

elseblock:
        ELSEIF expr THEN stmt_list elseblock { $$ = globalP->opr3(ELSEIF, $2, $4, $5); }
        | ELSE stmt_list                    { $$ = globalP->opr1(ELSE, $2); }
        |                                   { $$ = NULL; }
        ;

funcstat:
		FUNCTION VARIABLE '(' arguments ')' stmt_list END         { globalP->set_name($2); globalP->set_arguments($4); $$ = $6; }
		;


arguments:
		VARIABLE ',' arguments			{ $$ = globalP->opr2(',', globalP->id($1), $3); }
		| VARIABLE						{ $$ = globalP->id($1); }
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
        | VARIABLE '(' arguments2 ')' { $$ = globalP->function(globalP->id($1), $3); }
        | VARIABLE              { $$ = globalP->id($1); }
        | '-' expr %prec UMINUS { $$ = globalP->opr1(UMINUS, $2); }
        | expr '+' expr         { $$ = globalP->opr2('+', $1, $3); }
        | expr '-' expr         { $$ = globalP->opr2('-', $1, $3); }
        | expr '*' expr         { $$ = globalP->opr2('*', $1, $3); }
        | expr '/' expr         { $$ = globalP->opr2('/', $1, $3); }
        | expr '<' expr         { $$ = globalP->opr2('<', $1, $3); }
        | expr '>' expr         { $$ = globalP->opr2('>', $1, $3); }
        | expr GE expr          { $$ = globalP->opr2(GE, $1, $3); }
        | expr LE expr          { $$ = globalP->opr2(LE, $1, $3); }
        | expr NE expr          { $$ = globalP->opr2(NE, $1, $3); }
        | expr EQ expr          { $$ = globalP->opr2(EQ, $1, $3); }
		| expr AND expr          { $$ = globalP->opr2(AND, $1, $3); }
		| expr OR expr          { $$ = globalP->opr2(OR, $1, $3); }
        | '(' expr ')'          { $$ = $2; }
		| MATH_COS '(' expr ')' { $$ = globalP->opr1(MATH_COS, $3); }
		| MATH_SIN '(' expr ')' { $$ = globalP->opr1(MATH_SIN, $3); }
		| MATH_EXP '(' expr ')' { $$ = globalP->opr1(MATH_EXP, $3); }
        | MATH_ABS '(' expr ')' { $$ = globalP->opr1(MATH_ABS, $3); }
        | MATH_LOG '(' expr ')' { $$ = globalP->opr1(MATH_LOG, $3); }
        | MATH_LOG10 '(' expr ')' { $$ = globalP->opr1(MATH_LOG10, $3); }
        | MATH_SQRT '(' expr ')' { $$ = globalP->opr1(MATH_SQRT, $3); }
        | MATH_FLOOR '(' expr ')' { $$ = globalP->opr1(MATH_FLOOR, $3); }
        | MATH_CEIL '(' expr ')' { $$ = globalP->opr1(MATH_CEIL, $3); }
        | MATH_POW '(' expr ',' expr ')' { $$ = globalP->opr2(MATH_CEIL, $3, $5); }
        | MATH_MIN '(' expr ',' expr ')' { $$ = globalP->opr2(MATH_MIN, $3, $5); }
        | MATH_MAX '(' expr ',' expr ')' { $$ = globalP->opr2(MATH_MAX, $3, $5); }
        | MATH_PI                  { $$ = globalP->opr0(MATH_PI); }
		
        ;

%%


void yyerror(const char *s)
{
	globalP->err << "error: " << s << "\n";
}

void yaccparse(const char*command, ug::pclass *p)
{	
	globalP = p;
	if(command!=NULL)
		yy_scan_string(command);
    yyparse();
}