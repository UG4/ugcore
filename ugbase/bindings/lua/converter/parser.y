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
%token WHILE IF THEN END LOCAL MATHSIN MATHCOS MATHEXP FUNCTION RETURN
%nonassoc ELSE

%left AND OR
%left GE LE EQ NE '>' '<'
%left '+' '-'
%left '*' '/'

%nonassoc UMINUS

%type <nPtr> stmt expr stmt_list funcstat arguments exprlist arguments2

%%

program:
        func                { }
        ;

func:
          funcstat			{ globalP->add($1); }
        ;


stmt:
          expr							{ $$ = $1; }
		| RETURN exprlist				{ $$ = globalP->opr1('R', 1, $2);}
		| LOCAL VARIABLE '=' expr		{ $$ = globalP->opr2('=', 2, globalP->id($2), $4); globalP->set_local($2); }
		| LOCAL VARIABLE				{ $$ = NULL; globalP->set_local($2); }        
		| VARIABLE '=' expr				{ $$ = globalP->opr2('=', 2, globalP->id($1), $3); globalP->assert_local($1); }
        | IF expr THEN stmt_list END	{ $$ = globalP->opr2(IF, 2, $2, $4); }        
		;
funcstat:
		FUNCTION VARIABLE '(' arguments ')' stmt_list END         { globalP->set_name($2); globalP->set_arguments($4); $$ = $6; }
		;


arguments:
		VARIABLE ',' arguments			{ $$ = globalP->opr2(',', 2, globalP->id($1), $3); }
		| VARIABLE						{ $$ = globalP->id($1); }
		;

arguments2:
		expr',' arguments2			{ $$ = globalP->opr2(',', 2, $1, $3); }
		| expr						{ $$ = $1; }
		;


exprlist:
		expr ',' exprlist			{ $$ = globalP->opr2(',', 2, $1, $3); }
		| expr						{ $$ = $1; }
		;

stmt_list:
          stmt                  { $$ = $1; }
        |  stmt stmt_list { $$ = globalP->opr2(';', 2, $1, $2); }
        ;

expr:
          YY_INTEGER               { $$ = globalP->con($1); }
        | VARIABLE '(' arguments2 ')' { $$ = globalP->function(globalP->id($1), $3); }
        | VARIABLE              { $$ = globalP->id($1); }
        | '-' expr %prec UMINUS { $$ = globalP->opr1(UMINUS, 1, $2); }
        | expr '+' expr         { $$ = globalP->opr2('+', 2, $1, $3); }
        | expr '-' expr         { $$ = globalP->opr2('-', 2, $1, $3); }
        | expr '*' expr         { $$ = globalP->opr2('*', 2, $1, $3); }
        | expr '/' expr         { $$ = globalP->opr2('/', 2, $1, $3); }
        | expr '<' expr         { $$ = globalP->opr2('<', 2, $1, $3); }
        | expr '>' expr         { $$ = globalP->opr2('>', 2, $1, $3); }
        | expr GE expr          { $$ = globalP->opr2(GE, 2, $1, $3); }
        | expr LE expr          { $$ = globalP->opr2(LE, 2, $1, $3); }
        | expr NE expr          { $$ = globalP->opr2(NE, 2, $1, $3); }
        | expr EQ expr          { $$ = globalP->opr2(EQ, 2, $1, $3); }
		| expr AND expr          { $$ = globalP->opr2(AND, 2, $1, $3); }
		| expr OR expr          { $$ = globalP->opr2(OR, 2, $1, $3); }
        | '(' expr ')'          { $$ = $2; }
		| MATHCOS '(' expr ')' { $$ = globalP->opr1(MATHCOS, 1, $3); }
		| MATHSIN '(' expr ')' { $$ = globalP->opr1(MATHSIN, 1, $3); }
		| MATHEXP '(' expr ')' { $$ = globalP->opr1(MATHEXP, 1, $3); }
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