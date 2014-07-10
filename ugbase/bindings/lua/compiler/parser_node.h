/*
 * \file	parser_node.cpp
 * \author	Martin Rupp
 *
 * Created on 20. November 2012, 10:16
 */


#ifndef PARSER_NODE_H
#define PARSER_NODE_H

typedef enum { typeCon=1, typeId=2, typeOpr=3 } nodeEnum;

/* constants */
typedef struct {
    double value;                  /* value of constant */
} conNodeType;

/* identifiers */
typedef struct {
    int i;                      /* subscript to sym array */
} idNodeType;

struct nodeType;

/* operators */
typedef struct {
    int oper;                   /* operator */
    int nops;                   /* number of operands */
    nodeType **op;	/* operands */
} oprNodeType;

struct nodeType
{

    nodeEnum type;              /* type of node */

    union {
        conNodeType con;        /* constants */
        idNodeType id;          /* identifiers */
        oprNodeType opr;        /* operators */
    };
};

#endif /*PARSER_NODE_H*/
