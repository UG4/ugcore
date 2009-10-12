/****
 * Copyright (c) 2008 Nicolas Tessore
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 * THE SOFTWARE.
 ****/

#ifndef TOKSTREAM_H_
#define TOKSTREAM_H_


/* Token stream data structure */
typedef struct tokstream tokstream;


/* Open a new tokstream from file */
tokstream* ts_open(const char* file);

/* Close a tokstream */
void ts_close(tokstream* ts);


/* Push the current tokstream state onto stack */
int ts_push(tokstream* ts);

/* Pop the current state from the stack */
int ts_pop(tokstream* ts);


/* Check if a tokstream is at EOF */
int ts_eof(const tokstream* ts);

/* Get file error flag of a tokstream */
int ts_error(const tokstream* ts);

/* Return current line number */
int ts_line(const tokstream* ts);

/* Return current character position */
int ts_char(const tokstream* ts);

/* Return the current token */
const char* ts_tok(const tokstream* ts);


/* Get the next token from stream */
const char* ts_get(tokstream* ts);

/* Unget current token */
int ts_unget(tokstream* ts);

/* Get rest of line from stream */
const char* ts_getline(tokstream* ts);


/* Skip over the next token */
int ts_skip(tokstream* ts);

/* Skip line in stream */
int ts_skipline(tokstream* ts);


/* Seek to token */
int ts_seek(tokstream* ts, const char* tok);

/* Seek to character */
const char* ts_seekc(tokstream* ts, char c);

/* Seek to any character from array */
const char* ts_seekca(tokstream* ts, const char* ca);


/* Set separator characters */
void ts_sep(tokstream* ts, const char* sep);

/* Set character as separator */
void ts_sep_on(tokstream* ts, char c);

/* Unset character as separator */
void ts_sep_off(tokstream* ts, char c);

/* Set delimiter characters */
void ts_delim(tokstream* ts, const char* delim);

/* Set character as delimiter */
void ts_delim_on(tokstream* ts, char c);

/* Unset character as delimiter */
void ts_delim_off(tokstream* ts, char c);


/* Set input buffer size for stream */
int ts_bufsiz(tokstream* ts, int size);

#endif
