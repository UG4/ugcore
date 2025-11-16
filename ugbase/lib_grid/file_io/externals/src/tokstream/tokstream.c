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

/**
 * \mainpage
 *
 * The <em>tokstream</em> library is a simple, flexible and fast tokenizer
 * written in C.
 *
 * It contains only one struct, tokstream, and a number
 * of associated functions, which are all prefixed with a <tt>ts_</tt> tag.
 *
 * \section building Building
 *
 * Since the library consists of nothing more than a pair of header and
 * implementation files, the easiest solution is to directly compile it with
 * your project.
 *
 * If you, however, prefer to compile <em>tokstream</em> as a library, refer
 * to the \ref folders "cmake folder".
 *
 * \section folders Folder structure
 *
 * You are currently seeing documentation from the <em>doc</em> folder of the
 * library. The <tt>Doxyfile</tt> for building documentation is also located
 * there.
 *
 * A <tt>CMakeLists.txt</tt> file for building the library with
 * <a href="http://www.cmake.org/" target="_blank">CMake</a> can be found
 * in the <em>cmake</em> folder.
 *
 * The sources are located in the <em>tokstream</em> folder.
 *
 * \section license License
 *
 * The <em>tokstream</em> library is released under the MIT License:
 *
 * \verbatim
Copyright (c) 2008 Nicolas Tessore

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.
\endverbatim
 */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

/**
 * \file tokstream.h
 *
 * \brief The tokstream library header
 *
 * This is the header file of the tokstream library. It contains a number of
 * functions, all prefixed with <strong>ts_</strong>, to operate on files as a
 * stream of tokens.
 */
#include "../../include/tokstream/tokstream.h"


/****
 * settings
 */

/* default read buffer size */
#define TS_BUFSIZ BUFSIZ

/* number of possible characters */
#define TS_CHARMAP_SIZE 256

/* character flag map data type */
typedef char ts_charmap[TS_CHARMAP_SIZE];


/****
 * tokstream structs
 */

struct ts_state
{
    ts_charmap sep;
    ts_charmap sep2;
    ts_charmap delim;

    int eof;
    int error;

    int buf_rev;

    char* cur;
    char* tok;

    long int pos;
    int line_no;
    int char_no;

    int tok_len;
    long int tok_pos;
    int tok_line_no;
    int tok_char_no;

    char* tok_buf;
};

/**
 * \struct tokstream tokstream.h
 *
 * \brief Token stream data structure
 *
 * Data structure representing a token stream.
 *
 * Like a FILE object, a tokstream object will always be created dynamically
 * by a call to ts_open(), and deleted by the according call to ts_close().
 *
 * The structure has no publicly accessible members.
 */
struct tokstream
{
    FILE* fp;

    char* file;

    char* buf;
    int buf_size;
    int buf_len;
    int buf_rev;

    struct ts_state* state;
    struct ts_state* stack;
    int stack_size;
};


/****
 * charmap operations
 */

#define ts_charmap_clr(map) memset(map, 0, sizeof(ts_charmap))
#define ts_charmap_cpy(dst, src) memcpy(dst, src, sizeof(ts_charmap))

#define ts_charmap_0(map, c) (map[(int)c] = 0)
#define ts_charmap_1(map, c) (map[(int)c] = 1)

#define ts_charmap_get(map, c) (map[(int)c])


/****
 * inlining macros
 */

/* check if buffer is valid, if not get new, return status */
#define ts_bad_buf(ts) ((!ts->state->cur || !(*ts->state->cur)) && ts_read(ts))

#define ts_issep(ts) ts_charmap_get(ts->state->sep, *ts->state->cur)
#define ts_cissep(ts, c) ts_charmap_get(ts->state->sep, c)
#define ts_isdelim(ts) ts_charmap_get(ts->state->delim, *ts->state->cur)
#define ts_cisdelim(ts, c) ts_charmap_get(ts->state->delim, c)

/* advance cursor */
#define ts_adv_cur(ts) \
    do { \
        /* count chars and lines */ \
        if(*ts->state->cur == '\n') \
        { \
            ++ts->state->line_no; \
            ts->state->char_no = 1; \
        } \
        else \
        { \
            ++ts->state->char_no; \
        } \
         \
        /* increase cursor and position */ \
        ++ts->state->cur; \
        ++ts->state->pos; \
    } while(0) \

/* expand token */
#define ts_exp_tok(ts) \
    do { \
        /* advance cursor */ \
        ts_adv_cur(ts); \
         \
        /* increase token length */ \
        ++ts->state->tok_len; \
    } while(0) \

/* copy token to buffer */
#define ts_copy_tok(ts) \
    do { \
        /* free old token */ \
        free(ts->state->tok_buf); \
         \
        /* allocate space and copy token */ \
        ts->state->tok_buf = strncpy(malloc(ts->state->tok_len+1), ts->state->tok, ts->state->tok_len); \
         \
        /* terminate token */ \
        ts->state->tok_buf[ts->state->tok_len] = '\0'; \
    } while(0) \


/* copy string with allocation */
#define ts_strdup(str) strcpy(malloc(strlen(str)+1), str)


/****
 * internal functions declaration
 */

/* initialize a new state */
void ts_state_init(struct ts_state* state);

/* copy new state */
void ts_state_copy(struct ts_state* dst, const struct ts_state* src);

/* clean an old state */
void ts_state_clean(struct ts_state* state);

/* read new buffer for tokstream */
int ts_read(tokstream* ts);

/* normalize buffer contents */
int ts_normalize(tokstream* ts);


/****
 * implementation
 */

/**
 * \relatesalso tokstream
 *
 * \brief Open a new tokstream from file
 *
 * This opens the file given in the argument and constructs a tokstream around
 * it.
 *
 * Before reading can begin, you have to set the separators and delimiters for
 * the tokstream by using ts_sep(), ts_delim() and the according on and
 * off functions.
 *
 * \sa ts_sep(), ts_sep_on(), ts_sep_off()
 * \sa ts_delim(), ts_delim_on(), ts_delim_off()
 *
 * \param file The filename to be opened.
 *
 * \returns Returns a pointer to the opened tokstream, or nullptr on error.
 */
tokstream* ts_open(const char* file)
{
    FILE* fp;
    tokstream* ts;

    /* open file and keep it open */
    fp = fopen(file, "rb");
    if(!fp)
        return nullptr;

    /* allocate new tokstream */
    ts = malloc(sizeof(tokstream));

    /* initialize */

    /* set file pointer */
    ts->fp = fp;

    /* copy filename */
    ts->file = ts_strdup(file);

    /* create buffer */
    ts->buf_size = TS_BUFSIZ;
    ts->buf = malloc(ts->buf_size);
    ts->buf_len = 0;

    /* start with 0 buffer revision */
    ts->buf_rev = 0;

    /* initialize state stack */
    ts->stack = malloc(sizeof(struct ts_state));
    ts->state = ts->stack;
    ts->stack_size = 1;

    /* initialize main state */
    ts_state_init(ts->state);

    /* set file status */
    ts->state->eof = feof(ts->fp);
    ts->state->error = ferror(ts->fp);

    return ts;
}

/**
 * \relatesalso tokstream
 *
 * \brief Close a tokstream
 *
 * Closes the file of the tokstream and frees all allocated memory.
 *
 * \param ts The tokstream to be closed.
 */
void ts_close(tokstream* ts)
{
    /* close file */
    fclose(ts->fp);

    /* clean all stacked states */
    while(ts->state >= ts->stack)
    {
        ts_state_clean(ts->state);
        --ts->state;
    }

    /* free memory */
    free(ts->stack);
    free(ts->buf);
    free(ts->file);
    free(ts);
}

/**
 * \relatesalso tokstream
 *
 * \brief Push the current tokstream state onto stack
 *
 * This saves all current information of the tokstream, but does not alter any
 * of its properties. Reading will continue just like before the call.
 *
 * Using ts_pop(), reading can later be continued exactly from the moment
 * ts_push() was called. Every setting in the tokstream is restored.
 *
 * It is possible to push multiple states on top of each other onto the stack.
 * The limit of this is defined by available memory.
 *
 * If pushing the state onto the stack fails (due to memory allocation), all of
 * the stack is lost. If possible, a single empty state is pushed onto the
 * stack, so that ts_ functions will not fail mysteriously.
 *
 * \sa ts_pop()
 *
 * \param ts The tokstream of which the state is to be saved to the stack.
 *
 * \returns In case of an error, a non-zero value is returned.
 */
int ts_push(tokstream* ts)
{
    /* no current state */
    ts->state = nullptr;

    /* resize stack */
    ts->stack = realloc(ts->stack, (ts->stack_size + 1) * sizeof(struct ts_state));

    /* check there was enough memory available */
    if(!ts->stack)
    {
        /* create fallout state */
        ts->stack = malloc(sizeof(struct ts_state));
        ts->state = ts->stack;

        /* clean fallout state if possible */
        if(ts->state)
            ts_state_init(ts->state);

        /* return error */
        return 1;
    }

    /* duplicate previous state */
    ts_state_copy(ts->stack + ts->stack_size + 1, ts->stack + ts->stack_size);

    /* increment stack counter */
    ++ts->stack_size;

    /* set current state */
    ts->state = ts->stack + ts->stack_size;

    /* success */
    return 0;
}

/**
 * \relatesalso tokstream
 *
 * \brief Pop the current state from the stack
 *
 * This returns the tokstream to a state which was previously pushed to the
 * stack with ts_push(). Reading will continue as if the calling of ts_push()
 * and subsequently ts_pop() never occurred.
 *
 * \sa ts_push()
 *
 * \param ts The tokstream of which a stacked state is to be restored.
 *
 * \returns On attempting to pop the last state from the stack, the function
 *          returns a non-zero value.
 */
int ts_pop(tokstream* ts)
{
    /* prevent stack underflow */
    if(ts->state == ts->stack)
        return 1;

    /* clean discarded state */
    ts_state_clean(ts->state);

    /* resize state stack */
    ts->stack = realloc(ts->stack, (ts->stack_size - 1) * sizeof(struct ts_state));

    /* decrement stack counter */
    --ts->stack_size;

    /* set current state */
    ts->state = ts->stack + ts->stack_size;

    /* check if buffer changed */
    if(ts->buf_rev > ts->state->buf_rev)
    {
        /* invalidate cursor and token */
        ts->state->cur = nullptr;
        ts->state->tok = nullptr;
    }

    /* success */
    return 0;
}

/**
 * \name Stream information
 *
 * Functions to get status information about a token stream.
 *
 * \{
 */

/**
 * \relatesalso tokstream
 *
 * \brief Check if a tokstream is at EOF
 *
 * Checks if the tokstream has reached end-of-file (EOF). Calls feof() for the
 * underlying file object.
 *
 * \param ts The tokstream to check for EOF.
 *
 * \returns Returns a non-zero value if file is at EOF.
 */
int ts_eof(const tokstream* ts)
{
    /* delay eof to end of buffer */
    return ts->state->eof && (!ts->state->cur || !(*ts->state->cur));
}

/**
 * \relatesalso tokstream
 *
 * \brief Get file error flag of a tokstream
 *
 * This function returns the value of ferror() from <em>after the last read
 * operation</em>. Because of buffering, this is asynchronuous with ts_get()
 * calls.
 *
 * \param ts The tokstream of which to get the file error flag.
 *
 * \returns Report the last value of ferror() for the tokstream's FILE* object.
 */
int ts_error(const tokstream* ts)
{
    /* delay errors to end of buffer */
    return ts->state->error && (!ts->state->cur || !(*ts->state->cur));
}

/**
 * \relatesalso tokstream
 *
 * \brief Return current line number
 *
 * This function gives the number of the line at which the current token was
 * read in the file. If there is no current token, the line number of the next
 * processed character is returned.
 *
 * \param ts The stream to return its line number.
 *
 * \returns Returns the line number of the current token.
 */
int ts_line(const tokstream* ts)
{
    if(!ts->state->tok)
        return ts->state->line_no;

    return ts->state->tok_line_no;
}

/**
 * \relatesalso tokstream
 *
 * \brief Return current character position
 *
 * This function gives the position of the character at which the current token
 * was read in its line. If there is no current token, the character position
 * of the next processed character is returned.
 *
 * \param ts The stream to return its character position.
 *
 * \returns Returns the character position of the current token.
 */
int ts_char(const tokstream* ts)
{
    if(!ts->state->tok)
        return ts->state->char_no;

    return ts->state->tok_char_no;
}

/**
 * \relatesalso tokstream
 *
 * \brief Return the current token
 *
 * This function returns the exact same token that was fetched and returned
 * by the last call to ts_get(). If no such call was ever made, it will return
 * nullptr.
 *
 * The string returned belongs to the tokstream object. It will be invalid
 * on the next call to ts_get().
 *
 * \sa ts_get()
 *
 * \param ts The stream to get the current token of.
 *
 * \returns A string containing the current token, or nullptr if no such token
 *          exists.
 */
const char* ts_tok(const tokstream* ts)
{
    /* check if token is valid */
    if(!ts->state->tok)
        return nullptr;

    /* return buffered token */
    return ts->state->tok_buf;
}

/**
 * \}
 */

/**
 * \name Token getters
 *
 * Functions to get tokens from the input stream, modifying the current token.
 *
 * \{
 */

/**
 * \relatesalso tokstream
 *
 * \brief Get the next token from stream
 *
 * Search the input stream for the next token, according to current separator
 * and delimiter settings.
 *
 * The string returned belongs to the tokstream object. It will be invalid
 * on the next call to ts_get().
 *
 * \note This fetches the <em>next</em> token from the stream. If you want to
 *       get the <em>current</em> token, use ts_tok().
 *
 * \sa ts_tok()
 *
 * \param ts The stream to get the token from.
 *
 * \returns A zero-terminated string containing the next token is returned, or
 *          nullptr in case of an error (ie. EOF occurred while getting token).
 */
const char* ts_get(tokstream* ts)
{
    /* check if buffer is good */
    if(ts_bad_buf(ts))
        return nullptr;

    /* seek beginning of token */
    while(ts_issep(ts))
    {
        /* advance cursor */
        ts_adv_cur(ts);

        /* check if buffer is still good */
        if(ts_bad_buf(ts))
            return nullptr;
    }

    /* tokenize string beginning from cursor */
    ts->state->tok = ts->state->cur;

    /* reset token length */
    ts->state->tok_len = 0;

    /* store position of token */
    ts->state->tok_pos = ts->state->pos;
    ts->state->tok_line_no = ts->state->line_no;
    ts->state->tok_char_no = ts->state->char_no;

    /* expand token */
    ts_exp_tok(ts);

    /* check if token is not a delimiter */
    if(!ts_cisdelim(ts, *ts->state->tok))
    {
        /* move cursor forward until separator or delimiter */
        while(!ts_issep(ts) && !ts_isdelim(ts))
        {
            /* expand token */
            ts_exp_tok(ts);

            /* buffer ends here, and so does token */
            if(!(*ts->state->cur))
                break;
        }
    }

    /* copy token to token buffer */
    ts_copy_tok(ts);

    /* return the token found, from buffer */
    return ts->state->tok_buf;
}

/**
 * \relatesalso tokstream
 *
 * \brief Unget current token
 *
 * Put the current token back into the stream, so that the next call to ts_get()
 * might again return it.
 *
 * This is useful when you need to changing separators and delimiters, as well
 * as when you need to peek at the next token.
 *
 * \note This might cause a buffer refresh.
 *
 * \param ts The stream to unget the token to.
 *
 * \returns If the previous token is no longer available and cannot be ungot,
 *          a non-zero value is returned.
 */
int ts_unget(tokstream* ts)
{
    /* check if there is a token in buffer */
    if(!ts->state->tok)
        return 1;

    /* set cursor to token */
    ts->state->cur = ts->state->tok;
    ts->state->pos = ts->state->tok_pos;
    ts->state->line_no = ts->state->tok_line_no;
    ts->state->char_no = ts->state->tok_char_no;

    /* no current token */
    free(ts->state->tok_buf);
    ts->state->tok_buf = nullptr;
    ts->stack->tok = nullptr;

    /* success */
    return 0;
}

/**
 * \relatesalso tokstream
 *
 * \brief Get rest of line from stream
 *
 * Stores the rest of the current line (everything until newline character)
 * as the current token and returns it.
 *
 * The string will begin with the first non-separator character. If there is no
 * non-separator character until the end of line, the returned string will be
 * empty.
 *
 * The newline will be consumed, the stream will be positioned at the beginning
 * of the next line.
 *
 * Subsequent calls to ts_get() will return the line as returned by this
 * function.
 *
 * \param ts The stream to get the line from.
 *
 * \returns Returns a string containing the line, or nullptr if an error occurred.
 */
const char* ts_getline(tokstream* ts)
{
    /* tokenize until newline */
    if(!ts_seekc(ts, '\n'))
        return nullptr;

    /* advance cursor past newline */
    ts_adv_cur(ts);

    /* return the line from token buffer */
    return ts->state->tok_buf;
}

/**
 * \}
 */

/**
 * \name Input skipping
 *
 * Functions to skip over parts of the input without modifying the current
 * stream status, ie. current token.
 *
 * \{
 */

/**
 * \relatesalso tokstream
 *
 * \brief Skip over the next token
 *
 * Using this function, the next token can be skipped without invalidating the
 * current token. This might be useful if the next token is already known, ie.
 * from a call to a seek function.
 *
 * \param ts The stream in which to skip a token.
 *
 * \returns Returns a non-zero value if an error occurred.
 */
int ts_skip(tokstream* ts)
{
    /* check if buffer is good */
    if(ts_bad_buf(ts))
        return 1;

    /* seek beginning of token */
    while(ts_issep(ts))
    {
        /* advance cursor */
        ts_adv_cur(ts);

        /* check if buffer is still good */
        if(ts_bad_buf(ts))
            return 1;
    }

    /* advance cursor */
    ts_adv_cur(ts);

    /* check if token is not a delimiter */
    if(!ts_cisdelim(ts, *ts->state->tok))
    {
        /* move cursor forward until separator or delimiter */
        while(!ts_issep(ts) && !ts_isdelim(ts))
        {
            /* advance cursor */
            ts_adv_cur(ts);

            /* buffer ends here, and so does token */
            if(!(*ts->state->cur))
                break;
        }
    }

    /* done */
    return 0;
}

/**
 * \relatesalso tokstream
 *
 * \brief Skip line in stream
 *
 * Discards the current line in the stream and sets the stream position to the
 * beginning of the next line.
 *
 * \note Invalidates current token.
 *
 * \param ts The stream in which to skip a line.
 *
 * \returns On error, a non-zero value is returned.
 */
int ts_skipline(tokstream* ts)
{
    /* invalidate token */
    ts->state->tok = nullptr;

    /* check if buffer is good */
    if(ts_bad_buf(ts))
        return 1;

    /* increment cursor until we find newline */
    while(*ts->state->cur != '\n')
    {
        /* advance cursor */
        ts_adv_cur(ts);

        /* make sure buffer is still filled */
        if(ts_bad_buf(ts))
            return 1;
    }

    /* advance past newline */
    ts_adv_cur(ts);

    /* success */
    return 0;
}

/**
 * \name Input seeking
 *
 * Seek to specific position in token stream.
 *
 * \{
 */

/**
 * \relatesalso tokstream
 *
 * \brief Seek to token
 *
 * The searched token will be the <em>current</em> token. The next call to
 * ts_get() will fetch a new token.
 *
 * \param ts The token stream to operate on.
 * \param tok The token to seek.
 *
 * \returns A non-zero value is returned to indicate the token was not found.
 */

int ts_seek(tokstream* ts, const char* tok)
{
	/* get tokens from ts until tok is found */
	do
	{
		/* check if current token is right */
		if(strcmp(ts->state->tok, tok) == 0)
			return 0;
	}
	while(ts_get(ts) != nullptr);

	/* token was not found */
	return 1;
}

/**
 * \relatesalso tokstream
 *
 * \brief Seek to character
 *
 * Stores the input until it encounters the \a c character and returns it as a
 * token.
 *
 * The string will begin with the first non-separator character. If there is no
 * non-separator character until the character \a c is found, the returned
 * string will be empty.
 *
 * The character \a c will not be consumed, it can be part of the next token.
 *
 * Subsequent calls to ts_get() will return the token as returned by this
 * function.
 *
 * \param ts The stream to get the token from.
 * \param c The character that ends the token.
 *
 * \returns Returns a string containing the token, or nullptr if an error occurred.
 */
const char* ts_seekc(tokstream* ts, char c)
{
    /* check if buffer is good */
    if(ts_bad_buf(ts))
        return nullptr;

    /* seek beginning of token */
    while(ts_issep(ts) && *ts->state->cur != c)
    {
        /* advance cursor */
        ts_adv_cur(ts);

        /* check if buffer is still good */
        if(ts_bad_buf(ts))
            return nullptr;
    }

    /* tokenize string beginning from cursor */
    ts->state->tok = ts->state->cur;

    /* reset token length */
    ts->state->tok_len = 0;

    /* store position of token */
    ts->state->tok_pos = ts->state->pos;
    ts->state->tok_line_no = ts->state->line_no;
    ts->state->tok_char_no = ts->state->char_no;

    /* move cursor forward until char is found */
    while(*ts->state->cur != c)
    {
        /* expand token */
        ts_exp_tok(ts);

        /* buffer ends here, and so does token */
        if(!(*ts->state->cur))
            break;
    }

    /* copy token to token buffer */
    ts_copy_tok(ts);

    /* return the token buffer */
    return ts->state->tok_buf;
}

/**
 * \relatesalso tokstream
 *
 * \brief Seek to any character from array
 *
 * Stores the input until it encounters any of the \a ca characters and return
 * it as a token.
 *
 * The string will begin with the first non-separator character. If there is no
 * non-separator character until a character in \a ca is found, the returned
 * string will be empty.
 *
 * The character from \a ca will not be consumed, it can be part of the next
 * token.
 *
 * Subsequent calls to ts_get() will return the token as returned by this
 * function.
 *
 * \param ts The stream to get the token from.
 * \param ca The characters that end the token.
 *
 * \returns Returns a string containing the token, or nullptr if an error occurred.
 */
const char* ts_seekca(tokstream* ts, const char* ca)
{
    /* check if buffer is good */
    if(ts_bad_buf(ts))
        return nullptr;

    /* seek beginning of token */
    while(ts_issep(ts) && !strchr(ca, *ts->state->cur))
    {
        /* advance cursor */
        ts_adv_cur(ts);

        /* check if buffer is still good */
        if(ts_bad_buf(ts))
            return nullptr;
    }

    /* tokenize string beginning from cursor */
    ts->state->tok = ts->state->cur;

    /* reset token length */
    ts->state->tok_len = 0;

    /* store position of token */
    ts->state->tok_pos = ts->state->pos;
    ts->state->tok_line_no = ts->state->line_no;
    ts->state->tok_char_no = ts->state->char_no;

    /* move cursor forward until char is found */
    while(!strchr(ca, *ts->state->cur))
    {
        /* expand token */
        ts_exp_tok(ts);

        /* buffer ends here, and so does token */
        if(!(*ts->state->cur))
            break;
    }

    /* copy token to token buffer */
    ts_copy_tok(ts);

    /* return the token buffer */
    return ts->state->tok_buf;
}

/**
 * \}
 */

/**
 * \name Separator and delimiter control
 *
 * Functions to set which characters act as separators and which act as
 * delimiters.
 *
 * \{
 */

/**
 * \relatesalso tokstream
 *
 * \brief Set separator characters
 */
void ts_sep(tokstream* ts, const char* sep)
{
    /* turn all separator flags off */
    ts_charmap_clr(ts->state->sep);

    /* set separators */
    for(; *sep; ++sep)
        ts_charmap_1(ts->state->sep, *sep);

    /* make backup of separator flags */
    ts_charmap_cpy(ts->state->sep2, ts->state->sep);

    /* renormalize buffer */
    ts_normalize(ts);
}

/**
 * \relatesalso tokstream
 *
 * \brief Set character as separator
 */
void ts_sep_on(tokstream* ts, char c)
{
    /* set separator */
    ts_charmap_1(ts->state->sep, c);
    ts_charmap_1(ts->state->sep2, c);

    /* renormalize buffer */
    ts_normalize(ts);
}

/**
 * \relatesalso tokstream
 *
 * \brief Unset character as separator
 */
void ts_sep_off(tokstream* ts, char c)
{
    /* unset separator */
    ts_charmap_0(ts->state->sep, c);
    ts_charmap_0(ts->state->sep2, c);

    /* renormalize buffer */
    ts_normalize(ts);
}

/**
 * \relatesalso tokstream
 *
 * \brief Set delimiter characters
 *
 * */
void ts_delim(tokstream* ts, const char* delim)
{
    /* turn all delimiter flags off */
    ts_charmap_clr(ts->state->delim);

    /* restore all separator flags */
    ts_charmap_cpy(ts->state->sep, ts->state->sep2);

    /* set delimiters */
    for(; *delim; ++delim)
    {
        /* remove sep flag */
        ts_charmap_0(ts->state->sep, *delim);

        /* set delim flag */
        ts_charmap_1(ts->state->delim, *delim);
    }

    /* renormalize buffer */
    ts_normalize(ts);
}

/**
 * \relatesalso tokstream
 *
 * \brief Set character as delimiter
 */
void ts_delim_on(tokstream* ts, char c)
{
    /* remove sep flag */
    ts_charmap_0(ts->state->sep, c);

    /* set delimiter */
    ts_charmap_1(ts->state->delim, c);

    /* renormalize buffer */
    ts_normalize(ts);
}

/**
 * \relatesalso tokstream
 *
 * \brief Unset character as delimiter
 */
void ts_delim_off(tokstream* ts, char c)
{
    /* unset delimiter */
    ts_charmap_0(ts->state->delim, c);

    /* restore sep flag */
    if(ts_charmap_get(ts->state->sep2, c))
        ts_charmap_1(ts->state->sep, c);

    /* renormalize buffer */
    ts_normalize(ts);
}

/**
 * \}
 */

/**
 * \relatesalso tokstream
 *
 * \brief Set input buffer size for stream
 */
int ts_bufsiz(tokstream* ts, int size)
{
    /* invalidate buffer for all states */
    ts->state->cur = nullptr;
    ts->state->tok = nullptr;
    ++ts->buf_rev;

    /* reallocate buffer */
    ts->buf = realloc(ts->buf, size);

    /* check if realloc failed */
    if(!ts->buf)
    {
        /* allocate old buffer size */
        ts->buf = malloc(ts->buf_size);

        /* error */
        return 1;
    }

    /* set size of buffer */
    ts->buf_size = size;

    /* success */
    return 0;
}


/****
 * internal functions
 */

void ts_state_init(struct ts_state* state)
{
    ts_charmap_clr(state->sep);
    ts_charmap_clr(state->sep2);
    ts_charmap_clr(state->delim);

    state->eof = 0;
    state->error = 0;

    state->buf_rev = 0;

    state->cur = nullptr;
    state->tok = nullptr;

    state->pos = 0;
    state->line_no = 1;
    state->char_no = 1;

    state->tok_len = 0;
    state->tok_pos = 0;
    state->tok_line_no = 1;
    state->tok_char_no = 1;

    state->tok_buf = nullptr;
}

void ts_state_copy(struct ts_state* dst, const struct ts_state* src)
{
    ts_charmap_cpy(dst->sep, src->sep);
    ts_charmap_cpy(dst->sep2, src->sep2);
    ts_charmap_cpy(dst->delim, src->delim);

    dst->eof = src->eof;
    dst->error = src->error;

    dst->buf_rev = src->buf_rev;

    dst->cur = src->cur;
    dst->tok = src->tok;

    dst->pos = src->pos;
    dst->line_no = src->line_no;
    dst->char_no = src->char_no;

    dst->tok_len = src->tok_len;
    dst->tok_pos = src->tok_pos;
    dst->tok_line_no = src->tok_line_no;
    dst->tok_char_no = src->tok_char_no;

    dst->tok_buf = src->tok_buf ? ts_strdup(src->tok_buf) : nullptr;
}

void ts_state_clean(struct ts_state* state)
{
    free(state->tok_buf);
}

int ts_read(tokstream* ts)
{
    int seek_err;

    /* check if file is at eof already */
    if(ts->state->eof)
        return 1;

    /* seek to file position */
    seek_err = fseek(ts->fp, ts->state->pos, SEEK_SET);

    /* update error and eof data */
    ts->state->eof = feof(ts->fp);
    ts->state->error = ferror(ts->fp);

    /* if there was a seek error, rather not give fp */
    if(seek_err)
        return 1;

    /* invalidate cursor and token */
    ts->state->cur = nullptr;
    ts->state->tok = nullptr;

    /* increase buffer revision */
    ++ts->buf_rev;
    ts->state->buf_rev = ts->buf_rev;

    /* get BUFSIZ chars from file to buffer */
    ts->buf_len = fread(ts->buf, 1, ts->buf_size-1, ts->fp);

    /* terminate buffer string */
    ts->buf[ts->buf_len] = '\0';

    /* set error indicators */
    ts->state->eof = feof(ts->fp);
    ts->state->error = ferror(ts->fp);

    /* break on error before updating tokstream */
    if(ts->state->error)
        return 1;

    /* set cursor to beginning of buffer */
    ts->state->cur = ts->buf;

    /* normalize tokstream */
    ts_normalize(ts);

    /* success */
    return 0;
}

int ts_normalize(tokstream* ts)
{
    int trim;

    /* test buffer */
    if(ts->buf_rev == 0)
        return 0;

    /* don't trim when at eof */
    trim = 0;
    if(!ts->state->eof)
    {
        /* trim token chars from end of buffer */
        char* back = ts->buf + ts->buf_len - 1;
        while(ts->buf_len > 0)
        {
            /* trim until separator or delimiter */
            if(ts_cissep(ts, *back) || ts_cisdelim(ts, *back))
                break;

            /* trim buffer */
            --back;
            --ts->buf_len;
        }

        /* check whether buffer was trimmed */
        if(*(back+1))
        {
            /* buffer was trimmed */
            trim = 1;

            /* terminate trimmed buffer */
            *(back+1) = '\0';
        }
    }

    if(trim)
    {
        /* buffer changed, update buffer revision */
        ++ts->buf_rev;
        ++ts->state->buf_rev;
    }

    /* report changes to buffer */
    return trim;
}
