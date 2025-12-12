/* linenoise.c -- guerrilla line editing library against the idea that a
* line editing lib needs to be 20,000 lines of C code.
*
* You can find the latest source code at:
*
* http://github.com/antirez/linenoise
*
* Does a number of crazy assumptions that happen to be true in 99.9999% of
* the 2010 UNIX computers around.
*
* Copyright (c) 2010, Salvatore Sanfilippo <antirez at gmail dot com>
* All rights reserved.
*
* Redistribution and use in source and binary forms, with or without
* modification, are permitted provided that the following conditions are met:
*
* * Redistributions of source code must retain the above copyright notice,
* this list of conditions and the following disclaimer.
* * Redistributions in binary form must reproduce the above copyright
* notice, this list of conditions and the following disclaimer in the
* documentation and/or other materials provided with the distribution.
* * Neither the name of Redis nor the names of its contributors may be used
* to endorse or promote products derived from this software without
* specific prior written permission.
*
* THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
* AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
* IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
* ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
* LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
* CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
* SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
* INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
* CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
* ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
* POSSIBILITY OF SUCH DAMAGE.
*
* References:
* - http://invisible-island.net/xterm/ctlseqs/ctlseqs.html
* - http://www.3waylabs.com/nw/WWW/products/wizcon/vt220.html
*
* Todo list:
* - Switch to gets() if $TERM is something we can't support.
* - Filter bogus Ctrl+<char> combinations.
* - Win32 support
*
* Bloat:
* - Completion?
* - History search like Ctrl+r in readline?
*
* List of escape sequences used by this program, we do everything just
* with three sequences. In order to be so cheap we may have some
* flickering effect with some slow terminal, but the lesser sequences
* the more compatible.
*
* CHA (Cursor Horizontal Absolute)
* Sequence: ESC [ n G
* Effect: moves cursor to column n
*
* EL (Erase Line)
* Sequence: ESC [ n K
* Effect: if n is 0 or missing, clear from cursor to end of line
* Effect: if n is 1, clear from beginning of line to cursor
* Effect: if n is 2, clear entire line
*
* CUF (CUrsor Forward)
* Sequence: ESC [ n C
* Effect: moves cursor forward of n chars
*
*/

/* I changed some variable names (new -> newChar) and added two casts in order
 * fix issues with c++ compilers (void* -> char**).
 * Sebastian Reiter (s.b.reiter@googlemail.com)
 */

/*
 * Added support for word completion (see linenoiseSetCompletionFunction)
 * Martin Rupp (m.rupp@gcsc.uni-frankfurt.de)
 */

#include <termios.h>
#include <unistd.h>
#include <cstdlib>
#include <cstdio>
#include <cerrno>
#include <cstring>
#include <sys/ioctl.h>

#include "linenoise.h"

#define LINENOISE_DEFAULT_HISTORY_MAX_LEN 100
#define LINENOISE_MAX_LINE 4096
static const char *unsupported_term[] = {"dumb","cons25",nullptr};

static termios orig_termios; /* in order to restore at exit */
static int rawmode = 0; /* for atexit() function to check if restore is needed*/
static int atexit_registered = 0; /* register atexit just 1 time */
static int history_max_len = LINENOISE_DEFAULT_HISTORY_MAX_LEN;
static int history_len = 0;
char **history = nullptr;

static CompletionFunctionPtr complete = nullptr;

// linenoiseSetCompletionFunction
/**
 * \brief Set a function to implement word completion
 * The function "complete" is only called when the cursor is at the last sign
 * of the buffer (pos = strlen(pos) and the user hits <tab>.
 * The complete function is of type
 * typedef int (*CompletionFunctionPtr)(char *buf, int len, int buflen, int iPrintCompletionList);
 * where buf is the current line buffer, line its length and buflen its maximal length
 * complete has to return the strlen of buf after modifications on it.
 * iPrintCompletionList is #(of times tab is hit)-1.
 * The line then gets redrawn with pos = strlen(buf).
 */
int linenoiseSetCompletionFunction(CompletionFunctionPtr new_complete)
{
	complete = new_complete;
	return 0;
}

void linenoiseAtExit();
int linenoiseHistoryAdd(const char *line);

static int isUnsupportedTerm() {
    char *term = getenv("TERM");
    int j;

    if (term == nullptr) return 0;
    for (j = 0; unsupported_term[j]; j++)
        if (!strcasecmp(term,unsupported_term[j])) return 1;
    return 0;
}



static void freeHistory() {
    if (history) {
        int j;

        for (j = 0; j < history_len; j++)
            free(history[j]);
        free(history);
    }
}

static int enableRawMode(int fd) {
    termios raw;

    if (!isatty(STDIN_FILENO)) goto fatal;
    if (!atexit_registered) {
        atexit(linenoiseAtExit);
        atexit_registered = 1;
    }
    if (tcgetattr(fd,&orig_termios) == -1) goto fatal;

    raw = orig_termios; /* modify the original mode */
    /* input modes: no break, no CR to NL, no parity check, no strip char,
* no start/stop output control. */
    raw.c_iflag &= ~(BRKINT | ICRNL | INPCK | ISTRIP | IXON);
    /* output modes - disable post processing */
    raw.c_oflag &= ~(OPOST);
    /* control modes - set 8 bit chars */
    raw.c_cflag |= (CS8);
    /* local modes - choing off, canonical off, no extended functions,
* no signal chars (^Z,^C) */
    raw.c_lflag &= ~(ECHO | ICANON | IEXTEN | ISIG);
    /* control chars - set return condition: min number of bytes and timer.
* We want read to return every single byte, without timeout. */
    raw.c_cc[VMIN] = 1; raw.c_cc[VTIME] = 0; /* 1 byte, no timer */

    /* put terminal in raw mode after flushing */
    if (tcsetattr(fd,TCSAFLUSH,&raw) < 0) goto fatal;
    rawmode = 1;
    return 0;

fatal:
    errno = ENOTTY;
    return -1;
}

static void disableRawMode(int fd) {
    /* Don't even check the return value as it's too late. */
    if (rawmode && tcsetattr(fd,TCSAFLUSH,&orig_termios) != -1)
        rawmode = 0;
}

/* At exit we'll try to fix the terminal to the initial conditions. */
void linenoiseAtExit() {
    disableRawMode(STDIN_FILENO);
    freeHistory();
}

static int getColumns() {
    winsize ws;

    if (ioctl(1, TIOCGWINSZ, &ws) == -1) return 80;
    return ws.ws_col;
}

static void refreshLine(int fd, const char *prompt, char *buf, size_t len, size_t pos, size_t cols) {
    char seq[64];
    size_t plen = strlen(prompt);
    
    while((plen+pos) >= cols) {
        buf++;
        len--;
        pos--;
    }
    while (plen+len > cols) {
        len--;
    }

    /* Cursor to left edge */
    snprintf(seq,64,"\x1b[0G");
    if (write(fd,seq,strlen(seq)) == -1) return;
    /* Write the prompt and the current buffer content */
    if (write(fd,prompt,strlen(prompt)) == -1) return;
    if (write(fd,buf,len) == -1) return;
    /* Erase to right */
    snprintf(seq,64,"\x1b[0K");
    if (write(fd,seq,strlen(seq)) == -1) return;
    /* Move cursor to original position. */
    snprintf(seq,64,"\x1b[0G\x1b[%dC", static_cast<int>(pos + plen));
    if (write(fd,seq,strlen(seq)) == -1) return;
}

static int linenoisePrompt(int fd, char *buf, size_t buflen, const char *prompt) {
    size_t plen = strlen(prompt);
    size_t pos = 0;
    size_t len = 0;
    size_t cols = getColumns();
    int history_index = 0;

    buf[0] = '\0';
    buflen--; /* Make sure there is always space for the nulterm */

    /* The latest history entry is always our current buffer, that
* initially is just an empty string. */
    linenoiseHistoryAdd("");
    
    if (write(fd,prompt,plen) == -1) return -1;
    while(true)
    {
        char c;
        char seq[2], seq2[2];

        int nread = read(fd, &c, 1);
        if (nread <= 0) return len;

        static int tabhitcount = 0;
        if(c == '\t' && pos == len && complete != nullptr)
        {
        	tabhitcount++;
        	disableRawMode(fd);
        	pos = len = complete(buf, len, buflen, tabhitcount-1);
        	enableRawMode(fd);
        	refreshLine(fd,prompt,buf,len,pos,cols);
        }
        else
        	tabhitcount = 0;

        switch(c) {
        case '\t': break;
        case 13: /* enter */
        case 4: /* ctrl-d */
            history_len--;
            free(history[history_len]);
            return (len == 0 && c == 4) ? -1 : static_cast<int>(len);
        case 3: /* ctrl-c */
            errno = EAGAIN;
            return -1;
        case 127: /* backspace */
        case 8: /* ctrl-h */
            if (pos > 0 && len > 0) {
                memmove(buf+pos-1,buf+pos,len-pos);
                pos--;
                len--;
                buf[len] = '\0';
                refreshLine(fd,prompt,buf,len,pos,cols);
            }
            break;
        case 20: /* ctrl-t */
            if (pos > 0 && pos < len) {
                int aux = buf[pos-1];
                buf[pos-1] = buf[pos];
                buf[pos] = aux;
                if (pos != len-1) pos++;
                refreshLine(fd,prompt,buf,len,pos,cols);
            }
            break;
        case 2: /* ctrl-b */
            goto left_arrow;
        case 6: /* ctrl-f */
            goto right_arrow;
        case 16: /* ctrl-p */
            seq[1] = 65;
            goto up_down_arrow;
        case 14: /* ctrl-n */
            seq[1] = 66;
            goto up_down_arrow;
            break;
        case 27: /* escape sequence */
            if (read(fd,seq,2) == -1) break;
            if (seq[0] == 91 && seq[1] == 68) {
left_arrow:
                /* left arrow */
                if (pos > 0) {
                    pos--;
                    refreshLine(fd,prompt,buf,len,pos,cols);
                }
            } else if (seq[0] == 91 && seq[1] == 67) {
right_arrow:
                /* right arrow */
                if (pos != len) {
                    pos++;
                    refreshLine(fd,prompt,buf,len,pos,cols);
                }
            } else if (seq[0] == 91 && (seq[1] == 65 || seq[1] == 66)) {
up_down_arrow:
                /* up and down arrow: history */
                if (history_len > 1) {
                    /* Update the current history entry before to
* overwrite it with tne next one. */
                    free(history[history_len-1-history_index]);
                    history[history_len-1-history_index] = strdup(buf);
                    /* Show the new entry */
                    history_index += (seq[1] == 65) ? 1 : -1;
                    if (history_index < 0) {
                        history_index = 0;
                        break;
                    } else if (history_index >= history_len) {
                        history_index = history_len-1;
                        break;
                    }
                    strncpy(buf,history[history_len-1-history_index],buflen);
                    buf[buflen] = '\0';
                    len = pos = strlen(buf);
                    refreshLine(fd,prompt,buf,len,pos,cols);
                }
            } else if (seq[0] == 91 && seq[1] > 48 && seq[1] < 55) {
                /* extended escape */
                if (read(fd,seq2,2) == -1) break;
                if (seq[1] == 51 && seq2[0] == 126) {
                    /* delete */
                    if (len > 0 && pos < len) {
                        memmove(buf+pos,buf+pos+1,len-pos-1);
                        len--;
                        buf[len] = '\0';
                        refreshLine(fd,prompt,buf,len,pos,cols);
                    }
                }
            }
            break;

        default:
            if (len < buflen) {
                if (len == pos) {
                    buf[pos] = c;
                    pos++;
                    len++;
                    buf[len] = '\0';
                    if (plen+len < cols) {
                        /* Avoid a full update of the line in the
* trivial case. */
                        if (write(fd,&c,1) == -1) return -1;
                    } else {
                        refreshLine(fd,prompt,buf,len,pos,cols);
                    }
                } else {
                    memmove(buf+pos+1,buf+pos,len-pos);
                    buf[pos] = c;
                    len++;
                    pos++;
                    buf[len] = '\0';
                    refreshLine(fd,prompt,buf,len,pos,cols);
                }
            }
            break;
        case 21: /* Ctrl+u, delete the whole line. */
            buf[0] = '\0';
            pos = len = 0;
            refreshLine(fd,prompt,buf,len,pos,cols);
            break;
        case 11: /* Ctrl+k, delete from current to end of line. */
            buf[pos] = '\0';
            len = pos;
            refreshLine(fd,prompt,buf,len,pos,cols);
            break;
        case 1: /* Ctrl+a, go to the start of the line */
            pos = 0;
            refreshLine(fd,prompt,buf,len,pos,cols);
            break;
        case 5: /* ctrl+e, go to the end of the line */
            pos = len;
            refreshLine(fd,prompt,buf,len,pos,cols);
            break;
        }
    }
    return len;
}

static int linenoiseRaw(char *buf, size_t buflen, const char *prompt) {
    int fd = STDIN_FILENO;
    int count;

    if (buflen == 0) {
        errno = EINVAL;
        return -1;
    }
    if (!isatty(STDIN_FILENO)) {
        if (fgets(buf, buflen, stdin) == nullptr) return -1;
        count = strlen(buf);
        if (count && buf[count-1] == '\n') {
            count--;
            buf[count] = '\0';
        }
    } else {
        if (enableRawMode(fd) == -1) return -1;
        count = linenoisePrompt(fd, buf, buflen, prompt);
        disableRawMode(fd);
        printf("\n");
    }
    return count;
}

char *linenoise(const char *prompt) {
    char buf[LINENOISE_MAX_LINE];

    if (isUnsupportedTerm()) {
        printf("%s",prompt);
        fflush(stdout);
        if (fgets(buf,LINENOISE_MAX_LINE,stdin) == nullptr) return nullptr;
        size_t len = strlen(buf);
        while(len && (buf[len-1] == '\n' || buf[len-1] == '\r')) {
            len--;
            buf[len] = '\0';
        }
        return strdup(buf);
    } else {
        int count = linenoiseRaw(buf,LINENOISE_MAX_LINE, prompt);
        if (count == -1) return nullptr;
        return strdup(buf);
    }
}

/* Using a circular buffer is smarter, but a bit more complex to handle. */
int linenoiseHistoryAdd(const char *line) {
    if (history_max_len == 0) return 0;
    if (history == nullptr) {
        history = (char**)malloc(sizeof(char*)*history_max_len);
        if (history == nullptr) return 0;
        memset(history,0,(sizeof(char*)*history_max_len));
    }
    char *linecopy = strdup(line);
    if (!linecopy) return 0;
    if (history_len == history_max_len) {
        free(history[0]);
        memmove(history,history+1,sizeof(char*)*(history_max_len-1));
        history_len--;
    }
    history[history_len] = linecopy;
    history_len++;
    return 1;
}

int linenoiseHistorySetMaxLen(int len) {
    if (len < 1) return 0;
    if (history) {
        int tocopy = history_len;

        char **newChar = (char **) malloc(sizeof(char *) * len);
        if (newChar == nullptr) return 0;
        if (len < tocopy) tocopy = len;
        memcpy(newChar,history+(history_max_len-tocopy), sizeof(char*)*tocopy);
        free(history);
        history = newChar;
    }
    history_max_len = len;
    if (history_len > history_max_len)
        history_len = history_max_len;
    return 1;
}

/* Save the history in the specified file. On success 0 is returned
* otherwise -1 is returned. */
int linenoiseHistorySave(char *filename) {
    FILE *fp = fopen(filename,"w");
    int j;
    
    if (fp == nullptr) return -1;
    for (j = 0; j < history_len; j++)
        fprintf(fp,"%s\n",history[j]);
    fclose(fp);
    return 0;
}

/* Load the history from the specified file. If the file does not exist
* zero is returned and no operation is performed.
*
* If the file exists and the operation succeeded 0 is returned, otherwise
* on error -1 is returned. */
int linenoiseHistoryLoad(char *filename) {
    FILE *fp = fopen(filename,"r");
    char buf[LINENOISE_MAX_LINE];
    
    if (fp == nullptr) return -1;

    while (fgets(buf,LINENOISE_MAX_LINE,fp) != nullptr) {
        char *p;
        
        p = strchr(buf,'\r');
        if (!p) p = strchr(buf,'\n');
        if (p) *p = '\0';
        linenoiseHistoryAdd(buf);
    }
    fclose(fp);
    return 0;
}
