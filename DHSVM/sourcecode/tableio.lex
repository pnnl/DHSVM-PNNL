%{
/* -------------------------------------------------------------
   file: tableio.lex
   ------------------------------------------------------------- */
/* -------------------------------------------------------------
   Battelle Memorial Institute
   Pacific Northwest Laboratory
   ------------------------------------------------------------- */
/* -------------------------------------------------------------
   Created January 22, 1996 by William A. Perkins
   $Id: tableio.lex,v 1.4 2003/07/01 21:26:33 olivier Exp $
   ------------------------------------------------------------- */


const char* TABLE_LEX_ID = "$Id: tableio.lex,v 1.4 2003/07/01 21:26:33 olivier Exp $";

#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include "tableio.h"
#include "errorhandler.h"
#include "settings.h"

static FILE *table_file = NULL;
static char *table_file_name = NULL;
static int table_lines = 1;
int table_errors = 0;
int table_warnings = 0;
static int table_field;
static int comment_caller, string_caller;

static char string_buffer[1024];

typedef enum {
  TABLE_END = -1,
  TABLE_RECORD = 0,
  TABLE_FIELD = 1
} TableStatus;

typedef struct {  
  int type;
  int integer;
  float real;
  char string[TABLE_MAX_FIELD_LEN];
} TableScanField;



/* -------------------------------------------------------------
   table_yy_lex 
   This is how the generated scanner is declared.  If an end of record
   is reached, zero is returned and field->type is set only if an
   unterminated string is found. At EOF, -1 is returned and
   field->type is set only for an unterminated string.  1 is returned,
   and field->type is set if a field is successfully read.
   ------------------------------------------------------------- */
#define YY_DECL static TableStatus table_yy_lex(TableScanField *field)

%}
%option noyywrap
%option nostdinit
%x comment_start record_start string_start


WHITE	[ \t]+
INTEGER [0-9]+
REAL 	[-\+]?[0-9]*\.[0-9]+([-+]?[EeDd][0-9]+)?
WORD 	[A-Za-z][A-Za-z0-9_]+
HASH    #

%%


<string_start>\\\"            { strncat(string_buffer, "\"", 1024); }
<string_start>\\              { strncat(string_buffer, "\\", 1024); }
<string_start>[^\\\"\n]+      { strncat(string_buffer, yytext, 1024); }
<string_start>\"              { 
  field->type = TABLE_STRING;
  strncpy(field->string, string_buffer, TABLE_MAX_FIELD_LEN);
  BEGIN(string_caller);
  table_field++; 
  return(TABLE_FIELD);
}

<string_start>\n              { 
  field->type = TABLE_STRING;
  strncpy(field->string, string_buffer, TABLE_MAX_FIELD_LEN);
  BEGIN(string_caller);
  table_field++; 
  error_handler(ERRHDL_WARNING, 
                "%s: line %d, field %d: unterminated string field",
                table_file_name, table_lines, table_field);
  table_warnings++;
  yyless(0);
  return (TABLE_FIELD);
}

<record_start>{INTEGER}       { 
  field->type = TABLE_INTEGER;
  field->integer = strtol(yytext, NULL, 10);
  ++table_field; 
  return (TABLE_FIELD);
}

<record_start>{REAL}          { 
  field->type = TABLE_REAL;
  field->real = strtod(yytext, NULL);
  ++table_field;
  return (TABLE_FIELD);
}

<record_start>{WORD}          { 
  field->type = TABLE_WORD;
  strncpy(field->string, yytext, TABLE_MAX_FIELD_LEN);
  ++table_field; 
  return (TABLE_FIELD);
}

<record_start>\"        { string_caller = YY_START; BEGIN(string_start); 
                          string_buffer[0] = '\0'; }

<record_start>{WHITE}+        { /* eat it up */ }
  
<record_start>\n              { 
  ++table_lines;
  BEGIN(INITIAL); 
  return (TABLE_RECORD);
}

<INITIAL>\"             { BEGIN(record_start); table_field = -1; yyless(0); }

<INITIAL>{INTEGER}|{REAL}|{WORD} { BEGIN(record_start); table_field = -1; yyless(0); }

<INITIAL,record_start>{HASH}+        { comment_caller = YY_START; BEGIN(comment_start); }

<comment_start>{HASH}+         { /* eat up any of these within the comment */ }
<comment_start>[^\n]+          { /* eat up comments text */ }
<comment_start>\n              { BEGIN(comment_caller); yyless(0); }

<*><<EOF>>      { 
  switch (YY_START) {
  case record_start:
  case comment_start:
    break;
  case string_start:
    field->type = TABLE_STRING;
    strncpy(field->string, string_buffer, TABLE_MAX_FIELD_LEN);
    table_warnings++;
    table_field++; 
    error_handler(ERRHDL_WARNING, 
                  "%s: line %d, field %d: unterminated string",
                  table_file_name, table_lines, table_field);
    break;
  default:
    table_field = -1;
    break;
  }
  BEGIN(INITIAL);
  return(TABLE_END);
}

{WHITE}         { /* eat up extra white space */ }
\n		{ table_lines++;}

%%

/* -------------------------------------------------------------
   table_open

   opens the specified file and sets up the scanner.  Zero is returned
   if all goes well
   ------------------------------------------------------------- */
int
table_open(const char *filename) 
{
  if (table_file != NULL) {
    error_handler(ERRHDL_ERROR,
                 "table_open: already reading another table: %s",
                 table_file_name);
    return (-1);
  } 
  if ((table_file = fopen(filename, "r")) == NULL) {
    error_handler(ERRHDL_ERROR,
                 "table_open: can not open file \"%s\": %s",
                 filename, strerror(errno));
    table_file = NULL;
    return (-1);
  }
  table_file_name = (char *)strdup(filename);
  table_lines = 1;
  table_errors = 0;  
  table_warnings = 0;
  yyrestart(table_file);
  return (0);
}
  

/* -------------------------------------------------------------
   handle_field
   ------------------------------------------------------------- */
static int
handle_field(TableScanField *field, int nfields, TableField *fld_defn) 
{
  int i;
  int err = 0;
  
  if (field->type == fld_defn[0].type) {
    switch (field->type) {

    case TABLE_INTEGER:

      fld_defn[0].value.integer = field->integer;
      fld_defn[0].read = TRUE;
      break;

    case TABLE_REAL:

      fld_defn[0].value.real = field->real;
      fld_defn[0].read = TRUE;
      break;

    case TABLE_STRING:

      strncpy(fld_defn[0].field, field->string, TABLE_MAX_FIELD_LEN);
      fld_defn[0].read = TRUE;
      break;

    case TABLE_WORD:

      for (i = 0; fld_defn[0].words[i] != NULL; i++) {
        if (!strncmp(fld_defn[0].words[i], 
                     field->string, TABLE_MAX_FIELD_LEN)) {
          fld_defn[0].value.integer = i;
          fld_defn[0].read = TRUE;
          break;
        }
      }
      if (!fld_defn[0].read) {
        err++;
        error_handler(ERRHDL_ERROR, 
                      "%s: line %d, field %d: unknown keyword \"%s\"",
                      (const char *)table_file_name, 
                      (int)table_lines, (int)table_field, 
                      (const char *)field->string);
      }
      break;

    default:

      error_handler(ERRHDL_FATAL, 
                    "handle_field: what is this status: %d?", 
                    field->type);
      break;
          
    }
      
  } else if (field->type == TABLE_INTEGER && 
             fld_defn[0].type == TABLE_REAL) {

    fld_defn[0].value.real = field->integer;
    fld_defn[0].read = TRUE;

  } else if (fld_defn[0].required) {

    err++;
    error_handler(ERRHDL_ERROR, 
                  "%s: line %d, field %d: type mismatch",
                  table_file_name, table_lines, table_field);

  } else {

    error_handler(ERRHDL_STATUS,
                  "%s: line %d, field %d: type mismatch for optional \"%s\"",
                  table_file_name, table_lines, table_field,
                  fld_defn[0].name);
    if (--nfields > 0) {
      err += handle_field(field, nfields, &fld_defn[1]);
    } else {
      error_handler(ERRHDL_ERROR, 
                    "%s: line %d, field %d: unable to match field",
                    table_file_name, table_lines, table_field);
    }

  }
  return (err);
}
  

/* -------------------------------------------------------------
   table_get_fields

   This is used to read a record from the file and decide if errors
   have occurred. Zero is returned if a record is read OK, -1 if the
   end of file was reached, a positive integer if errors occured.
   This assumes table_open has been called.
   ------------------------------------------------------------- */
int
table_get_fields(int nfields, TableField *fld_defn) 
{
  int i, required;
  int err = 0;
  char done = FALSE;
  TableStatus status = TABLE_END;
  TableScanField field;

  error_handler(ERRHDL_DEBUG, 
                "%s: line %d: table_scan_record: looking for %d fields",
                table_file_name, table_lines, nfields);

  for (i = 0, required = 0; i < nfields; i++) {
    if (fld_defn[i].required) required++;
    fld_defn[i].read = FALSE;
    fld_defn[i].value.real = 0.0;
    strcpy(fld_defn[i].field, "\0");
  }
  
  while (!done) { 
    switch (status = table_yy_lex(&field)) {

    case TABLE_END:
    case TABLE_RECORD:
      if (table_field >= 0) {
        for (i = 0; i < nfields; i++) {
          if (fld_defn[i].required && !fld_defn[i].read) {
            err++;
            error_handler(ERRHDL_ERROR,
                          "%s: line %d: \"%s\" (field %d) required but not read",
                          table_file_name, table_lines,
                          fld_defn[i].name, i + 1);
          }
        }
      }
      done = TRUE;
      break;

    case TABLE_FIELD:

      if (table_field >= nfields) {
        err++;
        error_handler(ERRHDL_ERROR,
                      "%s: line %d, field %d: too many fields, expected %d",
                      table_file_name, table_lines, table_field, nfields);
        continue;
      }

      err += handle_field(&field, nfields - table_field, &fld_defn[table_field]);
      break;

    default:

      error_handler(ERRHDL_FATAL, 
                    "table_get_fields: what is this status: %d?", 
                    status);
      break;
          
    }
  }
  
  table_errors += err;

  if (status == TABLE_END)
    return (-1);

  return (err);
}

/* -------------------------------------------------------------
   table_lineno
   ------------------------------------------------------------- */
int
table_lineno(void)
{
  return table_lines;
}



/* -------------------------------------------------------------
   table_close
   ------------------------------------------------------------- */
void
table_close(void) 
{
  if (fclose(table_file) != 0) {
    error_handler(ERRHDL_WARNING,
                  "table_close: %s: error closing",
                  table_file_name, strerror(errno));
  }
  table_file = NULL;
  if (table_file_name != NULL) {
    free(table_file_name);
    table_file_name = NULL;
  }
  
}


#ifdef TEST_MAIN

typedef enum { ONE, TWO, THREE } RecordWord;
typedef struct _record_ {
  int ID;
  float value;
  RecordWord item;
  char *comment;
  struct _record_ *next;
} Record;

/* -------------------------------------------------------------
   Main Program
   ------------------------------------------------------------- */
int
main(int argc, char **argv)
{
  int status, i, rec, done;
  static char *crown_words[4] = { 
    "ONE", "TWO", "THREE", NULL 
  };
  static const int nfields = 4;
  static TableField fldesc[4] = {
    {
      "ID",
      TABLE_INTEGER,
      TRUE, FALSE,
      0.0, "", NULL
    },
    {
      "Channel Width",
      TABLE_REAL,
      TRUE, FALSE,
      0.0, "", NULL
    },
    {
      "Road Crown Type",
      TABLE_WORD,
      FALSE, FALSE,
      0.0, "", crown_words
    },
    {
      "Comment",
      TABLE_STRING,
      FALSE, FALSE,
      0.0, "", NULL
    }
  };
  static const char *file = "table_test.dat";

  error_handler_init(argv[0], NULL, ERRHDL_DEBUG);
  
  if (table_open(file) != 0) {
    error_handler(ERRHDL_ERROR, 
                  "channel_read_classes: unable to open file \"%s\": %s",
                  file, strerror(errno));
    exit(1);
  }

  printf("Successfully read stuff:\n\n");
    
  rec = 0;
  done = FALSE;
  while (!done) {
    done = ((status = table_get_fields(nfields, fldesc)) < 0);
    printf("Record %d:\n", ++rec);
    for (i = 0; i < nfields; i++) {
      printf("\t%20s:\t", fldesc[i].name);
      if (fldesc[i].read) {
        switch (fldesc[i].type) {
        case TABLE_INTEGER:
          printf("%d (INTEGER)\n", fldesc[i].value.integer);
          break;
        case TABLE_REAL:
          printf("%g (REAL)\n", (double)fldesc[i].value.real);
          break;
        case TABLE_WORD:
          printf("%s (WORD, index = %d)\n", 
                 fldesc[i].words[fldesc[i].value.integer],
                 fldesc[i].value.integer);
          break;
        case TABLE_STRING:
          printf("%s (STRING)\n", fldesc[i].field);
          break;
        default:
          error_handler(ERRHDL_ERROR,
                        "main: what is this type: %d?",
                        fldesc[i].type);
          break;
        }
      } else {
        printf("not read\n");
      }
    }
  }
  
  error_handler(ERRHDL_MESSAGE,
                "%s: %d errors, %d warnings",
                file, table_errors, table_warnings);
  error_handler_done();
  exit(0);
}



#endif
