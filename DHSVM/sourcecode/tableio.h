/* -------------------------------------------------------------
   file: tableio.h
   ------------------------------------------------------------- */
/* -------------------------------------------------------------
   Battelle Memorial Institute
   Pacific Northwest Laboratory
   ------------------------------------------------------------- */
/* -------------------------------------------------------------
   Created January 12, 1996 by  William A Perkins
   Last Change: Thu Feb 15 16:02:18 1996 by  William A Perkins <perk@yama.pnl.gov>
   ------------------------------------------------------------- */

/* $Id: tableio.h,v 1.4 2003/07/01 21:26:33 olivier Exp $ */

#ifndef _tableio_h_
#define _tableio_h_

#define TABLE_MAX_FIELD_LEN 128

typedef enum { 
  TABLE_INTEGER = 1, 
  TABLE_REAL = 2, 
  TABLE_STRING = 3, 
  TABLE_WORD = 4 
} TableFieldType;

typedef struct {
  const char *name;
  TableFieldType type;
  char required;
  char read;
  union {
    int integer;
    float real;
  } value;
  char field[TABLE_MAX_FIELD_LEN];
  char **words;
} TableField;

/* -------------------------------------------------------------
   externally available variables
   ------------------------------------------------------------- */

extern int table_errors;
extern int table_warnings;

/* -------------------------------------------------------------
   externally available functions
   ------------------------------------------------------------- */
int table_open(const char *filename);
int table_get_fields(int num, TableField *fields);
int table_lineno(void);
void table_close(void);

#endif
