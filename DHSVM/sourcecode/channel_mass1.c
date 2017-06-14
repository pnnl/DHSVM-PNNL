/*
 * SUMMARY:      channel_mass1.c
 * USAGE:        Part of DHSVM
 *
 * AUTHOR:       William A. Perkins
 * ORG:          Pacific NW National Laboratory
 * E-MAIL:       william.perkins@pnnl.gov
 * ORIG-DATE:    June 2017
 * DESCRIPTION:  Reads a DHSVM stream network and makes a MASS1 network. 
 *
 * DESCRIP-END.cd
 * FUNCTIONS:    
 * LAST CHANGE: 2017-06-14 07:41:06 d3g096
 * COMMENTS:
 */


#include <stdlib.h>
#include <unistd.h>
#include <stdio.h>
#include <libgen.h>
#include <string.h>
#include <sys/param.h>
#include <math.h>

#include "errorhandler.h"
#include "channel.h"

/******************************************************************************/
/*                             channel_compute_elevation                      */
/******************************************************************************/
void
channel_compute_elevation(Channel *network, float elev0)
{
  Channel *current = NULL;
  int max_order, o;

  max_order = 0;
  for (current = network; current != NULL; current = current->next) {
    if (current->order > max_order) max_order = current->order;
  }

  error_handler(ERRHDL_DEBUG, "computing channel elevations (maxorder = %d)",
                max_order);

  for (o = max_order; o > 0; --o) {
    for (current = network; current != NULL; current = current->next) {
      if (current->order == o) {
        if (current->outlet) {
          current->outlet_elevation = current->outlet->inlet_elevation;
        } else {
          current->outlet_elevation = elev0;
        }
        current->inlet_elevation = current->outlet_elevation + 
          current->length*current->slope;
      }
    }
  }
}

/******************************************************************************/
/*                           mass1_write_sections                             */
/******************************************************************************/
void
mass1_write_sections(const char *outname, ChannelClass *classes)
{
  ChannelClass *current;
  char outfile[MAXPATHLEN];
  FILE *out;
  float width;
  
  sprintf(outfile, "%ssection.dat", outname);
  if ((out = fopen(outfile, "wt")) == NULL) {
    error_handler(ERRHDL_ERROR, "cannot open section file \"%s\"",
                  outfile);
    return;
  }

  error_handler(ERRHDL_DEBUG, "writing MASS1 cross sections to \"%s\"",
                outfile);

  for (current = classes; current != NULL; current = current->next) {
    fprintf(out, "%d     1\n%.2f /\n", current->id, current->width);
  }

  fclose(out);
}

/******************************************************************************/
/*                          mass1_write_links                                 */
/******************************************************************************/
void
mass1_write_links(const char *outname, Channel *network, const float spacing)
{
  char outfile[MAXPATHLEN];
  FILE *out;
  Channel *current;

  sprintf(outfile, "%slink.dat", outname);
  if ((out = fopen(outfile, "wt")) == NULL) {
    error_handler(ERRHDL_ERROR, "cannot open link file \"%s\"",
                  outfile);
    return;
  }

  error_handler(ERRHDL_DEBUG, "writing MASS1 link information to \"%s\"",
                outfile);

  for (current = network; current != NULL; current = current->next) {
    int npts;
    npts = (int)truncf(current->length/spacing + 0.5);
    if (npts < 2) npts = 2;
    fprintf(out, "%5d", current->id); 
    fprintf(out, " %5d", 2);              /* input option */
    fprintf(out, " %5d", npts);           /* number of points */
    fprintf(out, " %5d", current->order); /* link order (unused) */
    fprintf(out, " %5d", 0);              /* upstream links (unused) */

                                          /* upstream BC */
    if (current->order > 1) {
      fprintf(out, " %5d", 0);
    } else {
      fprintf(out, " %5d", 1);
    }
                                          /* downstream BC */
    if (current->outlet == NULL) {
      fprintf(out, " %5d", 2);
    } else {
      fprintf(out, " %5d", 0);
    }
    
    fprintf(out, " %5d", 0);              /* temperature BC */
    fprintf(out, " %5d", current->id);    /* met zone */
                                          /* lateral inflow, TDG, temp BC */
    fprintf(out, " %5d %5d %5d", current->id, 0, current->id);
    fprintf(out, " %5.1f", 3.5);          /* LPI coefficient */

    fprintf(out, " /\n");

    if (current->outlet == NULL) {
      fprintf(out, " %5d", 0);
    } else {
      fprintf(out, "%5d", current->outlet->id); 
    }
    fprintf(out, "%72.72s /\n", " ");
  }
  fclose(out);
}

/******************************************************************************/
/*                           mass1_write_points                               */
/******************************************************************************/
void
mass1_write_points(const char *outname, Channel *network, const float spacing)
{
  char outfile[MAXPATHLEN];
  FILE *out;
  Channel *current;

  sprintf(outfile, "%spoint.dat", outname);
  if ((out = fopen(outfile, "wt")) == NULL) {
    error_handler(ERRHDL_ERROR, "cannot open point file \"%s\"",
                  outfile);
    return;
  }

  error_handler(ERRHDL_DEBUG, "writing MASS1 point information to \"%s\"",
                outfile);

  for (current = network; current != NULL; current = current->next) {
    int npts;
    npts = (int)truncf(current->length/spacing + 0.5);
    if (npts < 2) npts = 2;

    fprintf(out, "%5d", current->id);                   /* link id */
    fprintf(out, " %10.2f", current->length);           /* link length */
    fprintf(out, " %10.2f", current->inlet_elevation);  /* upstream elevation */
    fprintf(out, " %10.2f", current->outlet_elevation); /* downstream elevation */
    fprintf(out, " %5d", current->class2->id);          /* section id */
    fprintf(out, " %10.4f", current->class2->friction); /* Manning's n */
    fprintf(out, " %10.1f", 300.0);                     /* longitudinal dispersion */
    fprintf(out, " %10.4f", 0.0);                       /* unused */
    fprintf(out, " /\n");
  }
  fclose(out);
}

/******************************************************************************/
/*                          mass1_write_bclists                               */
/******************************************************************************/
void
mass1_write_bclists(const char *outname, Channel *network)
{
  char outfile[MAXPATHLEN];
  FILE *out;
  Channel *current;

  sprintf(outfile, "%linkbc.dat", outname);
  if ((out = fopen(outfile, "wt")) == NULL) {
    error_handler(ERRHDL_ERROR, "cannot open link BC file \"%s\"",
                  outfile);
    return;
  }
  fclose(out);

  sprintf(outfile, "%lateral.dat", outname);
  if ((out = fopen(outfile, "wt")) == NULL) {
    error_handler(ERRHDL_ERROR, "cannot open lateral inflow file \"%s\"",
                  outfile);
    return;
  }
  fclose(out);
}


/******************************************************************************/
/*                            Main Program                                    */
/******************************************************************************/
int
main(int argc, char **argv)
{
  char *program;
  char ch;
  char class_file[MAXPATHLEN], network_file[MAXPATHLEN];
  char outname[MAXPATHLEN];
  ChannelClass *classes, *c;
  Channel *network, *l;
  int n, ierr;
  int maxid;
  float spacing;
  float elev0;

  const char usage[] = "usage: [-v] [-s spacing] [-n roughness] [-o basename] %s class.dat network.dat";

  program = basename(argv[0]);

  error_handler_init(program, NULL, ERRHDL_MESSAGE);

  strncpy(outname, "", MAXPATHLEN);
  spacing = 250.0;
  elev0 = 0.0;
  ierr = 0;

  while ((ch = getopt(argc, argv, "vs:o:")) != -1) {
    switch (ch) {
    case 'v': 
      error_handler_init(program, NULL, ERRHDL_DEBUG);
      break;
    case 's':
      spacing = strtof(optarg, NULL);
      if (spacing <= 0.0) {
        error_handler(ERRHDL_ERROR, "spacing \"%s\" not understood", optarg);
        ierr++;
      }
      break;
    case 'o':
      strncpy(outname, optarg, MAXPATHLEN);
      break;
    default:
      error_handler(ERRHDL_FATAL, usage, program);
      exit(3);
    }
  }
  argc -= optind;
  argv += optind;

  if (argc < 2) {
    error_handler(ERRHDL_FATAL, usage, program);
    exit(2);
  } else {
    strncpy(class_file, argv[0], MAXPATHLEN);
    strncpy(network_file, argv[1], MAXPATHLEN);
  }

  error_handler(ERRHDL_DEBUG, "nominal section spacing = %.1f", spacing);
  error_handler(ERRHDL_DEBUG, "reading channel classes from %s...",
                class_file);
  classes = channel_read_classes(class_file, 0);
  c = classes;
  n = 0;
  while (c != NULL) {
    ++n;
    c = c->next;
  }
  error_handler(ERRHDL_DEBUG, "%d channel classes read.", n);


  error_handler(ERRHDL_DEBUG, "reading channel segments from %s...",
                network_file);
  network = channel_read_network(network_file, classes, &maxid);
  l = network;
  n = 0;
  while (l != NULL) {
      ++n;
      l = l->next;
  }
  error_handler(ERRHDL_DEBUG, "%d channel segments read.", n);

  channel_compute_elevation(network, elev0);

  mass1_write_sections(outname, classes);
  mass1_write_links(outname, network, spacing);
  mass1_write_points(outname, network, spacing);

  channel_free_network(network);
  channel_free_classes(classes);
  
  error_handler_done();
}
