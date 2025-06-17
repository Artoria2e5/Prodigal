/*******************************************************************************
    PRODIGAL (PROkaryotic DynamIc Programming Genefinding ALgorithm)
    Copyright (C) 2007-2016 University of Tennessee / UT-Battelle

    Code Author:  Doug Hyatt

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
*******************************************************************************/

#ifndef _TRAIN_H
#define _TRAIN_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* FIXME: This is a giant 558 KB struct which is awful for locality. Many bytes
   can be saved by replacing double with float, and a couple more by squeezing
   the two ints together. That would however change the on-disk format which is
   a big bother. (Either you break compatibility or you write a conversion
   routine.)
   
   Just appending should not break things given some care. */
struct _training {
  double gc;                    /* GC Content */
  int trans_table;              /* For compatibility: NCBI trans table ID */
  double st_wt;                 /* Start weight */
  double bias[3];               /* GC frame bias for each of the 3 positions */
  double type_wt[3];            /* Weights for ATG vs GTG vs TTG */
  int uses_sd;                  /* 0 if doesn't use SD motif, 1 if it does */
  double rbs_wt[28];            /* Set of weights for RBS scores */
  double ups_comp[32][4];       /* Base composition weights for non-RBS-distance
                                   motifs.  0-1 are the -1/-2 position, 2-31 are
                                   the -15 to -44 positions.  Second array is
                                   the base A,C,T,G,etc. */
  double mot_wt[4][4][4096];    /* Weights for upstream motifs.  First index is
                                   the motif length (3-6), the second is the
                                   spacer distance (0 = 5-10bp, 1 = 3-4bp, 2 =
                                   11-12bp, 3 = 13-15bp), and the last is the
                                   numerical value of the motif (ranging from 0
                                   to 4095 for 6-mers, less for shorter
                                   motifs) */
  double no_mot;                /* Weight for the case of no motif */
  double gene_dc[4096];         /* Coding statistics for the genome */
  int    version;               /* Future-proofing */
  char   table[64];             /* Translation table as trinuc()-to-AA mapping */ 
};

/* The old struct for figuring out where the original file should end.
   Can't just subtract by 64 due to struct layout magic */
struct _training_v0 {
  double gc;
  int trans_table;
  double st_wt;
  double bias[3];
  double type_wt[3];
  int uses_sd;
  double rbs_wt[28];
  double ups_comp[32][4];
  double mot_wt[4][4][4096];
  double no_mot;
  double gene_dc[4096];
};

void meta_to_v1(struct _training *);
void v0_to_v1(struct _training_v0 *old);
int write_training_file(char *, const struct _training *);
int read_training_file(char *, struct _training *);

#endif
