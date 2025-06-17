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
#include <stdint.h>

/* FIXME: This is a giant 558,392 B struct which is awful for locality.
   There is no squeezing it in a 32 KB L1D cache anyways though.
   
   The log prob weights are not very sensitive to precision AND has a
   guaranteed range limit, so turning them into Fixed16 should be fine.
   Same for GC%. Assuming we do that with bias (always sums to 3 and non-neg?)
   too we get 139,662 B.

   mot_wt can be changed to a [4][4096 + 1024 + 256 + 64], bringing us
   down to 52,110 B... Huh, worth a try some day!
   */
struct _training
{
  uint8_t version;           /* Future-proofing, now = 1 (nucmer, table) */
  uint8_t trans_table;       /* NCBI trans table ID, default = 11, custom = -1 */
  uint8_t uses_sd;           /* 0 if doesn't use SD motif, 1 if it does */
  uint8_t padding0;
  char table[64];            /* Translation table as trinuc()-to-AA mapping */
  double gc;                 /* GC Content */
  double bias[3];            /* GC frame bias for each of the 3 positions */
  double st_wt;              /* Start weight */
  double type_wt[3];         /* Weights for ATG vs GTG vs TTG */
  double rbs_wt[28];         /* Set of weights for RBS scores */
  double ups_comp[32][4];    /* Base composition weights for non-RBS-distance
                                motifs.  0-1 are the -1/-2 position, 2-31 are
                                the -15 to -44 positions.  Second array is
                                the base A,C,T,G,etc. */
  double mot_wt[4][4][4096]; /* Weights for upstream motifs.  First index is
                                the motif length (3-6), the second is the
                                spacer distance (0 = 5-10bp, 1 = 3-4bp, 2 =
                                11-12bp, 3 = 13-15bp), and the last is the
                                numerical value of the motif (ranging from 0
                                to 4095 for 6-mers, less for shorter
                                motifs) */
  double no_mot;             /* Weight for the case of no motif */
  double gene_dc[4096];      /* Coding statistics for the genome */
};

/* The old struct for reading old files. */
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

#define MAX(a, b) ((a) > (b) ? (a) : (b))
#define MAX_TRAINING_SIZE MAX(sizeof(struct _training), sizeof(struct _training_v0))

void meta_to_v1(struct _training *);
void v0_to_v1(struct _training_v0 *old);
int write_training_file(char *, const struct _training *);
int read_training_file(char *, struct _training *);

#endif
