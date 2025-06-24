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

/* FIXME: This is a 208,240 B struct which is awful for locality.
   There is no squeezing it in a 32 KB L1D cache, but now that it's in
   L2 we've saved a lot of L3 miss cycles (about 20!).
   
   Possible further shaves:
    float    -> 104,152B
    fixed16  -> 52,110B   (our uses here are not very sensitive to precision
                           and have a guaranteed range limit: log prob weights,
                           GC%, nonnegative bias summing to 3, etc.)
   */
#define NMOTIF (64 + 256 + 1024 + 4096)
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
  double mot_wt[4][NMOTIF];  /* Weights for upstream motifs.  First index is
                                the spacer distance (0 = 5-10bp, 1 = 3-4bp, 2 =
                                11-12bp, 3 = 13-15bp), and the last is the
                                numerical value of the motif plus a length-
                                dependent offset (0-63 for 3-mer, 64-319
                                for 4-mer, 320-1344 for 5-mer, 1345-5439
                                for 6-mer). */
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
  double mot_wt[4][4][4096]; /* Old style weights for upstream motifs.
                                First index is length+3
                                Second is spacer
                                Third is numeric value of motif */
  double no_mot;
  double gene_dc[4096];
};

#define MAX(a, b) ((a) > (b) ? (a) : (b))
#define MAX_TRAINING_SIZE MAX(sizeof(struct _training), \
                              sizeof(struct _training_v0))

/* convert mot_wt[n][x][idx] to mot_wt[x][i] and back;
  n \in [0..3], idx \in [0..4095], i \in [0..5439] */
inline static uint16_t mot_idx_v0_to_v1(uint16_t n, uint16_t idx) {
  uint16_t mot_wt_02[] = {0, 64, 320, 1344};
  return mot_wt_02[n] + idx;
}
inline static uint32_t mot_idx_v1_to_v0(uint16_t i) {
  uint16_t mot_wt_02[] = {0, 64, 320, 1344};
  uint16_t n = (i >= 64) + (i >= 320) + (i >= 1344);
  return (n << 16) | (i - mot_wt_02[n]);
}

/* Convert in situ knowing that the new struct is smaller */
void v0_to_v1(struct _training_v0 *old);
int write_training_file(char *, const struct _training *);
int read_training_file(char *, struct _training *);

#endif
