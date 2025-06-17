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

#include "training.h"
#include "table.h"
#include "bitmap.h"
#include <stdio.h>


/* Updates meta training data to v1 (flip 4096, table) */
void meta_to_v1(struct _training *tinf) {
  double old_mot_wt[4][4][4096];
  double old_gene_dc[4096];

  memcpy(old_mot_wt, tinf->mot_wt, sizeof(tinf->mot_wt));
  memcpy(old_gene_dc, tinf->gene_dc, sizeof(tinf->gene_dc));
  for (unsigned i = 0; i < 4096; i++) {
    unsigned old_i = range_preserving_bitrev(i, 12);
    for (unsigned j = 0; j < 4; j++) {
      for (unsigned k = 0; k < 4; k++) {
        tinf->mot_wt[j][k][i] = old_mot_wt[j][k][old_i];
      }
    }
    tinf->gene_dc[i] = old_gene_dc[old_i];
  }

  tinf->version = 1;
  id_to_table(tinf->table, tinf->trans_table);
}

#define COPY(f)  tinf.f = old->f;
#define CLONE(f) memcpy(tinf.f, old->f, sizeof(old->f));
/* Updates on-disk training data to v1 (copy fields then meta_to_v1()) */
void v0_to_v1(struct _training_v0 *old) {
  struct _training tinf;
  COPY(gc);
  COPY(trans_table);
  COPY(st_wt);
  CLONE(bias);
  CLONE(type_wt);
  COPY(uses_sd);
  CLONE(rbs_wt);
  CLONE(ups_comp);
  CLONE(mot_wt);
  COPY(no_mot);
  CLONE(gene_dc);
  tinf.version = 0;
  memset(tinf.table, '\0', sizeof(tinf.table));
  meta_to_v1(&tinf);
  memcpy(old, &tinf, sizeof(tinf));
}
#undef COPY
#undef CLONE

/* Reads a training file to use for gene prediction */
int read_training_file(char *fn, struct _training *tinf) {
  size_t rv;
  FILE *fh;
  int ret = 0;
  fh = fopen(fn, "rb");
  if(fh == NULL) return 1;
  rv = fread(tinf, 1, MAX_TRAINING_SIZE, fh);
  switch (rv) {
    case sizeof(struct _training):
      break;
    case sizeof(struct _training_v0):
      v0_to_v1((struct _training_v0 *) tinf);
      break;
    default:
      fprintf(stderr, "Read %zu bytes (should be v1 %zu or v0 %zu)\n",
        rv, sizeof(struct _training), sizeof(struct _training_v0));
      ret = -1; // Unexpected size
      break;
  }
  fclose(fh);
  return ret;
}

/* Writes a training file to use for a later run of gene prediction */
int write_training_file(char *fn, const struct _training *tinf) {
  size_t rv;
  FILE *fh;
  fh = fopen(fn, "wb");
  if(fh == NULL) return -1;
  rv = fwrite(tinf, sizeof(struct _training), 1, fh);
  fclose(fh);
  if(rv != 1) return -1;
  return 0;
}
