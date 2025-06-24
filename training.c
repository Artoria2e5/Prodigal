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
#include <math.h>
#include <stdint.h>
#include <float.h>


#define COPY(f)  tinf.f = old->f;
#define CLONE(f) for (unsigned i = 0; i < sizeof(tinf.f)/sizeof(tinf.f[0]); i++) { \
  tinf.f[i] = old->f[i]; \
}
/* Updates on-disk training data to v1 */
void v0_to_v1(struct _training_v0 *old) {
  struct _training tinf;
  COPY(gc);
  COPY(trans_table);
  COPY(st_wt);
  CLONE(bias);
  CLONE(type_wt);
  COPY(uses_sd);
  CLONE(rbs_wt);
  CLONE(ups_comp[0]);
  CLONE(ups_comp[1]);
  CLONE(ups_comp[2]);
  CLONE(ups_comp[3]);
  COPY(no_mot);
  for (unsigned i = 0; i < 4096; i++) {
    unsigned old_i = range_preserving_bitrev(i, 12);
    tinf.gene_dc[i] = old->gene_dc[old_i];
  }
  for (unsigned i = 0; i < 64 + 256 + 1024 + 4096; i++) {
    unsigned idx = mot_idx_v1_to_v0(i);
    unsigned n = idx >> 16;
    idx = idx & 0xFFFF;
    for (unsigned j = 0; j < 4; j++) {
      tinf.mot_wt[j][i] = old->mot_wt[n][j][idx];
    }
  }

  tinf.version = 1;
  id_to_table(tinf.table, old->trans_table);
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

static void fprint_fp(FILE *fh, double val) {
  int digits = val == rint(val) ? 0 : DBL_DIG;
  char fmt[16];
  snprintf(fmt, sizeof(fmt), "%%.%dg", digits);
  fprintf(fh, fmt, val);
}

void fprint_fps(FILE *fh, const double *vals, size_t n) {
  int allzero = 1;
  for (size_t i = 0; i < n; i++) {
    if (vals[i] != 0.0) {
      allzero = 0;
      break;
    }
  }
  if (allzero) {
    fprintf(fh, "0");
    return;
  }
  for (size_t i = 0; i < n; i++) {
    fprint_fp(fh, vals[i]);
    if (i < n - 1) {
      fprintf(fh, ", ");
    }
  }
}

void dump_training_as_c_literal(const struct _training *tinf, FILE *fh) {
  fprintf(fh, "{\n");
  fprintf(fh, "  .version = %u,\n", tinf->version);
  fprintf(fh, "  .trans_table = %u,\n", tinf->trans_table);
  fprintf(fh, "  .uses_sd = %u,\n", tinf->uses_sd);
  fprintf(fh, "  .padding0 = %u,\n", tinf->padding0);

  fprintf(fh, "  .table = \"");
  for (int i = 0; i < 64; i++) {
    unsigned char c = tinf->table[i];
    if (c == '\"' || c == '\\') {
      fprintf(fh, "\\%c", c); // escape quotes and backslashes
    } else if (c >= 0x20 && c < 0x7f) {
      fputc(c, fh);
    } else {
      fprintf(fh, "\\x%02x", c);
    }
  }
  fprintf(fh, "\",\n");

  fprintf(fh, "  .gc = "); fprint_fp(fh, tinf->gc); fprintf(fh, ",\n");

  fprintf(fh, "  .bias = { ");
  fprint_fps(fh, tinf->bias, 3);
  fprintf(fh, " },\n");

  fprintf(fh, "  .st_wt = "); fprint_fp(fh, tinf->st_wt); fprintf(fh, ",\n");

  fprintf(fh, "  .type_wt = { ");
  fprint_fps(fh, tinf->type_wt, 3);
  fprintf(fh, " },\n");

  fprintf(fh, "  .rbs_wt = { ");
  fprint_fps(fh, tinf->rbs_wt, 28);
  fprintf(fh, " },\n");

  fprintf(fh, "  .ups_comp = {\n");
  for (int i = 0; i < 32; i++) {
    fprintf(fh, "    { ");
    fprint_fps(fh, tinf->ups_comp[i], 4);
    fprintf(fh, " }%s%c", (i < 31) ? "," : "", (i % 8 == 7) ? '\n' : ' ');
  }
  fprintf(fh, "  },\n");

  fprintf(fh, "  .mot_wt = {\n");
  for (int i = 0; i < 4; i++) {
    fprintf(fh, "    { ");
    fprint_fps(fh, tinf->mot_wt[i], NMOTIF);
    fprintf(fh, " }%s\n", (i < 3) ? "," : "");
  }
  fprintf(fh, "  },\n");

  fprintf(fh, "  .no_mot = "); fprint_fp(fh, tinf->no_mot); fprintf(fh, ",\n");

  fprintf(fh, "  .gene_dc = { ");
  fprint_fps(fh, tinf->gene_dc, 4096);
  fprintf(fh, " }\n");
  fprintf(fh, "}");
}

#ifdef TRAINING_UTIL
#include <stdlib.h>
#include "training_baked.h"
int main(int argc, char *argv[]) {
  printf("include \"training_baked.h\"\n"); 
  printf("struct _training meta_pretrain[NMETA] ={\n");
  for (int i = 0; i < NMETA; i++) {
    v0_to_v1(&meta_pretrain[i]);
    dump_training_as_c_literal((struct _training*)&meta_pretrain[i], stdout);
    if (i < NMETA - 1) {
      printf(",\n");
    } else {
      printf("\n");
    }
  }
  printf("};\n");
}
#endif
