/*******************************************************************************
  PRODIGAL (PROkaryotic DynamIc Programming Genefinding ALgorithm)
  Copyright (C) 2007-2016 University of Tennessee / UT-Battelle

  table.h: routines and data structures for parsing NCBI-style 64-character
  genetic code strings

  Code Author:  Mingye Wang

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

#include "table.h"
#include <string.h> /* memcpy */
#include <stdio.h>

/* Part 1. Convert the NCBI amino-acid order from/to trinuc() order. */
/* 1.2 N_*: ncbi codon order: tt{t,c,a,g}, tc{t,c,a,g}, ... */
#define N_T 0
#define N_C 1
#define N_A 2
#define N_G 3
#define N_3(x, y, z) ((x) << 4 | (y) << 2 | (z))

const char ncbi_to_trinuc[64] = {
  P_3(P_T, P_T, P_T), P_3(P_T, P_T, P_C), P_3(P_T, P_T, P_A), P_3(P_T, P_T, P_G),
  P_3(P_T, P_C, P_T), P_3(P_T, P_C, P_C), P_3(P_T, P_C, P_A), P_3(P_T, P_C, P_G),
  P_3(P_T, P_A, P_T), P_3(P_T, P_A, P_C), P_3(P_T, P_A, P_A), P_3(P_T, P_A, P_G),
  P_3(P_T, P_G, P_T), P_3(P_T, P_G, P_C), P_3(P_T, P_G, P_A), P_3(P_T, P_G, P_G),
  P_3(P_C, P_T, P_T), P_3(P_C, P_T, P_C), P_3(P_C, P_T, P_A), P_3(P_C, P_T, P_G),
  P_3(P_C, P_C, P_T), P_3(P_C, P_C, P_C), P_3(P_C, P_C, P_A), P_3(P_C, P_C, P_G),
  P_3(P_C, P_A, P_T), P_3(P_C, P_A, P_C), P_3(P_C, P_A, P_A), P_3(P_C, P_A, P_G),
  P_3(P_C, P_G, P_T), P_3(P_C, P_G, P_C), P_3(P_C, P_G, P_A), P_3(P_C, P_G, P_G),
  P_3(P_A, P_T, P_T), P_3(P_A, P_T, P_C), P_3(P_A, P_T, P_A), P_3(P_A, P_T, P_G),
  P_3(P_A, P_C, P_T), P_3(P_A, P_C, P_C), P_3(P_A, P_C, P_A), P_3(P_A, P_C, P_G),
  P_3(P_A, P_A, P_T), P_3(P_A, P_A, P_C), P_3(P_A, P_A, P_A), P_3(P_A, P_A, P_G),
  P_3(P_A, P_G, P_T), P_3(P_A, P_G, P_C), P_3(P_A, P_G, P_A), P_3(P_A, P_G, P_G),
  P_3(P_G, P_T, P_T), P_3(P_G, P_T, P_C), P_3(P_G, P_T, P_A), P_3(P_G, P_T, P_G),
  P_3(P_G, P_C, P_T), P_3(P_G, P_C, P_C), P_3(P_G, P_C, P_A), P_3(P_G, P_C, P_G),
  P_3(P_G, P_A, P_T), P_3(P_G, P_A, P_C), P_3(P_G, P_A, P_A), P_3(P_G, P_A, P_G),
  P_3(P_G, P_G, P_T), P_3(P_G, P_G, P_C), P_3(P_G, P_G, P_A), P_3(P_G, P_G, P_G)
};

/* 1.1 P_*: trinuc()/is_a() order: {a,g,c,t}aa, {a,g,c,t}ga, ... */
const char trinuc_to_ncbi[64] = {
  N_3(N_A, N_A, N_A), N_3(N_G, N_A, N_A), N_3(N_C, N_A, N_A), N_3(N_T, N_A, N_A),
  N_3(N_A, N_G, N_A), N_3(N_G, N_G, N_A), N_3(N_C, N_G, N_A), N_3(N_T, N_G, N_A),
  N_3(N_A, N_C, N_A), N_3(N_G, N_C, N_A), N_3(N_C, N_C, N_A), N_3(N_T, N_C, N_A),
  N_3(N_A, N_T, N_A), N_3(N_G, N_T, N_A), N_3(N_C, N_T, N_A), N_3(N_T, N_T, N_A),
  N_3(N_A, N_A, N_G), N_3(N_G, N_A, N_G), N_3(N_C, N_A, N_G), N_3(N_T, N_A, N_G),
  N_3(N_A, N_G, N_G), N_3(N_G, N_G, N_G), N_3(N_C, N_G, N_G), N_3(N_T, N_G, N_G),
  N_3(N_A, N_C, N_G), N_3(N_G, N_C, N_G), N_3(N_C, N_C, N_G), N_3(N_T, N_C, N_G),
  N_3(N_A, N_T, N_G), N_3(N_G, N_T, N_G), N_3(N_C, N_T, N_G), N_3(N_T, N_T, N_G),
  N_3(N_A, N_A, N_C), N_3(N_G, N_A, N_C), N_3(N_C, N_A, N_C), N_3(N_T, N_A, N_C),
  N_3(N_A, N_G, N_C), N_3(N_G, N_G, N_C), N_3(N_C, N_G, N_C), N_3(N_T, N_G, N_C),
  N_3(N_A, N_C, N_C), N_3(N_G, N_C, N_C), N_3(N_C, N_C, N_C), N_3(N_T, N_C, N_C),
  N_3(N_A, N_T, N_C), N_3(N_G, N_T, N_C), N_3(N_C, N_T, N_C), N_3(N_T, N_T, N_C),
  N_3(N_A, N_A, N_T), N_3(N_G, N_A, N_T), N_3(N_C, N_A, N_T), N_3(N_T, N_A, N_T),
  N_3(N_A, N_G, N_T), N_3(N_G, N_G, N_T), N_3(N_C, N_G, N_T), N_3(N_T, N_G, N_T),
  N_3(N_A, N_C, N_T), N_3(N_G, N_C, N_T), N_3(N_C, N_C, N_T), N_3(N_T, N_C, N_T),
  N_3(N_A, N_T, N_T), N_3(N_G, N_T, N_T), N_3(N_C, N_T, N_T), N_3(N_T, N_T, N_T)
};

/* We only do ATG, GTG, TTG start, and it's only ever Met */
const char *sncbieaa_template = 
  "---M---------------M---------------M----------------------------";
inline static char sncbieaa_get(const char sncbieaa[65], int i) {
  if (sncbieaa[i] != '*' && sncbieaa[i] != '-')
    return sncbieaa_template[i];
  else
    return sncbieaa[i];
}

void eaa_to_table(char out[64], const char ncbieaa[65], const char sncbieaa[65]) {
  char temp_table[64];
  memcpy(temp_table, ncbieaa, 64);
  for (int i = 0; i < 64; i++) {
    /* The content of the table is the same as ncbieaa, except we encode the data
       from sncbieaa in bits. */
    if (sncbieaa[i] == '*' && ncbieaa[i] != '*')
        temp_table[i] |= 0x80; /* Only set the special bit for optional stops */
    else if (sncbieaa_get(sncbieaa, i) != '-')
        temp_table[i] |= 0x20; 

    out[(unsigned) ncbi_to_trinuc[i]] = temp_table[i];
  }
}

void table_to_eaa(const char table[64], char ncbieaa[65], char sncbieaa[65]) {
  for (int i = 0; i < 64; i++) {
    unsigned trinuc = trinuc_to_ncbi[i];
    char amino = table[trinuc];
    ncbieaa[i] = amino &~ 0xA0;
    if (amino & 0x80 || ncbieaa[i] == '*')
      sncbieaa[i] = '*';
    else if (amino & 0x20)
      sncbieaa[i] = 'M';
    else
      sncbieaa[i] = '-';
  }
  ncbieaa[64] = '\0';
  sncbieaa[64] = '\0';
}

/* Part 2. Predefined tables, which is a mix of NCBI and other data.
   Someday there will be a script to do this automatically from gc.prt
    */
const char predefined_tables[34][2][65] = {
  /* 0. Unused, fill with Veronika Kivenson "34" (11 + Pyl) */
  {
    "FFLLSSSSYY*OCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
    "---M------*---*----M---------------M----------------------------"
  },
  /* 1. Standard */
  {
    "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
    "---M------**--*----M---------------M----------------------------"
  },
  /* 2. Vertebrate Mito */
  {
    "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSS**VVVVAAAADDEEGGGG",
    "----------**--------------------MMMM----------**---M------------"
  },
  /* 3. Yeast Mito */
  {
    "FFLLSSSSYY**CCWWTTTTPPPPHHQQRRRRIIMMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
    "----------**----------------------MM---------------M------------",
  },
  /* 4. Mold Mito; Protozoan Mito; Coelenterate Mito; Mycoplasma; Spiroplasma */
  {
    "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
    "--MM------**-------M------------MMMM---------------M------------"
  },
  /* 5. Invertebrate Mito */
  {
    "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSSSSVVVVAAAADDEEGGGG",
    "---M------**--------------------MMMM---------------M------------"
  },
  /* 6. Ciliate Nuc; Dasycladacean Nuc; Hexamita Nuc */
  {
    "FFLLSSSSYYQQCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
    "--------------*--------------------M----------------------------"
  },
  /* 7. Unused, fill with 11 + Sec for fun (SECIS says otherwise) */
  {
    "FFLLSSSSYY**CCUWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
    "---M------**-------M---------------M----------------------------"
  },
  /* 8. Unused, fill with 11 + Sec + Pyl for fun */
  {
    "FFLLSSSSYY*OCCUWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
    "---M------*--------M---------------M----------------------------"
  },
  /* 9. Echinoderm Mito; Flatworm Mito */
  {
    "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNNKSSSSVVVVAAAADDEEGGGG",
    "----------**-----------------------M---------------M------------"
  },
  /* 10. Euplotid Nuc */
  {
    "FFLLSSSSYY**CCCWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
    "----------**-----------------------M----------------------------"
  },
  /* 11. Bacterial, Archaeal and Plant Plastid (OUR DEFAULT) */
  {
    "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
    "---M------**--*----M---------------M----------------------------"
  },
  /* 12. Alt Yeast Nuc */
  {
    "FFLLSSSSYY**CC*WLLLSPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
    "----------**--*----M---------------M----------------------------"
  },
  /* 13. Ascidian Mito */
  {
    "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSSGGVVVVAAAADDEEGGGG",
    "---M------**----------------------MM---------------M------------"
  },
  /* 14. Alt Flatworm Mito */
  {
    "FFLLSSSSYYY*CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNNKSSSSVVVVAAAADDEEGGGG",
    "-----------*-----------------------M----------------------------"
  },
  /* 15. Blepharisma Macronuc */
  {
    "FFLLSSSSYY*QCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
    "----------*---*--------------------M----------------------------"
  },
  /* 16. Chlorophycean Mito */
  {
    "FFLLSSSSYY*LCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
    "----------*---*--------------------M----------------------------"
  },
  /* 17. Unused, fill with Eddy Enterosoma/UBA4682 */
  {
    "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRMVVVVAAAADDEEGGGG",
    "---M------**--*----M---------------M----------------------------"
  },
  /* 18. Unused, fill with Eddy Peptacetobacter */
  {
    "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRQIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
    "---M------**--*----M---------------M----------------------------"
  },
  /* 19. Unused, fill with Eddy Anaerococcus, Onthovivens/UBA4855 */
  {
    "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRWIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
    "---M------**--*----M---------------M----------------------------"
  },
  /* 20. Unused, fill with Eddy Absconditabacterales */
  {
    "FFLLSSSSYY**CCGWLLLLPPPPHHQQRRWWIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
    "---M------**-----------------------M---------------M------------"
  },
  /* 21. Trematode Mito */
  {
    "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNNKSSSSVVVVAAAADDEEGGGG",
    "----------**-----------------------M---------------M------------"
  },
  /* 22. Scenedesmus obliquus Mito */
  {
    "FFLLSS*SYY*LCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
    "------*---*---*--------------------M----------------------------"
  },
  /* 23. Thraustochytrium Mito */
  {
    "FF*LSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
    "--*-------**--*-----------------M--M---------------M------------"
  },
  /* 24. Rhabdopleuridae Mito */
  {
    "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSSKVVVVAAAADDEEGGGG",
    "---M------**-------M---------------M---------------M------------"
  },
  /* 25. Candidate Division SR1 and Gracilibacteria */
  {
    "FFLLSSSSYY**CCGWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
    "---M------**-----------------------M---------------M------------"
  },
  /* 26. Pachysolen tannophilus Nuc */
  {
    "FFLLSSSSYY**CC*WLLLAPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
    "----------**--*----M---------------M----------------------------"
  },
  /* 27. Karyorelict Nuc */
  {
    "FFLLSSSSYYQQCCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
    "--------------*--------------------M----------------------------"
  },
  /* 28. Condylostoma Nuc */
  {
    "FFLLSSSSYYQQCCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
    "----------**--*--------------------M----------------------------"
  },
  /* 29. Mesodinium Nuc */
  {
    "FFLLSSSSYYYYCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
    "--------------*--------------------M----------------------------"
  },
  /* 30. Peritrich Nuc */
  {
    "FFLLSSSSYYEECC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
    "--------------*--------------------M----------------------------"
  },
  /* 31. Blastocrithidia Nuc */
  {
    "FFLLSSSSYYEECCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
    "----------**-----------------------M----------------------------"
  },
  /* 32. Balanophoraceae Plastid */
  {
    "FFLLSSSSYY*WCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
    "---M------*---*----M---------------M----------------------------"
  },
  /* 33. Cephalodiscidae Mito */
  {
    "FFLLSSSSYYY*CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSSKVVVVAAAADDEEGGGG",
    "---M-------*-------M---------------M---------------M------------"
  }
};

int id_to_table(char out[64], int id) {
  if (id < 0 || id > MAXTABLE) {
    return -1;
  }
  eaa_to_table(out, predefined_tables[id][0], predefined_tables[id][1]);
  return 0;
}

int sncbieaacmp(const char a[65], const char b[65]) {
  for (int i = 0; i < 64; i++) {
    if (sncbieaa_get(a, i) != sncbieaa_get(b, i)) {
      return 1;
    }
  }
  return 0;
}

int table_to_id(const char table[64]) {
  char ncbieaa[65], sncbieaa[65];
  table_to_eaa(table, ncbieaa, sncbieaa);
  for (int i = 0; i <= MAXTABLE; i++) {
    if (strcmp(ncbieaa, predefined_tables[i][0]) == 0 &&
        sncbieaacmp(sncbieaa, predefined_tables[i][1]) == 0) {
      return i;
    }
  }
  return -1;
}

/* Part 3. Handling user input */
/* Unfortunately * also has 0x20 set, so we use the unset version... */
const char *alphabet = "ACDEFGHIKLMNPQRSTVWYUO\x0A";
const char *salphabet = "ACDEFGHIKLMNPQRSTVWYUO*-";
int cmdline_eaa_to_table(char table[64], const char input[130]) {
  int len = strlen(input);
  int input_is_table = 0;
  char ncbieaa[65], sncbieaa[65];
  if (len != 64 && len != 129) {
    return -1;
  }

  for (int i = 0; i < 64; i++) {
    if (strchr(alphabet, (input[i] & (~0xA0))) == NULL)
      return -2;
    if (input[i] != '*' && (input[i] & 0xA0) != 0)
      input_is_table = 1;
  }
  if (input_is_table) {
    memcpy(table, input, 64);
    return 0;
  }

  memcpy(ncbieaa, input, 64);
  ncbieaa[64] = '\0';
  if (len == 64) {
    memcpy(sncbieaa, sncbieaa_template, 64);
    sncbieaa[64] = '\0';
    for (int i = 0; i < 64; i++) {
      if (ncbieaa[i] == '*') {
        sncbieaa[i] = '*';
      }
    }
  } else {
    memcpy(sncbieaa, input + 65, 64);
    sncbieaa[64] = '\0';
    for (int i = 0; i < 64; i++) {
      if (strchr(salphabet, sncbieaa[i]) == NULL)
        return -3;
    }
  }

  eaa_to_table(table, ncbieaa, sncbieaa);
  return 0;
}