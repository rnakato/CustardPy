#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <zlib.h>

#define max(a, b) ((a) > (b)) ? (a) :(b)
#define min(a, b) ((a) < (b)) ? (a) :(b)

#define STR_LEN 10000
#define ELEM_NUM 1000
#define STRUCT_ARRAY_MAX 10000
#define WOBBLE1 4
#define WOBBLE2 4

typedef enum{
  FILE_MODE_READ,
  FILE_MODE_WRITE,
  FILE_MODE_WB,
  FILE_MODE_A
} File_Mode;

typedef struct {
  char str[STR_LEN];
} Elem;

typedef struct {
  int pos1, pos2;
  int x, y;
  char str[STR_LEN];
  char tile[3*STR_LEN+10];
} Line;


void check_optimal(int count_dup, Line *read, int *dups, int *optdups, FILE *OUT_nodupname, FILE *OUT_optname, FILE *OUT_dupname);
void scan_dupreads(int count_dup, Line *read, int *dups, int *optdups, FILE *OUT_nodupname, FILE *OUT_optname, FILE *OUT_dupname);
void add_read(Line *read, int count_dup, char *str, int pos1, int pos2, Elem *clm);
int ParseLine(char *str, Elem clm[]);
int ParseLine_arbit(char *str, Elem clm[], char token);
void *my_calloc(size_t n, size_t s, char *name);
void *my_realloc(void *p, size_t s, char *name);
FILE *my_fopen(char *filename, File_Mode mode);

int tooclose(Line *read, int i, int j) {
  int diff1 = abs(read[i].pos1 - read[j].pos1);
  int diff2 = abs(read[i].pos2 - read[j].pos2);
  if (diff1 <= WOBBLE1 && diff2 <= WOBBLE2) {
      return 1;
  }
/*  if ((diff1 <= WOBBLE1 && diff2 <= WOBBLE2) || (diff2 <= WOBBLE1 && diff1 <= WOBBLE2)) {
      return 1;
  }*/
  return 0;
}

int optcheck(Line *read, int i, int j) {
  if (!strcmp(read[i].tile, read[j].tile)
      && abs(read[i].x - read[j].x) < 50
      && abs(read[i].y - read[j].y) < 50) {
    return 1;
  }
  return 0;
}


int main(int argc, char *argv[])
{
  FILE *IN=NULL;
  struct gzFile_s *gzIN=NULL;

  if(argc<=2) {
    printf("Usage: Juicer_remove_duplicate <merged_sort> <output dir>\n");
    printf("   <merged_sort>: Input file (merged_sort.txt(.gz))\n");
    printf("   <output dir>: Output directory)\n");
    exit(0);
  }

  char *inputfile = argv[1];
  int zipped=0;
  if (strstr(inputfile, ".gz")) zipped=1;

  char *str = (char *)my_calloc(STR_LEN, sizeof(char), "str");
  Elem *clm = (Elem *)my_calloc(ELEM_NUM, sizeof(Elem), "clm");

  if (zipped) {
    if ((gzIN = gzopen(inputfile, "r"))==NULL) {
      fprintf(stderr,"[E] Cannot open <%s>.\n", inputfile);
      exit(1);
    }
  } else {
    if ((IN = fopen(inputfile, "r")) == NULL) {
      fprintf(stderr,"[E] cannot open %s.\n", inputfile);
      exit(1);
    }
  }

  int count_dup=0;
  int nummax=STRUCT_ARRAY_MAX;
  Line *read = (Line *)my_calloc(nummax, sizeof(Line), "read");
/*  for(int i=0; i<nummax; i++) {
    read[i].str = (char *)my_calloc(STR_LEN, sizeof(char), "read");
    read[i].tile = (char *)my_calloc(STR_LEN, sizeof(char), "read");
  }*/

  int nummax_dups=STRUCT_ARRAY_MAX;
  int *dups = (int *)my_calloc(nummax_dups, sizeof(int), "dups");
  int *optdups = (int *)my_calloc(nummax_dups, sizeof(int), "optdups");

  char *chr1 = (char *)my_calloc(1000, sizeof(char), "chr1");
  char *chr2 = (char *)my_calloc(1000, sizeof(char), "chr2");
  char *pre_chr1 = (char *)my_calloc(1000, sizeof(char), "pre_chr1");
  char *pre_chr2 = (char *)my_calloc(1000, sizeof(char), "pre_chr2");

  int pre_strand1 = 0;
  int pre_pos1 = 0;
  int pre_index1 = 0;
  int pre_strand2 = 0;
  int pre_index2 = 0;

	// names of output files
  char *outdir = argv[2];
  char *dupname = (char *)my_calloc(10000, sizeof(char), "dupname");
  char *nodupname = (char *)my_calloc(10000, sizeof(char), "nodupname");
  char *optname = (char *)my_calloc(10000, sizeof(char), "optname");
  sprintf(dupname, "%s/dups.txt", outdir);
  sprintf(nodupname, "%s/merged_nodups.txt", outdir);
  sprintf(optname, "%s/optdups.txt", outdir);
  FILE *OUT_dupname = my_fopen(dupname, FILE_MODE_WRITE);
  FILE *OUT_nodupname = my_fopen(nodupname, FILE_MODE_WRITE);
  FILE *OUT_optname = my_fopen(optname, FILE_MODE_WRITE);

  while(1){
    char *c=NULL;
    if (zipped) c = gzgets(gzIN, str, STR_LEN);
    else        c = fgets(str, STR_LEN, IN);
    if(!c) break;

    if(str[0]=='\n') continue;
    int nclm = ParseLine(str, clm);

    int strand1 = atoi(clm[0].str);
    strcpy(chr1, clm[1].str);
    int pos1 = atoi(clm[2].str);
    int index1 = atoi(clm[3].str);

    int strand2 = atoi(clm[4].str);
    strcpy(chr2, clm[5].str);
    int pos2 = atoi(clm[6].str);
    int index2 = atoi(clm[7].str);

/*    printf("strand1: %d, chr1: %s, pos1: %d, index1: %d, strand2: %d, chr2: %s, pos2: %d, index2: %d\n",
    strand1, chr1, pos1, index1, strand2, chr2, pos2, index2);
    printf("pre_strand1: %d, pre_chr1: %s, pre_pos1: %d, pre_index1: %d, pre_strand2: %d, pre_chr2: %s, pre_index2: %d\n",
    pre_strand1, pre_chr1, pre_pos1, pre_index1, pre_strand2, pre_chr2, pre_index2);
*/
    // Duplicate
    if (strand1 == pre_strand1
        && strand2 == pre_strand2
        && index1  == pre_index1
        && index2  == pre_index2
        && abs(pos1 - pre_pos1) <= WOBBLE1
        && !strcmp(chr1, pre_chr1)
        && !strcmp(chr2, pre_chr2)) {
      add_read(read, count_dup, str, pos1, pos2, clm);
      count_dup++;
      if (count_dup >= nummax) {
        nummax += STRUCT_ARRAY_MAX;
        read = (Line *)my_realloc(read, nummax*sizeof(Line), "read");
      }
      if (count_dup >= nummax_dups) {
        nummax_dups += STRUCT_ARRAY_MAX;
        dups = (int *)my_realloc(dups, nummax_dups*sizeof(int), "dups");
        optdups = (int *)my_realloc(optdups, nummax_dups*sizeof(int), "optdups");
      }
    } else {
      scan_dupreads(count_dup, read, dups, optdups, OUT_nodupname, OUT_optname, OUT_dupname);

      // reset all the potential duplicate array variables
      count_dup = 0;
      strcpy(read[count_dup].str, str);
      add_read(read, count_dup, str, pos1, pos2, clm);
      count_dup++;
    }

    pre_strand1 = strand1;
    pre_pos1    = pos1;
    pre_index1  = index1;
    pre_strand2 = strand2;
    pre_index2  = index2;
    strcpy(pre_chr1, chr1);
    strcpy(pre_chr2, chr2);
  }

	if (count_dup == 1) {
    fprintf(OUT_nodupname, "%s", read[0].str);
	} else if (count_dup > 2) {
    check_optimal(count_dup, read, dups, optdups, OUT_nodupname, OUT_optname, OUT_dupname);
  }

  if (zipped) gzclose(gzIN);
  else fclose(IN);
  fclose(OUT_dupname);
  fclose(OUT_nodupname);
  fclose(OUT_optname);

  free(read);
  free(dups);
  free(optdups);
  free(str);
  free(clm);

  return 0;
}


void check_optimal(int count_dup, Line *read, int *dups, int *optdups, FILE *OUT_nodupname, FILE *OUT_optname, FILE *OUT_dupname)
{
  for (int j=0; j<count_dup; j++) {
    dups[j] = 0;
    optdups[j] = 0;
  }
  for (int j=0; j<count_dup; j++) {
    // (no daisy-chaining)
    if (dups[j] || optdups[j]) continue;

    // if reads are not already marked duplicate
    for (int k=j+1; k<count_dup; k++) {
      if (tooclose(read, j, k)) {
        if (optcheck(read, j, k)) {
          optdups[k] = 1;
        } else dups[k] = 1;
      }
      if (abs(read[j].pos1 - read[k].pos1) > WOBBLE1) break;
    }
  }

  // print dups out to dup file, non-dups out to non-dup file
  for (int j=0; j<count_dup; j++) {
    if (dups[j])         fprintf(OUT_dupname, "%s", read[j].str);
    else if (optdups[j]) fprintf(OUT_optname, "%s", read[j].str);
    else                 fprintf(OUT_nodupname, "%s", read[j].str);
  }
}

void scan_dupreads(int count_dup, Line *read, int *dups, int *optdups, FILE *OUT_nodupname, FILE *OUT_optname, FILE *OUT_dupname)
{

  if (count_dup==1) {
    fprintf(OUT_nodupname, "%s", read[0].str);
  } else if (count_dup==2) {
/*  printf("2 count_dup: %d %s", count_dup, read[1].str);
  printf("2 %s %s\n", read[0].tile, read[1].tile);
  printf("2 %d %d\n", read[0].x, read[1].x);
  printf("2 %d %d\n", read[0].y, read[1].y);*/
    fprintf(OUT_nodupname, "%s", read[0].str);

    if (tooclose(read, 0, 1)) {
      if (optcheck(read, 0, 1)) {
        fprintf(OUT_optname, "%s", read[1].str);
      }
      else fprintf(OUT_dupname, "%s", read[1].str);
    }
    else {
//  printf("3 count_dup: %d\n", count_dup);
      fprintf(OUT_nodupname, "%s", read[1].str);
    }
  } else if (count_dup > 2) {
    check_optimal(count_dup, read, dups, optdups, OUT_nodupname, OUT_optname, OUT_dupname);
  }
}

void add_read(Line *read, int count_dup, char *str, int pos1, int pos2, Elem *clm) {
  strcpy(read[count_dup].str, str);
  read[count_dup].pos1 = pos1;
  read[count_dup].pos2 = pos2;

//  printf("test1, count_dup %d\n", count_dup);
  Elem clm2[100];
  int nclm2 = ParseLine_arbit(clm[14].str, clm2, ':');
// printf("nclm2=%d, strlen=%ld\n", nclm2, strlen(clm[14].str));
  if (nclm2 > 100) {
    printf("Error: too many clm2.\n");
    exit(1);
  }
  if (nclm2 > 1) {
    sprintf(read[count_dup].tile, "%s%s%s", clm2[2].str, clm2[3].str, clm2[4].str);
    read[count_dup].x = atoi(clm2[5].str);
    char yStr[20];
    sscanf(clm2[6].str, "%[^'/']", yStr);
    read[count_dup].y = atoi(yStr);
  } else {
    sprintf(read[count_dup].tile, "%d", count_dup);
  }
//  printf("tile: %s, x: %d, y: %d, str: %s, pos1: %d, pos2: %d count_dup %d\n",read[count_dup].tile,
 // read[count_dup].x, read[count_dup].y, read[count_dup].str, read[count_dup].pos1, read[count_dup].pos2, count_dup);
}

int ParseLine(char *str, Elem clm[])
{
  int i, j=0, num=0, len=strlen(str);
  char *strtemp = (char *)my_calloc(len+1, sizeof(char), "strtemp");
  for(i=0; i<len; i++){
    if(str[i]=='\0' || str[i]=='\n'){
      strtemp[j]='\0';
      strcpy(clm[num].str, strtemp);
      return ++num;
    }
    if(str[i]==' '){
      strtemp[j]='\0';
      strcpy(clm[num].str, strtemp);
      num++;
      if(num >= ELEM_NUM){
        fprintf(stderr, "error: too many columns: %s", str);
        exit(0);
      }
      j=0;
    }else{
      strtemp[j]=str[i];
      j++;
    }
  }
  free(strtemp);
  return num;
}

int ParseLine_arbit(char *str, Elem clm[], char token){
  int i, j=0, num=0, len=strlen(str);
  char *strtemp = (char *)my_calloc(len+1, sizeof(char), "strtemp");
  for(i=0; i<len; i++){
    if(str[i]=='\0'){
      strtemp[j]='\0';
      strcpy(clm[num].str, strtemp);
      free(strtemp);
      return ++num;
    }
    if(str[i]==token){
      strtemp[j]='\0';
      strcpy(clm[num].str, strtemp);
      num++;
      j=0;
    }else{
      strtemp[j]=str[i];
      j++;
    }
  }
  free(strtemp);
  return num;
}

void *my_calloc(size_t n, size_t s, char *name)
{
  void *p;
  p = calloc(n,s);
  if(!p){
    fprintf(stderr,"[E]failed calloc: %s\n", name);
    exit(1);
  }
  return p;
}

void *my_realloc(void *p, size_t s, char *name)
{
  p = realloc(p,s);
  if(!p){
    fprintf(stderr,"[E]failed realloc: %s\n", name);
    exit(1);
  }
  return p;
}

FILE *my_fopen(char *filename, File_Mode mode)
{
  FILE *IN=NULL;
  switch(mode){
  case FILE_MODE_READ:  IN = fopen(filename, "r"); break;
  case FILE_MODE_WRITE: IN = fopen(filename, "w"); break;
  case FILE_MODE_WB:    IN = fopen(filename, "wb"); break;
  case FILE_MODE_A:     IN = fopen(filename, "a+"); break;
  default: fprintf(stderr,"[E] Invalid File_Mode: %d.\n", mode); exit(1);
  }
  if(!IN){
    fprintf(stderr,"[E] Cannot open <%s>.\n", filename);
    exit(1);
  }
  return IN;
}