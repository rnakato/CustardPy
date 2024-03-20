#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <zlib.h>
#include <math.h>

#define max(a, b) ((a) > (b)) ? (a) :(b)

#define MAPQTHRE_DEFAULT 30
#define ARRAYNUM 100
#define STR_LEN 10000
#define ELEM_NUM 256

typedef struct{
  char str[10240];
} Elem;

int ParseLine(char *str, Elem clm[]);
void *my_calloc(size_t n, size_t s, char *name);

int main(int argc, char *argv[])
{
  int i;
  FILE *IN;
  struct gzFile_s *gzIN=NULL;
  Elem clm[ELEM_NUM];

  if(argc<=1) {
    printf("Usage: distance_vs_count.Juicer.log <file> <MAPQ>\n");
    printf("       <file>:    Input file  (merged_nodups.txt)\n");
    printf("       <MAPQ>:    MAPQ threshold (default: %d)\n", MAPQTHRE_DEFAULT);
    exit(0);
  }

  char *inputfile = argv[1];
  int zipped=0;
  if(strstr(inputfile, ".gz")) zipped=1;

  int thre_mapq = MAPQTHRE_DEFAULT;
  if(argc>=3 && atoi(argv[2]) > 0) {
    thre_mapq = atoi(argv[2]);
  }

  char *str = (char *)my_calloc(STR_LEN, sizeof(char), "str");
  long *array = (long *)my_calloc(ARRAYNUM, sizeof(long), "array");
  unsigned long *array_maxval = (unsigned long *)my_calloc(ARRAYNUM, sizeof(unsigned long), "array_maxval");
  for(i=0; i<ARRAYNUM; i++) array_maxval[i] = pow(10, (i+1)*0.1);

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

  // 0 str1
  // 1 chr1
  // 2 pos1
  // 3 frag1
  // 4 str2
  // 5 chr2
  // 6 pos2
  // 7 frag2
  // 8 MAPQ1
  // 9 CIGAR1
  // 10 seq1
  // 11 MAPQ2
  // 12 CIGAR2
  // 13 seq2
  // 14 readname1
  // 15 readname2

  long nread=0;
//  long nread_mapq_either=0;
//  long nread_mapq_both=0;
  long max_distance=0;

  unsigned long maxval = pow(10, ARRAYNUM*0.1);
//  printf("maxval = %lu\n", maxval);

  while(1){
    char *c=NULL;
    if (zipped) c = gzgets(gzIN, str, STR_LEN);
    else        c = fgets(str, STR_LEN, IN);
    if(!c) break;

    if(str[0]=='\n') continue;
    int nclm = ParseLine(str, clm);
    if(nclm < 7) continue;
    if(strcmp(clm[1].str, clm[5].str)) continue;
    nread++;

    int mapq1 = atoi(clm[8].str);
    int mapq2 = atoi(clm[11].str);
//    if(mapq1 >= thre_mapq && mapq2 >= thre_mapq) nread_mapq_both++;
//    else if(mapq1 >= thre_mapq || mapq2 >= thre_mapq) nread_mapq_either++;
    if(mapq1 < thre_mapq || mapq2 < thre_mapq) continue;

    long start = atol(clm[2].str);
    long end = atol(clm[6].str);
    unsigned long length = abs(end - start);
    if(length >= maxval) continue;

    if (max_distance < length) max_distance = length;
//    printf("length=%lu\n", length);
    for(i=0; i<ARRAYNUM; i++) {
      if(length < array_maxval[i]){
        array[i]++;
        break;
      }
    }
  }

  for (i=15; i<ARRAYNUM; i++) {
    unsigned long s = pow(10,i*0.1);
    unsigned long e = pow(10,(i+1)*0.1)-1;
    if (s > max_distance) break;
    printf("%lu - %lu | %ld\n", s, e, array[i]);
  }

  free(str);
  free(array);
  free(array_maxval);
  return 0;
}


int ParseLine(char *str, Elem clm[])
{
  int i, j=0, num=0, len=strlen(str);
  char *strtemp = (char *)my_calloc(len, sizeof(char), "ParseLine");
  for(i=0; i<=len; i++){
    if(str[i]=='\0' || str[i]=='\n'){
      strtemp[j]='\0';
      strcpy(clm[num].str, strtemp);
      free(strtemp);
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
