#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <zlib.h>

#define max(a, b) ((a) > (b)) ? (a) :(b)

#define WINSIZE_DEFAULT 10000
#define MAPQTHRE_DEFAULT 30
#define ARRAYNUM 25000
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
    printf("Usage: distance_vs_count.Juicer <file> <winsize> <MAPQ>\n");
    printf("       <file>:    Input file  (merged_nodups.txt.gz)\n");
    printf("       <winsize>: window size (default: %d)\n", WINSIZE_DEFAULT);
    printf("       <MAPQ>:    MAPQ threshold (default: %d)\n", MAPQTHRE_DEFAULT);
    exit(0);
  }

  char *inputfile = argv[1];
  int zipped=0;
  if(strstr(inputfile, ".gz")) zipped=1;

  int winsize = WINSIZE_DEFAULT;
  if(argc>=3 && atoi(argv[2]) > 0) {
    winsize = atoi(argv[2]);
  }
  int thre_mapq = MAPQTHRE_DEFAULT;
  if(argc>=4 && atoi(argv[3]) > 0) {
    thre_mapq = atoi(argv[3]);
  }

  char *str = (char *)my_calloc(STR_LEN, sizeof(char), "str");
  int *array = (int *)my_calloc(ARRAYNUM, sizeof(int), "array");

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

  int nread=0;
  int nread_mapq_either=0;
  int nread_mapq_both=0;

  while(1){
    char *c=NULL;
    if (zipped) c = gzgets(gzIN, str, STR_LEN);
    else        c = fgets(str, STR_LEN, IN);
    if(!c) break;

    if(str[0]=='\n') continue;
    int nclm = ParseLine(str, clm);
    if(nclm < 7) continue;
    if(strcmp(clm[1].str, clm[5].str)) continue;

    int mapq1 = atoi(clm[8].str);
    int mapq2 = atoi(clm[11].str);
    nread++;
    if(mapq1 >= thre_mapq && mapq2 >= thre_mapq) nread_mapq_both++;
    else if(mapq1 >= thre_mapq || mapq2 >= thre_mapq) nread_mapq_either++;

    int start = atoi(clm[2].str);
    int end = atoi(clm[6].str);
    int length = abs(end - start);
    if(length >= winsize*ARRAYNUM) continue;
    array[length/winsize]++;
  }

  for (i=0; i<ARRAYNUM; i++) {
    printf("%d - %d | %d\n", winsize*i, winsize*(i+1)-1, array[i]);
  }

  //  printf("nread: %d, nread_either: %d, nread_both: %d\n", nread, nread_mapq_either, nread_mapq_both);

  free(str);
  free(array);
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
