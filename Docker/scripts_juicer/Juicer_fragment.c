#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <zlib.h>

#define max(a, b) ((a) > (b)) ? (a) :(b)
#define min(a, b) ((a) < (b)) ? (a) :(b)

#define STR_LEN 1000000
#define ELEM_NUM 10000000
#define MAPQ_ARRAYNUM 201
#define MAX_CHROMOSOMES 100
#define NUM_DISTCOUNT 2001

typedef enum{
  FILE_MODE_READ,
  FILE_MODE_WRITE,
  FILE_MODE_WB,
  FILE_MODE_A
} File_Mode;

typedef struct {
  char str[10240];
} Elem;

typedef struct {
  char name[100];
  long *pos;
  int npos;
} Chromosome;

int read_restrictionsitefile(Chromosome *chr, char *restrictionsitefile);
long my_bsearch(long pos, char *chrname, Chromosome chr[MAX_CHROMOSOMES], int nchr);
int ParseLine(char *str, Elem clm[]);
int ParseLine_tab(char *str, Elem clm[]);
void *my_calloc(size_t n, size_t s, char *name);
FILE *my_fopen(char *filename, File_Mode mode);


int main(int argc, char *argv[])
{
  FILE *IN=NULL;

  if(argc<=2) {
    printf("Usage: Juicer_fragment <norm.txt> <site file>\n");
    printf("       <norm.txt>:  Input file  (*_norm.txt, output of 'chimeric_blacklist.awk')\n");
    printf("       <site file>: Restriction site file  (<build>_<enzyme>.txt)\n");
    exit(0);
  }

  char *inputfile = argv[1];
  char *restrictionsitefile = argv[2];

  Chromosome chr[MAX_CHROMOSOMES];
  int nchr = read_restrictionsitefile(chr, restrictionsitefile);

  if ((IN = fopen(inputfile, "r")) == NULL) {
    fprintf(stderr,"[E] cannot open %s.\n", inputfile);
    exit(1);
  }

  char *str = (char *)my_calloc(STR_LEN, sizeof(char), "str");
  Elem *clm = (Elem *)my_calloc(ELEM_NUM, sizeof(Elem), "clm");

  while (fgets(str, STR_LEN, IN) != NULL) {
    if(str[0]=='\n') continue;
    int nclm = ParseLine_tab(str, clm);

    long pos1 = atoi(clm[2].str);
    long pos2 = atoi(clm[5].str);
    long index1 = my_bsearch(pos1, clm[1].str, chr, nchr);
    long index2 = my_bsearch(pos2, clm[4].str, chr, nchr);

    printf("%s %s %s %ld %s %s %s %ld ", clm[0].str, clm[1].str, clm[2].str, index1, clm[3].str, clm[4].str, clm[5].str, index2);
    for (int i=6; i<nclm; i++) printf("%s ", clm[i].str);
    printf("\n");
  }

  free(str);
  free(clm);
  return 0;
}


int read_restrictionsitefile(Chromosome *chr, char *restrictionsitefile)
{
  FILE *IN;
  if ((IN = fopen(restrictionsitefile, "r")) == NULL) {
    fprintf(stderr,"[E] cannot open %s.\n", restrictionsitefile);
    exit(1);
  }
  int nchr=0;
  char *str = (char *)my_calloc(STR_LEN, sizeof(char), "str");
  Elem *clm = (Elem *)my_calloc(ELEM_NUM, sizeof(Elem), "clm");

  while (fgets(str, STR_LEN, IN) != NULL) {
    if (str[0]=='\n') continue;

    int nclm = ParseLine(str, clm);
    strcpy(chr[nchr].name, clm[0].str);
    chr[nchr].pos = (long *)my_calloc(nclm-1, sizeof(long), "chr.pos");
    for (int i=1; i<nclm; i++) {
      chr[nchr].pos[i-1] = atol(clm[i].str);
    }
    chr[nchr].npos = nclm-1;

//   printf("%s, nclm%d\n", chr[nchr].name,  chr[nchr].npos);
    nchr++;
    if (nchr >= MAX_CHROMOSOMES) {
      fprintf(stderr, "error: too many chromosomes > %d\n", MAX_CHROMOSOMES);
      exit(0);
    }
  }
  fclose(IN);

  free(str);
  free(clm);

  return nchr;
}

long my_bsearch(long pos, char *chrname, Chromosome chr[MAX_CHROMOSOMES], int nchr)
{
  int chromIndex = -1;
  for (int i=0; i<nchr; i++) {
      if (strcmp(chrname, chr[i].name) == 0){
          chromIndex = i;
          break;
      }
  }

  long l = 0; // Lower bound
  long u = chr[chromIndex].npos; // Upper bound
  long index; 

  while (l <= u) {
    index = (l + u) / 2;
    if      (chr[chromIndex].pos[index] < pos) l = index + 1;
    else if (chr[chromIndex].pos[index] > pos) u = index - 1;
    else return index+1; 
  }

  return l; 
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

int ParseLine_tab(char *str, Elem clm[])
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
    if(str[i]=='\t'){
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