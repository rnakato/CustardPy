#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <zlib.h>

#define max(a, b) ((a) > (b)) ? (a) :(b)

#define WINSIZE_DEFAULT 10000
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
  char *inputfile = argv[1];
  int zipped=0;
  if(strstr(inputfile, ".gz")) zipped=1;

  int winsize = WINSIZE_DEFAULT;
  if(argc>=3 && atoi(argv[2]) > 0) {
    winsize = atoi(argv[2]);
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
  
  //  while ((fgets(str, STR_LEN, IN))!=NULL) { 
  //   if(str[0]=='\n') continue;
  while(1){
    char *c=NULL;
    if (zipped) c = gzgets(gzIN, str, STR_LEN);
    else        c = fgets(str, STR_LEN, IN);
    if(!c) break;
    
    if(str[0]=='\n') continue;
    int nclm = ParseLine(str, clm);
    if(nclm < 6) continue;
    if(strcmp(clm[1].str, clm[4].str)) continue;
    
    int start = atoi(clm[2].str);
    int end = atoi(clm[5].str);
    int length = end - start;
    if(length >= winsize*ARRAYNUM) continue;
    array[length/winsize]++;
  }

  for (i=0; i<ARRAYNUM; i++) {
    printf("%d - %d | %d\n", winsize*i, winsize*(i+1)-1, array[i]);
  }

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
