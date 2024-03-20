#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <zlib.h>
#include <ctype.h>

#define max(a, b) ((a) > (b)) ? (a) :(b)
#define min(a, b) ((a) < (b)) ? (a) :(b)

#define STR_LEN 10000
#define ELEM_NUM 10000
#define MAPQ_ARRAYNUM 201
#define MAX_CHROMOSOMES 100
#define NUM_DISTCOUNT 2001
#define MAX_RECORD_PARTS 12
#define MAX_PART_LENGTH 256

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
  int str;
  int pos;
  int m;
  int mapped;
  int read;
  char chr[1000];
  char cigarstr[1000];
  char seq[5000];
  char name[1000];
} SAMline;

// Static variables
int tottot=-1;
int count=1;
int count_abnorm=-1, count_reg=0, count_unmapped=0, count_norm=0;

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

int ParseLine(const char *str, Elem clm[])
{
  int i, j=0, num=0, len=strlen(str);
  char *strtemp = (char *)my_calloc(len+1, sizeof(char), "str");
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

int less_than(SAMline *line, int i, int j) {
  int s1 = line[i].str;
  char *c1 = line[i].chr;
  int p1 = line[i].pos;
  int s2 = line[j].str;
  char *c2 = line[j].chr;
  int p2 = line[j].pos;

  // Compare chromosomes
  int cmp = strcmp(c1, c2);
  if (cmp < 0) return 1;
  if (cmp > 0) return 0;

  // Chromosomes are equal, compare strands
  if (s1 < s2) return 1;
  if (s1 > s2) return 0;

  // Strands are equal, compare positions
  if (p1 < p2) return 1;
  if (p1 > p2) return 0;

  // All are equal, doesn't matter
  return 1;
}

void print_data(SAMline *line, int read1, int read2, FILE* OUT_fname) {
  if (less_than(line, read1, read2)) {
    // ideally we'll get rid of printing out cigar string at some point
/*    printf("11 read1:%d read2:%d %d\t%s\t%d\t%d\t%s\t%d\t%d\t%s\t%s\t%d\t%s\t%s\t%s\t%s\n", read1, read2,
      line[read1].str, line[read1].chr, line[read1].pos, line[read2].str, line[read2].chr, line[read2].pos,
      line[read1].m, line[read1].cigarstr, line[read1].seq, line[read2].m, line[read2].cigarstr, line[read2].seq,
      line[read1].name, line[read2].name);*/
    fprintf(OUT_fname, "%d\t%s\t%d\t%d\t%s\t%d\t%d\t%s\t%s\t%d\t%s\t%s\t%s\t%s\n",
      line[read1].str, line[read1].chr, line[read1].pos, line[read2].str, line[read2].chr, line[read2].pos,
      line[read1].m, line[read1].cigarstr, line[read1].seq, line[read2].m, line[read2].cigarstr, line[read2].seq,
      line[read1].name, line[read2].name);
  } else {
/*    printf("22 read1:%d read2:%d %d\t%s\t%d\t%d\t%s\t%d\t%d\t%s\t%s\t%d\t%s\t%s\t%s\t%s\n", read1, read2,
      line[read1].str, line[read1].chr, line[read1].pos, line[read2].str, line[read2].chr, line[read2].pos,
      line[read1].m, line[read1].cigarstr, line[read1].seq, line[read2].m, line[read2].cigarstr, line[read2].seq,
      line[read1].name, line[read2].name);*/
    fprintf(OUT_fname, "%d\t%s\t%d\t%d\t%s\t%d\t%d\t%s\t%s\t%d\t%s\t%s\t%s\t%s\n",
      line[read2].str, line[read2].chr, line[read2].pos, line[read1].str, line[read1].chr, line[read1].pos,
      line[read2].m, line[read2].cigarstr, line[read2].seq, line[read1].m, line[read1].cigarstr, line[read1].seq,
      line[read2].name, line[read1].name);
  }
}

int calculate_adjusted_length(const char* cigar) {
    int seqlength = 0;
    int num = 0; 
    int length = strlen(cigar); 

    for (int i=0; i<length; i++) {
      if (isdigit(cigar[i])) { // 現在の文字が数字なら、numに追加
        num = num * 10 + (cigar[i] - '0');
      } else {  // 文字がCIGAR操作（M, D, N, X, =, S, H）に対応している場合
        switch (cigar[i]) {
          case 'M':
          case 'D':
          case 'N':
          case 'X':
          case '=':
            seqlength += num;
            break;
          case 'S':
          case 'H':
            // 文字列の末尾にある場合のみ加算
            if (i == length - 1) seqlength += num;
            break;
        }
        num = 0;
      }
    }
    return seqlength;
}


int safe_atoi(const char *str) {
  char *endptr;
  int value = strtol(str, &endptr, 10);
  if (*endptr) {
    fprintf(stderr, "Conversion error occurred: %s\n", str);
    exit(EXIT_FAILURE);
  }
  return value;
}

void extract_read_info(const char* record, SAMline *line, int j, Elem *clm) {
  int nclm = ParseLine(record, clm);
//  printf("extract_read_info record %s nclm %d\n", record, nclm);

  char* readname = strtok(clm[0].str, "/");
  char* nextPart = strtok(NULL, "/");

  if (nextPart) {
    line[j].read = safe_atoi(nextPart);
  } else {
    int flag = atoi(clm[1].str);
    line[j].read = (flag & 64) > 0 ? 1 : 0;
  }
}

int check_pattern(const char *cigar) {
    int length = strlen(cigar);
    if (length == 0) return 0;
    if (!isdigit(cigar[0])) return 0;

    for (int i=1; i<length; i++) {
      if (!isdigit(cigar[i])) {
        if (cigar[i] == 'S' || cigar[i] == 'H') return 1;
        else return 0;
      }
    }
    return 0;
}

int adjust_position(int strand, const char* cigar, int pos) {
  if (strand == 16) {
    int adjusted_length = calculate_adjusted_length(cigar);
    pos += adjusted_length - 1;
//    printf("3 pos %d %d\n", pos, adjusted_length);
  } else if (strand == 0 && (check_pattern(cigar))) {
    int leading_bases = 0;
    sscanf(cigar, "%d", &leading_bases);
//    printf("str %d, cigar %s, leading_bases %d\n", strand, cigar, leading_bases);
    pos -= leading_bases;
    if (pos <= 0) pos = 1;
//    printf("1 pos %d\n", pos);
  }

  return pos;
}

void extract_and_assign_data(const char *record, SAMline *line, int j, Elem *clm) {

  int nclm = ParseLine(record, clm);
//  printf("record %s nclm %d\n", record, nclm);

  int sv = atoi(clm[1].str);
  strcpy(line[j].name, clm[0].str);
  line[j].str = sv & 16;
  strcpy(line[j].chr, clm[2].str);
  line[j].pos = safe_atoi(clm[3].str);
  line[j].m = safe_atoi(clm[4].str);
  strcpy(line[j].cigarstr, clm[5].str);
  strcpy(line[j].seq, clm[9].str);

  // blacklist - if 3rd bit set (=4) it means unmapped
  line[j].mapped = (sv & 4) == 0;

//  if(!strcmp(line[j].name, "M00336:181:000000000-A29H6:1:1101:16729:11813")) printf("AAAAA name %s sv %d str %d chr %s pos %d m %d cigarstr %s seq %s\n", line[j].name, sv, line[j].str, line[j].chr, line[j].pos, line[j].m, line[j].cigarstr, line[j].seq);
  line[j].pos = adjust_position(line[j].str, line[j].cigarstr, line[j].pos);

/* if(!strcmp(line[j].name, "M00336:181:000000000-A29H6:1:1101:16729:11813")) {
    printf("AAAAA name %s sv %d str %d chr %s pos %d m %d cigarstr %s seq %s\n", line[j].name, sv, line[j].str, line[j].chr, line[j].pos, line[j].m, line[j].cigarstr, line[j].seq);
  }*/
}

void check_mapped_read1and2(SAMline *line, int read1, int read2, int count, Elem *c, int *count_norm, int *count_unmapped,
                            FILE *OUT_fname1, FILE *OUT_fname3) {
  if (line[read1].mapped && line[read2].mapped) {
    (*count_norm)++;
    print_data(line, read1, read2, OUT_fname1);
	} else {
    for (int i=1; i<=count; i++) fprintf(OUT_fname3, "%s", c[i].str);
    (*count_unmapped)++;
	}
}

long calculate_distance(SAMline *line, int i, int j) {
  if(strcmp(line[i].chr, line[j].chr)) {
    return 1000000000;
  } else {
    return abs(line[i].pos - line[j].pos);
  }
}

// Function to analyze chromatin interaction data and classify read pairs
void analyze_chromatin_interactions(SAMline *line, Elem *c,
                                  FILE *OUT_fname1, FILE *OUT_fname2, FILE *OUT_fname3,
                                  int count, int *count_unmapped, int *count_abnorm, int *count_norm) {
  long dist[50];
  int read1 = 0, read2 = 0;

  for (int i=1; i<=4; i++) {
    for (int j=i+1; j<=4; j++) {
      dist[i*10 + j] = calculate_distance(line, i, j);
    }
  }

  if ((dist[13] < 1000 && dist[24] < 1000) || (dist[14] < 1000 && dist[23] < 1000)) {
    read1 = 1; read2 = 2;
  } else if (dist[12] < 1000 && dist[34] < 1000) {
    read1 = 1; read2 = 3;
  } else {
    read1 = 0;
  }


  if (read1) {
    check_mapped_read1and2(line, read1, read2, count, c, count_norm, count_unmapped, OUT_fname1, OUT_fname3);
  } else {
    // Chimeric read with the 4 ends > 1KB apart
    for (int i=1; i<=count; i++) fprintf(OUT_fname2, "%s", c[i].str);
    (*count_abnorm)++;
  }
}

// Function to analyze distance metrics and classify read pairs
void classify_reads_based_on_distance(SAMline *line, Elem *c, FILE *OUT_fname1, FILE *OUT_fname2, FILE *OUT_fname3, int count, int* count_unmapped, int* count_abnorm, int* count_norm) {
  long dist[3];
  int read1=0, read2=0;

  dist[0] = calculate_distance(line, 1, 2);
  dist[1] = calculate_distance(line, 2, 3);
  dist[2] = calculate_distance(line, 1, 3);
  long mindist = min(dist[0], (min(dist[1], dist[2])));

/*  printf("line[1].chr %s, line[2].chr %s, line[3].chr %s, \n", line[1].chr, line[2].chr, line[3].chr);
  printf("line[1].pos %d, line[2].pos %d, line[3].pos %d\n", line[1].pos, line[2].pos, line[3].pos);
  printf("line[1].mapped %d, line[2].mapped %d, line[3].mapped %d\n", line[1].mapped, line[2].mapped, line[3].mapped);
  printf("line[1].read %d, line[2].read %d, line[3].read %d\n", line[1].read, line[2].read, line[3].read);
  printf("dist[0] %ld dist[1] %ld dist[2] %ld  min %ld\n", dist[0], dist[1], dist[2], mindist);
*/
  if (mindist < 1000) {
    if (line[1].read == line[2].read) {
      read2 = 3;
      read1 = dist[2] > dist[1] ? 1 : 2;
    } else if (line[1].read == line[3].read) {
      read2 = 2;
      read1 = dist[0] > dist[1] ? 1 : 3;
    } else if (line[2].read == line[3].read) {
      read2 = 1;
      read1 = dist[0] > dist[2] ? 2 : 3;
    } else {
      fprintf(stderr, "reads strange\n");
      exit(EXIT_FAILURE);
    }
 //   printf("read1 %d read2 %d line[read1].mapped %d, line[read2].mapped %d, chr[read1] %s chr[read2] %s\n", read1, read2, line[read1].mapped, line[read2].mapped, line[read1].chr, line[read2].chr);


/*    printf("1122 read1:%d read2:%d %d\t%s\t%d\t read2 %d\t%s\t%d\t%d\t%s\t%s\t%d\t%s\t%s\t%s\t%s\n", read1, read2,
      line[read1].str, line[read1].chr, line[read1].pos, line[read2].str, line[read2].chr, line[read2].pos,
      line[read1].m, line[read1].cigarstr, line[read1].seq, line[read2].m, line[read2].cigarstr, line[read2].seq,
      line[read1].name, line[read2].name);*/
    // Process based on mapping status
    check_mapped_read1and2(line, read1, read2, count, c, count_norm, count_unmapped, OUT_fname1, OUT_fname3);
  } else { // Process as chimeric read
    for (int i=1; i<=count; i++) fprintf(OUT_fname2, "%s", c[i].str);
    (*count_abnorm)++;
  }
}

void func(char a[1024], char prev[1024], char *str, Elem *c, SAMline *line, Elem *clm,
          FILE *OUT_fname1, FILE *OUT_fname2, FILE *OUT_fname3) {
  int i, j;
  tottot++;
  
//  if(strstr(prev, "M00336:181:000000000-A29H6:1:1101:16729:11813")) printf("M00336:181:000000000-A29H6:1:1101:16729:11813 cout %d prev %s\n", count, prev);

  if (count==3 || count==4) {
    // chimeric read
    for (j=1; j<=count; j++) {
      extract_read_info(c[j].str, line, j, clm);
      extract_and_assign_data(c[j].str, line, j, clm);
    }

    if (count == 4) {
      analyze_chromatin_interactions(line, c, OUT_fname1, OUT_fname2, OUT_fname3, count, &count_unmapped, &count_abnorm, &count_norm);
    } else { // count == 3
      classify_reads_based_on_distance(line, c,  OUT_fname1, OUT_fname2, OUT_fname3, count, &count_unmapped, &count_abnorm, &count_norm);
    }
  }
  else if (count > 3) {  // chimeric read > 3, too many to deal with
    for (i=1; i<=count; i++) fprintf(OUT_fname2, "%s", c[i].str);
    count_abnorm++;
  }
  else if (count == 2) { // code here should be same as above, but it's a "normal" read
    for (i=1, j=0; i<=count; i++, j++) {
//      printf("c[%d].str %s\n", i, c[i].str);
      extract_and_assign_data(c[i].str, line, j, clm);
    }

    if (line[0].mapped && line[1].mapped) {
      count_reg++;
      print_data(line, 0, 1, OUT_fname1);
    } else {
      for (i=1; i<=count; i++) fprintf(OUT_fname3, "%s", c[i].str);
      count_unmapped++;
    }
  }
  else if (count == 1) { // this actually shouldn't happen, but it happens with alternate aligners on occasion
    count_abnorm++;
    for (i=1; i<=count; i++) fprintf(OUT_fname2, "%s", c[i].str);
  }

  count=1;
  strcpy(prev, a);
  strcpy(c[count].str, str);
}


int main(int argc, char *argv[])
{
  if(argc<=2) {
    printf("Usage: Juicer_chimeric_blacklist <*.SAM> <output prefix>\n");
    printf("       <*.SAM>: Input SAM file\n");
    printf("       <output prefix>: Prefix of the output file)\n");
    exit(0);
  }

  char *inputfile = argv[1];
  char *out_prefix = argv[2];

  FILE *IN = my_fopen(inputfile, FILE_MODE_READ);

  char fname1[10000];
  char fname2[10000];
  char fname3[10000];
  sprintf(fname1, "%s_norm.txt", out_prefix);
  sprintf(fname2, "%s_abnorm.sam", out_prefix);
  sprintf(fname3, "%s_unmapped.sam", out_prefix);
  FILE *OUT_fname1 = my_fopen(fname1, FILE_MODE_WRITE);
  FILE *OUT_fname2 = my_fopen(fname2, FILE_MODE_WRITE);
  FILE *OUT_fname3 = my_fopen(fname3, FILE_MODE_WRITE);

  char *token=NULL;
  char a[STR_LEN],  prev[STR_LEN];
  char *str = (char *)my_calloc(STR_LEN, sizeof(char), "str");
  Elem *c = (Elem *)my_calloc(ELEM_NUM, sizeof(Elem), "c");
  SAMline line[10];
  Elem *clm = (Elem *)my_calloc(ELEM_NUM, sizeof(Elem), "clm");

  Elem *clmtemp = (Elem *)my_calloc(ELEM_NUM, sizeof(Elem), "extract_and_assign_data_clm");

  while((fgets(str, STR_LEN, IN))!=NULL) {
    if(str[0]=='\n') continue;
    if(str[0]=='@') continue;

    int nclm = ParseLine(str, clm);
    token = strtok(clm[0].str, "/");
    if (token != NULL) strcpy(a, token);
    else strcpy(a, clm[0].str);

//    printf("a  %s token %s clm[0].str %s\n", a, token, clm[0].str);

    if (!strcmp(a, prev)) {
      count++;
      strcpy(c[count].str, str);
//      printf("count %d %s\n", count, c[count].str);
      continue;
    }

    func(a, prev, str, c, line, clmtemp, OUT_fname1, OUT_fname2, OUT_fname3);
  }
  func(a, prev, str, c, line, clmtemp, OUT_fname1, OUT_fname2, OUT_fname3);

  char fname4[10010];
  sprintf(fname4, "%s.res.txt", fname1);
  FILE *OUT_fname4 = my_fopen(fname4, FILE_MODE_A);
  fprintf(OUT_fname4, "%d %d %d %d %d\n", tottot, count_unmapped, count_reg, count_norm, count_abnorm);

  fclose(IN);
  fclose(OUT_fname1);
  fclose(OUT_fname2);
  fclose(OUT_fname3);
  fclose(OUT_fname4);

//  printf("str\n");
  free(str);
//  printf("c\n");
  free(c);
//  printf("clm\n");
  free(clm);
//  printf("ckntemp\n");
  free(clmtemp);

  return 0;
}
