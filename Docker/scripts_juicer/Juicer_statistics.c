#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <zlib.h>

#define max(a, b) ((a) > (b)) ? (a) :(b)
#define min(a, b) ((a) < (b)) ? (a) :(b)

#define STR_LEN 100000000
#define ELEM_NUM 1000000
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
int calculate_dist_from_restrictionsite(Chromosome chr[MAX_CHROMOSOMES], int nchr,
                int strand, const char* chrname, int position, int index, int report,
                long *five_prime_end, long *three_prime_end);
void write_stats(char *output_prefix, long nreads, long intra_fragment, long under_mapq,
                  long total_current, long n_ligation, long five_prime_end, long three_prime_end,
                  long large, long small, long ninter, long nintra, long right, long left, long inner, long outer);
void write_histogram_data(char *output_prefix, long *distCount, long *mapQ, long *mapQ_intra, long *mapQ_inter,
                          int nbin, long *bins, long innerM[][2], long outerM[][2], long rightM[][2], long leftM[][2]);
int my_bsearch(int x, long *a, int n);
char* format_with_commas(long value);
char* define_ligation(char *site);
char* get_dangling_junction(const char *ligation_junction);
int ParseLine(char *str, Elem clm[]);
void *my_calloc(size_t n, size_t s, char *name);
FILE *my_fopen(char *filename, File_Mode mode);


int main(int argc, char *argv[])
{
  FILE *IN=NULL;
  struct gzFile_s *gzIN=NULL;

  if(argc<=5) {
    printf("Usage: Juicer_statistics <merged_nodups> <site file> <enzyme> <MAPQ> <output prefix>\n");
    printf("       <merged_nodups>:    Input file  (merged_nodups.txt(.gz))\n");
    printf("       <site file>:    restriction site file  (<build>_<enzyme>.txt)\n");
    printf("       <enzyme>:  enzyme\n");
    printf("       <MAPQ>:    MAPQ threshold\n");
    printf("       <output prefix>: Prefix of the output file (e.g., 'inter_30'))\n");
    exit(0);
  }

  char *inputfile = argv[1];
  int zipped=0;
  if(strstr(inputfile, ".gz")) zipped=1;

  char *restrictionsitefile = argv[2];

  char *enzyme = argv[3];
  char *ligation = define_ligation(enzyme);
  char *dangling_junction = get_dangling_junction(ligation);

 // printf("enzyme: %s, ligation: %s, dangling_junction: %s\n",enzyme, ligation, dangling_junction);
 // exit(0);

  int thre_mapq = atoi(argv[4]);
  char *output_prefix = argv[5];

  char *str = (char *)my_calloc(STR_LEN, sizeof(char), "str");
  Elem *clm = (Elem *)my_calloc(ELEM_NUM, sizeof(Elem), "clm");

  Chromosome chr[MAX_CHROMOSOMES];
  int nchr = read_restrictionsitefile(chr, restrictionsitefile);

/*  for (int i=0; i<nchr; i++) {
    printf("%s: %ld, %d\n", chr[i].name, chr[i].pos[0], chr[i].npos);
  }*/

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

  long nreads=0;
  long intra_fragment=0;
  long under_mapq=0;
  long total_current=0;
  long n_ligation=0;

  int dangling=0;

  long nintra=0, ninter=0;
  long right=0, left=0, inner=0, outer=0;
  long very_small=0, small=0, large=0;

  long five_prime_end=0, three_prime_end=0;

  long bins[] = {10,12,15,19,23,28,35,43,53,66,81,100,123,152,187,231,285,351,433,534,658,811,1000,1233,1520,1874,2310,2848,3511,4329,5337,6579,8111,10000,12328,15199,18738,23101,28480,35112,43288,53367,65793,81113,100000,123285,151991,187382,231013,284804,351119,432876,533670,657933,811131,1000000,1232847,1519911,1873817,2310130,2848036,3511192,4328761,5336699,6579332,8111308,10000000,12328467,15199111,18738174,23101297,28480359,35111917,43287613,53366992,65793322,81113083,100000000,123284674,151991108,187381742,231012970,284803587,351119173,432876128,533669923,657933225,811130831,1000000000,1232846739,1519911083,1873817423,2310129700,2848035868,3511191734,4328761281,5336699231,6579332247,8111308308,10000000000};
  int nbin = sizeof(bins) / sizeof(bins[0]);
  long distCount[NUM_DISTCOUNT];
  for (int i=0; i<NUM_DISTCOUNT; i++) distCount[i] = 0;

  long rightM[nbin][2], leftM[nbin][2], innerM[nbin][2], outerM[nbin][2];
  for (int i=0; i<nbin; i++) {
    rightM[i][0] = leftM[i][0] = innerM[i][0] = outerM[i][0] = bins[i];
    rightM[i][1] = leftM[i][1] = innerM[i][1] = outerM[i][1] = 0;
  }
  long mapQ[MAPQ_ARRAYNUM], mapQ_inter[MAPQ_ARRAYNUM], mapQ_intra[MAPQ_ARRAYNUM];
  for (int i=0; i<MAPQ_ARRAYNUM; i++) {
    mapQ[i] = mapQ_inter[i] = mapQ_intra[i] = 0;
  }

  while(1){
    char *c=NULL;
    if (zipped) c = gzgets(gzIN, str, STR_LEN);
    else        c = fgets(str, STR_LEN, IN);
    if(!c) break;

    if(str[0]=='\n') continue;
    int nclm = ParseLine(str, clm);
    nreads++;

  // 0 strand1
  // 1 chr1
  // 2 pos1
  // 3 frag1
  // 4 strand2
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

    int strand1 = atoi(clm[0].str);
    char chr1[100];
    strcpy(chr1, clm[1].str);
    int pos1 = atoi(clm[2].str);
    int index1 = atoi(clm[3].str);

    int strand2 = atoi(clm[4].str);
    char chr2[100];
    strcpy(chr2, clm[5].str);
    int pos2 = atoi(clm[6].str);
    int index2 = atoi(clm[7].str);

    int mapq1 = atoi(clm[8].str);
    int mapq2 = atoi(clm[11].str);

    int is_hiccontact=1;
    if(!strcmp(chr1, chr2) && (index1 == index2)) {
      intra_fragment++;
      is_hiccontact=0;
    } else if (nclm > 8) {
      if(mapq1 < thre_mapq || mapq2 < thre_mapq) {
        under_mapq++;
        is_hiccontact=0;
      }
    }

    // Hi-C contact
    if(!is_hiccontact) continue;

    total_current++;
    int pos_dist = abs(pos1 - pos2);
    int hist_dist = my_bsearch(pos_dist, bins, nbin);

    // one part of read pair has unligated end
    int is_dangling=0;
    char *foundAt10 = strstr(clm[10].str, dangling_junction);
    char *foundAt13 = strstr(clm[13].str, dangling_junction);
    if (nclm > 8 && ((foundAt10 != NULL && foundAt10 == clm[10].str) || (foundAt13 != NULL && foundAt13 == clm[13].str))) {
      dangling++;
      is_dangling=1;
    }

    if(!strcmp(chr1, chr2)) { // intrachromosomal
      nintra++;

      if (strand1 == strand2) {
        if (strand1 == 0) {
          if (pos_dist >= 20000) right++;
          rightM[hist_dist][1]++;
        } else {
          if (pos_dist >= 20000) left++;
          leftM[hist_dist][1]++;
        }
      } else { // strand1 != strand2
        if (strand1 == 0) {
          if (pos1 < pos2) {
            if (pos_dist >= 20000) inner++;
            innerM[hist_dist][1]++;
          }
          else {
            if (pos_dist >= 20000) outer++;
            outerM[hist_dist][1]++;
          }
        }
        else {
          if (pos1 < pos2) {
            if (pos_dist >= 20000) outer++;
            outerM[hist_dist][1]++;
          } else {
            if (pos_dist >= 20000) inner++;
            innerM[hist_dist][1]++;
          }
        }
      }
      // intra reads less than 20 kbp apart
      if (pos_dist < 10) very_small++;
      else if (pos_dist < 20000) small++;
      else large++;

    } else {  // interchromosomal
      ninter++;
    }

    if (nclm > 8) {
      int mapq_val = min(mapq1, mapq2);
      if (mapq_val < MAPQ_ARRAYNUM) {
        mapQ[mapq_val]++;
        if (!strcmp(chr1, chr2)) mapQ_intra[mapq_val]++;
        else mapQ_inter[mapq_val]++;
      }
      // read pair contains ligation junction
      if (strstr(clm[10].str, ligation) != NULL || strstr(clm[13].str, ligation) != NULL) {
        n_ligation++;
      }
    }

    // determine distance from nearest restriction site, add to histogram
    int report = (strcmp(chr1, chr2)) || (pos_dist >= 20000);
//    printf("chr1 %s chr2 %s strand1 %d strand2 %d pos1 %d pos2 %d index1 %d index2 %d report %d\n", chr1, chr2, strand1, strand2, pos1, pos2, index1, index2, report);
    int dist = calculate_dist_from_restrictionsite(chr, nchr, strand1, chr1, pos1, index1, report, &five_prime_end, &three_prime_end);
    if (dist < NUM_DISTCOUNT) distCount[dist]++;
    dist = calculate_dist_from_restrictionsite(chr, nchr, strand2, chr2, pos2, index2, report, &five_prime_end, &three_prime_end);
    if (dist < NUM_DISTCOUNT) distCount[dist]++;

    if (is_dangling) {
      if (!strncmp(clm[10].str, dangling_junction, strlen(dangling_junction))) {
        dist = calculate_dist_from_restrictionsite(chr, nchr, strand1, chr1, pos1, index1, 1, &five_prime_end, &three_prime_end);
      }
      else { //	$record[13] =~ m/^$dangling_junction/)
        dist = calculate_dist_from_restrictionsite(chr, nchr, strand2, chr2, pos2, index2, 1, &five_prime_end, &three_prime_end);
      }
    }
  }

  /// Report the statistics
  write_stats(output_prefix, nreads, intra_fragment, under_mapq, total_current, n_ligation, five_prime_end, three_prime_end, large, small, ninter, nintra, right, left, inner, outer);

  /// Write the histogram data to a file
  write_histogram_data(output_prefix, distCount, mapQ, mapQ_intra, mapQ_inter, nbin, bins, innerM, outerM, rightM, leftM);

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

int calculate_dist_from_restrictionsite(Chromosome chr[MAX_CHROMOSOMES], int nchr,
                int strand, const char* chrname, int position, int index, int report,
                long *five_prime_end, long *three_prime_end)
{
  int chromIndex = -1;
  for (int i=0; i<nchr; i++) {
      if (strcmp(chrname, chr[i].name) == 0){
          chromIndex = i;
          break;
      }
  }

  long dist1=0, dist2=0;
  if (!index) dist1 = position;
  else {
      dist1 = labs(position - chr[chromIndex].pos[index-1]);
  }
  dist2 = labs(position - chr[chromIndex].pos[index]);

  int retval = (dist1 <= dist2) ? dist1 : dist2;

/*  printf("chromIndexname %s chranme %s position %d index %d dist1 %ld dist2 %ld retval %d\n",
  chr[chromIndex].name, chrname, position, index, dist1, dist2, retval);
  exit(0);*/

  if (retval == dist1 && report) {
      strand == 0 ? (*five_prime_end)++ : (*three_prime_end)++;
  } else if (retval == dist2 && report) {
      strand == 16 ? (*five_prime_end)++ : (*three_prime_end)++;
  }

  return retval;
}

int get_nreads_seq(char *stats_file, long *nreads_seq)
{
  if (access(stats_file, F_OK)) return 0;

  char line[1024];
  int seq = 0;
  long reads = 0;
  FILE *IN = fopen(stats_file, "r");
  if (IN == NULL) {
      perror("Error opening file");
      return 0;
  }

  while (fgets(line, sizeof(line), IN)) {
    if (!strstr(line, "Sequenced")) continue;

    char *label = strtok(line, ":");
    char *readsStr = strtok(NULL, ":");
    char *ptr;

    if (readsStr != NULL) {
      char cleanNumberStr[100] = {0};
      ptr = readsStr;

      int i=0;
      while (*ptr) {
          if (*ptr != ',' && *ptr != ' ') {
              cleanNumberStr[i++] = *ptr;
          }
          ptr++;
      }
      reads = atol(cleanNumberStr);
      seq = 1;
    }
    break;
  }
  fclose(IN);

  *nreads_seq = reads;
//  printf("Reads: %ld\n", *nreads_seq);
  return seq;
}

void  write_stats(char *output_prefix, long nreads, long intra_fragment, long under_mapq,
                  long total_current, long n_ligation, long five_prime_end, long three_prime_end,
                  long large, long small, long ninter, long nintra,
                  long right, long left, long inner, long outer)
{
  char stats_file[10000];
  sprintf(stats_file, "%s.txt", output_prefix);

  // Read the number of sequenced reads from the previous stats file
  long nreads_seq = 0;
  int isseq = get_nreads_seq(stats_file, &nreads_seq);
//  printf("Reads: %ld\n", nreads_seq);

  // Write the stats to the file
  FILE *OUT = my_fopen(stats_file, FILE_MODE_A);

  char* formatted = format_with_commas(intra_fragment);
  if(isseq) {
    fprintf(OUT, "Intra-fragment Reads: %s (%.2f%% / %.2f%%)\n", formatted, intra_fragment*100/(double)nreads_seq, intra_fragment*100/(double)nreads);
  } else {
    fprintf(OUT, "Intra-fragment Reads: %s(%.2f%%)\n", formatted, intra_fragment*100/(double)nreads);
  }
  formatted = format_with_commas(under_mapq);
  if(isseq) {
    fprintf(OUT, "Below MAPQ Threshold: %s (%.2f%% / %.2f%%)\n", formatted, under_mapq*100/(double)nreads_seq, under_mapq*100/(double)nreads);
  } else {
    fprintf(OUT, "Below MAPQ Threshold: %s(%.2f%%)\n", formatted, under_mapq*100/(double)nreads);
  }
  formatted = format_with_commas(total_current);
  if(isseq) {
    fprintf(OUT, "Hi-C Contacts: %s (%.2f%% / %.2f%%)\n", formatted, total_current*100/(double)nreads_seq, total_current*100/(double)nreads);
  } else {
    fprintf(OUT, "Hi-C Contacts: %s(%.2f%%)\n", formatted, total_current*100/(double)nreads);
  }

  formatted = format_with_commas(n_ligation);
  if(isseq) {
    fprintf(OUT, " Ligation Motif Present: %s (%.2f%% / %.2f%%)\n", formatted, n_ligation*100/(double)nreads_seq, n_ligation*100/(double)nreads);
  } else {
    fprintf(OUT, " Ligation Motif Present: %s (%.2f%%)\n", formatted, n_ligation*100/(double)nreads);
  }
//    fprintf(OUT, "5 %d 3 %d\n", five_prime_end ,three_prime_end);
  if (five_prime_end + three_prime_end > 0) {
    long sum = five_prime_end + three_prime_end;
    float f1 = three_prime_end*100/(double)sum;
    float f2 = five_prime_end*100/(double)sum;
    fprintf(OUT, " 3' Bias (Long Range): %.1f%% - %.1f%%\n", f1, f2);
  }
  else {
    fprintf(OUT, " 3' Bias (Long Range): 0%% - 0%%\n");
  }
  if (large > 0) {
    float fleft  = left*100/(double)large;
    float finner = inner*100/(double)large;
    float fouter = outer*100/(double)large;
    float fright = right*100/(double)large;
    fprintf(OUT, " Pair Type %%(L-I-O-R): %.1f%% - %.1f%% - %.1f%% - %.1f%%\n", fleft, finner, fouter, fright);
  }
  else {
    fprintf(OUT, " Pair Type %%(L-I-O-R): 0%% - 0%% - 0%% - 0%%\n");
  }

  formatted = format_with_commas(ninter);
  if(isseq) {
    fprintf(OUT, "Inter-chromosomal: %s (%.2f%% / %.2f%%)\n", formatted, ninter*100/(double)nreads_seq, ninter*100/(double)nreads);
  } else {
    fprintf(OUT, "Inter-chromosomal: %s (%.2f%%)\n", formatted, ninter*100/(double)nreads);
  }
  formatted = format_with_commas(nintra);
  if(isseq) {
    fprintf(OUT, "Intra-chromosomal: %s (%.2f%% / %.2f%%)\n", formatted, nintra*100/(double)nreads_seq, nintra*100/(double)nreads);
  } else {
    fprintf(OUT, "Intra-chromosomal: %s (%.2f%%)\n", formatted, nintra*100/(double)nreads);
  }
  formatted = format_with_commas(small);
  if(isseq) {
    fprintf(OUT, "Short Range (<20Kb): %s (%.2f%% / %.2f%%)\n", formatted, small*100/(double)nreads_seq, small*100/(double)nreads);
  } else {
    fprintf(OUT, "Short Range (<20Kb): %s (%.2f%%)\n", formatted, small*100/(double)nreads);
  }
  formatted = format_with_commas(large);
  if(isseq) {
    fprintf(OUT, "Long Range (>20Kb): %s (%.2f%% / %.2f%%)\n", formatted, large*100/(double)nreads_seq, large*100/(double)nreads);
  } else {
    fprintf(OUT, "Long Range (>20Kb): %s (%.2f%%)\n", formatted, large*100/(double)nreads);
  }
  fclose(OUT);
}

void write_histogram_data(char *output_prefix, long *distCount,
                          long *mapQ, long *mapQ_intra, long *mapQ_inter,
                          int nbin, long *bins, long innerM[][2], long outerM[][2], long rightM[][2], long leftM[][2])
{
  char histm_file[10000];
  sprintf(histm_file, "%s_hists.m", output_prefix);
  FILE *OUT = my_fopen(histm_file, FILE_MODE_WRITE);

  fprintf(OUT, "A = [\n");
  for (int i=1; i<NUM_DISTCOUNT; i++) {
    fprintf(OUT, "%ld ", distCount[i]);
  }
  fprintf(OUT, "\n];\n");

  fprintf(OUT, "B = [\n");
  for (int i=0; i<MAPQ_ARRAYNUM; i++) {
    fprintf(OUT, "%ld %ld %ld\n ", mapQ[i], mapQ_intra[i], mapQ_inter[i]);
  }
  fprintf(OUT, "\n];\n");

  fprintf(OUT, "D = [\n");
  for (int i=0; i<nbin; i++) {
//    printf("%ld %ld %ld %ld\n", innerM[i][1], outerM[i][1], rightM[i][1], leftM[i][1]);
    fprintf(OUT, "%ld %ld %ld %ld\n", innerM[i][1], outerM[i][1], rightM[i][1], leftM[i][1]);
  }
  fprintf(OUT, "\n];");
  fprintf(OUT, "x = [\n");
  for (int i=0; i<nbin; i++) fprintf(OUT, "%ld ", bins[i]);
  fprintf(OUT, "\n];\n");

  fclose(OUT);
}


// Binary search function as previously defined
int my_bsearch(int x, long *a, int n)
{
    int l = 0; // Lower bound
    int u = n - 1; // Upper bound
    int i; // Index of probe

    while (l <= u) {
      i = (l + u) / 2; // Midpoint of l and u
      if (a[i] < x) l = i + 1;
      else if (a[i] > x) u = i - 1;
      else return i; // x found at index i
    }

    return l; // x not found, return upper bound
}

char* format_with_commas(long value)
{
    // Convert the integer to a string
    char buffer[50]; // Ensure this buffer is large enough to hold the number
    sprintf(buffer, "%ld", value);

    // Calculate the length of the formatted string
    int length = strlen(buffer);
    int commaCount = (length - 1) / 3; // Number of commas to insert
    int formattedLength = length + commaCount;

    // Allocate memory for the formatted string
    char* formatted = (char*)malloc(formattedLength + 1); // +1 for the null terminator
    if (formatted == NULL) return NULL; // Check for allocation failure

    // Insert commas every three digits from the right
    int j = formattedLength; // Index for the formatted string
    formatted[j--] = '\0'; // Null-terminate the string
    for (int i = length - 1, commaPos = 0; i >= 0; --i) {
        formatted[j--] = buffer[i];
        if (++commaPos == 3 && i > 0) {
            formatted[j--] = ',';
            commaPos = 0;
        }
    }

    return formatted;
}

char* define_ligation(char *site)
{
    if (strcmp(site, "HindIII") == 0) {
        return "AAGCTAGCTT";
    } else if (strcmp(site, "DpnII") == 0 || strcmp(site, "MboI") == 0) {
        return "GATCGATC";
    } else if (strcmp(site, "NcoI") == 0) {
        return "CCATGCATGG";
    } else if (strcmp(site, "Arima") == 0) {
        return "(GAATAATC|GAATACTC|GAATAGTC|GAATATTC|GAATGATC|GACTAATC|GACTACTC|GACTAGTC|GACTATTC|GACTGATC|GAGTAATC|GAGTACTC|GAGTAGTC|GAGTATTC|GAGTGATC|GATCAATC|GATCACTC|GATCAGTC|GATCATTC|GATCGATC|GATTAATC|GATTACTC|GATTAGTC|GATTATTC|GATTGATC)";
    } else {
        if (strcmp(site, "none") != 0) { // Log message for unrecognized enzymes
            printf("%s not listed as recognized enzyme. Using %s as site file\n", site, site);
            printf("Ligation junction is undefined\n");
        }
        return "XXXX";
    }
}

char* get_dangling_junction(const char *ligation_junction)
{
    int length = strlen(ligation_junction);
    int midpoint = length / 2;

    // Allocate memory for the dangling junction. Add 1 for the null terminator.
    char *dangling_junction = (char *)malloc((length - midpoint + 1) * sizeof(char));
    if (dangling_junction == NULL) {
        // Return NULL if memory allocation fails
        return NULL;
    }

    // Copy the second half of the ligation junction to the dangling junction
    strcpy(dangling_junction, ligation_junction + midpoint);

    return dangling_junction;
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