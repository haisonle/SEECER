#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <ctype.h>
#include <unistd.h>

char letter[] = {'A', 'T', 'G', 'C'};
#define BUFFER_LENGTH 1000

#define MIN_LENGTH 35
bool subN(char* read) {

  int nN = 0;
  int nonNsince = 0;

  int i = 0;
  while (*read != '\n') {
    if (toupper(*read) == 'N') {
      nN++;
      nonNsince = 3;
      *read = letter[rand() % 4];
    } else {
      nonNsince--;
      if (nonNsince < 0) {
	nonNsince = 0;
	nN = 0;
      }
    }

    if (nN > 3) {
      if (i > MIN_LENGTH) {
	*read = '\n';
	*(read+1) = '\0';
	return true;
      } else {
	return false;
      }
    }
    ++read;
    ++i;
  }

  return true;
}

void parseFileNames(FILE **f1, FILE **f2, char* fn) {
  char* f = fn;

  while (*f != ',' && *f != '\0') {
    f++;
  }
  if (*f == '\0') {
      fprintf(stderr, "Invalid arguments\n");
      exit(1);
  } else {
      *f = '\0';
      *f1 = fopen(fn, "r");
      *f2 = fopen(f+1, "w");
  }

}


bool GetReads(char* id, char* seq, FILE* f) {
    bool eof = fgets(id, BUFFER_LENGTH, f) != NULL;

    eof &= (fgets(seq, BUFFER_LENGTH, f) != NULL);
    
    if (eof && id[0] == '+') {
	// fastaq
	eof &= (fgets(id, BUFFER_LENGTH, f) != NULL);
	eof &= (fgets(seq, BUFFER_LENGTH, f) != NULL);
    }

    return eof;
}

int main(int argc, char* argv[]) {

  srand ( time(NULL) );
    
  FILE *f1 = NULL;
  FILE *f2 = NULL;
  FILE *f3 = NULL;
  FILE *f4 = NULL;

  if (argc == 0) {
      exit(1);
  }

  if (argc > 1) {
      parseFileNames(&f1, &f3, argv[1]);
  }
  if (argc > 2) {
      parseFileNames(&f2, &f4, argv[2]);
  }

  char read_n1[BUFFER_LENGTH];
  char read_n2[BUFFER_LENGTH];

  char read1[BUFFER_LENGTH];
  char read2[BUFFER_LENGTH];
  
  int ridx = 0;

  while (GetReads(read_n1, read1,  f1)
	 && ( !f2 || GetReads(read_n2, read2, f2) )
	 ) {

      //printf("%s", read_n1);

      if (subN(read1) && (!f2 || subN(read2)) ) {
      /*
      fprintf(f3, ">%d\n%s", ridx, read1);
      fprintf(f4, ">%d\n%s", ridx, read2);
      */
	  fprintf(f3, ">%d\n%s", ridx, read1);
	  if (f2) {
	      fprintf(f4, ">%d\n%s", ridx, read2);
	  }
      }

    ridx++;
  }

  if (f1) fclose(f1);
  if (f2) fclose(f2);
  if (f3) fclose(f3);
  if (f4) fclose(f4);

  return 0;
}
