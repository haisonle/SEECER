#include <stdio.h>
#include <string.h>
#include <unistd.h>

#include <vector>
using namespace std;

char* getread(int& rid, char* annot, char* seq, FILE* f2, FILE* f3) {
  char a;
  fscanf(f2, "%c%d\n%[^\n]\n", &a, &rid, seq);
  fscanf(f3, "%c%[^\n]\n%[^\n]\n", &a, annot, seq);
  char* ret = annot;
  while (*ret != '\0' && *ret != '<') {
      ret++;
  }
  if (*ret == '<') ret++;

  return ret;
}

void parseFileNames(vector<FILE*>& files, char* fn, const char* opt) {
  char* f = fn;
  bool stop = false;
  while (1) {
    fn++;
    if (*fn == ',' || *fn == '\0') {
      if (*fn == '\0') {
	stop = true;
      }
      *fn = '\0';
      files.push_back(fopen(f, opt));
      f = fn+1;
      if (stop) {
	return;
      }
    }
  }
}

void closeFiles(vector<FILE*>& files) {
    for (int i = 0; i < (int) files.size(); ++i) {
    fclose(files[i]);
  }
}

int main(int argc, char* argv[]) {

  vector<FILE*> origs;
  vector<FILE*> preprocessed;
  vector<FILE*> output;

  FILE* corrected = fopen(argv[optind], "r");
  parseFileNames(origs, argv[optind + 1], "r");
  parseFileNames(preprocessed, argv[optind + 2], "r");
  parseFileNames(output, argv[optind + 3], "w");

  if (origs.size() == preprocessed.size() 
      && origs.size() == output.size()) {

    int n = origs.size();
  
    char id[256];
    char seq[256];
    char seqc[256];
    char annot[256];
    int rid;
    char a;

    int ridx = 0;
    do {

	for (int fi = 0; fi < n; ++fi) {
	    char* an = getread(rid, annot, seqc, preprocessed[fi], corrected);

	    int rridx = ridx;
	    bool b = false;
	    do {
	      fscanf(origs[fi], "%c%[^\n]\n%[^\n]\n", &a, id, seq);

		if (a == '+') {
		  fscanf(origs[fi], "%c%[^\n]\n%[^\n]\n", &a, id, seq);
		}

		b = (rid == rridx);
		if (!b) {
		  fprintf(output[fi], ">%s\n%s\n", id, seq);
		}
		rridx++;
	    } while (!b);

	    fprintf(output[fi], ">%s%s\n%s\n", an, id, seqc);
	}

	ridx = rid+1;

    } while (!feof(corrected));
    
    for (int fi = 0; fi < n; ++fi) {
      while (!feof(origs[fi])) {
	fscanf(origs[fi], "%c%[^\n]\n%[^\n]\n", &a, id, seq);
	if (a != '+') {
	    fprintf(output[fi], ">%s\n%s\n", id, seq);
	}

      }
    }

  } else {
    printf("Invalid Inputs!\n");
  }

  closeFiles(origs);
  closeFiles(preprocessed);
  closeFiles(output);
}
