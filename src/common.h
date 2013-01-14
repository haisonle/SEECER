/*
Copyright (C) 2012  Hai-Son Le (haisonle@gmail.com)

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU Affero General Public License as
published by the Free Software Foundation, either version 3 of the
License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU Affero General Public License for more details.

You should have received a copy of the GNU Affero General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/
#ifndef COMMON_H_
#define COMMON_H_

#include <omp.h>

#include <deque>
#include <iostream>

#include <stdint.h>
#include <seqan/store.h>
#include <seqan/sequence.h>

#include <vector>

#define MAX_THREADS 4

#define INDELS 2
#define TOTAL_OFFSETS 2*INDELS+1

// shoud match Dna of seqan
typedef enum dnaletters {
    DNAA ,DNAC,DNAG,DNAT, NDNA
} DnaLetters;

typedef struct _error{
    int pos;
    char a, b;
    
    void DebugPrint() const {
	std::cerr << "Pos " << pos
		  << " " << a << " <= " << b;
    }

} Error;

typedef enum hmm_hidden_state {
    STATEM = 0,
    STATED, STATEI, N_HSTATES
} HMMHiddenStates;

typedef struct _hmm_parameters {

    double Init[N_HSTATES];
    double Trans[N_HSTATES][N_HSTATES];

    double pseudo;

    double alpha;
    int m;

    int k;

    double entropy_th;
    double cluster_llh_th;
    double emit_delta;

    // options
    bool do_em;
    bool restrict_failures;
    int failure_th;
    bool reuse_reads;
    bool do_read_thinning;
    int rthgTh;

    double sim_coefficient;
    double sim_th;
    int max_core_errors;
    int max_pre_corrected_errors;
    int max_corrections_per_read;

} HMMParameters;

using namespace seqan;

char Complement(char c);
bool DiscardRead(const DnaString& read);
bool DiscardKmer(const char* read);

typedef seqan::FragmentStore<> THMMFragStore;

typedef long unsigned ulong;


#define MIN(X,Y) ((X) < (Y) ? (X) : (Y) )
#define MAX(X,Y) ((X) > (Y) ? (X) : (Y) )

typedef struct position_info {
    int core_emission_offset;
    int emission_length;
    int offset;
    
    position_info()
	: core_emission_offset(0), emission_length(0), offset(0) {
    }
    void Reset() {
	core_emission_offset = 0;
	offset = 0;
	emission_length = 0;
    }

} PositionInfo;

class ReadThread {
public:
    std::vector<int> reads;
    std::vector<int> estAlign;
    std::vector<int> multi;
    std::vector<DnaString> RComplement;
    void clear() {
	reads.clear();
	estAlign.clear();
	multi.clear();
	RComplement.clear();
    }

    void GetReadString(THMMFragStore& fragStore) {
	for (std::vector<int>::const_iterator it = reads.begin();
	     it != reads.end();
	     ++it) {
	    if (*it < 0) {
		RComplement.push_back(fragStore.readSeqStore[-*it-1]);
		reverseComplement(RComplement.back());
	    } else {
		RComplement.push_back(fragStore.readSeqStore[*it]);
	    }
	}
    }

    ~ReadThread() {
	
    }
};

typedef struct count_ {
  double val[4];
  double entropy;

  void Assign(double v[5], double entropy) {
    val[0] = v[0];
    val[1] = v[1];
    val[2] = v[2];
    val[3] = v[3];
    this->entropy = entropy;
  }
} CoreCount;

typedef struct cluster_ {
    DnaString core;
    ReadThread rthread;
    
    void clear() {
	rthread.clear();
	::clear(core);
    }

  // for printing out
  std::deque<CoreCount> core_prob;

} Cluster;


#endif
