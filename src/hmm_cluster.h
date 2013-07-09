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

#ifndef HMM_CLUSTER_H
#define HMM_CLUSTER_H


#include "common.h"
#include "read_finder.h"
#include "stats_keeper.h"

#include <iostream>
#include <ext/hash_map>
#include <limits>
#include <math.h>
#include <stdio.h>
#include <seqan/sequence.h>
#include <seqan/index.h>
#include <vector>

#include <gsl/gsl_rng.h>

using namespace seqan;

#define INF std::numeric_limits<double>::infinity()

typedef struct transmission {

    double Trans[N_HSTATES][N_HSTATES];
    double v[N_HSTATES][N_HSTATES];
    
    double scale, scale_norm;
    int k;
    
    HMMParameters* params;

    void reset() {
	for (int i = 0 ; i < N_HSTATES; ++i) {
	    for (int j = 0 ; j < N_HSTATES; ++j) {
		Trans[i][j] = params->Trans[i][j];
		v[i][j] = 0;
	    }
	}
	
	scale = pow(2, - params->alpha);
	scale_norm = 1 - scale;
	scale = scale / scale_norm;
	k = 0;
	
    }

    transmission(HMMParameters& params)
	: params(&params) {
	reset();
    }

    void addStats(int i, int j, double val) {
	v[i][j] += scale * val;
    }

    void update() {
	
	for (int i = 0 ; i < N_HSTATES; ++i) {
	    double s = 0;
	    for (int j = 0 ; j < N_HSTATES; ++j) {
		s += v[i][j];
	    }
	    for (int j = 0 ; j < N_HSTATES; ++j) {
		Trans[i][j] = log(v[i][j] / s);
	    }
	}

	++k;
	scale = pow(2 + k, - params->alpha);
	scale_norm *= (1-scale);
	scale = scale / scale_norm;

    }

    


} Transmission;


typedef struct Emission {

    /* DEBUG */
    double c[NDNA + 1];
    /* DEBUG */

    double v[NDNA];
    double s[NDNA];
    double scale, scale_norm;
    int count;
    int k;
    double diff;
    double entropy;
    Dna letter;


    static int m;
    static double alpha;

    Emission(double pseudo) {
	init(pseudo);
    }
    
    void init(double pseudo) {
	scale = pow(2, -alpha);
	scale_norm = 1 - scale;
	scale = scale / scale_norm;
	count = 0;
	k = 0;

	for (int i = 0; i < NDNA; ++i) {
	    v[i] = pseudo;
	    s[i] = -log((double)NDNA);
	    c[i] = 0;
	}

	c[NDNA] = 0;
	
	diff = 10000;
	entropy = 1000;
    }
    

    void update() {
	double t = 0;
	double m = -INF;

	for (int i = 0; i < NDNA; ++i) {
	    assert(v[i] >= 0);
	    t += v[i];
	}
	diff = 0;
	entropy = 0;
	for (int i = 0; i < NDNA; ++i) {
	    diff += fabs(v[i]/t - exp(s[i]));
	    s[i] = log(v[i] / t);
	    
	    entropy += -s[i] * v[i] / t;
	    if (s[i] > m) {
		m = s[i];
		letter = i;
	    }
	}
	
    }

    void addStats(int i, double val, double pseudo) {
	v[i] += val * scale;
	count++;
	if (count == m) {
	    for (int t = 0; t < NDNA; ++t) {
		v[t] += scale * (pseudo); //// ******
		assert(v[i] >= 0);
	    }

	    ++k;
	    scale = pow(2 + k, -alpha);
	    scale_norm *= (1-scale);
	    scale = scale / scale_norm;
	    count = 0;

	    update();
	}
    }

    void estimateFromAlignment(double pseudo) {
	entropy = 0;
	double m = -INF;
	for (int t = 0; t < NDNA; ++t) {
	    s[t] = log(c[t]+pseudo) - log(c[NDNA]+NDNA*pseudo);
	    entropy += -s[t] *
		(c[t] + pseudo) / (c[NDNA]+NDNA*pseudo);
	    if (s[t] > m) {
		m = s[t];
		letter = t;
	    }
	}
    }

    double getDiff() {
	double t = diff;
	diff = 0;
	return t;
    }

    void resetStats(double pseudo) {
	for (int i = 0; i < NDNA; ++i) {
	    v[i] = pseudo;
	    s[i] = -log((double)NDNA);
	    c[i] = 0;
	}
	c[NDNA] = 0;
    }
} Emission;

typedef struct _hstateval {
    double v[TOTAL_OFFSETS][N_HSTATES];
} Val;

typedef struct _traceelem {
    HMMHiddenStates v[TOTAL_OFFSETS][N_HSTATES];
} TraceElem;

typedef __gnu_cxx::hash_map<int, int> EPFLMap;

class HMMCluster {

private:

    THMMFragStore& fragStore;
    ReadFinder& finder;
    StatsKeeper& stats_keeper;
    
    HMMParameters *param;

    Transmission trans;

    // Emission buffer!!!
    const int length_emission_buffer;
    const int length_emission_left_buffer;
    const int length_emission_right_buffer;

    Emission* emission_buffer;
    Emission* emission_left_buffer;
    Emission* emission_right_buffer;

    int it_emission_buffer;
    int it_emission_left_buffer;
    int it_emission_right_buffer;

    int max_num_threads;

    double alpha;
    int  m;
    int max_read_length;
    size_t num_reads;
    
    PositionInfo main_pos, left_pos, right_pos;

    int g_sub_errors;
    int g_ins_errors;
    int g_del_errors;

    int g_stats_gathered;
    int g_stats_accepted;

    DnaString main_core;
    DnaString right_core;
    DnaString left_core;


    Val* alpha_array;
    Val* beta_array;
    TraceElem* trace_array;

    bool interesting_c;

    std::ostream* corrected_pos_os;

    int ToAllEmissionIndex(int i, int o);
    int ToEmissionBufferIndex(int i, int o,
			      int buffer_iter,
			      int buffer_length);
    double emit(Emission* buffer,
		int buffer_iter,
		int buffer_length,
		const DnaString& read,
		int i,
		int o);
    Dna getLetter(Emission* buffer,
		   int buffer_iter,
		   int buffer_length,
		   int i,
		   int o);

    int PaddingBases(int align, int length);
    int PaddingBasesLeftOnly(int align, int length);
    int PaddingBasesRightOnly(int align, int length);

    void updateEmit(Emission* buffer,
		    int buffer_iter,
		    int buffer_length,
		    const DnaString& read,
		    int i,
		    int o,
		    double v);

    void ForwardBackward(Emission* buffer,
			 int buffer_iter,
			 int buffer_length,
			 const DnaString& read,
			 int estAlignment);
    void ViterbiFillMatrix(Emission* buffer,
			   int buffer_iter,
			   int buffer_length,
			   const DnaString& read,
			   int estAlignment);
    bool Viterbi(DnaString& core,
		 PositionInfo& pos,
		 Emission* buffer,
		 int buffer_iter,
		 int buffer_length,
		 const DnaString& read,
		 int estAlignment,
		 DnaString& corrected,
		 std::vector<Error>& corrections,
		 std::ostream* os,
		 int& sub_errors,
		 int& ins_errors,
		 int& del_errors,
		 int read_id);

    bool BuildErrorProfilesHelper(DnaString& core,
				  PositionInfo& pos,
				  Emission* buffer,
				  int buffer_iter,
				  int buffer_length,
				  const DnaString& read,
				  int estAlignment,
				  int read_id,
				  EPFLMap& ep_map,
				  std::vector<std::vector<int> >& error_profiles,
				  std::ostream* os);

    void PrintEmissions(Emission* buffer,
			int buffer_iter,
			int buffer_length); 
    void PrintEmissionEntropy();
    double GetEntropyEmissions(Emission* buffer, int i);
    void GatherInitialReads(int cluster_id,
			    const DnaString& seed,
			    ReadThread& rthread,
			    Emission* buffer,
			    int buffer_iter,
			    int buffer_length);
    void GatherReads(int cluster_id,
		     const DnaString& core,
		     PositionInfo& pos,
		     int idx[], int idx_l,
		     ReadThread& rthread,
		     Emission* buffer,
		     int buffer_iter,
		     int buffer_length);

    void EMLearning(gsl_rng *& r,
		    DnaString& core,
		    PositionInfo& pos,
		    Emission* buffer,
		    int& buffer_iter,
		    int buffer_length,
		    ReadThread& rthread,
		    bool fleft,
		    bool fright);
    
    int ExtendLeft(DnaString& core,
		   Cluster& output,
		   PositionInfo& pos,
		   Emission* buffer,
		   int buffer_iter,
		   int buffer_length,
		   std::ostream* os);
    int ExtendRight(DnaString& core,
		    Cluster& output,
		    PositionInfo& pos,
		    Emission* buffer,
		    int buffer_iter,
		    int buffer_length,
		    std::ostream* os);
    int AssignReads(int cluster_id,
		    DnaString& core,
		    PositionInfo& pos,
		    Emission* buffer,
		    int buffer_iter,
		    int buffer_length,
		    ReadThread& rthread,
		    Cluster& output,
		    std::ostream* os);
    void ResetStatistics(Emission* buffer,
			 int buffer_iter,
			 int buffer_length);
    void ResetStatisticsLeftOnly();
    void ResetStatisticsRightOnly();
    void PopulateEmission(const ReadThread& rthread,
			  DnaString& core,
			  PositionInfo& pos,
			  Emission* buffer,
			  int& buffer_iter,
			  int buffer_length,
			  bool fleft,
			  bool fright,
			  std::ostream* os);
    void ModifyEmission(const ReadThread& rthread,
			DnaString& core,
			PositionInfo& pos,
			Emission* buffer,
			int& buffer_iter,
			int buffer_length,
			bool fleft,
			bool fright,
			std::ostream* os,
			int start = -1,
			int end = -1);
    double CheckDebugEntropyEmissions(Emission* buffer, int i);
    double ErrorProfilesSim(const std::vector<int>& x,
			    const std::vector<int>& y);
    int FilterReads(ReadThread& rthread,
		    DnaString& core,
		    PositionInfo& pos,
		    Emission* buffer,
		    int buffer_iter,
		    int buffer_length,
		    std::ostream* os);

public:
    HMMCluster(HMMParameters *param,
	       ReadFinder& finder,
	       StatsKeeper& stats_keeper,
	       int max_read_length,
	       double alpha, int m,
	       THMMFragStore& fragStore);
    ~HMMCluster();

    bool BuildCluster(const DnaString& baseread,
		      Cluster& output,
		      int readid,
		      const char* debugPath);

    void Reset();
    void PrintStats();
    void ReportPerf();
    void PrintReads(const ReadThread& rthread,
		    std::ostream* os,
		    int start = -1,
		    int end = -1) const;
    bool EstimateErrors(Emission* buffer,
			int buffer_iter,
			int buffer_length,
			const DnaString& read,
			int estAlignment);

    void setCorrectedPoOStream(std::ostream* os) {
      corrected_pos_os = os;
    }
};


#endif
