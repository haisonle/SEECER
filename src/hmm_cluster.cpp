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

#include "hmm_cluster.h"

#include <math.h>
#include <unistd.h>

#include <algorithm>
#include <fstream>
#include <iostream>
#include <sstream>

#include <gsl/gsl_randist.h>
#include <gsl/gsl_permutation.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_sort_vector.h>

using namespace seqan;
using namespace std;

#ifndef DEBUG
#define DEBUG 0
#endif

#ifndef _DEBUG_XX_
#define _DEBUG_XX_ 0
#endif

double sumexplog(double a, double b) {
    assert ( a <= 0 && b <= 0);

    // calculate log(e^a + e^b)
    if (std::isinf(-a))
	return b;
    if (std::isinf(-b))
	return a;

    if (b > a) {
	double t = a;
	a = b;
	b = t;
    }

    double t =  a + log( 1 + exp(b-a));
    assert (t <= 0);

    return t;
}

#define BUFFER_LENGTH 3
#define BUFFER_SEGMENT (max_read_length + INDELS)

int Emission::m;
double Emission::alpha;

HMMCluster::HMMCluster(HMMParameters *param,
		       ReadFinder& finder,
		       StatsKeeper& stats_keeper,
		       int max_read_length,
		       double alpha, int m,
		       THMMFragStore& fragStore)
    : fragStore(fragStore),
      finder(finder),
      stats_keeper(stats_keeper),
      param(param),
      //
      trans(*param),
      length_emission_buffer(4*BUFFER_LENGTH * BUFFER_SEGMENT),
      length_emission_left_buffer(BUFFER_LENGTH * BUFFER_SEGMENT),
      length_emission_right_buffer(BUFFER_LENGTH * BUFFER_SEGMENT),
      emission_buffer((Emission*) malloc(length_emission_buffer * sizeof(Emission))),
      emission_left_buffer(NULL),
      emission_right_buffer(NULL),
      it_emission_buffer(2*BUFFER_LENGTH * BUFFER_SEGMENT - 1),
      //
      alpha(alpha), m(m),
      max_read_length(max_read_length),
      num_reads(length(fragStore.readSeqStore)),
      g_sub_errors(0),
      g_ins_errors(0),
      g_del_errors(0),
      g_stats_gathered(0),
      g_stats_accepted(0),
      alpha_array((Val*) malloc((max_read_length+1) * sizeof(Val))),
      beta_array((Val*) malloc((max_read_length+1) * sizeof(Val))),
      trace_array((TraceElem*) malloc((max_read_length+1) * sizeof(TraceElem))),
      interesting_c(false) {


    for (int i = 0; i < TOTAL_OFFSETS; ++i) {
	for (int j = 0; j < N_HSTATES; ++j) {
	    alpha_array[0].v[i][j] = -INF;
	}
    }
    alpha_array[0].v[INDELS][STATEM] = 0;

    double t = 0;
    for (int i = 1; i <= INDELS; ++i) {
	t += trans.Trans[STATEM][STATED];
	alpha_array[0].v[INDELS-i][STATEM] = t;
    }
    t = 0;
    for (int i = 1; i <= INDELS; ++i) {
	t += trans.Trans[STATEM][STATEI];
	alpha_array[0].v[INDELS+i][STATEM] = t;
    }

}

HMMCluster::~HMMCluster() {

    free(alpha_array);
    free(beta_array);
    free(trace_array);

    free(emission_buffer);
}

void HMMCluster::Reset() {
    it_emission_buffer = 2*BUFFER_LENGTH * BUFFER_SEGMENT - 1;
    main_pos.Reset();
}

inline int HMMCluster::ToAllEmissionIndex(int i, int o) {
    int t = i - 1 + o - INDELS;
    assert (t >= 0);

    return t;
}

inline int HMMCluster::ToEmissionBufferIndex(int i, int o,
					     int buffer_iter,
					     int buffer_length) {
    int t = ToAllEmissionIndex(i,o) + buffer_iter;
    assert(t>=0);
    return t % buffer_length;
}

inline double HMMCluster::emit(Emission* buffer,
			       int buffer_iter,
			       int buffer_length,
			       const DnaString& read,
			       int i,
			       int o) {
    double t = buffer[ToEmissionBufferIndex(i, o, buffer_iter, buffer_length)]
	.s[(int)read[i-1]];
    assert(t < 0);
    return t;
}

inline Dna HMMCluster::getLetter(Emission* buffer,
				  int buffer_iter,
				  int buffer_length,
				  int i,
				  int o) {
    return buffer[ToEmissionBufferIndex(i,o, buffer_iter, buffer_length)]
	.letter;
}

inline void HMMCluster::updateEmit(Emission* buffer,
				   int buffer_iter,
				   int buffer_length,
				   const DnaString& read,
				   int i,
				   int o,
				   double v) {
    assert( v >= 0 );
    buffer[ToEmissionBufferIndex(i,o,buffer_iter, buffer_length)]
	.addStats(read[i - 1], v, param->pseudo);

}

int HMMCluster::PaddingBasesLeftOnly(int align, int length) {

    align += left_pos.offset;

    // prepend
    if (align < INDELS) {
	left_pos.offset += INDELS - align;
	left_pos.core_emission_offset += INDELS - align;
	left_pos.emission_length += INDELS - align;

	it_emission_left_buffer -= INDELS - align;
	it_emission_buffer -= INDELS - align;

	for (int i = 0; i < INDELS - align; ++i) {
	    emission_left_buffer[(i + it_emission_left_buffer
				  + length_emission_left_buffer)
			    % length_emission_left_buffer].init(param->pseudo);
	}

	it_emission_left_buffer =
	    (it_emission_left_buffer + length_emission_left_buffer)
	    % length_emission_left_buffer;

	align = INDELS;
    }

    // append
    int added = align + length  + INDELS - left_pos.emission_length;
    
    assert(added <= 0);
    
    return align;
}

int HMMCluster::PaddingBasesRightOnly(int align, int length) {


    align += right_pos.offset;

    // prepend
    assert(align > INDELS);

    // append
    int added = align + length  + INDELS - right_pos.emission_length;
    if (added > 0) {
      for (int i = 0; i < added; ++i) {
        emission_right_buffer[(right_pos.emission_length
			       + i + it_emission_right_buffer)
			      % length_emission_right_buffer]
	    .init(param->pseudo);
      }
      right_pos.emission_length += added;
    }
    
    return align;
}


int HMMCluster::PaddingBases(int align, int length) {

    align += main_pos.offset;

    // prepend
    if (align < INDELS) {
	main_pos.offset += INDELS - align;
	main_pos.core_emission_offset += INDELS - align;
	main_pos.emission_length += INDELS - align;
	it_emission_buffer -= INDELS - align;

	for (int i = 0; i < INDELS - align; ++i) {
	    emission_buffer[i + it_emission_buffer].init(param->pseudo);
	}

	assert(it_emission_buffer >= 0);

	align = INDELS;
    }


    // append
    int added = align + length  + INDELS - main_pos.emission_length;
    if (added > 0) {
	for (int i = 0; i < added; ++i) {
	    emission_buffer[main_pos.emission_length + i + it_emission_buffer]
		.init(param->pseudo);
	}
	main_pos.emission_length += added;
    }
    
    return align;
}

/*
  offset: -INDELS ... 0 ... INDELS
  offset = o - INDELS; o = 0 .. (2*INDELS + 1)

  i is the index of the read.
  j is the index of  the consensus sequence.

   j - i = o - INDELS = offset

=> j = i + o - INDELS

 */

void HMMCluster::ForwardBackward(Emission* buffer,
				 int buffer_iter,
				 int buffer_length,
				 const DnaString& read,
				 int estAlignment) {
    SEQAN_PROTIMESTART(fw);

#if DEBUG
    std::cerr << "FB: " << read << "(" << length(read)
	      << ") offset versus base: " << estAlignment << std::endl;
#endif

    int len = length(read);

    for (int i = 0; i < TOTAL_OFFSETS; ++i) {
	for (int j = 0; j < N_HSTATES; ++j) {
	    beta_array[len].v[i][j] =  0;
	}
    }

    for (int i = 1; i <= len; ++i) {

	for (int o = 0; o < TOTAL_OFFSETS; ++o) {
	    alpha_array[i].v[o][STATEM] =
		emit(buffer, buffer_iter, buffer_length,
		     read, i,  estAlignment + o) +
		sumexplog(sumexplog(trans.Trans[STATEM][STATEM]
				    + alpha_array[i-1].v[o][STATEM],
				    trans.Trans[STATED][STATEM]
				    + alpha_array[i-1].v[o][STATED]),
			  trans.Trans[STATEI][STATEM]
			  + alpha_array[i-1].v[o][STATEI]);
	    
	    if (o > 0) {
		alpha_array[i].v[o][STATED] =
		    sumexplog(sumexplog(trans.Trans[STATEM][STATED]
					+ alpha_array[i].v[o-1][STATEM],
					trans.Trans[STATED][STATED]
					+ alpha_array[i].v[o-1][STATED]),
			      trans.Trans[STATEI][STATED]
			      + alpha_array[i].v[o-1][STATEI]);
	    } else {
		alpha_array[i].v[o][STATED] = -INF;
	    }
	    
	    if (o < TOTAL_OFFSETS - 1) {
		alpha_array[i].v[o][STATEI] = // uniform emission probabilities
		    log(0.25) + 
		    sumexplog(sumexplog(trans.Trans[STATEM][STATEI]
					+ alpha_array[i-1].v[o+1][STATEM],
					trans.Trans[STATED][STATEI]
					+ alpha_array[i-1].v[o+1][STATED]),
			      trans.Trans[STATEI][STATEI]
			      + alpha_array[i-1].v[o+1][STATEI]);
	    } else {
		alpha_array[i].v[o][STATEI] = -INF;
	    }
	}

    }

    for (int i = len-1; i > 0; --i) {

	for (int o = TOTAL_OFFSETS - 1; o >= 0; --o) {
	    beta_array[i].v[o][STATEM] =
		emit(buffer, buffer_iter, buffer_length,
		     read, i + 1, estAlignment + o) +
		trans.Trans[STATEM][STATEM] + beta_array[i+1].v[o][STATEM];

	    if (o < TOTAL_OFFSETS-1)
		beta_array[i].v[o][STATEM] =
		    sumexplog(beta_array[i].v[o][STATEM],
			      trans.Trans[STATEM][STATED]
			      + beta_array[i].v[o+1][STATED]);
	    if (o > 0)
		beta_array[i].v[o][STATEM] =
		    sumexplog(beta_array[i].v[o][STATEM], log(0.25) +
			      trans.Trans[STATEM][STATEI] +
			      beta_array[i+1].v[o-1][STATEI]);

	    beta_array[i].v[o][STATED] =
		emit(buffer, buffer_iter, buffer_length, read,
		     i + 1, estAlignment + o) +
		trans.Trans[STATED][STATEM] + beta_array[i+1].v[o][STATEM];

	    if (o < TOTAL_OFFSETS-1)
		beta_array[i].v[o][STATED] =
		    sumexplog(beta_array[i].v[o][STATED],
			      trans.Trans[STATED][STATED]
			      + beta_array[i].v[o+1][STATED]);
	    
	    if (o > 0)
		beta_array[i].v[o][STATED] =
		    sumexplog(beta_array[i].v[o][STATED], log(0.25) +
			      trans.Trans[STATED][STATEI]
			      + beta_array[i+1].v[o-1][STATEI]);
	    
	    beta_array[i].v[o][STATEI] = 
		emit(buffer, buffer_iter, buffer_length,
		     read, i + 1, estAlignment + o) + 
		trans.Trans[STATEI][STATEM] + beta_array[i+1].v[o][STATEM];

	    if (o < TOTAL_OFFSETS-1)
		beta_array[i].v[o][STATEI] =
		    sumexplog(beta_array[i].v[o][STATEI],
			      trans.Trans[STATEI][STATED]
			      + beta_array[i].v[o+1][STATED]);
	    if (o > 0)
		beta_array[i].v[o][STATEI] =
		    sumexplog(beta_array[i].v[o][STATEI],
			      log(0.25) + trans.Trans[STATEI][STATEI]
			      + beta_array[i+1].v[o-1][STATEI]);
	    
	}
	
	
    }

#if DEBUG
    
    // Print debug
    for (int i = 0; i <= len; ++i) {
	std::cerr << i << ":\t";
	for (int o = 0 ; o< TOTAL_OFFSETS; ++o) {
	    std::cerr << "( " << exp(alpha_array[i].v[o][STATEM])
		      << " " << exp(alpha_array[i].v[o][STATEI])
		      << " " << exp(alpha_array[i].v[o][STATED]) << ")\t";
	}
	std::cerr << std::endl;
    }
    std::cerr << "=========================" << std::endl;

    for (int i = 0; i <= len; ++i) {
	std::cerr << i << ":\t";
	for (int o = 0 ; o < TOTAL_OFFSETS; ++o) {
	    std::cerr << "( " << exp(beta_array[i].v[o][STATEM])
		      << " " << exp(beta_array[i].v[o][STATEI])
		      << " " << exp(beta_array[i].v[o][STATED])
		      << ")\t";
	}
	std::cerr << std::endl;
    }

#endif

    double totalp;
    double ttotalp;
    
    for (int a = 1; a <= len; ++a) {
#if DEBUG
	std::cerr << a << ":\t";
#endif
	ttotalp = -INF;
	for (int i = 0; i < TOTAL_OFFSETS; ++i) {
	    totalp = -INF;
	    for (int j = 0; j < N_HSTATES; ++j) {
		if (j != STATED)
		    totalp = sumexplog(totalp,
				       alpha_array[a].v[i][j]
				       + beta_array[a].v[i][j]);
	    }
#if DEBUG
	    std::cerr << totalp << "\t";
#endif
	    ttotalp = sumexplog(totalp, ttotalp);
	}
#if DEBUG
	std::cerr << " = " << ttotalp << " ~ " << exp(ttotalp) << std::endl;
#endif
    }

    // Updating parameters
    for (int i = 1; i <= len; ++i) {

	for (int o = 0; o < TOTAL_OFFSETS; ++o) {

	    double v = exp(alpha_array[i].v[o][STATEM]
			   + beta_array[i].v[o][STATEM] - ttotalp);
#if DEBUG
	    if ( v > 1 ) {
		std::cerr << alpha_array[i].v[o][STATEM]
		    + beta_array[i].v[o][STATEM]
			  << " " << ttotalp << std::endl;
	    }
#endif
	    updateEmit(buffer, buffer_iter, buffer_length, 
		       read, i, estAlignment + o, v);
	}

    }

#if DEBUG
    std::cerr << "Forward takes " << SEQAN_PROTIMEDIFF(fw) << std::endl;
#endif

#if 0 // updating the transmissions
    for (int i = 1; i < len; ++i) {
	for (int o = 0; o < TOTAL_OFFSETS; ++o) {
	    
	    double e = emit(buffer, buffer_iter, buffer_length,
			    read, i,  estAlignment + o);

	    trans.addStats(STATEM, STATEM,
			   exp(e + trans.Trans[STATEM][STATEM] 
			       + alpha_array[i-1].v[o][STATEM]
			       + beta_array[i].v[o][STATEM] - ttotalp));
	    trans.addStats(STATED, STATEM,
			   exp(e + trans.Trans[STATED][STATEM] 
			       + alpha_array[i-1].v[o][STATED]
			       + beta_array[i].v[o][STATEM] - ttotalp));
	    trans.addStats(STATEI, STATEM,
			   exp(e + trans.Trans[STATEI][STATEM] 
			       + alpha_array[i-1].v[o][STATEI]
			       + beta_array[i].v[o][STATEM] - ttotalp));
	    
	    if (o > 0) {
		trans.addStats(STATEM, STATED,
			       exp(trans.Trans[STATEM][STATED]
				   + alpha_array[i].v[o-1][STATEM]
				   + beta_array[i].v[o][STATED] - ttotalp));
		trans.addStats(STATED, STATED,
			       exp(trans.Trans[STATED][STATED]
				   + alpha_array[i].v[o-1][STATED]
				   + beta_array[i].v[o][STATED] - ttotalp));
		trans.addStats(STATEI, STATED,
			       exp(trans.Trans[STATEI][STATED]
				   + alpha_array[i].v[o-1][STATEI]
				   + beta_array[i].v[o][STATED] - ttotalp));
	    }
	    
	    if (o < TOTAL_OFFSETS - 1) {
		trans.addStats(STATEM, STATEI,
			       0.25 * exp(trans.Trans[STATEM][STATEI]
					  + alpha_array[i-1].v[o+1][STATEM]
					  + beta_array[i].v[o][STATEI]
					  - ttotalp));
		trans.addStats(STATED, STATEI,
			       0.25 * exp(trans.Trans[STATED][STATEI]
					  + alpha_array[i-1].v[o+1][STATED]
					  + beta_array[i].v[o][STATEI]
					  - ttotalp));
		trans.addStats(STATEI, STATEI,
			       0.25 * exp(trans.Trans[STATEI][STATEI]
					  + alpha_array[i-1].v[o+1][STATEI]
					  + beta_array[i].v[o][STATEI]
					  - ttotalp));
	    }
	}
    }
    
#endif
    
}

double max_helper(double matchstate,
		  double deletionstate,
		  double insertionstate,
		  HMMHiddenStates* trace) {
    if (matchstate > deletionstate) {
	if (matchstate > insertionstate) {
	    *trace = STATEM;
	    return matchstate;
	} else {
	    *trace = STATEI;
	    return insertionstate;
	}
    } else {
	if (deletionstate > insertionstate) {
	    *trace = STATED;
	    return deletionstate;
	} else  {
	    *trace = STATEI;
	    return insertionstate;
	}
	
    }
}

#define VITERBI_DEBUG 1

void HMMCluster::ViterbiFillMatrix(Emission* buffer,
				   int buffer_iter,
				   int buffer_length,
				   const DnaString& read,
				   int estAlignment) {
    
    int len = length(read);
    
    for (int i = 1; i <= len; ++i) {
	
	for (int o = 0; o < TOTAL_OFFSETS; ++o) {
	    alpha_array[i].v[o][STATEM] =
		emit(buffer, buffer_iter, buffer_length, read,
		     i,  estAlignment + o) +
		max_helper(trans.Trans[STATEM][STATEM]
			   + alpha_array[i-1].v[o][STATEM],
			   trans.Trans[STATED][STATEM]
			   + alpha_array[i-1].v[o][STATED],
			   trans.Trans[STATEI][STATEM]
			   + alpha_array[i-1].v[o][STATEI],
			   &trace_array[i].v[o][STATEM]);

	    if (o > 0) {
		alpha_array[i].v[o][STATED] =
		    max_helper(trans.Trans[STATEM][STATED]
			       + alpha_array[i].v[o-1][STATEM],
			       trans.Trans[STATED][STATED]
			       + alpha_array[i].v[o-1][STATED],
			       trans.Trans[STATEI][STATED]
			       + alpha_array[i].v[o-1][STATEI],
			       &trace_array[i].v[o][STATED]);
	    } else {
		alpha_array[i].v[o][STATED] = -INF;
	    }
	    

	    if (o < TOTAL_OFFSETS - 1) {
		alpha_array[i].v[o][STATEI] =
		    log(0.25) +
		    max_helper(trans.Trans[STATEM][STATEI]
			       + alpha_array[i-1].v[o+1][STATEM],
			       trans.Trans[STATED][STATEI]
			       + alpha_array[i-1].v[o+1][STATED],
			       trans.Trans[STATEI][STATEI]
			       + alpha_array[i-1].v[o+1][STATEI],
			       &trace_array[i].v[o][STATEI]);
	    } else {
		alpha_array[i].v[o][STATEI] = -INF;
	    }
	    
	}
    }
    
#if DEBUG
    // Print debug
    for (int i = 0; i <= len; ++i) {
	std::cerr << i << ":\t";
	for (int o = 0 ; o< TOTAL_OFFSETS; ++o) {
	    std::cerr << "( " << exp(alpha_array[i].v[o][STATEM])
		      << " " << exp(alpha_array[i].v[o][STATEI])
		      << " " << exp(alpha_array[i].v[o][STATED]) << ")\t";
	}
	std::cerr << std::endl;
    }
    std::cerr << "=========================" << std::endl;
#endif

}


bool HMMCluster::Viterbi(DnaString& core,
			 PositionInfo& pos,
			 Emission* buffer,
			 int buffer_iter,
			 int buffer_length,
			 const DnaString& read,
			 int estAlignment,
			 DnaString& corrected,
			 vector<Error>& corrections,
			 ostream* os,
			 int& sub_errors,
			 int& ins_errors,
			 int& del_errors,
			 int read_id) {

#if DEBUG
    std::cerr << "Viterbi: " << read
	      << "(" << length(read)
	      << ") offset versus base: "
	      << estAlignment << std::endl;
#endif

    int len = length(read);

    ViterbiFillMatrix(buffer, buffer_iter,
		      buffer_length, read, estAlignment);
    
    
    double max_v = -INF; int b_state = -1; int b_o = -1;
    for (int o = 0 ; o < TOTAL_OFFSETS; ++o) {
	for (int s = 0; s < N_HSTATES; ++s) {
	    if (alpha_array[len].v[o][s] > max_v) {
		    b_state = s;
		    b_o = o;
		    max_v = alpha_array[len].v[o][s];
		}
	}
    }

    // Tracing back and printing
    stringstream ss(stringstream::in | stringstream::out);
    vector<int> position;


    // boundary reads should be remained for next rounds
    bool fitted = ToAllEmissionIndex(len, estAlignment + b_o)
	< pos.core_emission_offset + static_cast<int>(length(core));

    if (!fitted) {
	return false;
    }

    int i = len;


    sub_errors = 0;
    ins_errors = 0;
    del_errors = 0;

#if CORRECTIONS_SAVING
    clear(corrections);
    Error e;
#endif

    while (i > 0) {

	assert(b_o >= 0 && b_o <= TOTAL_OFFSETS);

	HMMHiddenStates nb_state = trace_array[i].v[b_o][b_state];

#if VITERBI_DEBUG
	string st;
	if (os) {
	    int ei = i - 1 + estAlignment + b_o - INDELS; 
	    stringstream sss(stringstream::out);
	    sss << ei;
	    st = sss.str();
	    reverse(st.begin(), st.end());
	}
#endif

	int idx = ToAllEmissionIndex(i, estAlignment + b_o);

	// boundary reads should be remained for next rounds
	if (idx < pos.core_emission_offset) {
	    return false;
	}

	if ( b_state == STATEM) {
	    assert(emit(buffer, buffer_iter, buffer_length, 
			read, i,  estAlignment + b_o)
		   + trans.Trans[nb_state][STATEM]
		   +  alpha_array[i-1].v[b_o][nb_state]  -
		   alpha_array[i].v[b_o][b_state] < 1e-10);

	    if (read[i - 1] != core[idx - pos.core_emission_offset]
		// only correct if the margin of correction is at least 0.1
		// (allows SNP)
		&& exp(buffer[(idx + buffer_iter)
			    % buffer_length].s[static_cast<int>(read[i - 1])])
		<= exp(buffer[(idx + buffer_iter) % buffer_length]
		       .s[static_cast<int>(core[idx - pos.core_emission_offset])])
		- param->emit_delta
		) {

#if CORRECTIONS_SAVING
		// substitution
		e.pos = i - 1;
		e.a = read[i - 1];
		e.b = core[idx - pos.core_emission_offset];

		corrections.push_back(e);
#endif

		sub_errors++;
#if VITERBI_DEBUG
		ss << "*";
#endif
		append(corrected, core[idx - pos.core_emission_offset]);
	    } else {
#if VITERBI_DEBUG
		//ss << read[i-1] << st << ")"<<b_o<<"("<< "\t";
		ss << core[idx - pos.core_emission_offset];//read[i-1];
#endif
		append(corrected, read[i - 1]);
	    }

	    --i;
	} else if (b_state == STATED) {
	    assert(trans.Trans[nb_state][STATED]
		   +  alpha_array[i].v[b_o-1][nb_state] - 
		   alpha_array[i].v[b_o][b_state] < 1e-10);
#if VITERBI_DEBUG
	    if (os) {
		//ss << "*" << st << ")"<<b_o<<"("<< "\t";
		ss << "-";
	    }
#endif
	    append(corrected, core[idx - pos.core_emission_offset]);
	
#if CORRECTIONS_SAVING    
	    e.pos = i;
	    e.a = '-';
	    e.b = core[idx - pos.core_emission_offset];

	    corrections.push_back(e);
#endif

	    del_errors++;

	    --b_o;

	} else {
	    assert(trans.Trans[nb_state][STATEI] + log(0.25) + 
		   +  alpha_array[i-1].v[b_o+1][nb_state] - 
		   alpha_array[i].v[b_o][b_state] < 1e-10);
#if VITERBI_DEBUG
	    if (os) {
		//ss << read[i-1];// << ")"<<b_o<<"("<< st;
	    }
#endif

#if CORRECTIONS_SAVING
	    e.pos = i - 1;
	    e.a = '+';
	    e.b = core[idx - pos.core_emission_offset];

	    corrections.push_back(e);
#endif

	    ins_errors++;

	    ++b_o;
	    --i;
	}
	b_state = nb_state;
    }

    max_v /= len;

#if VITERBI_DEBUG
    if (os) {
	if (estAlignment + b_o > INDELS + 1) {
	    int limit = ToEmissionBufferIndex(0, estAlignment + b_o,
					      buffer_iter, buffer_length);
	    if (limit < buffer_iter) {
		limit += buffer_length;
	    }
	    for (int o = buffer_iter; o <= limit; ++o) {
		ss << " ";//"\t";
	    }
	}
	string s = ss.str();
	reverse(s.begin(), s.end());
	
	*os << s << " " << read_id << " " <<  max_v;
	if (fitted && max_v > param->cluster_llh_th) {
	    *os << "* " << sub_errors + ins_errors + del_errors;
	}
	*os << std::endl;
    }
    
#endif

    //    std::cerr << max_v << std::endl;


    DnaStringReverse b(corrected);
    corrected = b;

    /*
    if (sub_errors + ins_errors + del_errors) {
	std::cerr << read << "\n" << corrected << "\n*****\n";// << corrections.size() << "\n";
    }
    */

#if _DEBUG_XX_
    
    /*
    if (fitted) {
      std::cerr << "#\t" << max_v << "\t"
		<< sub_errors + ins_errors + del_errors << std::endl;
    }
    */

    if (os) {
      *os << max_v * len << "\t" << max_v << "\t"
	  << sub_errors + ins_errors + del_errors << std::endl;
    }
#endif

    return fitted && max_v > param->cluster_llh_th;
}

bool HMMCluster::EstimateErrors(Emission* buffer,
				int buffer_iter,
				int buffer_length,
				const DnaString& read,
				int estAlignment) {
    int len = length(read);

    /*
    std::cerr << "Estimate errors " << read
	      << " " << estAlignment << std::endl;
    */

    ViterbiFillMatrix(buffer, buffer_iter,
		      buffer_length, read, estAlignment);
    
    
    double max_v = -INF; int b_state = -1; int b_o = -1;
    for (int o = 0 ; o < TOTAL_OFFSETS; ++o) {
	for (int s = 0; s < N_HSTATES; ++s) {
	    if (alpha_array[len].v[o][s] > max_v) {
		    b_state = s;
		    b_o = o;
		    max_v = alpha_array[len].v[o][s];
		}
	}
    }



    int i = len;
    int errors = 0;

    while (i > 0) {

	assert(b_o >= 0 && b_o <= TOTAL_OFFSETS);

	HMMHiddenStates nb_state = trace_array[i].v[b_o][b_state];

	if ( b_state == STATEM) {
	    assert(emit(buffer, buffer_iter, buffer_length, 
			read, i,  estAlignment + b_o)
		   + trans.Trans[nb_state][STATEM]
		   +  alpha_array[i-1].v[b_o][nb_state]  -
		   alpha_array[i].v[b_o][b_state] < 1e-10);

	    if (read[i - 1] != getLetter(buffer, buffer_iter, buffer_length, 
					 i,  estAlignment + b_o)) {
		errors++;
	    }

	    --i;
	} else if (b_state == STATED) {
	    assert(trans.Trans[nb_state][STATED]
		   +  alpha_array[i].v[b_o-1][nb_state] - 
		   alpha_array[i].v[b_o][b_state] < 1e-10);
	    errors++;
	    --b_o;

	} else {
	    assert(trans.Trans[nb_state][STATEI] + log(0.25) + 
		   +  alpha_array[i-1].v[b_o+1][nb_state] - 
		   alpha_array[i].v[b_o][b_state] < 1e-10);

	    errors++;
	    ++b_o;
	    --i;
	}
	b_state = nb_state;
    }
    
    return errors < param->max_core_errors;
}

// Return false if the read should be removed!
bool HMMCluster::BuildErrorProfilesHelper(DnaString& core,
					  PositionInfo& pos,
					  Emission* buffer,
					  int buffer_iter,
					  int buffer_length,
					  const DnaString& read,
					  int estAlignment,
					  int read_id,
					  EPFLMap& ep_map,
					  vector<vector<int> >& error_profiles,
					  ostream* os) {
    bool ret = true;
    int len = length(read);

    ViterbiFillMatrix(buffer, buffer_iter,
		      buffer_length, read, estAlignment);
    
    
    double max_v = -INF; int b_state = -1; int b_o = -1;
    for (int o = 0 ; o < TOTAL_OFFSETS; ++o) {
	for (int s = 0; s < N_HSTATES; ++s) {
	    if (alpha_array[len].v[o][s] > max_v) {
		    b_state = s;
		    b_o = o;
		    max_v = alpha_array[len].v[o][s];
		}
	}
    }

#if _DEBUG_XX_
    if (os) {
	*os << read_id << ": ";
    }
#endif
    int i = len;

    int core_errors = 0;
    int total_errors = 0;

    while (i > 0) {


	bool error = false;
	assert(b_o >= 0 && b_o <= TOTAL_OFFSETS);

	HMMHiddenStates nb_state = trace_array[i].v[b_o][b_state];
	int idx = ToAllEmissionIndex(i, estAlignment + b_o);


	if ( b_state == STATEM) {
	    assert(emit(buffer, buffer_iter, buffer_length, 
			read, i,  estAlignment + b_o)
		   + trans.Trans[nb_state][STATEM]
		   +  alpha_array[i-1].v[b_o][nb_state]  -
		   alpha_array[i].v[b_o][b_state] < 1e-10);

	    if (read[i - 1] != getLetter(buffer, buffer_iter, buffer_length, 
					 i,  estAlignment + b_o)) {
		error = true;
	    }

	    --i;
	} else if (b_state == STATED) {
	    assert(trans.Trans[nb_state][STATED]
		   +  alpha_array[i].v[b_o-1][nb_state] - 
		   alpha_array[i].v[b_o][b_state] < 1e-10);

	    error = true;

	    --b_o;

	} else {
	    assert(trans.Trans[nb_state][STATEI] + log(0.25) + 
		   +  alpha_array[i-1].v[b_o+1][nb_state] - 
		   alpha_array[i].v[b_o][b_state] < 1e-10);

	    error = true;
	    ++b_o;
	    --i;
	}
	b_state = nb_state;

	// dont care about extra possible indels
	if (error
	    && idx >= INDELS && idx < pos.emission_length - INDELS) {

	    idx -= pos.core_emission_offset;
	    if (idx < static_cast<int>(length(core)) && idx >= 0) {
		core_errors++;
		if (core_errors > param->max_core_errors) {
		    ret = false;
		    goto quit;
		}
	    } 
	    
	    {
		
		int error_profiles_idx;
		if (ep_map.find(idx) == ep_map.end()) {
		    ep_map.insert(pair<int, int>(idx, ep_map.size()));
		    error_profiles_idx = error_profiles.size();
		    error_profiles.resize(error_profiles.size() + 1);
		} else {
		    error_profiles_idx = ep_map[idx];
		}

		if (error_profiles[error_profiles_idx].empty() ||
		    // Because of possible indels
		    error_profiles[error_profiles_idx].back() != read_id) {
#if _DEBUG_XX_
		    if (os) {
			*os << "[" << idx << " " << error_profiles_idx << "] ";
		    }
#endif
		    error_profiles[error_profiles_idx].push_back(read_id);
		}
	    }

	    ++total_errors;

	}
    }

 quit:
#if _DEBUG_XX_
    if (os) {
	*os << std::endl;
    }
#endif
    return ret && total_errors < param->max_pre_corrected_errors;
}

double HMMCluster::ErrorProfilesSim(const vector<int>& x,
				    const vector<int>& y) {
    double sim = 0;
    vector<int>::const_iterator itx = x.begin();
    vector<int>::const_iterator ity = y.begin();
    while (itx != x.end()
	   && ity != y.end()) {
	if (*itx == *ity) {
	    ++sim;
	    ++itx;++ity;
	} else if (*itx < *ity) {
	    ++itx;
	} else {
	    ++ity;
	}
    }
    
    return sim;
}

int HMMCluster::FilterReads(ReadThread& rthread,
			    DnaString& core,
			    PositionInfo& pos,
			    Emission* buffer,
			    int buffer_iter,
			    int buffer_length,
			    ostream* os) {
    
    int count = 0;
    bool* mask = static_cast<bool*>(calloc(rthread.reads.size(), sizeof(bool)));
    
    EPFLMap ep_map;
    vector<vector<int> > error_profiles;
    // Building Profiles
    
    for (size_t i = 0; i < rthread.reads.size(); ++i) {
	DnaString& read = rthread.RComplement[i];
	
	if (!BuildErrorProfilesHelper(core, pos,
				      buffer, buffer_iter, buffer_length,
				      read, rthread.estAlign[i] + pos.offset, i,
				      ep_map,
				      error_profiles,
				      os)) {
#if _DEBUG_XX_
	    if (os) {
		*os << "Read " << i
		    << " has so many errors!" << std::endl;
	    }
#endif
	    count++;
	    mask[i] = true;
	}
    }
    
    SEQAN_PROTIMESTART(FAA);
	    
    // If we have enough reads to do TODO
    if (rthread.reads.size() - count > 2 && error_profiles.size() > 1) {
	
       // Remove error profiles with less than 2 reads
      int i = 0;
      int j = error_profiles.size() - 1;

      while (j > i) {
	  while (i < j && error_profiles[i].size() > 3) {
	      ++i;
	  }
	  while (i < j && error_profiles[j].size() <= 3) {
	      --j;
	  }
	  
	  if (i < j) {
	      error_profiles[i].swap(error_profiles[j]);
	      ++i; --j;
	  }
      }
      
      if (i == j && error_profiles[i].size() > 3) {
	++i;
      }

    int size = i;
    assert(size < pos.emission_length);
    error_profiles.resize(size);
	// Calculate the similarity
#if DEBUG
	std::cerr << "Errorneous bases: " << size << " of "
		  << pos.emission_length - length(core) << " new bases"
		  << std::endl;
#endif
	if (size > 1) {
	    gsl_matrix *A = gsl_matrix_alloc(size, size);
	    gsl_matrix *V = gsl_matrix_alloc(size, size);
	    gsl_vector *E = gsl_vector_alloc(size);
	    double *D = static_cast<double*>(calloc(size, sizeof(double)));
	    
	    // constructing the normalized Laplacian
	    vector<vector<int> >::const_iterator it1 = error_profiles.begin();
	    for (int i = 0; i < size - 1; ++i) {
		vector<vector<int> >::const_iterator it2 = it1 + 1;
		
		for (int j = i + 1; j < size; ++j) {
		    double v = ErrorProfilesSim(*it1, *it2);
		    v = 1 / (1 + exp( - (v-3)));
		    D[i]+=v;
		    D[j]+=v;
		    gsl_matrix_set(A, i, j, -v);
		    //gsl_matrix_set(A, j, i, -v);
		    ++it2;
		}
		++it1;
	    }
#if _DEBUG_XX_
	    if (os) {
		*os << "------------------------" << std::endl;
		for (int i = 0; i < size; ++i) {
		    for (int j = 0; j < size; ++j) {
			double v;
			if (i < j) {
			    v = -gsl_matrix_get(A, i, j);
			} else if (i > j) {
			    v = -gsl_matrix_get(A, j, i);
			} else {
			    v = 1;
			}
			*os << v << " ";
		    }
		    *os << std::endl;
		}
	    }
#endif	    
	    for (int i = 0; i < size; ++i) {
		gsl_matrix_set(A, i, i, D[i] / ( D[i] + 1));
	    }
	    
	    for (int i = 0; i < size; ++i) {
		D[i] = sqrt(D[i] + 1);
	    }
	    for (int i = 0; i < size - 1; ++i) {
		for (int j = i + 1; j < size; ++j) {
		    double v = gsl_matrix_get(A, i,j);
		    v /= D[i] * D[j];
		    gsl_matrix_set(A, i, j, v);
		    gsl_matrix_set(A, j, i, v);
		}
	    }

	    /*	    
	    if (os) {
		*os << "------------------------" << std::endl;
		for (int i = 0; i < size; ++i) {
		    for (int j = 0; j < size; ++j) {
			double v = gsl_matrix_get(A, i, j);
			*os << v << " ";
		    }
		    *os << std::endl;
		}
	    }
	    */

	    // Solving the eigenvalue problem
	    gsl_eigen_symmv_workspace * ws = gsl_eigen_symmv_alloc(size);
	    gsl_eigen_symmv(A, E, V, ws);
	    	    
	    /*
	    if (os) {
		*os << "------------------------" << std::endl;
		for (int i = 0; i < size; ++i) {
		    for (int j = 0; j < size; ++j) {
			*os << gsl_matrix_get(V, i, j ) << " ";
		    }
		    *os << std::endl;
		}
		*os << "------------------------" << std::endl;
		for (int i = 0; i < size; ++i) {
		    *os << gsl_vector_get(E, i) << " ";
		}
		*os << std::endl;
	    }
	    */

	    // sorting the eigens
	    gsl_permutation * perm = gsl_permutation_alloc(size);
	    if (gsl_sort_vector_index(perm, E) != 0) {
		abort();
	    }
	    // finding the number of clusters
	    double ce = gsl_vector_get(E, gsl_permutation_get(perm, 0));
	    double max = -1;
	    int nclusters = size;
	    for (int i = 1; i < size; ++i) {
		double v = gsl_vector_get(E, gsl_permutation_get(perm, i));
		if (v - ce > max) {
		    max = v - ce;
		    nclusters = i;
		}
		ce = v;
	    }

	    /*
	    if (nclusters == 1) {
		nclusters = 2;
	    }
	    */

#if _DEBUG_XX_
	    if (os) {
		*os << " N Clusters " << nclusters << std::endl;
	    }
#endif

	    // if (nclusters > 0)
            // if (nclusters > 1 && nclusters < size)
	    {
		// Constructing X^t
		gsl_matrix *X = gsl_matrix_alloc(nclusters, size);
		for (int i = 0; i < nclusters; ++i) {
		    gsl_vector_const_view col =
			gsl_matrix_const_column(V, gsl_permutation_get(perm, i));
		    gsl_matrix_set_row(X, i, &col.vector);
		}
#if _DEBUG_XX_
		if (os) {
		    *os << "------------------------" << std::endl;
		    for (int i = 0; i < nclusters; ++i) {
			for (int j = 0; j < size; ++j) {
			*os << gsl_matrix_get(X, i, j ) << " ";
			}
			*os << std::endl;
		    }
		}
#endif
		// QR decomposition
		int signum;
		gsl_vector *QRws = gsl_vector_alloc(size);
		// X = [R11 R12]
		gsl_vector *tau = gsl_vector_alloc(nclusters);
		if (gsl_linalg_QRPT_decomp(X, tau, perm, &signum, QRws) != 0) {
		    abort();
		}
		// Setting the lower part of X to 0 - clearing out Q
		for (int j = 0; j < nclusters; ++j) {
		    for (int i = j+1; i < nclusters; ++i) {
			gsl_matrix_set(X, i, j, 0 );
		    }
		}
		/*
		if (os) {
		    *os << "------------------------" << std::endl;
		    for (int i = 0; i < nclusters; ++i) {
			for (int j = 0; j < size; ++j) {
			    *os << gsl_matrix_get(X, i, j ) << " ";
			}
			*os << std::endl;
		    }
		    *os << "------------------------" << std::endl;
		    for (int j = 0; j < size; ++j) {
			*os << gsl_permutation_get(perm, j ) << " ";
		    }
		    *os << std::endl;
		}
		*/

		// Using LU decomposition with U = R11
		gsl_matrix_const_view R1 =
		    gsl_matrix_const_submatrix(X, 0, 0, nclusters, nclusters);
		gsl_permutation * permI = gsl_permutation_calloc(nclusters);
		vector<int>* sets = new vector<int>[nclusters];
		for (int i = 0; i < nclusters; i++) {
		    int c = gsl_permutation_get(perm, i);
		    sets[i].push_back(c);
		}
		for (int i = nclusters; i < size; i++) {
		    int c = gsl_permutation_get(perm, i);
		    gsl_vector_const_view col =
			gsl_matrix_const_column(X, i);
		    gsl_linalg_LU_solve(&R1.matrix, permI, &col.vector, tau);
		    // finding the maximum in absolute values
		    double m = abs(gsl_vector_get(tau, 0));
		    int mi = 0;
		    for (int j = 1; j < nclusters; j++) {
			double v = abs(gsl_vector_get(tau, j));
			if (v  > m) {
			    mi = j;
			    m = v;
			}
		    }
		    if (m > 0.5) {
			sets[mi].push_back(c);
		    }
		}
		
		int* errors_count =
		    static_cast<int*>(calloc(rthread.reads.size(), sizeof(int)));
		
		for (int i = 0; i < nclusters; ++i) {
#if _DEBUG_XX_
		    if (os) {
			*os << i << ": ";
			for (vector<int>::const_iterator it = sets[i].begin();
			     it != sets[i].end();
			     ++it) {
			    *os << *it << " ";
			}
			*os << std::endl;
		    }
#endif
		    if (sets[i].size() >= 2) {
			for (vector<int>::const_iterator it = sets[i].begin();
			     it != sets[i].end();
			     ++it) {
			    for (vector<int>::const_iterator itx =
				     error_profiles[*it].begin();
				 itx != error_profiles[*it].end(); ++itx) {
				errors_count[*itx]++;
			    }
			}
		    }
		    int th = MIN(3, 0.5 * sets[i].size());
		    for (int i = 0; i < static_cast<int>(rthread.reads.size());
			 ++i) {
			if (errors_count[i] > th) {
			    if (!mask[i]) {
				/*
				if (os) {
				    *os << "Remove "
				        << rthread.reads[i]
					<<std::endl;
				}
				*/
				mask[i] = true;
				++count;

			    }
			}
			errors_count[i] = 0;
		    }
		}
		free(errors_count);
		
		delete[] sets;
		gsl_vector_free(QRws);
		gsl_vector_free(tau);
		gsl_matrix_free(X);

		gsl_permutation_free(permI);
	    }

	    gsl_permutation_free(perm);
	    gsl_eigen_symmv_free(ws);
	    gsl_matrix_free(A);
	    gsl_matrix_free(V);
	    gsl_vector_free(E);
	    free(D);
	}

	if (os) {
	    *os << "###\t" << size << "\t" << rthread.reads.size()
		<< "\t" << count << "\t" << SEQAN_PROTIMEDIFF(FAA) << std::endl;
	}
	

    }
    
#if DEBUG
    std::cerr << "Filter reads take " << SEQAN_PROTIMEDIFF(FAA) <<
	std::endl;
#endif
    
#if _DEBUG_XX_
    if (os) {
	*os << "Remove " << count << " reads "
	    << rthread.reads.size() - count << " reads left" << std::endl;
	
	for (int j = 0; j < pos.core_emission_offset; ++j) {
	    int idx = (j + buffer_iter) % buffer_length;
	    *os << Dna(buffer[idx].letter);
	}
	*os << "|||||";
	for (int j = pos.core_emission_offset + length(core);
	     j < pos.emission_length;
	     ++j) {
	    int idx = (j + buffer_iter) % buffer_length;
	    *os << Dna(buffer[idx].letter);
	}
	*os << std::endl;
	
	*os << "BEFORE THINNING" << std::endl;
    }
#endif

    PrintReads(rthread, os);

    // bring all kept reads to the front, filtered reads
    // to the back of the vectors
    if (count > 0) {
	
	int i = 0;
	int j = rthread.reads.size() - 1;
	while (j > i) {
	    while (j > i && !mask[i]) {
		i++;
	    }
	    while (j > i && mask[j]) {
		j--;
	    }
	    if (j > i) {
		// swapping;
		int t;
		t = rthread.reads[j];
		rthread.reads[j] = rthread.reads[i];
		rthread.reads[i] = t;
		t = rthread.estAlign[j];
		rthread.estAlign[j] = rthread.estAlign[i];
		rthread.estAlign[i] = t;
		t = rthread.multi[j];
		rthread.multi[j] = rthread.multi[i];
		rthread.multi[i] = t;
		swap(rthread.RComplement[i], rthread.RComplement[j]);
		i++;j--;
	    }
	}
	
	if (i == j && !mask[i]) {
	    ++i;
	}

	assert(count + i == static_cast<int>(rthread.reads.size()));

	// reporting this read as failure
	if (param->restrict_failures) {
	  for (int j = rthread.reads.size() - 1; j >= ((int) rthread.reads.size()) - count; --j)
	    stats_keeper.ReportFailure(rthread.reads[j]);
	}		


#if _DEBUG_XX_
	if (os) {
	    *os << "THINNING" << std::endl;
	    PrintReads(rthread, os, 0, rthread.reads.size() - count);
	    /*
	    *os << "Discarded" << std::endl;
	    PrintReads(rthread, os, rthread.reads.size() - count);
	    */
	}
#endif
	// Adjusting the offset incase the core is shifted
	// This happens when the baseread is filtered,
	// only the first iteration!!!
	
	if (length(core) == 0) {
	    /*	    
	    std::cerr << "*********************** READJUST *******************"
		      << std::endl;
	    */
	    int m = rthread.estAlign[0];
	    for (int i = 1; i < static_cast<int>(rthread.estAlign.size())
		     - count; ++i) {
		if (rthread.estAlign[i] > m) {
		    m = rthread.estAlign[i];
		}
	    }

	    
	    for (int i = 0; i < static_cast<int>(rthread.estAlign.size());
		 ++i) {
		rthread.estAlign[i] -=  m;
	    }

	    pos.offset += m;
	    pos.core_emission_offset =  pos.offset;

	    
	}
	
    }
    
    free(mask);
    
    return count;
}

void HMMCluster::GatherInitialReads(int cluster_id,
				    const DnaString& seed,
				    ReadThread& rthread,
				    Emission* buffer,
				    int buffer_iter,
				    int buffer_length) {
    if (static_cast<int>(length(seed)) < param->k) {
	return;
    }

    int idx[] = {0, length(seed)};
    finder.GetReads(cluster_id,
		    seed, idx, 2, 0, // so that we dont care extension
		    rthread,
		    main_pos,
		    this, buffer, buffer_iter, buffer_length);
}

void HMMCluster::GatherReads(int cluster_id,
			     const DnaString& core,
			     PositionInfo& pos,
			     int idx[], int idx_l,
			     ReadThread& rthread,
			     Emission* buffer,
			     int buffer_iter,
			     int buffer_length) {
    int t = (int) length(core);
    assert(t >= 0);
    
    finder.GetReads(cluster_id,
		    core, idx, idx_l, length(core),
		    rthread, pos,
		    this, buffer, buffer_iter, buffer_length);

    // PrintReads(rthread, &std::cerr);
}

#define MAX_ITER 100

void HMMCluster::EMLearning(gsl_rng *& r,
			    DnaString& core,
			    PositionInfo& pos,
			    Emission* buffer,
			    int& buffer_iter,
			    int buffer_length,
			    ReadThread& rthread,
			    bool fleft,
			    bool fright) {
    if (rthread.estAlign.size() <= 0) {
	std::cerr << "Empty list of reads" << std::endl;
	return;
    }
    
    gsl_permutation * p = gsl_permutation_alloc(rthread.estAlign.size());
    
    // learning!!!
    
    gsl_permutation_init(p);
    
    SEQAN_PROTIMESTART(hmmlearning);
    
    int mm = 0;
    int iter = 0;
    
    for (; iter < MAX_ITER; ++iter) {

	
	gsl_ran_shuffle(r, p->data, rthread.estAlign.size(), sizeof(size_t));
	
	// step-wise EM
	int iter_read = 0;
	
	while (iter_read < static_cast<int>(rthread.estAlign.size())) {
	    int pi = gsl_permutation_get(p, iter_read);
	    DnaString& read = rthread.RComplement[pi];

	    //std::cerr << rthread.reads[pi] << " " << read << std::endl;
	    ForwardBackward(buffer,
			    buffer_iter,
			    buffer_length,
			    read, rthread.estAlign[pi] + pos.offset);
	    
	    
#if DEBUG
	    std::cerr << "offset:" << pos.offset
		      << " core_emission_offset:" << pos.core_emission_offset
		      << std::endl;
#endif
	    
	    iter_read++;

	    // Updating	    
	    if (mm % 5 == 0) {
		// update
		// need to be smarter, not overupdating		    
		double diff = 0;
		int c = pos.emission_length - length(core);
		
		if (fleft) {
		    for (int j = 0; j < pos.core_emission_offset; ++j) {
			int idx = (j + buffer_iter) % buffer_length;
			diff += buffer[idx].getDiff();
		    }
		}
		if (fright) {
		    for (int j = pos.core_emission_offset + length(core);
			 j < pos.emission_length;
			 ++j) {
			int idx = (j + buffer_iter) % buffer_length;
			diff += buffer[idx].getDiff();
		    }
		}
		
		// std::cerr << mm << " " << diff << " " << c << std::endl;
		if (iter > 0 && diff / c < 0.01) {
		    iter = MAX_ITER;
		    break;
		}		
	    }
	    
	    mm++;
	    /*
	      if (mm % param->m == 0) {
	      trans.update();
	      }
	    */
	    
	}
	
    }
#if DEBUG	
    std::cerr << "Learning takes " << SEQAN_PROTIMEDIFF(hmmlearning)
	      << " seconds (" << mm << ")" << std::endl;
#endif
    gsl_permutation_free(p);
    
}

int HMMCluster::ExtendLeft(DnaString& core,
				PositionInfo& pos,
				Emission* buffer,
				int buffer_iter,
				int buffer_length,
				ostream* os) {
    // Pre-Core
    int o = pos.core_emission_offset - 1;
    double entropy;
    DnaString a;
    Dna s;

    /*
    if (os) {
	for (int i = 0; i < o; ++i) {
	    int c = buffer[(i + buffer_iter) % buffer_length].c[NDNA];
	    for (int k = 0; k < NDNA; ++k) {
		*os << exp(buffer[(i + buffer_iter) % buffer_length].s[k])
		    << " "
		    << (double) buffer[(i + buffer_iter) % buffer_length].c[k] / c
		    << " ";
	    }
	    *os << CheckDebugEntropyEmissions(buffer, (i + buffer_iter)
	           % buffer_length) << " ";
	}
    }
    */

    while (o >= 0 
	   && (entropy = buffer[(o + buffer_iter)
				% buffer_length].entropy)
	   <= param->entropy_th) {
	a += buffer[(o + buffer_iter) % buffer_length].letter;
	--o;pos.core_emission_offset--;
	/*
	if (os) {
	  // *os << s << "(" << entropy <<") ";
	  *os << "100\t" << entropy << "\t100" << std::endl;
	}
	*/
	
    }
    /*
    if (os) {
	*os << std::endl;
    }
    */
    DnaStringReverse b(a);
    DnaString c = b;
    append(c, core);
    core = c;

    return length(b);
}

int HMMCluster::ExtendRight(DnaString& core,
				 PositionInfo& pos,
				 Emission* buffer,
				 int buffer_iter,
				 int buffer_length,
				 ostream* os) {
double entropy;
    Dna s;
    int ext = 0;
    // Post-Core
    /*
    if (os) {
	for (int i = pos.core_emission_offset + length(core);
	     i < pos.emission_length; ++i) {
	    int c = buffer[(i + buffer_iter) % buffer_length].c[NDNA];
	    for (int k = 0; k < NDNA; ++k) {
		*os << exp(buffer[(i + buffer_iter) % buffer_length].s[k])
		    << " "
		    << (double) buffer[(i + buffer_iter) % buffer_length].c[k] / c
		    << " ";
	    }
	    *os << CheckDebugEntropyEmissions(buffer, (i + buffer_iter)
	           % buffer_length) << " ";
	}
    }
    */
    int i = pos.core_emission_offset + length(core);
    while (i < pos.emission_length
	   && (entropy = buffer[(i + buffer_iter)
				% buffer_length].entropy)
	   <= param->entropy_th) {
	core += buffer[(i + buffer_iter) % buffer_length].letter;
	i++;
	ext++;
	/*
	  if (os) {
	    // *os << s << "(" << entropy <<") ";
	    *os << "100\t" << entropy << "\t100" << std::endl;
	  }
	*/
    }
    /*
    if (os) {
	*os << std::endl;
    }
    */

    return ext;
}

int HMMCluster::AssignReads(int cluster_id,
			    DnaString& core,
			    PositionInfo& pos,
			    Emission* buffer,
			    int buffer_iter,
			    int buffer_length,
			    ReadThread& rthread,
			    Cluster& output,
			    ostream* os) {

     int stats_accepted = 0;
     bool collided = false;

     // Accepting reads that belongs to this clusters wholely	
     SEQAN_PROTIMESTART(viterbi);
     
     for (int i = 0; i < static_cast<int>(rthread.estAlign.size()); ++i) {
 #if DEBUG
	 std::cerr << rthread.reads[i] << "\t";
 #endif
	 vector<Error> corrections;
	 DnaString& read = rthread.RComplement[i];

	 DnaString c;
	 int sub_errors, del_errors, ins_errors;

	 // std::cerr << read << " ";
	 if (Viterbi(core, pos, buffer, buffer_iter,
		     buffer_length, read, rthread.estAlign[i]
		     + pos.offset, c, corrections, os,
		     sub_errors, ins_errors, del_errors, rthread.reads[i])
	     && (sub_errors + ins_errors + del_errors
		 < param->max_corrections_per_read)){
	     // This read belongs to this cluster!!!
	     // Saving later

	     if (stats_keeper.RemoveReads(cluster_id, rthread.reads[i],
					  c, corrections,
					  sub_errors + ins_errors + del_errors
					  > 0)
		 ) {
	       stats_accepted++;
	       g_sub_errors += sub_errors;
	       g_ins_errors += ins_errors;
	       g_del_errors += del_errors;
	     } else {
	       collided = true;
	     }
	     

	     // saving to the cluster
	     output.rthread.reads.push_back(rthread.reads[i]);
	     output.rthread.estAlign.push_back(rthread.estAlign[i]);
	     output.rthread.multi.push_back(rthread.multi[i]);

#if _DEBUG_XX_
	     if (stats_keeper.InterestingRead(rthread.reads[i])) {
		 interesting_c = true;
	     }
#endif
	 } else {
	     // reporting this read as failure
	     if (param->restrict_failures
		 && sub_errors + ins_errors + del_errors
		 < param->max_corrections_per_read) {
		 stats_keeper.ReportFailure(rthread.reads[i]);
	     }
	 }
     }

#if DEBUG
     std::cerr << "Viterbi takes " << SEQAN_PROTIMEDIFF(viterbi)
	       << " seconds." << std::endl;

     std::cerr << "Total: " << rthread.estAlign.size()
	       << " Accepted: " << stats_accepted;
     std::cerr << std::endl;
#endif

     g_stats_gathered += rthread.estAlign.size();
     g_stats_accepted += stats_accepted;

     return collided ? 0 : stats_accepted;
 }


 void HMMCluster::ResetStatistics(Emission* buffer,
				  int buffer_iter,
				  int buffer_length) {
     // reset statistics
     for (int i = 0; i < main_pos.core_emission_offset; ++i) {
	 buffer[i + buffer_iter].resetStats(param->pseudo);
     }

     for (int i = main_pos.core_emission_offset + length(main_core);
	  i < main_pos.emission_length; ++i) {
	 buffer[i + buffer_iter].resetStats(param->pseudo);
     }
 }

 void HMMCluster::ResetStatisticsLeftOnly() {
     for (int i = 0; i < left_pos.core_emission_offset; ++i) {
	 emission_left_buffer[(i + it_emission_left_buffer)
			      % length_emission_left_buffer]
	     .resetStats(param->pseudo);
     }
 }

 void HMMCluster::ResetStatisticsRightOnly() {
     for (int i = right_pos.core_emission_offset + length(right_core);
	  i < right_pos.emission_length; ++i) {
	 emission_right_buffer[(i + it_emission_right_buffer)
			       % length_emission_right_buffer]
	     .resetStats(param->pseudo);
     }
 }

bool HMMCluster::BuildCluster(const DnaString& baseread,
			      Cluster& output,
			      int readid,
			      const char* debugPath) {

     /*
     {
	 for (int i = 0 ; i < N_HSTATES; ++i) {
	    for (int j = 0 ; j < N_HSTATES; ++j) {
		std::cerr << exp(trans.Trans[i][j]) << "\t";
	    }
	    std::cerr << std::endl;
	}
	
     }
     */

     /* ***** */
     // Reset();

     stringstream* outfile = NULL;

     if (debugPath) {
	 outfile = new stringstream(stringstream::out);
     }

#if DEBUG
     std::cerr << "Base read: " << baseread << std::endl;
#endif

     ReadThread rthread;
     
     // Seed the reads with the baseread
     // HACK: need to be before offset/core_emission_offset is modified
     GatherInitialReads(readid, baseread, rthread,
			emission_buffer, it_emission_buffer,
			length_emission_buffer);
     rthread.GetReadString(fragStore);
     
     main_core = "";

     // initial init
     /* we shuould not bias to the baseread!!!
     PaddingBases(0, length(baseread));
     for (int i = 0; i < length(baseread); ++i) {
	 emissions[offset + i].v[baseread[i]] += 1;
	 emissions[offset + i].update();
     }
     */


     const gsl_rng_type * T;
     gsl_rng * r;
     gsl_rng_env_setup();
     T = gsl_rng_default;
     r = gsl_rng_alloc (T);

     int round_counter = 0;
     int extended_left_length = 0;
     int extended_right_length = 0;
     int assigned;

     bool split_extend = false;

     int nbases;

     do {
	 assigned = 0;
	 round_counter++;
	 assert(main_pos.core_emission_offset <= main_pos.offset);

 #if DEBUG
	 SEQAN_PROTIMESTART(round);
	 std::cerr << "Reads:" << rthread.reads.size() << std::endl;
	 // Output reads
	 PrintReads(rthread, outfile);
 #endif

	 for (int i = 0; i < static_cast<int>(rthread.estAlign.size()); ++i) {
	     int len = length(rthread.RComplement[i]);
	     PaddingBases(rthread.estAlign[i], len);
	 }

	 bool do_thinning = param->do_read_thinning
	     && (static_cast<int>(rthread.reads.size()) > param->rthgTh);

	 if (do_thinning) {
	     PopulateEmission(rthread, main_core, main_pos,
			      emission_buffer, it_emission_buffer,
			      length_emission_buffer,
			      true, true, outfile);
	     
	     int removed;
	     removed = FilterReads(rthread, main_core, main_pos,
				   emission_buffer, it_emission_buffer,
				   length_emission_buffer, outfile);
	     
	     if (removed > 0) {
		 int s = rthread.reads.size() - removed;
		 ModifyEmission(rthread, main_core, main_pos,
				    emission_buffer, it_emission_buffer,
				length_emission_buffer,
				true, true, outfile, s);
		 rthread.reads.resize(s);
		 rthread.estAlign.resize(s);
		 rthread.multi.resize(s);
		 rthread.RComplement.resize(s);
	     }
	     
	 }
	 
	 if (param->do_em) {
	     EMLearning(r, main_core, main_pos, emission_buffer, it_emission_buffer, length_emission_buffer,
			rthread, true, true);
	 } else if (!do_thinning) {
	     PopulateEmission(rthread, main_core, main_pos,
			      emission_buffer, it_emission_buffer, length_emission_buffer,
			      true, true, outfile);
	 }

	 //PrintEmissions(emission_buffer, it_emission_buffer, length_emission_buffer);

#if DEBUG
	 PrintEmissionEntropy();
 #endif

	 if (rthread.reads.size()  == 0) {
	     break;
	 }

	 // do not allow changing core
	 // CheckCoreCoherence();

	 extended_left_length =
	     ExtendLeft(main_core, main_pos, emission_buffer,
			it_emission_buffer, length_emission_buffer, outfile);
	 extended_right_length =
	     ExtendRight(main_core, main_pos, emission_buffer,
			 it_emission_buffer, length_emission_buffer, outfile);

	 if (outfile) {
	     *outfile << std::endl;
	 }

	 assigned = AssignReads(readid, main_core, main_pos,
				emission_buffer, it_emission_buffer,
				length_emission_buffer, rthread, output, outfile);

	 ResetStatistics(emission_buffer, it_emission_buffer,
			 length_emission_buffer);

#if DEBUG
	 if (outfile) {
	     *outfile << "Get more reads pre-core: "
		      << extended_left_length << std::endl;
	     *outfile << "Get more reads post-core: "
		      << extended_right_length << ","
		      << length(main_core) << std::endl;

	     *outfile << "============= CORE: "
		      << main_core << "(" << length(main_core)
		      << ")" << std::endl;
	     *outfile << "============= EMISSION: "
		      << main_pos.core_emission_offset
		      << " + " << length(main_core) 
		       << " + " << main_pos.emission_length
		 - (main_pos.core_emission_offset + length(main_core))
		      << " = " << main_pos.emission_length << std::endl;
	 }
 #endif
	 
	 SEQAN_PROTIMESTART(gatherRead);
	 // Gathering reads
	 rthread.clear();

	 if (assigned == 0) {
	   break;
	 }

	 if (static_cast<int>(length(main_core)) > extended_left_length
	     + extended_right_length + 2 * BUFFER_SEGMENT) {
	   split_extend = true;
	   break;
	 }

	 int idx[] = {0,
		      MIN(extended_left_length + param->k - 1,
			  static_cast<int>(length(main_core))),
		      MAX(0, static_cast<int>(length(main_core))
			  - extended_right_length - param->k + 1),
		      length(main_core)};

	 if (idx[2] <= idx[1]) {
	   idx[1] = idx[3];
	   GatherReads(readid, main_core, main_pos, idx, 2, rthread,
		       emission_buffer, it_emission_buffer,
		       length_emission_buffer);
	 } else {
	   GatherReads(readid, main_core, main_pos, idx, 4, rthread,
		       emission_buffer, it_emission_buffer,
		       length_emission_buffer);
	 }
	 rthread.GetReadString(fragStore);
	   
 #if DEBUG
	 std::cerr << "Gathering reads takes " << SEQAN_PROTIMEDIFF(gatherRead)
		   << " seconds." << std::endl;
	 std::cerr << "Round: " << round_counter << " takes "
		   << SEQAN_PROTIMEDIFF(round) << " seconds." << std::endl;
 #endif

     } while (rthread.reads.size() > 0); // end iteratively build the consensus



     if (split_extend) {

	 left_core = main_core;
	 right_core = main_core;
	 left_pos = main_pos;
	 right_pos = main_pos;

	 // Setting new buffer
	 // Have to set here before we used up the buffer
	 it_emission_left_buffer =
	     -(main_pos.core_emission_offset + extended_left_length
	       + BUFFER_SEGMENT) +  length_emission_left_buffer;
	 emission_left_buffer
	     = emission_buffer + it_emission_buffer - it_emission_left_buffer;

	 it_emission_right_buffer =
	     -( main_pos.core_emission_offset + length(main_core)
		- extended_right_length - BUFFER_SEGMENT);
	 emission_right_buffer
	     = emission_buffer + it_emission_buffer  - it_emission_right_buffer;

	 if (extended_left_length > 0) {
#if _DEBUG_XX_
	     if (outfile) {
		 *outfile << "LEFT EXTENSION!!!!" << std::endl;
		 }
#endif
	     SEQAN_PROTIMESTART(left_extension);
	     // Extending left only
	     rthread.clear();

	     int idx[] = {0, extended_left_length + param->k - 1};
	     GatherReads(readid, left_core, left_pos, idx, 2, rthread,
			 emission_left_buffer, it_emission_left_buffer,
			 length_emission_left_buffer);
	     rthread.GetReadString(fragStore);
	     
	     round_counter = 0;
	     while (rthread.reads.size() > 0) {
		 round_counter++;
		 assert(left_pos.core_emission_offset <= left_pos.offset);
		 
#if _DEBUG_XX_
		 SEQAN_PROTIMESTART(round);
		 if (outfile) {
		     *outfile << "Extending Left Reads:"
			      << rthread.reads.size() << std::endl;
		 }
		 PrintReads(rthread, outfile);
#endif
		 
		 
		 for (int i = 0; i < static_cast<int>(rthread.estAlign.size());
		      ++i) {
		     int len = length(rthread.RComplement[i]);
		     
		     PaddingBasesLeftOnly(rthread.estAlign[i], len);
		 }
		 
		 /*
		 PrintEmissions(emission_left_buffer,
				it_emission_left_buffer,
				length_emission_left_buffer);
		 */
		 bool do_thinning = param->do_read_thinning
		     && (static_cast<int>(rthread.reads.size())
			 > param->rthgTh);
		 
		 if (do_thinning) {
		     PopulateEmission(rthread, left_core, left_pos,
				      emission_left_buffer,
				      it_emission_left_buffer,
				      length_emission_left_buffer,
				      true, false, outfile);
		     int removed
			 = FilterReads(rthread, left_core, left_pos,
				       emission_left_buffer,
				       it_emission_left_buffer,
				       length_emission_left_buffer, outfile);
		     if (removed > 0) {
			 int s = rthread.reads.size() - removed;
			 ModifyEmission(rthread, left_core, left_pos,
					emission_left_buffer,
					it_emission_left_buffer,
					length_emission_left_buffer,
					true, false, outfile, s);
			 rthread.reads.resize(s);
			 rthread.estAlign.resize(s);
			 rthread.multi.resize(s);
			 rthread.RComplement.resize(s);
		     }
		     
		 }
		     
		 if (param->do_em) {
		     EMLearning(r, left_core, left_pos,
				emission_left_buffer, it_emission_left_buffer,
				length_emission_left_buffer,
				rthread, true, false);
		 } else if (!do_thinning) {
		     PopulateEmission(rthread, left_core, left_pos,
				      emission_left_buffer,
				      it_emission_left_buffer,
				      length_emission_left_buffer,
				      true, false, outfile);
		 }
		 
		 if (rthread.reads.size()  == 0) {
		     break;
		 }
		 
		 /*
		 PrintEmissions(emission_left_buffer,
				it_emission_left_buffer,
				length_emission_left_buffer);
		 */


		 
		 extended_left_length =
		     ExtendLeft(left_core, left_pos, emission_left_buffer,
				it_emission_left_buffer,
				length_emission_left_buffer, outfile);
		 if (outfile) {
		     *outfile << std::endl;
		 }
		 
		 int assigned =
		     AssignReads(readid, left_core, left_pos,
				 emission_left_buffer, it_emission_left_buffer,
				 length_emission_left_buffer, rthread,
				 output, outfile);
		 
		 // Reset
		 ResetStatisticsLeftOnly();
		 //  PrintEmissions(emission_left_buffer, it_emission_left_buffer,
		 //                 length_emission_left_buffer);
		 
 #if DEBUG
		 std::cerr << "Get more reads pre-core: "
			   << extended_left_length << std::endl;
		 
		 std::cerr << "============= CORE: "
			   << left_core << "(" << length(left_core)
			   << ")" << std::endl;
		 std::cerr << "============= EMISSION: "
			   << left_pos.core_emission_offset
			   << " + " << length(left_core) 
			   << " + " << left_pos.emission_length - (left_pos.core_emission_offset + length(left_core))
			   << " = " << left_pos.emission_length << std::endl;
#endif
		 
		 if (assigned == 0)
		     break;
		 
		 // Gathering reads
		 rthread.clear();
		 
		 assert(extended_left_length + param->k - 1 > 0);
		 
		 int idx[] = {0, extended_left_length + param->k - 1};
		 GatherReads(readid, left_core, left_pos, idx, 2, rthread,
			     emission_left_buffer, it_emission_left_buffer,
			     length_emission_left_buffer);
		 
		 rthread.GetReadString(fragStore);
#if DEBUG
		 std::cerr << "Round: " << round_counter
			   << " takes " << SEQAN_PROTIMEDIFF(round)
			   << " seconds." << std::endl;
#endif

	     } // end extending left
#if DEBUG
	     std::cerr << "Left extension takes " << SEQAN_PROTIMEDIFF(left_extension) << std::endl;
#endif
	 }
	 
	 if (extended_right_length > 0) {
#if _DEBUG_XX_
	     if (outfile) {
		 *outfile << "RIGHT EXTENSION!!!!" << std::endl;
	     }
#endif
	     SEQAN_PROTIMESTART(right_extension);
	     
	     rthread.clear();
	     
	     int idx[] = { MAX(0, ((int)length(right_core))
			       - extended_right_length - param->k + 1),
			  length(right_core)};
	     GatherReads(readid, right_core, right_pos,
			 idx, 2, rthread,
			 emission_right_buffer, it_emission_right_buffer,
			 length_emission_right_buffer);
	     
	     // PrintEmissions(emission_right_buffer, it_emission_right_buffer,
	     //		    length_emission_right_buffer);

	     rthread.GetReadString(fragStore);
	     
	     round_counter = 0;
	     while (rthread.reads.size() > 0) {
		 round_counter++;
		 assert(right_pos.core_emission_offset <= right_pos.offset);
		 
#if _DEBUG_XX_
		 SEQAN_PROTIMESTART(round);
		 if (outfile) {
		     *outfile << "Extending Right Reads:"
			      << rthread.reads.size() << std::endl;
		 }
		 PrintReads(rthread, outfile);
#endif


		 for (int i = 0; i < static_cast<int>(rthread.estAlign.size());
		      ++i) {
		     int len = length(rthread.RComplement[i]);
		     
		     PaddingBasesRightOnly(rthread.estAlign[i], len);
		 }
		 
		 bool do_thinning = param->do_read_thinning
		     && (static_cast<int>(rthread.reads.size())
			 > param->rthgTh);
		 /*	 
		 PrintEmissions(emission_right_buffer,
				it_emission_right_buffer,
				length_emission_right_buffer);
		 */
		 if (do_thinning) {
		     PopulateEmission(rthread, right_core, right_pos,
				      emission_right_buffer,
				      it_emission_right_buffer,
				      length_emission_right_buffer,
				      false, true, outfile);
		     int removed  =
			 FilterReads(rthread, right_core, right_pos,
				     emission_right_buffer,
				     it_emission_right_buffer,
				     length_emission_right_buffer, outfile);
		     if (removed > 0) {
			 int s = rthread.reads.size() - removed;
			 ModifyEmission(rthread, right_core, right_pos,
					emission_right_buffer,
					it_emission_right_buffer,
					length_emission_right_buffer,
					false, true, outfile, s);
			 rthread.reads.resize(s);
			 rthread.estAlign.resize(s);
			 rthread.multi.resize(s);
			 rthread.RComplement.resize(s);
		     }
		     
		 }
		 
		 if (param->do_em) {
		     EMLearning(r, right_core, right_pos,
				emission_right_buffer,
				it_emission_right_buffer,
				length_emission_right_buffer,
				rthread, false, true);
		 } else if (!do_thinning) {
		     PopulateEmission(rthread, right_core, right_pos,
				      emission_right_buffer,
				      it_emission_right_buffer,
				      length_emission_right_buffer,
				      false, true, outfile);
		 }
		 
		 if (rthread.reads.size()  == 0) {
		     break;
		 }
		 
		 // PrintEmissions(emission_right_buffer, it_emission_right_buffer,
		 //		length_emission_right_buffer);
		 

		 extended_right_length =
		     ExtendRight(right_core, right_pos,
				 emission_right_buffer,
				 it_emission_right_buffer,
				 length_emission_right_buffer, outfile);
		 if (outfile) {
		     *outfile << std::endl;
		 }
		 
		 int assigned =
		     AssignReads(readid, right_core, right_pos,
				 emission_right_buffer,
				 it_emission_right_buffer,
				 length_emission_right_buffer,
				 rthread, output, outfile);
		 
		     // Reset
		 ResetStatisticsRightOnly();
		 /*
		 PrintEmissions(emission_right_buffer,
				it_emission_right_buffer,
				length_emission_right_buffer);
		 */
#if DEBUG
		 std::cerr << "Get more reads post-core: "
			   << extended_right_length << ","
			   << length(right_core) << std::endl;
		 
		     std::cerr << "============= CORE: "
			       << right_core << "("
			       << length(right_core) << ")" << std::endl;
		     std::cerr << "============= EMISSION: "
			       << right_pos.core_emission_offset
			       << " + " << length(right_core) 
			       << " + " << right_pos.emission_length
			 - (right_pos.core_emission_offset + length(right_core))
			       << " = " << right_pos.emission_length << std::endl;
#endif

		     if (assigned == 0)
			 break;
		     
		     // Gathering reads
		     rthread.clear();
		     
		     assert(length(right_core) - extended_right_length
			    - param->k + 1 > 0);
		     
		     int idx[] = {((int) length(right_core))
				  - extended_right_length - param->k + 1,
				  length(right_core)};
		     GatherReads(readid, right_core, right_pos, idx, 2, rthread,
				 emission_right_buffer, it_emission_right_buffer,
				 length_emission_right_buffer);
		     rthread.GetReadString(fragStore);
#if DEBUG
		     std::cerr << "Round: " << round_counter
			       << " takes " << SEQAN_PROTIMEDIFF(round)
			       << " seconds." << std::endl;
 #endif

	     } // end extending right
#if DEBUG
	     std::cerr << "Right extension takes "
		       << SEQAN_PROTIMEDIFF(right_extension) << std::endl;
#endif
	 }
	 
	 nbases = ((int) length(left_core)) + length(right_core) - length(main_core);
	 // Saving the core
	 output.core = left_core;
	 append(output.core, infix(right_core, length(main_core), length(right_core)));
     } else {
	 nbases = length(main_core);
	 output.core = main_core;
     }


     /*
     std::cerr << "#PERF#: " << SEQAN_PROTIMEDIFF(assembly)
     << " " << nbases
	       << " " << g_stats_accepted
	       << " " << g_sub_errors
	       << " " << g_del_errors
	       << " " << g_ins_errors
	       << std::endl;
     */



     // free mem
     {
	 gsl_rng_free (r);
	 
	 //# if _DEBUG_XX_
	 if (outfile) {
	     if (interesting_c /* (double) g_stats_accepted / g_stats_gathered < 0.1*/) {
		 char fn2[256];
		 sprintf(fn2, "%s/db_%d_%d_%d_%d",
			 debugPath, nbases, g_stats_gathered,
			 g_stats_accepted, readid);
		 FILE* file = fopen(fn2, "w");
		 string str = outfile->str();
		 fprintf(file, "%s", str.c_str());
		 fclose(file);
	     }
	     delete outfile;
	 }
     }

     /*     
	    for (int i = 0 ; i < N_HSTATES; ++i) {
	    for (int j = 0 ; j < N_HSTATES; ++j) {
	    std::cerr << exp(trans.Trans[i][j]) << "\t";
	    }
	    std::cerr << std::endl;
	}
     */

     return interesting_c;
 }

double HMMCluster::GetEntropyEmissions(Emission* buffer, int i) {
    double e = 0;
    for (int j = 0; j < NDNA; ++j) {
	e += -buffer[i].s[j] * exp(buffer[i].s[j]);
     }
    
    return e;
}

double HMMCluster::CheckDebugEntropyEmissions(Emission* buffer, int i) {
    return buffer[i].entropy;
}

void HMMCluster::PrintEmissionEntropy() {
    std::cerr << " ******************************* " << std::endl;
    std::cerr << setprecision(3);
    for (int i = 0; i < main_pos.emission_length; ++i) {
	std::cerr << i << ":" << GetEntropyEmissions(emission_buffer, i + it_emission_buffer) << "\t";
    }
    std::cerr << "\n ******************************* " << std::endl;

}

 void HMMCluster::PrintEmissions(Emission* buffer,
				 int buffer_iter,
				 int buffer_length) {
     std::cerr << " ******************************* " << std::endl;
     std::cerr << setprecision(3);

     std::cerr << "*\t";
     for (int i = 0; i < buffer_length; ++i) {
	 std::cerr << i << "(" << (i + buffer_iter) % buffer_length << ")\t";
     }
     std::cerr << std::endl;

     for (int k = 0; k < NDNA; ++k) {
	 std::cerr << Dna(k) << "\t";
	 for (int i = 0; i < buffer_length; ++i) {
	     std::cerr << exp(buffer[(i + buffer_iter) % buffer_length].s[k]) << "\t";
	 }
	 std::cerr << std::endl;

     }
     std::cerr << " ******************************* " << std::endl;
 }

 void HMMCluster::PrintStats() {
     std::cerr << " Total Accepted: " << g_stats_accepted
	       << " Total Sub Errors: " << g_sub_errors
	       << " Total Ins Errors: " << g_ins_errors
	       << " Total Del Errors: " << g_del_errors
	       << std::endl;
 }

 void HMMCluster::PrintReads(const ReadThread& rthread,
			     ostream* os,
			     int start,
			     int end) const {
     if (start < 0) {
	 start = 0;
     }
     if ( end < 0) {
	 end = rthread.reads.size();
     }
#if _DEBUG_XX_
     if (os) {
	 int m = INT_MAX;
	 for (int i = start; i < end; ++i) {
	     if (rthread.estAlign[i] < m) {
		 m = rthread.estAlign[i];
	     }
	 }

	 for (int i = start; i < end; ++i) {
	     
	   DnaString read;
	   if (rthread.reads[i] >= 0) {
	     read = fragStore.readSeqStore[rthread.reads[i]];
	   } else {
	     read = fragStore.readSeqStore[-rthread.reads[i]-1];
	     reverseComplement(read);
	   }


	     for (int k = m; k < rthread.estAlign[i]; k++) {
		 *os << " ";
	     }
	     *os << read << " " << rthread.reads[i] << " "
		 << rthread.estAlign[i] << " " 
		 << rthread.multi[i] << std::endl;
	 }
     }
#endif
 }

// saving time by buffer this?
inline double MultiWeight(int m) {
    // std::cerr << 1.00 / (1 + exp( -(m-3))) << " " ;
    return 1.00 / (1 + exp( -2*(m-3)));
}

void HMMCluster::PopulateEmission(const ReadThread& rthread,
				  DnaString& core,
				  PositionInfo& pos,
				  Emission* buffer,
				  int& buffer_iter,
				  int buffer_length,
				  bool fleft,
				  bool fright,
				  ostream* os) {
    int min_pos = pos.emission_length;
    int max_pos = 0;
    for (int i = 0; i < static_cast<int>(rthread.reads.size()); ++i) {
	int alignment = rthread.estAlign[i] + pos.offset;
	const DnaString& read = rthread.RComplement[i];
	double v = MultiWeight(rthread.multi[i]);
	    
	for (int j = 0; j < static_cast<int>(length(read)); ++j) {
	    int idx = ToAllEmissionIndex(j + 1, alignment + INDELS);
	    min_pos = MIN(min_pos, idx);
	    max_pos = MAX(max_pos, idx);

	    idx = (idx + buffer_iter) % buffer_length;
	    buffer[idx].c[(int)read[j]] += v;
	    buffer[idx].c[NDNA]+=v;
	}
    }

    /* Update the estimation too */
    for (int j = min_pos; j <= max_pos; ++j) {
      int idx = (j + buffer_iter) % buffer_length;
      buffer[idx].estimateFromAlignment(param->pseudo);
    }
}


void HMMCluster::ModifyEmission(const ReadThread& rthread,
				DnaString& core,
				PositionInfo& pos,
				Emission* buffer,
				int& buffer_iter,
				int buffer_length,
				bool fleft,
				bool fright,
				ostream* os,
				int start,
				int end) {
    if (start < 0) {
	start = 0;
    }
    if ( end < 0) {
	end = rthread.reads.size();
    }
    
    int min_pos = pos.emission_length;
    int max_pos = 0;

    for (int i = start; i < end; ++i) {
        int alignment = rthread.estAlign[i] + pos.offset;
	const DnaString& read = rthread.RComplement[i];

	double v = MultiWeight(rthread.multi[i]);

	for (int j = 0; j < static_cast<int>(length(read)); ++j) {
	    int idx = ToAllEmissionIndex(j + 1, alignment + INDELS);
	    min_pos = MIN(min_pos, idx);
	    max_pos = MAX(max_pos, idx);

	    idx = (idx + buffer_iter) % buffer_length;
	    buffer[idx].c[(int)read[j]] -= v;
	    buffer[idx].c[NDNA] -= v;
	}
    }

    for (int j = min_pos; j <= max_pos; ++j) {
      int idx = (j + buffer_iter) % buffer_length;
      buffer[idx].estimateFromAlignment(param->pseudo);
    }



}

void HMMCluster::ReportPerf() {
    stats_keeper.UpdatePerf(g_sub_errors, g_ins_errors, g_del_errors);
}
