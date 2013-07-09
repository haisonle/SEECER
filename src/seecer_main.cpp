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

//#define SEQAN_PROFILE		// enable time measuring


#define SEQAN_PARALLEL

#include <cstdio>
#include <fstream>
#include <getopt.h>
#include <iostream>
#include <set>
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <vector>

#include <omp.h>
#include <gsl/gsl_permutation.h>
#include <gsl/gsl_randist.h>

#include "common.h"
#include "hashmap_read_finder.h"
#include "hmm_cluster.h"
#include "smart_hashmap_read_finder.h"
#include "stats_keeper.h"


#ifndef DEBUG
#define DEBUG 0
#endif

using namespace seqan;

bool ReadFastaFiles(THMMFragStore& fragStore,
		    char* fn1, char* fn2) {
    char* fn1_i = fn1;
    char* fn2_i = fn2;
    bool br = false;
    while (!br) {
	while (*fn1 != ',' && *fn1 != '\0') {
	    fn1++;
	}
	if (*fn1 == '\0') {
	    br = true;
	} else {
	    *fn1 ='\0';
	}
	if (fn2) {
	    while (*fn2 != ',' && *fn2 != '\0') {
		fn2++;
	    }
	    *fn2 ='\0';
	    if (!loadReads(fragStore, fn1_i, fn2_i)) {
		return false;
	    }
	    fn2_i = ++fn2;
	} else {
	    if (!loadReads(fragStore, fn1_i)) {
		return false;
	    }
	}
	fn1_i = ++fn1;
    }

    std::cerr << "Total Reads " << length(fragStore.readSeqStore)
	      << std::endl;

    return true;
}


void help_msg() {
    const char* msg =
	"SEECER: SEquencing Error CorrEction for Rna-Seq data\n"
	"seecer [options] read1 [read2]\n"
	"--------------------------------------------\n"
	" read1, read2: are Fasta/Fastaq files.\n"
	"        If only read1 is provided, the reads are considered singles.\n"
        "        Otherwise, read1 and read2 are paired-end reads.\n"
        " *** Important ***:\n"
        " Reads should not contain Ns. Please use the provided run_seecer.sh \n"
        " script to handle Ns.\n"
	"--------------------------------------------\n"
	"Options:\n"
	" --kmer <k> : specify a different K value (default = 17).\n"
	" --kmerCount <f> : specify the file containing kmer counts. This file\n"
	"        is produced by JELLYFISH, we provided a Bash script to generate this file (run_seecer.sh).\n"
        "        If the parameter is not set, SEECER will count kmers by itself\n"
	"        (slower and memory-inefficient).\n"
	" --clusterLLH <e> : specify a different log likelihood threshold (default = -1).\n"
	" --entropy <e> : specify a different entropy threshold (default = 0.6).\n"
	" --help, -h : this help message.\n";
    fprintf(stderr, "%s\n", msg);
}

void writeContig(const HMMParameters& param, std::ofstream &f,
		 const DnaString& contig, Cluster& res) {
  if (f.is_open() && static_cast<int>(length(contig)) > param.k) {
    f << length(contig) << "\t" << contig << "\t" << res.rthread.reads.size() << "\t";

    for (std::deque<CoreCount>::iterator it = res.core_prob.begin();
	 it != res.core_prob.end(); ++it) {
      f << it->val[0] << "\t";
      f << it->val[1] << "\t";
      f << it->val[2] << "\t";
      f << it->val[3] << "\t";
      f << it->entropy << "\t";
    }
    
    f << std::endl;
  }
}

int correct_errors(int argc, char * argv[])
{
    HMMParameters param;
    Emission::m = 20;
    Emission::alpha = 0.5;

    param.m = Emission::m;
    param.alpha = Emission::alpha;

    param.pseudo = 0.01;
    
    param.Init[STATEM] = log(0.8);
    param.Init[STATEI] = log(0.1);
    param.Init[STATED] = log(0.1);
	    
    param.Trans[STATEM][STATEM] = log(0.94);
    param.Trans[STATEM][STATEI] = log(0.03);
    param.Trans[STATEM][STATED] = log(0.03);
    
    param.Trans[STATEI][STATEM] = log(0.95);
    param.Trans[STATEI][STATEI] = log(0.05);
    param.Trans[STATEI][STATED] = log(0.0);
	    
    param.Trans[STATED][STATEM] = log(0.95);
    param.Trans[STATED][STATEI] = log(0.0);
    param.Trans[STATED][STATED] = log(0.05);

    param.k = 17;

    param.entropy_th = 0.6;
    param.cluster_llh_th = -1;
    param.do_em = false;
    param.restrict_failures = true;
    param.reuse_reads = true;
    param.do_read_thinning = true;
    // Threadshold to do read-thinning
    param.rthgTh = 0;


    param.sim_coefficient = 10;
    param.sim_th = 1e-2;
    param.max_core_errors = 5;//5
    param.max_pre_corrected_errors = 10;//10
    param.max_corrections_per_read = 5;//5
    param.failure_th = 3;
    param.emit_delta = 0.1;

    const char *kmerFn = NULL;
    const char *debugPath = NULL;
    const char *interestingFn = NULL;
    const char *stats_suffix = NULL;
    const char *ctrlFn = NULL;
    const char *contigFn = NULL;
    const char *corPosFn = NULL;

    const char *outputFn = static_cast<const char*>("corrected_reads.fa");
    int startSeed = -1;
    int endSeed = -1;

    /* Parsing options */
    int c;
    opterr = 0;

    std::ofstream fContig;
    std::ofstream fCorPos;

    while (1) {
	static struct option long_options[] =
	    {
		{"verbose", no_argument, 0, 'v'},
		{"help", no_argument, 0, 'h'},
		{"doEM", no_argument, 0, 'm'},
		{"do-read-thinning", no_argument, 0, 't'},
		{"rthgTh", required_argument, 0, '6'},
		{"no-reuse-reads", no_argument, 0, 'r'},
		{"allow-failures", no_argument, 0, 'f'},
		{"kmerCount", required_argument, 0, 'k'},
		{"kmer", required_argument, 0, '5'},
		{"entropy", required_argument, 0, 'E'},
		{"clusterLLH", required_argument, 0, 'C'},
		{"debugPath", required_argument, 0, 'd'},
		{"output", required_argument, 0, 'o'},
		{"maxCorrections", required_argument, 0, 'M'},
		{"start", required_argument, 0, '0'},
		{"end", required_argument, 0, '9'},
		{"pseudocount", required_argument, 0, '1'},
		{"failureTh", required_argument, 0, '2'},
		{"eDelta", required_argument, 0, '3'},
		{"dbReads", required_argument, 0, '4'},
		{"stats_suffix", required_argument, 0, 'S'},
		{"ctrlFn", required_argument, 0, 'c'},
		{"contigFn", required_argument, 0, 'G'},
		{"corPosFn", required_argument, 0, 'P'},
		{0, 0, 0, 0}
	    };
	int option_index = 0;
	
	c = getopt_long (argc, argv, "vk:E:C:d:o:p:h",
			 long_options, &option_index);
	/* Detect the end of the options. */
           if (c == -1)
	       break;

	switch (c) {
	case 0:
	    /* If this option set a flag, do nothing else now. */
	    if (long_options[option_index].flag != 0)
		break;
	    printf ("option %s", long_options[option_index].name);
	    if (optarg)
		printf (" with arg %s", optarg);
	    printf ("\n");
	    break;
	case 'v':
	    break;
	case '0':
	    startSeed = atoi(optarg);
	    break;
	case '9':
	    endSeed = atoi(optarg);
	    break;
	case '1':
	    param.pseudo = atof(optarg);
	    break;
	case '2':
	    param.restrict_failures = true;
	    param.failure_th = MIN(0xff, atoi(optarg));
	    break;
	case '3':
	    param.emit_delta = atof(optarg);
	    break;
	case '4':
	    interestingFn = optarg;
	    break;
	case '5':
	    param.k = atoi(optarg);
	    break;
	case '6':
	    param.rthgTh = atoi(optarg);
	    break;
	case 'm':
	    param.do_em = true;
	    break;
	case 'M':
	    param.max_corrections_per_read = atoi(optarg);
	    break;
	case 't':
	    param.do_read_thinning = true;
	    break;
	case 'r':
	    param.reuse_reads = false;
	    break;
	case 'f':
	    param.restrict_failures = false;
	    break;
	case 'k':
	    kmerFn = optarg;
	    break;
	case 'E':
	    param.entropy_th = atof(optarg);
	    break;
	case 'C':
	    param.cluster_llh_th = atof(optarg);
	    break;
	case 'd':
	    debugPath = optarg;
	    break;
	case 'o':
	    outputFn = optarg;
	    break;
	case 'G':
	    contigFn = optarg;
	    fContig.open(contigFn);
	    if (!fContig.is_open()) {
	      std::cerr << "Config file is invalid!!!" << std::endl;
	      return -1;
	    }
	    break;
	case 'P':
	    corPosFn = optarg;
	    fCorPos.open(corPosFn);
	    if (!fCorPos.is_open()) {
	      std::cerr << "File to store corrected positions is invalid!!!" << std::endl;
	      return -1;
	    }
	    break;
	case 'c':
	    ctrlFn = optarg;
	    break;
	case 'p':
	    // Do not set more than 8
	    // omp_set_num_threads(MIN(8, atoi(optarg)));
	    omp_set_num_threads(atoi(optarg));
	    break;
	case 'S':
	    stats_suffix = optarg;
            break;
	case 'h':
	    help_msg();
	    exit(0);
	    break;
	default:
	    printf("Too many arguments\n");
	    return -1; //abort()
	}
	
    }

    // making sure the output is valid first
    std::ofstream output(outputFn);
    
    if (!output.is_open()) {
	std::cerr << "Output file is invalid!!!" << std::endl;
	return -1;
    }

    
    THMMFragStore fragStore;
    
    if (optind == argc - 1) {
	ReadFastaFiles(fragStore, argv[optind], 0);
    } else if (optind == argc - 2) {
	ReadFastaFiles(fragStore, argv[optind], argv[optind+1]);
    } else {
	// print help
	help_msg();
	return -1;
    }

#if 0

    std::ifstream res("corPos.txt");
    char s[256];
    while (!res.eof()) {
      res.getline(s, 256);
      int a,b;
      char c,d;
      
      sscanf(s, "%d\t%d\t%c\t%c", &a, &b, &c, &d);
      
      assert((unsigned int) a < length(fragStore.readSeqStore));

      switch (c) {
      case 'M':
	fprintf(stderr, "%d %d M %c%c\n", a, b, (char) fragStore.readSeqStore[a][b], d);
	break;
      case 'I':
	fprintf(stderr, "%d %d I %c\n", a, b, d);
	break;
      case 'D':
	fprintf(stderr, "%d %d D %c\n", a, b, (char) fragStore.readSeqStore[a][b-1]);
	break;
      default:
	break;
      }


      

    }
    


    res.close();

#endif

    StatsKeeper stats_keeper(&param, fragStore);

    ReadFinder* finder;

    if (kmerFn) {
	finder = new QGramSmartHashMapReadFinder(kmerFn, fragStore, stats_keeper, param.k);
    } else {
	finder = new QGramHashMapReadFinder(fragStore, stats_keeper, param.k);	
    }


#if DEBUG
    if (interestingFn) {
	stats_keeper.ReadIReads(interestingFn);
    }
#endif

    SEQAN_PROTIMESTART(constructQgramExt);

    std::cerr << "done. (" << SEQAN_PROTIMEDIFF(constructQgramExt) << " seconds)" << std::endl;

    String<DnaQ> seq;

    int skipped = 0;

    if (startSeed < 0) {
	startSeed = 0;
    }
    if (endSeed < 0) {
	endSeed = length(fragStore.readSeqStore);
    }

    std::cerr << "Total " << length(fragStore.readSeqStore) << " reads." << std::endl;


    int n_clusters = 0;

    SEQAN_PROTIMESTART(total_correction);

    bool abort = false;

#pragma omp parallel
    {
#pragma omp for
	for (int ridx = startSeed; ridx < endSeed; ++ridx) {
	    if (!abort) {
	      HMMCluster cluster(&param, *finder, stats_keeper,
				 finder->GetMaximumReadLength(),
				 Emission::alpha, Emission::m, fragStore);

	      if (corPosFn) {
		cluster.setCorrectedPoOStream(&fCorPos);
	      }

	      if (stats_keeper.FreeRead(ridx) && !DiscardRead(fragStore.readSeqStore[ridx])) {
		
		Cluster res;
		
		bool interesting_c =
		  cluster.BuildCluster(fragStore.readSeqStore[ridx], res, ridx, debugPath);
		cluster.ReportPerf();
#pragma omp critical
		{
		  writeContig(param, fContig, res.core, res);
		}
		
		n_clusters++;
		
		
		// dummy
		interesting_c = interesting_c;

#ifdef SEQAN_PROFILE
		SEQAN_PROTIMESTART(correction);
#endif
		    
#if DEBUG  
		    {
#pragma omp atomic
			// outputing core
			if (debugPath && res.rthread.reads.size() > 0
			    && interesting_c) {
			    char fn[256];
			    sprintf(fn, "%s/db_clust_%d", debugPath, ridx);
			    ofstream os(fn);
			    os << res.core << std::endl
			       << "--------------------------------"
			       << std::endl;
			    cluster.PrintReads(res.rthread, &os);
			    os.close();
			}
		    }
		
#endif
   
	    } else {
	        if (ridx % 1000 == 0) {
		    char otime[256];
		    time_t tim = time(NULL);
		    char *s = ctime_r(&tim, otime);
		    s[strlen(s) - 1] = '\0'; 
		    
		    // std::cerr << "Baseread " << ridx << std::endl;
		    std::cerr << s << " (" << tim << ") \n Assigned " 
			      << stats_keeper.NumProcessedReads() << " reads, " << stats_keeper.NumCollidedReads() << " collisions("
			      << stats_keeper.NumFailures() << "), ";
		    std::cerr << "Processed " << n_clusters+skipped << "/" << endSeed - startSeed
			      << " seeds = " << (double) (n_clusters+skipped) / (endSeed - startSeed) * 100
			      << " % complete" << std::endl;
		    
		    if (ctrlFn && fileExists(ctrlFn)) {
		      abort = true;
                      #pragma omp flush (abort)
		    }

		}
#pragma omp atomic
		skipped++;
	      }


	    }
	}
    }
    
    stats_keeper.PrintStats();

    std::cerr << "Total Reads: " << length(fragStore.readSeqStore)
	      << std::endl;

    //cluster.PrintStats();
#ifdef SEQAN_PROFILE
    std::cerr << "*** Total time required for execution: " <<
	SEQAN_PROTIMEDIFF(total_correction) << " seconds." << std::endl;
#endif

    std::cerr << "*** Total " << n_clusters << " contigs." << std::endl;
    
    /*
    for (unsigned long i = 0; i < all_cores.size(); ++i) {
	std::cerr << all_cores[i] << std::endl;
    }
    */

    const DnaString* s;
    
    for (unsigned long t = 0; t < length(fragStore.readSeqStore); ++t) {
	if (stats_keeper.ReadAccepted(t)) {
	    s = stats_keeper.GetReads(t);
	    if (s) {
		output << ">" << fragStore.readNameStore[t] << "<Corrected" << std::endl;
		output << *s << std::endl;
	    } else {
		output << ">" << fragStore.readNameStore[t] << "<Assigned" << std::endl;
		output << fragStore.readSeqStore[t] << std::endl;
	    }
	} else {
	    output << ">" << fragStore.readNameStore[t] << std::endl;
	    output << fragStore.readSeqStore[t] << std::endl;
	}
    }
    
    if (stats_suffix) {
	char fn[256];
	sprintf(fn, "nretrievals_%s", stats_suffix);
	finder->PrintStats(fn);
	sprintf(fn, "ncollisions_%s", stats_suffix);
	stats_keeper.PrintStats(fn);
    }
    
    delete finder;

    fContig.close();
    fCorPos.close();

    return 0;
}

int main(int argc, char * argv[]) {
    //omp_set_nested(1);
    //omp_set_dynamic(false);
    
    if (omp_get_nested()) {
	std::cerr << "NESTED THREADS supported" << std::endl;
    }
    omp_set_num_threads(MIN(8, omp_get_max_threads()));
    correct_errors(argc, argv);
}
