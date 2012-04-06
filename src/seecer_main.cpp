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

int correct_errors(int argc, char * argv[])
{
    HMMParameters param;
    Emission::m = 20;
    Emission::alpha = 0.5;

    param.m = Emission::m;
    param.alpha = Emission::alpha;

    param.pseudo = 1.005;
    
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

    param.k = 20;

    param.entropy_th = 1;
    param.cluster_llh_th = -1;
    param.do_em = false;
    param.restrict_failures = false;
    param.reuse_reads = true;
    param.do_read_thinning = false;
    param.rthgTh = 10;


    param.sim_coefficient = 10;
    param.sim_th = 1e-2;
    param.max_core_errors = 5;
    param.max_pre_corrected_errors = 10;
    param.max_corrections_per_read = 5;
    param.failure_th = 10;
    param.emit_delta = 0.00;

    char *kmerFn = NULL;
    char *debugPath = NULL;
    char *interestingFn = NULL;
    char *stats_suffix = NULL;

    const char *outputFn = static_cast<const char*>("corrected_reads.fa");
    int startSeed = -1;
    int endSeed = -1;

    /* Parsing options */
    int c;
    opterr = 0;

    while (1) {
	static struct option long_options[] =
	    {
		{"verbose", no_argument, 0, 'v'},
		{"doEM", no_argument, 0, 'm'},
		{"do-read-thinning", no_argument, 0, 't'},
		{"rthgTh", required_argument, 0, '6'},
		{"no-reuse-reads", no_argument, 0, 'r'},
		{"no-failure", no_argument, 0, 'f'},
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
		{0, 0, 0, 0}
	    };
	int option_index = 0;
	
	c = getopt_long (argc, argv, "vk:E:C:d:o:p:",
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
	    param.failure_th = atoi(optarg);
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
	    param.restrict_failures = true;
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
	case 'p':
	    omp_set_num_threads(atoi(optarg));
	    break;
	case 'S':
	    stats_suffix = optarg;
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
	printf("Too many arguments\n");
	return -1;
    }

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

    int n_clusters = 0;

    SEQAN_PROTIMESTART(total_correction);

    for (int rr = 0; rr < 1; rr++) {
#pragma omp parallel
    {
#pragma omp for
	for (int ridx = startSeed; ridx < endSeed; ++ridx) {
	    HMMCluster cluster(&param, *finder, stats_keeper,
			       finder->GetMaximumReadLength(),
			       Emission::alpha, Emission::m, fragStore);
	    if (!DiscardRead(fragStore.readSeqStore[ridx]) && stats_keeper.FreeRead(ridx)) {
		if (ridx % 1000 == 0) {
		    char otime[256];
		    time_t tim = time(NULL);
		    char *s = ctime_r(&tim, otime);
		    s[strlen(s) - 1] = '\0'; 
		    
		    std::cerr << "Baseread " << ridx << std::endl;
		    std::cerr << s << " " << tim << " Processed " << stats_keeper.NumProcessedReads() << " " << stats_keeper.NumCollidedReads() << " "
			      << n_clusters+skipped << "/" << endSeed - startSeed
			      << " = " << (double) (n_clusters+skipped) / (endSeed - startSeed) << " " << std::endl;
		    
		}
		
		Cluster res;
		
		SEQAN_PROTIMESTART(correction);
		
		    bool interesting_c =
			cluster.BuildCluster(fragStore.readSeqStore[ridx], res, ridx, debugPath);
		    cluster.ReportPerf();
		    
#if DEBUG
		    //std::cerr << "Time required for execution: " << SEQAN_PROTIMEDIFF(correction) << " seconds." << std::endl;
#endif
		    
#if DEBUG  
		    {
#pragma omp atomic
			n_clusters++;
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

#pragma omp atomic
		skipped++;
		}
	}
	
    }
    
    std::cerr << "============== Round " << rr << " ============"
	      << std::endl;
    stats_keeper.PrintStats();
    
    }



    std::cerr << "Total Reads: " << length(fragStore.readSeqStore)
	      << std::endl;
    std::cerr << "Skipped: " << skipped << std::endl;

    //cluster.PrintStats();

    std::cerr << "*** Total time required for execution: " << SEQAN_PROTIMEDIFF(total_correction) << " seconds." << std::endl;
    std::cerr << "*** " << n_clusters << " cores." << std::endl;
    
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
		output << ">" << fragStore.readNameStore[t] << "<#%" << std::endl;
		output << *s << std::endl;
	    } else {
		output << ">" << fragStore.readNameStore[t] << "<#" << std::endl;
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

    return 0;
}

int main(int argc, char * argv[]) {
    //omp_set_nested(1);
    //omp_set_dynamic(false);
    
    if (omp_get_nested()) {
	std::cerr << "NESTED THREADS supported" << std::endl;
    }
    
    correct_errors(argc, argv);
}
