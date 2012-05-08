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

#include "stats_keeper.h"
#include <iostream>
#include <stdio.h>

#include <omp.h>

StatsKeeper::StatsKeeper(HMMParameters *param, THMMFragStore& fragStore)
    : param(param),
      fragStore(fragStore),
      readCount(length(fragStore.readSeqStore)),
      g_hs(static_cast<int*>(calloc(readCount, sizeof(int)))),
      g_failures(static_cast<unsigned char*>(calloc(readCount, sizeof(unsigned char)))),
      g_ncollisions(static_cast<int*>(calloc(readCount, sizeof(int)))),
      creads_index(static_cast<int*>(calloc(readCount, sizeof(int)))),
      nAcceptedReads(0),
      nCollidedReads(0),
      nCollidedMatchedReads(0),
      nCorrectedReads(0),
      g_sub_errors(0),
      g_ins_errors(0),
      g_del_errors(0) {
}

StatsKeeper::~StatsKeeper() {
    free(g_hs);
    free(g_failures);
}

void StatsKeeper::AcceptNewRead(bool corrected) {
#pragma omp atomic
    nAcceptedReads++;
    if (corrected) {
#pragma omp atomic
	nCorrectedReads++;
    }
    
}

void StatsKeeper::CollidedRead(bool matched) {
#pragma omp atomic
    nCollidedReads++;
    if (!matched) {
#pragma omp atomic
	nCollidedMatchedReads++;
    }
}

void StatsKeeper::PrintStats() const {
    std::cerr << "Accepted: " << nAcceptedReads
	      << " Corrected: " << nCorrectedReads
	      << std::endl
	      << "CollidedReads: " << nCollidedReads
	      << " BadCollided: " << nCollidedMatchedReads
	      << std::endl;
    std::cerr << "Substitutions: " << g_sub_errors
	      << " Deletions: " << g_del_errors
	      << " Insertions: " << g_ins_errors
	      << std::endl;
    std::cout << nAcceptedReads << " "
	      << nCorrectedReads << std::endl;
}

bool StatsKeeper::InterestingRead(int read_id) {
    if (read_id < 0) {
	read_id = - read_id - 1;
    }

    return interesting_reads.find(read_id) != interesting_reads.end();
}

bool StatsKeeper::RemoveReads(int cluster_id,
			      int read_index,
			      const DnaString& corrected,
			      std::vector<Error>& corrections,
			      bool fC) {

    bool ret = true;    
    DnaString c = corrected;
    if (read_index < 0 ) {
	read_index = -read_index - 1;
	reverseComplement(c);
    } else {
    }

#pragma omp critical (STATS_KEEPER)
{
    if (g_hs[read_index]) {
        bool matched;
	if (creads_index[read_index]) {
	    matched = corrected_reads[creads_index[read_index] - 1] == c;
	} else {
	    matched = fragStore.readSeqStore[read_index] == c;

	}
	g_ncollisions[read_index]++;
	CollidedRead(matched);
	ret = false;
    } else {

	g_hs[read_index] = cluster_id + 1;
	if (c != fragStore.readSeqStore[read_index]) {
	    corrected_reads.push_back(c);
	    creads_index[read_index] = corrected_reads.size();
	    AcceptNewRead(true);
	} else {
	    AcceptNewRead(false);
	}
    }
 }

 return ret;
}

void StatsKeeper::ReportFailure(int read_index) {
    read_index = read_index < 0 ? -read_index - 1 : read_index;

#pragma omp critical
    {
    if (g_failures[read_index] < 0xff) {
	g_failures[read_index]++;

    }
    }
}

bool StatsKeeper::ReadAccepted(int read_index) {
    return g_hs[read_index] != 0;
}

bool StatsKeeper::ReadUsefull(int cluster_id, int read_index) {
    bool ret;
#pragma omp critical (STATS_KEEPER)
    {
	ret = (g_failures[read_index] < param->failure_th) &&
	    param->reuse_reads ? g_hs[read_index] != cluster_id+1 : g_hs[read_index] == 0;
    }
    return ret;
}


bool StatsKeeper::FreeRead(int read_index) {
    bool ret;
#pragma omp critical (STATS_KEEPER)
{
    read_index = read_index < 0 ? -read_index - 1 : read_index;

    ret = g_hs[read_index] == 0
	&& g_failures[read_index] < param->failure_th;
}

    return ret;
}

void StatsKeeper::UpdatePerf(int sub_errors, int ins_errors, int del_errors) {
    g_sub_errors += sub_errors;
    g_ins_errors += ins_errors;
    g_del_errors += del_errors;
}

void StatsKeeper::ReadIReads(const char* fn) {

    FILE* f = fopen(fn, "r");
    int r = 0;
    while (!feof(f)) {
	fscanf(f, "%d\n", &r);
	interesting_reads.insert(r);
    }
    fclose(f);

}

const DnaString* StatsKeeper::GetReads(int read_index) {
    if (creads_index[read_index]) {
	return &corrected_reads[creads_index[read_index] - 1];
    } else {
        return NULL;
    }

}

void StatsKeeper::PrintStats(const char* fn) {
    FILE* f = fopen(fn, "w");

    for (unsigned i = 0; i < length(fragStore.readSeqStore); ++i) {
        fprintf(f, "%u %d\n", i, g_ncollisions[i]);
    }    

    fclose(f);
}
