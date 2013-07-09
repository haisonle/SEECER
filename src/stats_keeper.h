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

#ifndef STATS_KEEPER_H_
#define STATS_KEEPER_H_

#include "common.h"

#include <ext/hash_map>
#include <set>
#include <vector>

#include <seqan/store.h>

using namespace seqan;

class StatsKeeper {

public:
    StatsKeeper(HMMParameters *param, THMMFragStore& fragStore);
    ~StatsKeeper();
    void AcceptNewRead(bool corrected);
    void CollidedRead(bool matched);
    void PrintStats() const;

    bool RemoveReads(int cluster_id,
		     int read_index,
		     const DnaString& corrected,
		     std::vector<Error>& corrections,
		     bool fC);
    void ReportFailure(int read_index);
    bool FreeRead(int read_index);

    bool ReadUsefull(int cluster_id, int read_index);
    bool ReadAccepted(int read_index);
    int GetClusterID(int read_index) {
	read_index = read_index >= readCount ? read_index - readCount : read_index;
	return g_hs[read_index];
    }

    int NumProcessedReads() {
      return nAcceptedReads;
    }

    int NumCollidedReads() {
      return nCollidedReads;
    }

    int NumFailures() {
      return g_nfailures;
    }

    void UpdatePerf(int sub_errors, int ins_errors, int del_errors);

    void ReadIReads(const char* fn);

    bool InterestingRead(int read_id);

    const DnaString* GetReads(int read_index);

    void PrintStats(const char* fn);
private:
    HMMParameters *param;
    THMMFragStore& fragStore;
    int readCount;
    int* g_hs;
    unsigned char* g_failures;
    int* g_ncollisions;
    int* creads_index;
    std::vector<DnaString> corrected_reads;

    unsigned int nAcceptedReads;
    unsigned int nCollidedReads;
    unsigned int nCollidedMatchedReads;
    unsigned int nCorrectedReads;
    unsigned int g_sub_errors;
    unsigned int g_ins_errors;
    unsigned int g_del_errors;

    
    unsigned int g_nfailures;
    // FOR DEBUG
    std::set<int> interesting_reads;
    
};

#endif
