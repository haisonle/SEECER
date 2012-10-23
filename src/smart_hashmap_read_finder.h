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


#ifndef SMART_HASHMAP_READ_FINDER_
#define SMART_HASHMAP_READ_FINDER_

#include <omp.h>

#include "common.h"
#include "read_finder.h"
#include <fstream>
#include <set>
#include <tr1/unordered_map>
#include <vector>

#define BIT_SHIFT 2

using namespace seqan;

typedef struct {
    int32_t max_count;
    int32_t count;
    uint32_t* ids;
    char* positions;
} s_map_struct;

class QGramSmartHashMapReadFinder : public ReadFinder {
public:
    QGramSmartHashMapReadFinder(const char* qgram_count_f,
				THMMFragStore& fragStore,
				StatsKeeper& stats_keeper,
				int k);
    virtual ~QGramSmartHashMapReadFinder();
    

    virtual void GetReads(int cluster_id,
			  const DnaString& core,
			  int idx[], int idx_l,
			  int core_len,
			  ReadThread& rthread,
			  PositionInfo& pos,
			  HMMCluster* hmmcluster,
			  Emission* buffer,
			  int buffer_iter,
			  unsigned buffer_length);

    virtual int GetMaximumReadLength() {
	return max_read_length;
    }



    virtual void PrintStats(const char* fn);

private:
    THMMFragStore& fragStore;
    StatsKeeper& stats_keeper;

    int rnaseq_k;

    ulong readCount;

    int* retrieved_count;

    int max_read_length;

    std::tr1::unordered_map<uint64_t, s_map_struct> grammap;

    uint64_t mask;
    
    void UpdateGram(char letter, uint64_t& gram, uint64_t& r_gram);
    void BuildIndex(const char* qgram_count_f);
    void AddReadToIndex(ulong i);
    DnaString GramToString(uint64_t gram);
    static const char codes[256];
    uint64_t GramToBinary(const char *in);
};

#endif
