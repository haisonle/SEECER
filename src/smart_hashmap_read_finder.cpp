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

#include "smart_hashmap_read_finder.h"
#include "hmm_cluster.h"

QGramSmartHashMapReadFinder::QGramSmartHashMapReadFinder(const char* qgram_count_f,
							 THMMFragStore& fragStore,
							 StatsKeeper& stats_keeper,
							 int k)
    : fragStore(fragStore),
      stats_keeper(stats_keeper),
      rnaseq_k(k),
      readCount(length(fragStore.readSeqStore)),
      retrieved_count(static_cast<int*>(calloc(readCount, sizeof(int)))),
      max_read_length(0),
      mask((1ULL << rnaseq_k * BIT_SHIFT) - 1) {
    SEQAN_PROTIMESTART(build_index);
    BuildIndex(qgram_count_f);
    std::cerr << "Building Index in " << SEQAN_PROTIMEDIFF(build_index) << " seconds." << std::endl;
}

QGramSmartHashMapReadFinder::~QGramSmartHashMapReadFinder() {
    
    
    for (std::tr1::unordered_map<uint64_t, s_map_struct>::iterator it = grammap.begin(); it != grammap.end(); ++it) {
	delete[] it->second.ids;
	delete[] it->second.positions;
    }

    free(retrieved_count);
}

void QGramSmartHashMapReadFinder::UpdateGram(char letter, uint64_t& gram, uint64_t& r_gram) {

    switch (letter) {
    case DNAA:
	gram = gram << BIT_SHIFT | DNAA;
	r_gram = r_gram >> BIT_SHIFT | ((uint64_t) DNAT << ((rnaseq_k - 1) * BIT_SHIFT ));
	break;
    case DNAT:
	gram = gram << BIT_SHIFT | DNAT;
	r_gram = r_gram >> BIT_SHIFT | ((uint64_t) DNAA << ((rnaseq_k - 1) * BIT_SHIFT ));
	break;
    case DNAG:
	gram = gram << BIT_SHIFT | DNAG;
	    r_gram = r_gram >> BIT_SHIFT | ((uint64_t) DNAC << ((rnaseq_k - 1) * BIT_SHIFT ));
	    break;
    case DNAC:
	gram = gram << BIT_SHIFT | DNAC;
	r_gram = r_gram >> BIT_SHIFT | ((uint64_t) DNAG << ((rnaseq_k - 1) * BIT_SHIFT ));
	break;
    default:
	gram = gram << BIT_SHIFT;
	r_gram = r_gram >> BIT_SHIFT;
	break;				    
    }

    gram &= mask;
    r_gram &= mask;
}

const char QGramSmartHashMapReadFinder::codes[256] = {
    -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -2, -3, -3, -3, -3, -3, 
    -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, 
    -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -1, -3, -3, 
    -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, 
    -3,  0, -1,  1, -1, -3, -3,  2, -1, -3, -3, -1, -3, -1, -1, -3, 
    -3, -3, -1, -1,  3, -3, -1, -1, -1, -1, -3, -3, -3, -3, -3, -3, 
    -3,  0, -1,  1, -1, -3, -3,  2, -1, -3, -3, -1, -3, -1, -1, -3, 
    -3, -3, -1, -1,  3, -3, -1, -1, -1, -1, -3, -3, -3, -3, -3, -3, 
    -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, 
    -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, 
    -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, 
    -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, 
    -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, 
    -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, 
    -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, 
    -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3
  };

uint64_t QGramSmartHashMapReadFinder::GramToBinary(const char *in) {
    uint64_t res = 0;
    for(int i = 0; i < rnaseq_k; i++) {
	const char c = codes[(int)*in++];
        if(c < 0)
	    return 0;
        res = (res << 2) | c;
    }
    return res;
}

DnaString QGramSmartHashMapReadFinder::GramToString(uint64_t gram) {
    DnaString t;

    uint64_t mask = (1 << BIT_SHIFT) - 1;

    for (int i = 0; i < rnaseq_k; ++i) {
	append(t, Dna(gram & mask));
	gram >>= BIT_SHIFT;
    }
    
    reverse(t);

    return t;
}

void QGramSmartHashMapReadFinder::AddReadToIndex(ulong i) {

    int lastN = -1;
    uint64_t gram = 0;
    uint64_t r_gram = 0;
    std::set<uint64_t> local_gram_set;

    // std::cerr << "Contruct gram set: " << s;

    for (int k = 0; k < static_cast<int>(length(fragStore.readSeqStore[i])); ++k) {
	char letter = (int) fragStore.readSeqStore[i][k];
	
	if (letter == NDNA) {
	    lastN = k; 
	}
	
	UpdateGram(letter, gram, r_gram);
	
	if (lastN <= (k - rnaseq_k)) {
	    // add the gram
	    uint64_t lookup = gram;
	    if (gram > r_gram) {
		lookup = r_gram;
	    }

	    if (local_gram_set.find(lookup) == local_gram_set.end()) {
		if(grammap.find(lookup) != grammap.end()) {

		    if (grammap[lookup].count >= grammap[lookup].max_count) {
			std::cerr << "ERROR!!!! " << GramToString(gram) << " "
				  << GramToString(r_gram) << " "
				  << grammap[lookup].count << std::endl;
			std::cerr << "Kmer count is not correct. Sorry we have to crash!"
				  << std::endl;
			exit(1);
		    }
		    
		    // Insert to the array
		    int idx = __sync_fetch_and_add(&grammap[lookup].count, 1);

		    assert(idx >=0);
		    grammap[lookup].ids[idx] = i;
		    if (lookup == gram) {
			grammap[lookup].positions[idx] = k + 1;
		    } else {
			grammap[lookup].positions[idx] = -(k + 1);
		    }		 
		}
		
		local_gram_set.insert(lookup);
	    }
	
	} // end getting all grams

    }

}


void QGramSmartHashMapReadFinder::BuildIndex(const char* qgram_count_f) {
    
    // Counting the number of qgrams

    int gram_size = 0;
    uint64_t total_count = 0;
    // allocating space

    std::cerr << "Allocating space" << std::endl;
    
    FILE* f = fopen(qgram_count_f, "r");
    char str[256];
    uint32_t count;

    while (!feof(f)) {
	if (gram_size % 100000 == 0) {
	    std::cerr << "Processed " << gram_size << " grams" << std::endl;
	}
	fscanf(f, "%s %u\n", str, &count);
	if (!DiscardKmer(str)) {
	  uint64_t hash = GramToBinary(str);
	  grammap[hash].max_count = count;
	  grammap[hash].count = 0;
	  grammap[hash].ids = new uint32_t[count];
	  grammap[hash].positions = new char[count];
	  gram_size++;
	  total_count += count;
	} else {
	  std::cerr << " SKIPPING " << str << std::endl;
	}
    }

    fclose(f);
    std::cerr << total_count << " of indices" << std::endl;

    std::cerr << "Total " << gram_size << " grams " << std::endl;
    std::cerr << total_count << " of indices" << std::endl;

    std::cerr << "Finding locations" << std::endl;

    // TODO: This can be parallelized    
#pragma omp parallel for
    for (int i = 0; i < static_cast<int>(length(fragStore.readSeqStore)); ++i) {
        if (i % 100000 == 0) {
	  std::cerr << "Getting qgrams locations of " << i << " reads" << std::endl;
        }

	if (!DiscardRead(fragStore.readSeqStore[i])) {
	    if (max_read_length < static_cast<int>(length(fragStore.readSeqStore[i]))) {
		max_read_length = static_cast<int>(length(fragStore.readSeqStore[i]));
	    }
	    
	    AddReadToIndex(i);
	}
      
    }

    uint64_t n_total_count = 0;
    for (std::tr1::unordered_map<uint64_t, s_map_struct>::iterator it = grammap.begin(); it != grammap.end(); ++it) {
	n_total_count += it->second.count;
    }
    if (n_total_count > total_count) {
      std::cerr << "Under-allocated memory!!!!" << std::endl;
      exit(1);
    }

    std::cerr << n_total_count << "/" << total_count <<" of indices" << std::endl;
}

void QGramSmartHashMapReadFinder::GetReads(int cluster_id,
					   const DnaString& core,
					   int idx[], int idx_l,
					   int core_len,
					   ReadThread& rthread,
					   PositionInfo& pos,
					   HMMCluster* hmmcluster,
					   Emission* buffer,
					   int buffer_iter,
					   unsigned buffer_length) {
    // std::cerr << "Get reads " << length(kmer) << std::endl;
    std::multiset<int> read_set;

    for (int i_idx = 0; i_idx < idx_l; i_idx+=2) {
	int start = idx[i_idx];
	int end = idx[i_idx + 1];
	uint64_t gram = 0;
	uint64_t r_gram = 0;
    
	start = MAX(0, start);
	end = MIN(end, static_cast<int>(length(core)));
	
	assert(static_cast<size_t>(core_len) <= length(core)
	       && static_cast<size_t>(end) <= length(core)
	       && static_cast<size_t>(start) < length(core)
	       && static_cast<size_t>(start) >= 0);
	
	int k = start;
	for (; k < start + rnaseq_k - 1; ++k) {
	    char letter = (int) core[k];
	    
	    UpdateGram(letter, gram, r_gram);
	}
	
	for (; k < end; ++k) {
	    assert(static_cast<size_t>(k) < length(core));
	    char letter = (int) core[k];
	
	    UpdateGram(letter, gram, r_gram);
	    
	    uint64_t lookup = gram;
	    // looking up the hashmap
	    if (gram > r_gram) {
		lookup = r_gram;
	    }
	    
	    
	    if (grammap.find(lookup) == grammap.end()) {
		
		// std::cerr << "Gram " << lookup << " " << GramToString(lookup) <<
		// " found nothing!!!" << std::endl;
		// assert(false);
		continue;
	    }
	
	    s_map_struct& read_ids = grammap[lookup];
	    
	    // std::cerr << "Looking up " << lookup << " " << GramToString(lookup) << " found " << read_ids.count << std::endl;
	    
	    for (int i = 0; i < read_ids.count; ++i) {
		int id = read_ids.ids[i];
	    
		//std::cerr << "Finding overlap " << fragStore.readSeqStore[id] << "(" << id << ")" << std::endl;
		int nid, loc;
		char read_k = abs(read_ids.positions[i]) - 1;
		
		if ((gram == lookup) == (read_ids.positions[i] > 0)) {
		    loc = k - read_k;
		    nid = id;
		    
		    // std::cerr << "Read " << fragStore.readSeqStore[id] << "(" << id << ") " 
		    //      << loc << std::endl;
		} else  {
		    loc = k - (rnaseq_k - 1) - (length(fragStore.readSeqStore[id]) - read_k - 1);
		    nid = -id - 1;
		    
		    // std::cerr << "Read " << fragStore.readSeqStore[id] << "(" << id << ") " 
		    //      << loc << std::endl;
		}
		
		if (stats_keeper.ReadUsefull(cluster_id, id)
		    && (loc < 0 
			|| (loc + static_cast<int>(length(fragStore.readSeqStore[id])) > core_len))
		    ) {
		    if (read_set.find(id) == read_set.end()) {
			bool overlapGood = true;
			
			if (core_len > 0) {
			    DnaString read = fragStore.readSeqStore[id];
			    if (static_cast<int>(id) != nid) {
				reverseComplement(read);
			    }
			    
			    int b = MAX(0, INDELS - loc);
			    int e = MIN(static_cast<int>(length(read)), core_len - loc - INDELS);
			    if (e > b) {
				read = infix(read, b, e);
				overlapGood = 
				hmmcluster->EstimateErrors(buffer, buffer_iter, buffer_length,
							   read, loc + b + pos.core_emission_offset);
			    }
			}
			
			if (overlapGood) {
			    retrieved_count[id]++;
			    rthread.reads.push_back(nid);
			    rthread.estAlign.push_back(loc + pos.core_emission_offset - pos.offset);
			    read_set.insert(id);
			}
		    } else {
			read_set.insert(id);
		    }
		}    
		
	    }
	} // end getting all grams
    }
    
    for (std::vector<int>::const_iterator it = rthread.reads.begin();
	 it != rthread.reads.end(); ++it) {
	if (*it > 0) {
	    rthread.multi.push_back(read_set.count(*it));
	} else {
	    rthread.multi.push_back(read_set.count(-*it - 1));
	}
    }
}

void QGramSmartHashMapReadFinder::PrintStats(const char* fn) {
    FILE* f = fopen(fn, "w");

    for (unsigned i = 0; i < length(fragStore.readSeqStore); ++i) {
	fprintf(f, "%u %d\n", i, retrieved_count[i]);
    }    

    fclose(f);
}
