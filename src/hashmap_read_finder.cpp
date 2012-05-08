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

#include "hashmap_read_finder.h"
#include "hmm_cluster.h"

QGramHashMapReadFinder::QGramHashMapReadFinder(THMMFragStore& fragStore,
					       StatsKeeper& stats_keeper,
					       int k)
    : fragStore(fragStore),
      stats_keeper(stats_keeper),
      rnaseq_k(k),
      readCount(length(fragStore.readSeqStore)),
      retrieved_count(static_cast<int*>(calloc(length(fragStore.readSeqStore), sizeof(int)))),
      extra_bits(-1),
      mask ((1ULL << rnaseq_k * BIT_SHIFT) - 1),
      max_read_length(0),
      grammap(NULL) {

    SEQAN_PROTIMESTART(build_index);
    BuildIndex();
    std::cerr << "Building Index in " << SEQAN_PROTIMEDIFF(build_index) << " seconds." << std::endl;
}

QGramHashMapReadFinder::~QGramHashMapReadFinder() {
    
    if (extra_bits > 0) {
	for (int ext = 0; ext < 1<<extra_bits; ++ext) {
	    
	  for (std::tr1::unordered_map<uint64_t, map_struct>::iterator it = grammap[ext].begin(); it != grammap[ext].end(); ++it) {
		delete[] it->second.ids;
		delete[] it->second.positions;
	    }
	}
    }

    free(retrieved_count);
    delete[] grammap;
}

void QGramHashMapReadFinder::UpdateGram(char letter, uint64_t& gram, uint64_t& r_gram) {

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

DnaString QGramHashMapReadFinder::GramToString(uint64_t gram) {
    DnaString t;

    uint64_t mask = (1 << BIT_SHIFT) - 1;

    for (int i = 0; i < rnaseq_k; ++i) {
	append(t, Dna(gram & mask));
	gram >>= BIT_SHIFT;
    }
    
    reverse(t);

    return t;
}

void QGramHashMapReadFinder::ConstructGramSet(const DnaString& s, std::set<uint64_t>& local_gram_set) {

    int lastN = -1;
    uint64_t gram = 0;
    uint64_t r_gram = 0;
    local_gram_set.clear();

    // std::cerr << "Contruct gram set: " << s;

    for (int k = 0; k < static_cast<int>(length(s)); ++k) {
	char letter = (int) s[k];
	
	if (letter == NDNA) {
	    lastN = k; 
	}
	
	UpdateGram(letter, gram, r_gram);
	
	if (lastN <= (k - rnaseq_k)) {
	    // add the gram
	    if (gram < r_gram) {
		local_gram_set.insert(gram);
	    } else {
		local_gram_set.insert(r_gram);
	    }
	}
	
    } // end getting all grams

    // std::cerr << " " << local_gram_set.size() << std::endl;

}

void QGramHashMapReadFinder::AddReadToIndex(ulong i, std::set<uint64_t>& local_gram_set) {

    int lastN = -1;
    uint64_t gram = 0;
    uint64_t r_gram = 0;
    local_gram_set.clear();
    int ext_bit = i >> 32;
    uint32_t idx = i & ((1UL << 32) - 1);

    // std::cerr << "Contruct gram set: " << s;

    for (int k = 0; k < static_cast<int>(length(fragStore.readSeqStore[i]));
	 ++k) {
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
#if DEBUG
		DnaString t = GramToString(lookup);
		
		std::cerr << "Gram " << lookup << " " << t;
		reverseComplement(t);
		std::cerr << " " << t << " " << fragStore.readSeqStore[i] << "(" << i << ")" << std::endl;
#endif
		//#pragma omp critical
		{
		grammap[ext_bit][lookup].ids[grammap[ext_bit][lookup].count] = idx;
		if (lookup == gram) {
		    grammap[ext_bit][lookup].positions[grammap[ext_bit][lookup].count++] = k + 1;
		} else {
		    grammap[ext_bit][lookup].positions[grammap[ext_bit][lookup].count++] = -(k + 1);
		}		 
		}

		local_gram_set.insert(lookup);
	    }
	
	} // end getting all grams

    }

}


void QGramHashMapReadFinder::BuildIndex() {
    
    // Counting the number of qgrams

    std::set<uint64_t> local_gram_set;

    extra_bits = MAX(ceil(log2(length(fragStore.readSeqStore))) - 32, 0);

    int gram_size = 0;

    std::cerr << "Extra # bits = " << extra_bits << std::endl;
    grammap = new std::tr1::unordered_map<uint64_t, map_struct>[1 << extra_bits];


    // TODO: This can be parallelized
    //#pragma omp parallel for
    for (ulong i = 0; i < length(fragStore.readSeqStore); ++i) {
	
	int ext_bit = i >> 32;

      if (i % 1000 == 0) {
	std::cerr << "Counting qgrams " << i << " reads" << std::endl;
      }

      if (!DiscardRead(fragStore.readSeqStore[i])) {
	  if (max_read_length < static_cast<int>(length(fragStore.readSeqStore[i]))) {
	      max_read_length = length(fragStore.readSeqStore[i]);
	  }
	  
	  ConstructGramSet(fragStore.readSeqStore[i], local_gram_set);
	  
	  // inserting into the hashmap
	  //#pragma omp critical
	  {
	      for (std::set<uint64_t>::const_iterator it = local_gram_set.begin(); it != local_gram_set.end(); ++it) {
		  if (grammap[ext_bit].find(*it) != grammap[ext_bit].end()) {
		      grammap[ext_bit][*it].count++;
		  } else {
		      grammap[ext_bit][*it].count = 1;
		  }
	      }
	      gram_size += local_gram_set.size();
	  }
      }
    }
    
    std::cerr << "Total " << gram_size << " grams " << std::endl;

    // allocating space

    std::cerr << "Allocating space" << std::endl;
    
    uint64_t total_count = 0;
    for (int ext = 0; ext < 1<<extra_bits; ++ext) {
      for (std::tr1::unordered_map<uint64_t, map_struct>::iterator it = grammap[ext].begin(); it != grammap[ext].end(); ++it) {
	    // std::cerr << "Gram: " << it->first << " " << it->second.count << std::endl;
	    it->second.ids = new uint32_t[it->second.count];
	    it->second.positions = new char[it->second.count];
	    total_count += it->second.count;

	    it->second.count = 0;
	}
    }
    // constructing the id lists

    std::cerr << total_count << " of indices" << std::endl;

    std::cerr << "Finding locations" << std::endl;
    // TODO: This can be parallelized    
    //#pragma omp parallel for
    for (ulong i = 0; i < length(fragStore.readSeqStore); ++i) {
      if (i % 1000 == 0) {
	std::cerr << "Getting qgrams locations of " << i << " reads" << std::endl;
      }

      if (!DiscardRead(fragStore.readSeqStore[i])) {
	  AddReadToIndex(i, local_gram_set);
      }
    }

}

void QGramHashMapReadFinder::GetReads(int cluster_id,
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

	int k = start;
	for (; k < start + rnaseq_k - 1; ++k) {
	    char letter = (int) core[k];
	    
	    UpdateGram(letter, gram, r_gram);
	}
	
	for (; k < end; ++k) {
	    assert(k < static_cast<int>(length(core)));
	    char letter = (int) core[k];
	    
	    UpdateGram(letter, gram, r_gram);
	    
	    uint64_t lookup = gram;
	    // looking up the hashmap
	    if (gram > r_gram) {
		lookup = r_gram;
	    }
	    
	    for (int ext = 0; ext < 1<<extra_bits; ++ext) {
		if (grammap[ext].find(lookup) == grammap[ext].end()) {
		    
		    // std::cerr << "Gram " << lookup << " " << GramToString(lookup) <<
		    // " found nothing!!!" << std::endl;
		    // assert(false);
		    continue;
		}
		
		map_struct& read_ids = grammap[ext][lookup];
		
		// std::cerr << "Looking up " << lookup << " " << GramToString(lookup) << " found " << read_ids.count << std::endl;
		
		for (ulong i = 0; i < read_ids.count; ++i) {
		    ulong id = read_ids.ids[i] | ((long) ext << 32);
		    
		    // std::cerr << "Finding overlap " << fragStore.readSeqStore[id] << "(" << id << ")" << std::endl;
		    
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
				int e = MIN(static_cast<int>(length(read)),
					    core_len - loc - INDELS);
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

void QGramHashMapReadFinder::PrintStats(const char* fn) {
    FILE* f = fopen(fn, "w");

    for (unsigned i = 0; i < length(fragStore.readSeqStore); ++i) {
	fprintf(f, "%d %d\n", i, retrieved_count[i]);
    }    

    fclose(f);
}
