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

#ifndef READ_FINDER_H_
#define READ_FINDER_H_

#include "common.h"
#include "stats_keeper.h"

#include <fstream>

#include <set>
#include <vector>

#include <seqan/sequence.h>
#include <seqan/file.h>
#include <seqan/index.h>
#include <seqan/store.h>


using namespace seqan;

class HMMCluster;
struct Emission;

class ReadFinder {

public:
    virtual ~ReadFinder() {};

    virtual void GetReads(int cluster_id,
			  const DnaString& core,
			  int idx[], int idx_l,
			  int core_len,
			  ReadThread& rthread,
			  PositionInfo& pos,
			  HMMCluster* hmmcluster,
			  Emission* buffer,
			  int buffer_iter,
			  unsigned buffer_length) = 0;

    virtual int GetMaximumReadLength() = 0;

    virtual void PrintStats(const char* fn) {
    }

};

#endif
