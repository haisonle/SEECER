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
#include "common.h"

char Complement(char c) {
    c = toupper(c);

    switch (c) {
    case 'A': return 'T';
    case 'T': return 'A';
    case 'G': return 'C';
    case 'C': return 'G';
    default: return c;
    }

}

bool DiscardRead(const DnaString& read) {
  int nA = 0;
  int nT = 0;
  for (unsigned i = 0; i < length(read); ++i) {
    
    if ((int) read[i] == NDNA)
      return true;

    if ((int) read[i] == DNAA)
      ++nA;
    if ((int) read[i] == DNAT)
      ++nT;
  }

  return (nA > length(read) * 0.7) || (nT > length(read) * 0.7);
}

bool DiscardKmer(const char* read) {
  int nA = 0;
  int nT = 0;
  int n = 0;
  while (*read != '\0') {

    if (*read == 'A')
      ++nA;
    if (*read == 'T')
      ++nT;
    n++; read++;
  }

  return ((nA > n * 0.9) || (nT > n * 0.9));
}
