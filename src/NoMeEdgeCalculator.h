/*Copyright 2017 Shounak Chakraborty
 * 
 * This file is part of ChromaClique.
 *
 * ChromaClique is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * ChromaClique is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with ChromaClique.  If not, see <http://www.gnu.org/licenses/>.
 *
 *
*/

#ifndef NOMEEDGECALCULATOR_H_
#define NOMEEDGECALCULATOR_H_


#include <set>

#include "AlignmentRecord.h"
#include "EdgeCalculator.h"

using namespace std;

class NoMeEdgeCalculator : public EdgeCalculator {
private:



   double nomeParam;
   std::unordered_map<int, double> switchRateFS;
   std::unordered_map<int, double> switchRateRS;



public:
   NoMeEdgeCalculator(double nomeParam, std::unordered_map<int, double> switchRateFS , std::unordered_map<int, double> switchRateRS);
   virtual ~NoMeEdgeCalculator();

   /** Decides whether an edge is to be drawn between the two given nodes. */
   virtual bool edgeBetween(const AlignmentRecord& ap1, const AlignmentRecord& ap2,int numGCAllowedPos,int ct) const;

   virtual void getPartnerLengthRange(const AlignmentRecord& ap, unsigned int* min, unsigned int* max) const;

   void setOverlapCliques(double d);
   double scoreSequences(string code1, std::string qual1, string code2, std::string qual2, int start1, int end1, int start2, int end2, bool isStrand1, bool isClique, int numGCAllowedPos, int ct) const;
   std::pair<std::string , int> fixIndels(string sequence,string reference) const ;
   std::pair<int , int> mixedReadDetermination(const AlignmentRecord & ap1, const AlignmentRecord & ap2) const;
   bool scoreReads(const AlignmentRecord& ap1, const AlignmentRecord& ap2, bool isClique, int numGCAllowedPos, int ct) const;
   bool determineStrand(const AlignmentRecord & ap) const;
   int pairedReadDetermination(const AlignmentRecord& ar1, const AlignmentRecord& ar2)const;

};
#endif /* NOMEEDGECALCULATOR_H_ */

