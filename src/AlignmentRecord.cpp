
/* Copyright 2012 Tobias Marschall
 *
 * This file is part of HaploClique.
 *
 * HaploClique is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * HaploClique is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with HaploClique.  If not, see <http://www.gnu.org/licenses/>.
 */

#include <vector>
#include <algorithm>
#include <math.h>
#include <ctype.h>
#include <map>
#include <tuple>
#include <array>

#include <boost/tokenizer.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/algorithm/string/split.hpp>
#include <boost/algorithm/string/compare.hpp>

#include "AlignmentRecord.h"
#include "Clique.h"

using namespace std;
using namespace boost;

namespace {
//helper functions to merge DNA sequences
int agreement(const char& qual1, const char& qual2){
    float prob1 = std::pow(10,(float)-qual1/10);
    float prob2 = std::pow(10,(float)-qual2/10);
    float posterior = (prob1*prob2/3)/(1-prob1-prob2+4*prob1*prob2/3);
    posterior = std::round(-10*log10(posterior));
    return posterior;
}

int disagreement(const char& qual1, const char& qual2){
    float prob1 = std::pow(10,(float)-qual1/10);
    float prob2 = std::pow(10,(float)-qual2/10);
    float posterior = ((prob1*(1-prob2/3))/(prob1+prob2-4*prob1*prob2/3));
    posterior = std::round(-10*log10(posterior));
    return posterior;
}

float phredProb(const char& qual){
    return std::pow(10, (double)(-qual)/10.0);
}
//error probabilites are precomputed for all phred scores
std::array<float, 127> compute_error_probs(){
    std::array<float, 127> result;
    for (int i = 33; i < result.size(); i++){
        result[i] = phredProb(i-33);
    }
    return result;
}
//new values for updated error probabilites are precomputed
std::array<std::array<int, 127>, 127> compute_error_agreement(){
    std::array<std::array<int, 127>, 127> result;
    for (int i = 33; i < result.size(); i++){
        for(int j = 33; j < result.size(); j++){
            result[i][j] = agreement(i-33,j-33);
        }
    }
    return result;
}

std::array<std::array<int, 127>, 127> compute_error_disagreement(){
    std::array<std::array<int, 127>, 127> result;
    for (int i = 33; i < result.size(); i++){
        for(int j = 33; j < result.size(); j++){
            result[i][j] = disagreement(i-33,j-33);
        }
    }
    return result;
}

std::array<float, 127> error_probs = compute_error_probs();
std::array<std::array<int, 127>, 127> error_agreement = compute_error_agreement();
std::array<std::array<int, 127>, 127> error_disagreement = compute_error_disagreement();
}




int phred_sum(const string& phred, char phred_base=33) {
    int result = 0;
    for (size_t i=0; i<phred.size(); ++i) {
        result += phred[i] - phred_base;
    }
    return result;
}

//called by getmergedDnaSequence() to create new CigarData
std::vector<BamTools::CigarOp> createCigar(std::string& nucigar){
    std::vector<BamTools::CigarOp> res;
    unsigned int counter = 1;
    unsigned int pos = 0;
    char c;
    while(pos < nucigar.size()){
        c = nucigar[pos];
        pos++;
        while(pos < nucigar.size()){
            if (c == nucigar[pos]){
                counter++;
                pos++;
            } else {
                res.push_back(BamTools::CigarOp(c,counter));
                counter = 1;
                break;
            }
        }
        if(pos == nucigar.size()){
            res.push_back(BamTools::CigarOp(c,counter));
        }
    }
    return res;
}

int computeOffset(const std::vector<char>& cigar){

    int offset = 0;
    for(auto& i : cigar){

        if (i == 'S' || i == 'H'){
            offset++;
        } else break;
    }
    return offset;
}

void computeSOffset(const std::vector<char>& cigar,int& c, int& q){
    for(auto& i : cigar){
        if (i == 'S'){
            c++;
            q++;
        } else if(i == 'H'){

            c++;
        } else break;
    }
}

int computeRevOffset(const std::vector<char>& cigar){
    std::vector<char> t = cigar;
    int offset = 0;
    for (std::vector<char>::reverse_iterator it = t.rbegin(); it != t.rend(); ++it){
        if (*it == 'S' || *it == 'H'){
            offset++;
        } else break;
    }
    return offset;
}

int computeRevSOffset(const std::vector<char>& cigar){
    std::vector<char> t = cigar;
    int offset = 0;
    for (std::vector<char>::reverse_iterator it = t.rbegin(); it != t.rend(); ++it){
        if (*it == 'S'){
            offset++;
        } else break;
    }
    return offset;
}

std::pair<char,char> computeEntry(const char& base1, const char& qual1, const char& base2, const char& qual2){
    std::pair<char,char> result;

    if (base1==base2){
        result.first = base1;
        result.second = std::min(error_agreement[qual1][qual2]+33,126);
    }
    else if (qual1>=qual2) {
        result.first = base1;
        result.second = error_disagreement[qual1][qual2]+33;
    } else {
        result.first = base2;
        result.second = error_disagreement[qual2][qual1]+33;
    }
    return result;
}
std::vector<string> CTandCTOTReadCoder(string reference,string sequence,vector<int> qualityList,char char1,char char2,char char3){

    string codedRead = "";
    string codedReadQualities = "";
    bool endChecker;
    int counter = 0;

    for(int i =0; i<sequence.size()-2 ; i++){
        counter ++;
        if((sequence.at(i)==char1) && (sequence.at(i+1)==char2) && (sequence.at(i+2)!=char1)){
            if(reference.at(i+1)!=char2 || reference.at(i)!=char1){
                continue;
            }

            if(counter == 1){
                codedRead += "O";
                codedReadQualities += static_cast<char> (qualityList[i+1]);
            }else{
                codedRead += std::to_string(counter-1) + "O";
                codedReadQualities += std::to_string(counter-1);
                codedReadQualities += static_cast<char> (qualityList[i+1]);
            }
            i += 1;
            counter = 0;
        }else if((sequence.at(i)==char1) && (sequence.at(i+1)==char3) && (sequence.at(i+2)!=char1)){
            if(reference.at(i+1)!=char2 || reference.at(i)!=char1){
                continue;
            }
            if(counter == 1){
                codedRead += "C";
                codedReadQualities += static_cast<char> (qualityList[i+1]);
            }else{
                codedRead += std::to_string(counter-1) + "C";
                codedReadQualities += std::to_string(counter-1);
                codedReadQualities += static_cast<char> (qualityList[i+1]);
            }

            i += 1;
            counter = 0;
        }
        //This portion is specifically for calculating the last three nucleotides in the read sequence
        if(i>=sequence.size()-3){
            i++;
            if(i==sequence.size()-2){
                endChecker = true;
                if((sequence.at(i)==char1) && (sequence.at(i+1)==char2) && (reference.at(i+2)!=char1)){
                    counter++;
                    if(reference.at(i+1)!=char2 || reference.at(i)!=char1){
                        counter-=1;
                        break;
                    }
                    if(counter == 1){
                        codedRead += "O";
                        codedReadQualities += static_cast<char> (qualityList[i+1]);
                    }else{
                        codedRead += std::to_string(counter-1) + "O";
                        codedReadQualities += std::to_string(counter-1);
                        codedReadQualities += static_cast<char> (qualityList[i+1]);
                    }
                    i += 1;
                    counter = -1;

                }else if((sequence.at(i)==char1) && (sequence.at(i+1)==char3) && (reference.at(i+2)!=char1)){
                    counter++;
                    if(reference.at(i+1)!=char2|| reference.at(i)!=char1){
                        counter-=1;
                        break;
                    }
                    if(counter == 1){
                        codedRead += "C";
                        codedReadQualities += static_cast<char> (qualityList[i+1]);
                    }else{
                        codedRead += std::to_string(counter-1) + "C";
                        codedReadQualities += std::to_string(counter-1);
                        codedReadQualities += static_cast<char> (qualityList[i+1]);
                    }

                    i += 1;
                    counter = -1;
                }

            }else if(i==sequence.size()-1){
                endChecker = false;
            }
        }

    }
    if(counter!=0 && counter!=-1){
        codedRead += std::to_string(counter+2);
        codedReadQualities += std::to_string(counter+2);
    }else if(counter!=-1){
        if(endChecker){
            codedRead += std::to_string(counter+2);
            codedReadQualities += std::to_string(counter+2);
        }else{
            codedRead += std::to_string(counter+1);
            codedReadQualities += std::to_string(counter+1);
        }
    }

    return {codedRead,codedReadQualities};

}
//this coding needs to be updated and corrected as the other has been already done
std::vector<string> OBandCTOBReadCoder(string reference,string sequence,vector<int> qualityList,char char1,char char2,char char3){

    string codedRead = "";
    string codedReadQualities = "";

    int counter = 0;
    for(int i =0; i<sequence.size()-2 ; i++){
        counter ++;
        if((sequence.at(i)!=char1) && (sequence.at(i+1)==char2) && (sequence.at(i+2)==char1)){
            if(counter == 1){
                codedRead += "1O";
                codedReadQualities += "1";
                codedReadQualities += static_cast<char> (qualityList[i+1]);
            }else{
                codedRead += std::to_string(counter) + "O";
                codedReadQualities += std::to_string(counter);
                codedReadQualities += static_cast<char> (qualityList[i+1]);
            }
            i += 2;
            counter = 0;
        }else if((sequence.at(i)!=char1) && (sequence.at(i+1)==char3) && (sequence.at(i+2)==char1)){
            if(reference.at(i+1)==char3){
                //counter+=1;
                continue;
            }
            if(counter == 1){
                codedRead += "1C";
                codedReadQualities += "1";
                codedReadQualities += static_cast<char> (qualityList[i+1]);
            }else{
                codedRead += std::to_string(counter) + "C";
                codedReadQualities += std::to_string(counter);
                codedReadQualities += static_cast<char> (qualityList[i+1]);
            }

            i += 2;
            counter = 0;
        }

    }
    if(counter!=0){
        codedRead += std::to_string(counter+2);
        codedReadQualities += std::to_string(counter+2);

    }
    return {codedRead,codedReadQualities};

}

//Creates the align sequence and reference sequences and list containing the qualities
std::tuple<string,string,vector<int>> AlignmentRecord::referenceString(const std::vector<char> cigarData, int startPosition, std::string& sequence, string  qualities, string & reference){

    string constructedReference = "";
    string alignSeq="";
    vector<int> qualityList;
    int referenceCounter = startPosition-computeOffset(cigarData);
    int sequenceCounter=0;//this is the quality counter as well
    string replaceString = "G";

    for ( int s = 0; s < cigarData.size(); s++) {
        if (cigarData[s] == 'I') {
            alignSeq+=sequence.at(sequenceCounter);
            qualityList.push_back(qualities.at(sequenceCounter)-33);
            sequenceCounter++;
            constructedReference +='I';
        }else if (cigarData[s] == 'D'){
            alignSeq+='D';
            qualityList.push_back(0);
            constructedReference += (reference)[referenceCounter];
            referenceCounter++;
        }else if (cigarData[s] == 'S'){
            constructedReference += (reference)[referenceCounter];
            //Correcting GC position errors in the reads by checking the reference
            if(sequenceCounter!=0 && (sequence.at(sequenceCounter)=='C' || sequence.at(sequenceCounter)=='T') && (alignSeq.at(alignSeq.size()-1)!='G' && constructedReference.at(alignSeq.size()-1)=='G')){
                alignSeq = alignSeq.substr(0,alignSeq.size()-1);
                alignSeq += "G";
                alignSeq += sequence.at(sequenceCounter);
                replaceString += sequence.at(sequenceCounter);
                sequence.replace(sequenceCounter-1,2,replaceString);
                replaceString = "G";
            }else{
                alignSeq+=sequence.at(sequenceCounter);
            }
            qualityList.push_back(qualities.at(sequenceCounter)-33);
            referenceCounter++;
            sequenceCounter++;
        }else if (cigarData[s] == 'H'){
            alignSeq+='H';
            qualityList.push_back(qualities.at(sequenceCounter)-33);
            constructedReference += 'H';
            sequenceCounter++;
            referenceCounter++;
        }else if (cigarData[s] == 'M') {
            constructedReference += (reference)[referenceCounter];
            //Correcting GC position errors in the reads by checking the reference
            if(sequenceCounter!=0 && (sequence.at(sequenceCounter)=='C' || sequence.at(sequenceCounter)=='T') && (alignSeq.at(alignSeq.size()-1)!='G' && constructedReference.at(alignSeq.size()-1)=='G')){
                alignSeq = alignSeq.substr(0,alignSeq.size()-1);
                alignSeq += "G";
                alignSeq += sequence.at(sequenceCounter);
                replaceString += sequence.at(sequenceCounter);
                sequence.replace(sequenceCounter-1,2,replaceString);
                replaceString = "G";
            }else{
                alignSeq+=sequence.at(sequenceCounter);
            }
            qualityList.push_back(qualities.at(sequenceCounter)-33);
            sequenceCounter++;
            referenceCounter++;
        }else if (cigarData[s] == 'N') {
            constructedReference += (reference)[referenceCounter];
            //Correcting GC position errors in the reads by checking the reference
            if(sequenceCounter!=0 && (sequence.at(sequenceCounter)=='C' || sequence.at(sequenceCounter)=='T') && (alignSeq.at(alignSeq.size()-1)!='G' && constructedReference.at(alignSeq.size()-1)=='G')){
                alignSeq = alignSeq.substr(0,alignSeq.size()-1);
                alignSeq += "G";
                alignSeq += sequence.at(sequenceCounter);
                replaceString += sequence.at(sequenceCounter);
                sequence.replace(sequenceCounter-1,2,replaceString);
                replaceString = "G";
            }else{
                alignSeq+=sequence.at(sequenceCounter);
            }
            qualityList.push_back(qualities.at(sequenceCounter)-33);
            sequenceCounter++;
            referenceCounter++;
        }else if (cigarData[s] == 'P') {
            constructedReference += (reference)[referenceCounter];
            //Correcting GC position errors in the reads by checking the reference
            if(sequenceCounter!=0 && (sequence.at(sequenceCounter)=='C' || sequence.at(sequenceCounter)=='T') && (alignSeq.at(alignSeq.size()-1)!='G' && constructedReference.at(alignSeq.size()-1)=='G')){
                alignSeq = alignSeq.substr(0,alignSeq.size()-1);
                alignSeq += "G";
                alignSeq += sequence.at(sequenceCounter);
                replaceString += sequence.at(sequenceCounter);
                sequence.replace(sequenceCounter-1,2,replaceString);
                replaceString = "G";
            }else{
                alignSeq+=sequence.at(sequenceCounter);
            }
            qualityList.push_back(qualities.at(sequenceCounter)-33);
            sequenceCounter++;
            referenceCounter++;
        }else {
            assert(false);
        }
    }
    constructedReference+=reference[referenceCounter];


    if(cigarData.size()!=constructedReference.size()-1){
        assert(false);
    }
    boost::algorithm::to_upper(constructedReference);
    boost::algorithm::to_upper(alignSeq);
    return std::make_tuple(constructedReference,alignSeq,qualityList);

}




bool determineStrand(string sequence, string reference) {
    int flag=0;

    for (int i = 0; i < sequence.size(); i++) {
        if(sequence.at(i)=='A'){
            if(reference.at(i)=='G'){
                flag = 1;
                break;
            }
        }
    }
    if(flag==0){
        return true;
    }
    return false;
}

//Returns the complement of a given sequence
string complement(string sequence)  {
    string complementString = "";
    for(char base : sequence) {
        if(base=='A'){
            complementString += 'T';
        }else if(base=='T'){
            complementString += 'A';
        }else if(base=='G'){
            complementString += 'C';
        }else if(base=='C'){
            complementString += 'G';
        }else {
            complementString += base;
        }
    }

    return complementString;
}

bool is_number(const std::string& s)
{
    return !s.empty() && std::find_if(s.begin(),
        s.end(), [](char c) { return !std::isdigit(c); }) == s.end();
}



AlignmentRecord::AlignmentRecord(const BamTools::BamAlignment& alignment, int readRef, vector<string>* rnm , string &mainReference) : readNameMap(rnm) {

    this->single_end = true;
    this->readNames.insert(readRef);
    this->name = alignment.Name;
    this->start1 = alignment.Position + 1;
    this->end1 = alignment.GetEndPosition();
    this->cigar1 = alignment.CigarData;
    this->uniqueSupport = 0;
    this->support = 0;

    this->phred_sum1 = phred_sum(alignment.Qualities);
    //Constructing the unrolled cigar strings
    for (const auto& it : cigar1) {
        for (unsigned int s = 0; s < it.Length; ++s) {
            this->cigar1_unrolled.push_back(it.Type);
        }
        if (it.Type == 'D') {
            this->length_incl_deletions1+=it.Length;
            if (it.Length > 1) {
                this->length_incl_longdeletions1+=it.Length;
            }
        }
    }

    //NoMe Additions
    //Determines if the alignment is forward or reverse stranded
    if(!alignment.IsPaired()){
    this->strand1 = !alignment.IsReverseStrand();
    }else{
        if(alignment.IsFirstMate()){
            if(!alignment.IsReverseStrand()){
                this->strand1=true;
            }else{
                this->strand1=false;
            }
        }else if(alignment.IsSecondMate()){
            if(alignment.IsReverseStrand()){
                this->strand1=true;
            }else{
                this->strand1=false;
            }
        }
    }
    this->refId=alignment.RefID;
    string querySequence = alignment.QueryBases;
    //Reference string returns the reference, the aligned sequence and the quality list
    auto ref1 = referenceString(this->getCigar1Unrolled(),alignment.Position,querySequence,alignment.Qualities,mainReference);
    this->sequence1 = ShortDnaSequence(querySequence, alignment.Qualities);

    this->reference_seq1 = std::get<0>(ref1);
    this->alignSequence1=std::get<1>(ref1);
    this->qualityList1 = std::get<2>(ref1);
    this->length_incl_deletions1 = this->sequence1.size();
    this->length_incl_longdeletions1 = this->sequence1.size();

    //End Nome additions



    this->cov_pos = this->coveredPositions();
}

AlignmentRecord::AlignmentRecord(unique_ptr<vector<const AlignmentRecord*>>& alignments, unsigned int clique_id, std::map<string, int> *appearanceMap) : cigar1_unrolled(), cigar2_unrolled() {
    //no longer majority vote, phred scores are updated according to Edgar et al.
    assert ((*alignments).size()>1);
    //get first AlignmentRecord
    auto& al1 = (*alignments)[0];
    this->start1 = al1->getStart1();
    this->end1 = al1->getEnd1();
    this->cigar1 = al1->getCigar1();
    this->cigar1_unrolled = al1->getCigar1Unrolled();
    this->sequence1 = al1->getSequence1();
    this->readNameMap = al1->readNameMap;
    this->readNames.insert(al1->readNames.begin(), al1->readNames.end());
    this->single_end = al1->isSingleEnd();
    this->reference_seq1 = al1->getReferenceSeq1();
    this->alignSequence1 = al1->getAlignSequence1();
    this->qualityList1 = al1->getQualityList1();
    this->uniqueSupport = 0;
    this->support = 0;
    this->refId = al1->getRefId();

    this->start2 = 0;
    this->end2 = 0;

    //NoMe Additions
    this->strand1 = al1->isStrand1();
    // End NoMe Additions
    if(al1->isPairedEnd()){
        this->start2 = al1->getStart2();
        this->end2 = al1->getEnd2();
        this->cigar2 = al1->getCigar2();
        this->cigar2_unrolled = al1->getCigar2Unrolled();
        this->sequence2 = al1->getSequence2();
        this->reference_seq2 = al1->getReferenceSeq2();
        this->alignSequence2 = al1->getAlignSequence2();
        this->qualityList2 = al1->getQualityList2();
    }



    if(al1->getName().find("Clique")==std::string::npos){//meaning that its not a clique
        this->consistingReads.push_back(al1->getName());
    }else{
        this->consistingReads.push_back(al1->getName());
        for(auto s : al1->getConsistingReads()){
            this->consistingReads.push_back(s);
        }
        this->consistingReads.push_back("End");
    }
    if(appearanceMap->find(al1->getName())!=appearanceMap->end()){
        appearanceMap->operator[](al1->getName()) = appearanceMap->operator[](al1->getName()) + 1;
    }else{
        appearanceMap->operator[](al1->getName()) = 1;
    }

    //merge recent AlignmentRecord with all other alignments of Clique
    for (unsigned int i = 1; i < (*alignments).size(); i++){
        auto& al = (*alignments)[i];

        if(appearanceMap->find(al->getName())!=appearanceMap->end()){
            appearanceMap->operator [](al->getName()) = appearanceMap->operator [](al->getName()) + 1;
        }else{
            appearanceMap->operator [](al->getName()) = 1;
        }
        if(al1->getName().find("Clique") == std::string::npos){

            this->consistingReads.push_back(al->getName());
        }else{
            this->consistingReads.push_back(al->getName());
            for(auto s : al->getConsistingReads()){
                this->consistingReads.push_back(s);
            }
            this->consistingReads.push_back("End");
        }


        if (this->single_end && al->isSingleEnd()){
            mergeAlignmentRecordsSingle(*al,1,1);
        }
        else if (!(this->single_end) && al->isPairedEnd()){

            mergeAlignmentRecordsPaired(*al);
        }
        else {
            mergeAlignmentRecordsMixed(*al);
        }
        this->readNames.insert(al->readNames.begin(),al->readNames.end());


    }
    if(al1->getName().find("Clique") == std::string::npos){

        this->consistingReads.push_back("End");
    }
    //update name of new Clique Superread
    this->name = "Clique_" + to_string(clique_id);

    if(this->isSingleEnd()){

        if(this->isStrand1()){
            std::vector<string> codesAndQuality1 = CTandCTOTReadCoder(this->reference_seq1,this->alignSequence1,this->qualityList1,'G','C','T');
            this->code1 = codesAndQuality1[0];
            this->gcQuality1 = codesAndQuality1[1];
        }else{
            std::vector<string> codesAndQuality1 = OBandCTOBReadCoder(this->reference_seq1,this->alignSequence1,this->qualityList1,'G','C','T');
            this->code1 = codesAndQuality1[0];
            this->gcQuality1 = codesAndQuality1[1];
        }
    }else{
        if(this->isStrand1()){
            std::vector<string> codesAndQuality1 = CTandCTOTReadCoder(this->reference_seq1,this->alignSequence1,this->qualityList1,'G','C','T');
            this->code1 = codesAndQuality1[0];
            this->gcQuality1 = codesAndQuality1[1];
            std::vector<string> codesAndQuality2 = CTandCTOTReadCoder(this->reference_seq2,this->alignSequence2,this->qualityList2,'G','C','T');
            this->code2 = codesAndQuality2[0];
            this->gcQuality2 = codesAndQuality2[1];
        }else{
            std::vector<string> codesAndQuality1 = OBandCTOBReadCoder(this->reference_seq1,this->alignSequence1,this->qualityList1,'G','C','T');
            this->code1 = codesAndQuality1[0];
            this->gcQuality1 = codesAndQuality1[1];
            std::vector<string> codesAndQuality2 = OBandCTOBReadCoder(this->reference_seq2,this->alignSequence2,this->qualityList2,'G','C','T');
            this->code2 = codesAndQuality2[0];
            this->gcQuality2 = codesAndQuality2[1];
        }
    }
    this->cov_pos=this->coveredPositions();

}

//combines reads to paired end reads (single end given overlapping paired ends) when readBamFile is run
void AlignmentRecord::pairWith(const BamTools::BamAlignment& alignment, string &mainReference) {

    if ((unsigned)(alignment.Position+1) > this->end1) {
        this->single_end = false;
        this->start2 = alignment.Position + 1;
        this->end2 = alignment.GetEndPosition();
        if (!(this->end2 > 0)){
            //cout << this->name << " end: " << this->end2 << endl;
        }
        this->cigar2 = alignment.CigarData;

        this->phred_sum2 = phred_sum(alignment.Qualities);


        for (const auto& it : cigar2) {
            for (unsigned int s = 0; s < it.Length; ++s) {
                this->cigar2_unrolled.push_back(it.Type);
            }
            if (it.Type == 'D') {
                this->length_incl_deletions2+=it.Length;
                if (it.Length > 1) {
                    this->length_incl_longdeletions2+=it.Length;
                }
            }
        }
        //NoMe Additions
        string querySequence = alignment.QueryBases;
        std::tuple<string,string,vector<int>> ref2 = referenceString(this->getCigar2Unrolled(),alignment.Position,querySequence,alignment.Qualities,mainReference);

        this->sequence2 = ShortDnaSequence(querySequence, alignment.Qualities);

        this->reference_seq2=std::get<0>(ref2);
        this->alignSequence2=std::get<1>(ref2);
        this->qualityList2 = std::get<2>(ref2);

        this->length_incl_deletions2 = this->sequence2.size();
        this->length_incl_longdeletions2 = this->sequence2.size();
        //End NoMe Additions
        this->cov_pos = this->coveredPositions();

    } else if ((unsigned)alignment.GetEndPosition() < this->start1) {
        this->single_end = false;
        this->start2 = this->start1;
        this->end2 = this->end1;
        this->cigar2 = this->cigar1;
        this->sequence2 = this->sequence1;
        this->reference_seq2 = this->reference_seq1;
        this->alignSequence2 = this->alignSequence1;
        this->code2 = this->code1;
        this->qualityList2 = this->qualityList1;
        this->phred_sum2 = this->phred_sum1;
        this->cigar2_unrolled = this->cigar1_unrolled;
        this->length_incl_deletions2 = this->length_incl_deletions1;
        this->length_incl_longdeletions2 = this->length_incl_longdeletions1;

        this->start1 = alignment.Position + 1;
        this->end1 = alignment.GetEndPosition();
        this->cigar1 = alignment.CigarData;

        this->phred_sum1 = phred_sum(alignment.Qualities);

        this->cigar1_unrolled.clear();
        for (const auto& it : cigar1) {
            for (unsigned int s = 0; s < it.Length; ++s) {
                this->cigar1_unrolled.push_back(it.Type);
            }
            if (it.Type == 'D') {
                this->length_incl_deletions1+=it.Length;
                if (it.Length > 1) {
                    this->length_incl_longdeletions1+=it.Length;
                }
            }
        }

        //NoMe Additions
        string querySequence = alignment.QueryBases;
        std::tuple<string,string,vector<int>> ref3 = referenceString(this->getCigar1Unrolled(),alignment.Position,querySequence,alignment.Qualities,mainReference);

        this->sequence1 = ShortDnaSequence(querySequence, alignment.Qualities);
        this->reference_seq1 = std::get<0>(ref3);
        this->alignSequence1 = std::get<1>(ref3);
        this->qualityList1 = std::get<2>(ref3);

        this->length_incl_deletions1 = this->sequence1.size();
        this->length_incl_longdeletions1 = this->sequence1.size();

        //End NoMe Additions
        this->cov_pos = this->coveredPositions();
        if(this->reference_seq1.size()-1!=this->cigar1_unrolled.size()){
            assert(false);
        }
    }//merging of overlapping paired ends to single end reads
    else {
        this->getMergedDnaSequence(alignment,mainReference);
    }

    //Reverse Complementing the pair 1 and two based on which strand they come from
    if(this->isStrand1()){
            //OT and CTOT reads, no need to complement these reads as they are already in the forward direction
            std::vector<string> codesAndQuality1 = CTandCTOTReadCoder(this->reference_seq1,this->alignSequence1,this->qualityList1,'G','C','T');
            this->code1 = codesAndQuality1[0];
            this->gcQuality1 = codesAndQuality1[1];
            vector<string> codesAndQuality2 = CTandCTOTReadCoder(this->reference_seq2,this->alignSequence2,this->qualityList2,'G','C','T');
            this->code2 = codesAndQuality2[0];
            this->gcQuality2 = codesAndQuality2[1];

    }else if(!this->strand1){
            //OB and CTOB reads. We need to complement these reads since we want the original reverse strand reads
            this->sequence1 = ShortDnaSequence(complement(this->sequence1.toString()), this->sequence1.qualityString());
            this->reference_seq1 = complement(this->reference_seq1);
            this->alignSequence1 = complement(this->alignSequence1);
            this->sequence2 = ShortDnaSequence(complement(this->sequence2.toString()), this->sequence2.qualityString());
            this->reference_seq2 = complement(this->reference_seq2);
            this->alignSequence2 = complement(this->alignSequence2);
            vector<string> codesAndQuality1 = OBandCTOBReadCoder(this->reference_seq1,this->alignSequence1,this->qualityList1,'G','C','T');
            this->code1 = codesAndQuality1[0];
            this->gcQuality1 = codesAndQuality1[1];
            vector<string> codesAndQuality2 = OBandCTOBReadCoder(this->reference_seq2,this->alignSequence2,this->qualityList2,'G','C','T');
            this->code2 = codesAndQuality2[0];
            this->gcQuality2 = codesAndQuality2[1];

    }
}



//computes vector for AlignmentRecord which contains information about the mapping position in the reference, the base and its phred score, the error probability, the position of the base in the original read and an annotation in which read the base occurs given paired end reads.
std::vector<AlignmentRecord::mapValue> AlignmentRecord::coveredPositions() const{
    std::vector<AlignmentRecord::mapValue> cov_positions;
    //position in ref
    int r = this->start1;
    //position in querybases / quality string of read
    int q = 0;
    for (unsigned int i = 0; i< this->cigar1_unrolled.size(); ++i){
        char c = this->cigar1_unrolled[i];
        switch(c){
            case 'M': {
                c = this->sequence1[q];
                char qual = this->sequence1.qualityChar(q);
                cov_positions.push_back({r,c,qual,error_probs[qual],q,0});
                //char d = this->sequence1[q];
                ++q;
                ++r;
                break;
            }
            case 'D': {
                ++r;
                break;
            }
            case 'S':
            case 'I': {
                ++q;
                break;
            }
            case 'H':
                break;
        }
    }
    //In case of a paired end read
    if (!this->single_end){
        assert(this->start1 <= this->start2);
            //position in ref
            r = this->start2;
            //position in query bases
            q = 0;
            //in case of overlapping paired end it is declared as a single end
            //if (r <= this->end1){
            //   this->single_end = true;
            //}
            for (unsigned int i = 0; i< this->cigar2_unrolled.size(); ++i){
                char c = this->cigar2_unrolled[i];
                switch(c){
                    case 'M': {
                        c = this->sequence2[q];
                        char qual = this->sequence2.qualityChar(q);
                        cov_positions.push_back({r,c,qual,error_probs[qual],q,1});
                        ++q;
                        ++r;
                        break;
                    }
                    case 'D': {
                        ++r;
                        break;
                    }
                    case 'S':
                    case 'I': {
                        ++q;
                        break;
                    }
                    case 'H':
                        break;
              }
           }
        }
    return cov_positions;
}

//NO Overlap Merge only for Bam alignments
void AlignmentRecord::noOverlapMergeBAM(std::string& dna, std::string& qualities, std::string& nucigar, int& c_pos, int& q_pos, int& ref_pos, int i) const{
    //cout<<"noOverlapMerge"<<endl;
    char c;
    const ShortDnaSequence* s = 0;
    if (i == 1){
        c = this->cigar1_unrolled[c_pos];
        s = &this->sequence1;
    } else {
        c = this->cigar2_unrolled[c_pos];
        s = &this->sequence2;
    }
    if (c == 'H'){
        ref_pos++;
        c_pos++;
    } else if (c == 'I') {
        dna += (*s)[q_pos];
        qualities += s->qualityChar(q_pos);
        nucigar += 'I';
        q_pos++;
        c_pos++;
    } else if (c == 'D') {
        nucigar += 'D';
        ref_pos++;
        c_pos++;
    } else if (c == 'S'){
        ref_pos++;
        q_pos++;
        c_pos++;
    } else if (c == 'M'){
        dna += (*s)[q_pos];
        qualities += s->qualityChar(q_pos);
        nucigar += c;
        ref_pos++;
        q_pos++;
        c_pos++;
    } else {
        assert(false);

    }
}


//creates merged DNA sequences and Cigar string out of overlapping paired end reads while reading in BAM Files
void AlignmentRecord::getMergedDnaSequence(const BamTools::BamAlignment& alignment,string & mainReference){
        std::string dna = "";
        std::string qualities = "";
        std::string nucigar = "";
        std::vector<char> cigar_temp_unrolled;
        //vector of CigarOp
        for (const auto& it : alignment.CigarData) {
            for (unsigned int s = 0; s < it.Length; ++s) {
                cigar_temp_unrolled.push_back(it.Type);
            }
        }
        //get starting position and ending position according to ref position, paying attention to clipped bases
        int offset_f1 = computeOffset(this->cigar1_unrolled);
        int offset_f2 = computeOffset(cigar_temp_unrolled);
        int offset_b1 = computeRevOffset(this->cigar1_unrolled);
        int offset_b2 = computeRevOffset(cigar_temp_unrolled);

        //updated ref position including clips
        int ref_s_pos1 = this->start1-offset_f1;
        int ref_e_pos1 = this->end1+offset_b1;
        int ref_s_pos2 = alignment.Position+1-offset_f2;
        int ref_e_pos2 = alignment.GetEndPosition()+offset_b2;
        //position in query sequences // phred scores
        int q_pos1 = 0;
        int q_pos2 = 0;
        //position in unrolled cigar vectors
        int c_pos1 = 0;
        int c_pos2 = 0;
        //4 cases of different overlaps
        //------------
        //     ------------
        if(ref_s_pos1 <= ref_s_pos2 && ref_e_pos1 <= ref_e_pos2){
            while(ref_s_pos1<ref_s_pos2){
                noOverlapMergeBAM(dna,qualities,nucigar,c_pos1,q_pos1,ref_s_pos1,1);
            }
            while(ref_s_pos1<=ref_e_pos1){
                overlapMerge(alignment,dna,qualities,nucigar,cigar_temp_unrolled,c_pos1,c_pos2,q_pos1,q_pos2,ref_s_pos1);
            }
            while(ref_s_pos1<=ref_e_pos2){
                noOverlapMerge(alignment, dna, qualities, nucigar, cigar_temp_unrolled, c_pos2, q_pos2, ref_s_pos1);
            }
        }//------------------------------
            //           ----------
         else if (ref_s_pos1 <= ref_s_pos2 && ref_e_pos1 >= ref_e_pos2){
            while(ref_s_pos1<ref_s_pos2){
                noOverlapMergeBAM(dna,qualities,nucigar,c_pos1,q_pos1,ref_s_pos1,1);
            }
            while(ref_s_pos1<=ref_e_pos2){
                overlapMerge(alignment,dna,qualities,nucigar,cigar_temp_unrolled,c_pos1,c_pos2,q_pos1,q_pos2,ref_s_pos1);
            }
            while(ref_s_pos1<=ref_e_pos1){
                noOverlapMergeBAM(dna,qualities,nucigar,c_pos1,q_pos1,ref_s_pos1,1);
            }
            //         ----------
            //--------------------------
        } else if (ref_s_pos1 >= ref_s_pos2 && ref_e_pos1 <= ref_e_pos2){
            while(ref_s_pos2<ref_s_pos1){
                noOverlapMerge(alignment, dna, qualities, nucigar, cigar_temp_unrolled, c_pos2, q_pos2, ref_s_pos2);
            }
            while(ref_s_pos2<=ref_e_pos1){
                overlapMerge(alignment,dna,qualities,nucigar,cigar_temp_unrolled,c_pos1,c_pos2,q_pos1,q_pos2,ref_s_pos2);
            }
            while(ref_s_pos2<=ref_e_pos2){
                noOverlapMerge(alignment, dna, qualities, nucigar, cigar_temp_unrolled, c_pos2, q_pos2, ref_s_pos2);
            }
            //           --------------------
            //---------------------
        } else {
            assert(ref_s_pos1 >= ref_s_pos2 && ref_e_pos1 >= ref_e_pos2);
            while(ref_s_pos2<ref_s_pos1){
                noOverlapMerge(alignment, dna, qualities, nucigar, cigar_temp_unrolled, c_pos2, q_pos2, ref_s_pos2);
            }
            while(ref_s_pos2<=ref_e_pos2){
                overlapMerge(alignment,dna,qualities,nucigar,cigar_temp_unrolled,c_pos1,c_pos2,q_pos1,q_pos2,ref_s_pos2);
            }
            while(ref_s_pos2<=ref_e_pos1){
                noOverlapMergeBAM(dna,qualities,nucigar,c_pos1,q_pos1,ref_s_pos2,1);
            }
        }
        this->start1 = std::min(this->start1,(unsigned int)alignment.Position+1);
        this->end1=std::max((unsigned int)alignment.GetEndPosition(),this->end1);
        this->single_end= true;
        this->cigar1 = createCigar(nucigar);

        this->phred_sum1=phred_sum(qualities);

        this->cigar1_unrolled.clear();
        for (char i : nucigar){
            this->cigar1_unrolled.push_back(i);
        }

        //NoMe Additions
        std::tuple<string,string,vector<int>> ref = referenceString(this->getCigar1Unrolled(),alignment.Position,dna,qualities,mainReference);
        if(!this->strand1){
           this->sequence1 = ShortDnaSequence(complement(dna), qualities);
           this->reference_seq1=complement(std::get<0>(ref));
           this->alignSequence1=complement(std::get<1>(ref));
            this->qualityList1 = std::get<2>(ref);
        }else{
           this->sequence1=ShortDnaSequence(dna,qualities);
           this->reference_seq1=std::get<0>(ref);
           this->alignSequence1=std::get<1>(ref);
            this->qualityList1 = std::get<2>(ref);
        }
        //fill out thr frs ffs thing here

        this->length_incl_deletions1 = this->sequence1.size();
        this->length_incl_longdeletions1 = this->sequence1.size();
        //End NoMe Additions
        //this->strandInfo = this->determineStrand(this->alignSequence1,this->reference_seq1);
        this->cov_pos = this->coveredPositions();
}

//helper function for getMergedDnaSequence
void AlignmentRecord::noOverlapMerge(const BamTools::BamAlignment& alignment, std::string& dna, std::string& qualities, std::string& nucigar, std::vector<char>& cigar_temp_unrolled, int& c_pos, int& q_pos, int& ref_pos) const{
    char c = cigar_temp_unrolled[c_pos];
    if (c == 'H'){
        ref_pos++;
        c_pos++;
    } else if (c == 'I') {
        dna += alignment.QueryBases[q_pos];
        qualities += alignment.Qualities[q_pos];
        nucigar += 'I';
        q_pos++;
        c_pos++;
    } else if (c == 'D') {
        nucigar += 'D';
        ref_pos++;
        c_pos++;
    } else if (c == 'S'){
        ref_pos++;
        q_pos++;
        c_pos++;
    } else if (c == 'M'){
        dna += alignment.QueryBases[q_pos];
        qualities += alignment.Qualities[q_pos];
        nucigar += c;
        ref_pos++;
        q_pos++;
        c_pos++;
    } else {
        assert(false);
    }
}

//helper function for getMergedDnaSequence, clipped bases are NOT contained in final sequence
void AlignmentRecord::overlapMerge(const BamTools::BamAlignment& alignment, std::string& dna, std::string& qualities, std::string& nucigar, std::vector<char>& cigar_temp_unrolled, int& c_pos1, int& c_pos2, int& q_pos1, int& q_pos2, int& ref_pos) const{
    char c1 = this->cigar1_unrolled[c_pos1];
    char c2 = cigar_temp_unrolled[c_pos2];
    if((c1 == 'M' && c2 == 'M') || (c1 == 'S' && c2 == 'S') || (c1 == 'I' && c2 == 'I')){
        if (c1 != 'S'){
            std::pair<char,char> resPair = computeEntry(this->sequence1[q_pos1],this->sequence1.qualityChar(q_pos1),alignment.QueryBases[q_pos2],alignment.Qualities[q_pos2]);
            dna += resPair.first;
            qualities += resPair.second;
            nucigar += c1;
        }
        if (c1 != 'I') ref_pos++;
        q_pos1++;
        q_pos2++;
        c_pos1++;
        c_pos2++;
    } else if ((c1 == 'D' && c2 == 'D') || (c1 == 'H' && c2 == 'H') || (c1 == 'D' && c2 == 'H') || (c1 == 'H' && c2 == 'D') || (c1 == 'D' && c2 == 'S') || (c1 == 'S' && c2 == 'D')){
        c_pos1++;
        c_pos2++;
        ref_pos++;
        if (c1 == 'D' || c2 == 'D'){
            nucigar += 'D';
        }
        if (c1 == 'S'){
            q_pos1++;
        } else if(c2 == 'S'){
            q_pos2++;
        }
    } else if ((c1 == 'M' && (c2 == 'D' || c2 == 'H' || c2 == 'S')) || ((c1 == 'D' || c1 == 'H' || c1 == 'S') && c2 == 'M') || (c1 == 'S' && c2 == 'H') || (c1 == 'H' && c2 == 'S')) {
        if (c1 == 'M'){
            nucigar += 'M';
            dna += this->sequence1[q_pos1];
            qualities += this->sequence1.qualityChar(q_pos1);
            ref_pos++;
            c_pos1++;
            c_pos2++;
            q_pos1++;
            if (c2 == 'S') q_pos2++;
        } else if (c2 == 'M'){
            nucigar += 'M';
            dna +=  alignment.QueryBases[q_pos2];
            qualities += alignment.Qualities[q_pos2];
            ref_pos++;
            c_pos1++;
            c_pos2++;
            q_pos2++;
            if (c1 == 'S') q_pos1++;
        } else if (c1 == 'S'){
            ref_pos++;
            c_pos1++;
            c_pos2++;
            q_pos1++;
        } else {
            ref_pos++;
            c_pos1++;
            c_pos2++;
            q_pos2++;
        }
    } else if (c1 == 'I' || c2 == 'I'){
        if(c1 == 'I'){
            nucigar += 'I';
            dna += this->sequence1[q_pos1];
            qualities += this->sequence1.qualityChar(q_pos1);
            c_pos1++;
            q_pos1++;
        } else {
            nucigar += 'I';
            dna +=  alignment.QueryBases[q_pos2];
            qualities += alignment.Qualities[q_pos2];
            c_pos2++;
            q_pos2++;
        }
    } else {
        assert(false);
    }
}

void AlignmentRecord::updateReference(std::string& ref,int i,int refFinalPos) const{
    std::string r = "";
    if (i == 1){
        r = this->getReferenceSeq1();
    } else {
        r = this->getReferenceSeq2();
    }
    ref = ref + r.at(refFinalPos);
}

//method is partly also used by mergeAlignmentRecordsMixed and mergeAlignmentRecordsPaired, i = ith cigar of this, j = ith cigar of AlignmentRecord ar
void AlignmentRecord::mergeAlignmentRecordsSingle(const AlignmentRecord& ar, int i, int j){
    std::string dna = "";
    std::string ref = "";
    std::string algn = "";
    std::string qualities = "";
    std::string nucigar = "";
    int offset_f1, offset_f2, offset_b1, offset_b2, ref_s_pos1, ref_e_pos1, ref_s_pos2,ref_e_pos2 = 0;
    //i and j determine which cigar strings / sequences are considered

    if (i == 1){
         //get starting position and ending position according to ref position, paying attention to clipped bases
         //updated ref position including clips
         offset_f1 = computeOffset(this->cigar1_unrolled);
         offset_b1 = computeRevOffset(this->cigar1_unrolled);
         ref_e_pos1 = this->end1+offset_b1;
         ref_s_pos1 = this->start1-offset_f1;
         if(j == 1){
             std::vector<char> cigar = ar.getCigar1Unrolled();
             offset_f2 = computeOffset(cigar);
             offset_b2 = computeRevOffset(cigar);
             ref_s_pos2 = ar.getStart1()-offset_f2;
             ref_e_pos2 = ar.getEnd1()+offset_b2;
         } else {
             std::vector<char> cigar = ar.getCigar2Unrolled();
             offset_f2 = computeOffset(cigar);
             offset_b2 = computeRevOffset(cigar);
             ref_s_pos2 = ar.getStart2()-offset_f2;
             ref_e_pos2 = ar.getEnd2()+offset_b2;
         }
    } else {
         //get starting position and ending position according to ref position, paying attention to clipped bases
         //updated ref position including clips
         offset_f1 = computeOffset(this->cigar2_unrolled);
         offset_b1 = computeRevOffset(this->cigar2_unrolled);
         ref_s_pos1 = this->start2-offset_f1;
         ref_e_pos1 = this->end2+offset_b1;
         if(j == 1){
             std::vector<char> cigar = ar.getCigar1Unrolled();
             offset_f2 = computeOffset(cigar);
             offset_b2 = computeRevOffset(cigar);
             ref_s_pos2 = ar.getStart1()-offset_f2;
             ref_e_pos2 = ar.getEnd1()+offset_b2;
         } else {
             std::vector<char> cigar = ar.getCigar2Unrolled();
             offset_f2 = computeOffset(cigar);
             offset_b2 = computeRevOffset(cigar);
             ref_s_pos2 = ar.getStart2()-offset_f2;
             ref_e_pos2 = ar.getEnd2()+offset_b2;
         }
    }

    //position in query sequences // phred scores
    int q_pos1 = 0;
    int q_pos2 = 0;
    //position in unrolled cigar vectors
    int c_pos1 = 0;
    int c_pos2 = 0;
    //Position in the reference strings
    int ref_rel_pos1 =0;
    int ref_rel_pos2 =0;
    //4 cases of different overlaps
    //------------
    //     ------------
    if(ref_s_pos1 <= ref_s_pos2 && ref_e_pos1 <= ref_e_pos2){
        while(ref_s_pos1<ref_s_pos2){
            noOverlapMerge(algn,ref,dna,qualities,nucigar,c_pos1,q_pos1,ref_s_pos1,ref_rel_pos1,i);
        }
        ref_rel_pos2 = ref_rel_pos1;
        while(ref_s_pos1<=ref_e_pos1){
            overlapMerge(algn,ref,ar,dna,qualities,nucigar,c_pos1,c_pos2,q_pos1,q_pos2,ref_s_pos1,ref_rel_pos1,i,j);
        }
        ref_rel_pos2 = ref_rel_pos1-ref_rel_pos2;
        while(ref_s_pos1<=ref_e_pos2){
            ar.noOverlapMerge(algn,ref,dna, qualities, nucigar, c_pos2, q_pos2, ref_s_pos1,ref_rel_pos2,j);
        }
        ar.updateReference(ref,j,ref_rel_pos2);

    }//------------------------------
        //           ----------
     else if (ref_s_pos1 <= ref_s_pos2 && ref_e_pos1 >= ref_e_pos2){
        while(ref_s_pos1<ref_s_pos2){
            noOverlapMerge(algn,ref,dna,qualities,nucigar,c_pos1,q_pos1,ref_s_pos1,ref_rel_pos1,i);
        }
        while(ref_s_pos1<=ref_e_pos2){
            overlapMerge(algn,ref,ar,dna,qualities,nucigar,c_pos1,c_pos2,q_pos1,q_pos2,ref_s_pos1,ref_rel_pos1,i,j);
        }
        while(ref_s_pos1<=ref_e_pos1){
            noOverlapMerge(algn,ref,dna,qualities,nucigar,c_pos1,q_pos1,ref_s_pos1,ref_rel_pos1,i);
        }
        updateReference(ref,i,ref_rel_pos1);
        //           ----------
        //--------------------------
    } else if (ref_s_pos1 >= ref_s_pos2 && ref_e_pos1 <= ref_e_pos2){
        while(ref_s_pos2<ref_s_pos1){
            ar.noOverlapMerge(algn,ref,dna, qualities, nucigar, c_pos2, q_pos2, ref_s_pos2,ref_rel_pos2,j);
        }
        while(ref_s_pos2<=ref_e_pos1){
            overlapMerge(algn,ref,ar,dna,qualities,nucigar,c_pos1,c_pos2,q_pos1,q_pos2,ref_s_pos2,ref_rel_pos1,i,j);
        }
        ref_rel_pos2 += ref_rel_pos1-1;
        while(ref_s_pos2<=ref_e_pos2){//Why is it less that equals here?Check with paper
            ar.noOverlapMerge(algn,ref,dna, qualities, nucigar, c_pos2, q_pos2, ref_s_pos2,ref_rel_pos2,j);
        }
        ar.updateReference(ref,j,ref_rel_pos2);
        //           --------------------
        //---------------------
    } else {
        assert(ref_s_pos1 >= ref_s_pos2 && ref_e_pos1 >= ref_e_pos2);
        while(ref_s_pos2<ref_s_pos1){
            ar.noOverlapMerge(algn,ref,dna, qualities, nucigar, c_pos2, q_pos2, ref_s_pos2,ref_rel_pos2,j);
        }
        while(ref_s_pos2<=ref_e_pos2){
            overlapMerge(algn,ref,ar,dna,qualities,nucigar, c_pos1,c_pos2,q_pos1,q_pos2,ref_s_pos2,ref_rel_pos1,i,j);
        }
        while(ref_s_pos2<=ref_e_pos1){
            noOverlapMerge(algn,ref,dna,qualities,nucigar, c_pos1,q_pos1,ref_s_pos2,ref_rel_pos1,i);
        }
        updateReference(ref,i,ref_rel_pos1);

    }
    if (i == 1){
        if(j == 1){
            this->start1 = std::min(this->start1,ar.getStart1());
            this->end1=std::max(ar.getEnd1(),this->end1);
            this->single_end= true;
            this->cigar1 = createCigar(nucigar);
            this->sequence1=ShortDnaSequence(dna,qualities);
            this->phred_sum1=phred_sum(qualities);
            this->length_incl_deletions1 = this->sequence1.size();
            this->length_incl_longdeletions1 = this->sequence1.size();
            this->cigar1_unrolled.clear();
            for (char i : nucigar){
                this->cigar1_unrolled.push_back(i);
            }
            this->reference_seq1 = ref;
            this->alignSequence1 = algn;
            this->qualityList1.assign(qualities.begin(),qualities.end());

        } else {
            this->start2 = std::min(this->start1,ar.getStart2());
            this->end2 = std::max(this->end1,ar.getEnd2());
            this->cigar2 = createCigar(nucigar);
            this->sequence2=ShortDnaSequence(dna,qualities);
            this->phred_sum2=phred_sum(qualities);
            this->length_incl_deletions2 = this->sequence2.size();
            this->length_incl_longdeletions2 = this->sequence2.size();
            this->cigar2_unrolled.clear();
            for (char i : nucigar){
                this->cigar2_unrolled.push_back(i);
            }
            this->reference_seq2 = ref;
            this->alignSequence2 = algn;
            this->qualityList2.assign(qualities.begin(),qualities.end());
        }
    } else {
        this->single_end= false;
        this->cigar2 = createCigar(nucigar);
        this->sequence2=ShortDnaSequence(dna,qualities);
        this->phred_sum2=phred_sum(qualities);
        this->length_incl_deletions2 = this->sequence2.size();
        this->length_incl_longdeletions2 = this->sequence2.size();
        this->cigar2_unrolled.clear();
        for (char i : nucigar){
            this->cigar2_unrolled.push_back(i);
        }
        if(j == 1){
            this->start2 = std::min(this->start2,ar.getStart1());
            this->end2=std::max(ar.getEnd1(),this->end2);

        }else {
            this->start2 = std::min(this->start2,ar.getStart2());
            this->end2=std::max(ar.getEnd2(),this->end2);
        }
        this->reference_seq2 = ref;
        this->alignSequence2 = algn;
        this->qualityList2.assign(qualities.begin(),qualities.end());
    }
}

void AlignmentRecord::mergeAlignmentRecordsPaired(const AlignmentRecord& ar){
    std::string dna = "";
    std::string qualities = "";
    std::string nucigar = "";
    std::string ref = "";
    std::string algn = "";
    //get starting position and ending position according to ref position, paying attention to clipped bases
    int offset_f1_c1 = computeOffset(this->cigar1_unrolled);
    int offset_f1_c2 = computeOffset(this->cigar2_unrolled);
    int offset_f2_c1 = computeOffset(ar.getCigar1Unrolled());
    int offset_f2_c2 = computeOffset(ar.getCigar2Unrolled());
    int offset_b1_c1 = computeRevOffset(this->cigar1_unrolled);
    int offset_b1_c2 = computeRevOffset(this->cigar2_unrolled);
    int offset_b2_c1 = computeRevOffset(ar.getCigar1Unrolled());
    int offset_b2_c2 = computeRevOffset(ar.getCigar2Unrolled());
    //updated ref position including clips
    int ref_s_pos1_c1 = this->start1-offset_f1_c1;
    int ref_e_pos1_c1 = this->end1+offset_b1_c1;
    int ref_s_pos1_c2 = this->start2-offset_f1_c2;
    int ref_e_pos1_c2 = this->end2+offset_b1_c2;
    int ref_s_pos2_c1 = ar.getStart1()-offset_f2_c1;
    int ref_e_pos2_c1 = ar.getEnd1()+offset_b2_c1;
    int ref_s_pos2_c2 = ar.getStart2()-offset_f2_c2;
    int ref_e_pos2_c2 = ar.getEnd2()+offset_b2_c2;
    //position in query sequences // phred scores
    int q_c1_pos1 = 0;
    int q_c2_pos1 = 0;
    int q_c1_pos2 = 0;
    int q_c2_pos2 = 0;
    //position in unrolled cigar vectors
    int c_c1_pos1 = 0;
    int c_c2_pos1 = 0;
    int c_c1_pos2 = 0;
    int c_c2_pos2 = 0;
    //Position in the reference strings
    int ref_rel_pos1 =0;
    int ref_rel_pos2 =0;
    //int i;
    // --------    |  -----------    <-this
    //   --------  |     ----------

    if(this->end1 < ar.getStart2() && this->start2 > ar.getEnd1()){
        mergeAlignmentRecordsSingle(ar,1,1);
        mergeAlignmentRecordsSingle(ar,2,2);
    }//----------    ------------   <-this
    //    ---------------   -----------
    else if(ref_s_pos1_c1 <= ref_s_pos2_c1 && this->end1 >= ar.getStart1() && this->start2 <= ar.getEnd1() && this->end2 >= ar.getStart2()){
        while(ref_s_pos1_c1 < ref_s_pos2_c1){
            noOverlapMerge(algn,ref,dna,qualities,nucigar,c_c1_pos1,q_c1_pos1,ref_s_pos1_c1,ref_rel_pos1,1);
        }
        ref_rel_pos2 = ref_rel_pos1;
        while(ref_s_pos1_c1<=this->end1){
            overlapMerge(algn,ref,ar,dna,qualities,nucigar,c_c1_pos1,c_c1_pos2,q_c1_pos1,q_c1_pos2,ref_s_pos1_c1,ref_rel_pos1,1,1);
        }
        ref_rel_pos2 = ref_rel_pos1 - ref_rel_pos2 ;
        while(ref_s_pos1_c1<this->start2){
            ar.noOverlapMerge(algn,ref,dna,qualities,nucigar,c_c1_pos2,q_c1_pos2,ref_s_pos1_c1,ref_rel_pos2,1);
        }
        ref_rel_pos1 =0;
        ref_rel_pos2 =0;
        computeSOffset(this->cigar2_unrolled,c_c2_pos1,q_c2_pos1);
        while(ref_s_pos1_c1<=ar.getEnd1()){
            overlapMerge(algn,ref,ar,dna,qualities,nucigar,c_c2_pos1,c_c1_pos2,q_c2_pos1,q_c1_pos2,ref_s_pos1_c1,ref_rel_pos1,2,1);
        }
        while(ref_s_pos1_c1<ar.getStart2()){
            noOverlapMerge(algn,ref,dna,qualities,nucigar,c_c2_pos1,q_c2_pos1,ref_s_pos1_c1,ref_rel_pos1,2);
        }
        ref_rel_pos2 = ref_rel_pos1;
        computeSOffset(ar.getCigar2Unrolled(),c_c2_pos2,q_c2_pos2);
        while(ref_s_pos1_c1<=ref_e_pos2_c2 && ref_s_pos1_c1 <= ref_e_pos1_c2){
            overlapMerge(algn,ref,ar,dna,qualities,nucigar,c_c2_pos1,c_c2_pos2,q_c2_pos1,q_c2_pos2,ref_s_pos1_c1,ref_rel_pos1,2,2);
        }
        ref_rel_pos2 = ref_rel_pos1-ref_rel_pos2;
        if(ref_s_pos1_c1-1 == ref_e_pos2_c2){
            while(ref_s_pos1_c1<=ref_e_pos1_c2){
                noOverlapMerge(algn,ref,dna,qualities,nucigar,c_c2_pos1,q_c2_pos1,ref_s_pos1_c1,ref_rel_pos1,2);
            }
            updateReference(ref,2,ref_rel_pos1);
        } else if (ref_s_pos1_c1-1 == ref_e_pos1_c2){
            while(ref_s_pos1_c1<=ref_e_pos2_c2){
                ar.noOverlapMerge(algn,ref,dna,qualities,nucigar,c_c2_pos2,q_c2_pos2,ref_s_pos1_c1,ref_rel_pos2,2);
            }
            ar.updateReference(ref,2,ref_rel_pos2);
        }
        this->start1=std::min(this->start1,ar.getStart1());
        this->end1=std::max(ar.getEnd2(),this->end2);
        this->single_end = true;
        this->cigar1=createCigar(nucigar);
        this->sequence1=ShortDnaSequence(dna,qualities);
        this->phred_sum1=phred_sum(qualities);
        this->length_incl_deletions1 = this->sequence1.size();
        this->length_incl_longdeletions1 = this->sequence1.size();
        this->cigar1_unrolled.clear();
        for (char i : nucigar){
            this->cigar1_unrolled.push_back(i);
        }
        this->reference_seq1 = ref;
        this->alignSequence1 = algn;
        this->qualityList1.assign(qualities.begin(),qualities.end());
   } //-----------       ----------- <-this
    //    -----------------               -----------
    else if(ref_s_pos1_c1 <= ref_s_pos2_c1 && this->end1 >=ar.getStart1() && this->start2 <= ar.getEnd1() && this->end2 < ar.getStart2()){
        while(ref_s_pos1_c1 < ref_s_pos2_c1){
            noOverlapMerge(algn,ref,dna,qualities,nucigar,c_c1_pos1,q_c1_pos1,ref_s_pos1_c1,ref_rel_pos1,1);
        }
        ref_rel_pos2 = ref_rel_pos1;
        while(ref_s_pos1_c1<=this->end1){
            overlapMerge(algn,ref,ar,dna,qualities,nucigar,c_c1_pos1,c_c1_pos2,q_c1_pos1,q_c1_pos2,ref_s_pos1_c1,ref_rel_pos1,1,1);
        }
        ref_rel_pos2 = ref_rel_pos1-ref_rel_pos2;
        while(ref_s_pos1_c1<this->start2){
            ar.noOverlapMerge(algn,ref,dna,qualities,nucigar,c_c1_pos2,q_c1_pos2,ref_s_pos1_c1,ref_rel_pos2,1);
        }
        ref_rel_pos1 = 0;
        computeSOffset(this->cigar2_unrolled,c_c2_pos1,q_c2_pos1);
        while(ref_s_pos1_c1<=ref_e_pos2_c1 && ref_s_pos1_c1<=ref_e_pos1_c2){
            overlapMerge(algn,ref,ar,dna,qualities,nucigar,c_c2_pos1,c_c1_pos2,q_c2_pos1,q_c1_pos2,ref_s_pos1_c1,ref_rel_pos1,2,1);
        }
        ref_rel_pos2 += ref_rel_pos1;
        if(ref_s_pos1_c1-1 == ref_e_pos2_c1 ){
            while(ref_s_pos1_c1<=ref_e_pos1_c2){
                noOverlapMerge(algn,ref,dna,qualities,nucigar,c_c2_pos1,q_c2_pos1,ref_s_pos1_c1,ref_rel_pos1,2);
            }
            updateReference(ref,2,ref_rel_pos1);
        } else if(ref_s_pos1_c1-1 == ref_e_pos1_c2){
            while(ref_s_pos1_c1<=ref_e_pos2_c1){
                ar.noOverlapMerge(algn,ref,dna,qualities,nucigar,c_c1_pos2,q_c1_pos2,ref_s_pos1_c1,ref_rel_pos2,1);
            }
            ar.updateReference(ref,1,ref_rel_pos2);
        }
        this->start1=std::min(this->start1,ar.getStart1());
        this->end1=std::max(ar.getEnd1(),this->end2);
        this->single_end = false;
        this->cigar1=createCigar(nucigar);
        this->sequence1=ShortDnaSequence(dna,qualities);
        this->phred_sum1=phred_sum(qualities);
        this->length_incl_deletions1 = this->sequence1.size();
        this->length_incl_longdeletions1 = this->sequence1.size();
        this->cigar1_unrolled.clear();
        for (char i : nucigar){
            this->cigar1_unrolled.push_back(i);
        }
        this->reference_seq1 = ref;
        this->alignSequence1 = algn;
        this->qualityList1.assign(qualities.begin(),qualities.end());
        this->start2=ar.getStart2();
        this->end2=ar.getEnd2();
        this->cigar2=ar.getCigar2();
        this->sequence2=ar.getSequence2();
        this->phred_sum2=ar.getPhredSum2();
        this->length_incl_deletions2 = ar.getSequence2().size();
        this->length_incl_longdeletions1 = ar.getSequence2().size();
        this->cigar2_unrolled =ar.getCigar2Unrolled();
        this->reference_seq2 = ar.getReferenceSeq2();
        this->alignSequence2 = ar.getAlignSequence2();
        this->qualityList2 = ar.getQualityList2();

   } //----------        ------------  <- this
    //                --------    ----------
    else if(ref_s_pos2_c1 <= ref_s_pos1_c2 && this->end1 < ar.getStart1() && this->start2 <= ar.getEnd1() && this->end2 >= ar.getStart2()){
        while(ref_s_pos2_c1 < ref_s_pos1_c2){
            ar.noOverlapMerge(algn,ref,dna,qualities,nucigar,c_c1_pos2,q_c1_pos2,ref_s_pos2_c1,ref_rel_pos2,1);
        }
        while(ref_s_pos2_c1 <= ar.getEnd1()){
            overlapMerge(algn,ref,ar,dna,qualities,nucigar,c_c2_pos1,c_c1_pos2,q_c2_pos1,q_c1_pos2,ref_s_pos2_c1,ref_rel_pos1,2,1);
        }
        while(ref_s_pos2_c1 < ar.getStart2()){
            noOverlapMerge(algn,ref,dna,qualities,nucigar,c_c2_pos1,q_c2_pos1,ref_s_pos2_c1,ref_rel_pos1,2);
        }
        ref_rel_pos2=ref_rel_pos1;
        computeSOffset(ar.getCigar2Unrolled(),c_c2_pos2,q_c2_pos2);
        while(ref_s_pos2_c1 <= ref_e_pos2_c2 && ref_s_pos2_c1 <= ref_e_pos1_c2){
            overlapMerge(algn,ref,ar,dna,qualities,nucigar,c_c2_pos1,c_c2_pos2,q_c2_pos1,q_c2_pos2,ref_s_pos2_c1,ref_rel_pos1,2,2);
        }
        ref_rel_pos2 = ref_rel_pos1 - ref_rel_pos2;
        if(ref_s_pos2_c1-1 == ref_e_pos2_c2){
            while(ref_s_pos2_c1<=ref_e_pos1_c2){
                noOverlapMerge(algn,ref,dna,qualities,nucigar,c_c2_pos1,q_c2_pos1,ref_s_pos2_c1,ref_rel_pos1,2);
            }
            updateReference(ref,2,ref_rel_pos1);
        } else if(ref_s_pos2_c1-1 == ref_e_pos1_c2){
            while(ref_s_pos2_c1<=ref_e_pos2_c2){
                ar.noOverlapMerge(algn,ref,dna,qualities,nucigar,c_c2_pos2,q_c2_pos2,ref_s_pos2_c1,ref_rel_pos2,2);
            }
            ar.updateReference(ref,2,ref_rel_pos2);
        }
        this->start2=std::min(this->start2,ar.getStart1());
        this->end2=std::max(ar.getEnd2(),this->end2);
        this->single_end = false;
        this->cigar2=createCigar(nucigar);
        this->sequence2=ShortDnaSequence(dna,qualities);
        this->phred_sum2=phred_sum(qualities);
        this->length_incl_deletions2 = this->sequence2.size();
        this->length_incl_longdeletions2 = this->sequence2.size();
        this->cigar2_unrolled.clear();
        for (char i : nucigar){
            this->cigar2_unrolled.push_back(i);
        }
        this->reference_seq2 = ref;
        this->alignSequence2 = algn;
        this->qualityList2.assign(qualities.begin(),qualities.end());


    } //--------      --------- <-this
      //                --  -----------
    else if(ref_s_pos1_c2<=ref_s_pos2_c1 && this->end1 < ar.getStart1() && this->start2<=ar.getEnd1() && this->end2 >= ar.getStart2() ){
        while(ref_s_pos1_c2<ref_s_pos2_c1){
            noOverlapMerge(algn,ref,dna,qualities,nucigar,c_c2_pos1,q_c2_pos1,ref_s_pos1_c2,ref_rel_pos1,2);
        }
        while(ref_s_pos1_c2<=ar.getEnd1()){
            overlapMerge(algn,ref,ar,dna,qualities,nucigar,c_c2_pos1,c_c1_pos2,q_c2_pos1,q_c1_pos2,ref_s_pos1_c2,ref_rel_pos1,2,1);
        }
        while(ref_s_pos1_c2<ar.getStart2()){
            noOverlapMerge(algn,ref,dna,qualities,nucigar,c_c2_pos1,q_c2_pos1,ref_s_pos1_c2,ref_rel_pos1,2);
        }
        ref_rel_pos2 = ref_rel_pos1;
        computeSOffset(ar.getCigar2Unrolled(),c_c2_pos2,q_c2_pos2);
        while(ref_s_pos1_c2<=ref_e_pos1_c2 && ref_s_pos1_c2<=ref_e_pos2_c2){
            overlapMerge(algn,ref,ar,dna,qualities,nucigar,c_c2_pos1,c_c2_pos2,q_c2_pos1,q_c2_pos2,ref_s_pos1_c2,ref_rel_pos1,2,2);
        }
        ref_rel_pos2=ref_rel_pos1-ref_rel_pos2;
        if(ref_s_pos1_c2-1==ref_e_pos1_c2){
            while(ref_s_pos1_c2<=ref_e_pos2_c2){
                ar.noOverlapMerge(algn,ref,dna,qualities,nucigar,c_c2_pos2,q_c2_pos2,ref_s_pos1_c2,ref_rel_pos2,2);
            }
            ar.updateReference(ref,2,ref_rel_pos2);
        } else if(ref_s_pos1_c2-1==ref_e_pos2_c2){
            while(ref_s_pos1_c2<=ref_e_pos1_c2){
                noOverlapMerge(algn,ref,dna,qualities,nucigar,c_c2_pos1,q_c2_pos1,ref_s_pos1_c2,ref_rel_pos1,2);
            }
            updateReference(ref,2,ref_rel_pos1);
        }
        this->start2=std::min(this->start2,ar.getStart1());
        this->end2=std::max(ar.getEnd2(),this->end2);
        this->single_end = false;
        this->cigar2=createCigar(nucigar);
        this->sequence2=ShortDnaSequence(dna,qualities);
        this->phred_sum2=phred_sum(qualities);
        this->length_incl_deletions2 = this->sequence2.size();
        this->length_incl_longdeletions2 = this->sequence2.size();
        this->cigar2_unrolled.clear();
        for (char i : nucigar){
            this->cigar2_unrolled.push_back(i);
        }
        this->reference_seq2 = ref;
        this->alignSequence2 = algn;
        this->qualityList2.assign(qualities.begin(),qualities.end());


    } // -------------     ----------- <-this
        //   ----  ---------------
    else if(ref_s_pos1_c1 <= ref_s_pos2_c1 &&this->start1 <=ar.getEnd1() && this->end1 >= ar.getStart2() && this->start2 <= ar.getEnd2()){\
        while(ref_s_pos1_c1<ref_s_pos2_c1){
            noOverlapMerge(algn,ref,dna,qualities,nucigar,c_c1_pos1,q_c1_pos1,ref_s_pos1_c1,ref_rel_pos1,1);
        }
        while(ref_s_pos1_c1<=ar.getEnd1()){
            overlapMerge(algn,ref,ar,dna,qualities,nucigar,c_c1_pos1,c_c1_pos2,q_c1_pos1,q_c1_pos2,ref_s_pos1_c1,ref_rel_pos1,1,1);
        }
        while(ref_s_pos1_c1<ar.getStart2()){
            noOverlapMerge(algn,ref,dna,qualities,nucigar,c_c1_pos1,q_c1_pos1,ref_s_pos1_c1,ref_rel_pos1,1);
        }
        ref_rel_pos2 = ref_rel_pos1;
        computeSOffset(ar.getCigar2Unrolled(),c_c2_pos2,q_c2_pos2);
        while(ref_s_pos1_c1<=this->end1){
            overlapMerge(algn,ref,ar,dna,qualities,nucigar,c_c1_pos1,c_c2_pos2,q_c1_pos1,q_c2_pos2,ref_s_pos1_c1,ref_rel_pos1,1,2);
        }
        ref_rel_pos2 = ref_rel_pos1-ref_rel_pos2;
        while(ref_s_pos1_c1<this->start2){
            ar.noOverlapMerge(algn,ref,dna,qualities,nucigar,c_c2_pos2,q_c2_pos2,ref_s_pos1_c1,ref_rel_pos2,2);
        }
        ref_rel_pos1 = 0;
        computeSOffset(this->cigar2_unrolled,c_c2_pos1,q_c2_pos1);
        while(ref_s_pos1_c1<=ref_e_pos2_c2 && ref_s_pos1_c1<=ref_e_pos1_c2){
            overlapMerge(algn,ref,ar,dna,qualities,nucigar,c_c2_pos1,c_c2_pos2,q_c2_pos1,q_c2_pos2,ref_s_pos1_c1,ref_rel_pos1,2,2);
        }
        ref_rel_pos2+=ref_rel_pos1;
        if(ref_s_pos1_c1-1==ref_e_pos2_c2){
            while(ref_s_pos1_c1<=ref_e_pos1_c2){
                noOverlapMerge(algn,ref,dna,qualities,nucigar,c_c2_pos1,q_c2_pos1,ref_s_pos1_c1,ref_rel_pos1,2);
            }
            updateReference(ref,2,ref_rel_pos1);
        } else if(ref_s_pos1_c1-1==ref_e_pos1_c2){
            while(ref_s_pos1_c1<=ref_e_pos2_c2){
                ar.noOverlapMerge(algn,ref,dna,qualities,nucigar,c_c2_pos2,q_c2_pos2,ref_s_pos1_c1,ref_rel_pos2,2);
            }
            ar.updateReference(ref,2,ref_rel_pos2);
        }
        this->start1=std::min(this->start1,ar.getStart1());
        this->end1=std::max(ar.getEnd2(),this->end2);
        this->single_end = true;
        this->cigar1=createCigar(nucigar);
        this->sequence1=ShortDnaSequence(dna,qualities);
        this->phred_sum1=phred_sum(qualities);
        this->length_incl_deletions1 = this->sequence1.size();
        this->length_incl_longdeletions1 = this->sequence1.size();
        this->cigar1_unrolled.clear();
        for (char i : nucigar){
            this->cigar1_unrolled.push_back(i);
        }
        this->reference_seq1 = ref;
        this->alignSequence1 = algn;
        this->qualityList1.assign(qualities.begin(),qualities.end());


    } //----------         -----------  <-this
     //  ----   -------
    else if(ref_s_pos1_c1 <= ref_s_pos2_c1 && this->start2 > ar.getEnd2() && this->end1 >= ar.getStart2() && this->start1 <= ar.getEnd1()){
        while(ref_s_pos1_c1 < ref_s_pos2_c1){
            noOverlapMerge(algn,ref,dna,qualities,nucigar,c_c1_pos1,q_c1_pos1,ref_s_pos1_c1,ref_rel_pos1,1);
        }
        while(ref_s_pos1_c1 <= ar.getEnd1()){
            overlapMerge(algn,ref,ar,dna,qualities,nucigar,c_c1_pos1,c_c1_pos2,q_c1_pos1,q_c1_pos2,ref_s_pos1_c1,ref_rel_pos1,1,1);
        }
        while(ref_s_pos1_c1<ar.getStart2()){
            noOverlapMerge(algn,ref,dna,qualities,nucigar,c_c1_pos1,q_c1_pos1,ref_s_pos1_c1,ref_rel_pos1,1);
        }
        ref_rel_pos2 = ref_rel_pos1;
        computeSOffset(ar.getCigar2Unrolled(),c_c2_pos2,q_c2_pos2);
        while(ref_s_pos1_c1<=ref_e_pos1_c1 && ref_s_pos1_c1<=ref_e_pos2_c2){
            overlapMerge(algn,ref,ar,dna,qualities,nucigar,c_c1_pos1,c_c2_pos2,q_c2_pos1,q_c2_pos2,ref_s_pos1_c1,ref_rel_pos1,1,2);
        }
        ref_rel_pos2 = ref_rel_pos1 - ref_rel_pos2;
        if(ref_s_pos1_c1-1==ref_e_pos1_c1){
            while(ref_s_pos1_c1<=ref_e_pos2_c2){
                ar.noOverlapMerge(algn,ref,dna,qualities,nucigar,c_c2_pos2,q_c2_pos2,ref_s_pos1_c1,ref_rel_pos2,2);
            }
            ar.updateReference(ref,2,ref_rel_pos2);
        } else if(ref_s_pos1_c1-1==ref_e_pos2_c2){
            while(ref_s_pos1_c1<=ref_e_pos1_c1){
                noOverlapMerge(algn,ref,dna,qualities,nucigar,c_c1_pos1,q_c1_pos1,ref_s_pos1_c1,ref_rel_pos1,1);
            }
            updateReference(ref,1,ref_rel_pos1);
        }
        this->start1=std::min(this->start1,ar.getStart1());
        this->end1=std::max(ar.getEnd2(),this->end1);
        this->single_end = false;
        this->cigar1=createCigar(nucigar);
        this->sequence1=ShortDnaSequence(dna,qualities);
        this->phred_sum1=phred_sum(qualities);
        this->length_incl_deletions1 = this->sequence1.size();
        this->length_incl_longdeletions1 = this->sequence1.size();
        this->cigar1_unrolled.clear();
        for (char i : nucigar){
            this->cigar1_unrolled.push_back(i);
        }
        this->reference_seq1 = ref;
        this->alignSequence1 = algn;
        this->qualityList1.assign(qualities.begin(),qualities.end());
        //std::copy(qualities.begin(), qualities.end(), std::back_inserter(this->qualityList1));


    } //   ------   -----------   <-this
    //----------------  ----------------
      else if(ref_s_pos2_c1 <= ref_s_pos1_c1 && this->end1 >= ar.getStart1() && this->start2 <= ar.getEnd1() && this->end2 >= ar.getStart2()){
        while(ref_s_pos2_c1 < ref_s_pos1_c1){
            ar.noOverlapMerge(algn,ref,dna,qualities,nucigar,c_c1_pos2,q_c1_pos2,ref_s_pos2_c1,ref_rel_pos2,1);
        }
        while(ref_s_pos2_c1 <=this->end1){
            overlapMerge(algn,ref,ar,dna,qualities,nucigar,c_c1_pos1,c_c1_pos2,q_c1_pos1,q_c1_pos2,ref_s_pos2_c1,ref_rel_pos1,1,1);
        }
        ref_rel_pos2+= ref_rel_pos1;
        while(ref_s_pos2_c1 < this->start2){
            ar.noOverlapMerge(algn,ref,dna,qualities,nucigar,c_c1_pos2,q_c1_pos2,ref_s_pos2_c1,ref_rel_pos2,1);
        }
        ref_rel_pos1 = 0;
        computeSOffset(this->cigar2_unrolled,c_c2_pos1,q_c2_pos1);
        while(ref_s_pos2_c1<=ar.getEnd1()){
            overlapMerge(algn,ref,ar,dna,qualities,nucigar,c_c2_pos1,c_c1_pos2,q_c2_pos1,q_c1_pos2,ref_s_pos2_c1,ref_rel_pos1,2,1);
        }
        while(ref_s_pos2_c1<ar.getStart2()){
            noOverlapMerge(algn,ref,dna,qualities,nucigar,c_c2_pos1,q_c2_pos1,ref_s_pos2_c1,ref_rel_pos1,2);
        }
        ref_rel_pos2 = ref_rel_pos1;
        computeSOffset(ar.getCigar2Unrolled(),c_c2_pos2,q_c2_pos2);
        while(ref_s_pos2_c1<=ref_e_pos1_c2 && ref_s_pos2_c1<=ref_e_pos2_c2){
            overlapMerge(algn,ref,ar,dna,qualities,nucigar,c_c2_pos1,c_c2_pos2,q_c2_pos1,q_c2_pos2,ref_s_pos2_c1,ref_rel_pos1,2,2);
        }
        ref_rel_pos2 = ref_rel_pos1 - ref_rel_pos2;
        if(ref_s_pos2_c1-1 == ref_e_pos1_c2){
            while(ref_s_pos2_c1<=ref_e_pos2_c2){
                ar.noOverlapMerge(algn,ref,dna,qualities,nucigar,c_c2_pos2,q_c2_pos2,ref_s_pos2_c1,ref_rel_pos2,2);
            }
            ar.updateReference(ref,2,ref_rel_pos2);
        } else if(ref_s_pos2_c1-1 == ref_e_pos2_c2){
            while(ref_s_pos2_c1<=ref_e_pos1_c2){
                noOverlapMerge(algn,ref,dna,qualities,nucigar,c_c2_pos1,q_c2_pos1,ref_s_pos2_c1,ref_rel_pos1,2);
            }
            updateReference(ref,2,ref_rel_pos1);
        }
        this->start1=std::min(this->start1,ar.getStart1());
        this->end1=std::max(ar.getEnd2(),this->end2);
        this->single_end = true;
        this->cigar1=createCigar(nucigar);
        this->sequence1=ShortDnaSequence(dna,qualities);
        this->phred_sum1=phred_sum(qualities);
        this->length_incl_deletions1 = this->sequence1.size();
        this->length_incl_longdeletions1 = this->sequence1.size();
        this->cigar1_unrolled.clear();
        for (char i : nucigar){
            this->cigar1_unrolled.push_back(i);
        }
        this->reference_seq1 = ref;
        this->alignSequence1 = algn;
        this->qualityList1.assign(qualities.begin(),qualities.end());


    } //      ----    --------- <-this
      //------------------          --------------
      else if(ref_s_pos2_c1 <= ref_s_pos1_c1 && this->end1>=ar.getStart1() && this->start2 <=ar.getEnd1() && this->end2 < ar.getStart2()){
        while(ref_s_pos2_c1<ref_s_pos1_c1){
            ar.noOverlapMerge(algn,ref,dna,qualities,nucigar,c_c1_pos2,q_c1_pos2,ref_s_pos2_c1,ref_rel_pos2,1);
        }
        while(ref_s_pos2_c1<=this->end1){
            overlapMerge(algn,ref,ar,dna,qualities,nucigar,c_c1_pos1,c_c1_pos2,q_c1_pos1,q_c1_pos2,ref_s_pos2_c1,ref_rel_pos1,1,1);
        }
        ref_rel_pos2 += ref_rel_pos1;
        while(ref_s_pos2_c1<this->start2){
            ar.noOverlapMerge(algn,ref,dna,qualities,nucigar,c_c1_pos2,q_c1_pos2,ref_s_pos2_c1,ref_rel_pos2,1);
        }
        ref_rel_pos1 = 0;
        computeSOffset(this->cigar2_unrolled,c_c2_pos1,q_c2_pos1);
        while(ref_s_pos2_c1<=ref_e_pos2_c1 && ref_s_pos2_c1<=ref_e_pos1_c2){
            overlapMerge(algn,ref,ar,dna,qualities,nucigar,c_c2_pos1,c_c1_pos2,q_c2_pos1,q_c1_pos2,ref_s_pos2_c1,ref_rel_pos1,2,1);
        }
        ref_rel_pos2 += ref_rel_pos1;
        if(ref_s_pos2_c1-1==ref_e_pos2_c1){
            while(ref_s_pos2_c1<=ref_e_pos1_c2){
                noOverlapMerge(algn,ref,dna,qualities,nucigar,c_c2_pos1,q_c2_pos1,ref_s_pos2_c1,ref_rel_pos1,2);
            }
            updateReference(ref,2,ref_rel_pos1);
        } else if(ref_s_pos2_c1-1 == ref_e_pos1_c2){
            while(ref_s_pos2_c1<=ref_e_pos2_c1){
                ar.noOverlapMerge(algn,ref,dna,qualities,nucigar,c_c1_pos2,q_c1_pos2,ref_s_pos2_c1,ref_rel_pos2,1);
            }
            ar.updateReference(ref,1,ref_rel_pos2);
        }
        this->start1=std::min(this->start1,ar.getStart1());
        this->end1=std::max(ar.getEnd1(),this->end2);
        this->single_end = false;
        this->cigar1=createCigar(nucigar);
        this->sequence1=ShortDnaSequence(dna,qualities);
        this->phred_sum1=phred_sum(qualities);
        this->length_incl_deletions1 = this->sequence1.size();
        this->length_incl_longdeletions1 = this->sequence1.size();
        this->cigar1_unrolled.clear();
        for (char i : nucigar){
            this->cigar1_unrolled.push_back(i);
        }
        this->reference_seq1 = ref;
        this->alignSequence1 = algn;
        this->qualityList1.assign(qualities.begin(),qualities.end());
        this->start2=ar.getStart2();
        this->end2=ar.getEnd2();
        this->cigar2=ar.getCigar2();
        this->sequence2=ar.getSequence2();
        this->phred_sum2=ar.getPhredSum2();
        this->length_incl_deletions2 = ar.getSequence2().size();
        this->length_incl_longdeletions1 = ar.getSequence2().size();
        this->cigar2_unrolled =ar.getCigar2Unrolled();
        this->reference_seq2 = ar.getReferenceSeq2();
        this->alignSequence2 = ar.getAlignSequence2();
        this->qualityList2 = ar.getQualityList2();



    } //                -------    ------- <-this
      //-------------       ------------
      else if(ref_s_pos1_c1 <= ref_s_pos2_c2 && this->start1 > ar.getEnd1() && this->end1>=ar.getStart2() && this->start2 <= ar.getEnd2()){
        while(ref_s_pos1_c1<ref_s_pos2_c2){
            noOverlapMerge(algn,ref,dna,qualities,nucigar,c_c1_pos1,q_c1_pos1,ref_s_pos1_c1,ref_rel_pos1,1);
        }
        ref_rel_pos2 = ref_rel_pos1;
        while(ref_s_pos1_c1<=this->end1){
            overlapMerge(algn,ref,ar,dna,qualities,nucigar,c_c1_pos1,c_c2_pos2,q_c1_pos1,q_c2_pos2,ref_s_pos1_c1,ref_rel_pos1,1,2);
        }
        ref_rel_pos2 = ref_rel_pos1 - ref_rel_pos2-1;
        while(ref_s_pos1_c1<this->start2){
            ar.noOverlapMerge(algn,ref,dna,qualities,nucigar,c_c2_pos2,q_c2_pos2,ref_s_pos1_c1,ref_rel_pos2,2);
        }
        ref_rel_pos1 = 0;
        computeSOffset(this->cigar2_unrolled,c_c2_pos1,q_c2_pos1);
        while(ref_s_pos1_c1<= ref_e_pos1_c2 && ref_s_pos1_c1<= ref_e_pos2_c2){
            overlapMerge(algn,ref,ar,dna,qualities,nucigar,c_c2_pos1,c_c2_pos2,q_c2_pos1,q_c2_pos2,ref_s_pos1_c1,ref_rel_pos1,2,2);
        }
        ref_rel_pos2 += ref_rel_pos1;
        if(ref_s_pos1_c1-1 == ref_e_pos1_c2){
            while(ref_s_pos1_c1<=ref_e_pos2_c2){
                ar.noOverlapMerge(algn,ref,dna,qualities,nucigar,c_c2_pos2,q_c2_pos2,ref_s_pos1_c1,ref_rel_pos2,2);
            }
            ar.updateReference(ref,2,ref_rel_pos2);
        } else if(ref_s_pos1_c1-1 == ref_e_pos2_c2){
            while(ref_s_pos1_c1<=ref_e_pos1_c2){
                noOverlapMerge(algn,ref,dna,qualities,nucigar,c_c2_pos1,q_c2_pos1,ref_s_pos1_c1,ref_rel_pos1,2);
            }
            updateReference(ref,2,ref_rel_pos1);
        }
        this->start1=ar.getStart1();
        this->end1=ar.getEnd1();
        this->cigar1=ar.getCigar1();
        this->sequence1=ar.getSequence1();
        this->phred_sum1=ar.getPhredSum1();
        this->length_incl_deletions1 = ar.getSequence1().size();
        this->length_incl_longdeletions1 = ar.getSequence1().size();
        this->cigar1_unrolled =ar.getCigar1Unrolled();
        this->reference_seq1 = ar.getReferenceSeq1();
        this->alignSequence1 = ar.getAlignSequence1();
        this->qualityList1 = ar.getQualityList1();
        this->start2=std::min(this->start1,ar.getStart2());
        this->end2=std::max(ar.getEnd2(),this->end2);
        this->single_end = false;
        this->cigar2=createCigar(nucigar);
        this->sequence2=ShortDnaSequence(dna,qualities);
        this->phred_sum2=phred_sum(qualities);
        this->length_incl_deletions2 = this->sequence2.size();
        this->length_incl_longdeletions2 = this->sequence2.size();
        this->cigar2_unrolled.clear();
        for (char i : nucigar){
            this->cigar2_unrolled.push_back(i);
        }
        this->reference_seq2 = ref;
        this->alignSequence2 = algn;
        this->qualityList2.assign(qualities.begin(),qualities.end());

    } //                   ---   -------- <-this
      //-------------     -----------
     else if(ref_s_pos2_c2<=ref_s_pos1_c1 && ar.getEnd1()<this->start1 && this->end1>=ar.getStart2() && this->start2 <= ar.getEnd2()){
        while(ref_s_pos2_c2<ref_s_pos1_c1){
            ar.noOverlapMerge(algn,ref,dna,qualities,nucigar,c_c2_pos2,q_c2_pos2,ref_s_pos2_c2,ref_rel_pos2,2);
        }
        while(ref_s_pos2_c2<=this->end1){
            overlapMerge(algn,ref,ar,dna,qualities,nucigar,c_c1_pos1,c_c2_pos2,q_c1_pos1,q_c2_pos2,ref_s_pos2_c2,ref_rel_pos1,1,2);
        }
        ref_rel_pos2+=ref_rel_pos1;
        while(ref_s_pos2_c2<this->start2){
            ar.noOverlapMerge(algn,ref,dna,qualities,nucigar,c_c2_pos2,q_c2_pos2,ref_s_pos2_c2,ref_rel_pos2,2);
        }
        ref_rel_pos1 = 0;
        computeSOffset(this->cigar2_unrolled,c_c2_pos1,q_c2_pos1);
        while(ref_s_pos2_c2<=ref_e_pos1_c2 && ref_s_pos2_c2<=ref_e_pos2_c2){
            overlapMerge(algn,ref,ar,dna,qualities,nucigar,c_c2_pos1,c_c2_pos2,q_c2_pos1,q_c2_pos2,ref_s_pos2_c2,ref_rel_pos1,2,2);
        }
        ref_rel_pos2 += ref_rel_pos1;
        if(ref_s_pos2_c2-1==ref_e_pos1_c2){
            while(ref_s_pos2_c2<=ref_e_pos2_c2){
                ar.noOverlapMerge(algn,ref,dna,qualities,nucigar,c_c2_pos2,q_c2_pos2,ref_s_pos2_c2,ref_rel_pos2,2);
            }
            ar.updateReference(ref,2,ref_rel_pos2);
        } else if(ref_s_pos2_c2-1==ref_e_pos2_c2){
            while(ref_s_pos2_c2<=ref_e_pos1_c2){
                noOverlapMerge(algn,ref,dna,qualities,nucigar,c_c2_pos1,q_c2_pos1,ref_s_pos2_c2,ref_rel_pos1,2);
            }
            updateReference(ref,2,ref_rel_pos1);
        }
        this->start1=ar.getStart1();
        this->end1=ar.getEnd1();
        this->cigar1=ar.getCigar1();
        this->sequence1=ar.getSequence1();
        this->phred_sum1=ar.getPhredSum1();
        this->length_incl_deletions1 = this->sequence1.size();
        this->length_incl_longdeletions1 = this->sequence1.size();
        this->cigar1_unrolled = ar.getCigar1Unrolled();
        this->reference_seq1 = ar.getReferenceSeq1();
        this->alignSequence1 = ar.getAlignSequence1();
        this->qualityList1 = ar.getQualityList1();
        this->start2=std::min(this->start1,ar.getStart2());
        this->end2=std::max(ar.getEnd2(),this->end2);
        this->single_end = false;
        this->cigar2=createCigar(nucigar);
        this->sequence2=ShortDnaSequence(dna,qualities);
        this->phred_sum2=phred_sum(qualities);
        this->length_incl_deletions2 = this->sequence2.size();
        this->length_incl_longdeletions2 = this->sequence2.size();
        this->cigar2_unrolled.clear();
        for (char i : nucigar){
            this->cigar2_unrolled.push_back(i);
        }
        this->reference_seq2 = ref;
        this->alignSequence2 = algn;
        this->qualityList2.assign(qualities.begin(),qualities.end());
    }
       //  --------    --------------  <-this
       //-----    ----------------
      else if(ref_s_pos2_c1 <= ref_s_pos1_c1 && this->start1 <= ar.getEnd1() && this->end1 >= ar.getStart2() && this->start2 <= ar.getEnd2()){
        while(ref_s_pos2_c1<ref_s_pos1_c1){
            ar.noOverlapMerge(algn,ref,dna,qualities,nucigar,c_c1_pos2,q_c1_pos2,ref_s_pos2_c1,ref_rel_pos2,1);
        }
        while(ref_s_pos2_c1<=ar.getEnd1()){
            overlapMerge(algn,ref,ar,dna,qualities,nucigar,c_c1_pos1,c_c1_pos2,q_c1_pos1,q_c1_pos2,ref_s_pos2_c1,ref_rel_pos1,1,1);
        }
        while(ref_s_pos2_c1<ar.getStart2()){
            noOverlapMerge(algn,ref,dna,qualities,nucigar,c_c1_pos1,q_c1_pos1,ref_s_pos2_c1,ref_rel_pos1,1);
        }
        ref_rel_pos2 = ref_rel_pos1;
        computeSOffset(ar.getCigar2Unrolled(),c_c2_pos2,q_c2_pos2);
        while(ref_s_pos2_c1<=this->end1){
            overlapMerge(algn,ref,ar,dna,qualities,nucigar,c_c1_pos1,c_c2_pos2,q_c1_pos1,q_c2_pos2,ref_s_pos2_c1,ref_rel_pos1,1,2);
        }
        ref_rel_pos2 = ref_rel_pos1-ref_rel_pos2;
        while(ref_s_pos2_c1<this->start2){
            ar.noOverlapMerge(algn,ref,dna,qualities,nucigar,c_c2_pos2,q_c2_pos2,ref_s_pos2_c1,ref_rel_pos2,2);
        }
        ref_rel_pos1 = 0;
        computeSOffset(this->cigar2_unrolled,c_c2_pos1,q_c2_pos1);
        while(ref_s_pos2_c1<=ref_e_pos2_c2 && ref_s_pos2_c1<=ref_e_pos1_c2){
            overlapMerge(algn,ref,ar,dna,qualities,nucigar,c_c2_pos1,c_c2_pos2,q_c2_pos1,q_c2_pos2,ref_s_pos2_c1,ref_rel_pos1,2,2);
        }
        ref_rel_pos2 += ref_rel_pos1;
        if(ref_s_pos2_c1-1==ref_e_pos2_c2){
            while(ref_s_pos2_c1<=ref_e_pos1_c2){
                   noOverlapMerge(algn,ref,dna,qualities,nucigar,c_c2_pos1,q_c2_pos1,ref_s_pos2_c1,ref_rel_pos1,2);
            }
            updateReference(ref,2,ref_rel_pos1);
        } else if(ref_s_pos2_c1-1==ref_e_pos1_c2){
            while(ref_s_pos2_c1<=ref_e_pos2_c2){
                ar.noOverlapMerge(algn,ref,dna,qualities,nucigar,c_c2_pos2,q_c2_pos2,ref_s_pos2_c1,ref_rel_pos2,2);
            }
            ar.updateReference(ref,2,ref_rel_pos2);
        }
        this->start1=std::min(this->start1,ar.getStart1());
        this->end1=std::max(ar.getEnd2(),this->end2);
        this->single_end = true;
        this->cigar1=createCigar(nucigar);
        this->sequence1=ShortDnaSequence(dna,qualities);
        this->phred_sum1=phred_sum(qualities);
        this->length_incl_deletions1 = this->sequence1.size();
        this->length_incl_longdeletions1 = this->sequence1.size();
        this->cigar1_unrolled.clear();
        for (char i : nucigar){
            this->cigar1_unrolled.push_back(i);
        }
        this->reference_seq1 = ref;
        this->alignSequence1 = algn;
        this->qualityList1.assign(qualities.begin(),qualities.end());
    } //   --------------                 -------------- <-this
      //----------     -------------
      else if(ref_s_pos2_c1<=ref_s_pos1_c1 && this->start2 > ar.getEnd2() && this->start1 <= ar.getEnd1() && this->end1 >= ar.getStart2()){
        while(ref_s_pos2_c1<ref_s_pos1_c1){
            ar.noOverlapMerge(algn,ref,dna,qualities,nucigar,c_c1_pos2,q_c1_pos2,ref_s_pos2_c1,ref_rel_pos2,1);
        }
        while(ref_s_pos2_c1<=ar.getEnd1()){
            overlapMerge(algn,ref,ar,dna,qualities,nucigar,c_c1_pos1,c_c1_pos2,q_c1_pos1,q_c1_pos2,ref_s_pos2_c1,ref_rel_pos1,1,1);
        }
        while(ref_s_pos2_c1<ar.getStart2()){
            noOverlapMerge(algn,ref,dna,qualities,nucigar,c_c1_pos1,q_c1_pos1,ref_s_pos2_c1,ref_rel_pos1,1);
        }
        ref_rel_pos2 = ref_rel_pos1;
        computeSOffset(ar.getCigar2Unrolled(),c_c2_pos2,q_c2_pos2);
        while(ref_s_pos2_c1<=ref_e_pos1_c1 && ref_s_pos2_c1<=ref_e_pos2_c2){
            overlapMerge(algn,ref,ar,dna,qualities,nucigar,c_c1_pos1,c_c2_pos2,q_c1_pos1,q_c2_pos2,ref_s_pos2_c1,ref_rel_pos1,1,2);
        }
        ref_rel_pos2 = ref_rel_pos1 - ref_rel_pos2-1;
        if(ref_s_pos2_c1-1==ref_e_pos1_c1){
            while(ref_s_pos2_c1<=ref_e_pos2_c2){
                ar.noOverlapMerge(algn,ref,dna,qualities,nucigar,c_c2_pos2,q_c2_pos2,ref_s_pos2_c1,ref_rel_pos2,2);
            }
            ar.updateReference(ref,2,ref_rel_pos2);
        } else if(ref_s_pos2_c1-1==ref_e_pos2_c2){
            while(ref_s_pos2_c1<=ref_e_pos1_c1){
                noOverlapMerge(algn,ref,dna,qualities,nucigar,c_c1_pos1,q_c1_pos1,ref_s_pos2_c1,ref_rel_pos1,1);
            }
            updateReference(ref,1,ref_rel_pos1);
        }
        this->start1=std::min(this->start1,ar.getStart1());
        this->end1=std::max(ar.getEnd2(),this->end1);
        this->single_end = false;
        this->cigar1=createCigar(nucigar);
        this->sequence1=ShortDnaSequence(dna,qualities);
        this->phred_sum1=phred_sum(qualities);
        this->length_incl_deletions1 = this->sequence1.size();
        this->length_incl_longdeletions1 = this->sequence1.size();
        this->cigar1_unrolled.clear();
        for (char i : nucigar){
            this->cigar1_unrolled.push_back(i);
        }
        this->reference_seq1 = ref;
        this->alignSequence1 = algn;
        this->qualityList1.assign(qualities.begin(),qualities.end());
    }
}

void AlignmentRecord::mergeAlignmentRecordsMixed(const AlignmentRecord& ar){
    if(ar.isSingleEnd()){
        std::string dna, qualities, nucigar = "";
        std::string ref = "";
        std::string algn = "";
        int offset_s_f, offset_s_b, offset_p_f1, offset_p_f2, offset_p_b1, offset_p_b2 = 0;
        offset_s_f = computeOffset(ar.getCigar1Unrolled());
        offset_s_b = computeRevOffset(ar.getCigar1Unrolled());
        offset_p_f1 = computeOffset(this->cigar1_unrolled);
        offset_p_f2 = computeOffset(this->cigar2_unrolled);
        offset_p_b1 = computeRevOffset(this->cigar1_unrolled);
        offset_p_b2 = computeRevOffset(this->cigar2_unrolled);
        int ref_s_pos1 = ar.getStart1()-offset_s_f;
        int ref_e_pos1 = ar.getEnd1()+offset_s_b;
        int ref_p_s_pos1 = this->start1-offset_p_f1;
        int ref_p_e_pos1 = this->end1+offset_p_b1;
        int ref_p_s_pos2 = this->start2-offset_p_f2;
        int ref_p_e_pos2 = this->end2+offset_p_b2;
        //position in query sequences // phred scores
        int q_pos1 = 0;
        int q_p_pos1 = 0;
        int q_p_pos2 = 0;
        //position in unrolled cigar vectors
        int c_pos1 = 0;
        int c_p_pos1 = 0;
        int c_p_pos2 = 0;
        //Position in the reference strings
        int ref_rel_pos1 =0;
        int ref_rel_pos2 =0;
        //int i;
        // ---------     -------- ->this (second read not changed)
        //----------
        if(ar.getEnd1() < this->start2){
            mergeAlignmentRecordsSingle(ar,1,1);
            this->single_end = false;
            assert(this->end1 < this->start2);
        } // -------     -------- ->this (second read not changed)
        //             ----------
        else if (ar.getStart1() > this->end1){
            mergeAlignmentRecordsSingle(ar,2,1);
            assert(this->end1 < this->start2);
        } //----------          -----------   ->this OR  -------       -----------
        //-------------------------------              -----------------------------
        else if(ref_s_pos1 <= ref_p_s_pos1){
            while(ref_s_pos1<ref_p_s_pos1){
                ar.noOverlapMerge(algn,ref,dna, qualities, nucigar, c_pos1, q_pos1,ref_s_pos1,ref_rel_pos2,1);
            }
            while(ref_s_pos1<=this->end1){
                overlapMerge(algn,ref,ar,dna,qualities,nucigar,c_p_pos1,c_pos1,q_p_pos1,q_pos1,ref_s_pos1,ref_rel_pos1,1,1);
            }
            ref_rel_pos2 += ref_rel_pos1;
            while(ref_s_pos1<this->start2){
                ar.noOverlapMerge(algn,ref,dna,qualities,nucigar,c_pos1,q_pos1,ref_s_pos1,ref_rel_pos2,1);
            }
            ref_rel_pos1 = 0;
            computeSOffset(this->cigar2_unrolled,c_p_pos2,q_p_pos2);
            while(ref_s_pos1<=ref_p_e_pos2 && ref_s_pos1 <= ref_e_pos1){
                overlapMerge(algn,ref,ar,dna,qualities,nucigar,c_p_pos2,c_pos1,q_p_pos2,q_pos1,ref_s_pos1,ref_rel_pos1,2,1);
            }
            ref_rel_pos2 += ref_rel_pos1;
            if(ref_s_pos1-1 == ref_p_e_pos2){
                while(ref_s_pos1<=ref_e_pos1){
                    ar.noOverlapMerge(algn,ref,dna,qualities,nucigar,c_pos1,q_pos1,ref_s_pos1,ref_rel_pos2,1);
                }
                ar.updateReference(ref,1,ref_rel_pos2);
            } else if (ref_s_pos1-1 == ref_e_pos1){
                while(ref_s_pos1<=ref_p_e_pos2){
                    noOverlapMerge(algn,ref,dna,qualities,nucigar,c_p_pos2,q_p_pos2,ref_s_pos1,ref_rel_pos1,2);
                }
                updateReference(ref,2,ref_rel_pos1);
            }
            this->start1=std::min(this->start1,ar.getStart1());
            this->end1=std::max(ar.getEnd1(),this->end2);
            this->single_end = true;
            this->cigar1=createCigar(nucigar);
            this->sequence1=ShortDnaSequence(dna,qualities);
            this->phred_sum1=phred_sum(qualities);
            this->length_incl_deletions1 = this->sequence1.size();
            this->length_incl_longdeletions1 = this->sequence1.size();
            this->cigar1_unrolled.clear();
            for (char i : nucigar){
                this->cigar1_unrolled.push_back(i);
            }
            this->reference_seq1 = ref;
            this->alignSequence1 = algn;
            this->qualityList1.assign(qualities.begin(),qualities.end());
        } //----------          ------------ ->this OR ----------       -----------
        //     -------------------------------            ----------------------
        else if(ref_s_pos1 >= ref_p_s_pos1){
            while(ref_p_s_pos1 < ref_s_pos1){
                noOverlapMerge(algn,ref,dna,qualities,nucigar, c_p_pos1,q_p_pos1,ref_p_s_pos1,ref_rel_pos1,1);
            }
            ref_rel_pos2 = ref_rel_pos1;
            while(ref_p_s_pos1<=this->end1){
                overlapMerge(algn,ref,ar,dna,qualities,nucigar,c_p_pos1,c_pos1,q_p_pos1,q_pos1,ref_p_s_pos1,ref_rel_pos1,1,1);
            }
            ref_rel_pos2 = ref_rel_pos1 - ref_rel_pos2 ;
            while(ref_p_s_pos1<this->start2){
                ar.noOverlapMerge(algn,ref,dna,qualities,nucigar,c_pos1,q_pos1,ref_p_s_pos1,ref_rel_pos2,1);
            }
            ref_rel_pos1 = 0;
            computeSOffset(this->cigar2_unrolled,c_p_pos2,q_p_pos2);
            while(ref_p_s_pos1<=ref_p_e_pos2 && ref_p_s_pos1 <= ref_e_pos1){
                overlapMerge(algn,ref,ar,dna,qualities,nucigar,c_p_pos2,c_pos1,q_p_pos2,q_pos1,ref_p_s_pos1,ref_rel_pos1,2,1);
            }
            ref_rel_pos2 += ref_rel_pos1-1;
            if(ref_p_s_pos1-1 == ref_p_e_pos2){
                while(ref_p_s_pos1<=ref_e_pos1){
                    ar.noOverlapMerge(algn,ref,dna,qualities,nucigar,c_pos1,q_pos1,ref_p_s_pos1,ref_rel_pos2,1);
                }
                ar.updateReference(ref,1,ref_rel_pos2);
            } else if (ref_p_s_pos1-1 == ref_e_pos1){
                while(ref_p_s_pos1<=ref_p_e_pos2){
                    noOverlapMerge(algn,ref,dna,qualities,nucigar,c_p_pos2,q_p_pos2,ref_p_s_pos1,ref_rel_pos1,2);
                }
                updateReference(ref,2,ref_rel_pos1);
            }
            this->start1=std::min(this->start1,ar.getStart1());
            this->end1=std::max(ar.getEnd1(),this->end2);
            this->single_end = true;
            this->cigar1=createCigar(nucigar);
            this->sequence1=ShortDnaSequence(dna,qualities);
            this->phred_sum1=phred_sum(qualities);
            this->length_incl_deletions1 = this->sequence1.size();
            this->length_incl_longdeletions1 = this->sequence1.size();
            this->cigar1_unrolled.clear();
            for (char i : nucigar){
                this->cigar1_unrolled.push_back(i);
            }
            this->reference_seq1 = ref;
            this->alignSequence1 = algn;
            this->qualityList1.assign(qualities.begin(),qualities.end());
        }
    }
    else if (ar.isPairedEnd()){
        std::string dna, qualities, nucigar = "";
        std::string ref = "";
        std::string algn = "";
        int offset_s_f, offset_s_b, offset_p_f1, offset_p_f2, offset_p_b1, offset_p_b2 = 0;
        //get starting position and ending position according to ref position, paying attention to clipped bases
        //updated ref position including clips
        offset_s_f = computeOffset(this->cigar1_unrolled);
        offset_s_b = computeRevOffset(this->cigar1_unrolled);
        offset_p_f1 = computeOffset(ar.getCigar1Unrolled());
        offset_p_f2 = computeOffset(ar.getCigar2Unrolled());
        offset_p_b1 = computeRevOffset(ar.getCigar1Unrolled());
        offset_p_b2 = computeRevOffset(ar.getCigar2Unrolled());
        int ref_s_pos1 = this->start1-offset_s_f;
        int ref_e_pos1 = this->end1+offset_s_b;
        int ref_p_s_pos1 = ar.getStart1()-offset_p_f1;
        //int ref_p_e_pos1 = ar.getEnd1()+offset_p_b1;
        //int ref_p_s_pos2 = ar.getStart2()-offset_p_f2;
        int ref_p_e_pos2 = ar.getEnd2()+offset_p_b2;
        //position in query sequences // phred scores
        int q_pos1 = 0;
        int q_p_pos1 = 0;
        int q_p_pos2 = 0;
        //position in unrolled cigar vectors
        int c_pos1 = 0;
        int c_p_pos1 = 0;
        int c_p_pos2 = 0;
        //Position in the reference strings
        int ref_rel_pos1 =0;
        int ref_rel_pos2 =0;
        // ---------  -----------          OR ---------- ------------
        // ---------               ->this                -------------
        if(this->end1 < ar.getStart2()){
            mergeAlignmentRecordsSingle(ar,1,1);
            this->start2 = ar.getStart2();
            this->end2=ar.getEnd2();
            this->single_end= false;
            this->cigar2 = ar.getCigar2();
            this->sequence2=ar.getSequence2();
            this->phred_sum2=ar.getPhredSum2();
            this->length_incl_deletions2 = ar.getLengthInclDeletions2();
            this->length_incl_longdeletions2 = ar.getLengthInclLongDeletions2();
            this->cigar2_unrolled = ar.getCigar2Unrolled();
            this->reference_seq2 = ar.getReferenceSeq2();
            this->alignSequence2 = ar.getAlignSequence2();
            this->qualityList2 = ar.getQualityList2();
        } else if (this->start1 > ar.getEnd1()){
            mergeAlignmentRecordsSingle(ar,1,2);
            this->start1= ar.getStart1();
            this->end1=ar.getEnd1();
            this->single_end= false;
            this->cigar1 = ar.getCigar1();
            this->sequence1=ar.getSequence1();
            this->phred_sum1=ar.getPhredSum1();
            this->length_incl_deletions1 = ar.getLengthInclDeletions1();
            this->length_incl_longdeletions1 = ar.getLengthInclLongDeletions1();
            this->cigar1_unrolled = ar.getCigar1Unrolled();
            this->reference_seq1 = ar.getReferenceSeq1();
            this->alignSequence1 = ar.getAlignSequence1();
            this->qualityList1 = ar.getQualityList1();
        }//----------          -----------        OR  -------       -----------
        //----------------------------    <-this    -----------------------------
        else if(ref_s_pos1 <= ref_p_s_pos1){
                while(ref_s_pos1<ref_p_s_pos1){
                    noOverlapMerge(algn,ref,dna, qualities, nucigar, c_pos1, q_pos1,ref_s_pos1,ref_rel_pos1,1);
                }
                while(ref_s_pos1<=ar.getEnd1()){
                    overlapMerge(algn,ref,ar,dna,qualities,nucigar,c_pos1,c_p_pos1,q_pos1,q_p_pos1,ref_s_pos1,ref_rel_pos1,1,1);
                }
                while(ref_s_pos1<ar.getStart2()){
                    noOverlapMerge(algn,ref,dna,qualities,nucigar,c_pos1,q_pos1,ref_s_pos1,ref_rel_pos1,1);
                }
                ref_rel_pos2 = ref_rel_pos1;
                computeSOffset(ar.getCigar2Unrolled(),c_p_pos2,q_p_pos2);

                while(ref_s_pos1<=ref_p_e_pos2 && ref_s_pos1 <=ref_e_pos1){
                    overlapMerge(algn,ref,ar,dna,qualities,nucigar,c_pos1,c_p_pos2,q_pos1,q_p_pos2,ref_s_pos1,ref_rel_pos1,1,2);
                }
                ref_rel_pos2 = ref_rel_pos1 -ref_rel_pos2;
                if(ref_s_pos1-1 == ref_p_e_pos2){
                    while(ref_s_pos1<=ref_e_pos1){
                        noOverlapMerge(algn,ref,dna,qualities,nucigar,c_pos1,q_pos1,ref_s_pos1,ref_rel_pos1,1);
                    }
                    updateReference(ref,1,ref_rel_pos1);
                } else if (ref_s_pos1-1 == ref_e_pos1){
                    while(ref_s_pos1<=ref_p_e_pos2){
                        ar.noOverlapMerge(algn,ref,dna,qualities,nucigar,c_p_pos2,q_p_pos2,ref_s_pos1,ref_rel_pos2,2);
                    }
                    ar.updateReference(ref,2,ref_rel_pos2);
                }
                this->start1=std::min(this->start1,ar.getStart1());
                this->end1=std::max(ar.getEnd2(),this->end1);
                this->single_end = true;
                this->cigar1=createCigar(nucigar);
                this->sequence1=ShortDnaSequence(dna,qualities);
                this->phred_sum1=phred_sum(qualities);
                this->length_incl_deletions1 = this->sequence1.size();
                this->length_incl_longdeletions1 = this->sequence1.size();
                this->cigar1_unrolled.clear();
                for (char i : nucigar){
                    this->cigar1_unrolled.push_back(i);
                }
                this->reference_seq1 = ref;
                this->alignSequence1 = algn;
                this->qualityList1.assign(qualities.begin(), qualities.end());

        }//----------          ------------        OR ----------       -----------
        //     -------------------------------   <-this   ----------------------
        else if(ref_s_pos1 >= ref_p_s_pos1){
            while(ref_p_s_pos1 < ref_s_pos1){
                ar.noOverlapMerge(algn,ref,dna,qualities,nucigar, c_p_pos1,q_p_pos1,ref_p_s_pos1,ref_rel_pos2,1);
            }
            while(ref_p_s_pos1<=ar.getEnd1()){
                overlapMerge(algn,ref,ar,dna,qualities,nucigar,c_pos1,c_p_pos1,q_pos1,q_p_pos1,ref_p_s_pos1,ref_rel_pos1,1,1);
            }

            while(ref_p_s_pos1<ar.getStart2()){
                noOverlapMerge(algn,ref,dna,qualities,nucigar,c_pos1,q_pos1,ref_p_s_pos1,ref_rel_pos1,1);
            }
            ref_rel_pos2 = ref_rel_pos1;
            computeSOffset(ar.getCigar2Unrolled(),c_p_pos2,q_p_pos2);
            while(ref_p_s_pos1<=ref_p_e_pos2 && ref_p_s_pos1 <= ref_e_pos1){
                overlapMerge(algn,ref,ar,dna,qualities,nucigar,c_pos1,c_p_pos2,q_pos1,q_p_pos2,ref_p_s_pos1,ref_rel_pos1,1,2);
            }
            ref_rel_pos2 = ref_rel_pos1 -1;
            if(ref_p_s_pos1-1 == ref_p_e_pos2){
                while(ref_p_s_pos1<=ref_e_pos1){
                    noOverlapMerge(algn,ref,dna,qualities,nucigar,c_pos1,q_pos1,ref_p_s_pos1,ref_rel_pos1,1);
                }
                updateReference(ref,1,ref_rel_pos1);
            } else if (ref_p_s_pos1-1 == ref_e_pos1){
                while(ref_p_s_pos1<=ref_p_e_pos2){
                    ar.noOverlapMerge(algn,ref,dna,qualities,nucigar,c_p_pos2,q_p_pos2,ref_p_s_pos1,ref_rel_pos2,2);
                }
                ar.updateReference(ref,2,ref_rel_pos2);
            }
            this->start1=std::min(this->start1,ar.getStart1());
            this->end1=std::max(ar.getEnd2(),this->end1);
            this->single_end = true;
            this->cigar1=createCigar(nucigar);
            this->sequence1=ShortDnaSequence(dna,qualities);
            this->phred_sum1=phred_sum(qualities);
            this->length_incl_deletions1 = this->sequence1.size();
            this->length_incl_longdeletions1 = this->sequence1.size();
            this->cigar1_unrolled.clear();
            for (char i : nucigar){
                this->cigar1_unrolled.push_back(i);
            }
            this->reference_seq1 = ref;
            this->alignSequence1 = algn;
            this->qualityList1.assign(qualities.begin(), qualities.end());
        }
    }
}

//helper functions for merging DNA Sequences to create combined Alignment Record
void AlignmentRecord::noOverlapMerge(std::string& algn,std::string& ref,std::string& dna, std::string& qualities, std::string& nucigar, int& c_pos, int& q_pos, int& ref_pos, int& ref_rel_pos, int i) const{
    char c;
    const ShortDnaSequence* s = 0;
    std::string r = "";
    if (i == 1){
        c = this->cigar1_unrolled[c_pos];
        s = &this->sequence1;
        r = this->getReferenceSeq1();
    } else {
        c = this->cigar2_unrolled[c_pos];
        s = &this->sequence2;
        r = this->getReferenceSeq2();
    }
    if (c == 'H'){
        ref_rel_pos++;
        ref_pos++;
        c_pos++;
    } else if (c == 'I') {
        dna += (*s)[q_pos];
        algn += (*s)[q_pos];
        qualities += s->qualityChar(q_pos);
        ref+= 'I';
        nucigar += 'I';
        q_pos++;
        c_pos++;
    } else if (c == 'D') {
        algn += 'D';
        ref+= r.at(ref_rel_pos);
        ref_rel_pos++;
        nucigar += 'D';
        ref_pos++;
        c_pos++;
    } else if (c == 'S'){
        ref_rel_pos++;
        ref_pos++;
        q_pos++;
        c_pos++;
    } else if (c == 'M'){
        dna += (*s)[q_pos];
        algn += (*s)[q_pos];
        qualities += s->qualityChar(q_pos);
        ref+= r.at(ref_rel_pos);

        nucigar += c;
        ref_rel_pos++;
        ref_pos++;
        q_pos++;
        c_pos++;
    } else {
        assert(false);

    }
}



//helper function for mergeAlignmentRecords function
void AlignmentRecord::overlapMerge(std::string& algn,std::string& ref,const AlignmentRecord& ar, std::string& dna, std::string& qualities, std::string& nucigar, int& c_pos1, int& c_pos2, int& q_pos1, int& q_pos2, int& ref_pos, int& ref_rel_pos, int i, int j) const{
    char c1, c2;
    const ShortDnaSequence* s1,* s2 = 0;
    std::string r = "";
    if (i == 1){
        c1 = this->cigar1_unrolled[c_pos1];
        s1 = &this->sequence1;
        r = this->getReferenceSeq1();
        if(j == 1){
            c2 = ar.getCigar1Unrolled()[c_pos2];
            s2 = &ar.getSequence1();
        } else {
            c2 = ar.getCigar2Unrolled()[c_pos2];
            s2 = &ar.getSequence2();
        }
    } else {
        c1 = this->cigar2_unrolled[c_pos1];

        s1 = &this->sequence2;
        r = this->getReferenceSeq2();
        if(j == 1){
            c2 = ar.getCigar1Unrolled()[c_pos2];
            s2 = &ar.getSequence1();
        } else {
            c2 = ar.getCigar2Unrolled()[c_pos2];
            s2 = &ar.getSequence2();
        }
    }
    if((c1 == 'M' && c2 == 'M') || (c1 == 'S' && c2 == 'S') || (c1 == 'I' && c2 == 'I')){

        if (c1 != 'S'){
            std::pair<char,char> resPair = computeEntry((*s1)[q_pos1],s1->qualityChar(q_pos1),(*s2)[q_pos2],s2->qualityChar(q_pos2));
            dna += resPair.first;
            algn += resPair.first;
            qualities += resPair.second;
            nucigar += c1;
            if(c1=='I'){
                ref+= 'I';
            }else{
                ref+= r.at(ref_rel_pos);
            }
        }
        if (c1 != 'I'){
            ref_rel_pos++;
            ref_pos++;
        }
        q_pos1++;
        q_pos2++;
        c_pos1++;
        c_pos2++;
    } else if ((c1 == 'D' && c2 == 'D') || (c1 == 'H' && c2 == 'H') || (c1 == 'D' && c2 == 'H') || (c1 == 'H' && c2 == 'D') || (c1 == 'D' && c2 == 'S') || (c1 == 'S' && c2 == 'D')){

        c_pos1++;
        c_pos2++;
        ref+= r.at(ref_rel_pos);
        ref_rel_pos++;
        ref_pos++;
        if (c1 == 'D' || c2 == 'D'){
            nucigar += 'D';
            algn +='D';
        }
        if (c1 == 'S'){
            q_pos1++;
        } else if(c2 == 'S'){
            q_pos2++;
        }
    } else if ((c1 == 'M' && (c2 == 'D' || c2 == 'H' || c2 == 'S')) || ((c1 == 'D' || c1 == 'H' || c1 == 'S') && c2 == 'M') || (c1 == 'S' && c2 == 'H') || (c1 == 'H' && c2 == 'S')) {

        if (c1 == 'M'){
            nucigar += 'M';
            dna += (*s1)[q_pos1];
            algn+=(*s1)[q_pos1];
            qualities += s1->qualityChar(q_pos1);
            ref+= r.at(ref_rel_pos);
            ref_rel_pos++;
            ref_pos++;
            c_pos1++;
            c_pos2++;
            q_pos1++;
            if (c2 == 'S') q_pos2++;
        } else if (c2 == 'M'){
            nucigar += 'M';
            dna +=  (*s2)[q_pos2];
            algn+= (*s2)[q_pos2];
            qualities += s2->qualityChar(q_pos2);
            ref+= r.at(ref_rel_pos);
            ref_rel_pos++;
            ref_pos++;
            c_pos1++;
            c_pos2++;
            q_pos2++;
            if (c1 == 'S') q_pos1++;
        } else if (c1 == 'S'){
            ref_rel_pos++;
            ref_pos++;
            c_pos1++;
            c_pos2++;
            q_pos1++;
        } else {
            ref_rel_pos++;
            ref_pos++;
            c_pos1++;
            c_pos2++;
            q_pos2++;
        }
    } else if (c1 == 'I' || c2 == 'I'){

        if(c1 == 'I'){
            nucigar += 'I';
            dna += (*s1)[q_pos1];
            algn+=(*s1)[q_pos1];
            qualities += s1->qualityChar(q_pos1);
            ref+= 'I';
            c_pos1++;
            q_pos1++;
        } else {
            nucigar += 'I';
            dna +=  (*s2)[q_pos2];
            algn+=(*s2)[q_pos2];
            qualities += s2->qualityChar(q_pos2);
            ref+= 'I';
            c_pos2++;
            q_pos2++;
        }
    } else{
        assert(false);
    }
}



size_t AlignmentRecord::intersectionLength(const AlignmentRecord& ap) const {
    assert(single_end == ap.single_end);
    int left = max(getIntervalStart(), ap.getIntervalStart());
    int right = min(getIntervalEnd(), ap.getIntervalEnd()) + 1;
    return max(0, right-left);
}

size_t AlignmentRecord::internalSegmentIntersectionLength(const AlignmentRecord& ap) const {
    int left = max(getInsertStart(), ap.getInsertStart());
    int right = min(getInsertEnd(), ap.getInsertEnd()) + 1;
    return max(0, right-left);
}

int AlignmentRecord::getPhredSum1() const {
    return phred_sum1;
}

int AlignmentRecord::getPhredSum2() const {
    assert(!single_end);
    return phred_sum2;
}

double AlignmentRecord::getProbability() const {
    return probability;
}

unsigned int AlignmentRecord::getEnd1() const {
    return end1;
}

unsigned int AlignmentRecord::getEnd2() const {
    return end2;
}

std::string AlignmentRecord::getName() const {
    return name;
}

unsigned int AlignmentRecord::getStart1() const {
    return start1;
}

unsigned int AlignmentRecord::getStart2() const {
    return start2;
}

const std::vector<BamTools::CigarOp>& AlignmentRecord::getCigar1() const {
    return cigar1;
}

const std::vector<BamTools::CigarOp>& AlignmentRecord::getCigar2() const {
    assert(!single_end);
    return cigar2;
}

const ShortDnaSequence& AlignmentRecord::getSequence1() const {
    return sequence1;
}

const ShortDnaSequence& AlignmentRecord::getSequence2() const {
    return sequence2;
}

unsigned int AlignmentRecord::getIntervalStart() const {
    return start1;
}

unsigned int AlignmentRecord::getIntervalEnd() const {
    if (single_end) {
        return end1;
    } else {
        return end2;
    }
}

unsigned int AlignmentRecord::getInsertStart() const {
    assert(!single_end);
    return end1 + 1;
}

unsigned int AlignmentRecord::getInsertEnd() const {
    assert(!single_end);
    return start2 - 1;
}

unsigned int AlignmentRecord::getInsertLength() const {
    assert(!single_end);
    return start2 - (end1 + 1);
}

alignment_id_t AlignmentRecord::getID() const {
    return id;
}

void AlignmentRecord::setID(alignment_id_t id) {
    this->id = id;
}

bool AlignmentRecord::isSingleEnd() const {
    return single_end;
}

bool AlignmentRecord::isPairedEnd() const {
    return !single_end;
}

std::vector<std::string> AlignmentRecord::getReadNames() const {
    vector<string> rnames;
    for (int i : this->readNames) {
        rnames.push_back((*readNameMap)[i]);
    }
    return rnames;
}

const std::vector<char>& AlignmentRecord::getCigar1Unrolled() const {
    return this->cigar1_unrolled;
}
const std::vector<char>& AlignmentRecord::getCigar2Unrolled() const {
    return this->cigar2_unrolled;
}
int AlignmentRecord::getLengthInclDeletions1() const {
    return this->length_incl_deletions1;
}
int AlignmentRecord::getLengthInclDeletions2() const {
    return this->length_incl_deletions2;
}
int AlignmentRecord::getLengthInclLongDeletions1() const {
    return this->length_incl_longdeletions1;
}
int AlignmentRecord::getLengthInclLongDeletions2() const {
    return this->length_incl_longdeletions2;
}

//NoMe Additions
bool AlignmentRecord::isStrand1() const{
    return this->strand1;
}


int AlignmentRecord::getRefId() const{
    return this->refId;
}

std::string AlignmentRecord::getReferenceSeq1() const{
    return this->reference_seq1;
}
std::string AlignmentRecord::getReferenceSeq2() const{
    return this->reference_seq2;
}

std::string AlignmentRecord::getCode1() const{
    return this->code1;
}
std::string AlignmentRecord::getCode2() const{
    return this->code2;
}

std::string AlignmentRecord::getGcQuality1() const{
    return this->gcQuality1;
}
std::string AlignmentRecord::getGcQuality2() const{
    return this->gcQuality2;
}


std::string AlignmentRecord::getAlignSequence1() const{
    return this->alignSequence1;
}


std::string AlignmentRecord::getAlignSequence2() const{
    return this->alignSequence2;
}

vector<int> AlignmentRecord::getQualityList1() const{
    return this->qualityList1;
}

vector<int> AlignmentRecord::getQualityList2() const{
    return this->qualityList2;
}

std::vector<string> AlignmentRecord::getConsistingReads()const{
    return this->consistingReads;
}

//End Nome Additions

double setProbabilities(std::deque<AlignmentRecord*>& reads) {
    double read_usage_ct = 0.0;
    double mean = 1.0 / reads.size();

    for(auto&& r : reads) {
        read_usage_ct += r->getReadCount();
    }

    if (not reads.empty()) {
        read_usage_ct = max(read_usage_ct, (double) reads[0]->readNameMap->size());
    }
    double stdev = 0.0;

    for (auto&& r : reads) {
        r->probability = r->getReadCount() / read_usage_ct;

        stdev += (r->probability - mean)*(r->probability - mean);
    }

    return sqrt(1.0 / (reads.size() - 1) * stdev);
}

void printConsistingReads(std::ostream& outfile, std::deque<AlignmentRecord*>& reads, double uniqueSupportFilter){

    for (auto&& r : reads) {
        if(uniqueSupportFilter>-1){
            if(r->uniqueSupport<=uniqueSupportFilter){
                continue;
            }
        }
        outfile<< r->name<<" : ";
        for(string readName : r->getConsistingReads()){
            outfile<<readName<<" , ";
        }
        outfile<<endl;
    }
}

void printReads(std::ostream& outfile, std::deque<AlignmentRecord*>& reads, int doc_haplotypes) {
    auto comp = [](AlignmentRecord* al1, AlignmentRecord* al2) { return al1->probability > al2->probability; };
    std::sort(reads.begin(), reads.end(), comp);

    outfile.precision(5);
    outfile << std::fixed;


    if (doc_haplotypes == 0){
        for (auto&& r : reads) {
            unsigned int abs_number_reads = r->getReadNames().size();
            outfile << ">" <<r->name;
            if (not r->single_end) outfile << "|paired";
            outfile << "|ht_freq:" << r->probability;
            outfile << "|start1:" << r->getStart1();
            outfile << "|end1:" << r->getEnd1();
            if (not r->single_end){
                outfile << "|start2:" << r->getStart2();
                outfile << "|end2:" << r->getEnd2();
            }
            outfile << "|#reads:" << abs_number_reads;
            outfile << "|#strand:" << r->isStrand1();
            outfile << endl;

            outfile << r->sequence1;
            if (not r->single_end) {
                for(unsigned int i = r->end1+1; i < r->start2; i++) {
                    outfile << "N";
                }
                outfile << r->sequence2;
            }
            outfile << endl;
        }
    } else if(doc_haplotypes == 5) {
        std::vector<std::string> names;
        for (auto&& r : reads) {
            names = r->getReadNames();
            int haplo1counter = 0;
            int haplo2counter = 0;
            int haplo3counter = 0;
            int haplo4counter = 0;
            int haplo5counter = 0;
            outfile << ">" << r->name;

            if(r->single_end){
                outfile << "|ht_freq:" << r->probability;
                outfile << "|start1:" << r->getStart1();
                outfile << "|end1:" << r->getEnd1();
            }
            else if(!r->single_end && (r->getEnd1()+1 < r->getStart2())){
                outfile << "|paired";
                outfile << "|ht_freq:" << r->probability;
                outfile << "|start1:" << r->getStart1();
                outfile << "|end1:" << r->getEnd1();
                outfile << "|start2:" << r->getStart2();
                outfile << "|end2:" << r->getEnd2();
            }
            else if (!r->single_end && (r->getEnd1()+1 == r->getStart2())){
                outfile << "|ht_freq:" << r->probability;
                outfile << "|start1:" << r->getStart1();
                outfile << "|end1:" << r->getEnd2();
            }
            for(auto& i: names){
                 if (i.find("mutant1") != std::string::npos){
                    haplo1counter++;
                 } else if (i.find("mutant2") != std::string::npos) {
                    haplo2counter++;
                 } else if (i.find("mutant3") != std::string::npos) {
                    haplo3counter++;
                 } else if (i.find("mutant4") != std::string::npos) {
                    haplo4counter++;
                 } else if (i.find("mutant5") != std::string::npos) {
                    haplo5counter++;
                 }
            }
            outfile << "|ht1:" << haplo1counter;
            outfile << "|ht2:" << haplo2counter;
            outfile << "|ht3:" << haplo3counter;
            outfile << "|ht4:" << haplo4counter;
            outfile << "|ht5:" << haplo5counter;
            outfile << endl;

            outfile << r->sequence1;

            if (not r->single_end) {
                for(unsigned int i = r->end1+1; i < r->start2; i++) {
                    outfile << "N";
                }
                outfile << r->sequence2;
            }
            outfile << endl;
        }
    } else if(doc_haplotypes == 2){
        std::vector<std::string> names;
        for (auto&& r : reads) {
            names = r->getReadNames();
            int haplo1counter = 0;
            int haplo2counter = 0;
            outfile << ">" << r->name;

            if(r->single_end){
                outfile << "|ht_freq:" << r->probability;
                outfile << "|start1:" << r->getStart1();
                outfile << "|end1:" << r->getEnd1();
            }
            else if(!r->single_end && (r->getEnd1()+1 < r->getStart2())){
                outfile << "|paired";
                outfile << "|ht_freq:" << r->probability;
                outfile << "|start1:" << r->getStart1();
                outfile << "|end1:" << r->getEnd1();
                outfile << "|start2:" << r->getStart2();
                outfile << "|end2:" << r->getEnd2();
            }
            else if (!r->single_end && (r->getEnd1()+1 == r->getStart2())){
                outfile << "|ht_freq:" << r->probability;
                outfile << "|start1:" << r->getStart1();
                outfile << "|end1:" << r->getEnd2();
            }
            for(auto& i: names){
                 if (i.find("normal") != std::string::npos){
                    haplo1counter++;
                 } else if (i.find("mutant") != std::string::npos) {
                    haplo2counter++;
                 }
            }
            outfile << "|ht1:" << haplo1counter;
            outfile << "|ht2:" << haplo2counter;
            outfile << endl;

            outfile << r->sequence1;

            if (not r->single_end) {
                for(unsigned int i = r->end1+1; i < r->start2; i++) {
                    outfile << "N";
                }
                outfile << r->sequence2;
            }
            outfile << endl;
        }
    }
}


void printBAM(std::string filename, std::deque<AlignmentRecord*>& reads ,BamTools::SamHeader& header, BamTools::RefVector& references,double uniqueSupportFilter){
    BamTools::BamAlignment al;
    filename = filename+".bam";
    BamTools::BamWriter writer;
    //get Header and get References
    if ( !writer.Open(filename, header, references) ) {
        cerr << "Could not open output BAM file" << endl;
        throw std::runtime_error("Couldn't open output Bamfile");
    return;
    }
    /*
        std::string Name;               // read name
        int32_t     Length;             // length of query sequence
        std::string QueryBases;         // 'original' sequence (contained in BAM file)
        std::string AlignedBases;       // 'aligned' sequence (QueryBases plus deletion, padding, clipping chars)
        std::string Qualities;          // FASTQ qualities (ASCII characters, not numeric values)
        std::string TagData;            // tag data (use provided methods to query/modify)
        int32_t     RefID;              // ID number for reference sequence
        int32_t     Position;           // position (0-based) where alignment starts
        uint16_t    Bin;                // BAM (standard) index bin number for this alignment
        uint16_t    MapQuality;         // mapping quality score
        uint32_t    AlignmentFlag;      // alignment bit-flag (use provided methods to query/modify)
        std::vector<CigarOp> CigarData; // CIGAR operations for this alignment
        int32_t     MateRefID;          // ID number for reference sequence where alignment's mate was aligned
        int32_t     MatePosition;       // position (0-based) where alignment's mate starts
        int32_t     InsertSize;         // mate-pair insert size
        std::string Filename;           // name of BAM file which this alignment comes from

    BamAlignment::BamAlignment(void)
    : Length(0)
    , RefID(-1)
    , Position(-1)
    , Bin(0)
    , MapQuality(0)
    , AlignmentFlag(0)
    , MateRefID(-1)
    , MatePosition(-1)
    , InsertSize(0)
    */
    // iterate through all alignments and write them to output
    int counterFS = 0;
    int counterRS = 0;
    for (auto&& r : reads){
        if(r->isStrand1()){
            counterFS += 1;
        }else{
            counterRS +=1;
        }
    }
    int num = 0;
    for (auto&& r : reads){
        if(uniqueSupportFilter>-1){
            if(r->uniqueSupport<=uniqueSupportFilter){
                num ++;
                continue;
            }
        }

        //set members of al
        if(r->single_end){
            if(r->isStrand1()){
                al.QueryBases = r->getSequence1().toString();
            }else{
                al.QueryBases =  complement(r->getSequence1().toString());
            }
            al.Name = r->getName();
            al.Length = r->getSequence1().size();

            //no aligned bases
            al.Qualities = r->getSequence1().qualityString();
            //no tag data
            al.RefID = r->getRefId() ;// Changed for NoMe
            al.Position = r->getStart1()-1;
            //no bin, map quality. alignment flag is set extra
            al.MateRefID = r->getRefId() ;// Changed for NoMe. Considers the mates mapped to the same reference
            al.MatePosition = r->getStart1()-1;
            al.InsertSize = 0;
            al.CigarData = r->getCigar1();
            al.Filename = filename;

            //set flags
            al.SetIsDuplicate(false);
            al.SetIsFailedQC(false);
            al.SetIsFirstMate(true);
            al.SetIsMapped(true);
            al.SetIsMateMapped(false);
            al.SetIsMateReverseStrand(r->isStrand1());
            al.SetIsPaired(false);
            al.SetIsPrimaryAlignment(true);
            al.SetIsProperPair(false);
            al.SetIsReverseStrand(!r->isStrand1());
            al.SetIsSecondMate(false);
            writer.SaveAlignment(al);
        } else {
            if(r->isStrand1()){
                al.QueryBases = r->getSequence1().toString();
            }else{
                al.QueryBases =  complement(r->getSequence1().toString());
            }
            al.Name = r->getName();
            al.Length = r->getSequence1().size();
            //al.QueryBases = r->getSequence1().toString();
            //no aligned bases
            al.Qualities = r->getSequence1().qualityString();
            //no tag data
            //given the case that reads were mapped against more than one reference, al.RefID has to be modified
            al.RefID = r->getRefId() ;
            al.Position = r->getStart1()-1;
            //no bin, map quality. alignment flag is set extra
            al.MateRefID = r->getRefId();
            al.CigarData = r->getCigar1();
            al.MatePosition = r->getStart2()-1;
            al.InsertSize =r->getEnd2()-r->getStart1()+1;
            al.Filename = filename;

            //set flags
            al.SetIsDuplicate(false);
            al.SetIsFailedQC(false);
            al.SetIsFirstMate(true);
            al.SetIsMapped(true);
            al.SetIsMateMapped(true);
            al.SetIsMateReverseStrand(r->isStrand1());
            al.SetIsPaired(true);
            al.SetIsPrimaryAlignment(true);
            al.SetIsProperPair(true);
            al.SetIsReverseStrand(!r->isStrand1());
            al.SetIsSecondMate(false);
            writer.SaveAlignment(al);
            if(r->isStrand1()){
                al.QueryBases = r->getSequence2().toString();
            }else{
                al.QueryBases =  complement(r->getSequence2().toString());
            }
            al.Name = r->getName();
            al.Length = r->getSequence2().size();
            //al.QueryBases = r->getSequence2().toString();
            //no aligned bases
            al.Qualities = r->getSequence2().qualityString();
            //no tag data
            //given the case that reads were mapped against more than one reference, al.RefID has to be modified
            al.RefID = r->getRefId();
            al.Position = r->getStart2()-1;
            //no bin, map quality. alignment flag is set extra.
            al.CigarData = r->getCigar2();
            al.MatePosition = r->getStart1()-1;
            al.MateRefID = r->getRefId();
            al.InsertSize =(-1)*(r->getEnd2()-r->getStart1()+1);
            al.Filename = filename;
            //set flags
            al.SetIsDuplicate(false);
            al.SetIsFailedQC(false);
            al.SetIsFirstMate(false);
            al.SetIsMapped(true);
            al.SetIsMateMapped(true);
            al.SetIsMateReverseStrand(!r->isStrand1());
            al.SetIsPaired(true);
            al.SetIsPrimaryAlignment(true);
            al.SetIsProperPair(true);
            al.SetIsReverseStrand(r->isStrand1());
            al.SetIsSecondMate(true);
            writer.SaveAlignment(al);
        }
        //write to output
    }
    writer.Close();
}
