/*
 *Copyright 2017    Shounak Chakraborty
 *
 * This is a part of the preprocessing step for ChromaClique.
 * ChromaClique is a program that is used to reconstruct chromatypes from NOMe-seq data.
 *
 *
 *GC switchfinder generates the GC switch profiles for a given bamfile
 *
 * Input :  1. BamFile
 *          2. OutputFile path
 *          3. Reference File (fasta)
 *
 * Output :
 *          Format :  distance | switch/switch+nonSwitch | switch+nonSwitch | switch | nonSwitch
 *          1. FSSwitches.txt
 *          2. RSSwitches.txt
 */

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <unordered_map>
#include <stdexcept>
#include <limits>
#include <cassert>
#include <ctime>
#include <algorithm>
#include <deque>
#include <math.h>
#include "docopt/docopt.h"

#include <boost/tokenizer.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/algorithm/string.hpp>

#include <api/BamReader.h>
#include <api/BamAlignment.h>
#include <api/BamWriter.h>
#include <forward_list>

using namespace std;
using namespace boost;


int computeOffset(const std::vector<char>& cigar){

    int offset = 0;
    for(auto& i : cigar){

        if (i == 'S' || i == 'H'){
            offset++;
        } else break;
    }
    return offset;
}

//Returns the reverse complement of a given sequence
string complement(string sequence)  {
    //cout<<"reverseComplement"<<endl;
    string reverseCompliment = "";
    for(char base : sequence) {
        if(base=='A'){
            reverseCompliment += 'T';
        }else if(base=='T'){
            reverseCompliment += 'A';
        }else if(base=='G'){
            reverseCompliment += 'C';
        }else if(base=='C'){
            reverseCompliment += 'G';
        }else {
            reverseCompliment += base;
        }
    }

    return reverseCompliment;
}



//Creates the align sequence and reference sequences
std::tuple<string,string,vector<int>> referenceString(const std::vector<char> cigarData, int startPosition, std::string sequence,string & reference,string  qualities){

    //cout<<"constructing reference"<<endl;
    string constructedReference = "";
    string alignSeq="";
    vector<int> qualityList;
    int referenceCounter = startPosition-computeOffset(cigarData);
    int sequenceCounter=0;//this is the quality counter as well

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
            alignSeq+=sequence.at(sequenceCounter);
            qualityList.push_back(qualities.at(sequenceCounter)-33);
            constructedReference += (reference)[referenceCounter];
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
            alignSeq+=sequence.at(sequenceCounter);
            qualityList.push_back(qualities.at(sequenceCounter)-33);
            sequenceCounter++;
            referenceCounter++;
        }else if (cigarData[s] == 'N') {
            constructedReference += (reference)[referenceCounter];
            alignSeq+=sequence.at(sequenceCounter);
            qualityList.push_back(qualities.at(sequenceCounter)-33);
            sequenceCounter++;
            referenceCounter++;
        }else if (cigarData[s] == 'P') {
            constructedReference += (reference)[referenceCounter];
            alignSeq+=sequence.at(sequenceCounter);
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

//Takes in a coded read and  populates the switch, non-switch maps
void frequencyFinder(string codedRead, std::map<int,std::pair<int,int>>& frequencyMap){
    bool status;
    string previousCount;
    int switchCount=0;
    for (int i = 0; i < codedRead.size(); ++i) {
        if(codedRead.at(i)=='C'){
            if(switchCount!=0 && status==true){
                if(previousCount==""){
                    if(frequencyMap.find(0)== frequencyMap.end()){
                        frequencyMap.insert(std::pair<int,std::pair<int,int>>(0,std::pair<int,int>(1,0)));
                    }else{
                        frequencyMap.at(0).first  += 1;
                    }
                }else{
                    if(frequencyMap.find(std::stoi(previousCount))== frequencyMap.end()){
                         frequencyMap.insert(std::pair<int,std::pair<int,int>>(std::stoi(previousCount),std::pair<int,int>(1,0)));
                    }else{
                        frequencyMap.at(std::stoi(previousCount)).first  += 1;
                    }
                }
            }else if(switchCount!=0 && status==false){
                if(previousCount==""){
                    if(frequencyMap.find(0)== frequencyMap.end()){
                        frequencyMap.insert(std::pair<int,std::pair<int,int>>(0,std::pair<int,int>(0,1)));
                    }else{
                        frequencyMap.at(0).second  += 1;
                    }
                }else{
                    if(frequencyMap.find(std::stoi(previousCount))== frequencyMap.end()){
                         frequencyMap.insert(std::pair<int,std::pair<int,int>>(std::stoi(previousCount),std::pair<int,int>(0,1)));
                    }else{
                        frequencyMap.at(std::stoi(previousCount)).second  += 1;
                    }
                }
            }
            status = false;
            previousCount = "";
            switchCount += 1;
        }else if(codedRead.at(i)=='O'){
            if(switchCount!=0 && status==false){
                if(previousCount==""){
                    if(frequencyMap.find(0)== frequencyMap.end()){
                        frequencyMap.insert(std::pair<int,std::pair<int,int>>(0,std::pair<int,int>(1,0)));
                    }else{
                        frequencyMap.at(0).first  += 1;
                    }
                }else{
                    if(frequencyMap.find(std::stoi(previousCount))== frequencyMap.end()){
                         frequencyMap.insert(std::pair<int,std::pair<int,int>>(std::stoi(previousCount),std::pair<int,int>(1,0)));
                    }else{
                        frequencyMap.at(std::stoi(previousCount)).first  += 1;
                    }
                }
            }else if(switchCount!=0 && status==true){
                if(previousCount==""){
                    if(frequencyMap.find(0)== frequencyMap.end()){
                        frequencyMap.insert(std::pair<int,std::pair<int,int>>(0,std::pair<int,int>(0,1)));
                    }else{
                        frequencyMap.at(0).second  += 1;
                    }
                }else{
                    if(frequencyMap.find(std::stoi(previousCount))== frequencyMap.end()){
                         frequencyMap.insert(std::pair<int,std::pair<int,int>>(std::stoi(previousCount),std::pair<int,int>(0,1)));
                    }else{
                        frequencyMap.at(std::stoi(previousCount)).second  += 1;
                    }
                }
            }
            status = true;
            previousCount = "";
            switchCount += 1;
        }else{
            previousCount += codedRead.at(i);
        }
    }
}

//Codes the reads based on GC and GT regions for the OTandCTOT reads
//Sends the coded reads to  the methos "frequencFinder" which populates the switchRate maps based on the distances in the codes.
std::vector<string> FFSAndFRSSwitchDetermination(string reference,string sequence,vector<int> qualityList,std::map<int,std::pair<int,int>>& frequencyMap,char char1,char char2,char char3){

    string codedRead = "";
    string codedReadQualities = "";
    bool endChecker;

    int counter = 0;
    for(int i =0; i<sequence.size()-2 ; i++){
        counter ++;
        if((sequence.at(i)==char1) && (sequence.at(i+1)==char2) && (sequence.at(i+2)!=char1)){
            if(reference.at(i+1)!=char2){
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
            if(reference.at(i+1)==char3){
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
                    if(reference.at(i+1)!=char2){
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
                    if(reference.at(i+1)==char3){
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


    if((codedRead.find("O") != string::npos)||(codedRead.find("C") != string::npos)){
        frequencyFinder(codedRead,frequencyMap);
    }
    return {codedRead,codedReadQualities};

}

//Codes the reads based on GC and GT regions for the OBandCTOB reads
//Sends the coded reads to  the methos "frequencFinder" which populates the switchRate maps based on the distances in the codes.
std::vector<string> RFSAndRRSSwitchDetermination(string reference,string sequence,vector<int> qualityList,std::map<int,std::pair<int,int>>& frequencyMap,char char1,char char2,char char3){

    string codedRead = "";
    string codedReadQualities = "";

    int counter = 0;
    for(int i =0; i<sequence.size()-2 ; i++){
        counter ++;
        if((sequence.at(i)!=char1) && (sequence.at(i+1)==char2) && (sequence.at(i+2)==char1)){
            if(counter == 1){
                codedRead += "2O";
                codedReadQualities += "2";
                codedReadQualities += static_cast<char> (qualityList[i+1]);
            }else{
                codedRead += std::to_string(counter+1) + "O";
                codedReadQualities += std::to_string(counter+1);
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
                codedRead += "2C";
                codedReadQualities += "2";
                codedReadQualities += static_cast<char> (qualityList[i+1]);
            }else{
                codedRead += std::to_string(counter+1) + "C";
                codedReadQualities += std::to_string(counter+1);
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
    if((codedRead.find("O") != string::npos)||(codedRead.find("C") != string::npos)){
        frequencyFinder(codedRead,frequencyMap);
    }
    return {codedRead,codedReadQualities};

}



void scoreIndividualReads(string code,string qualCode,std::map<int,int> & readScoresFrequency,std::map<int,int> & switchScoresFrequency,std::map<int,int> & qualScoresFrequency,std::vector<pair<int,int>> & indTotalScoresFrequency,std::map<int,std::pair<int,int>> &switchMap){
    cout<<code<<endl;
    char prev='\0';
    char next='\0';
    int numCGpos = 0;
    string prevCount="";
    int swScore = 0;
    int qualScore=0;
    int switchScore = 0;
    int qualityScore = 0;
    int totalScore = 0;

    for(int i=0; i<code.size(); i++){
        if(code.at(i)=='C' || code.at(i)=='O'){
            if(numCGpos == 0){
                numCGpos +=1;
                prev=code.at(i);
                qualScore += qualCode.at(i);
                prevCount = "";
                continue;
            }
            next = code.at(i);
            if(prev==next){
                swScore+= -10 * log((switchMap.at(std::stoi(prevCount)).first/(double)(switchMap.at(std::stoi(prevCount)).first+switchMap.at(std::stoi(prevCount)).second))/2);
            }else{
                swScore+= -10 * log((switchMap.at(std::stoi(prevCount)).second/(double)(switchMap.at(std::stoi(prevCount)).first+switchMap.at(std::stoi(prevCount)).second))/2);
            }
            prevCount = "";
            prev = code.at(i);
            numCGpos += 1;
            qualScore += qualCode.at(i);
        }else{
            prevCount+=code.at(i);
        }
    }
   if(swScore !=0){
       switchScore = swScore/(numCGpos-1);
   }

   if(qualScore !=0){
       qualityScore = qualScore/numCGpos;
   }

   totalScore = switchScore + qualityScore;

   if(qualScore>0 && switchScore>0){
       indTotalScoresFrequency.push_back(std::pair<int,int>(qualityScore,switchScore));
   }


    if(totalScore != 0){
        if(readScoresFrequency.find(totalScore)== readScoresFrequency.end()){
            readScoresFrequency.insert(std::pair<int,int>(totalScore,1));
        }else{
            readScoresFrequency.at(totalScore)  += 1;
        }
    }

    if(switchScore != 0){
        if(switchScoresFrequency.find(switchScore)== switchScoresFrequency.end()){
            switchScoresFrequency.insert(std::pair<int,int>(switchScore,1));
        }else{
            switchScoresFrequency.at(switchScore)  += 1;
        }
    }

    if(qualityScore != 0){
        if(qualScoresFrequency.find(qualityScore)== qualScoresFrequency.end()){
            qualScoresFrequency.insert(std::pair<int,int>(qualityScore,1));
        }else{
            qualScoresFrequency.at(qualityScore)  += 1;
        }
    }


}

//Takes in a bamfile path and an outputfile path and a reference string.
//Checks whether an alignment is from the OT/CTOT or OB/CTOR strand ans codes the reads accordingly
// Creates two maps(distance,switchrate) and writes them out to file
void findSwitches(string filename , string outfile, string & reference) {
    BamTools::BamReader bamreader;
    if (not bamreader.Open(filename)) {
        cerr << bamreader.GetFilename() << endl;
        throw std::runtime_error("Couldn't open Bamfile");
    }
    BamTools::BamAlignment alignment;
    std::map<int,std::pair<int,int>> FSSwitchMap;// map containg switch and nonswitch values for FS reads
    std::map<int,std::pair<int,int>> RSSwitchMap;// map containg switch and nonswitch values for RS reads
    vector<string> FSCodes;
    vector<string> FSQualCodes;
    std::map<int,int> FSReadScoresFrequency;// map contating the score frequencies for FS reads
    std::map<int,int> FSSwitchScoresFrequency;// map contating the score frequencies for FS reads
    std::map<int,int> FSQualScoresFrequency;// map contating the score frequencies for FS reads
    std::vector<pair<int,int>> indFSTotalScores;// vector contating the individual score for FS reads
    vector<string> RSCodes;
    vector<string> RSQualCodes;
    std::map<int,int> RSReadScoresFrequency;// map contating the score frequencies for FS reads
    std::map<int,int> RSSwitchScoresFrequency;// map contating the score frequencies for FS reads
    std::map<int,int> RSQualScoresFrequency;// map contating the score frequencies for FS reads
    std::vector<pair<int,int>> indRSTotalScores;// map contating the individual score for RS reads

    while (bamreader.GetNextAlignment(alignment)){
        std::vector<char> cigar;
        for (const auto& it : alignment.CigarData) {
            for (unsigned int s = 0; s < it.Length; ++s) {
                cigar.push_back(it.Type);
            }
        }
        std::tuple<string,string,vector<int>> ref = referenceString(cigar,alignment.Position,alignment.QueryBases,reference,alignment.Qualities);
        string constructedReference = std::get<0>(ref);
        string alignSequence = std::get<1>(ref);
        vector<int> qualityList = std::get<2>(ref);


        bool strand1;
        if(!alignment.IsPaired()){
        strand1 = !alignment.IsReverseStrand();
        }else{
            if(alignment.IsFirstMate()){
                if(!alignment.IsReverseStrand()){
                    strand1=true;
                }else{
                    strand1=false;
                }
            }else if(alignment.IsSecondMate()){
                if(!alignment.IsReverseStrand()){
                    strand1=false;
                }else{
                    strand1=true;
                }
            }
        }

        // All the reads in a bam file are stored in the forward direction. So for OT and CTOT we dont chenge anything since they are forward reads
        // For OB and CTOB we take the compliment of both since we want both to be aligned to the reverse strand
        if(strand1){//OT nad CTOT
            vector<string> codes= FFSAndFRSSwitchDetermination(constructedReference,alignSequence,qualityList,FSSwitchMap,'G','C','T');
            FSCodes.push_back(codes[0]);
            FSQualCodes.push_back(codes[1]);
        }else if(!strand1){//OB and CTOB
            vector<string> codes= RFSAndRRSSwitchDetermination(complement(constructedReference),complement(alignSequence),qualityList,RSSwitchMap,'G','C','T');
            RSCodes.push_back(codes[0]);
            RSQualCodes.push_back(codes[1]);
        }


    }

    //Scoring forward strand reads
    for(int i = 0; i<FSCodes.size(); i++){
        scoreIndividualReads(FSCodes[i],FSQualCodes[i],FSReadScoresFrequency,FSSwitchScoresFrequency,FSQualScoresFrequency,indFSTotalScores,FSSwitchMap);
    }
    //Scoring reverse strand reads
    for(int i = 0; i<RSCodes.size(); i++){
        scoreIndividualReads(RSCodes[i],RSQualCodes[i],RSReadScoresFrequency,RSSwitchScoresFrequency,RSQualScoresFrequency,indRSTotalScores,RSSwitchMap);
    }



    ofstream FSSwitchFile;
    ofstream RSSwitchFile;

    FSSwitchFile.open (outfile+"FSSwitches.txt");
    RSSwitchFile.open (outfile+"RSSwitches.txt");

    //Outfile(switchrates) format: distance | switch/switch+nonSwitch | switch+nonSwitch | switch | nonSwitch

    std::map<int,std::pair<int,int>>::iterator it = FSSwitchMap.begin();
    for (it=FSSwitchMap.begin(); it!=FSSwitchMap.end(); ++it){
        double stat;
        stat = it->second.first/(double)(it->second.second+it->second.first);
        FSSwitchFile << it->first << "    " << stat <<"    "<<it->second.second+it->second.first<<"    "<<it->second.first<<"    "<<it->second.second<< '\n';

    }


    std::map<int,std::pair<int,int>>::iterator it3 = RSSwitchMap.begin();
    for (it3=RSSwitchMap.begin(); it3!=RSSwitchMap.end(); ++it3){
        double stat;
        stat = it3->second.first/(double)(it3->second.second+it3->second.first);
        RSSwitchFile << it3->first << "    " << stat <<"    "<<it3->second.second+it3->second.first<<"    "<<it3->second.first<<"    "<<it3->second.second<< '\n';

    }

}

static string readReferenceFile(std::string referencePath){
    //std::vector<string> reference;
    std::string reference = "";
    std::ifstream input(referencePath);
    std::string line="";
    if(!input.good()){
            //std::cerr << "Error opening "<<reference_path<< "  . Could Not Open Reference File." << std::endl;
            //return reference;
    }
    while( std::getline( input, line ).good() ){
        if( line.empty() || line[0] == '>' ){ // Identifier marker
            continue;
        }else if(!line.empty() || line[0] != '>'){
            reference+=line;
            line = "";
        }
    }
    boost::algorithm::to_upper(reference);
    return reference;
}

// Arguments: BamFile | OutputFilePath | ReferenceFile
int main(int argc, char* argv[]) {
   
    cout << "Starting "<< argv[1] << endl;
    string reference = readReferenceFile(argv[3]);
    findSwitches(argv[1],argv[2],reference);
    //scoreIndividualReads();
    cout << "Ended"<< argv[1] << endl;


}
