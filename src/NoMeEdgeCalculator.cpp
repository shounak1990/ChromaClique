/* Author: Shounak Chakraborty
 */

#include <math.h>
#include <boost/math/distributions/normal.hpp>
#include <stdlib.h>
#include <map>
#include <set>
#include <algorithm>
#include <vector>
#include <boost/tokenizer.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/algorithm/string.hpp>
#include <string>

#include "NoMeEdgeCalculator.h"
#include "NewEdgeCalculator.h"
#include "hmm.h"

 using namespace std;
 using namespace boost;


NoMeEdgeCalculator::NoMeEdgeCalculator(double nomeParam, std::unordered_map<int, double> switchRateFS, std::unordered_map<int, double> switchRateRS) {
    this->nomeParam = nomeParam;
    this->switchRateFS = switchRateFS;
    this->switchRateRS = switchRateRS;
 }

NoMeEdgeCalculator::~NoMeEdgeCalculator() {
}


std::pair<int , int> NoMeEdgeCalculator::mixedReadDetermination(const AlignmentRecord & ap1, const AlignmentRecord & ap2) const{

    //ap1 is always the read on top in the figures
    if(ap1.isSingleEnd()){
        if((ap1.getStart1()>ap2.getEnd1())&&(ap1.getEnd1()<ap2.getStart2())){
            //             ------
            //   ---------          ----------
            return std::make_pair(0, 0);
        }else if(ap1.getEnd1()<ap2.getStart1()){
            //  ------
            //              ---------       ----------
            return std::make_pair(0, 0);
        }else if(ap1.getStart1()>ap2.getEnd2()){
            //                              ------
            //  ---------       ----------
            return std::make_pair(0, 0);
        }else if((ap1.getStart1()<ap2.getEnd1())&&(ap1.getEnd1()<ap2.getStart2())){
            //   ----------
            //  --------       -----------
            return std::make_pair(1, 1);
        }else if((ap1.getStart1()>ap2.getEnd1())&&(ap1.getEnd1()>ap2.getStart2())){
            //               ------------
            //   -------        ------------
            return std::make_pair(1, 2);
        }else if((ap1.getStart1()<ap2.getEnd1())&&(ap1.getEnd1()>ap2.getStart2())){
            //       --------------
            //   ---------    ----------
            return std::make_pair(1, 3);
        }

    }else if(ap1.isPairedEnd()){
        if((ap2.getStart1()>ap1.getEnd1())&&(ap2.getEnd1()<ap1.getStart2())){
            //  ----------             ------------
            //             ---------
            return std::make_pair(0, 0);
        }else if(ap2.getStart1()>ap1.getEnd2()){
            //  ----------      ------------
            //                                   ---------
            return std::make_pair(0, 0);
        }else if(ap2.getEnd1()<ap1.getStart1()){
            //                     ----------      ------------
            //    ---------
            return std::make_pair(0, 0);
        }else if((ap2.getStart1()<ap1.getEnd1())&&(ap2.getEnd1()<ap1.getStart2())){
            //   ----------    -------------
            //  --------
            return std::make_pair(1, 1);
        }else if((ap2.getStart1()>ap1.getEnd1())&&(ap2.getEnd1()>ap1.getStart2())){
            // ------------       ------------
            //                  ------------
            return std::make_pair(2, 1);
        }else if((ap2.getStart1()<ap1.getEnd2())&&(ap2.getEnd1()>ap1.getStart2())){
            //  ----------      ------------
            //      -------------------
            return std::make_pair(1, 4);
        }
    }

}


int NoMeEdgeCalculator::pairedReadDetermination(const AlignmentRecord& ar1, const AlignmentRecord& ar2)const{

    // --------    |  -----------    <-ar1
    //   --------  |     ----------

    if(ar1.getEnd1() < ar2.getStart2() && ar1.getStart2() > ar2.getEnd1()){
        return 1;
    }//----------    ------------   <-ar1
    //    ---------------   -----------
    else if(ar1.getStart1() <= ar2.getStart2() && ar1.getEnd1() >= ar2.getStart1() && ar1.getStart2() <= ar2.getEnd1() && ar1.getEnd2() >= ar2.getStart2()){
        return 2;

   } //-----------       ----------- <-ar1
    //    -----------------               -----------
    else if(ar1.getStart1() <= ar2.getStart2() && ar1.getEnd1() >=ar2.getStart1() && ar1.getStart2() <= ar2.getEnd1() && ar1.getEnd2() < ar2.getStart2()){
        return 3;
   } //----------        ------------  <- ar1
    //                --------    ----------
    else if(ar2.getStart1() <= ar1.getStart2() && ar1.getEnd1() < ar2.getStart1() && ar1.getStart2() <= ar2.getEnd1() && ar1.getEnd2() >= ar2.getStart2()){
        return 4;
    } //--------      --------- <-ar1
      //                --  -----------
    else if(ar1.getStart2() <= ar2.getStart1() && ar1.getEnd1() < ar2.getStart1() && ar1.getStart2() <= ar2.getEnd1() && ar1.getEnd2() >= ar2.getStart2() ){
        return 5;
    } // -------------     ----------- <-ar1
        //   ----  ---------------
    else if(ar1.getStart1() <= ar2.getStart1() &&ar1.getStart1() <= ar2.getEnd1() && ar1.getEnd1() >= ar2.getStart2() && ar1.getStart2() <= ar2.getEnd2()){\
        return 6;
    } //----------         -----------  <-ar1
     //  ----   -------
    else if(ar1.getStart1() <= ar2.getStart1() && ar1.getStart2() > ar2.getEnd2() && ar1.getEnd1() >= ar2.getStart2() && ar1.getStart1() <= ar2.getEnd1()){
        return 7;
    } //   ------   -----------   <-ar1
    //----------------  ----------------
      else if(ar2.getStart1() <= ar1.getStart1() && ar1.getEnd1() >= ar2.getStart1() && ar1.getStart2() <= ar2.getEnd1() && ar1.getEnd2() >= ar2.getStart2()){
        return 8;
    } //      ----    --------- <-ar1
      //------------------          --------------
      else if(ar2.getStart1() <= ar1.getStart1() && ar1.getEnd1() >= ar2.getStart1() && ar1.getStart2() <= ar2.getEnd1() && ar1.getEnd2() < ar2.getStart2()){
        return 9;

    } //                -------    ------- <-ar1
      //-------------       ------------
      else if(ar1.getStart1() <= ar2.getStart2() && ar1.getStart1() > ar2.getEnd1() && ar1.getEnd1() >= ar2.getStart2() && ar1.getStart2() <= ar2.getEnd2()){
        return 10;
    } //                   ---   -------- <-ar1
      //-------------     -----------
     else if(ar2.getStart2() <= ar1.getStart1() && ar2.getEnd1()<ar1.getStart1() && ar1.getEnd1() >= ar2.getStart2() && ar1.getStart2() <= ar2.getEnd2()){
        return 11;
    }
       //  --------    ------------  <-ar1
       //-----    ----------------
      else if(ar2.getStart1() <= ar1.getStart1() && ar1.getStart1() <= ar2.getEnd1() && ar1.getEnd1() >= ar2.getStart2() && ar1.getStart2() <= ar2.getEnd2()){
        return 12;
    } //   --------------                 -------------- <-ar1
      //----------     -------------
      else if(ar2.getStart1() <= ar1.getStart1() && ar1.getStart2() > ar2.getEnd2() && ar1.getStart1() <= ar2.getEnd1() && ar1.getEnd1() >= ar2.getStart2()){
        return 13;
    }
}



bool NoMeEdgeCalculator::edgeBetween(const AlignmentRecord & ap1, const AlignmentRecord & ap2, int numGCAllowedPos, int ct) const{
    if (ap1.isStrand1()&& !ap2.isStrand1()) {
        return false;
    }else if (!ap1.isStrand1()&& ap2.isStrand1()) {
        return false;
    }
    bool isClique = false;
    if((ap1.getName().find("Clique")!=std::string::npos)||(ap2.getName().find("Clique")!=std::string::npos)){
        isClique = true;
    }


      //For Ideal Graph. Assumes that the reads have chroma1 or chroma2 written in their names
//    if((ap1.getName().find("chroma1") != std::string::npos && ap2.getName().find("chroma1") != std::string::npos) || (ap1.getName().find("chroma2") != std::string::npos && ap2.getName().find("chroma2") != std::string::npos)){
//        if(((ap1.getCode1().find("C") != std::string::npos) || (ap1.getCode1().find("O") != std::string::npos)) || ((ap1.getCode2().find("C") != std::string::npos) || (ap1.getCode2().find("O") != std::string::npos))){
//            if(((ap2.getCode1().find("C") != std::string::npos) || (ap2.getCode1().find("O") != std::string::npos)) || ((ap2.getCode2().find("C") != std::string::npos) || (ap2.getCode2().find("O") != std::string::npos))){
//                if(scoreReads(ap1,ap2, isClique,numGCAllowedPos,ct)){
//                    return true;
//                }else{
//                    return false;
//                }
//            }else{
//                return false;
//            }
//        }else{
//            return false;
//        }
//    }else{
//        return false;
//    }

    //For normal edge criterion
    if(((ap1.getCode1().find("C") != std::string::npos) || (ap1.getCode1().find("O") != std::string::npos)) || ((ap1.getCode2().find("C") != std::string::npos) || (ap1.getCode2().find("O") != std::string::npos))){
        if(((ap2.getCode1().find("C") != std::string::npos) || (ap2.getCode1().find("O") != std::string::npos)) || ((ap2.getCode2().find("C") != std::string::npos) || (ap2.getCode2().find("O") != std::string::npos))){
            cout<<"---------------------------------"<<endl;
            cout<<ap1.getName()<<endl;
            cout<<ap2.getName()<<endl;
            if(scoreReads(ap1,ap2, isClique,numGCAllowedPos,ct)){
                return true;
            }else{
                return false;
            }
        }else{
            return false;
        }
    }else{
        return false;
    }
}

void NoMeEdgeCalculator::getPartnerLengthRange(const AlignmentRecord& ap, unsigned int* min, unsigned int* max) const {
    assert(false);
}

void NoMeEdgeCalculator::setOverlapCliques(double d){
    assert(false);
}


void trimCodes(string & code, string & qual, int overlap, int & rightNonOverlapGCCount, bool hasRightOverlap){

    int i = 0;
    string prevCount="0";
    int distance = 0;
    string returnCode = "";
    string returnQual = "";
    string rightNonOverlap = "";
    while(i<code.size()){
        if(code.at(i)=='O' || code.at(i)=='C'){
            distance += stoi(prevCount);
            distance += 2;
            if(distance>=overlap){
                if(hasRightOverlap){
                    rightNonOverlap = code.substr(i,(code.size()-i));
                    int noOs= std::count(rightNonOverlap.begin(), rightNonOverlap.end(), 'O');
                    int noCs= std::count(rightNonOverlap.begin(), rightNonOverlap.end(), 'C');
                    rightNonOverlapGCCount = noOs + noCs;
                }
                break;
            }
            returnCode += prevCount;
            returnCode += code.at(i);
            returnQual += prevCount;
            returnQual += qual.at(i);
            prevCount = "0";
        }else{
            if(prevCount=="0"){
                prevCount = code.at(i);
            }else{
                prevCount += code.at(i);
            }
        }
        i++;
    }
    code = returnCode;
    qual = returnQual;
}




//Return 1 if there should be an edge: i.e all the GC positions match or the ChromatypeProb > nomeParam
//Return -1 if there should be no edge: i.e all the GC positions do not match and the ChromatypeProb < nomeParam
//Return 0 if there are no GC params
double NoMeEdgeCalculator::scoreSequences(std::string code1,std::string qual1, std::string code2,std::string qual2 , int start1, int end1, int start2, int end2 , bool isStrand1,bool isClique,int numGCAllowedPos, int ct) const{

    if((code1.find("C") != std::string::npos) || (code1.find("O") != std::string::npos)){
        if((code2.find("C") != std::string::npos) || (code2.find("O") != std::string::npos)){

        }else{
            return -1;
        }
    }else{
        return -1;
    }

    int overlap= 0;
    int s1=0;
    int s2=0;
    if(start1<start2){
        s2 = 0;
        s1 = start2-start1;
    }else if (start1>start2) {
        s1=0;
        s2=start1-start2;
    }else if (start1==start2) {
        s1=0;
        s2=0;
    }

    if(end1<end2){
        if(start1>start2){
            //          ----------------
           //     -------------------------------
            overlap=end1-start1+1;
        }else if (start1<start2) {
            //     --------------
            //        -----------------
            overlap=end1-start2+1;
        }else if (start1==start2) {
            //     -------------
            //     -----------------------
            overlap=end1-start1+1;
        }

    }else if (end1>end2) {
        if(start1>start2){
            //               ------------------
            //        -----------------
            overlap=end2-start1+1;
        }else if (start1<start2) {
            //     --------------------
            //        -----------
            overlap=end2-start2;
        }else if (start1==start2) {
            //        ------------------
            //        -----------
            overlap=end2-start1+1;
        }
    }else if (end1==end2) {
        if(start1<start2){
            //        ------------------
            //             -------------
            overlap=end1-start2+1;
        }else if (start1>start2) {
            //               -----------
            //        ------------------
            overlap=end1-start1+1;
        }else if (start1 == start2) {
            //        ------------------
            //        ------------------
            overlap=end1-start1+1 ;
        }
    }

    if(overlap<0){
        return -1;
    }


    string compareString1="";
    string compareQual1="";
    string compareString2="";
    string compareQual2="";
    int leftNonOverlapGCCount = -1;
    int rightNonOverlapGCCount = 0;

    if(s1==s2){
        leftNonOverlapGCCount = 0;
        compareString1 = code1;
        compareQual1 = qual1;
        compareString2 = code2;
        compareQual2 = qual2;
    }else if(s2==0){
        compareString2 = code2;
        compareQual2 = qual2;
        int distance = 0;
        int i = 0;
        string prevCountCode = "0";
        string prevCountQual = "0";
        while(i<code1.size() && i<qual1.size()){
            if(code1.at(i)=='O'||code1.at(i)=='C'){
                leftNonOverlapGCCount++;
                if(prevCountCode!=prevCountQual){
                    return -1;
                }
                distance = distance + std::stoi(prevCountCode);
                if(distance==s1){
                    compareString1 += code1.substr(i);
                    compareQual1 += qual1.substr(i);
                    break;
                }else if(distance>s1){
                   compareString1 += std::to_string(distance - s1);
                   compareQual1 += std::to_string(distance - s1);
                   compareString1 += code1.substr(i);
                   compareQual1 += qual1.substr(i);
                   break;
                }

                distance = distance + 2;
                prevCountCode = "0";
                prevCountQual = "0";
                i++;
            }else{
                prevCountCode = prevCountCode + code1.at(i);
                prevCountQual = prevCountQual + code1.at(i);
                i++;
            }
        }
    }else if(s1==0){
        compareString1 = code1;
        compareQual1 = qual1;
        int distance = 0;
        int i = 0;
        string prevCountCode = "0";
        string prevCountQual = "0";
        while(i<code2.size() && i<qual2.size()){
            if(code2.at(i)=='O'||code2.at(i)=='C'){
                leftNonOverlapGCCount++;
                if(prevCountCode!=prevCountQual){
                    return -1;
                }
                distance = distance + std::stoi(prevCountCode);
                if(distance==s2){
                    compareString2 += code2.substr(i);
                    compareQual2 += qual2.substr(i);
                    break;
                }else if(distance>s2){
                   compareString2 += std::to_string(distance - s2);
                   compareQual2 += std::to_string(distance - s2);
                   compareString2 += code2.substr(i);
                   compareQual2 += qual2.substr(i);
                   break;
                }

                distance = distance + 2;
                prevCountCode = "0";
                prevCountQual = "0";
                i++;
            }else{
                prevCountCode = prevCountCode + code2.at(i);
                prevCountQual = prevCountQual + code2.at(i);
                i++;
            }
        }
    }

    if(compareString1.size()==0){
        compareString1 += std::to_string(overlap);
        compareQual1 += std::to_string(overlap);
    }else if(compareString2.size()==0){
        compareString2 += std::to_string(overlap);
        compareQual2 += std::to_string(overlap);
    }
      if(compareString1.size()>compareString2.size()){
        trimCodes(compareString1,compareQual1,overlap,rightNonOverlapGCCount,true);
        trimCodes(compareString2,compareQual2,overlap,rightNonOverlapGCCount,false);
    }else if(compareString1.size()<compareString2.size()){
        trimCodes(compareString1,compareQual1,overlap,rightNonOverlapGCCount,false);
        trimCodes(compareString2,compareQual2,overlap,rightNonOverlapGCCount,true);
    }else{
        trimCodes(compareString1,compareQual1,overlap,rightNonOverlapGCCount,false);
        trimCodes(compareString2,compareQual2,overlap,rightNonOverlapGCCount,false);
    }

    int numberOfPositions1= std::count(compareString1.begin(), compareString1.end(), 'O') + std::count(compareString1.begin(), compareString1.end(), 'C');
    int numberOfPositions2= std::count(compareString2.begin(), compareString2.end(), 'O') + std::count(compareString2.begin(), compareString2.end(), 'C');
    if(numberOfPositions1!=numberOfPositions2) return -1;

    if(numberOfPositions1<numGCAllowedPos){
        return -1;
    }


    if(compareString1.size()!=compareString2.size()){
        return -1;
    }


    //Needed for checking the non overlapping sections
//    if(nonOverlap>0){
//        if(leftNonOverlapGCCount>nonOverlap || rightNonOverlapGCCount>nonOverlap){
//            return -1;
//        }
//    }

    std::unordered_map<int, double> switchRateMap;
    if(isStrand1){
        switchRateMap = this->switchRateFS;
    }else{
        switchRateMap = this->switchRateRS;
    }
    std::auto_ptr<hmm> hmmProbability(new hmm(switchRateMap,compareString1,compareQual1,compareString2,compareQual2,numGCAllowedPos));

   double score = hmmProbability->calculateProbability();
   double prob = nomeParam;

   cout<<"code1: "<<compareString1<<endl;
   cout<<"code2: "<<compareString2<<endl;
   cout<<"score: "<<score<<endl;
   cout<<"---------------------------------"<<endl;
   if(score == -1000){//meaning that there is something wrong with the codes
       return -1;
   }else{
//       //For ideal Graph
//       //specificPositionErrorMap[numberOfPositions1] = std::make_pair((specificPositionErrorMap[numberOfPositions1].first+1),(specificPositionErrorMap[numberOfPositions1].second));
//       return 1;
       //For normal edge Criterion
       if(score > prob){
           return 1;
       }else{
           return -1;
       }
   }
}


bool NoMeEdgeCalculator::scoreReads(const AlignmentRecord& ap1, const AlignmentRecord& ap2, bool isClique,int numGCAllowedPos, int ct) const{

    double score =0;

    double scorePair1 = 0;
    double scorePair2 = 0;


    bool pe1 = ap1.isPairedEnd();
    bool pe2 = ap2.isPairedEnd();

    // special cases of paired end reads for which no edge is allowed
    if(pe1 && pe2){

        unsigned int e11= ap1.getEnd1();
        unsigned int e12= ap1.getEnd2();
        unsigned int s11= ap1.getStart1();
        unsigned int s12= ap1.getStart2();
        unsigned int e21= ap2.getEnd1();
        unsigned int e22= ap2.getEnd2();
        unsigned int s21= ap2.getStart1();
        unsigned int s22= ap2.getStart2();
        //--------   ---------
        //           ---------   -------
        if(e11 < s21 && e12 < s12){
            if((ap1.getName().find("Clique") != string::npos)&&(ap2.getName().find("Clique") != string::npos)){
            }

            return false;
        } //         ---------    ---------
        //-------    ---------
        else if(e21 < s11 && e22 < s12){
            if((ap1.getName().find("Clique") != string::npos)&&(ap2.getName().find("Clique") != string::npos)){
            }
            return false;
        } // ------------     ----------
        //   -----------                  ---------
        else if(e21 < s12 && e12 < s22){
            if((ap1.getName().find("Clique") != string::npos)&&(ap2.getName().find("Clique") != string::npos)){
            }
            return false;
        }//-----------            ----------
        //-----------  ----------
        else if(e11 < s22 && e22 < s12){
            if((ap1.getName().find("Clique") != string::npos)&&(ap2.getName().find("Clique") != string::npos)){
            }
            return false;
        }//          --------- ---------
        //----------           ---------
        else if (e21 < s11 && e11 < s22){
            if((ap1.getName().find("Clique") != string::npos)&&(ap2.getName().find("Clique") != string::npos)){
            }
            return false;
        }//--------             --------
        //           --------   --------
        else if(e11 < s21 && e21 < s12){
            if((ap1.getName().find("Clique") != string::npos)&&(ap2.getName().find("Clique") != string::npos)){
            }
            return false;
        }
    }


    //Scoring for single ended cliques
     if(ap1.isSingleEnd() && ap2.isSingleEnd() ){
         score = scoreSequences(ap1.getCode1(),ap1.getGcQuality1(),ap2.getCode1(),ap2.getGcQuality1(),ap1.getStart1(),ap1.getEnd1(),ap2.getStart1(),ap2.getEnd1() , ap1.isStrand1(),isClique,numGCAllowedPos, ct);
         if(score==1)return true;
         else return false;

     }



    //For mixed single and paired end cliques
    if((ap1.isSingleEnd() && pe2)||(pe1 && ap2.isSingleEnd())){
        std::pair<int , int> result =mixedReadDetermination(ap1,ap2);
        int first = result.first;
        int second = result.second;
        if(first==0||second==0){
            return false;
        }
        if(first==1){
            if(second==1){
                scorePair1 = scoreSequences(ap1.getCode1(),ap1.getGcQuality1(),ap2.getCode1(),ap2.getGcQuality1(),ap1.getStart1(),ap1.getEnd1(),ap2.getStart1(),ap2.getEnd1() , ap1.isStrand1(),isClique,numGCAllowedPos,ct);
                scorePair2 = 0;
            }else if(second == 2){
                scorePair2 = scoreSequences(ap1.getCode1(),ap1.getGcQuality1(),ap2.getCode2(),ap2.getGcQuality2(),ap1.getStart1(),ap1.getEnd1(),ap2.getStart2(),ap2.getEnd2() , ap1.isStrand1(),isClique,numGCAllowedPos,ct);
                scorePair1 = 0;
            }else if(second == 3){
                scorePair1 = scoreSequences(ap1.getCode1(),ap1.getGcQuality1(),ap2.getCode1(),ap2.getGcQuality1(),ap1.getStart1(),ap1.getEnd1(),ap2.getStart1(),ap2.getEnd1() , ap1.isStrand1(),isClique,numGCAllowedPos,ct);
                if(scorePair1 == -1) return false;
                //For the second Mate
                scorePair2 = scoreSequences(ap1.getCode1(),ap1.getGcQuality1(),ap2.getCode2(),ap2.getGcQuality2(),ap1.getStart1(),ap1.getEnd1(),ap2.getStart2(),ap2.getEnd2() , ap1.isStrand1(),isClique,numGCAllowedPos,ct);
                if(scorePair2 == -1) return false;
            }else if(second == 4){
                scorePair1 = scoreSequences(ap1.getCode1(),ap1.getGcQuality1(),ap2.getCode1(),ap2.getGcQuality1(),ap1.getStart1(),ap1.getEnd1(),ap2.getStart1(),ap2.getEnd1() , ap1.isStrand1(),isClique,numGCAllowedPos,ct);
                if(scorePair1 == -1) return false;
                //For the second Mate
                scorePair2=scoreSequences(ap1.getCode2(),ap1.getGcQuality2(),ap2.getCode1(),ap2.getGcQuality1(),ap1.getStart2(),ap1.getEnd2(),ap2.getStart1(),ap2.getEnd1() , ap1.isStrand1(),isClique,numGCAllowedPos,ct);
                if(scorePair2 == -1) return false;
            }
        }else{
            if(second == 1){
                scorePair2 = scoreSequences(ap1.getCode2(),ap1.getGcQuality2(),ap2.getCode1(),ap2.getGcQuality1(),ap1.getStart2(),ap1.getEnd2(),ap2.getStart1(),ap2.getEnd1() , ap1.isStrand1(),isClique,numGCAllowedPos,ct);
                scorePair1 = 0;
            }
        }

        if(scorePair1+scorePair2>=2)return true;
        else return false;


    }



     double tempScore = 0;
    //Scoring sequences for paired end cliques
    if(pe1&&pe2){
        int returnVal = pairedReadDetermination(ap1,ap2);
        if(returnVal == 1 ){
          tempScore = scoreSequences(ap1.getCode1(),ap1.getGcQuality1(),ap2.getCode1(),ap2.getGcQuality1(),ap1.getStart1(),ap1.getEnd1(),ap2.getStart1(),ap2.getEnd1() , ap1.isStrand1(),isClique,numGCAllowedPos,ct);
          if (tempScore < 0) return false;
          score += tempScore;
          tempScore = scoreSequences(ap1.getCode2(),ap1.getGcQuality2(),ap2.getCode2(),ap2.getGcQuality2(),ap1.getStart2(),ap1.getEnd2(),ap2.getStart2(),ap2.getEnd2() , ap1.isStrand1(),isClique,numGCAllowedPos,ct);
          if (tempScore < 0) return false;
          score += tempScore;
        }else if(returnVal == 2 || returnVal == 8){
            tempScore = scoreSequences(ap1.getCode1(),ap1.getGcQuality1(),ap2.getCode1(),ap2.getGcQuality1(),ap1.getStart1(),ap1.getEnd1(),ap2.getStart1(),ap2.getEnd1() , ap1.isStrand1(),isClique,numGCAllowedPos,ct);
            if (tempScore < 0) return false;
            score += tempScore;
            tempScore = scoreSequences(ap1.getCode2(),ap1.getGcQuality2(),ap2.getCode1(),ap2.getGcQuality1(),ap1.getStart2(),ap1.getEnd2(),ap2.getStart1(),ap2.getEnd1() , ap1.isStrand1(),isClique,numGCAllowedPos,ct);
            if (tempScore < 0) return false;
            score += tempScore;
            tempScore=scoreSequences(ap1.getCode2(),ap1.getGcQuality2(),ap2.getCode2(),ap2.getGcQuality2(),ap1.getStart2(),ap1.getEnd2(),ap2.getStart2(),ap2.getEnd2() , ap1.isStrand1(),isClique,numGCAllowedPos,ct);
            if (tempScore < 0) return false;
            score += tempScore;
        }else if(returnVal == 3 || returnVal == 9){
            tempScore = scoreSequences(ap1.getCode1(),ap1.getGcQuality1(),ap2.getCode1(),ap2.getGcQuality1(),ap1.getStart1(),ap1.getEnd1(),ap2.getStart1(),ap2.getEnd1() , ap1.isStrand1(),isClique,numGCAllowedPos,ct);
            if (tempScore < 0) return false;
            score += tempScore;
            tempScore = scoreSequences(ap1.getCode2(),ap1.getGcQuality2(),ap2.getCode1(),ap2.getGcQuality1(),ap1.getStart2(),ap1.getEnd2(),ap2.getStart1(),ap2.getEnd1() , ap1.isStrand1(),isClique,numGCAllowedPos,ct);
            if (tempScore < 0) return false;
            score += tempScore;
        }else if(returnVal == 4 || returnVal == 5){
            tempScore = scoreSequences(ap1.getCode2(),ap1.getGcQuality2(),ap2.getCode1(),ap2.getGcQuality1(),ap1.getStart2(),ap1.getEnd2(),ap2.getStart1(),ap2.getEnd1() , ap1.isStrand1(),isClique,numGCAllowedPos,ct);
            if (tempScore < 0) return false;
            score += tempScore;
            tempScore=scoreSequences(ap1.getCode2(),ap1.getGcQuality2(),ap2.getCode2(),ap2.getGcQuality2(),ap1.getStart2(),ap1.getEnd2(),ap2.getStart2(),ap2.getEnd2() , ap1.isStrand1(),isClique,numGCAllowedPos,ct);
            if (tempScore < 0) return false;
            score += tempScore;
        }else if(returnVal == 6){
            tempScore = scoreSequences(ap1.getCode1(),ap1.getGcQuality1(),ap2.getCode1(),ap2.getGcQuality1(),ap1.getStart1(),ap1.getEnd1(),ap2.getStart1(),ap2.getEnd1() , ap1.isStrand1(),isClique,numGCAllowedPos,ct);
            if (tempScore < 0) return false;
            score += tempScore;
            tempScore = scoreSequences(ap1.getCode1(),ap1.getGcQuality1(),ap2.getCode2(),ap2.getGcQuality2(),ap1.getStart1(),ap1.getEnd1(),ap2.getStart2(),ap2.getEnd2() , ap1.isStrand1(),isClique,numGCAllowedPos,ct);
            if (tempScore < 0) return false;
            score += tempScore;
            tempScore=scoreSequences(ap1.getCode2(),ap1.getGcQuality2(),ap2.getCode2(),ap2.getGcQuality2(),ap1.getStart2(),ap1.getEnd2(),ap2.getStart2(),ap2.getEnd2() , ap1.isStrand1(),isClique,numGCAllowedPos,ct);
            if (tempScore < 0) return false;
            score += tempScore;
        }else if(returnVal == 7 || returnVal == 13){
            tempScore = scoreSequences(ap1.getCode1(),ap1.getGcQuality1(),ap2.getCode1(),ap2.getGcQuality1(),ap1.getStart1(),ap1.getEnd1(),ap2.getStart1(),ap2.getEnd1() , ap1.isStrand1(),isClique,numGCAllowedPos,ct);
            if (tempScore < 0) return false;
            score += tempScore;
            //For the second Mate
            tempScore=scoreSequences(ap1.getCode1(),ap1.getGcQuality1(),ap2.getCode2(),ap2.getGcQuality2(),ap1.getStart1(),ap1.getEnd1(),ap2.getStart2(),ap2.getEnd2() , ap1.isStrand1(),isClique,numGCAllowedPos,ct);
            if (tempScore < 0) return false;
            score += tempScore;
        }else if(returnVal == 10 || returnVal == 11){
            tempScore = scoreSequences(ap1.getCode1(),ap1.getGcQuality1(),ap2.getCode2(),ap2.getGcQuality2(),ap1.getStart1(),ap1.getEnd1(),ap2.getStart2(),ap2.getEnd2() , ap1.isStrand1(),isClique,numGCAllowedPos,ct);
            if (tempScore < 0) return false;
            score += tempScore;
            //For the second Mate
            tempScore = scoreSequences(ap1.getCode2(),ap1.getGcQuality2(),ap2.getCode2(),ap2.getGcQuality2(),ap1.getStart2(),ap1.getEnd2(),ap2.getStart2(),ap2.getEnd2() , ap1.isStrand1(),isClique,numGCAllowedPos,ct);
            if (tempScore < 0) return false;
            score += tempScore;
        }else if(returnVal == 12){
            tempScore = scoreSequences(ap1.getCode1(),ap1.getGcQuality1(),ap2.getCode1(),ap2.getGcQuality1(),ap1.getStart1(),ap1.getEnd1(),ap2.getStart1(),ap2.getEnd1() , ap1.isStrand1(),isClique,numGCAllowedPos,ct);
            if (tempScore < 0) return false;
            score += tempScore;
            tempScore = scoreSequences(ap1.getCode1(),ap1.getGcQuality1(),ap2.getCode2(),ap2.getGcQuality2(),ap1.getStart1(),ap1.getEnd1(),ap2.getStart2(),ap2.getEnd2() , ap1.isStrand1(),isClique,numGCAllowedPos,ct);
            if (tempScore < 0) return false;
            score += tempScore;
            tempScore=scoreSequences(ap1.getCode2(),ap1.getGcQuality2(),ap2.getCode2(),ap2.getGcQuality2(),ap1.getStart2(),ap1.getEnd2(),ap2.getStart2(),ap2.getEnd2() , ap1.isStrand1(),isClique,numGCAllowedPos,ct);
            if (tempScore < 0) return false;
            score += tempScore;
        }

        if(score>=2)return true;
        else return false;
    }



}

