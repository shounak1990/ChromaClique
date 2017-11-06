/* Author: Shounak Chakraborty
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
#include<hmm.h>

using namespace std;


hmm::hmm(std::unordered_map<int, double> switchRate,string code1, string qual1, string code2, string qual2,int minGCPositions){
    this->switchRateMap = switchRate;
    this->code1 = code1;
    this->code2 = code2;
    this-> qual1 = qual1;
    this-> qual2 = qual2;
    this->minGCPositions = minGCPositions;
}

bool hmm::createMaps(){
    int i = 0;
    int numGCPos = 0;
    string previousCount1 = "0";
    string previousCount2 = "0";
    //cout<<"here1"<<endl;
    while(i<this->code1.size() && i<this->qual1.size() && i<this->code2.size() && i<this->qual2.size()){
        //cout<<"here2: "<<i<<endl;
        if(this->code1.at(i)=='O'||this->code1.at(i)=='C'){
            if(this->code2.at(i)=='O'||this->code2.at(i)=='C'){
                if(previousCount1!=previousCount2){
                    return false;
                }
                if(numGCPos==0){
                    distanceMap[numGCPos+1] = (0);
                    observedChroma1[numGCPos+1] = (this->code1.at(i));
                    this->qualityMap1[numGCPos+1] = (pow(10,(-this->qual1.at(i)/(double)10)));
                    //qualityMap1.push_back(qual1.at(i));
                    observedChroma2[numGCPos+1] = (this->code2.at(i));
                    qualityMap2[numGCPos+1] = (pow(10,(-this->qual2.at(i)/(double)10)));
                    //qualityMap2.push_back(qual2.at(i));
                    previousCount1 = "0";
                    previousCount2 = "0";
                    numGCPos += 1;
                    i++;
                    continue;
                }
                observedChroma1[numGCPos+1] = (this->code1.at(i));
                qualityMap1[numGCPos+1] = (pow(10,(-this->qual1.at(i)/(double)10)));
                //qualityMap1.push_back(qual1.at(i));
                observedChroma2[numGCPos+1] = (this->code2.at(i));
                qualityMap2[numGCPos+1] = (pow(10,(-this->qual2.at(i)/(double)10)));
                //qualityMap2.push_back(qual2.at(i));
                distanceMap[numGCPos+1] = std::stoi(previousCount1);
                numGCPos += 1;
                previousCount1 = "0";
                previousCount2 = "0";
            }else{
                return false;
            }
        }else{
            if(previousCount1=="0"||previousCount2=="0"){
                previousCount1 = code1.at(i);
                previousCount2 = code2.at(i);
            }else{
                previousCount1 = previousCount1 + code1.at(i);
                previousCount2 = previousCount2 + code2.at(i);
            }
        }
        i++;
    }
    T = numGCPos;
    if(T<this->minGCPositions){
        return false;
    }
    return true;
}

double hmm::A(char i, char j, int t){
    if(i==j){
        return (1-(this->switchRateMap[this->distanceMap[t]]/2));
    }else{
        return (this->switchRateMap[this->distanceMap[t]]/2);
    }
}

double hmm::B1(char original, char observed, int t){
    if(original==observed){
        return (1-this->qualityMap1[t]);
    }else{
        return (this->qualityMap1[t]);
    }
}

double hmm::B2(char original, char observed, int t){
    if(original==observed){
        return (1-this->qualityMap2[t]);
    }else{
        return (this->qualityMap2[t]);
    }
}

double hmm::calculateProbability(){
    if(!this->createMaps()){
        return -1000;
    }
    this->CProb[1] = 0.5 * B1('C',this->observedChroma1[1],1) * B2('C',this->observedChroma2[1],1);
    this->OProb[1] = 0.5 * B1('O',this->observedChroma1[1],1) * B2('O',this->observedChroma2[1],1);
//    cout<<"Initial: "<<endl;
//    cout<<"Cprob[1]: "<<this->CProb[1]<<endl;
//    cout<<"Oprob[1]: "<<this->CProb[1]<<endl;
//    cout<<"----------------------------"<<endl;
//    double scaleFactor = 1/((double)(this->CProb[1]+this->OProb[1]));
    this->CProb[1] = this->CProb[1];
    this->OProb[1] = this->OProb[1];
//    double scaleFactorLog = log10(scaleFactor);

    for(int t = 1; t<T; t++){
        this->CProb[t+1] = B1('C',this->observedChroma1[t+1],t+1) * B2('C',this->observedChroma2[t+1],t+1) * ( A('C','C',t+1)*this->CProb[t] + A('O','C',t+1)*this->OProb[t] );
        this->OProb[t+1] = B1('O',this->observedChroma1[t+1],t+1) * B2('O',this->observedChroma2[t+1],t+1) * ( A('C','O',t+1)*this->CProb[t] + A('O','O',t+1)*this->OProb[t] );
//        cout<<"Iteration: "<<t+1<<endl;
//        cout<<"Cprob["<<t+1<<"]: "<<this->CProb[t+1]<<endl;
//        cout<<"Oprob["<<t+1<<"]: "<<this->OProb[t+1]<<endl;
//        cout<<"----------------------------"<<endl;
//        scaleFactor = 1/((double)(this->CProb[t+1]+this->OProb[t+1]));
//        this->CProb[t+1] = this->CProb[t+1]*(double)scaleFactor;
//        this->OProb[t+1] =this->OProb[t+1]*(double)scaleFactor;

//        scaleFactorLog += log10(scaleFactor);
    }

    double prob = this->CProb[T]+this->OProb[T];
    return prob/*pow10(-1*scaleFactorLog)*/;

}
