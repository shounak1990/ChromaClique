/* Author: Shounak Chakraborty
 */

#ifndef HMM_H_
#define HMM_H_

#include <iostream>
#include <fstream>
#include <vector>
#include <list>
#include <algorithm>
#include <map>



class hmm {
private:
    std::unordered_map<int, double> switchRateMap;
    std::unordered_map<int, int> distanceMap;
    std::unordered_map<int, double> qualityMap1;
    std::unordered_map<int, double> qualityMap2;
    std::unordered_map<int, char> observedChroma1;
    std::unordered_map<int, char> observedChroma2;
    std::unordered_map<int, double> CProb;
    std::unordered_map<int, double> OProb;
    int T;
    std::string code1;
    std::string qual1;
    std::string code2;
    std::string qual2;
    int minGCPositions;
public:
    hmm(std::unordered_map<int, double> switchRate,std::string code1, std::string qual1, std::string code2, std::string qual2,int minGCPositions);
    bool createMaps();
    double A(char i, char j, int t);
    double B1(char original, char observed, int t);
    double B2(char original, char observed, int t);
    double calculateProbability();
};

#endif //HMM_H_
