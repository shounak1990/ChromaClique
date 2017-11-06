/*bamFileChecker- Checks if two adjacent reads are first and second mate of a read
 *                Ensure that the number of first mates is equal to the number of second mates
 *
 *Requires the bamfile to be sorted by position -> samtools sort bamFile.bam > sortedBamFile.bam
 *
 * Arguements: 1. Mode: "f" if you want to find the max depth and "c" if you want to output another file which has been clipped to a required coverage
 *             2. Input Bam file (Sorted by position)
 *             3. Output Bam file (Only if the mode is "c")
 *             4. Max coverage required (Only if the mode is "c")
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
#include <queue>

#include "docopt/docopt.h"

#include <boost/tokenizer.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/algorithm/string.hpp>

#include <api/BamReader.h>
#include <api/BamAlignment.h>
#include <api/BamWriter.h>
#include <api/algorithms/Sort.h>
#include <forward_list>

using namespace std;
using namespace boost;


void clipCoverage(string filename , string outfile, int maxDepth) {
    BamTools::BamReader bamreader;
    if (not bamreader.Open(filename)) {
        cerr << bamreader.GetFilename() << endl;
        throw std::runtime_error("Couldn't open Bamfile");
    }
    BamTools::BamAlignment alignment;
    BamTools::BamWriter writer;
    if (not writer.Open(outfile,bamreader.GetConstSamHeader(),bamreader.GetReferenceData())) {
            throw std::runtime_error("Couldn't open out Bamfile");
    }

    int currentStart =0;
    int currentEnd =0;


    std::queue<std::pair<int,int>> positions;
    bool firstAlignment = true;

     while (bamreader.GetNextAlignment(alignment)){

        //Working with the first alignment
        if(firstAlignment){
            positions.push(std::make_pair((alignment.Position+1),alignment.GetEndPosition()));
            firstAlignment = false;
            continue;
        }

        //For all alignments other than the first alignment

        //Assigning the positions
        currentStart = (alignment.Position+1);
        currentEnd = alignment.GetEndPosition();

//        cout<<"current start1: "<<currentStart<<" end1: "<<currentEnd<<endl;

        while(!positions.empty()){
            //cout<<"start1: "<<firstMatePositions.front().first<<" end1: "<<firstMatePositions.front().second<<endl;
            int endPair1 = positions.front().second;
            if(currentStart>endPair1 ){
                //cout<<"pop"<<endl;
                positions.pop();
                continue;
            }
            break;
        }

        positions.push(std::make_pair(currentStart,currentEnd));

        if((positions.size()+1)<maxDepth){
//            maxDepth = firstMatePositions.size();
            writer.SaveAlignment(alignment);
        }
    }

    writer.Close();
}


void countCoverage(string filename) {
    BamTools::BamReader bamreader;
    if (not bamreader.Open(filename)) {
        cerr << bamreader.GetFilename() << endl;
        throw std::runtime_error("Couldn't open Bamfile");
    }
    BamTools::BamAlignment alignment;

    int currentStart =0;
    int currentEnd =0;
    int maxDepth =0;
    std::queue<std::pair<int,int>> positions;
    bool firstAlignment = true;

     while (bamreader.GetNextAlignment(alignment)){

        //Working with the first alignment
        if(firstAlignment){
            positions.push(std::make_pair((alignment.Position+1),alignment.GetEndPosition()));
            maxDepth = 1;
            firstAlignment = false;
            continue;
        }

        //For all alignments other than the first alignment

        //Assigning the positions
        currentStart = (alignment.Position+1);
        currentEnd = alignment.GetEndPosition();

        //cout<<"current start1: "<<currentStartPair1<<" end1: "<<currentEndPair1<<endl;

        while(!positions.empty()){
            //cout<<"start1: "<<firstMatePositions.front().first<<" end1: "<<firstMatePositions.front().second<<endl;
            int endPair1 = positions.front().second;
            if(currentStart>endPair1 ){
                //cout<<"pop"<<endl;
                positions.pop();
                continue;
            }
            break;
        }

        positions.push(std::make_pair(currentStart,currentEnd));

        if((positions.size())>maxDepth){
            maxDepth = positions.size();

        }
    }
     cout<<"Maximum depth: "<<maxDepth<<endl;

}



int main(int argc, char* argv[]) {


    cout << "Starting Coverage Clipper"<<endl;
    if(std::string(argv[1]) == "--c"){
        cout<<"Clipping coverage: "<<argv[2]<<endl;
        clipCoverage(argv[2],argv[3],stoi(argv[4]));
    }else if(std::string(argv[1]) == "--f"){
        cout<<"Counting max depth: "<<argv[2]<<endl;
        countCoverage(argv[2]);
    }

    cout << "Ended Coverage Clipper"<< endl;


}
