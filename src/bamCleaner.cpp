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
 *
 *bamFileCleaner- Checks if all the reads are
 *          proper pairs
 *          not duplicates
 *          both mates mapped to the same chromosome
 *          not failed QC reads
 *          both mates are mapped
 *
 * Arguements: Input Bam file
 *             Output Bam file
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

#include "docopt/docopt.h"

#include <boost/tokenizer.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/algorithm/string.hpp>

#include <api/BamReader.h>
#include <api/BamAlignment.h>
#include <api/BamWriter.h>

using namespace std;
using namespace boost;



void clipCoverage(string filename, string outfile) {
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
    while (bamreader.GetNextAlignment(alignment)){
        string referenceID= bamreader.GetReferenceData().at(alignment.RefID).RefName;
        string mateReferenceID= bamreader.GetReferenceData().at(alignment.MateRefID).RefName;
        if (!alignment.IsProperPair() || alignment.IsDuplicate() || !(referenceID==mateReferenceID) || alignment.IsFailedQC() || !(alignment.IsMapped()&&alignment.IsMateMapped())) {
            continue;
        }
        writer.SaveAlignment(alignment);
    }
    writer.Close();
}


int main(int argc, char* argv[]) {
    cout << "Starting "<< argv[1] << endl;
    clipCoverage(argv[1],argv[2]);
    cout << "Ended"<< argv[1] << endl;


}
