//
// Changed from GraphTools library
// Date: 2023-08-09
// Author: Huang Fei <huangfei@genomics.cn>
// All rights reserved.
//


#pragma once

#include <algorithm>
#include <string>
#include <vector>

namespace varViewer
{

struct opts_s {
    opts_s(){
        ref = NULL;
        bam = NULL;
        var = NULL;
        type = NULL;
        list = NULL;
        prefix = NULL;
        map = NULL;
        bed =NULL;
        gt = NULL;
        width = 1300;
        lane  = 10.0;
        cov   = 5000;
        per  = 1000;
        mqual = 20;
        bqual = 20;
        flank = 500;
        dedup = false;
    }

    char *ref; //reference sequence
    char *bam; //bam file
    char *var; //variant file
    char *type;//variant type
    char *list; //bam list file for CNV reference
    char *prefix; // output prefix
    char *bed;   // gene bed file for CNV
    char *map;   //mappability
    char *gt;  //variant genotype information
    int width;  //plot width
    float lane;  // lane height
    int  cov;  //CNV max cover
    int  per; // coverage depth window
    int mqual;  //mapping quality
    int bqual;//base quality
    int flank; // flank size
    bool dedup; // remove duplication
    int mean_cov=0;
};

using StringPair = std::pair<std::string, std::string>;
std::vector<std::string> splitStringByDelimiter(const std::string& str, char sep = ' ');
std::vector<std::string> splitStringByWhitespace(const std::string& str);
static inline std::string reverseString(std::string seq) {
    std::reverse(seq.begin(), seq.end());
    return seq;
}
std::string reverseComplement(std::string seq);
char complementBase(char base);
bool checkIfReferenceSequence(const std::string& sequence);
bool checkIfNucleotideReferenceSequence(const std::string& sequence);
}
