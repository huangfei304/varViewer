//
// Changed from GraphTools library
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//      http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.
//

#include "utils.hh"

#include <algorithm>
#include <sstream>
#include <stdexcept>
#include <unordered_map>

using std::string;
using std::unordered_map;
using std::vector;

namespace varViewer
{

vector<string> splitStringByDelimiter(const std::string& str, char sep) {
    vector<string> words;
    std::stringstream sstream(str);
    string word;
    while (std::getline(sstream, word, sep))
    {
        words.push_back(word);
    }

    return words;
}

vector<string> splitStringByWhitespace(const std::string& str) {
    vector<string> words;
    std::stringstream sstream(str);
    string word;
    while (sstream >> word)
    {
        words.push_back(word);
    }
    return words;
}

char complementBase(char base) {
    switch (base) {
    case 'A':
        return 'T';
    case 'a':
        return 't';
    case 'C':
        return 'G';
    case 'c':
        return 'g';
    case 'G':
        return 'C';
    case 'g':
        return 'c';
    case 'T':
        return 'A';
    case 't':
        return 'a';
    case 'R':
        return 'Y';
    case 'Y':
        return 'R';
    case 'K':
        return 'M';
    case 'M':
        return 'K';
    case 'S':
        return 'S';
    case 'W':
        return 'W';
    case 'B':
        return 'V';
    case 'D':
        return 'H';
    case 'H':
        return 'D';
    case 'V':
        return 'B';
    case ' ':
        return ' ';
    default:
        return 'N';
    }
}

string reverseComplement(string seq) {
    std::transform(seq.begin(), seq.end(), seq.begin(), complementBase);
    std::reverse(seq.begin(), seq.end());
    return seq;
}

const unordered_map<char, string> kSymbolExpansion
    = { { 'A', "A" },   { 'C', "C" },   { 'T', "T" },    { 'G', "G" },  { 'R', "AG" },  { 'Y', "CT" },
        { 'K', "GT" },  { 'M', "AC" },  { 'S', "CG" },   { 'W', "AT" }, { 'B', "CGT" }, { 'D', "AGT" },
        { 'H', "ACT" }, { 'V', "ACG" }, { 'N', "ACGT" }, { 'X', "X" } };

static bool checkIfNucleotideReferenceSymbol(char symbol) {
    return (symbol == 'A') || (symbol == 'C') || (symbol == 'T') || (symbol == 'G');
}

static bool hasExpandableSymbols(const std::string& s) {
    static const struct ExpandableSyms
    {
        ExpandableSyms()
        {
            for (const auto& sym : kSymbolExpansion)
            {
                if (sym.second.size() > 1)
                {
                    value.push_back(sym.first);
                }
            }
        }
        string value;
    } expandableSyms;
    return s.find_first_of(expandableSyms.value) != string::npos;
}

bool checkIfNucleotideReferenceSequence(const std::string& sequence) {
    for (char symbol : sequence)
    {
        if (!checkIfNucleotideReferenceSymbol(symbol))
        {
            return false;
        }
    }
    return true;
}

static bool checkIfReferenceSymbol(char symbol) { return kSymbolExpansion.find(symbol) != kSymbolExpansion.end(); }

bool checkIfReferenceSequence(const std::string& sequence) {
    for (char symbol : sequence) {
        if (!checkIfReferenceSymbol(symbol)) {
            return false;
        }
    }
    return true;
}
}
