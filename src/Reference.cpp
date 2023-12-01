//
// varViewer
// Date: 2023-08-10
// Author: Huang Fei <huangfei@genomics.cn>
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.
//

#include <algorithm>
#include <memory>
#include <stdexcept>
#include "Reference.hh"

using std::string;
using std::to_string;
using std::vector;

namespace varViewer
{
void Reference::chr_idx_size() {
    for (int contigIndex = 0; contigIndex != faidx_nseq(htsFastaIndexPtr_); ++contigIndex) {
        const char* sequenceName = faidx_iseq(htsFastaIndexPtr_, contigIndex);
        int sequenceLength = faidx_seq_len(htsFastaIndexPtr_, sequenceName);
        chr_size_.emplace_back(sequenceLength);
        chr_index_.emplace(std::make_pair(sequenceName, contigIndex));
    }
}
string Reference::getSequence(const string& chrName, int start, int end) const {
    int extractedLength;
    char* sequencePtr = faidx_fetch_seq(htsFastaIndexPtr_, chrName.c_str(), start - 1, end, &extractedLength);

    if (!sequencePtr || extractedLength < 0 || extractedLength < end - start) {
        const string encoding(chrName + ":" + to_string(start) + "-" + to_string(end));
        const string message = "Unable to extract " + encoding + " from genome sequence";
        throw std::logic_error(message);
    }

    string sequence(sequencePtr);
    free(sequencePtr);
    std::transform(sequence.begin(), sequence.end(), sequence.begin(), ::toupper);

    return sequence;
}
string Reference::getSequence(const string& chrName) const {
    const int contigIndex = chr_index_.at(chrName);
    return getSequence(chrName, 0, chr_size_[contigIndex]);
}
}
