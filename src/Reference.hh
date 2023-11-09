//
// varViewer
// Date: 2023-08-10
//
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

#pragma once

#include <iostream>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>

// Include the fai class from samtools
#include "htslib/faidx.h"

namespace varViewer
{
class Reference {
public:
    Reference(const faidx_t* faidx_): htsFastaIndexPtr_(faidx_){
        chr_idx_size();
    };

    void chr_idx_size();
    std::string getSequence(const std::string& chrName, int64_t start, int64_t end) const;
    std::string getSequence(const std::string& chrName) const;

private:
    const faidx_t* htsFastaIndexPtr_;
    std::vector<int64_t> chr_size_;
    std::unordered_map<std::string, int32_t> chr_index_;
};
}