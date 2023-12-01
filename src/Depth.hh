//
// varViewer
// Date: 2023-09-05
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
#include <map>
#include <vector>
#include <stdexcept>

extern "C"
{
    #include "htslib/faidx.h"
    #include "htslib/hts.h"
    #include "htslib/sam.h"
}
using std::string;
using std::vector;
namespace varViewer {
struct region {
    region(int start, int end)
        : start_(start)
        , end_(end)
    {
    }

    int start_;
    int end_;
};

class Depth {
  public:
    Depth(int flank, int mapqual, std::string bed_file,std::string bam_list, std::string bam_file)
    : flank_(flank), mapqual_(mapqual), bed_file_(bed_file), bam_list_(bam_list)
    , bam_file_(bam_file)
    {
        readBed();
        all_bam.push_back(bam_file_);
        add_other_bam();
        calDepth();
        normalized_factor();
    }
    Depth(std::string bed_file,std::string bam_list, std::string bam_file)
    : bed_file_(bed_file)
    , bam_list_(bam_list)
    , bam_file_(bam_file)
    {
        readBed();
        all_bam.push_back(bam_file_);
        add_other_bam();
        calDepth();
        normalized_factor();
    }

    void readBed();
    void add_other_bam();
    void calDepth();
    void normalized_factor();
    std::vector<float> get_factors(){return factors_;}
    std::vector<int> get_total_depth(){return total_depth;}
    std::vector<std::string> get_bam_list(){return all_bam; }

  private:
    int flank_ = 500;
    int mapqual_ = 20;
    string bed_file_;
    string bam_list_;
    string bam_file_;
    samFile *input_sf=NULL;
    bam_hdr_t *hdr=NULL;
    hts_itr_t *iter=NULL;
    hts_idx_t *index=NULL;
    int  filter_ = BAM_FSECONDARY | BAM_FQCFAIL | BAM_FDUP;
    std::map<std::string, vector<region>> regions;
    std::vector<std::string> all_bam;
    std::vector<int> total_depth;
    std::vector<float> factors_;
};
}