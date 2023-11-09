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

#pragma once

#include <map>
#include <memory>
#include <string>
#include <vector>
#include "LinearAlign.hh"
#include "Reference.hh"
#include "Variant.hh"
#include "Coverage.hh"

extern "C"
{
    #include "htslib/faidx.h"
    #include "htslib/hts.h"
    #include "htslib/sam.h"
}

namespace varViewer
{
enum class FusionType {
    // first: first break point
    // first_left: first break point left part
    FLSR=0,   // first_left + second_right
    SLFR,   // second_left + first_right
    FLrevSL, // first_left + revComp second_left
    revSRFR, // revComp second_right + first_right
    //revFRSR, // revComp first_right + second_right, equal to revSRFR
    //SLrevFL  // second_left + revComp first_left, equal to FLrevSL
};
struct Read {
    Read()
    {
        score=0, strand=true,count=0;
    }
    Read(int score,bool strand, int count)
    : score(score)
    , strand(strand)
    , count(count)
    {
    }
    Read(std::list<std::string>bases, int score, bool strand, int count, Alignment align)
        : bases(std::move(bases))
        , score(score)
        , strand(std::move(strand))
        , count(std::move(count))
        , align(std::move(align))
    {
    }
    std::list<std::string> bases; // ' ': same as reference, quality< 20 lower, quality>20: upper.
    int score;
    bool strand;
    int count;
    Alignment align;
};
struct Frag {
    Frag(Read& read_, Read& mate_)
        : read(read_)
        , mate(mate_)
    {
    }
    Read read;
    Read mate;
};
struct Fusion {
    Fusion(std::vector<Frag>& AlignReads, FusionType fusion_type)
    : AlignReads_(AlignReads)
    , fusion_type_(fusion_type)
    {
    }

    std::vector<Frag> AlignReads_;
    FusionType fusion_type_;
}
class SV {
    public:
    SV(int flank, Reference& seq, VariantSpecific& var, int mapqual, htsFile* htsFilePtr, bam_hdr_t* htsHeaderPtr, hts_idx_t* htsIndexPtr)
    : flank_(flank), seq_(seq), var_(var), mapqual_(mapqual)
    , htsFilePtr(htsFilePtr)
    , htsHeaderPtr(htsHeaderPtr)
    , htsIndexPtr(htsIndexPtr)
    {
        std::cerr<< "Start FusionType" << std::endl;
        AlignsFusion();
        sort_AlignReads();
        update_offset();
        handle_PE_overlap();
    }

    void AlignsFusion();
    void sort_AlignReads();
    void handle_PE_overlap();
    void min_start_offset();
    void max_end_offset();
    void update_offset(){
        min_start_offset();
        max_end_offset();
    }

    std::vector<Fusion> fusion_;
    private:
        int flank_ = 0;
        int64_t refer_start_ = 0; // the most left position
        VariantSpecific& var_;
        Reference& seq_;
        int mapqual_= 20;
        htsFile* htsFilePtr=nullptr;
        bam_hdr_t* htsHeaderPtr=nullptr;
        hts_idx_t* htsIndexPtr=nullptr;
        int filter_ = BAM_FSECONDARY | BAM_FQCFAIL | BAM_FDUP;
        int start_offset_ = 0;  // index from refer_start_
        int end_offset_ = 0;   // length from start_offset_
};
std::ostream& operator<<(std::ostream& os, const Read& read);
std::ostream& operator<<(std::ostream& os, Aligns& align);
}