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
    Read(int score,bool strand, int count)
    : score(score)
    , strand(strand)
    , count(count)
    {
        score=0, strand=true, count=0;
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
class Aligns {
    public:
    Aligns(int flank, Reference& seq, VariantSpecific& var, Coverage& cover, int mapqual, htsFile* htsFilePtr, bam_hdr_t* htsHeaderPtr, hts_idx_t* htsIndexPtr)
    : flank_(flank), seq_(seq), var_(var), cover_(cover), mapqual_(mapqual)
    , htsFilePtr(htsFilePtr)
    , htsHeaderPtr(htsHeaderPtr)
    , htsIndexPtr(htsIndexPtr)
    {
        set_filter_flag();
        if( var_.var_type == VariantType::SV ){
            std::cerr<<"Start FusionType"<<std::endl;
            AlignsFusion();
            //std::cerr<<"FusionType:"<<std::to_string(int(fusion_type_))<<std::endl;
            //AlignsFusion(pe_record);
        }else {
            refer_start_ = (var_.genomic_pos.front().start_ > flank_) ? (var_.genomic_pos.front().start_ - flank_ ) : 1;
            if( var_.var_type == VariantType::CNV || var_.var_type == VariantType::LOSS || var_.var_type == VariantType::DUP || var_.var_type == VariantType::LOSSDUP ){
                std::cerr<<"CNV analysis"<<std::endl;
                AlignCNV();
            }else{
                AlignSNV();
            }
        }
        sort_AlignReads();
        update_offset();
        handle_PE_overlap();
    }

    void AlignsFusion();
    void AlignCNV();
    void reverseRead(Read& read_);
    void add_coverage(int, int, int, int, int);
    void set_filter_flag();
    void print_Frag();
    void sort_AlignReads();
    bool check_dup(Read& read, Read& mate,bool add=true);
    void handle_PE_overlap();
    std::list<std::string> splitAtBasePosition(Alignment& align, std::list<std::string>& str, size_t pos, bool first=false);
    int chr_Name_Idx(const std::string& chrName);
    void AlignSNV();
    void min_start_offset();
    void max_end_offset();
    void update_offset(){
        min_start_offset();
        max_end_offset();
    }
    int get_read_count(){ return AlignReads.size();}
    int get_start_offset(){ return start_offset_; } // index from refer_start_
    int get_end_offset() { return end_offset_; } // index from start_offset_
    std::string get_chr(){ return var_.genomic_pos.front().chrName;}
    std::string get_end_chr(){return var_.genomic_pos.back().chrName;}
    int64_t get_start() {return var_.genomic_pos.front().start_; }
    int64_t get_end() { return var_.genomic_pos.back().end_; }
    int get_refer_start() {return refer_start_; }
    std::string get_ref_base() {return var_.ref;}
    std::string get_alt_base(){return var_.allele;}
    std::string get_mutant_info(){
        return var_.genomic_pos.front().chrName+":"+std::to_string(var_.genomic_pos.front().start_)+":"+var_.ref+">"+var_.allele;
    }
    FusionType get_fusion_type(){return fusion_type_;}
    //std::string getSequence(){ seq_.getSequence(chr_); }
    Reference get_reference(){ return seq_; }
    //std::string getSequence(std::string chrName, int64_t ss, int64_t ee){ seq_.getSequence(chrName,ss, ee); }
    int get_flank(){return flank_;}
    VariantSpecific get_VariantSpecific(){return var_;}

    std::vector<Frag> AlignReads;
    private:
        int flank_ = 0;
        int64_t refer_start_ = 0; // the most left position
        VariantSpecific& var_;
        Coverage& cover_;
        Reference& seq_;
        int mapqual_=20;
        htsFile* htsFilePtr=nullptr;
        bam_hdr_t* htsHeaderPtr=nullptr;
        hts_idx_t* htsIndexPtr=nullptr;
        int filter_ ;
        int start_offset_ = 0;  // index from refer_start_
        int end_offset_ = 0;   // length from start_offset_
        FusionType fusion_type_;
};
std::ostream& operator<<(std::ostream& os, const Read& read);
std::ostream& operator<<(std::ostream& os, Aligns& align);
}