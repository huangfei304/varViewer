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
#include <fstream>
#include "LinearAlign.hh"
#include "Reference.hh"
#include "Variant.hh"
//#include "Coverage.hh"
#include "Depth.hh"

extern "C"
{
    #include "htslib/faidx.h"
    #include "htslib/hts.h"
    #include "htslib/sam.h"
}

namespace varViewer
{
struct CNVinfo {
    CNVinfo(int start_offset, std::vector<int> depth )
    : start_on_chr(start_offset)
    , depth_(depth)
    {
    }
    std::vector<int> depth_;
    int start_on_chr;
};
struct GTinfo {
    GTinfo(int pos, std::string infor)
    : pos_(pos)
    , infor_(infor)
    {
    }
    int pos_;
    std::string infor_;
};
class CNV {
    public:
    CNV(int flank, VariantSpecific& var, int mapqual, std::vector<std::string> bam_list, std::vector<float>factors, std::string mapability, Reference& seq,std::string var_file)
    //CNV(int flank, VariantSpecific& var, int mapqual, Depth& depth, std::string mapability, Reference& seq)
    : flank_(flank), var_(var), mapqual_(mapqual)
    , bam_list_(bam_list)
    , factors_(factors)
    //, depth_(depth)
    ,mapability_(mapability)
    ,seq_(seq)
    ,var_file_(var_file)
    {
        //std::cerr<<"CNV analysis begin"<<std::endl;
        refer_start_ = (var_.genomic_pos.front().start_ > flank_) ? (var_.genomic_pos.front().start_ - flank_ ) : 1;
        //bam_list_ = depth_.get_bam_list();
        //factors_ = depth_.get_factors();
        //for(int i = 0; i <factors_.size(); ++i ){
        //    std::cerr<<"bam file:"<<bam_list_[i]<<std::endl;
        //    std::cerr<<"factor:"<<std::to_string(factors_[i])<<std::endl;
       //}
        readBam();
        normalizedDepth();
        calCoverage();
        averageDepth();
        readMapability();
        exonGC();
        mxDepth();
        readGT();
    }
    CNV(VariantSpecific& var, std::string mapability, Reference& seq)
    //CNV(VariantSpecific& var, Depth& depth, std::string mapability, Reference& seq)
    : var_(var)
    //, depth_(depth)
    , mapability_(mapability)
    , seq_(seq)
    {
        //std::cerr<<"CNV analysis begin"<<std::endl;
        refer_start_ = (var_.genomic_pos.front().start_ > flank_) ? (var_.genomic_pos.front().start_ - flank_ ) : 1;
        readBam();
        normalizedDepth();
        calCoverage();
        averageDepth();
        readMapability();
        exonGC();
        mxDepth();
        readGT();
    }

    void readBam();
    void normalizedDepth();
    void calCoverage();
    void averageDepth();
    void readMapability();
    void exonGC();
    void mxDepth();
    void readGT();
    //CNVinfo get_CNVinfo(int i){ return cnv_info_[i];}
    int get_CNVinfo_size() {return cnv_info_.size();}
    VariantSpecific get_VariantSpecific(){return var_;}
    std::string get_chr(){return var_.genomic_pos.front().chrName;}
    std::vector<float> get_gc(){ return gc_;}
    std::vector<float> get_cover(int i=1){if( i==1 ){return cover_1;}else{return cover_30;}};
    std::vector<float> get_mapability(){return map_score_;}
    int get_ref_start() {return refer_start_;}
    int get_cnv_size() {return (refer_end_ - refer_start_+1);}
    int get_min_depth() { return min_depth_;}
    int get_max_depth() { return max_depth_;}

    std::vector<CNVinfo> cnv_info_;
    std::vector<GTinfo> gt_info_;
    private:
        int flank_ = 500;
        int64_t refer_start_ = 0; // the most left position
        int64_t refer_end_ = 0; // the most right position
        VariantSpecific& var_;
        Reference& seq_;
        int mapqual_= 20;
        samFile *input_sf=NULL;
        bam_hdr_t *hdr=NULL;
        //hts_itr_t *iter=NULL;
        hts_idx_t *index=NULL;
        int filter_ = BAM_FSECONDARY | BAM_FQCFAIL | BAM_FDUP;
        //Depth& depth_;
        std::vector<std::string> bam_list_;
        std::vector<float> factors_;
        std::string mapability_="";
        int dcutoff_ = 30; // depth cutoff
        std::vector<float> gc_;
        std::vector<float> cover_1;
        std::vector<float> cover_30;
        std::vector<float> map_score_;
        int min_depth_ = 10000;
        int max_depth_ = 0;
        std::string var_file_ =""; // variant genotye file

};
}