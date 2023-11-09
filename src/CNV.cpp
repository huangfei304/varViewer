//
// varViewer
// Date: 2023-08-09
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

#include <vector>
#include <boost/algorithm/string.hpp>
#include <iostream>
#include <sstream>
#include <zlib.h>
#include "Aligns.hh"
#include "utils.hh"
#include "CNV.hh"
extern "C"
{
#include "htslib/faidx.h"
#include "htslib/hts.h"
#include "htslib/sam.h"
}

namespace varViewer
{
using std::string;
using std::vector;
using std::to_string;
using std::logic_error;

void CNV::readBam() {
    const string chr_ = var_.genomic_pos.front().chrName;
    const int64_t start_ = var_.genomic_pos.front().start_;
    const int64_t end_ = var_.genomic_pos.back().end_;
    const char *mode = "rb";
    //std::vector<std::string> bam_list_ = depth_.get_bam_list();
    for(int i = 0; i <bam_list_.size(); i++) {
        if ((input_sf = sam_open(bam_list_[i].c_str(), mode)) == NULL) {
            fprintf(stderr, "Failed to open input file %s.\n", bam_list_[i].c_str());
            exit(EXIT_FAILURE);
        }
        if ((hdr = sam_hdr_read(input_sf)) == NULL) {
            fprintf(stderr, "Failed to read header for input file %s.\n", bam_list_[i].c_str());
            exit(EXIT_FAILURE);
        }
        if ((index = sam_index_load(input_sf, bam_list_[i].c_str())) == NULL) {
            fprintf(stderr, "Failed to read index for input file %s.\n", bam_list_[i].c_str());
            exit(EXIT_FAILURE);
        }
        fprintf(stderr, "[CNV] Start read Bam file: %s.\n", bam_list_[i].c_str());
        const int chr_idx = sam_hdr_name2tid(hdr,chr_.c_str());
        int32_t chr_len = (int32_t)hdr->target_len[chr_idx];
        refer_end_ = ( (end_ + flank_) > chr_len ) ? chr_len : (end_+flank_);
        hts_itr_t *iter = sam_itr_queryi(index, chr_idx, refer_start_, refer_end_);
        if (iter == nullptr) {
            string region= chr_ + ":"+ to_string(refer_start_) + "-" + to_string(refer_end_);
            throw logic_error("Failed to extract reads from the specified region:"+region);
        }

        std::vector<int> cover_;
        for(int j=0;j<(refer_end_ - refer_start_+1);j++){
            cover_.push_back(0);
        }
        bam1_t *rec = bam_init1();
        while((sam_itr_next(input_sf,iter,rec)) >=0 ) {
            if( rec->core.flag & BAM_FUNMAP ) continue;
            if( rec->core.flag & filter_ ) continue;
            if( rec->core.qual < (uint8_t) mapqual_ ) continue;

            int start,end,refer_start_offset,c_len,c_op,rd_pos,ref_pos;
            uint32_t *cigar;

            start = rec->core.pos;
            refer_start_offset = start - refer_start_;
            if ( refer_start_offset < 0 ) continue;
            //fprintf(stderr, "read name: %s.\n", bam_get_qname(rec));
            cigar = bam_get_cigar(rec); // CIGAR
            ref_pos=rd_pos=0;
            for (int j = 0; j < rec->core.n_cigar; ++j) {
                c_op = bam_cigar_op(cigar[j]);
                c_len = bam_cigar_oplen(cigar[j]);
                switch (c_op) {
                    case BAM_CMATCH: // M
                    case BAM_CEQUAL: // =
                    case BAM_CDIFF: // X
                        for(int x=0;x<c_len;x++) cover_[refer_start_offset+ref_pos+x]++;
                        break;
                    case BAM_CDEL: // D
                            cover_[refer_start_offset+ref_pos]++;
                        break;
                    /** case BAM_CHARD_CLIP:
                    *    cigar_ = "H";
                    *    break;
                    * case BAM_CPAD:
                    *    cigar_ = "P";
                    *    break;
                    **/
                    case BAM_CINS: // CNV not consider I
                    case BAM_CREF_SKIP: // CNV not consider N
                    case BAM_CSOFT_CLIP: // CNV not consider S
                    default:
                        break;
                }
                if((bam_cigar_type(cigar[j]) & 0x02) !=0 ){
                    ref_pos += c_len;
                }
            }
        }
        bam_destroy1(rec);
        hts_itr_destroy(iter);

        //fprintf(stderr,"AAAAAAAAAAAAAAAAAAAAAA\n");
        //std::cerr<<"refer_start:"<<std::to_string(refer_start_)<<",max:"<<std::to_string(max)<<std::endl;
        //for(int k=0;k<cover_.size();k++){
        //    std::cerr<<"Index:"<<std::to_string(k)<<",depth:"<<cover_[i]<<std::endl;
        //}
        CNVinfo this_info(refer_start_, cover_);
        cnv_info_.push_back(this_info);

        fprintf(stderr, "[CNV] Finish read Bam file: %s.\n", bam_list_[i].c_str());
    }
}
void CNV::normalizedDepth(){// normalized reference sample depth
    std::cerr<<"[CNV] Start normalized depth..."<<std::endl;
    int size = (refer_end_ - refer_start_ + 1);
    for(int i=0;i<bam_list_.size();++i){
        float this_factor = factors_[i];
        for(int j=0;j<size;++j){
            int dp = cnv_info_[i].depth_[j];
            if( dp > 0 ){
                float this_dp = (dp * this_factor);
                if( this_dp < 1 ){
                    dp = 1;
                }else{
                    dp = int(this_dp+0.5);
                }
            }
            cnv_info_[i].depth_[j] = dp;
        }
    }
    std::cerr<<"[CNV] End normalized depth..."<<std::endl;
}
void CNV::calCoverage(){// calculate coverage
    std::cerr<<"[CNV] Start calculate coverage..."<<std::endl;
    for(auto it=var_.genomic_pos.begin();it!=var_.genomic_pos.end();++it){
        int64_t start = it->start_;
        int64_t end = it->end_;
        int start_offset, end_offset, exon_size, this_cover, this_dcover;
        start_offset = start - refer_start_;
        end_offset = end - refer_start_;
        exon_size = (end - start + 1);
        this_cover = 0;
        this_dcover = 0;
        for(int k=start_offset;k<=end_offset;++k){
            if( cnv_info_[0].depth_[k] > 0 ) this_cover++;
            if( cnv_info_[0].depth_[k] > dcutoff_ ) this_dcover++;
        }
        cover_1.push_back(1.0*(this_cover/exon_size));
        cover_30.push_back(1.0*(this_dcover/exon_size));
    }
    //for(int i=0;i<cover_1.size();i++){
    //    std::cerr<<"Index: "<<std::to_string(i)<<", Coverage_1: "<<std::to_string(cover_1[i])<<", Coverage_30: "<<std::to_string(cover_30[i])<<std::endl;
    //}
    std::cerr<<"[CNV] End calculate coverage..."<<std::endl;
}
void CNV::averageDepth(){// reference average depth
   std::cerr<<"[CNV] Start reference average depth..."<<std::endl;
   int size = (refer_end_ - refer_start_ + 1);
   int num = bam_list_.size()-1;
   std::vector<int> this_cover;
   for(int i=0;i<size;++i){
        int total = 0;
        for(int j = 1; j <=num; ++j ){
            total += cnv_info_[j].depth_[i];
        }
        total = int(1.0*total/num+0.5);
        this_cover.push_back(total);
   }
   CNVinfo this_info(refer_start_,this_cover);
   cnv_info_.push_back(this_info);
   std::cerr<<"[CNV] End reference average depth..."<<std::endl;
}
void CNV::readMapability(){
    if ( mapability_ !="" ){
        std::cerr<<"[CNV] Start read Mapability file:"+mapability_<<std::endl;
        gzFile fp=NULL;
        std::string line;
        std::string chr = var_.genomic_pos.front().chrName;
        std::vector<float> this_map((refer_end_ - refer_start_ + 1),0);
        if(Z_NULL==(fp=gzopen(mapability_.c_str(),"r"))){
            fprintf(stderr,"Error opening file: %s\n",mapability_);
            exit(0);
        }

    const int LENS = 1024;
    char buf[LENS];
    const char *delims = "\t \n";
     while(gzgets(fp,buf,LENS)!=NULL){
        char* map_chr = strdup(strtok(buf,delims));
        std::string chr_s(map_chr, strlen(map_chr));
        int map_start = atoi(strtok(NULL,delims));
        int map_end = atoi(strtok(NULL,delims));
        float map_score = atof(strtok(NULL,delims));
        if ( chr_s != chr ) continue;
        if( (map_start >= refer_start_ && map_start <= refer_end_) || (map_end >= refer_start_ && map_end <= refer_end_) ){
            int min = (map_start < refer_start_ ) ? refer_start_ : map_start;
            int max = (map_end < refer_end_ ) ?  map_end : refer_end_;
            int len = max - min+1;
            int idx = min - refer_start_;
            for(int i=0;i<len;++i){
                this_map[i+idx] = map_score;
            }
        }
     }
     gzclose(fp);
     map_score_ = this_map;
     std::cerr<<"[CNV] Finish read Mapability file:"+mapability_<<std::endl;
    }
}
void CNV::readGT(){
    if ( var_file_ !="" ){
        std::cerr<<"[CNV] Start read variant file:"+var_file_<<std::endl;
        std::string line;
        std::string chr = var_.genomic_pos.front().chrName;
        std::ifstream varFile(var_file_);
        if (varFile.is_open()) {
            while (getline(varFile, line)){
                std::vector<std::string> pieces;
                boost::split(pieces, line, boost::is_any_of("\t"));
                //vector<string> pieces = splitStringByDelimiter(line,'\t');
                if (pieces[0] == "#Chr" || pieces[0] != chr ) continue;
                int64_t pos;
                std::stringstream sstream;
                sstream<<pieces[1];sstream>>pos;
                if( pos >= refer_start_ && pos <= refer_end_ ){
                    std::string info = pieces[0]+":"+pieces[1]+":"+pieces[2]+">"+pieces[3]+":"+pieces[4];
                    GTinfo this_gt(pos, info);
                    gt_info_.push_back(this_gt);
                }
            }
            varFile.close();
            std::cerr<<"[CNV] Finish read variant file:"+var_file_<<std::endl;
        } else{
            throw std::runtime_error("Unable to open variantion file: " + var_file_);
        }
    }

    for (int i=0; i<gt_info_.size();++i){
        std::string flag = "intron";
        for(auto it=var_.genomic_pos.begin();it!=var_.genomic_pos.end();++it){
            int start = it->start_;
            int end = it->end_;
            if( gt_info_[i].pos_ < start ) break;
            if (gt_info_[i].pos_ > end ) continue;
            flag = "exon";
        }
        gt_info_[i].infor_ += ":"+flag;
    }

    //for(auto gt_: gt_info_){
    //    std::cerr<<"pos:"<<std::to_string(gt_.pos_)<<",infor:"<<gt_.infor_<<std::endl;
    //}
}
void CNV::exonGC(){
     std::cerr<<"[CNV] Start calculate the GC of the CNV exon."<<std::endl;
    const string chr_ = var_.genomic_pos.front().chrName;
    string ref_seq = seq_.getSequence(chr_);
    std::transform(ref_seq.begin(), ref_seq.end(), ref_seq.begin(), ::toupper);
    std::vector<float> gc_rate;

    for(auto it=var_.genomic_pos.begin();it!=var_.genomic_pos.end();++it){
        int64_t start = it->start_;
        int64_t end = it->end_;
        int total = (end - start + 1);
        int gc_count = 0;
        for (int i=start;i<=end; ++i){
            if( ref_seq[i]== 'G' || ref_seq[i] == 'C' ){
                gc_count++;
            }
        }
        float this_gc_rate = 1.0 * gc_count / total;
        gc_rate.push_back(this_gc_rate);
    }
    gc_ = gc_rate;
    std::cerr<<"[CNV] Finish calculate the GC of the CNV exon."<<std::endl;
}
void CNV::mxDepth(){
    std::cerr<<"[CNV] Start min and max exon Depth..."<<std::endl;
    for(int i=0;i<=bam_list_.size();++i){
        int min = 1000000;
        int max = 0;
        // just the test and the reference set mean depth
        if( i !=0 && i != bam_list_.size() ) continue;
        for(auto it=var_.genomic_pos.begin();it!=var_.genomic_pos.end();++it){
            int64_t start = it->start_;
            int64_t end = it->end_;
            int start_offset = start - refer_start_;
            int end_offset = end - refer_start_;
            for(int k=start_offset;k<=end_offset;++k){
                if( cnv_info_[i].depth_[k] > max ) max = cnv_info_[i].depth_[k];
                if( cnv_info_[i].depth_[k] < min ) min = cnv_info_[i].depth_[k];
            }
        }
        if( min_depth_ > min ) min_depth_ = min;
        if( max_depth_ < max ) max_depth_ = max;
    }

    if( min_depth_ >= max_depth_ ){
        std::cerr<<"The exon min depth is more than max depth, please check..."<<std::endl;
        exit(EXIT_FAILURE);
   }
   std::cerr<<"[CNV] min depth in all sample: "<<std::to_string(min_depth_)<<std::endl;
   std::cerr<<"[CNV] max depth in all sample: "<<std::to_string(max_depth_)<<std::endl;
   std::cerr<<"[CNV] End min and max exon Depth..."<<std::endl;
}
}

