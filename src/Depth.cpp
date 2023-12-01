// varViewer
// Date: 2023-09-05
// Author: Huang Fei <huangfei@genomics.cn>
//

#include "Depth.hh"
#include <sstream>
#include <fstream>
#include <iostream>
#include <boost/algorithm/string.hpp>
#include <map>
#include <vector>
#include <stdexcept>

namespace varViewer {

void Depth::readBed(){
    std::string line;
    std::ifstream bedFile(bed_file_);
    if (bedFile.is_open()) {
        while (getline(bedFile, line)){
            std::vector<std::string> pieces;
            boost::split(pieces, line, boost::is_any_of("\t"));
            //if( chr_ != "all" && chr_ != pieces[0] ) continue;
            int sstart,eend;
            std::stringstream sstream,estream;
            sstream<<pieces[1];sstream>>sstart;
            estream<<pieces[2];estream>>eend;

            region this_region(sstart, eend);
            std::map<std::string,std::vector<region>>::iterator iter=regions.find(pieces[0]);
            if( iter != regions.end()){
                iter->second.push_back(this_region);
            }else{
                std::vector<region> tmp_regions;
                tmp_regions.push_back(this_region);
                regions.emplace(pieces[0], tmp_regions);
            }
        }
        bedFile.close();
    }else{
        throw std::logic_error("Unable to open BED file: " + bed_file_);
    }
    std::cerr<<"Finish read BED file: "+bed_file_<<std::endl;
}
void Depth::add_other_bam(){
    std::string line;
    std::ifstream bamList(bam_list_);
    if (bamList.is_open()) {
        while (getline(bamList, line)){
            all_bam.push_back(line);
        }
        bamList.close();
    }else{
        throw std::logic_error("Unable to open bam list file: " + bam_list_);
    }
}
void Depth::normalized_factor(){
    int depth = total_depth[0];
    for(int i=0;i<total_depth.size();++i){
        factors_.push_back((1.00 * depth/total_depth[i]));
        std::cerr<<"Index: "<<std::to_string(i)<<", factor:"<<std::to_string(factors_[i])<<std::endl;
    }
}
void Depth::calDepth(){
    const char *mode = "rb";
    for(auto it = all_bam.begin(); it !=all_bam.end(); ++it ){
        if ((input_sf = sam_open((*it).c_str(), mode)) == NULL) {
            fprintf(stderr, "Failed to open input file \"%s\".\n", (*it).c_str());
            exit(EXIT_FAILURE);
        }
        if ((hdr = sam_hdr_read(input_sf)) == NULL) {
            fprintf(stderr, "Failed to read header for input file \"%s\".\n", (*it).c_str());
            exit(EXIT_FAILURE);
        }
        if ((index = sam_index_load(input_sf, (*it).c_str())) == NULL) {
            fprintf(stderr, "Failed to read index for input file \"%s\".\n", (*it).c_str());
            exit(EXIT_FAILURE);
        }
        fprintf(stderr, "[Depth] Start read Bam file: %s.\n", (*it).c_str());
        int this_total_depth=0;
        for(auto iter=regions.begin();iter != regions.end();++iter){
            for(auto i=iter->second.begin();i !=iter->second.end();++i){
                int start = i->start_;
                int end = i->end_;

                int this_depth = 0;
                int chr_idx = sam_hdr_name2tid(hdr,iter->first.c_str());
                bam1_t *rec = bam_init1();
                hts_itr_t *rec_iter = sam_itr_queryi(index, chr_idx, start-flank_, end+flank_);
                if (rec_iter == nullptr) {
                    string region= iter->first + ":"+ std::to_string(start) + "-" + std::to_string(end);
                    throw std::logic_error("Failed to extract reads from the specified region:"+region);
                }
                while((sam_itr_next(input_sf,rec_iter,rec)) >=0 ) {
                    if( rec->core.flag & BAM_FUNMAP ) continue;
                    if( rec->core.flag & filter_ ) continue;
                    if( rec->core.qual < (uint8_t) mapqual_ ) continue;

                    int pos = rec->core.pos;
                    int pos_end = pos;
                    uint32_t *cigar = bam_get_cigar(rec); // CIGAR
                    for (int i = 0; i < rec->core.n_cigar; ++i) {
                        int c_op = bam_cigar_op(cigar[i]);
                        int c_len = bam_cigar_oplen(cigar[i]);
                        if((bam_cigar_type(cigar[i]) & 0x02) !=0 ){
                            pos_end += c_len;
                        }
                    }
                    if( pos <= end  && pos_end >= start ){
                        int min = (pos < start ) ? start : pos;
                        int max = (pos_end < end ) ?  pos_end : end;
                        int len = max - min+1;
                        this_depth += len;
                    }
                }
                bam_destroy1(rec);
                hts_itr_destroy(rec_iter);
                this_total_depth += this_depth;
            }
        }
        fprintf(stderr, "[Depth] Finish read Bam file: %s.\n", (*it).c_str());
        total_depth.push_back(this_total_depth);
    }
}
}
