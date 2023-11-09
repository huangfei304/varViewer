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
#include "Aligns.hh"
#include "utils.hh"
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

void Aligns::reverseRead(Read& this_read){
    std::list<string> this_bases;
     for (std::list<string>::reverse_iterator it = this_read.bases.rbegin(); it != this_read.bases.rend();++it) {
        string this_base = reverseComplement(*it);
        this_bases.push_back(this_base);
     }
     this_read.bases = this_bases;
     this_read.align.reverseOperation();
     this_read.strand = (!this_read.strand);
}
bool softclip_match(string first, string second, int len, int mismatch, bool right=false){
    int mis = 0;
    int ss_len = (second.size()-1);
    int i, j;
    if( right ){
        for(i=(len-1),j=ss_len;i>=0;--i,--j){
            if( first[i] != second[i]) mis++;
        }
    }else{
        for(i=0;i<len;++i){
            if(first[i] != second[i]) mis++;
        }
    }
    return (mis <= mismatch);
}
void Aligns::set_filter_flag(){
    /* Secondary read or QC failure */
    filter_ = BAM_FSECONDARY | BAM_FQCFAIL;
    if( var_.var_type == VariantType::LOSS || var_.var_type == VariantType::DUP || var_.var_type == VariantType::LOSSDUP ){
         /* CNV Duplicate reads */
        filter_ = BAM_FSECONDARY | BAM_FQCFAIL | BAM_FDUP;
    }
}
void Aligns::sort_AlignReads (){
    std::sort(AlignReads.begin(), AlignReads.end(), [](Frag const& a, Frag const& b) -> bool {
        if( a.read.align.referenceStart()==b.read.align.referenceStart()){
            return ((a.read.align.referenceStart()+a.read.align.referenceLength())< (b.read.align.referenceStart() + b.read.align.referenceLength()));
        }else{
            return(a.read.align.referenceStart() < b.read.align.referenceStart());
        }
    });
}
void Aligns::print_Frag(){
    int total=0;
    for(int i=0;i<AlignReads.size();++i){
        std::cerr<<"First read: "<<AlignReads[i].read.align<<", Strand:"<<AlignReads[i].read.strand<<", Count:"<<AlignReads[i].read.count<<std::endl;

        for(auto j=AlignReads[i].read.bases.begin();j != AlignReads[i].read.bases.end(); ++j){
            std::cerr<<*j;
        }
        std::cerr<<std::endl;
        total += AlignReads[i].read.count;
        if( AlignReads[i].mate.count > 0){
            std::cerr<<"Mate read: "<<AlignReads[i].mate.align<<", Strand:"<<AlignReads[i].mate.strand<<", Count:"<<AlignReads[i].mate.count<<std::endl;
        }
    }
    std::cerr<<"Total read:"<<std::to_string(total)<<std::endl;
}
std::list<string> Aligns::splitAtBasePosition(Alignment& align, std::list<string>& bases, size_t reference_position, bool first) {
    const size_t end_of_reference_positions = align.referenceStart() + align.referenceLength();
    if (reference_position == 0 || end_of_reference_positions <= reference_position) {
        std::ostringstream os;
        os <<align;
        throw logic_error("Cannot split " + os.str() + " for bases, at reference position " + to_string(reference_position));
    }
    if( align.size() != bases.size()){
        throw logic_error("The Alignment and bases list size is not equal.\n");
    }

    size_t first_unused_position = align.referenceStart();
    std::list<Operation> operations_ = align.operations();
    std::list<Operation>::const_iterator operation_it = operations_.begin();
    std::list<string>::const_iterator base_it = bases.begin();
    while (operation_it != operations_.end()) {
        const size_t first_unused_position_after_applying_operation
            = first_unused_position + operation_it->referenceLength();
        if (first_unused_position_after_applying_operation <= reference_position){
            ++operation_it;
            ++base_it;
            first_unused_position = first_unused_position_after_applying_operation;
        }else{
            break;
        }
    }

    if (operation_it == operations_.end()) {
        throw logic_error("The must be some error in  splitAtReferencePosition");
    }

    if (first_unused_position == reference_position) {
        std::list<std::string> suffix_bases;
        if( first ){
            suffix_bases.splice(suffix_bases.begin(), bases, bases.begin(), --base_it);
        }else{
            suffix_bases.splice(suffix_bases.begin(), bases, base_it, bases.end());
        }
        return suffix_bases;
    } else {
        const size_t first_piece_reference_length = reference_position - first_unused_position;
        std::list<std::string> suffix_bases;
        if( first ){
            string first_bases = base_it->substr(0, first_piece_reference_length);
            suffix_bases.splice(suffix_bases.begin(), bases, bases.begin(), --base_it);
            suffix_bases.push_back(first_bases);
        }else{
            //std::cerr<<"pos:"<<std::to_string(first_piece_reference_length)<<", all:"<<*base_it<<"#\n";
            string first_bases = base_it->substr(first_piece_reference_length);
            suffix_bases.splice(suffix_bases.begin(), bases, ++base_it, bases.end());
            suffix_bases.push_front(first_bases);
        }
        return suffix_bases;
    }
}

void Aligns::handle_PE_overlap(){
    for (auto iter = AlignReads.begin(); iter < AlignReads.end(); ++iter) {
        if( iter->mate.count == 0) continue;
        int read_end = (int)(iter->read.align.referenceStart() + iter->read.align.referenceLength());
        int mate_start = (int)(iter->mate.align.referenceStart());
        int mate_end = (int)(iter->mate.align.referenceStart() + iter->mate.align.referenceLength());
        if( mate_end <= read_end){
            Read mate_read(0, true, 0);
            iter->mate = mate_read;
        }else if( mate_start <= read_end ){
            //std::cerr<<"read:"<<iter->read<<std::endl<<"mate:"<<iter->mate<<std::endl;
            //std::cerr<<"mate_start:"<<std::to_string(mate_start)<<",read_end:"<<std::to_string(read_end)<<std::endl;
            iter->mate.bases = splitAtBasePosition(iter->mate.align, iter->mate.bases, read_end);
            iter->mate.align = iter->mate.align.splitAtReferencePosition(read_end);
        }
    }
}
bool Aligns::check_dup(Read& read, Read& mate,bool add){
    if( mate.count>0 && read.align.referenceStart() > mate.align.referenceStart()){
        throw logic_error("The must be some error: the bam maybe not sort.\n");
    }
    if (AlignReads.size() ==0  || !add ){
        AlignReads.emplace_back(Frag(read, mate));
        return true;
    }
    if( mate.count > 0 ){// PE reads
        for(int i=0;i<AlignReads.size();++i){
            Read read_ = AlignReads[i].read;
            Read mate_ =AlignReads[i].mate;
            if( mate_.count > 0 &&  read_.align == read.align && mate_.align == mate.align ){
                AlignReads[i].read.count++;
                AlignReads[i].mate.count++;
                return true;
            }
        }
    }else{// SE reads
        for(int i=0;i<AlignReads.size();++i){
            Read read_ = AlignReads[i].read;
            Read mate_ =AlignReads[i].mate;
            if( mate_.count==0 && read_.align == read.align ){
                AlignReads[i].read.count++;
                return true;
            }
        }
    }
    AlignReads.emplace_back(Frag(read, mate));
    return true;
}
int Aligns::chr_Name_Idx(const string& chrName){
    for(int32_t i=0;i< htsHeaderPtr->n_targets;++i){
        if( chrName.c_str() == htsHeaderPtr->target_name[i] ){
            return i;
        }
    }
    return -1;
}
void Aligns::min_start_offset(){
    if (AlignReads.size()==0){
        string region= var_.genomic_pos.front().chrName +":"+to_string(var_.genomic_pos.front().start_)+"-"+to_string(var_.genomic_pos.back().end_);
        throw logic_error("No reads mapped to this site: "+region);
    }
    int min_offset = std::numeric_limits<int>::max();
    for(int i=0;i<AlignReads.size();++i){
        Read tmp_read = AlignReads[i].read;
        if (min_offset > tmp_read.align.referenceStart() ){
            min_offset = tmp_read.align.referenceStart();
        }
    }
    if( min_offset == std::numeric_limits<int>::max() || min_offset < 0 ){
        throw logic_error("The min offset is too big or less than 0.");
    }

    for(int i=0;i<AlignReads.size();++i){
        AlignReads[i].read.align.set_reference_start(AlignReads[i].read.align.referenceStart()-min_offset);
        if( AlignReads[i].mate.count>0){
            AlignReads[i].mate.align.set_reference_start(AlignReads[i].mate.align.referenceStart()-min_offset);
        }
       //AlignReads[i].align = AlignReads[i].align.splitAtReferencePosition(min_offset);
    }
    start_offset_ = min_offset;
}
void Aligns::max_end_offset(){
    if (AlignReads.size()==0){
        string region= var_.genomic_pos.front().chrName +":"+to_string(var_.genomic_pos.front().start_)+"-"+to_string(var_.genomic_pos.back().end_);
        throw logic_error("No reads mapped to this site: "+region);
    }
    int tmp_offset = 0;
    for(int i=0;i<AlignReads.size();++i){
        int end = (int)(AlignReads[i].read.align.referenceStart() + AlignReads[i].read.align.referenceLength());
        if( tmp_offset < end ){
            tmp_offset = end;
        }
        if( AlignReads[i].mate.count>0 ){
            end = (int)(AlignReads[i].mate.align.referenceStart()+AlignReads[i].mate.align.referenceLength());
            if( tmp_offset < end ){
                tmp_offset = end;
            }
        }
    }
    end_offset_ = tmp_offset + start_offset_;
}

void Aligns::AlignSNV() {
    const string chr_ = var_.genomic_pos.front().chrName;
    const int64_t start_ = var_.genomic_pos.front().start_;
    const int64_t end_ = var_.genomic_pos.back().end_;
    VariantType type_  =  var_.var_type;
    const string ref_ = var_.ref;
    const string alt_ = var_.allele;
    const int chr_idx = sam_hdr_name2tid(htsHeaderPtr,chr_.c_str());
    int32_t chr_len = (int32_t)htsHeaderPtr->target_len[chr_idx];
    const int32_t tmp_end = ( (end_ + flank_) > chr_len ) ? chr_len : (end_+flank_);
    hts_itr_t *iter = sam_itr_queryi(htsIndexPtr, chr_idx, refer_start_, tmp_end);
    if (iter == nullptr) {
        string region= chr_ + ":"+ to_string(start_) + "-" + to_string(end_);
        throw logic_error("Failed to extract reads from the specified region:"+region);
    }
    bam1_t *rec = bam_init1();
    std::map<string, Frag>pairedCache;
    string ref_seq = seq_.getSequence(chr_);
    //std::cerr<<ref_seq<<std::endl;
    std::transform(ref_seq.begin(), ref_seq.end(), ref_seq.begin(), ::toupper);
    while((sam_itr_next(htsFilePtr,iter,rec)) >=0 ) {
        if( rec->core.flag & BAM_FUNMAP ) continue;
        if( rec->core.flag & filter_ ) continue;
        if( rec->core.qual < (uint8_t) mapqual_ ) continue;

        bool strand;
        string bases,quals,read_name,this_cigar;
        int start,end,readLength,refer_start_offset,c_len,c_op,rd_pos, ref_pos;
        uint8_t *htsSeqPtr,*htsQualsPtr;
        uint32_t *cigar;
        std::list<string> seqInfo;
        std::list<Operation> operations;

        start = rec->core.pos;
        refer_start_offset = start - refer_start_;
        if ( refer_start_offset < 0 ) continue;
        read_name = bam_get_qname(rec);
        strand = (bam_is_rev(rec) == 0);
        htsSeqPtr = bam_get_seq(rec);// read sequence
        htsQualsPtr = bam_get_qual(rec);// read quality
        cigar = bam_get_cigar(rec); // CIGAR
        readLength = rec->core.l_qseq;
        ref_pos=rd_pos=0;
        string coverage_ins,coverage_del="";
        coverage_ins="";
        for (int i = 0; i < rec->core.n_cigar; ++i) {
            c_op = bam_cigar_op(cigar[i]);
            c_len = bam_cigar_oplen(cigar[i]);
            string cigar_base(c_len,' ');
            switch (c_op) {
                case BAM_CMATCH: // M
                    for(int j=0;j<c_len;++j){
                        if( seq_nt16_str[bam_seqi(htsSeqPtr, (rd_pos+j))] != ref_seq[start+ref_pos+j]){
                            if(htsQualsPtr[rd_pos+j] >= mapqual_ ){
                                cigar_base[j] = seq_nt16_str[bam_seqi(htsSeqPtr, (rd_pos+j))];
                            }else{
                                cigar_base[j] = tolower(seq_nt16_str[bam_seqi(htsSeqPtr, (rd_pos+j))]);
                            }
                            cover_.add_coverage(refer_start_offset+ref_pos+j, VAR_VR);
                        }else{
                            cover_.add_coverage(refer_start_offset+ref_pos+j, VAR_RR);
                        }
                    }
                    operations.push_back(Operation(OperationType::kMatch, (uint32_t)c_len));
                    seqInfo.push_back(cigar_base);
                    break;
                case BAM_CINS: // I
                    // CNV not consider Insertion
                    if( type_ == VariantType::LOSS || type_ == VariantType::DUP || type_ == VariantType::LOSSDUP ) break;
                    operations.push_back(Operation(OperationType::kInsertionToRef, (uint32_t)c_len));
                    for(int j=0;j<c_len;++j){
                        if(htsQualsPtr[rd_pos+j] >= mapqual_ ){
                            cigar_base[j] = seq_nt16_str[bam_seqi(htsSeqPtr, (rd_pos+j))];
                        }else{
                            cigar_base[j] = tolower(seq_nt16_str[bam_seqi(htsSeqPtr, (rd_pos+j))]);
                        }
                    }
                    seqInfo.push_back(cigar_base);
                    if( VariantType::INS == type_ && (refer_start_offset+ref_pos)==flank_ && alt_ == cigar_base ){
                        cover_.add_coverage(refer_start_offset+ref_pos, VAR_VR);
                    }else{
                        cover_.add_coverage(refer_start_offset+ref_pos, VAR_RR);
                    }
                    if(VariantType::INDEL == type_ && refer_start_offset+ref_pos == flank_ ) coverage_ins = cigar_base;
                    break;
                case BAM_CDEL: // D
                    operations.push_back(Operation(OperationType::kDeletionFromRef, (uint32_t)c_len));
                    for(int j=0;j<c_len;++j){
                        cigar_base = ref_seq[start+ref_pos+j];
                    }
                    if( VariantType::DEL == type_ && (refer_start_offset+ref_pos)==flank_ && alt_ == cigar_base ){
                        cover_.add_coverage(refer_start_offset+ref_pos, VAR_VR);
                    }else{
                        cover_.add_coverage(refer_start_offset+ref_pos, VAR_RR);
                    }
                    if(VariantType::INDEL == type_ && refer_start_offset+ref_pos == flank_ ) coverage_del = cigar_base;
                    seqInfo.push_back(cigar_base);
                    break;
                case BAM_CREF_SKIP: // N
                    operations.push_back(Operation(OperationType::kMissingBases, (uint32_t)c_len));
                    for(int j=0;j<c_len;++j){
                        cigar_base = ref_seq[start+ref_pos+j];
                    }
                    seqInfo.push_back(cigar_base);
                    break;
                case BAM_CSOFT_CLIP: // S
                    // CNV not consider SOFT_CLIP
                    if( type_ == VariantType::LOSS || type_ == VariantType::DUP || type_ == VariantType::LOSSDUP ) break;
                    operations.push_back(Operation(OperationType::kSoftclip, (uint32_t)c_len));
                    for(int j=0;j<c_len;++j){
                        cigar_base[j] = seq_nt16_str[bam_seqi(htsSeqPtr, (rd_pos+j))];
                    }
                    seqInfo.push_back(cigar_base);
                    break;
                /** case BAM_CHARD_CLIP:
                *    cigar_ = "H";
                *    break;
                * case BAM_CPAD:
                *    cigar_ = "P";
                *    break;
                **/
                case BAM_CEQUAL: // =
                    operations.push_back(Operation(OperationType::kMatch, (uint32_t)c_len));
                    seqInfo.push_back(cigar_base);
                    cover_.add_coverage_range(refer_start_offset+ref_pos, VAR_RR, c_len);
                    break;
                case BAM_CDIFF: // X
                    operations.push_back(Operation(OperationType::kMismatch, (uint32_t)c_len));
                    for(int j=0;j<c_len;++j){
                        if(htsQualsPtr[rd_pos+j] >= mapqual_ ){
                            cigar_base[j] = seq_nt16_str[bam_seqi(htsSeqPtr, (rd_pos+j))];
                            cover_.add_coverage(refer_start_offset+ref_pos+j, VAR_VR);
                        }else{
                            cigar_base[j] = tolower(seq_nt16_str[bam_seqi(htsSeqPtr, (rd_pos+j))]);
                        }
                    }
                    cover_.add_coverage_range(refer_start_offset+ref_pos, VAR_VR, c_len);
                    seqInfo.push_back(cigar_base);
                    break;
                default:
                    break;
            }
            if ((bam_cigar_type(cigar[i]) & 0x01) != 0 ){
                rd_pos += c_len;
            }
            if((bam_cigar_type(cigar[i]) & 0x02) !=0 ){
                ref_pos += c_len;
            }
        }
        //just for insdel
        if( type_ == VariantType::INDEL && (alt_== (coverage_del+coverage_ins) || alt_== (coverage_ins+coverage_del)) ){
            cover_.add_coverage(flank_, VAR_VR);
        }
        Alignment this_align(refer_start_offset, operations);
        end = start + this_align.referenceLength() - 1;
        if( start > end_ || end < start_ ) continue;
        // Store read and alignment
        Read this_read(seqInfo, rec->core.qual, strand, 1, this_align);

        /**
        * std::cerr<<"Read:"<<read_name<<",start:"<<std::to_string(start)<<",end:"<<std::to_string(end)<<std::endl;
        * std::cerr<<"offset:"<<std::to_string(refer_start_offset)<<",Alignment:"<<this_align<<std::endl;
        * for(auto j=seqInfo.begin();j != seqInfo.end(); ++j){
        *    std::cerr<<*j;
        * }
        * std::cerr<<std::endl;
        **/
        if (pairedCache.find(read_name) == pairedCache.end()) {
            Read mate_read(0, true, 0);
            Frag tmp_Frag(this_read, mate_read);
            pairedCache.emplace(read_name, tmp_Frag);
        } else {
            pairedCache.at(read_name).mate = this_read;
        }
    }
    bam_destroy1(rec);
    hts_itr_destroy(iter);
    for (auto i = pairedCache.begin(); i != pairedCache.end(); ++i) {
        //std::cerr<< i->first << " " << i->second.read <<i->second.mate << endl;
        check_dup(i->second.read,i->second.mate);
    }
    if( AlignReads.empty()){
        throw std::logic_error("No read mapped to this region.");
    }
}

//add_coverage(start_,end_, start+ref_pos, (start+ref_pos+c_len-1), (refer_start_offset+ref_pos));
void Aligns::add_coverage(int start, int end, int ss, int ee, int idx){
    if( ss <= end  && ee >= start ){
        int min = (ss < start ) ? start : ss;
        int max = (ee < end ) ?  ee : end;
        int len = max - min+1;
        if( len > 1 ){
            cover_.add_coverage_range(idx, VAR_RR, len);
        }else{
            cover_.add_coverage(idx, VAR_RR);
        }
    }
}
void Aligns::AlignCNV() {
    const string chr_ = var_.genomic_pos.front().chrName;
    const int chr_idx = sam_hdr_name2tid(htsHeaderPtr,chr_.c_str());
    int32_t chr_len = (int32_t)htsHeaderPtr->target_len[chr_idx];
    string ref_seq = seq_.getSequence(chr_);
    std::transform(ref_seq.begin(), ref_seq.end(), ref_seq.begin(), ::toupper);
    std::map<string, Frag>pairedCache;

    const int64_t start_ = var_.genomic_pos.front().start_;
    const int64_t end_ = var_.genomic_pos.back().end_;
    const int32_t tmp_end = ( (end_ + flank_) > chr_len ) ? chr_len : (end_+flank_);
    const int32_t tmp_start = (start_ > flank_) ? start_-flank_ : 1;
    hts_itr_t *iter = sam_itr_queryi(htsIndexPtr, chr_idx, tmp_start, tmp_end);
    if (iter == nullptr) {
        string region= chr_ + ":"+ to_string(start_) + "-" + to_string(end_);
        throw logic_error("Failed to extract reads from the specified region:"+region);
    }
    bam1_t *rec = bam_init1();
    while((sam_itr_next(htsFilePtr,iter,rec)) >=0 ) {
        if( rec->core.flag & BAM_FUNMAP ) continue;
        if( rec->core.flag & filter_ ) continue;
        if( rec->core.qual < (uint8_t) mapqual_ ) continue;

        bool strand;
        string bases,quals,read_name,this_cigar;
        int start,end,readLength,refer_start_offset,c_len,c_op,rd_pos, ref_pos;
        uint8_t *htsSeqPtr,*htsQualsPtr;
        uint32_t *cigar;
        std::list<string> seqInfo;
        std::list<Operation> operations;

        start = rec->core.pos;
        refer_start_offset = start - refer_start_;
        if ( refer_start_offset < 0 ) continue;
        read_name = bam_get_qname(rec);
        strand = (bam_is_rev(rec) == 0);
        htsSeqPtr = bam_get_seq(rec);// read sequence
        htsQualsPtr = bam_get_qual(rec);// read quality
        cigar = bam_get_cigar(rec); // CIGAR
        readLength = rec->core.l_qseq;
        ref_pos=rd_pos=0;
        for (int i = 0; i < rec->core.n_cigar; ++i) {
            c_op = bam_cigar_op(cigar[i]);
            c_len = bam_cigar_oplen(cigar[i]);
            string cigar_base(c_len,' ');
            switch (c_op) {
                case BAM_CMATCH: // M
                    for(int j=0;j<c_len;++j){
                        if( seq_nt16_str[bam_seqi(htsSeqPtr, (rd_pos+j))] != ref_seq[start+ref_pos+j]){
                            if(htsQualsPtr[rd_pos+j] >= mapqual_ ){
                                cigar_base[j] = seq_nt16_str[bam_seqi(htsSeqPtr, (rd_pos+j))];
                            }else{
                                cigar_base[j] = tolower(seq_nt16_str[bam_seqi(htsSeqPtr, (rd_pos+j))]);
                            }
                        }
                    }
                    add_coverage(start_,end_, start+ref_pos, (start+ref_pos+c_len-1), (refer_start_offset+ref_pos));
                    operations.push_back(Operation(OperationType::kMatch, (uint32_t)c_len));
                    seqInfo.push_back(cigar_base);
                    break;
                case BAM_CINS: // I
                    // CNV not consider Insertion
                    operations.push_back(Operation(OperationType::kInsertionToRef, (uint32_t)c_len));
                    //for(int j=0;j<c_len;++j){
                    //    if(htsQualsPtr[rd_pos+j] >= mapqual_ ){
                    //        cigar_base[j] = seq_nt16_str[bam_seqi(htsSeqPtr, (rd_pos+j))];
                    //    }else{
                    //        cigar_base[j] = tolower(seq_nt16_str[bam_seqi(htsSeqPtr, (rd_pos+j))]);
                    //    }
                    //}
                    seqInfo.push_back(cigar_base);
                    break;
                case BAM_CDEL: // D
                    operations.push_back(Operation(OperationType::kDeletionFromRef, (uint32_t)c_len));
                    for(int j=0;j<c_len;++j){
                        cigar_base[j] = ref_seq[start+ref_pos+j];
                    }
                    seqInfo.push_back(cigar_base);
                    break;
                case BAM_CREF_SKIP: // N
                    operations.push_back(Operation(OperationType::kMissingBases, (uint32_t)c_len));
                    //for(int j=0;j<c_len;++j){
                    //    cigar_base = ref_seq[start+ref_pos+j];
                    //}
                    seqInfo.push_back(cigar_base);
                    break;
                case BAM_CSOFT_CLIP: // S
                    // CNV not consider SOFT_CLIP
                    operations.push_back(Operation(OperationType::kSoftclip, (uint32_t)c_len));
                    seqInfo.push_back(cigar_base);
                    break;
                case BAM_CEQUAL: // =
                    operations.push_back(Operation(OperationType::kMatch, (uint32_t)c_len));
                    seqInfo.push_back(cigar_base);
                    add_coverage(start_,end_, start+ref_pos, (start+ref_pos+c_len-1), (refer_start_offset+ref_pos));
                    break;
                case BAM_CDIFF: // X
                    operations.push_back(Operation(OperationType::kMismatch, (uint32_t)c_len));
                    for(int j=0;j<c_len;++j){
                        if(htsQualsPtr[rd_pos+j] >= mapqual_ ){
                            cigar_base[j] = seq_nt16_str[bam_seqi(htsSeqPtr, (rd_pos+j))];
                        }else{
                            cigar_base[j] = tolower(seq_nt16_str[bam_seqi(htsSeqPtr, (rd_pos+j))]);
                        }
                    }
                    add_coverage(start_,end_, start+ref_pos, (start+ref_pos+c_len-1), (refer_start_offset+ref_pos));
                    seqInfo.push_back(cigar_base);
                    break;
                default:
                    break;
            }
            if ((bam_cigar_type(cigar[i]) & 0x01) != 0 ){
                rd_pos += c_len;
            }
            if((bam_cigar_type(cigar[i]) & 0x02) !=0 ){
                ref_pos += c_len;
            }
        }
        Alignment this_align(refer_start_offset, operations);
        end = start + this_align.referenceLength() - 1;
        if( start > end_ || end < start_ ) continue;
        // Store read and alignment
        Read this_read(seqInfo, rec->core.qual, strand, 1, this_align);

        if (pairedCache.find(read_name) == pairedCache.end()) {
            Read mate_read(0, true, 0);
            Frag tmp_Frag(this_read, mate_read);
            pairedCache.emplace(read_name, tmp_Frag);
        } else {
            pairedCache.at(read_name).mate = this_read;
        }
    }
    bam_destroy1(rec);
    hts_itr_destroy(iter);

    for (auto i = pairedCache.begin(); i != pairedCache.end(); ++i) {
        //std::cerr<< i->first << " " << i->second.read <<i->second.mate << endl;
        check_dup(i->second.read,i->second.mate);
    }
    if( AlignReads.empty()){
        throw std::logic_error("No read mapped to this region.");
    }
}
void Aligns::AlignsFusion() {
    const string first_chr = var_.genomic_pos.front().chrName;
    const int64_t first_bp = var_.genomic_pos.front().start_;
    const string second_chr = var_.genomic_pos.back().chrName;
    const int64_t second_bp = var_.genomic_pos.back().start_;

    string ref_seq = seq_.getSequence(first_chr);
    std::transform(ref_seq.begin(), ref_seq.end(), ref_seq.begin(), ::toupper);
    string first_left_seq=ref_seq.substr(first_bp-flank_, flank_);
    string rev_first_left = reverseComplement(first_left_seq);
    string first_right_seq=ref_seq.substr(first_bp, flank_);
    string rev_first_right = reverseComplement(first_right_seq);
    string first_all = ref_seq.substr(first_bp-flank_, 2*flank_);
    //std::cerr<<"first_left:"<<first_left_seq<<std::endl;
    //std::cerr<<"rev_first_left:"<<rev_first_left<<std::endl;
    //std::cerr<<"first_right:"<<first_right_seq<<std::endl;
    //std::cerr<<"rev_first_right:"<<rev_first_right<<std::endl;
    //std::cerr<<"first_all:"<<first_all<<std::endl;
    if( second_chr != first_chr) ref_seq=seq_.getSequence(second_chr);
    string second_left_seq = ref_seq.substr(second_bp-flank_, flank_);
    string rev_second_left = reverseComplement(second_left_seq);
    string second_right_seq = ref_seq.substr(second_bp, flank_);
    string rev_second_right = reverseComplement(second_right_seq);
    string second_all = ref_seq.substr(second_bp-flank_, 2*flank_);
    //std::cerr<<"second_left:"<<second_left_seq<<std::endl;
    //std::cerr<<"rev_second_left:"<<rev_second_left<<std::endl;
    //std::cerr<<"second_right:"<<second_right_seq<<std::endl;
    //std::cerr<<"rev_second_right:"<<rev_second_right<<std::endl;
    //std::cerr<<"second_all:"<<second_all<<std::endl;

    //record the FL-SR, SL-FR, FL-revSL, revSR-FR, revFR-SR, SL-revFL read count;
    // revFR-SR equal to revSR-FR, SL-revFL equal to FL-revSL
    // 0-index for softclip reads
    // 1-index for PE reads
    int direct[2][6]={{0,0,0,0},{0,0,0,0}};
    std::map<string, Frag>pairedCache;
    std::map<string, Frag>fusionCache;

    // First break point
    int left_chr_idx = sam_hdr_name2tid(htsHeaderPtr,first_chr.c_str());
    int right_chr_idx = sam_hdr_name2tid(htsHeaderPtr, second_chr.c_str());
    hts_itr_t *iter = sam_itr_queryi(htsIndexPtr, left_chr_idx, first_bp-flank_-1, first_bp+flank_);
    if (iter == nullptr) {
        string region= first_chr + ":"+ to_string(first_bp-flank_) + "-" + to_string(first_bp+flank_);
        throw logic_error("Failed to extract reads from the specified region:"+region);
    }

    //std::cerr<<"Read first bam: chr:"<<first_chr<<",breakpoint:"<<first_bp<<std::endl;
    bam1_t *rec = bam_init1();
    ref_seq = seq_.getSequence(first_chr);
    while((sam_itr_next(htsFilePtr,iter,rec)) >=0 ) {
        if( rec->core.flag & BAM_FUNMAP ) continue;
        if( rec->core.flag & BAM_FDUP ) continue;
        if( rec->core.flag & filter_ ) continue;
        //if (bam_aux2i(bam_aux_get(b,"NM")) >3) continue;

        uint32_t *cigar, *end_cigar;
        int start, ref_pos, rd_pos,c_op, c_len;
        uint8_t *htsSeqPtr, *htsQualsPtr;
        string read_seq,read_name;
        bool strand,softclip_read,fusion_pe_read,fusion_read;
        std::list<string> seqInfo;
        std::list<Operation> operations;

        htsSeqPtr = bam_get_seq(rec);// read sequence
        htsQualsPtr = bam_get_qual(rec);// read quality
        //for (int i = 0; i < rec->core.l_qseq; i++) {
	    //    read_seq += seq_nt16_str[bam_seqi(htsSeqPtr, i)];
        //}
        //std::cerr<<"read_name:"<<bam_get_qname(rec)<<",seq:"<<read_seq<<std::endl;

        start = rec->core.pos; // 0-based
        read_name = bam_get_qname(rec);
        strand = (bam_is_rev(rec)==0);
        cigar = bam_get_cigar(rec); // CIGAR
        end_cigar = cigar+rec->core.n_cigar-1;
        softclip_read = false;
        fusion_pe_read = false;
        fusion_read = false;
        end_cigar = cigar+rec->core.n_cigar-1;
        if(bam_cigar_op(*cigar) == BAM_CSOFT_CLIP || bam_cigar_op(*end_cigar) == BAM_CSOFT_CLIP ){
            softclip_read = true;
        }else if(rec->core.mtid == right_chr_idx && rec->core.qual >= (uint8_t) mapqual_ ){
            if( rec->core.mpos >= (second_bp - flank_) && rec->core.mpos <= (second_bp + flank_) ){
                fusion_pe_read = true;
            }
        }
        if( !softclip_read && !fusion_pe_read ) continue;
        //std::cerr<<"read name:"<<read_name<<",pos:"<<std::to_string(start)<<std::endl;
        ref_pos=rd_pos=0;
        for(int i=0; i< rec->core.n_cigar;++i){
            c_op = bam_cigar_op(cigar[i]);
            c_len = bam_cigar_oplen(cigar[i]);
            string cigar_base(c_len, ' ');
            switch(c_op){
                case BAM_CMATCH: // M
                    for(int j=0;j<c_len;++j){
                        if( seq_nt16_str[bam_seqi(htsSeqPtr, (rd_pos+j))]!= ref_seq[start+ref_pos+j]){
                            if(htsQualsPtr[rd_pos+j] >= mapqual_ ){
                                cigar_base[j] = seq_nt16_str[bam_seqi(htsSeqPtr, (rd_pos+j))];
                            }else{
                                cigar_base[j] = tolower(seq_nt16_str[bam_seqi(htsSeqPtr, (rd_pos+j))]);
                            }
                        }
                    }
                    operations.push_back(Operation(OperationType::kMatch, (uint32_t)c_len));
                    seqInfo.push_back(cigar_base);
                    break;
                case BAM_CINS: // I
                    operations.push_back(Operation(OperationType::kInsertionToRef, (uint32_t)c_len));
                    for(int j=0;j<c_len;++j){
                        if(htsQualsPtr[rd_pos+j] >= mapqual_ ){
                            cigar_base[j] = seq_nt16_str[bam_seqi(htsSeqPtr, (rd_pos+j))];
                        }else{
                            cigar_base[j] = tolower(seq_nt16_str[bam_seqi(htsSeqPtr, (rd_pos+j))]);
                        }
                    }
                    seqInfo.push_back(cigar_base);
                    break;
                case BAM_CDEL: // D
                    operations.push_back(Operation(OperationType::kDeletionFromRef, (uint32_t)c_len));
                    for(int j=0;j<c_len;++j){
                        cigar_base = ref_seq[start+ref_pos+j];
                    }
                    seqInfo.push_back(cigar_base);
                    break;
                case BAM_CREF_SKIP: // N
                    operations.push_back(Operation(OperationType::kMissingBases, (uint32_t)c_len));
                    for(int j=0;j<c_len;++j){
                        cigar_base = ref_seq[start+ref_pos+j];
                    }
                    seqInfo.push_back(cigar_base);
                    break;
                case BAM_CSOFT_CLIP: // S
                    operations.push_back(Operation(OperationType::kSoftclip, (uint32_t)c_len));
                    for(int j=0;j<c_len;++j){
                        cigar_base[j] = seq_nt16_str[bam_seqi(htsSeqPtr, (rd_pos+j))];
                    }
                    seqInfo.push_back(cigar_base);
                    if( i == 0 && (start+ref_pos) == first_bp ){
                        //S-first --> SL-FR
                        if( softclip_match(cigar_base, second_left_seq, c_len , 1, true) ){
                            if( c_len > 10 ) direct[0][1]++;
                            fusion_read = true;
                        //S-first ---> revSR-FR
                        }else if( softclip_match(cigar_base, rev_second_right, c_len, 1, true)){
                            if(c_len > 10 ) direct[0][3]++;
                            fusion_read = true;
                        }
                    }else if(i == (rec->core.n_cigar - 1) && (start+ref_pos) == first_bp ){
    		            // first-S ---> FL-SR
		                if( softclip_match(cigar_base, second_right_seq, c_len, 1) ){
		                    if(c_len > 10 ) direct[0][0]++;
		                    fusion_read = true;
   		                //first-S ---> FL-revSL
		                }else if( softclip_match(cigar_base, rev_second_left, c_len, 1)){
		                    if( c_len > 10 ) direct[0][2]++;
		                    fusion_read = true;
		                }
		            }
		            break;
                case BAM_CEQUAL: // =
                    operations.push_back(Operation(OperationType::kMatch, (uint32_t)c_len));
                    seqInfo.push_back(cigar_base);
                    break;
                case BAM_CDIFF: // X
                    operations.push_back(Operation(OperationType::kMismatch, (uint32_t)c_len));
                    for(int j=0;j<c_len;++j){
                        if(htsQualsPtr[rd_pos+j] >= mapqual_ ){
                            cigar_base[j] = seq_nt16_str[bam_seqi(htsSeqPtr, (rd_pos+j))];
                        }else{
                            cigar_base[j] = tolower(seq_nt16_str[bam_seqi(htsSeqPtr, (rd_pos+j))]);
                        }
                    }
                    seqInfo.push_back(cigar_base);
                    break;
                default:
                    break;
            }
            if ((bam_cigar_type(cigar[i]) & 0x01) != 0 && c_op != BAM_CBACK ){
                rd_pos += c_len;
            }
            if((bam_cigar_type(cigar[i]) & 0x02) !=0 && c_op != BAM_CBACK){
                ref_pos += c_len;
            }
        }
        Alignment this_align(start, operations);
        Read this_read(seqInfo, rec->core.qual, strand, 1, this_align);
        Read mate_read(0, true, 0);
        Frag tmp_Frag(this_read, mate_read);
        if( fusion_read ){
            fusionCache.emplace(read_name, tmp_Frag);
        }else if( fusion_pe_read ){
            pairedCache.emplace(read_name, tmp_Frag);
        }
    }

    //std::cerr<<"Read second bam..."<<std::endl;
    //Second break point
    iter = sam_itr_queryi(htsIndexPtr, right_chr_idx, second_bp-flank_-1, second_bp+flank_);
    if (iter == nullptr) {
        string region= second_chr + ":"+ to_string(second_bp-flank_) + "-" + to_string(second_bp+flank_);
        throw logic_error("Failed to extract reads from the specified region:"+region);
    }
    ref_seq = seq_.getSequence(second_chr);
    while((sam_itr_next(htsFilePtr,iter,rec)) >=0 ) {
        if( rec->core.flag & BAM_FUNMAP ) continue;
        if( rec->core.flag & BAM_FDUP ) continue;
        if( rec->core.flag & filter_ ) continue;
        //if (bam_aux2i(bam_aux_get(b,"NM")) >3) continue;

        uint32_t *cigar, *end_cigar;
        int start, ref_pos, rd_pos,c_op,c_len;
        uint8_t *htsSeqPtr, *htsQualsPtr;
        string read_seq,read_name;
        bool strand,fusion_read,fusion_pe_read, softclip_read;
        std::list<string> seqInfo;
        std::list<Operation> operations;

        htsSeqPtr = bam_get_seq(rec);// read sequence
        htsQualsPtr = bam_get_qual(rec);// read quality
        //for (int i = 0; i < rec->core.l_qseq; i++) {
	    //    read_seq += seq_nt16_str[bam_seqi(htsSeqPtr, i)];
        //}
        //std::cerr<<"read_name:"<<bam_get_qname(rec)<<",seq:"<<read_seq<<std::endl;
        start = rec->core.pos;
        read_name = bam_get_qname(rec);
        strand = (bam_is_rev(rec)==0);
        cigar = bam_get_cigar(rec); // CIGAR
        end_cigar = cigar+rec->core.n_cigar-1;
        softclip_read = false;
        fusion_read = false;
        fusion_pe_read = false;
        end_cigar = cigar+rec->core.n_cigar-1;
        if( fusionCache.find(read_name) != fusionCache.end() || bam_cigar_op(*cigar) == BAM_CSOFT_CLIP || bam_cigar_op(*end_cigar) == BAM_CSOFT_CLIP){
            softclip_read = true;
        }else if(rec->core.mtid == left_chr_idx && rec->core.qual >= (uint8_t)mapqual_){
            if(pairedCache.find(read_name) != pairedCache.end()){
                fusion_pe_read = true;
            }
        }
        if( !softclip_read && !fusion_pe_read ) continue;
        ref_pos=rd_pos=0;
        for(int i=0; i< rec->core.n_cigar;++i){
            c_op = bam_cigar_op(cigar[i]);
            c_len = bam_cigar_oplen(cigar[i]);
            string cigar_base(c_len, ' ');
            switch(c_op){
                case BAM_CMATCH: // M
                    for(int j=0;j<c_len;++j){
                        if( seq_nt16_str[bam_seqi(htsSeqPtr, (rd_pos+j))] !=  ref_seq[start+ref_pos+j]){
                            if(htsQualsPtr[rd_pos+j] >= mapqual_ ){
                                cigar_base[j] = seq_nt16_str[bam_seqi(htsSeqPtr, (rd_pos+j))];
                            }else{
                                cigar_base[j] = tolower(seq_nt16_str[bam_seqi(htsSeqPtr, (rd_pos+j))]);
                            }
                        }
                    }
                    operations.push_back(Operation(OperationType::kMatch, (uint32_t)c_len));
                    seqInfo.push_back(cigar_base);
                    break;
                case BAM_CINS: // I
                    operations.push_back(Operation(OperationType::kInsertionToRef, (uint32_t)c_len));
                    for(int j=0;j<c_len;++j){
                        if(htsQualsPtr[rd_pos+j] >= mapqual_ ){
                            cigar_base[j] = seq_nt16_str[bam_seqi(htsSeqPtr, (rd_pos+j))];
                        }else{
                            cigar_base[j] = tolower(seq_nt16_str[bam_seqi(htsSeqPtr, (rd_pos+j))]);
                        }
                    }
                    seqInfo.push_back(cigar_base);
                    break;
                case BAM_CDEL: // D
                    operations.push_back(Operation(OperationType::kDeletionFromRef, (uint32_t)c_len));
                    for(int j=0;j<c_len;++j){
                        cigar_base = ref_seq[start+ref_pos+j];
                    }
                    seqInfo.push_back(cigar_base);
                    break;
                case BAM_CREF_SKIP: // N
                    operations.push_back(Operation(OperationType::kMissingBases, (uint32_t)c_len));
                    for(int j=0;j<c_len;++j){
                        cigar_base = ref_seq[start+ref_pos+j];
                    }
                    seqInfo.push_back(cigar_base);
                    break;
                case BAM_CSOFT_CLIP: // S
                    operations.push_back(Operation(OperationType::kSoftclip, (uint32_t)c_len));
                    for(int j=0;j<c_len;++j){
                        cigar_base[j] = seq_nt16_str[bam_seqi(htsSeqPtr, (rd_pos+j))];
                    }
                    seqInfo.push_back(cigar_base);
                    if( i == 0 && start == second_bp ){
                        //S-second --> FL-SR
                        if( softclip_match(cigar_base, first_left_seq, c_len , 1, true) ){
                            if ( c_len > 10 ) direct[0][0]++;
                                fusion_read = true;
                        //S-second ---> revFR-SR, equal to revSR-FR
                        }else if( softclip_match(cigar_base, rev_first_right, c_len, 1, true)){
                            if( c_len > 10 ) direct[0][3]++;
                            fusion_read = true;
                        }
                        //std::cerr<<"second Left: "<<cigar_base<<",read_name:"<<bam_get_qname(rec)<<std::endl;
                    } else if( i == (rec->core.n_cigar - 1) ){
		                if( (start+ref_pos) == second_bp ) {
		                    // second-S ---> SL-FR
		                    if( softclip_match(cigar_base, first_right_seq, c_len, 1) ){
		                        if( c_len > 10 ) direct[0][1]++;
                                fusion_read = true;
		                    // second-S ---> SL-revFL, equal to FL-revSL
		                    }else if( softclip_match(cigar_base, rev_first_left, c_len, 1)){
		                        if( c_len > 10 ) direct[0][2]++;
		                        fusion_read = true;
		                    }
		                    //std::cerr<<"second Right: "<<cigar_base<<",read_name:"<<bam_get_qname(rec)<<std::endl;
                        }
                    }
                    break;
                case BAM_CEQUAL: // =
                    operations.push_back(Operation(OperationType::kMatch, (uint32_t)c_len));
                    seqInfo.push_back(cigar_base);
                    break;
                case BAM_CDIFF: // X
                    operations.push_back(Operation(OperationType::kMismatch, (uint32_t)c_len));
                    for(int j=0;j<c_len;++j){
                        if(htsQualsPtr[rd_pos+j] >= mapqual_ ){
                            cigar_base[j] = seq_nt16_str[bam_seqi(htsSeqPtr, (rd_pos+j))];
                        }else{
                            cigar_base[j] = tolower(seq_nt16_str[bam_seqi(htsSeqPtr, (rd_pos+j))]);
                        }
                    }
                    seqInfo.push_back(cigar_base);
                    break;
                default:
                    break;
            }
            if ((bam_cigar_type(cigar[i]) & 0x01) != 0 && c_op != BAM_CBACK ){
                rd_pos += c_len;
            }
            if((bam_cigar_type(cigar[i]) & 0x02) !=0 && c_op != BAM_CBACK){
                ref_pos += c_len;
            }
        }
        Alignment this_align(start, operations);
        Read this_read(seqInfo, rec->core.qual, strand, 1, this_align);
        if( fusion_read ){
            if( fusionCache.find(read_name) != fusionCache.end()){
                fusionCache.at(read_name).mate = this_read;
            }else if( pairedCache.find(read_name) != pairedCache.end()){
                Read mate_read = pairedCache.at(read_name).read;
                Frag tmp_Frag(mate_read,this_read);
                fusionCache.emplace(read_name, tmp_Frag);
            }else{
                Read mate_read(0, true, 0);
                Frag tmp_Frag(mate_read, this_read);
                fusionCache.emplace(read_name, tmp_Frag);
            }
        }else if( fusion_pe_read ){
            if( fusionCache.find(read_name) != fusionCache.end()){
                fusionCache.at(read_name).mate = this_read;
            }else{
                if (pairedCache.find(read_name) != pairedCache.end()) {
                    pairedCache.at(read_name).mate = this_read;
                }
            }
        }
    }
    bam_destroy1(rec);
    hts_itr_destroy(iter);
    ref_seq = "";

    // filtered paired reads
    for (auto i = pairedCache.begin(); i != pairedCache.end(); ++i) {
        if( i->second.mate.count > 0 ){
            int first_read_start = i->second.read.align.referenceStart();
            int first_read_end = first_read_start + i->second.read.align.referenceLength();
            int second_read_start = i->second.mate.align.referenceStart();
            int second_read_end = second_read_end + i->second.mate.align.referenceLength();
            bool add = false;
            if( first_read_start > first_bp-flank_ && first_read_end < first_bp ) {
                if( second_read_start > second_bp-flank_ && second_read_end < second_bp){
                    direct[1][2]++;
                    add = true;
                }else if( second_read_start > second_bp && second_read_end < second_bp + flank_ ){
                    direct[1][0]++;
                    add = true;
                }
            }else if( first_read_start > first_bp && first_read_end < first_bp + flank_ ){
                if( second_read_start > second_bp-flank_ && second_read_end <= second_bp ){
                    direct[1][1]++;
                    add = true;
                }else if( second_read_start >= second_bp && second_read_end < second_bp + flank_ ){
                    direct[1][3]++;
                    add = true;
                }
            }
            if ( add ){
                Frag tmp_Frag(pairedCache.at(i->first).read, pairedCache.at(i->first).mate);
                fusionCache.emplace(i->first, tmp_Frag);
            }
        }
    }
	/**
	* samFile *bam_out;
	* std::string tmp_str="fusion_softclip_pe_read.sam";
	* char* outfile=(char*)calloc(1024,sizeof(char));
	* sprintf(outfile, "%s", tmp_str.c_str());
	* bam_out = sam_open(outfile, "w");
	* sam_close(bam_out);
	*/
	//std::cerr<<"Output for reads:"<<std::endl;
    //for (auto i = pairedCache.begin(); i != pairedCache.end(); ++i) {
        //std::cerr<< "first reads:"<<i->first << ",pos:"<<std::to_string(i->second.read.align.referenceStart())<<","<<i->second.read<<std::endl;
        //std::cerr<< "second reads:"<<i->first<<",pos:"<<std::to_string(i->second.mate.align.referenceStart())<<","<<i->second.mate<< std::endl;
    //}
    string fusion_direct[4]={"First_Left->Second_Right","Second_Left->First_Right","First_Left->revSecond_Left","revSecond_Right->First_Right"};
    for(int i=0;i<4;++i){
        //if(direct[0][i] == direct[0][idx_sc]) count_sc++;
        std::cerr<<"Index: "<<std::to_string(i)<<"\t"<<std::to_string(direct[0][i])<<"\t"<<std::to_string(direct[1][i])<<"\t"<<fusion_direct[i]<<std::endl;
    }
    int idx = 0;
    int count = 0;
    for(int i=0;i<4;++i){
        int total = (direct[0][i] + direct[1][i]);
        if( total > (direct[0][idx]+direct[1][idx]))  idx = i;
        if( total > 5) count++;
    }
    if( count > 1 ){
        idx = 0;
        for(int i=0;i<4;++i){
            if( direct[0][i] >0 && direct[1][i] >0){
                int total = (direct[0][i]+direct[1][i]);
                if(total > (direct[0][idx]+direct[1][idx])) idx = i;
            }
        }
        //string region= first_chr + ":"+ to_string(first_bp) + "-" + second_chr + ":" + to_string(second_bp);
        //throw logic_error("This fusion direction is not uniq:"+region+"\n");
    }
    //if( count > 1){
    //    throw logic_error("This fusion direction is not uniq:"+region+"\n");
    //}
    //if( direct[idx] < 3 ){
    //    std::cerr<<"The support softclip reads is less 3. in fusion: "<<region<<std::endl;
    //}
    pairedCache.clear();
    fusion_type_ = FusionType(idx);
    if( fusion_type_ == FusionType::FLSR){// first_left + second_right
            refer_start_ = first_bp - flank_;
            int right_refer_start_ = second_bp;
        for (auto i = fusionCache.begin(); i != fusionCache.end(); ++i) {
            if( i->second.read.count > 0 ){
                int refer_start_offset = i->second.read.align.referenceStart() - refer_start_;
                if( refer_start_offset < 0 ) std::cerr<<"FLSR-lllllllllllllllllll\n";
                i->second.read.align.set_reference_start(refer_start_offset);
            }
            if(i->second.mate.count > 0){
                int refer_start_offset = i->second.mate.align.referenceStart() - right_refer_start_;
                if( refer_start_offset < 0 ){
                    std::cerr<<"FLSR-rrrrrrrrrrrrrrrrrrrr\n";
                }
                i->second.mate.align.set_reference_start(refer_start_offset+flank_);
            }
            check_dup(i->second.read,i->second.mate, false);
        }
    }else if( fusion_type_ == FusionType::SLFR){// second_left + first_right
        refer_start_ = second_bp - flank_;
        int right_refer_start_ = first_bp;
        for (auto i = fusionCache.begin(); i != fusionCache.end(); ++i) {
            if( i->second.mate.count > 0 ){
                int refer_start_offset = i->second.mate.align.referenceStart() - refer_start_;
                if( refer_start_offset < 0 ) std::cerr<<"SLFR-lllllllllllllllllll\n";
                i->second.mate.align.set_reference_start(refer_start_offset);
            }
            if(i->second.read.count > 0){
                int refer_start_offset = i->second.read.align.referenceStart() - right_refer_start_;
                if( refer_start_offset < 0 ){
                    std::cerr<<"SLFR-rrrrrrrrrrrrrrr\n";
                }
                i->second.read.align.set_reference_start(refer_start_offset+flank_);
            }
            check_dup(i->second.mate,i->second.read, false);
        }
    }else if( fusion_type_ == FusionType::FLrevSL){// first_left + revComp second_left
        refer_start_ = first_bp - flank_;
        int right_refer_start_ = second_bp;
        for (auto i = fusionCache.begin(); i != fusionCache.end(); ++i) {
            if( i->second.read.count >0 && i->second.mate.count >0 ){
                int first_read_start = i->second.read.align.referenceStart();
                int first_read_end = first_read_end + i->second.read.align.referenceLength();
                int second_read_start = i->second.mate.align.referenceStart();
                int second_read_end = second_read_end + i->second.mate.align.referenceLength();
                if( first_read_start >= first_bp || second_read_start >= second_bp  ){
                    std::cerr<<" Paired-End read: "<<i->first<<" is support on different direct."<<std::endl;
                    continue;
                }
            }
            if( i->second.read.count > 0 ){
                int refer_start_offset = i->second.read.align.referenceStart() - refer_start_+1;
                if( refer_start_offset < 0 ){
                    string infor="FLrevSL:read_name:"+i->first+",left_read start:"+std::to_string(i->second.read.align.referenceStart())+",length:"+std::to_string(i->second.read.align.referenceLength())+"\n";
                    throw logic_error(infor);
                }
                i->second.read.align.set_reference_start(refer_start_offset);
            }
            if(i->second.mate.count > 0){
                //std::cerr<<"BeforeReverse, read_name:"<<i->first<<",right_read start:"<<std::to_string(i->second.mate.align.referenceStart())<<",length:"<<std::to_string(i->second.read.align.referenceLength())<<"\n";
                reverseRead(i->second.mate);
                //std::cerr<<"AfterReverse, read_name:"<<i->first<<",right_read start:"<<std::to_string(i->second.mate.align.referenceStart())<<",length:"<<std::to_string(i->second.read.align.referenceLength())<<"\n";
                int refer_start_offset = right_refer_start_ - i->second.mate.align.referenceStart() +1;
                if( refer_start_offset < 0 ){
                    string infor="FLrevSL:read_name:"+i->first+",right_read start:"+std::to_string(i->second.mate.align.referenceStart())+",length:"+std::to_string(i->second.read.align.referenceLength())+"\n";
                    throw logic_error(infor);
                }
                i->second.mate.align.set_reference_start(refer_start_offset+flank_);
            }
            if( i->second.read.count > 0 ){
                check_dup(i->second.read,i->second.mate, false);
            }else{
                check_dup(i->second.mate,i->second.read, false);
            }
        }
    }else{//revSRFR, revComp second_right + first_right
        refer_start_ = second_bp - flank_;
        int right_refer_start_ = first_bp;
        for (auto i = fusionCache.begin(); i != fusionCache.end(); ++i) {
            if( i->second.mate.count > 0 ){
                int refer_start_offset = i->second.mate.align.referenceStart() - refer_start_;
                if( refer_start_offset < 0 ) std::cerr<<"revSRFR-lllllllllllllllllll\n";
                i->second.mate.align.set_reference_start(refer_start_offset);
            }
            if(i->second.read.count > 0){
                int refer_start_offset = i->second.read.align.referenceStart() - right_refer_start_;
                if( refer_start_offset < 0 ){
                    std::cerr<<"revSRFR-rrrrrrrrrrrrrrr\n";
                }
                i->second.read.align.set_reference_start(refer_start_offset+flank_);
            }
            check_dup(i->second.mate,i->second.read, false);
        }
    }
}
std::ostream& operator<<(std::ostream& os, const Read& read) {
    std::list<Operation> operations_ = read.align.operations();
    std::list<Operation>::const_iterator operation_it = operations_.begin();
    std::list<string>::const_iterator base_it = read.bases.begin();
    while (operation_it != operations_.end()) {
        os<<"Operation: "<<*operation_it<<", Sequence:"<<*base_it<<",";
        ++operation_it;
        ++base_it;
    }
    os<<"Score:"<<read.score<<",Count:"<<read.count<<",Strand:"<<read.strand<<"\n";
    return os;
}
std::ostream& operator<<(std::ostream& os, Aligns& align) {
    os<<"refer_start_:"<<align.get_refer_start()<<",start_offset_:"<<align.get_start_offset()<<",end_offset_:"<<align.get_end_offset()<<"\n";
    os<<"Read infor:\n";
    for(int i=0;i<align.AlignReads.size();i++){
        os<<align.AlignReads[i].read;
        if( align.AlignReads[i].mate.count > 0) os<<align.AlignReads[i].mate;
    }
    return os;
}
}

