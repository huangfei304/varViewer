// varViewer
// Date: 2023-08-09
// Author: Huang Fei <huangfei@genomics.cn>
//

#pragma once

#include <iostream>
#include <memory>
#include <string>
#include <vector>
#include <boost/optional.hpp>
#include "utils.hh"

namespace varViewer
{
struct GenomicRegion {
    GenomicRegion(std::string chr, int start, int end)
        : chrName(chr)
        , start_(start)
        , end_(end)
    {
    }

    bool operator==(const GenomicRegion& other) const {
        return (other.chrName == chrName && other.start_ == start_ && other.end_ == end_);
    }

    std::string chrName;
    //int64_t start_;
    //int64_t end_;
    int start_;
    int end_;
};
enum class VariantType {
    SNP,
    INS, // insertion
    DEL, // deletion
    INDEL, // Insertion and deletion
    CNV,
    LOSS, // CNV Deletion
    DUP, // CNV Duplication
    LOSSDUP, //CNV Deletion and Duplicaiton in transcript
    SV
};

struct VariantSpecific {
    VariantSpecific(VariantType vartype, std::vector<GenomicRegion> gpos, std::string ref_, std::string allele_)
    : var_type(vartype)
    , genomic_pos(gpos)
    , ref(ref_)
    , allele(allele_)
    {
            infor=genomic_pos.front().chrName+"_"+std::to_string(genomic_pos.front().start_)+"_"+std::to_string(genomic_pos.back().end_);
            exon_size = 0;
            for(int i=0;i<genomic_pos.size();++i){
                exon_size += (genomic_pos[i].end_ - genomic_pos[i].start_ + 1);
            }
    }
    VariantSpecific(VariantType vartype, std::vector<GenomicRegion> gpos, std::string ref_, std::string allele_,std::string strand_)
    : var_type(vartype)
    , genomic_pos(gpos)
    , ref(ref_)
    , allele(allele_)
    ,strand(strand_)
    {
        infor=genomic_pos.front().chrName+"_"+std::to_string(genomic_pos.front().start_)+"_"+std::to_string(genomic_pos.back().end_);
        exon_size = 0;
        for(int i=0;i<genomic_pos.size();++i){
            exon_size += (genomic_pos[i].end_ - genomic_pos[i].start_ + 1);
        }
    }
    // for CNV
    VariantSpecific(VariantType vartype,  std::string tran, std::vector<GenomicRegion> gpos, std::string ref_, std::string allele_, std::string strand_, int exon_start, int exon_end)
    : var_type(vartype)
    , transcript(tran)
    , genomic_pos(gpos)
    , ref(ref_)
    , allele(allele_)
    ,strand(strand_)
    ,start(exon_start)
    ,end(exon_end)
    {
        exon_size = 0;
        for(int i=0;i<genomic_pos.size();++i){
            exon_size += (genomic_pos[i].end_ - genomic_pos[i].start_ + 1);
        }
    }
    bool operator==(const VariantSpecific& other) const {
        if(var_type == other.var_type && ref == other.ref && allele == other.allele ){
            if(other.genomic_pos.size()==genomic_pos.size() ){
                for(int i=0;i<genomic_pos.size();++i){
                    if( !(genomic_pos[i] == other.genomic_pos[i])) return false;
                }
                return true;
            }else{
                return false;
            }
        }else{
            return false;
        }
    }
    ~VariantSpecific(){
        if( !genomic_pos.empty()){
            std::vector<GenomicRegion>().swap(genomic_pos);
        }
    }
    VariantType var_type;
    std::vector<GenomicRegion> genomic_pos;
    std::string transcript="NA";
    std::string ref;
    std::string allele;
    std::string infor;
    int exon_size;
    std::string strand="+";
    int start=0; // CNV start exon number
    int end =0; // CNV end exon number
};
class Variant {
public:
    Variant(std::string variant_file, int chr_index, int start_index, int end_index, int ref_index, int allele_index, int var_type_index)
        : variant_file_(variant_file),chr_index_(chr_index),start_index_(start_index)
        ,end_index_(end_index),ref_index_(ref_index), allele_index_(allele_index)
        , var_type_index_(var_type_index){
            readVariant();
    }
    Variant(std::string variant_file):variant_file_(variant_file){
        chr_index_ = 0;
        start_index_ = 1;
        end_index_ = 2;
        ref_index_ = 3;
        allele_index_ = 4;
        var_type_index_ = 5;
        readVariant();
    }
    Variant(std::string variant_file, std::string type): variant_file_(variant_file), type_(type){
        if( type_ == "SNV" ){
            chr_index_ = 0;
            start_index_ = 1;
            end_index_ = 2;
            ref_index_ = 3;
            allele_index_ = 4;
            var_type_index_ = 5;
            readVariant();
        }
    }
    VariantType getVarType(const std::string tmp_str);
    VariantSpecific Index( const int idx) { return variant_spec[idx];}
    void readVariant();
    void readCNV(std::string bed_file);
    void readSV();
    const int size() const { return variant_spec.size();}
    bool check_exist(const VariantSpecific& other);

private:
    std::string variant_file_;
    std::vector<VariantSpecific> variant_spec;
    std::string type_;
    int chr_index_;
    int start_index_;
    int end_index_;
    int ref_index_;
    int allele_index_;
    int var_type_index_;
};
std::string decodeVarType(VariantType var_type);
}