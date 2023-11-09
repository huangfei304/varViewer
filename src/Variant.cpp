//
// varViewer
// Date: 2023-08-09.
// Author: Huang Fei <huangfei@genomics.cn>
//

#include <fstream>
#include <sstream>
#include <stdexcept>
#include <boost/algorithm/string.hpp>
#include <map>
#include "utils.hh"
#include "Variant.hh"

namespace varViewer
{
VariantType Variant::getVarType(const std::string varStr_) {
    std::string tmp_str = varStr_;
    std::transform(tmp_str.begin(), tmp_str.end(), tmp_str.begin(), ::toupper);
    if( tmp_str == "SNV" ){
        return VariantType::SNP;
    }else if(tmp_str == "INS" ){
        return VariantType::INS;
    } else if(tmp_str == "DEL"){
        return VariantType::DEL;
    }else if( tmp_str == "DUP"){
        return VariantType::DUP;
    }else if( tmp_str =="INDEL" || tmp_str=="DELINS"){
        return VariantType::INDEL;
    }else if( tmp_str == "CNV"){
        return VariantType::CNV;
    }else if( tmp_str == "LOSS" ){
        return VariantType::LOSS;
    }else if( tmp_str == "LOSSDUP" ){
        return VariantType::LOSSDUP;
    }else if(tmp_str == "SV" ){
        return VariantType::SV;
    }else{
        throw std::logic_error("Encountered invalid variant type: " + varStr_);
    }
}
bool Variant::check_exist(const VariantSpecific& other){
    if( variant_spec.size()== 0 ){
        return false;
    }
    for(int i=0; i<variant_spec.size(); ++i) {
        if( variant_spec[i] == other ){
            return true;
        }
    }
    return false;
}
void Variant::readVariant(){
    std::string line;
    std::ifstream varFile(variant_file_);
    if (varFile.is_open()) {
        while (getline(varFile, line)){
            std::vector<std::string> pieces;
            boost::split(pieces, line, boost::is_any_of("\t"));
            //vector<string> pieces = splitStringByDelimiter(line,'\t');
            if (pieces[0] == "Scale") continue;
            VariantType var_type = getVarType(pieces[var_type_index_]);
            int64_t sstart,eend;
            std::stringstream sstream,estream;
            sstream<<pieces[start_index_];sstream>>sstart;
            estream<<pieces[end_index_];estream>>eend;
            GenomicRegion genomic_region(pieces[chr_index_], sstart, eend);
            VariantSpecific single_var(var_type, {genomic_region}, pieces[ref_index_], pieces[allele_index_]);
            single_var.infor=pieces[chr_index_]+"_"+pieces[start_index_]+"_"+pieces[end_index_];
            if( !check_exist(single_var) ){
                variant_spec.push_back(single_var);
            }
        }
        varFile.close();
    }else{
        throw std::runtime_error("Unable to open SNV file: " + variant_file_);
    }
}
void Variant::readSV(){
    std::string line;
    std::ifstream varFile(variant_file_);
    if (varFile.is_open()) {
        while (getline(varFile, line)){
            std::vector<std::string> pieces;
            boost::split(pieces, line, boost::is_any_of("\t"));
            //vector<string> pieces = splitStringByDelimiter(line,'\t');
            if (pieces[0].substr(0,5) == "Judge") continue;

            std::vector<std::string> left_bp;
            std::vector<std::string> right_bp;
            //std::vector<GenomicRegion> grs;
            boost::split(left_bp, pieces[3], boost::is_any_of(":"));
            boost::split(right_bp, pieces[4], boost::is_any_of(":"));
            //VariantType var_type = getVarType(pieces[var_type_index_]);
            int64_t sstart,eend;
            std::stringstream sstream,estream;
            sstream<<left_bp[1];sstream>>sstart;
            estream<<right_bp[1];estream>>eend;
            GenomicRegion gr_left(left_bp[0], sstart, sstart);
            GenomicRegion gr_right(right_bp[0], eend, eend);
            VariantSpecific single_var(VariantType::SV, {gr_left, gr_right}, "", "");
            single_var.infor=pieces[1]+"_"+pieces[2];
            if( !check_exist(single_var) ){
                variant_spec.push_back(single_var);
            }
        }
        varFile.close();
    }else{
        throw std::runtime_error("Unable to open SV file: " + variant_file_);
    }
}
void Variant::readCNV(std::string bed_file){
    std::string line;

    //format: chr<tab>exon_start<tab>exon_end<tab>transcript<tab>strand
    std::map<std::string,std::vector<GenomicRegion>> tran_range;
    std::map<std::string, std::string> tran_strand;
    std::ifstream bedFile(bed_file);
    if (bedFile.is_open()) {
        while (getline(bedFile, line)){
            std::vector<std::string> pieces;
            boost::split(pieces, line, boost::is_any_of("\t"));

            int64_t sstart,eend;
            std::stringstream sstream,estream;
            sstream<<pieces[1];sstream>>sstart;
            estream<<pieces[2];estream>>eend;

            std::map<std::string,std::vector<GenomicRegion>>::iterator iter=tran_range.find(pieces[3]);
            if( iter != tran_range.end()){
                GenomicRegion genomic_region(pieces[0],sstart,eend);
                tran_range.at(pieces[3]).push_back(genomic_region);
            }else{
                GenomicRegion genomic_region(pieces[0], sstart, eend);
                std::vector<GenomicRegion> gps {genomic_region};
                //tran_range.insert(make_pair(i->first, gps));
                tran_range.emplace(pieces[3], gps);
                tran_strand.emplace(pieces[3], pieces[4]);
            }
        }
        bedFile.close();
    }else{
        throw std::runtime_error("Unable to open CNV bedfile: " + bed_file);
    }

    //format: chr<tab>varType<tab>transcript<tab>exon_start<tab>exon_end
    std::ifstream varFile(variant_file_);
    if (varFile.is_open()) {
        while (getline(varFile, line)){
            std::vector<std::string> pieces;
            boost::split(pieces, line, boost::is_any_of("\t"));
            //vector<string> pieces = splitStringByDelimiter(line,'\t');
            if (pieces[0] == "#Gene") continue;
            VariantType var_type = getVarType("CNV");
            std::transform(pieces[1].begin(), pieces[1].end(), pieces[1].begin(), ::toupper);
            std::string tmp_str = pieces[1].substr(0,3);
            if (pieces[1] == "COPY NUMBER LOSS" || tmp_str == "DEL" ){
                var_type = getVarType("LOSS");
            }else{
                var_type = getVarType("DUP");
            }

            int sstart,eend,raw_start,raw_end;
            std::stringstream sstream,estream;
            sstream<<pieces[3];sstream>>sstart;
            estream<<pieces[4];estream>>eend;
            if( tran_range.find(pieces[2]) == tran_range.end()){
                std::cerr<<"The gene:"<<pieces[2]<<" is not in the BED file, please check."<<std::endl;
                exit(1);
            }
            int exon_num = tran_range.at(pieces[2]).size();
            if( sstart > exon_num || eend > exon_num ){
                std::cerr<<"The gene:"<<pieces[2]<<" exon number ("<<std::to_string(sstart)<<"-"<<std::to_string(eend)<<")in CNV is large than total exon number:"<<std::to_string(exon_num)<<std::endl;
                exit(1);
            }
            raw_start=sstart;
            raw_end = eend;
            if ( tran_strand.at(pieces[2]) == "-" ){
                sstart = (exon_num - sstart + 1);
                eend = (exon_num - eend + 1);
            }
            if( sstart > eend ){
                int tmp = sstart;
                sstart = eend;
                eend = tmp;
            }
            VariantSpecific single_var(var_type, pieces[2], tran_range.at(pieces[2]), "", "", tran_strand.at(pieces[2]), sstart, eend);
            variant_spec.push_back(single_var);
            //chrName_geneName_EXstart_EXend_DEL
            std::string infor = pieces[0]+"_"+pieces[2]+"_EX"+std::to_string(raw_start)+"_EX"+std::to_string(raw_end);
            if( var_type == VariantType::DEL ){
                infor += "_DEL";
            }else if( var_type == VariantType::DUP){
                infor +="_DUP";
            }else if( var_type == VariantType::LOSS ){
                infor +="_DEL";
            }else{
                infor +="_CNV";
            }
            variant_spec.back().infor=infor;
        }
        varFile.close();
    }else{
        throw std::runtime_error("Unable to open CNV file: " + variant_file_);
    }
    //for(auto iter=transcript.begin();iter !=transcript.end(); ++iter){
    //    std::cerr<<"transcript:"<<iter->first<<std::endl;
    //}
//    for (auto iter = tran_range.begin(); iter != tran_range.end(); ++iter) {
//        VariantSpecific single_var(transcript.at(iter->first), iter->second, "", "",tran_strand.at(iter->first));
//        if( !check_exist(single_var) ){
//            variant_spec.push_back(single_var);
//        }
//    }
}
std::string decodeVarType(VariantType var_type)
{
    switch (var_type)
    {
    case VariantType::SNP:
        return "SNP";
    case VariantType::INS:
        return "INS";
    case VariantType::DEL:
        return "DEL";
    case VariantType::INDEL:
        return "INDEL";
    case VariantType::CNV:
        return "CNV";
    case VariantType::LOSS:
        return "DEL";
    case VariantType::DUP:
        return "DUP";
    default:
        std::cerr<<" unknown Varuant Type"<<std::endl;
        exit(1);
    }
}
}

