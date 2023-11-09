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

#include "LanePlot.hh"
#include <stdexcept>
#include <tuple>
#include <boost/optional.hpp>
#include <algorithm>
#include <iostream>

using boost::optional;
using std::list;
using std::string;
using std::unordered_map;
using std::vector;
using std::logic_error;

namespace varViewer
{
// TODO: Rename segment to displaySegment?
static Segment getSegment(Frag& frag, const std::string type, const double opacity=0.7) {
    Read this_read = frag.read;
    Read mate_read = frag.mate;
    std::list<Operation> operations_ = this_read.align.operations();
    int start = this_read.align.referenceStart();
    int end = this_read.align.referenceStart() + this_read.align.referenceLength();
    vector<Feature> features;
    //+: #FFB6C1(lightpink)   #8470FF(lightslateblue)
    //"#8da0cb"   "#fc8d62"
    string fill = (this_read.strand ) ? "#ffb6c1" : "#8470ff";
	list<Operation>::const_iterator operation_it = operations_.begin();
	list<std::string>::const_iterator base_it = this_read.bases.begin();
	if( operation_it->type()==OperationType::kSoftclip){
        start -= operation_it->length();
    }
    //std::cerr<<"Read start:"<<std::to_string(start)<<",end:"<<std::to_string(end)<<std::endl;
    //std::cerr<<"Read:"<<this_read<<std::endl;
    while (operation_it != operations_.end()) {
        const int opLength = operation_it->length();
        if (operation_it->type() == OperationType::kMatch) {
            Feature this_feature(FeatureType::kRect, opLength, fill, "none",1.0);
            this_feature.label = *base_it;
            features.push_back(this_feature);
        }else if (operation_it->type() == OperationType::kMismatch) {
            Feature this_feature(FeatureType::kRect, opLength, "#cdcdcd", "none",1.0);
            this_feature.label = *base_it;
            features.push_back(this_feature);
        } else if (operation_it->type() == OperationType::kDeletionFromRef){
            Feature this_feature(FeatureType::kLine, opLength, "none", "black");
			features.push_back(this_feature);
        } else if (operation_it->type() == OperationType::kSoftclip) {
            Feature this_feature(FeatureType::kRect, opLength, "none", "none", 1.0);
            this_feature.label = *base_it;
			features.push_back(this_feature);
        } else if (operation_it->type()  == OperationType::kInsertionToRef) {
            Feature this_feature(FeatureType::kVerticalLine, 0, "none", "black",1.0);
			features.push_back(this_feature);
        } else {// kMissingBases
            throw logic_error(" operation: "+operation_it->generateCigar() +" is not right.");
        }
        ++operation_it;
        ++base_it;
    }
    if( mate_read.count > 0 ){
        operations_ = mate_read.align.operations();
        operation_it = operations_.begin();
        base_it = mate_read.bases.begin();
        fill = (mate_read.strand ) ? "#ffb6c1" : "#8470ff";
        int mate_start = mate_read.align.referenceStart();
        if( operation_it->type()==OperationType::kSoftclip){
            mate_start -= operation_it->length();
        }
        int insert_len = mate_start-1 - end;
        if( insert_len > 0 ){
            Feature this_feature(FeatureType::kLine, insert_len, "#E5E5E5", "#E5E5E5");
            features.push_back(this_feature);
        }
        //std::cerr<<"mate-start:"<<std::to_string(mate_start)<<",insert_len:"<<std::to_string(insert_len)<<std::endl;
        //std::cerr<<mate_read<<std::endl;
        while (operation_it != operations_.end()) {
            const int opLength = operation_it->length();
            if (operation_it->type() == OperationType::kMatch) {
                Feature this_feature(FeatureType::kRect, opLength, fill, "none",1.0);
                this_feature.label = *base_it;
                features.push_back(this_feature);
            }else if (operation_it->type() == OperationType::kMismatch) {
                Feature this_feature(FeatureType::kRect, opLength, "#cdcdcd", "none",1.0);
                this_feature.label = *base_it;
                features.push_back(this_feature);
            } else if (operation_it->type() == OperationType::kDeletionFromRef){
                Feature this_feature(FeatureType::kLine, opLength, "none", "black");
			    features.push_back(this_feature);
            } else if (operation_it->type() == OperationType::kSoftclip) {
                Feature this_feature(FeatureType::kRect, opLength, "none", "none",1.0);
                this_feature.label = *base_it;
			    features.push_back(this_feature);
            } else if (operation_it->type()  == OperationType::kInsertionToRef) {
                Feature this_feature(FeatureType::kVerticalLine, 0, "none", "black");
			    features.push_back(this_feature);
            } else {// kMissingBases
                throw logic_error(" operation: "+operation_it->generateCigar() +" is not right.");
            }
            ++operation_it;
            ++base_it;
        }
    }
    if( type == "SNV" && this_read.count > 1 ){
        Feature this_feature(FeatureType::kLabel, 4, "none", "none",1.0);
        this_feature.label = "x"+std::to_string(this_read.count);
        features.push_back(this_feature);
    }
    Segment this_segment(start, features, opacity);
    return this_segment;
}

static std::vector<Segment> RefLane(Aligns& this_align, std::string type) {
    vector<Segment> seg_vec;
    vector<Feature> features;
    if( type == "SNV"){
        string fill="url(#OrangeWhiteOrange)";
        int var_length = this_align.get_end() - this_align.get_start() + 1;
        int var_start = this_align.get_start();
        int leftflank_length = this_align.get_start() - (this_align.get_refer_start() + this_align.get_start_offset()) - 1;
        int leftflank_start = this_align.get_start() - leftflank_length;
        int rightflank_start = this_align.get_end()+1;
        int rightflank_length = this_align.get_end_offset() -this_align.get_start_offset() - leftflank_length - var_length;
        std::string ntseq = this_align.get_reference().getSequence(this_align.get_chr(), leftflank_start, (leftflank_start+leftflank_length+var_length+rightflank_length));
        std::transform(ntseq.begin(), ntseq.end(), ntseq.begin(), ::toupper);
        //std::cerr<<"left length:"<<std::to_string(leftflank_length)<<",right length:"<<std::to_string(rightflank_length)<<std::endl;
        //std::cerr<<"sequence:"<<ntseq<<std::endl;
        // leftflank
        fill="url(#BlueWhiteBlue)";
        features.emplace_back(FeatureType::kRectWithLeftBreak, leftflank_length, fill, "black");
        features.back().label = ntseq.substr(0, leftflank_length);
        //variant
        fill="url(#OrangeWhiteOrange)";
        features.emplace_back(FeatureType::kRect, var_length, fill, "black",1.0);
        features.back().label = ntseq.substr(leftflank_length,var_length);
        // rightflank
        fill="url(#BlueWhiteBlue)";
        features.emplace_back(FeatureType::kRectWithRightBreak, rightflank_length, fill, "black");
        features.back().label = ntseq.substr(leftflank_length+var_length,rightflank_length);
    }else{//SV
         FusionType fusion_type = this_align.get_fusion_type();
         int leftflank_length, leftflank_start, rightflank_start, rightflank_length,left_bp,right_bp, left_most;
         left_bp = int(this_align.get_start());
         right_bp = int(this_align.get_end());
         left_most = this_align.get_refer_start();
         if( fusion_type == FusionType::FLSR ){
            leftflank_length = left_bp - (left_most + this_align.get_start_offset()) - 1;
            leftflank_start = left_bp - leftflank_length;
            rightflank_length = left_most + this_align.get_end_offset() - left_bp;
            rightflank_start = right_bp;
         }else if( fusion_type == FusionType::FLrevSL){
            leftflank_length = left_bp - (left_most + this_align.get_start_offset());
            leftflank_start = left_bp - leftflank_length - 1;
            rightflank_length = left_most + this_align.get_end_offset() - left_bp;
            rightflank_start = right_bp - rightflank_length - 1;
         }else if( fusion_type == FusionType::SLFR ){
         }else if( fusion_type == FusionType::revSRFR){
         }

        std::cerr<<"left chr:"<<this_align.get_chr()<<",start:"<<std::to_string(leftflank_start)<<",length:"<<std::to_string(leftflank_length)<<std::endl;
        std::cerr<<"right chr:"<<this_align.get_end_chr()<<",start:"<<std::to_string(rightflank_start)<<",length:"<<std::to_string(rightflank_length)<<std::endl;
        std::string left_seq = this_align.get_reference().getSequence(this_align.get_chr(), leftflank_start, (leftflank_start+leftflank_length));
        std::transform(left_seq.begin(), left_seq.end(), left_seq.begin(), ::toupper);
        if(this_align.get_fusion_type() == FusionType::revSRFR){
                left_seq = reverseComplement(left_seq);
        }

        // leftflank
        string fill="url(#BlueWhiteBlue)";
        features.emplace_back(FeatureType::kRectWithLeftBreak, leftflank_length, fill, "black");
        features.back().label = left_seq;
        // rightflank
        std::string right_seq = this_align.get_reference().getSequence(this_align.get_end_chr(), rightflank_start, (rightflank_start+rightflank_length));
        std::transform(right_seq.begin(), right_seq.end(), right_seq.begin(), ::toupper);
        if(this_align.get_fusion_type() == FusionType::FLrevSL){
            right_seq = reverseComplement(right_seq);
        }
        fill="url(#OrangeWhiteOrange)";
        features.emplace_back(FeatureType::kRectWithRightBreak, rightflank_length, fill, "black");
        features.back().label = right_seq;
    }

    Segment flank_segment(0, features, 1.0);
    seg_vec.push_back(flank_segment);
    return seg_vec;
}
static std::vector<Segment> ChrCoordLane(int start, int end, int first_exon_start, int offset=100){
    vector<Segment> seg_vec;
    vector<Feature> features;
    //divide into 7 parts
    float part = (end+offset - (start-offset)+1)/7;
    for (int i=0; i<6; ++i ) {
        int index = int(0.5+(first_exon_start-offset) + (i+1)*part);
        features.emplace_back(FeatureType::kLineSegment, part, "black", "black");
        features.back().label = std::to_string(index);
    }

    features.emplace_back(FeatureType::kLine, part, "black", "black", 0.0);
    Segment flank_segment(0, features, 1.0);
    seg_vec.push_back(flank_segment);
    return seg_vec;
}
static std::vector<Segment> TransLane(std::vector<GenomicRegion>&gregion, std::string strand, int offset=100, int exon_start=1, int exon_end=1){
    vector<Segment> seg_vec;
    vector<Feature> features;
    int start = gregion.front().start_ - offset;
    int end = gregion.back().end_ + offset;

    int exon_num = gregion.size();
    if( offset > 0 ){
        if( exon_num > 1 ){
            features.emplace_back(FeatureType::kLine, offset, "black","black",0.5);
        }else{
            if( strand == "+" ){
                features.emplace_back(FeatureType::kArrows, offset, "black","black",0.5, 0.0, 1.0);
            }else{
                features.emplace_back(FeatureType::kArrows, offset, "black","black",0.5,1.0, 0.0);
            }
        }
    }
    for (int i = 0; i< gregion.size(); ++i) {
        int size = gregion[i].end_ - gregion[i].start_ + 1;
        if( (i+1) >= exon_start && (i+1) <= exon_end ){
            features.emplace_back(FeatureType::kRect, size, "red", "red",1.0);
        }else{
            features.emplace_back(FeatureType::kRect, size, "black", "black",1.0);
        }
        if( i < (gregion.size()-1) ){
            int next_size = gregion[i+1].start_ - gregion[i].end_;
            if( next_size > (end-start+1)/50 ){
                if( strand == "+"){
                    features.emplace_back(FeatureType::kArrows, next_size, "black","black",0.5, 0.0, 1.0);
                }else{
                    features.emplace_back(FeatureType::kArrows, next_size, "black","black",0.5, 1.0, 0.0);
                }
            }else{
                features.emplace_back(FeatureType::kLine, next_size, "black","black",0.5);
            }
        }
    }
    if( offset > 0){
        if(exon_num > 1){
            features.emplace_back(FeatureType::kLine, offset, "black","black",0.5);
        }else{
           if( strand == "+" ){
                features.emplace_back(FeatureType::kArrows, offset, "black","black",0.5, 0.0, 1.0);
            }else{
                features.emplace_back(FeatureType::kArrows, offset, "black","black",0.5,1.0, 0.0);
            }
        }
    }

    Segment flank_segment(0, features, 1.0);
    seg_vec.push_back(flank_segment);
    return seg_vec;
}

static std::vector<Segment> CoverageLane(std::vector<GenomicRegion>&gregion, std::vector<int> this_cover, int left_most, int offset, int min_depth, int max_depth, std::string color="red", float opacity=0.6){
    vector<Segment> seg_vec;
    vector<Feature> features;
    int deta_depth = max_depth - min_depth;
    vector<float> depth;
    /*
    * if( offset >0 ){
    *        features.emplace_back(FeatureType::kLine, offset, "black", "black",0.0);
    * }

    * for (int i = 0; i< gregion.size(); ++i) {
    *    int size = int(gregion[i].end_ - gregion[i].start_ +1);
    *    int start_offset = gregion[i].start_ - left_most;
    *    int end_offset = gregion[i].end_ - left_most;
    *    vector<float> depth;
    *    for(int k=start_offset;k<=end_offset;++k){
    *        depth.push_back((this_cover[k]-min_depth)/deta_depth);
    *    }

    *    features.emplace_back(FeatureType::kPath, size, "none",color);
    *    features.back().depth = depth;
    *    if( i < (gregion.size()-1) ){
    *        int next_size = gregion[i+1].start_ - gregion[i].end_;
    *        features.emplace_back(FeatureType::kLine, next_size, "black","black",0.0);
    *    }
    * }

    * if( offset > 0){
    *    features.emplace_back(FeatureType::kLine, offset, "black","black",0.0);
    * }
    */

    /*
    * for(int i=0;i<offset;i++){
    *      depth.push_back(0.0);
    * }
    * for (int i = 0; i< gregion.size(); ++i) {
    *    int size = int(gregion[i].end_ - gregion[i].start_ +1);
    *    int start_offset = gregion[i].start_ - left_most;
    *    int end_offset = gregion[i].end_ - left_most;

    *    for(int k=start_offset;k<=end_offset;++k){
    *        float value = 1.0*(this_cover[k] - min_depth)/deta_depth;
    *        depth.push_back(value);
    *    }
    *   if( i < (gregion.size()-1) ){
    *        int next_size = gregion[i+1].start_ - gregion[i].end_;
    *        for(int j=0;j<next_size;j++){
    *            depth.push_back(0.0);
    *        }
    *   }
    *}

    *for(int i=0;i<offset;i++){
    *    depth.push_back(0.0);
    * }
    */

    int extend = 40;
    for(int i=0;i<(offset-extend);i++){
        depth.push_back(0.0);
    }
    for (int i = 0; i< gregion.size(); ++i) {
        int size = int(gregion[i].end_ - gregion[i].start_ +1);
        int start_offset = gregion[i].start_ - left_most;
        int end_offset = gregion[i].end_ - left_most;
        if( i == 0 ){ // first one
            start_offset -= extend;
        }else{
            int next_size = gregion[i].start_ - gregion[i-1].end_;
            if( next_size >=2*extend ){
                start_offset -= extend;
            }else{
                start_offset -= (next_size/2);
            }
        }
        if( i == (gregion.size()-1) ){
            end_offset += extend;
        }else{
            int next_size = gregion[i+1].start_ - gregion[i].end_;
            if(next_size >=2*extend){
                end_offset += extend;
            }else{
                end_offset +=(next_size/2);
            }
        }
        for(int k=start_offset;k<=end_offset;++k){
            float value = (this_cover[k] > min_depth) ? 1.0* (this_cover[k] - min_depth) : 0.0;
            if( value > deta_depth ) value = deta_depth; // may be extend depth > max ???
            //std::cerr<<"Index: "<<std::to_string(k)<<",depth:"<<std::to_string(this_cover[k])<<std::endl;
            depth.push_back(float(value/deta_depth));
        }
        if( i < (gregion.size()-1) ){
            int next_size = gregion[i+1].start_ - gregion[i].end_;
            if( next_size > 2*extend ){
                for(int j=0;j<(next_size-2*extend);j++){
                    depth.push_back(0.0);
                }
            }
        }
    }
    for(int i=0;i<(offset-extend);i++){
        depth.push_back(0.0);
    }

    features.emplace_back(FeatureType::kPath, depth.size(), "none",color);
    features.back().depth = depth;
    Segment cover_segment(0, features, opacity);
    seg_vec.push_back(cover_segment);
    return seg_vec;
}
static std::vector<Segment> RDratioLane(std::vector<GenomicRegion>&gregion, int index, std::vector<int> treat_cover, std::vector<int> ref_cover, int left_most, int offset,std::string color="black"){
    std::vector<Segment> seg_vec;
    std::vector<Feature> features;

    int dist = offset;
    for (int i = 0; i< gregion.size(); ++i) {
        int size = int(gregion[i].end_ - gregion[i].start_ +1);
        int start_offset = gregion[i].start_ - left_most;
        int end_offset = gregion[i].end_ - left_most;
        if( i == index ){
            int treat_total=0;
            int ref_total=0;
            for(int j=start_offset;j<=end_offset;++j){
                treat_total += treat_cover[j];
                ref_total += ref_cover[j];
            }
            std::string rate="0.00";
            if( ref_total > 0 ){
                float ratio = 1.0*treat_total/ref_total;
                //ratio = float(int(100.0*ratio +0.5))/100;
                rate=std::to_string(ratio);
            }else{
                std::cerr<<"The region ["<<std::to_string(gregion[i].start_)<<","<<std::to_string(gregion[i].end_)<<"] in reference set depth is 0, please check."<<std::endl;
                exit(1);
            }
            if( index % 3 == 0 ){
                if( index == 0 ){
                    features.emplace_back(FeatureType::kLabel, 10, color, color, 1.0);
                }else{
                    features.emplace_back(FeatureType::kLabel, 40, color, color, 1.0);
                }
            }else if( index % 3 == 1 ){
                features.emplace_back(FeatureType::kLabel, 40, color, color, 0.9);
            }else{
                features.emplace_back(FeatureType::kLabel, 40, color, color, 0.8);
            }
            features.back().label=rate.substr(0,4);
            Segment rate_segment(dist+size/2, features, 1.0);
            seg_vec.push_back(rate_segment);
            break;
        }
        if( (i+1) < gregion.size() ){
            dist += (gregion[i+1].start_ - gregion[i].start_);
        }
    }
    return seg_vec;
}

static std::vector<Segment> GCLane(std::vector<GenomicRegion>&gregion, std::vector<float>cover,int offset, string color="black"){
    std::vector<Segment> seg_vec;
    vector<Feature> features;
    int start = gregion.front().start_ - offset;
    int end = gregion.back().end_ + offset;
//    if( offset > 0 ){
//        features.emplace_back(FeatureType::kLine, offset, "black","black",0.0);
//    }
//    for (int i = 0; i< gregion.size(); ++i) {
//        int size = int(gregion[i].end_ - gregion[i].start_ +1);
//        features.emplace_back(FeatureType::kLine, size+20, color,color, cover[i]);
//        if( i < (gregion.size()-1) ){
//            int next_size = gregion[i+1].start_ - gregion[i].end_;
//            features.emplace_back(FeatureType::kLine, next_size-20, "black","black",0.0);
//        }
//    }
//    if( offset > 0){
//        features.emplace_back(FeatureType::kLine, offset, "black","black",0.0);
//    }

    int extend = 20;
    if( offset > extend ){
        features.emplace_back(FeatureType::kLine, (offset-extend),"black", "black",0.0);
    }

    for(int i=0; i<gregion.size();++i){
        int size = int(gregion[i].end_ - gregion[i].start_ + 1);
        if(i==0){
            size += extend;
        }else{
            int next_size = gregion[i].start_ - gregion[i-1].end_;
            if(next_size >= 2*extend){
                size += extend;
            }else{
                size += next_size/2;
            }
        }
        if( i == (gregion.size()-1)){
            size += extend;
        }else{
            int next_size = gregion[i+1].start_ - gregion[i].end_;
            if(next_size >= 2*extend ){
                size += extend;
            }else{
                size += (next_size/2);
            }
        }

        features.emplace_back(FeatureType::kLine, size, color,color, cover[i]);

        if( i < (gregion.size()-1) ){
            int next_size = gregion[i+1].start_ - gregion[i].end_;
            if( next_size > 2*extend ){
                features.emplace_back(FeatureType::kLine, (next_size-2*extend), "black","black",0.0);
            }
        }
    }

    if( offset > extend ){
        features.emplace_back(FeatureType::kLine, (offset-extend),"black", "black",0.0);
    }

    Segment flank_segment(0, features, 1);
    seg_vec.push_back(flank_segment);

    return seg_vec;
}
static std::vector<Segment> CoverLane(std::vector<GenomicRegion>&gregion, std::vector<float>cover,std::vector<float>cover30, int offset, string color="black",string color3="blue",float min=0.9,float max=1.0){
    std::vector<Segment> seg_vec;
    vector<Feature> features;
    int start = gregion.front().start_ - offset;
    int end = gregion.back().end_ + offset;
    float deta = max - min;
//    if( offset > 0 ){
//        features.emplace_back(FeatureType::kLine, offset, "black","black",0.0);
//    }
//    for (int i = 0; i< gregion.size(); ++i) {
//        int size = int((gregion[i].end_ - gregion[i].start_ + 1+1)/2);
//        float value_1 = (cover[i] > min) ? (cover[i] - min ) : 0.0;
//        features.emplace_back(FeatureType::kRect, size, color, color, value_1/deta);
//        float value_30 = (cover30[i] > min) ? (cover30[i] - min ) : 0.0;
//        features.emplace_back(FeatureType::kRect, size, color3, color3, value_30/deta);
//        if( i < (gregion.size()-1) ){
//            int next_size = gregion[i+1].start_ - gregion[i].end_;
//            features.emplace_back(FeatureType::kLine, next_size, "black","black",0.0);
//        }
//    }
//    if( offset > 0){
//        features.emplace_back(FeatureType::kLine, offset, "black","black",0.0);
//    }

    int extend = 40;
    if( offset > extend ){
        features.emplace_back(FeatureType::kLine, (offset-extend),"black", "black",0.0);
    }
    for(int i=0; i<gregion.size();++i){
        int size = int(gregion[i].end_ - gregion[i].start_ + 1);
        if(i==0){
            size += extend;
        }else{
            int next_size = gregion[i].start_ - gregion[i-1].end_;
            if(next_size >= 2*extend){
                size += extend;
            }else{
                size += next_size/2;
            }
        }
        if( i == (gregion.size()-1)){
            size += extend;
        }else{
            int next_size = gregion[i+1].start_ - gregion[i].end_;
            if(next_size >= 2*extend ){
                size += extend;
            }else{
                size += (next_size/2);
            }
        }
        float value_1 = (cover[i] > min) ? (cover[i] - min ) : 0.0;
        features.emplace_back(FeatureType::kRect, size/2, color, color, value_1/deta);
        float value_30 = (cover30[i] > min) ? (cover30[i] - min ) : 0.0;
        features.emplace_back(FeatureType::kRect, size/2, color3, color3, value_30/deta);

        if( i < (gregion.size()-1) ){
            int next_size = gregion[i+1].start_ - gregion[i].end_;
            if( next_size > 2*extend ){
                features.emplace_back(FeatureType::kLine, (next_size-2*extend), "black","black",0.0);
            }
        }
    }
    if( offset > extend ){
        features.emplace_back(FeatureType::kLine, (offset-extend),"black", "black",0.0);
    }

    Segment flank_segment(0, features, 1);
    seg_vec.push_back(flank_segment);

    return seg_vec;
}
static std::vector<Segment> CoordLane(int start, int end, std::string label="xxx", int pad=100, std::string color="black", float min=0, float max=1.0){
    vector<Segment> seg_vec;
    vector<Feature> features;
    features.emplace_back(FeatureType::kCoord, (end-start+1+2*pad), color, color, 0.0, min, max);
    features.back().label=label;
    Segment flank_segment(0, features, 1.0);
    seg_vec.push_back(flank_segment);

    return seg_vec;
}
static std::vector<Segment> MapLane(std::vector<GenomicRegion>&gregion, std::vector<float>map, int left_most, int offset, string color="red" ){
    int start = gregion.front().start_ - offset;
    int end = gregion.back().end_ + offset;
    vector<Segment> seg_vec;
    vector<Feature> features;
    vector<float>depth;

    int per =1000;
    float part = 1.0;
    if( (end-start+1) > 2*per ){
        part = (end-start+1.0)/per;
    } else {
        per = (end-start+1);
    }
    for(int i=0;i<per;++i) {
        int ss = int(start + i*part+0.5);
        int ee = int(start+ (i+1)*part+0.5);
        float sum = 0;
        for(int pos=ss;pos<ee;++pos){
            int idx = pos - left_most;
            sum += map[idx];
        }
        float value = sum/(ee-ss);
        depth.push_back(value);
    }
    features.emplace_back(FeatureType::kPath, (end-start+1), "none",color);
    features.back().depth = depth;
    //features.back().label = "Mapability";
    Segment flank_segment(0, features, 1.0);
    seg_vec.push_back(flank_segment);

    return seg_vec;
}
static std::vector<Segment> CABANALane(std::vector<GenomicRegion>&gregion, int left_most, std::vector<int>treat_depth, std::vector<int>ref_depth, int offset, string color="red",float opacity=0.6){
    vector<Segment> seg_vec;
    vector<Feature> features;
    vector<float> nrd;//the normalized read depth

    int size = offset*2 + (gregion.back().end_ - gregion.front().start_ + 1);
    for (int i = 0; i< gregion.size(); ++i) {
        int size = int(gregion[i].end_ - gregion[i].start_ +1);
        int start_offset = gregion[i].start_ - left_most;
        int end_offset = gregion[i].end_ - left_most;
        for(int k=start_offset;k<=end_offset;++k){
            float value=0.0;
            if( ref_depth[k] == 0 ){
                if( treat_depth[k]==0) value = 0.0;
                else value = 1.0;
            }else{
                value = 1.0*(treat_depth[k] - ref_depth[k])/ref_depth[k];
            }
            if( value > 1.0 ) value = 1.0;
            else if(value < -1.0) value = -1.0;
            nrd.push_back(value);
        }
    }

    features.emplace_back(FeatureType::kPath, size, "none", color, 0.5, -1.0, 1.0);
    features.back().depth = nrd;
    Segment cover_segment(0, features, opacity);
    seg_vec.push_back(cover_segment);
    return seg_vec;
}

LanePlot generateLanePlot(Aligns& align, opts_s& opts, Coverage& cover) {
    //std::cerr<<"GenerateLanPlot start..."<<std::endl;
    //std::vector<Frag> align_read = align.AlignReads;
    LanePlot this_laneplot;
    std::string type="SNV";
    float lane_height = (opts.lane >15.0) ? opts.lane : 15.0;
    if( strcmp(opts.type, "SNV") ==0 ){
        this_laneplot.push_back(Lane(lane_height, {RefLane(align,type)}));
    }else if(strcmp(opts.type, "SV")==0){
        type="SV";
        this_laneplot.push_back(Lane(lane_height, {RefLane(align,type)}));
    }
    //std::cerr<<"Finished RefLane."<<std::endl;
   //list<Segment> segments;
    for (int idx = 0; idx != align.AlignReads.size(); ++idx) {
        //std::cerr<<"Index: "<<std::to_string(idx)<<std::endl;
        auto this_segment = getSegment(align.AlignReads[idx],type);
        // merge non-overlapped segment
        //std::cerr<<"merge non-overlapped segment"<<std::endl;
        //std::cerr<<"Idx:"<<std::to_string(idx)<<", this_segment start:"<<std::to_string(this_segment.start)<<",end:"<<std::to_string(this_segment.end)<<std::endl;
        bool flag = true;
        for(int i=2;i<this_laneplot.size();++i){
            if( this_laneplot[i].end < this_segment.start){
                //std::cerr<<"laneplot index:"<<std::to_string(i)<<",start:"<<std::to_string(this_laneplot[i].start)<<",end:"<<std::to_string(this_laneplot[i].end);
                this_laneplot[i].multi_segments[0].push_back(this_segment);
                this_laneplot[i].end = this_segment.end;
                //std::cerr<<",end: "<<std::to_string(this_laneplot[i].end)<<std::endl;
                flag = false;
                break;
            }
        }
        std::vector<Segment> seg_vec;
        seg_vec.push_back(this_segment);
        if( flag ) this_laneplot.push_back(Lane(opts.lane,{seg_vec}));
        //segments.push_back(this_segment);
    }
    //this_laneplot.push_back(Lane(10, segments));
    return this_laneplot;
}
LanePlot generateLanePlot(CNV& cnv, opts_s& opts, int pad) {
    LanePlot this_laneplot;
    VariantSpecific this_var = cnv.get_VariantSpecific();
    int start = this_var.genomic_pos.front().start_;
    int end = this_var.genomic_pos.back().end_;
    int this_var_size = (end - start + 1);
    int this_var_exon_size = this_var.exon_size;

    std::cerr<<"[LanePlot] Chromosome Coord..."<<std::endl;
    this_laneplot.push_back(Lane(opts.lane-5, {ChrCoordLane(opts.flank, ((end-start+1)+opts.flank), start, pad)}));
    std::cerr<<"[LanePlot] Trancript Exon information..."<<std::endl;
    this_laneplot.push_back(Lane(opts.lane*2, {TransLane(this_var.genomic_pos, this_var.strand, pad, this_var.start, this_var.end )}));

    std::cerr<<"[LanePlot] GC Coverage Mappability..."<<std::endl;
    std::vector<std::vector<Segment>> multi_seg;
    //std::vector<std::vector<Segment>> multi_seg_cover;

    //GC, Mapability
    if( !cnv.get_mapability().empty() ){
        multi_seg.push_back(MapLane(this_var.genomic_pos, cnv.get_mapability(), cnv.get_ref_start(), pad, "grey"));
    }
    multi_seg.push_back(GCLane(this_var.genomic_pos, cnv.get_gc(), pad, "red"));
    multi_seg.push_back(CoordLane(opts.flank, ((end-start+1)+opts.flank), "GC/Mappability", pad, "black",0.0,100.0));
    this_laneplot.push_back(Lane(opts.lane*10,multi_seg));

    multi_seg.clear();
    //Coord, Coverage(>0,>30)
    float base = 95.0;
    multi_seg.push_back(CoverLane(this_var.genomic_pos, cnv.get_cover(), cnv.get_cover(30),pad, "red","blue", base/100, 1.00));
    //multi_seg.push_back(CoverGCLane(this_var.genomic_pos, cnv.get_cover(30), pad, "orange"));
    multi_seg.push_back(CoordLane(opts.flank, ((end-start+1)+opts.flank), "Coverage", pad, "black",base, 100.0));
    this_laneplot.push_back(Lane(opts.lane*10,multi_seg));

    std::cerr<<"[LanePlot] Sample Depth..."<<std::endl;
    // Normalized depth of reference and sample
    int cnv_size = cnv.cnv_info_.size();
    int min = cnv.get_min_depth();
    int max = cnv.get_max_depth();
    multi_seg.clear();
    //std::vector<std::vector<Segment>> multi_seg_depth;
    std::string color="grey";
    for(int i=1;i<cnv_size;i++){
        if( i == cnv_size-1){// just the mean depth of the reference set.
            color="blue";
            CNVinfo this_cnv_info = cnv.cnv_info_.at(i);
            std::vector<int> this_cnv = this_cnv_info.depth_;
            //for(int j=0;j<this_cnv_info.depth_.size();j++){
            //    std::cerr<<"Index:"<<std::to_string(i)<<","<<std::to_string(j)<<",Depth:"<<std::to_string(this_cnv_info.depth_[j])<<std::endl;
            //}
            vector<Segment> this_segs = CoverageLane(this_var.genomic_pos, this_cnv_info.depth_, cnv.get_ref_start(), pad, min, max, color,1.0);
            multi_seg.push_back(this_segs);
        }
    }

    //the analysis sample
    color="red";
    CNVinfo this_cnv_info = cnv.cnv_info_.at(0);
    vector<Segment> this_segs = CoverageLane(this_var.genomic_pos, this_cnv_info.depth_, cnv.get_ref_start(), pad, min, max, color,0.7);
    multi_seg.push_back(this_segs);

    // exon reads depth ratio
    for(int i=0;i<this_var.genomic_pos.size();++i){
        vector<Segment> rd_segs=RDratioLane(this_var.genomic_pos, i, cnv.cnv_info_.at(0).depth_, cnv.cnv_info_.at(cnv_size-1).depth_,cnv.get_ref_start(),pad,"black");
        multi_seg.push_back(rd_segs);
    }

    std::cerr<<"[LanePlot] Sample Depth Coord..."<<std::endl;
    // Insert Coverage Coord
    //multi_seg_depth.emplace(multi_seg_depth.begin(), CoordLane(start, end, "Normalized Depth", pad, "black", float(min), float(max)));
    multi_seg.push_back(CoordLane(start, end, "Normalized Depth", pad, "black", float(min), float(max)));
    this_laneplot.push_back(Lane(opts.lane*15, multi_seg));

    std::cerr<<"[LanePlot] Sample Variant Information..."<<std::endl;
    int laneplot_size = this_laneplot.size();
    if( !cnv.gt_info_.empty() ){
        int cnv_start = this_var.genomic_pos.at(this_var.start-1).start_;
        int cnv_end = this_var.genomic_pos.at(this_var.end-1).end_;
        for (int idx = 0; idx < cnv.gt_info_.size(); ++idx) {
            //std::cerr<<"pos: "<<std::to_string(cnv.gt_info_[idx].pos_)<<", infor: "<<cnv.gt_info_[idx].infor_<<std::endl;
            vector<Feature> features;
            std::vector<Segment> seg_vec;
            if( cnv.gt_info_[idx].pos_ >= cnv_start && cnv.gt_info_[idx].pos_ <= cnv_end ){
                features.emplace_back(FeatureType::kDot, 50, "orange", "orange", 1.0);// rate == 0.5 for black
            }else{
                features.emplace_back(FeatureType::kDot, 50, "orange", "orange", 0.5); // rate == 1.0 for red
            }
            //Feature this_feature(FeatureType::kDot, 50, "orange", "orange",0.5);
            //this_feature.label = cnv.gt_info_[idx].infor_;
            //features.push_back(this_feature);
            features.back().label = cnv.gt_info_[idx].infor_;
            Segment this_segment(pad+(cnv.gt_info_[idx].pos_ - start), features, 0.7);

            bool flag = true;
            int dist = int ((end - start)/7);
            for(int i=laneplot_size;i<this_laneplot.size();++i){
                if( (this_laneplot[i].end + dist) < this_segment.start){
                    this_laneplot[i].multi_segments.push_back({this_segment});
                    this_laneplot[i].end = this_segment.end;
                    flag = false;
                    break;
                }
            }
            seg_vec.push_back(this_segment);
            if( flag ) this_laneplot.push_back(Lane(opts.lane,{seg_vec}));
        }
    }

    //Copy-number Analysis by BAse-level NormAlization (CABANA)
    std::cerr<<"[LanePlot] CABANA..."<<std::endl;
    //std::vector<std::vector<Segment>> multi_seg_cabana;
    multi_seg.clear();
    for(int i=0;i<(cnv_size-1);i++){
        if( i == 0 ){// just the mean depth of the reference set.
            color="red";
        }else{
            color="grey";
        }
        //std::cerr<<"Index:"<<std::to_string(i)<<" CABANALane"<<std::endl;
        vector<Segment> cabana_segs = CABANALane(this_var.genomic_pos, cnv.get_ref_start(), cnv.cnv_info_.at(i).depth_, cnv.cnv_info_.at(cnv_size-1).depth_, pad, color,1.0);
        multi_seg.push_back(cabana_segs);
    }

    // Insert CABANA Coord
    int exon_length=0;
    int total_size = pad*2 + this_var_size;
    float large_factor = 1.0*total_size / this_var_exon_size;
    for (int i = 0; i< this_var.genomic_pos.size(); ++i) {
        int this_exon_size = int(0.5+large_factor*(this_var.genomic_pos[i].end_ - this_var.genomic_pos[i].start_ + 1));
        vector<Feature> features;
        std::vector<Segment> seg_vec;
        features.emplace_back(FeatureType::kVerticalLine, 20, "blue", "blue", 1.0); // dashed line
        //features.emplace_back(FeatureType::kVerticalLine, 20, "black", "black", 0.5); // dashed line
        Segment this_segment(exon_length+this_exon_size, features, 0.6);
        seg_vec.push_back(this_segment);
        multi_seg.push_back(seg_vec);
        exon_length += this_exon_size;
    }
    multi_seg.push_back(CoordLane(1, total_size, "Exon normalized Depth", 0, "black", -1.0, 1.0));
    this_laneplot.push_back(Lane(opts.lane*15, multi_seg));

    //features.emplace_back(FeatureType::kLabel, 40, color, color, 0.9);
    //this_laneplot[i].multi_segments.push_back({this_segment});
    multi_seg.clear();
    exon_length=0;
    for (int i = 0; i< this_var.genomic_pos.size(); ++i) {
        vector<Feature> features;
        int this_exon_size = int(0.5+large_factor*(this_var.genomic_pos[i].end_ - this_var.genomic_pos[i].start_+1));
        string label = "EX1";
        string color="black";
        if( this_var.strand == "+" ){
            label = "EX"+std::to_string(i+1);

        }else{
            label = "EX"+std::to_string(this_var.genomic_pos.size()-i);
        }
        if( (i+1) >= this_var.start && (i+1) <= this_var.end ){
            color="red";
        }
        features.emplace_back(FeatureType::kLabel, 40, color, color, 0.9,-45, 1.0);
        features.back().label = label;
        Segment this_segment(exon_length + this_exon_size/2, features, 0.6);
        multi_seg.push_back({this_segment});
        exon_length += this_exon_size;
    }
    this_laneplot.push_back(Lane(opts.lane,multi_seg));


    return this_laneplot;
}
}