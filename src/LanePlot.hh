// varViewer
// Date: 2023-08-09
// Author: Huang Fei <huangfei@genomics.cn>
//

#pragma once

#include <iostream>
#include <list>
#include <memory>
#include <string>
#include <unordered_map>
#include <vector>
#include <boost/optional.hpp>
#include "utils.hh"
#include "Coverage.hh"
#include "Aligns.hh"
#include "Variant.hh"
#include "CNV.hh"

namespace varViewer
{
enum class FeatureType
{
    kRect, // Bar
    kRectWithLeftBreak,
    kRectWithRightBreak,
    kLine,
    kLineSegment,
    kArrows,
    kVerticalLine,
    kLabel,
    kCoord,
    kPath,  // line with path
    kDot,  // dot
};

struct Feature
{
    Feature(FeatureType type, int length, std::string fill, std::string stroke)
        : type(type)
        , length(length)
        , fill(std::move(fill))
        , stroke(std::move(stroke))
    {
    }
    Feature(FeatureType type, int length, std::string fill, std::string stroke, float rate_)
        : type(type)
        , length(length)
        , fill(std::move(fill))
        , stroke(std::move(stroke))
        ,rate(rate_)
    {
    }
    Feature(FeatureType type, int length, std::string fill, std::string stroke, float rate_,float min_, float max_)
        : type(type)
        , length(length)
        , fill(std::move(fill))
        , stroke(std::move(stroke))
        ,rate(rate_)
        ,min(min_)
        ,max(max_)
    {
    }
    ~Feature(){
        if (!depth.empty()){
            std::vector<float>().swap(depth);
        }
    }
    FeatureType type;
    int length;
    boost::optional<std::string> label;
    std::vector<float> depth; // for kPath
    std::string fill;
    std::string stroke;
    float rate = 0.5; // rate of the height from the top, 1 for drawVerticalLine dashedline
    float min = 0.0; // Coverage min depth, or 1 for arrows left direct, kLabel rotate
    float max = 1.0; // Coverage max depth, or 1 for arrows right direct, kLabel align direct, 0 for start, 1,middle
};

//Single or Pair-End reads in the same segment
struct Segment
{
    Segment(int start, std::vector<Feature> features, double opacity)
        : start(start)
        , features(std::move(features))
        , opacity(opacity)
    {
        end = start;
        for (const auto& feature : this->features) {
            end += feature.length;
        }
    }
    ~Segment(){
        if (!features.empty()){
            std::vector<Feature>().swap(features);
        }
    }
    int start;
    int end;
    std::vector<Feature> features;
    double opacity;
};

// Different Pair-End reads in different segment
struct Lane
{
    Lane(float height, std::vector<std::vector<Segment>> segments)
        : height(height)
        , multi_segments(segments)
    {
        start = multi_segments[0].front().start;
        end = multi_segments[0].back().end;
    }
    ~Lane(){
        if ( !multi_segments.empty() ){
            for(int i=0;i<multi_segments.size();++i ){
                std::vector<Segment>().swap(multi_segments[i]);
            }
        }
    }
    float height;
    int start;
    int end;
    //std::vector<Segment> segments; just for one type draw
    std::vector<std::vector<Segment>> multi_segments;
};
using LanePlot = std::vector<Lane>;
LanePlot generateLanePlot(Aligns & align, opts_s& opts, Coverage& cover);
LanePlot generateLanePlot(CNV & cnv, opts_s & opts, int pad);

}