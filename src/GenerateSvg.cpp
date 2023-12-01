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

#include "GenerateSvg.hh"
#include <fstream>
#include <list>
#include <stdexcept>
#include <iostream>

using boost::optional;
using std::list;
using std::ofstream;
using std::string;
using std::to_string;
using std::vector;

namespace varViewer
{

static void drawRect(ofstream& out, float x, float y, float width, float height, const string& fill, const string& stroke, double opacity, float rate=1.0) {
    out << "<rect x=\"" << x << "\" y=\"" << (y+ height*(1-rate)) << "\"";
    out << " width=\"" << width << "\" height=\"" << height*rate << "\"";
    out << " stroke=\"" << stroke << "\"";
    out << " fill=\"" << fill << "\"";
    out << " opacity=\"" << opacity << "\"";
    out << " />\n";
}

static void drawRectWithLeftBreak(ofstream& out, float x, float y, float width, float height, const string& fill, const string& stroke) {
    out << "<path d=\"";
    out << "M " << x << " " << y << " ";
    out << "h " << width << " ";
    out << "v " << height << " ";
    out << "h -" << width << " ";
    out << "l " << 5 << " " << -height / 4 << " ";
    out << "l " << -5 << " " << -height / 4 << " ";
    out << "l " << 5 << " " << -height / 4 << " ";
    out << "Z\"";

    out << " fill=\"" << fill << "\"";
    out << " stroke=\"" << stroke << "\"";

    out << "/>\n";
}

static void drawRectWithRightBreak(ofstream& out, float x, float y, float width, float height, const string& fill, const string& stroke) {
    out << "<path d=\"";
    out << "M " << x + width << " " << y << " ";
    out << "h " << -width << " ";
    out << "v " << height << " ";
    out << "h " << width << " ";
    out << "l " << -5 << " " << -height / 4 << " ";
    out << "l " << 5 << " " << -height / 4 << " ";
    out << "l " << -5 << " " << -height / 4 << " ";
    out << "Z\"";

    out << " fill=\"" << fill << "\"";
    out << " stroke=\"" << stroke << "\"";

    out << "/>\n";
}

static void drawLine(ofstream& out, float x, float y, float width, float height, const string& stroke, float rate=0.5, float opacity=0.4, bool dash=false ) {
    int inter=4;
    out << "<line ";
    out << "x1=\"" << x << "\" y1=\"" << y +  height * (1-rate)<< "\" ";
    out << "x2=\"" << x + width << "\" y2=\"" << y + height * (1-rate) << "\" ";
    if( dash ){
        out<< "style=\"stroke:"<<stroke<<"; stroke-dasharray:"<<inter<<" "<<inter<<"; opacity:"<<opacity<<";\"";
    }else{
        out << "stroke=\"" << stroke << "\" ";
    }
    out << "/>\n";
}

static void drawLineSegment(ofstream& out, float x, float y, float width, float height, const string& stroke) {
    out << "<line ";
    out << "x1=\"" << x << "\" y1=\"" << (y + height) << "\" ";
    out << "x2=\"" << x + width << "\" y2=\"" << (y + height) << "\" ";
    out << "stroke=\"" << stroke << "\" ";
    out << "/>\n";
    out << "<line ";
    out << "x1=\"" << x + width << "\" y1=\"" << (y + height  - 5) << "\" ";
    out << "x2=\"" << x + width << "\" y2=\"" << (y + height) << "\" ";
    out << "stroke=\"" << stroke << "\" ";
    out<<"/>\n";
}

static void drawLetter(ofstream& out, float x, float y, float width, float height, char letter,float fontsize=11,std::string type="SNV",float rate=0.5) {
    string color = "black";
    if (letter == 'A' || letter == 'a') {
        color = "#FF6347";
    } else if (letter == 'T' || letter =='t') {
        color = "#FCA100";
    } else if (letter == 'C' || letter =='c') {
        color = "#393939";
    } else if (letter == 'G' || letter =='g' ) {
        color = "#2F8734";
    }
    if(type == "SNV" || type == "SV"){
        out << "<text x=\"" << x + width / 2 << "\" y=\"" << (y + height*(1.0-rate)) << "\"";
        out << " dy=\"0.25em\""; // dominant-baseline="middle"
        out << " text-anchor=\"middle\"";
        out << " font-family=\"monospace\"";
        out << " font-size=\""<< fontsize <<"px\"";
        out << " stroke=\"" << color << "\"";
        out << ">";
        out << letter << "</text>";
        out <<"\n";
    }else{ // color line
         out << "<line ";
        out << "x1=\"" << x << "\" y1=\"" << (y + height *(1-rate)) << "\" ";
        out << "x2=\"" << x + width << "\" y2=\"" << (y + height*(1-rate)) << "\" ";
        out << "stroke=\"" << color << "\" ";
        out << "/>\n";
    }
}

static void drawLabel(ofstream& out, float x, float y, float width, float height, const string& text, float fontsize=11, float rate=0.5, std::string align="start", int rotate=0,std::string color="black") {
    //string color = "#393939";
    out << "<text x=\"" << (x-width/2) << "\" y=\"" << y + height*(1.0-rate) << "\"";
    out << " dy=\"0.25em\""; // dominant-baseline="middle"
    //if( rotate != 0 ){
    //    out << " text-anchor=\"middle\"";
    //}else{
    //    out << " text-anchor=\"start\"";
    //}
    out <<" text-anchor=\""<<align<<"\"";
    out << " font-family=\"monospace\"";
    out << " font-size=\""<<fontsize<<"\"";
    out << " stroke=\"" << color << "\"";
    if( rotate != 0) out << " transform=\"rotate(" <<std::to_string(rotate)<<","<<std::to_string(x-width/2)<<","<<std::to_string(y+height*(1.0-rate))<<")\"";
    out << ">";
    out << text << "</text>";
    out <<"\n";

    //<text x="50" y="24" font-size="28" font-family="Times New Roman" fill="black" transform="rotate(angel, x, y)">love life</text>
}

static void drawText(ofstream& out, float x, float y, float width, float height, const string& text, float fontsize, std::string type="SNV") {
    const float letterWidth = width / text.length();
    for (int letterIndex = 0; letterIndex != text.length(); ++letterIndex) {
        if (text[letterIndex]==' ') continue;
        drawLetter(out, x + letterWidth * letterIndex, y, letterWidth, height, text[letterIndex], (fontsize+2.0), type);
    }
}

static void
drawArrows(ofstream& out, float x, float y, float width, float height, const string& stroke, const string& direct,  float rate, const optional<string>& text) {
    // rate = 0.5
    out << "<line x1=\"" << x << "\" y1=\"" << (y + height *(1-rate)) << "\"";
    out << " x2=\"" << x + width << "\" y2=\"" << (y + height *(1-rate)) << "\"";
    out << " stroke=\"" << stroke << "\"";
    if( direct == "left" ){
            out << " marker-start=\"url(#arrow)\" />\n";
    }else if( direct == "right" ){
            out << " marker-end=\"url(#arrow)\" />\n";
    }else{
        out << " marker-start=\"url(#arrow)\" marker-end=\"url(#arrow)\" />\n";
    }

    if (text) {
        out << "<text x=\"" << x + width / 2 << "\" y=\"" << y + height*(1-rate) << "\"";
        out << " dy=\"0.25em\"";
        out << " text-anchor=\"middle\" font-family=\"monospace\" font-size=\"13px\"";
        out << " style=\"stroke:white; stroke-width:1.0em\" >";
        out << *text << "</text>\n";

        out << "<text x=\"" << x + width / 2 << "\" y=\"" << y + height*(1-rate) << "\"";
        out << " dy=\"0.25em\"";
        out << " text-anchor=\"middle\" font-family=\"monospace\" font-size=\"13px\"";
        out << " font-weight=\"lighter\" stroke=\"" << stroke << "\" >";
        out << *text << "</text>\n";
    }
}

static void drawVerticalLine(ofstream& out, float x, float y, float width, float height, const string& stroke, bool stroke_width=true, bool dash=false) {
    out << "<line ";
    out << "x1=\"" << x << "\" y1=\"" << y << "\" ";
    out << "x2=\"" << x << "\" y2=\"" << y + height << "\" ";
    if( dash ){
        if( stroke_width ){
            out<< "style=\"stroke:"<<stroke<<"; stroke-dasharray:4 4; opacity: 0.4; strok-width:3px;\"";
        }else{
            out<< "style=\"stroke:"<<stroke<<"; stroke-dasharray:4 4; opacity: 0.4; \"";
        }
    }else{
        if (stroke_width){
            out << " style=\"stroke:" << stroke << "; stroke-width:3px\"";
        }else{
            out << "stroke=\"" << stroke << "\" ";
        }
    }
    out << "/>\n";
}
static void drawDot(ofstream& out, float x, float y, float width, float height, const string& fill, const string& stroke){
    float r = (width > height ) ? width/2 : height/2;
    out<<"<circle ";
    out<<" cx=\"" << x << "\" cy=\"" << y << "\" ";
    out<<" r=\"" << r <<"\" ";
    out<<" fill=\"" << fill <<"\" ";
    out<<" stroke=\""<< stroke<<"\" ";
    out<<"></circle>\n";
}

static void drawPath(ofstream& out, float x, float y, float width, float height, const string& fill, const string& stroke, const std::vector<float> depth,int min=0,int max=1, float opacity=0.6) {
    float x_bin = width / depth.size();
    float y_bin = height / (max-min);

    out << "<path d=\"";
    out << "M " << x << " " << y+ height << " L ";
    //for(int i=0;i<depth.size();++i){
    //    float x_index = x + x_bin*i;
    //    float y_index = y + (max - depth[i])*y_bin;
    //    out<<" "<< x_index <<" "<< y_index;
    //}
    float prev_x = x;
    float prev_y = y+(max - depth[0])*y_bin;
    out << " "<< prev_x <<" "<<prev_y;
    for(int i=1;i<depth.size();++i){
        float x_index = x + x_bin*i;
        float y_index = y + (max - depth[i])*y_bin;
        if( y_index != prev_y ){
            if( x_index > prev_x ){
                out<<" "<< prev_x <<" "<< prev_y;
            }
            out<<" "<< x_index <<" "<< y_index;
        }
        prev_x = x_index;
        prev_y = y_index;
    }
    out<<" "<< prev_x <<" "<< prev_y;

    float x_last = x + x_bin * (depth.size() -1);
    out<<" "<< x_last <<" "<< y+height;
    out << "\"";
    out << " fill=\"" << fill << "\"";
    out << " stroke=\"" << stroke << "\"";
    out << " opacity=\"" << opacity << "\"";
    out << "/>\n";
}

static void drawCoord(ofstream& out, float x, float y, float width, float height, string label="xxx", int ymin=0, int ymax=1, float fontsize=8){
    // x-axis
    drawLine(out, x, y, width, height,"black", 0.0);
    // y-axis
    drawVerticalLine(out, x, y, width, height, "black",false);
    // y-label
    if( ymin == -1 && ymax == 1){
        drawLabel(out, x-10, y+height-fontsize/2, 15, fontsize+2, "-1.0", fontsize);
        drawLabel(out,x-10, y-fontsize/2, 15, fontsize+2, "1.0", fontsize);
    }else{
        drawLabel(out, x-10, y+height-fontsize/2, 15, fontsize+2, std::to_string(ymin), fontsize);
        drawLabel(out,x-10, y-fontsize/2, 15, fontsize+2, std::to_string(ymax), fontsize);
    }

    // draw five part
    float per = (ymax-ymin)/4.0;
    for(int i=1;i<4;i++){
        //float y_coord = y+height * (1 - i/4);
        drawLine(out, x, y, width, height, "grey", i/4.0, 0.4, true);
        if( i == 2 ){
            float value = ymin+per*i;
            if( ymax > 1.0 ){
                std::string yvalue = "97.5";
                if( value != 97.5 ){
                    yvalue = std::to_string(int(value+0.5));
                }
                drawLabel(out, x-10, y+height*i/4.0-fontsize/2, 15, fontsize+2, yvalue, fontsize);
            }else{
                if( ymin == -1 && ymax == 1 ){
                    drawLabel(out, x-10, y+height*i/4.0-fontsize/2, 15, fontsize+2, "0.0", fontsize);
                }else{
                    drawLabel(out, x-10, y+height*i/4.0-fontsize/2, 15, fontsize+2, std::to_string(value), fontsize);
                }
            }
        }
    }
    if( label !="xxx"){
        drawLabel(out, x-18, y+height/2, 15, fontsize+2, label, fontsize+2, 0.5, "middle", -90);
    }
}

void drawLane(ofstream& out, float baseWidth, int xPosStart, float yPos, std::string type, const Lane& lane) {
    for (const auto& multi_segment: lane.multi_segments){
        for (const auto& segment : multi_segment) {
            float xPos = xPosStart + segment.start * baseWidth;
            for (const auto& feature : segment.features) {
                const float featureWidth = feature.length * baseWidth;
                    if (feature.type == FeatureType::kRect) {
                    drawRect(out, xPos, yPos, featureWidth, lane.height, feature.fill, feature.stroke, segment.opacity, feature.rate);
                    if (feature.label){
                        drawText(out, xPos, yPos, featureWidth, lane.height, *feature.label, baseWidth, type);
                    }
                } else if (feature.type == FeatureType::kLabel) {
                    if(feature.label){
                        if( feature.max == 0.0 ){
                            drawLabel(out, xPos+featureWidth/2, yPos, featureWidth, lane.height, *feature.label, 11,  feature.rate, "start", int(feature.min), feature.stroke);
                        }else{
                            drawLabel(out, xPos+featureWidth/2, yPos, featureWidth, lane.height, *feature.label, 11,  feature.rate, "middle", int(feature.min), feature.stroke);
                        }
                    }
                } else if (feature.type == FeatureType::kRectWithLeftBreak) {
                    drawRectWithLeftBreak(out, xPos, yPos, featureWidth, lane.height, feature.fill, feature.stroke);
                    if (feature.label) {
                        drawText(out, xPos, yPos, featureWidth, lane.height, *feature.label, baseWidth);
                    }
                } else if (feature.type == FeatureType::kRectWithRightBreak) {
                    drawRectWithRightBreak(out, xPos, yPos, featureWidth, lane.height, feature.fill, feature.stroke);
                    if (feature.label){
                        drawText(out, xPos, yPos, featureWidth, lane.height, *feature.label,baseWidth);
                    }
                } else if (feature.type == FeatureType::kLine) {
                    drawLine(out, xPos, yPos, featureWidth, lane.height, feature.stroke,feature.rate);
                } else if (feature.type == FeatureType::kArrows) {
                    std::string direct="right";
                    if( feature.min == 1.0 && feature.max == 1.0 ) {
                        direct = "both";
                    }else if( feature.min == 1.0 ){
                        direct = "left";
                    }
                    drawArrows(out, xPos, yPos, featureWidth, lane.height, feature.stroke, direct, feature.rate, feature.label);
                } else if (feature.type == FeatureType::kVerticalLine) {
                    if( feature.rate == 1.0 ){
                        drawVerticalLine(out, xPos, yPos, featureWidth, lane.height, feature.stroke, true, true);
                    }else{
                        drawVerticalLine(out, xPos, yPos, featureWidth, lane.height, feature.stroke);
                    }
                }else if( feature.type == FeatureType::kLineSegment){
                    drawLineSegment(out, xPos, yPos, featureWidth, lane.height, feature.stroke);
                    if(feature.label){
                        float fontsize=10;
                        float label_length = 4*fontsize;
                        drawLabel(out, xPos+featureWidth, yPos-15, label_length, 10, *feature.label);
                    }
                }else if( feature.type == FeatureType::kCoord ){
                    drawCoord(out, xPos, yPos, featureWidth, lane.height, *feature.label, feature.min, feature.max);
                }else if( feature.type == FeatureType::kPath ){
                    drawPath(out, xPos, yPos, featureWidth, lane.height, feature.fill, feature.stroke, feature.depth, int(feature.min), int(feature.max), segment.opacity);
                }else if( feature.type == FeatureType::kDot ){
                    drawDot(out, xPos, yPos, featureWidth, lane.height, feature.fill, feature.stroke);
                    if (feature.label){
                        if( feature.rate == 1.0 ){
                            drawLabel(out, xPos-20, yPos+10, 60, 10, *feature.label, 9, 0.5, "start", 0, "red");
                        }else{
                            drawLabel(out, xPos-20, yPos+10, 60, 10, *feature.label, 9);
                        }
                    }
                } else {
                    throw std::runtime_error("Encountered feature of unknown type");
                }
                xPos += featureWidth;
            }
        }
    }
}

static void printHeader(int width, int height, ofstream& out) {
    //svgFile << "<svg width=\"100%\" height=\"100%\" xmlns=\"http://www.w3.org/2000/svg\" viewBox=\"0 0 "<< plotWidth <<" "<< plotHeight<<"\">\n";
    out << "<svg width=\"" << width << "\" height=\"" << height << "\" xmlns=\"http://www.w3.org/2000/svg\">\n";
    out << "<defs>\n"
        "    <linearGradient id=\"BlueWhiteBlue\" x1=\"0%\" y1=\"0%\" x2=\"0%\" y2=\"100%\">\n"
        "      <stop offset=\"0%\" style=\"stop-color:#8da0cb;stop-opacity:0.8\" />\n"
        "      <stop offset=\"50%\" style=\"stop-color:#8da0cb;stop-opacity:0.1\" />\n"
        "      <stop offset=\"100%\" style=\"stop-color:#8da0cb;stop-opacity:0.8\" />\n"
        "    </linearGradient>\n"
        "    <linearGradient id=\"OrangeWhiteOrange\" x1=\"0%\" y1=\"0%\" x2=\"0%\" y2=\"100%\">\n"
        "     <stop offset=\"0%\" style=\"stop-color:#fc8d62;stop-opacity:0.8\" />\n"
        "      <stop offset=\"50%\" style=\"stop-color:#fc8d62;stop-opacity:0.1\" />\n"
        "      <stop offset=\"100%\" style=\"stop-color:#fc8d62;stop-opacity:0.8\" />\n"
        "    </linearGradient>\n"
        "    <linearGradient id=\"GreenWhiteGreen\" x1=\"0%\" y1=\"0%\" x2=\"0%\" y2=\"100%\">\n"
        "      <stop offset=\"0%\" style=\"stop-color:#66c2a5;stop-opacity:0.8\" />\n"
        "      <stop offset=\"50%\" style=\"stop-color:#66c2a5;stop-opacity:0.1\" />\n"
        "      <stop offset=\"100%\" style=\"stop-color:#66c2a5;stop-opacity:0.8\" />\n"
        "    </linearGradient>\n"
        "    <marker id=\"arrow\" viewBox=\"0 0 10 10\" refX=\"9\" refY=\"5\"\n"
        "        markerWidth=\"6\" markerHeight=\"6\"\n"
        "        orient=\"auto-start-reverse\">\n"
        "      <path d=\"M 0 0 L 10 5 L 0 10 z\" />\n"
        "    </marker>"
        "</defs>\n";
}

void generateSvg(Aligns & alignReads, const string& outputPath, Coverage& cover, opts_s& opts) {
    const int kPlotPadX = 20;
    const int kPlotPadY = 20;
    const int coverHeight = 30;
    float kSpacingBetweenLanes = opts.lane/3.0;
    //if( kSpacingBetweenLanes < 1 ) kSpacingBetweenLanes = 1;

    std::cerr<<"Start to SVG"<<std::endl;
    LanePlot lane_plot = generateLanePlot(alignReads, opts, cover);
    int read_num = lane_plot.size();
    int base_num = alignReads.get_end_offset() - alignReads.get_start_offset();
    //const int plotWidth = kBaseWidth * (base_num+4) +  2 * kPlotPadX;
    float kBaseWidth = (opts.width - 2.0*kPlotPadX)/(base_num+4);
    if (strcmp(opts.type, "SV")==0){
        kBaseWidth = 6;
        opts.width = kBaseWidth*(base_num+4) + 2 * kPlotPadX;
    }
    //plot height for SNV
    //std::cerr<<"Read num:"<<std::to_string(read_num)<<",lane and space:"<<std::to_string(opts.lane+kSpacingBetweenLanes)<<std::endl;
    int plotHeight = int(0.5+(1+read_num) * (opts.lane+kSpacingBetweenLanes) + 3 * kPlotPadY);
    //std::cerr<<"Read num:"<<std::to_string(read_num)<<",plotHeight:"<<std::to_string(plotHeight)<<std::endl;
    float yPos = 2*kPlotPadY;
    ofstream svgFile(outputPath);
    if (svgFile.is_open()) {
        printHeader(opts.width, plotHeight, svgFile);
        for(const auto& this_lane: lane_plot){
            drawLane(svgFile, kBaseWidth, kPlotPadX, yPos, opts.type, this_lane);
            yPos += this_lane.height + kSpacingBetweenLanes;
        }

        if(strcmp(opts.type, "SNV") ==0 ){
            int leftflank_length = alignReads.get_start() - (alignReads.get_refer_start() + alignReads.get_start_offset()) -1;
            int var_length = alignReads.get_end() - alignReads.get_start() + 1;
            drawRect(svgFile, kPlotPadX+kBaseWidth*leftflank_length, 2*kPlotPadY, var_length*kBaseWidth, plotHeight-3*kPlotPadY, "#E5E5E5", "none", 0.5);

            // add variant information
            std::string mutant_info = alignReads.get_mutant_info();
            drawLabel(svgFile, kPlotPadX+kBaseWidth*mutant_info.size()/2, kPlotPadY, kBaseWidth*mutant_info.size(), kPlotPadY, mutant_info, 12);
        }
        svgFile << "</svg>" << std::endl;
        std::vector<Lane>().swap(lane_plot);
    } else {
        throw std::runtime_error("Unable to open " + outputPath);
    }
}
void generateSvg(CNV& cnv, string outputPath, opts_s& opts) {
    const int kPlotPadX = 50;
    const int kPlotPadY = 20;
    //const int coverHeight = 30;
    float kSpacingBetweenLanes = opts.lane;
    if( kSpacingBetweenLanes < 20 ) kSpacingBetweenLanes = 20;

    std::cerr<<"[gSvg] Start to CNV SVG"<<std::endl;
    LanePlot lane_plot = generateLanePlot(cnv, opts, opts.flank);
    int total_lane_height = 0;
    for(int i=0;i<lane_plot.size(); ++i){
        //std::cerr<<"Index:"<<std::to_string(i)<<",heigth:"<<std::to_string(lane_plot[i].height)<<std::endl;
        total_lane_height += lane_plot[i].height;
    }
    int base_num = cnv.get_cnv_size();
    //const int plotWidth = kBaseWidth * (base_num+4) +  2 * kPlotPadX;
    float kBaseWidth = (opts.width - 2.0*kPlotPadX)/(base_num+4);
    //std::cerr<<"total lane height:"<<std::to_string(total_lane_height)<<",size:"<<std::to_string(lane_plot.size())<<",spacing:"<<std::to_string(kSpacingBetweenLanes);
    int plotHeight = (total_lane_height + int((lane_plot.size()-1)*kSpacingBetweenLanes) + 3 * kPlotPadY + opts.lane);
    //std::cerr<<",plotHeight:"<<std::to_string(plotHeight)<<std::endl;
    float yPos = 2*kPlotPadY;

    VariantSpecific this_var = cnv.get_VariantSpecific();
    const string outfile = outputPath+"_"+this_var.infor+".svg";
    ofstream svgFile(outfile);
    if (svgFile.is_open()) {
        printHeader(opts.width+20, plotHeight, svgFile);
        int index = 0;
        for(const auto& this_lane: lane_plot){
            //std::cerr<<"yPos:"<<std::to_string(yPos)<<std::endl;
            drawLane(svgFile, kBaseWidth, kPlotPadX, yPos, opts.type, this_lane);
            yPos += this_lane.height + kSpacingBetweenLanes;
            if( index == (lane_plot.size()-3) ) yPos += opts.lane;
            index++;
        }
        drawLabel(svgFile, kPlotPadX, kPlotPadY-12, kBaseWidth*40, kPlotPadY, this_var.infor, 12);

        int extend = 40;
        int var_before = opts.flank + (this_var.genomic_pos.at(this_var.start-1).start_ - this_var.genomic_pos.at(0).start_)-extend;
        int var_length = (this_var.genomic_pos.at(this_var.end-1).end_ - this_var.genomic_pos.at(this_var.start-1).start_)+2*extend;
        drawRect(svgFile, kPlotPadX+kBaseWidth*var_before, int(2.5*kPlotPadY), var_length*kBaseWidth, plotHeight-(2*kPlotPadY+17*opts.lane+2.5*kSpacingBetweenLanes), "#E5E5E5", "none", 0.5);
        svgFile << "</svg>" << std::endl;
    } else {
        throw std::runtime_error("Unable to open " + outputPath);
    }
}
}

