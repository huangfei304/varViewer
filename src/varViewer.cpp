#include <algorithm>
#include <cstring>
#include <fstream>
#include <getopt.h>
#include <iostream>
#include <sstream>

#include "Reference.hh"
#include "Variant.hh"
#include "Aligns.hh"
#include "GenerateSvg.hh"
#include "htslib/faidx.h"
#include "htslib/sam.h"
#include "htslib/hts.h"
#include "Coverage.hh"
#include "utils.hh"
#include "Depth.hh"

namespace varViewer {
#define VARVIEWER_NAME    "varViewer"
#define VARVIEWER_VERSION "0.3"

static const char *help = "\
Required arguments:\n\
    -r, --ref REF           Reference genome in FASTA format\n\
    -b, --bam IN            Sorted and indexed input BAM or CRAM file\n\
    -i, --in VAR            Variant information file\n\
    -t, --type STR          Variant Type for analysis: SNV, CNV or SV\n\
    -p, --prefix PFX        Output file prefix\n\
\n\
Options:\n\
    -n, --list STR          bam list file for CNV reference, one bam per line\n\
    -k, --map  STR          mapability file for CNV, gz format\n\
    -v, --var  STR          variant genotype information for CNV, format[chr,start,ref,alt,info]\n\
    -a, --bed  STR           transcript (gene) file( BED format) for CNV analysis\n\n\
    -w, --width INT          plot width, [1300]\n\
    -l, --lane  FLT          height per lane [10.0]\n\n\
    -s, --mis   FLT          Maxmium mismatch length rate in softclip [0.1]\n\
    -c, --cov   INT          Maxmium reads coverage for CNV [5000]\n\
    //-e, --per   INT          plot point size for CNV [1000]\n\
    -m, --mqual INT          Minimum read mapping quality for filtering,[20]\n\
    -q, --bqual INT          Minimum base quality for lower case [20]\n\
    -f, --flank INT          Maxmium Flank bases of the variant to plot, must be more than read length. [500]\n\
    -d, --dup                Remove duplication reads\n\
    -h, --help               Show this help\n\
";



void usage(){
    std::cerr << VARVIEWER_NAME << " v" << VARVIEWER_VERSION << std::endl
              << std::endl
              << help << std::endl;
}

enum input_format { INPUT_SAM = 0, INPUT_BAM, INPUT_CRAM };
#define MODE_LEN 3
#define FMT_LEN  5

int varViewer_main(int argc, char **argv)
{
    opts_s opts;

    // bam parameters
    char mode[MODE_LEN], fmt[FMT_LEN];
    enum input_format format;
    samFile *input_sf=NULL;
    bam_hdr_t *hdr=NULL;
    hts_itr_t *iter=NULL;
    hts_idx_t *index=NULL;
    faidx_t *_fai=NULL;

    memset(mode, '\0', MODE_LEN);
    memset(fmt, '\0', FMT_LEN);

    int c, optidx;
    static struct option long_options[] = {
        {"ref", 1, nullptr, 0},                     // r
        {"bam", 1, nullptr, 0},                     // b
        {"in", 1, nullptr, 0},                      // i
        {"type", 1, nullptr, 0},                    // t
        {"list", 1, nullptr, 0},                    // n
        {"prefix", 1, nullptr, 0},                  // p
        {"bed", 1, nullptr, 0},                    // a
        {"map", 1, nullptr, 0},                    // k
        {"var",1, nullptr, 0},                     // v
        {"width", 1, nullptr, 0},                   // w
        {"lane", 1, nullptr, 0},                    // l
        {"cov", 1, nullptr, 0},                     // c
        {"per", 1, nullptr, 0},                     // e
        {"mqual", 1, nullptr, 0},                   // m
        {"bqual", 1, nullptr, 0},                   // q
        {"flank", 1, nullptr, 0},                   // f
        {"dup", 0, nullptr,  0},                    //d
        {"help", 0, nullptr, 0},                    // h
        // hidden options
        //{"color", 0, nullptr, 0},                 // Z
        {nullptr, 0, nullptr, 0}};
    const char *short_options = "r:b:i:t:n:p:a:k:v:w:l:c:e:m:q:f:dh";
    const char *shorter_options = "rbitnpakvwlcemqfdh";

    if (argc == 1) {
        usage();
        return EXIT_SUCCESS;
    }

    while ((c = getopt_long(argc, argv, short_options, long_options, &optidx)) != -1) {
        if (c == 0) {
            c = shorter_options[optidx];
        }

        switch (c) {
        case 'r': // ref
            opts.ref = optarg;
            break;
        case 'b': // bam
            opts.bam = optarg;
            break;
        case 'i': // infile
            opts.var = optarg;
            break;
        case 't': // type
            opts.type = optarg;
            break;
        case 'n': // bam list
            opts.list = optarg;
            break;
        case 'p': // prefix
            opts.prefix = optarg;
            break;
        case 'a':
            opts.bed = optarg;
            break;
        case 'k':
            opts.map = optarg;
            break;
        case 'v': // variant information
            opts.gt = optarg;
            break;
        case 'w': // plot width
            opts.width = std::atoi(optarg);
            break;
        case 'l': // lane height
            opts.lane = std::atof(optarg);
            break;
        case 'c': //max cover
            opts.cov = std::atoi(optarg);
            break;
        case 'e': //per
            opts.per = std::atoi(optarg);
            break;
        case 'm': // mapping quality
            opts.mqual = std::atoi(optarg);
            break;
        case 'q': // base quality
            opts.bqual = std::atoi(optarg);
            break;
        case 'f': // flank
            opts.flank = std::atoi(optarg);
            break;
        case 'd':
            opts.dedup = true;
            break;
        case 'h': // help
            usage();
            return EXIT_SUCCESS;
        default:
            return EXIT_FAILURE;
        }
    }

    if (opts.bam != nullptr) {
        fprintf(stderr,"Input alignment file: %s\n", opts.bam);
    } else {
        fprintf(stderr, "No input alignment file given\n");
        exit(EXIT_FAILURE);
    }
    if (opts.ref != nullptr) {
        fprintf(stderr, "Reference file: %s\n", opts.ref);
    } else {
        fprintf(stderr, "No reference file given\n");
        exit(EXIT_FAILURE);
    }
    if (opts.var != nullptr) {
        fprintf(stderr, "Variant file: %s\n", opts.var);
    } else {
        fprintf(stderr, "No variant file given\n");
        exit(EXIT_FAILURE);
    }
    if (opts.mqual < 0 || opts.mqual > 60) {
        fprintf(stderr, "Minimum mapping quality must be in the range [0, 60]\n");
        exit(EXIT_FAILURE);
    }
    if (opts.bqual < 0 || opts.bqual > 60) {
        fprintf(stderr, "Minimum base quality must be in the range [0, 60]\n");
        exit(EXIT_FAILURE);
    }
    if (opts.prefix == nullptr) {
        fprintf(stderr, "No filename prefix given\n");
        exit(EXIT_FAILURE);
    }
    if(opts.type != nullptr ){
        if (strcmp(opts.type, "SNV") !=0 && strcmp(opts.type, "CNV") !=0 && strcmp(opts.type, "SV") !=0 ) {
            fprintf(stderr, "Unrecognized input type [SNV,CNV,SV].\n");
            exit(EXIT_FAILURE);
        }
    }else{
        fprintf(stderr,"The -t,--type is not set.\n");
    }

    format = INPUT_SAM;
    size_t fn_len = strlen(opts.bam);
    if (fn_len >= 4 && strcmp(opts.bam + fn_len - 4, ".bam") == 0) {
        format = INPUT_BAM;
    } else if (fn_len >= 5 && strcmp(opts.bam + fn_len - 5, ".cram") == 0) {
        format = INPUT_CRAM;
    } else if (!(fn_len >= 4 && strcmp(opts.bam + fn_len - 4, ".sam") == 0)) {
         fprintf(stderr, "Input filename has an unrecognized file extension.\n");
         exit(EXIT_FAILURE);
    }

    switch (format) {
        case INPUT_BAM:
            strncpy(mode, "rb", 3);
            break;
        case INPUT_CRAM:
            strncpy(mode, "rc", 3);
            break;
        case INPUT_SAM:
        default:
            strncpy(mode, "r", 2);
            break;
    }

    /* Open input file */
    if ((input_sf = sam_open(opts.bam, mode)) == NULL) {
        fprintf(stderr, "Failed to open input file \"%s\".\n", opts.bam);
        exit(EXIT_FAILURE);
    }
    if (format == INPUT_CRAM) {
        size_t rec_buff_size = strlen(opts.ref) + 5;
        char* ref_buff = (char*)calloc(rec_buff_size, sizeof(char));
        if ( ref_buff==NULL ){
            fprintf(stderr, "Failed to locate memory.\n");
            exit(EXIT_FAILURE);
        }
        strncpy(ref_buff, opts.ref, rec_buff_size);
        strncat(ref_buff, ".fai", strlen(".fai"));

        if (hts_set_fai_filename(input_sf, ref_buff) != 0) {
            fprintf(stderr, "hts_set_fai_filename() failed for \"%s\".", ref_buff);
            free(ref_buff);
            exit(EXIT_FAILURE);
        }
        free(ref_buff);
    }
    if ((hdr = sam_hdr_read(input_sf)) == NULL) {
        fprintf(stderr, "Failed to read header for input file \"%s\".", opts.bam);
        exit(EXIT_FAILURE);
    }
    if ((index = sam_index_load(input_sf, opts.bam)) == NULL) {
        fprintf(stderr, "Failed to read index for input file \"%s\".", opts.bam);
        exit(EXIT_FAILURE);
    }
    //fprintf(stderr, "Reading Reference file:%s.\n",opts.ref);
    if ( (_fai = fai_load(opts.ref)) == NULL ) {
        fprintf(stderr, "Fail to read reference file \"%s\".", opts.ref);
        exit(EXIT_FAILURE);
    }
    if( strcmp(opts.type, "CNV") ==0 ){
        if( opts.list == nullptr ){
            fprintf(stderr ,"The CNV analysis referenece bam file is not set.");
            exit(EXIT_FAILURE);
        }
        if( opts.list == nullptr ){
            fprintf(stderr, "The mapability file is not set.");
            exit(EXIT_FAILURE);
        }
    }

    Reference reference(_fai);
    Variant varInfo(opts.var,opts.type);
    if( strcmp(opts.type, "CNV") ==0 ){
        varInfo.readCNV(opts.bed);
    }else if(strcmp(opts.type, "SV")==0){
        varInfo.readSV();
    }
    int total = varInfo.size();
    if( total == 0 ){
        fprintf(stderr, "No variant was detected for file \"%s\", Exit....\n", opts.var);
        exit(EXIT_FAILURE);
    }else{
        fprintf(stderr,"Total detected %d variants.\n", total);
    }
    for(int i=0;i<total;++i){
        VariantSpecific this_var = varInfo.Index(i);
        std::string chr = this_var.genomic_pos.front().chrName;
        const int chr_idx = sam_hdr_name2tid(hdr,chr.c_str());
        int32_t chr_len = (int32_t)hdr->target_len[chr_idx];

        int64_t start = (this_var.genomic_pos.front().start_ < opts.flank) ? 1 : (this_var.genomic_pos.front().start_ - opts.flank);
        int64_t end = (this_var.genomic_pos.back().end_ + opts.flank > chr_len) ? chr_len : (this_var.genomic_pos.back().end_ + opts.flank);
        int64_t var_length = end - start + 1;
        if( strcmp(opts.type, "SV")==0 || strcmp(opts.type, "SNV")==0 ){
            if( strcmp(opts.type, "SV")==0 ) var_length = 2*(opts.flank+1);
            fprintf(stderr,"var_length:%d\n", var_length);
            Coverage this_cover((int)var_length);
            fprintf(stderr, "Deal with: %s:%lld-%lld:%s>%s....\n", chr.c_str(), this_var.genomic_pos.front().start_, this_var.genomic_pos.back().end_, this_var.ref.c_str(),this_var.allele.c_str());
            Aligns this_align(opts.flank, reference, this_var, this_cover, opts.mqual, input_sf, hdr, index);
            std::cerr<<"start_offset:"<<std::to_string(this_align.get_start_offset())<<",end_offset:"<<std::to_string(this_align.get_end_offset())<<std::endl;
            //std::cerr<<this_cover<<std::endl;
            size_t fn_len = strlen(opts.prefix)+strlen(this_var.infor.c_str())+15;
	        char* outfile=(char*)calloc(fn_len,sizeof(char));
	   	    sprintf(outfile, "%s_%s_%s.svg", opts.prefix, opts.type, this_var.infor.c_str());
            generateSvg(this_align, outfile, this_cover, opts);
        }else if( strcmp(opts.type, "CNV")==0){
            Depth this_depth(opts.bed, opts.list, opts.bam);
            //CNV this_cnv(opts.flank, this_var, opts.mqual, this_depth, opts.map, reference);
            CNV this_cnv(opts.flank, this_var, opts.mqual, this_depth.get_bam_list(), this_depth.get_factors(), opts.map, reference, opts.gt);
            std::cerr<<"cnv start_offset:"<<std::to_string(this_cnv.get_ref_start())<<",size:"<<std::to_string(this_cnv.get_cnv_size())<<std::endl;
            size_t fn_len = strlen(opts.prefix)+strlen(opts.type)+1;
	        char* outfile=(char*)calloc(fn_len,sizeof(char));
	   	    sprintf(outfile, "%s_%s", opts.prefix, opts.type);
            generateSvg(this_cnv, outfile, opts);
        }
    }

    fai_destroy(_fai);
    sam_hdr_destroy(hdr);
    sam_close(input_sf);

    std::cerr << "All finished." << std::endl;
    return EXIT_SUCCESS;
  }
} // end namespace varViewer

int main(int argc, char **argv)
{
    return varViewer::varViewer_main(argc, argv);
}
