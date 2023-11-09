// varViewer
// Date: 2023-08-24
// Author: Huang Fei <huangfei@genomics.cn>
//
#include <cstring>
#include "Coverage.hh"

namespace varViewer {

Coverage::Coverage()
    : _covs(nullptr)
    , _len(0)
{
}

Coverage::Coverage(size_t len)
    : _covs(new int[len][COVS_ALIGN])
    , _len(len)
{
    this->reset();
}

Coverage::~Coverage(){
    delete[] _covs;
}

void Coverage::reset() {
    std::memset(_covs, 0, COVS_ALIGN * _len * sizeof(int));
}

void Coverage::add_coverage(int pos, covs_idx idx) {
    size_t p = (size_t)pos;
    if ( p < _len ) {
        ++_covs[p][idx];
    }
}

void Coverage::add_coverage(int pos, covs_idx idx, int add_cov) {
    size_t p = (size_t)pos;
    if( p < _len ){
        _covs[p][idx] += add_cov;
    }
}

void Coverage::add_coverage_range(int pos, covs_idx idx, int len){
    size_t p = (size_t)pos;
    size_t end = p + (size_t)len;
    if (end > _len) {
        end = _len;
    }
    while (p < end ) {
        ++_covs[p][idx];
        ++p;
    }
}

void Coverage::add_coverage_range(int pos, covs_idx idx, int len, int add_cov) {
    size_t p = (size_t)pos;
    size_t end = p + (size_t)len;
    if (end > _len) {
        end = _len;
    }
    while (p < end) {
        _covs[p][idx] += add_cov;
        ++p;
    }
}

void Coverage::set_coverage(int pos, covs_idx idx, int cov) {
    size_t p = (size_t)pos;
    if( p < _len ){
        _covs[p][idx] = cov;
    }
}

int Coverage::get_pos_coverage(int pos){
    size_t p = (size_t) pos;
    int total=0;
    if(p >= 0 && p < _len ){
        for(int idx=0; idx<COVS_ALIGN; ++idx) {
            total += _covs[p][idx];
        }
    }else{
        throw std::logic_error("pos: "+std::to_string(pos)+" is out of range (0,"+std::to_string(_len)+")");
    }
    return total;
}
int Coverage::get_coverage(int pos, covs_idx idx) {
    size_t p = (size_t)pos;
    return (p < _len) ? _covs[p][idx] : 0;
}

int Coverage::get_max_coverage() {
    int total=0;
    int p=0;
    while(p < _len ){
        if( total < get_pos_coverage(p)){
            total = get_pos_coverage(p);
        }
        ++p;
    }
    return total;
}

std::ostream& operator<<(std::ostream& os, Coverage& cover) {
    os<<"Index\tSNP_VR\tSNP_RR\tINS_VR\tINS_RR\tDEL_VR\tDEL_RR\tINSDEL_VR\tTotal\n";
    int p = 0;
    const int size = cover.get_len();
    while(p < size){
        os<<std::to_string(p);
        for(int idx=0; idx<COVS_ALIGN; ++idx) {
            os<<"\t"<<std::to_string(cover.get_coverage(p, (covs_idx)idx));
        }
        os<<std::to_string(cover.get_pos_coverage(p))<<"\n";
        ++p;
    }
    return os;
}

}
