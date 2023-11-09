//
// varViewer
//Date: 2023-08-09
// Author: Huang Fei <huangfei@genomics.cn>

#pragma once

#include <list>
#include <string>

#include "Operation.hh"

namespace varViewer
{
class Alignment
{
public:
    using size_type = size_t;
    using const_iterator = std::list<Operation>::const_iterator;

    Alignment(){
        reference_start_=0;
    }
    Alignment(int32_t reference_start, std::list<Operation> operations)
        : reference_start_(reference_start)
        , operations_(std::move(operations))
    {}
    Alignment(uint32_t reference_start, const std::string& cigar);
    const std::list<Operation>& operations() const { return operations_; }
    size_type numOperations() const { return operations_.size(); }
    uint32_t queryLength() const;
    uint32_t referenceLength() const;
    uint32_t referenceStart() const { return reference_start_; }

    const_iterator begin() const { return operations_.begin(); }
    const_iterator end() const { return operations_.end(); }
    const Operation& front() const { return operations_.front(); }
    const Operation& back() const { return operations_.back(); }
    size_type size() const { return operations_.size(); }
    bool operator==(const Alignment& other) const
    {
        return operations_ == other.operations_ && reference_start_ == other.reference_start_;
    }
    bool operator<(const Alignment& other) const;
    std::string generateCigar() const;
    Alignment splitAtReferencePosition(size_t reference_position, bool first=false);
    int32_t set_reference_start(int32_t start_){ reference_start_ = start_; }
    void reverseOperation();

protected:
    void decodeCigar(const std::string& encoding);

private:
    int32_t reference_start_ = 0;
    std::list<Operation> operations_;
};

std::ostream& operator<<(std::ostream& os, const Alignment& alignment);
}
