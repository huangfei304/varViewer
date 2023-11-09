// varViewer
// Date: 2023-08-09
// Author: Huang Fei <huangfei@genomics.cn>
//

#include "LinearAlign.hh"

#include <algorithm>
#include <iostream>
#include <map>
#include <sstream>
#include <stdexcept>
#include <vector>

using std::list;
using std::logic_error;
using std::map;
using std::string;
using std::to_string;
using std::vector;

namespace varViewer
{

Alignment::Alignment(uint32_t reference_start, const string& cigar)
    : reference_start_(reference_start)
{
    decodeCigar(cigar);
}

void Alignment::decodeCigar(const string& cigar)
{
    string length_encoding;
    for (char c : cigar) {
        if (isalpha(c) != 0)
        {
            uint32_t operation_length = std::stoi(length_encoding);
            OperationType operation_type = decodeOperationType(c);
            operations_.emplace_back(operation_type, operation_length);
            length_encoding.clear();
        }else{
            if (isdigit(c) == 0){
                throw logic_error(cigar + " is malformed CIGAR string");
            }
            length_encoding += c;
        }
    }
}

uint32_t Alignment::queryLength() const
{
    int32_t query_span = 0;
    for (const auto& operation : operations_) {
        query_span += operation.queryLength();
    }
    return query_span;
}

uint32_t Alignment::referenceLength() const
{
    int32_t reference_span = 0;
    for (const auto& operation : operations_){
        reference_span += operation.referenceLength();
    }
    return reference_span;
}

string Alignment::generateCigar() const
{
    string cigar_string;
    for (const auto& operation : operations_) {
        cigar_string += operation.generateCigar();
    }
    return cigar_string;
}

void Alignment::reverseOperation(){
     reference_start_ += referenceLength();
     operations_.reverse();
}

Alignment Alignment::splitAtReferencePosition(size_t reference_position, bool first) {
    const size_t end_of_reference_positions = referenceStart() + referenceLength();
    if (reference_position == 0 || end_of_reference_positions <= reference_position) {
        std::ostringstream os;
        os << *this;
        throw logic_error("Cannot split " + os.str() + " at reference position " + to_string(reference_position));
    }

    size_t first_unused_position = reference_start_;
    list<Operation>::const_iterator operation_it = operations_.begin();
    while (operation_it != operations_.end()) {
        const size_t first_unused_position_after_applying_operation
            = first_unused_position + operation_it->referenceLength();
        if (first_unused_position_after_applying_operation <= reference_position){
            ++operation_it;
            first_unused_position = first_unused_position_after_applying_operation;
        }else{
            break;
        }
    }

    if (operation_it == operations_.end()) {
        throw logic_error("The must be some error in  splitAtReferencePosition");
    }

    if (first_unused_position == reference_position) {
        list<Operation> suffix_operations;
        if( first ){
            suffix_operations.splice(suffix_operations.begin(), operations_, operations_.begin(), --operation_it);
        }else{
            suffix_operations.splice(suffix_operations.begin(), operations_, operation_it, operations_.end());
        }
        return Alignment(first_unused_position, suffix_operations);
    } else {
        const size_t first_piece_reference_length = reference_position - first_unused_position;
        list<Operation> suffix_operations;
        if( first ){
            suffix_operations.splice(suffix_operations.begin(), operations_, operations_.begin(), --operation_it);
            Operation suffix_operation(operation_it->type(), first_piece_reference_length);
            suffix_operations.push_back(suffix_operation);
        }else{
            uint32_t suffix_reference_length = operation_it->referenceLength() - first_piece_reference_length;
            Operation suffix_operation(operation_it->type(), suffix_reference_length);
            suffix_operations.splice(suffix_operations.begin(), operations_, ++operation_it, operations_.end());
            suffix_operations.push_front(suffix_operation);
        }
        return Alignment(reference_position, suffix_operations);
    }
}

bool Alignment::operator<(const Alignment& other) const
{
    if (reference_start_ != other.reference_start_) {
        return reference_start_ < other.reference_start_;
    }

    return operations_ < other.operations_;
}

std::ostream& operator<<(std::ostream& os, const Alignment& alignment)
{
    os << "Ref start: " << alignment.referenceStart() << ", ";
    for (const Operation& operation : alignment) {
        os << operation;
    }

    return os;
}
}
