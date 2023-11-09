// varViewer
// Date: 2023-08-09
// Author: Huang Fei <huangfei@genomics.cn>
//

#include "Operation.hh"
#include <algorithm>
#include <cassert>
#include <iostream>
#include <map>
#include <sstream>
#include <stdexcept>

using std::logic_error;
using std::map;
using std::string;
using std::to_string;

namespace varViewer
{

string Operation::generateCigar() const {
    string cigar_string = to_string(length_);
    std::ostringstream os;
    os << type_;
    cigar_string += os.str();

    return cigar_string;
}

OperationType decodeOperationType(char type_encoding)
{
    switch (type_encoding)
    {
    case '=':
        return OperationType::kMatch;
    case 'M':
        return OperationType::kMatch;
    case 'N':
        return OperationType::kMissingBases;
    case 'X':
        return OperationType::kMismatch;
    case 'I':
        return OperationType::kInsertionToRef;
    case 'D':
        return OperationType::kDeletionFromRef;
    case 'S':
        return OperationType::kSoftclip;
    default:
        throw logic_error(to_string(type_encoding) + " is unknown CIGAR operation");
    }
}

Operation::Operation(string cigar)
{
    type_ = decodeOperationType(cigar.back());
    cigar.pop_back();
    length_ = std::stoi(cigar);
}

uint32_t Operation::referenceLength() const
{
    switch (type_)
    {
    case OperationType::kMatch:
    case OperationType::kMismatch:
    case OperationType::kMissingBases:
    case OperationType::kDeletionFromRef:
        return length_;
    default:
        return 0;
    }
}

uint32_t Operation::queryLength() const
{
    switch (type_)
    {
    case OperationType::kMatch:
    case OperationType::kMismatch:
    case OperationType::kMissingBases:
    case OperationType::kInsertionToRef:
    case OperationType::kSoftclip:
        return length_;
    default:
        return 0;
    }
}

bool Operation::operator<(const Operation& other) const
{
    if (type_ != other.type_){
        return type_ < other.type_;
    }

    return length_ < other.length_;
}
std::ostream& operator<<(std::ostream& os, OperationType operation_type)
{
    switch (operation_type)
    {
    case OperationType::kMatch:
        os << 'M';
        break;
    case OperationType::kMismatch:
        os << 'X';
        break;
    case OperationType::kInsertionToRef:
        os << 'I';
        break;
    case OperationType::kDeletionFromRef:
        os << 'D';
        break;
    case OperationType::kSoftclip:
        os << 'S';
        break;
    case OperationType::kMissingBases:
        os << 'N';
    }

    return os;
}

std::ostream& operator<<(std::ostream& os, const Operation& operation)
{
    os << operation.length() << operation.type();
    return os;
}
}
