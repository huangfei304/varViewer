// varViewer
// Date: 2023-08-09
// Author: Huang Fei <huangfei@genomics.cn>
//


#pragma once

#include <cstdint>
#include <string>

namespace varViewer
{
enum class OperationType
{
    kMatch,
    kMismatch,
    kInsertionToRef,
    kDeletionFromRef,
    kSoftclip,
    kMissingBases
};

// Represents a single alignment operation
class Operation
{
public:
    Operation(OperationType type, uint32_t length)
        : type_(type)
        , length_(length)
    {
    }
    explicit Operation(std::string cigar);

    OperationType type() const { return type_; }
    uint32_t length() const { return length_; }
    uint32_t referenceLength() const;
    uint32_t queryLength() const;

    bool operator==(const Operation& other) const { return type_ == other.type_ && length_ == other.length_; }
    bool operator<(const Operation& other) const;

    std::string generateCigar() const;

private:
    OperationType type_;
    uint32_t length_;
};

OperationType decodeOperationType(char type_encoding);
std::ostream& operator<<(std::ostream& os, OperationType operation_type);
std::ostream& operator<<(std::ostream& os, const Operation& operation);
}
