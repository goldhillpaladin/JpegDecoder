#include <algorithm>
#include <cassert>

#include "bit_operation.h"

uint32_t GetBigEndianInteger(const uint8_t* data, int len) {
    if (len <= 0 || 4 < len)
        return 0;

    uint32_t ret = 0;
    for (int i = 0; i < len; ++i) {
        ret <<= 8;
        ret |= data[i];
    }
    return ret;
}

uint8_t High4bit(uint8_t byte) {
    return byte >> 4;
}

uint8_t Low4bit(uint8_t byte) {
    return byte & 0xF;
}

uint8_t NBitsInOneByte(uint8_t byte, int start, int end) {
    assert(start <= end && start >= 0 && end <= 7);
    uint8_t ret = byte << start;
    ret >>= (start + 7 - end);
    return ret;
}

uint32_t GetBitsToInteger(const uint8_t* p, int start, int count) {
    int i = start / 8; // 第一个字节的下标
    int sub_start = start % 8; // 当前字节中要被截取的 bit 位的起始下标
    int sub_count = 0; // 当前字节中要被获得的 bit 位数
    int sub_end = 0;
    uint32_t ret = 0;
    while (count > 0) {
        sub_end = std::min(7, sub_start + count - 1);
        sub_count = sub_end - sub_start + 1;

        ret <<= sub_count;
        ret |= NBitsInOneByte(p[i], sub_start, sub_end);

        ++i;
        count -= sub_count;
        sub_start = 0;
    }
    return ret;
}

std::vector<bool> GetBits(const uint8_t* p, int start, int count) {
    std::vector<bool> bits; bits.reserve(count);
    int end = start + count - 1;
    for (int bit_id = start; bit_id <= end; ++bit_id)
        bits.push_back(GetBit(p, bit_id));
    return bits;
}

bool GetBit(const uint8_t* p, int bit_id) {
    return p[bit_id / 8] & (1 << (7 - bit_id % 8));
}

std::vector<bool> EncodeInNBits(uint16_t code, uint8_t size) {
    assert(size <= 16 && 0 < size);
    std::vector<bool> ret(size, false);
    for (uint8_t i = 0; i < size; ++i) {
        if (code & (1 << i))
            ret[size - 1 - i] = true;
    }
    return ret;
}

uint16_t NBitsToInteger(const std::vector<bool>& bits) {
    assert(0 < bits.size() && bits.size() <= 16);
    uint16_t value = 0;
    for (bool bit : bits) {
        value <<= 1;
        if (bit)
            value |= 1;
    }
    return value;
}