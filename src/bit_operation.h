#pragma once

#include <cstdint>
#include <vector>

// 按大端模式读取数据
uint32_t GetBigEndianInteger(const uint8_t* data, int len);

// 获取高 4 位
uint8_t High4bit(uint8_t byte);

// 获取低 4 位
uint8_t Low4bit(uint8_t byte);

// 获取字节内的连续某几位
uint8_t NBitsInOneByte(uint8_t byte, int start, int end);

// 获取以第 start 位 bit 开始的 count 个 bit，并保存为32位无符号整型数
uint32_t GetBitsToInteger(const uint8_t* p, int start, int count);

// 获取以第 start 位 bit 开始的 count 个 bit
std::vector<bool> GetBits(const uint8_t* p, int start, int count);

// 获取第 n 个 bit
bool GetBit(const uint8_t* p, int bit_id);

// 用 size 个 bit 去编码 code
std::vector<bool> EncodeInNBits(uint16_t code, uint8_t size);

// 从 bit 集合获取一个无符号 16 位数
uint16_t NBitsToInteger(const std::vector<bool>& bits);