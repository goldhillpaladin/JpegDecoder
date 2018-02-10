#pragma once

#include <vector>

#include "jpeg_data_structrue.h"
#include "jpeg_marker_codes.h"

std::vector<FrameHeader::ComponentParam> DecodeFrameComponentParams(const uint8_t* data, int Nf);

FrameHeader DecodeFrameHeader(const uint8_t* data);

std::vector<ScanHeader::ComponentParam> DecodeScanComponentParams(const uint8_t* data, int Ns);

ScanHeader DecodeScanHeader(const uint8_t* data);

DefQuantTable DecodeDQT(const uint8_t* data);

DefHuffmanTable DecodeDefHuffTable(const uint8_t* data);

void ParseJpeg(const std::string& file_path);