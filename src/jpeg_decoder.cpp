#include <cassert>
#include <cmath>
#include <iostream>
#include <memory>
#include <unordered_map>

#include "bit_operation.h"
#include "jpeg_decoder.h"

#pragma warning(disable : 4996)

const float PI = 3.14159265f;
const float PI_OVER_SIXTEEN = PI / 16;

const uint8_t* data_begin_;
const uint8_t* data_end_;
std::vector<std::vector<uint8_t>> APPn_data_vec_;
std::vector<DefQuantTable> DQT_vec_;
FrameHeader SOF0_;
std::vector<DefHuffmanTable> DHT_vec_;
std::vector<DefHuffmanTable::HuffmanTable> huffman_tables_;
ScanHeader SOS_;
std::vector<uint8_t> entropy_code_;

std::vector<FrameHeader::ComponentParam> DecodeFrameComponentParams(const uint8_t* data, int Nf) {
    std::vector<FrameHeader::ComponentParam> com_params(Nf);
    int len = Nf * 3;
    for (int i = 0; i < len; i += 3) {
        FrameHeader::ComponentParam& com_param = com_params[i / 3];
        com_param.Ci = data[i];
        com_param.Hi = data[i + 1] >> 4;
        com_param.Vi = data[i + 1] & 0xF;
        com_param.Tqi = data[i + 2];
    }
    return com_params;
}

FrameHeader DecodeFrameHeader(const uint8_t* data) {
    FrameHeader header;

    header.Lf = (((uint16_t)data[0]) << 8) | data[1];
    size_t i = 2;
    header.P = data[i]; ++i;
    header.Y = (((uint16_t)data[i]) << 8) | data[i + 1]; i += 2;
    header.X = (((uint16_t)data[i]) << 8) | data[i + 1]; i += 2;
    header.Nf = data[i]; ++i;
    header.com_params = DecodeFrameComponentParams(data + i, header.Nf);

    return header;
}

std::vector<ScanHeader::ComponentParam> DecodeScanComponentParams(const uint8_t* data, int Ns) {
    std::vector<ScanHeader::ComponentParam> com_params(Ns);
    int len = Ns * 2;
    for (int i = 0; i < len; i += 2) {
        ScanHeader::ComponentParam& param = com_params[i / 2];
        param.Csj = data[i];
        param.Tdj = data[i + 1] >> 4;
        param.Taj = data[i + 1] & 0xF;
    }
    return com_params;
}

ScanHeader DecodeScanHeader(const uint8_t* data) {
    ScanHeader header;

    size_t i = 0;
    header.Ls = (((uint16_t)data[i]) << 8) | data[i + 1]; i += 2;
    header.Ns = data[i]; ++i;
    header.com_params = DecodeScanComponentParams(data + i, header.Ns); i += 2 * header.Ns;
    header.Ss = data[i]; ++i;
    header.Se = data[i]; ++i;
    header.Ah = data[i] >> 4;
    header.Al = data[i] & 0xF;

    return header;
}

DefQuantTable DecodeDQT(const uint8_t* data) {
    DefQuantTable dqt;

    size_t i = 0;
    dqt.Lq = (((uint16_t)data[i]) << 8) | data[i + 1]; i += 2;
    while (i < dqt.Lq) {
        DefQuantTable::QuantTable table;
        table.Pq = data[i] >> 4;
        table.Tq = data[i] & 0xF; ++i;
        for (int count = 0; count < 64; ++count) {
            if (table.Pq == 0)
                table.Qk[count] = data[i];
            else
                table.Qk[count] = (((uint16_t)data[i]) << 8) | data[i + 1];
            i += table.Pq + 1;
        }
        dqt.tables.push_back(table);
    }

    return dqt;
}

DefHuffmanTable DecodeDefHuffTable(const uint8_t* data) {
    DefHuffmanTable dht;

    size_t i = 0;
    dht.Lh = GetBigEndianInteger(data, 2); i += 2;
    while (i < dht.Lh) {
        DefHuffmanTable::HuffmanTable table;

        table.Tc = High4bit(data[i]);
        table.Th = Low4bit(data[i]); ++i;
        for (int id = 0; id < 16; ++id)
            table.Li[id] = data[i++];
        for (int id = 0; id < 16; ++id) {
            uint8_t Li = table.Li[id];
            auto& Vi = table.Vij[id];
            Vi.reserve(Li);
            for (int count = 0; count < Li; ++count) {
                Vi.push_back(data[i++]);
            }
        }

        dht.tables.push_back(table);
    }

    return dht;
}

std::vector<uint8_t> FileToVector(const std::string& file_path) {
    std::vector<uint8_t> bytes;

    FILE* fp = fopen(file_path.c_str(), "rb");
    if (!fp) {
        std::cout << "The file is corrupt" << std::endl;
        return bytes;
    }

    fseek(fp, 0, SEEK_END);
    size_t file_size = ftell(fp);
    fseek(fp, 0, SEEK_SET);
    bytes.reserve(file_size);
    uint8_t buf[1024] = { 0 };
    for (int count = (int)file_size; count > 0; count -= 1024) {
        if (count >= 1024) {
            fread(buf, 1, 1024, fp);
            bytes.insert(bytes.end(), buf, buf + 1024);
        }
        else {
            fread(buf, 1, count, fp);
            bytes.insert(bytes.end(), buf, buf + count);
        }
    }
    fclose(fp);
    fp = nullptr;
    return bytes;
}

const uint8_t* FindSOI(const uint8_t* data, size_t data_len) {
    const uint8_t* p = data;
    while (p - data < (int)data_len) {
        if (*p == 0xFF && *(p + 1) == 0xD8)
            return p;
        ++p;
    }
    return nullptr;
}

std::vector<std::vector<uint8_t>> GetAPPnData(const uint8_t*& data) {
    std::vector<std::vector<uint8_t>> apps;
    const uint8_t* p = data;
    while (p[0] == 0xFF && (p[1] & 0xF0) == 0xE0) {
        uint8_t marker_code[2] = { p[0], p[1] }; p += 2;
        uint16_t La = GetBigEndianInteger(p, 2);
        std::vector<uint8_t> app; app.reserve(2 + La);
        app.insert(app.end(), marker_code, marker_code + 2);
        app.insert(app.end(), p, p + La);
        data += app.size();
        p = data;
        apps.push_back(std::move(app));
    }
    return apps;
}

DefQuantTable GetDQT(const uint8_t*& data) {
    data += 2;
    DefQuantTable table = DecodeDQT(data);
    data += table.Lq;
    return table;
}

std::vector<DefQuantTable::QuantTable> GetQuantTables(const std::vector<DefQuantTable>& dqt_vec) {
    std::vector<DefQuantTable::QuantTable> tables;
    for (auto& dqt : dqt_vec) {
        for (auto& table : dqt.tables)
            tables.push_back(table);
    }
    return tables;
}

FrameHeader GetSOF(const uint8_t*& data) {
    data += 2;
    FrameHeader header = DecodeFrameHeader(data);
    data += header.Lf;
    return header;
}

DefHuffmanTable GetDHT(const uint8_t*& data) {
    data += 2;
    DefHuffmanTable table = DecodeDefHuffTable(data);
    data += table.Lh;
    return table;
}

ScanHeader GetSOS(const uint8_t*& data) {
    data += 2;
    ScanHeader header = DecodeScanHeader(data);
    data += header.Ls;
    return header;
}

std::vector<uint8_t> GetEntropyCode(const uint8_t*& data) {
    std::vector<uint8_t> ret;

    for (;;) {
        if (*data == 0xFF) {
            if (*(data + 1) != 0)
                break;
            if (*(data + 1) == 0) {
                ret.push_back(*data);
                data += 2;
            }
        }
        else {
            ret.push_back(*data);
            ++data;
        }
    }

    return ret;
}

/*const uint8_t* ParseMarkerSegment(const uint8_t* data) {
    const uint8_t* p = data;
    while (p < data_end) {
        if (*p == 0xFF && *(p + 1) != 0)
            break;
    }
    if (p == data_end)
        return nullptr;

    switch (*p) {
    case 0xE0:
    case 0xE1:
    case 0xE2:
    case 0xE3:
    case 0xE4:
    case 0xE5:
    case 0xE6:
    case 0xE7:
    case 0xE8:
    case 0xE9:
    case 0xEA:
    case 0xEB:
    case 0xEC:
    case 0xED:
    case 0xEE:
    case 0xEF:
        APPn_data = GetAPPnData(p);
        break;
    case:
    }
    
    return p;
}*/

// Generate size table
std::vector<uint8_t> GenerateHuffmanSize(const DefHuffmanTable::HuffmanTable& huff_table) {
    std::vector<uint8_t> HuffSize;
    for (int len = 1; len <= 16; ++len) {
        if (huff_table.Li[len - 1] == 0)
            continue;
        for (int i = 0; i < huff_table.Li[len - 1]; ++i) {
            HuffSize.push_back(len);
        }
    }
    HuffSize.push_back(0);
    int last_k = (int)HuffSize.size() - 1;
    return HuffSize;
}

// Generate code table
std::vector<uint16_t> GenerateHuffmanCodeTable(const std::vector<uint8_t>& HuffSize) {
    std::vector<uint16_t> HuffCode(HuffSize.size(), 0);
    uint8_t value = 0;
    uint16_t code = 0;
    int size = HuffSize[0];
    for (;;) {
        HuffCode[value] = code;
        ++code;
        ++value;
        if (HuffSize[value] != size) {
            if (HuffSize[value] == 0)
                break;
            do {
                code <<= 1;
                ++size;
            } while (HuffSize[value] != size);
        }
    }
    return HuffCode;
}

std::unordered_map<uint8_t, std::unordered_map<uint8_t, uint16_t>> OrderHuffmanCodeTable(const std::vector<uint16_t>& unordered_huff_codes, const DefHuffmanTable::HuffmanTable& huff_table) {
    std::unordered_map<uint8_t, std::unordered_map<uint8_t, uint16_t>> ret;
    int id = 0;
    for (int i = 0; i < (int)huff_table.Vij.size(); ++i) {
        for (uint8_t Vij : huff_table.Vij[i]) {
            ret[i+1][Vij] = unordered_huff_codes[id++];
        }
    }
    return ret;
}

std::shared_ptr<HuffmanNode> ConstructHuffmanTree(const std::unordered_map<uint8_t, std::unordered_map<uint8_t, uint16_t>>& HuffCode,
                                                  const std::vector<uint8_t>& HuffSize) {
    std::shared_ptr<HuffmanNode> root(new HuffmanNode);

    auto InsertIntoHuffmanTree = [](HuffmanNode* root, uint8_t val, const std::vector<bool>& bits) {
        HuffmanNode* node = root;
        for (bool bit : bits) {
            if (!bit) {
                if (!node->left)
                    node->left = new HuffmanNode;
                node = node->left;
            }
            else {
                if (!node->right)
                    node->right = new HuffmanNode;
                node = node->right;
            }
        }
        node->val = val;
        assert(!node->left && !node->right);
    };

    for (const auto& size : HuffCode) {
        uint8_t bit_count = size.first;
        for (auto& val_code : size.second) {
            uint8_t val = val_code.first;
            std::vector<bool> bits = EncodeInNBits(val_code.second, bit_count);
            InsertIntoHuffmanTree(root.get(), val, bits);
        }
    }

    return root;
}

std::shared_ptr<HuffmanNode> HuffmanTableToHuffmanTree(const DefHuffmanTable::HuffmanTable& table) {
    std::vector<uint8_t> HuffSize = GenerateHuffmanSize(table);
    std::vector<uint16_t> originalHuffCode = GenerateHuffmanCodeTable(HuffSize);
    auto HuffCode = OrderHuffmanCodeTable(originalHuffCode, table);
    auto huff_root = ConstructHuffmanTree(HuffCode, HuffSize);
    return huff_root;
}

// bits 是否位于 Huffman 树中
bool ExistHuffmanCode(HuffmanNode* huff_root, const std::vector<bool>& bits) {
    HuffmanNode* node = huff_root;
    int i = 0;
    while (i < (int)bits.size() && (node->left || node->right)) {
        bool bit = bits[i++];
        node = bit ? node->right : node->left;
    }
    return !node->left && !node->right;
}

uint8_t DecodeHuffman(HuffmanNode* huff_root, const uint8_t* p, int& bit_id) {
    HuffmanNode* node = huff_root;
    while (node->left || node->right) {
        bool bit = GetBit(p, bit_id++);
        node = bit ? node->right : node->left;
    }
    return node->val;
}

std::vector<bool> IntegerToRLEBits(int16_t value) {
    if (value == 0)
        return {};

    auto GetBitCount = [](uint16_t value) -> int {
        int count = 0;
        while (value) {
            ++count;
            value >>= 1;
        }
        return count;
    };

    std::vector<bool> bits;
    if (value > 0) {
        int bit_count = GetBitCount((uint16_t)value);
        bits = EncodeInNBits((uint16_t)value, bit_count);
    }
    else {
        int bit_count = GetBitCount((uint16_t)abs(value));
        bits = EncodeInNBits((uint16_t)(value - 1), bit_count);
    }
    return bits;
}

int16_t RLEBitsToInteger(const std::vector<bool>& bits) {
    if (bits.empty())
        return 0;

    int16_t value = 0;
    if (bits[0]) {
        value = NBitsToInteger(bits);
    }
    else {
        auto tmp_bits = bits;
        for (int i = 0; i < (int)tmp_bits.size(); ++i)
            tmp_bits[i] = !tmp_bits[i];
        value = -NBitsToInteger(tmp_bits);
    }
    return value;
}

int16_t RLEBitsToInteger(const uint8_t* p, int bit_id, int bit_count) {
    std::vector<bool> bits = GetBits(p, bit_id, bit_count);
    return RLEBitsToInteger(bits);
}

JpegBlock ZigzagToMatrix(const std::array<int16_t, 64>& zigzag) {
    auto coords = ZigzigOrderOfMatrix(8);
    JpegBlock matrix;
    for (int i = 0; i < 64; ++i) {
        int x = coords[i].first;
        int y = coords[i].second;
        matrix[x][y] = zigzag[i];
    }
    return matrix;
}

JpegBlock MCUToMatrix(const std::vector<RLE>& mcu) {
    std::array<int16_t, 64> zigzag;
    int i = 0;
    for (auto& rle : mcu) {
        if (rle.runlength == 0 && rle.size == 0) { // EOB
            while (i < 64)
                zigzag[i++] = 0;
            break;
        }
        int pre_zero_count = rle.runlength;
        for (int j = 0; j < rle.runlength; ++j) {
            zigzag[i++] = 0;
        }
        zigzag[i++] = rle.amplitude;
    }
    return ZigzagToMatrix(zigzag);
}

JpegBlock EntropyCodeToQuantizedDCTMatrix(const uint8_t* code, int& bit_id, 
    std::shared_ptr<HuffmanNode> dc_root, std::shared_ptr<HuffmanNode> ac_root) {
    std::vector<RLE> mcu;
    int element_count = 0;

    // DC 系数
    RLE rle;
    uint8_t runlen_size = DecodeHuffman(dc_root.get(), code, bit_id);
    rle.runlength = High4bit(runlen_size);
    rle.size = Low4bit(runlen_size);
    rle.amplitude = RLEBitsToInteger(code, bit_id, rle.size); bit_id += rle.size;
    mcu.push_back(rle);
    ++element_count;

    // Y 通道的 AC 系数
    while (element_count < 64) {
        runlen_size = DecodeHuffman(ac_root.get(), code, bit_id);
        rle.runlength = High4bit(runlen_size);
        rle.size = Low4bit(runlen_size);
        rle.amplitude = RLEBitsToInteger(code, bit_id, rle.size); bit_id += rle.size;
        mcu.push_back(rle);
        element_count += rle.size + 1;
        if (runlen_size == 0) // EOB
            break;
    }
    JpegBlock block = MCUToMatrix(mcu);
    return block;
}

// 反量化
JpegBlock DeQuantize(const JpegBlock& quantized, const DefQuantTable::QuantTable& table) {
    JpegBlock origin;
    
    JpegBlock quant_table = ZigzagToMatrix(table.Qk);
    for (int i = 0; i < 8; ++i) {
        for (int j = 0; j < 8; ++j) {
            origin[i][j] = quantized[i][j] * quant_table[i][j];
        }
    }

    return origin;
}

// 由 DIFF 转换成真实的 DC 值
void AmendQuantizedDCTMatrix(std::vector<UnitBlock>& unit_blocks) {
    for (int component = 0; component < 3; ++component) {
        int16_t prev_dc = 0;
        for (int i = 0; i < (int)unit_blocks.size(); ++i) {
            unit_blocks[i][component][0][0] = unit_blocks[i][component][0][0] + prev_dc;
            prev_dc = unit_blocks[i][component][0][0];
        }
    }
}

double f1(uint8_t x, uint8_t y, const JpegBlock& dct_matrix) {
    auto Alpha = [](uint8_t u) -> double {
        return (u == 0) ? 0.70710678 : 1.0f;
    };
    double ret = 0;
    double tmp_x = ((double)(2 * x + 1))*PI_OVER_SIXTEEN;
    double tmp_y = ((double)(2 * y + 1))*PI_OVER_SIXTEEN;
    for (int u = 0; u < 8; ++u) {
        for (int v = 0; v < 8; ++v) {
            double tmp = Alpha(u);
            tmp = Alpha(v);
            tmp = dct_matrix[u][v];
            tmp = std::cos(tmp_x*u);
            tmp = std::cos(tmp_y*u);
            double dd = Alpha(u) * Alpha(v) * dct_matrix[u][v] * std::cos(tmp_x*u) * std::cos(tmp_y*v);
            ret += Alpha(u) * Alpha(v) * dct_matrix[u][v] * std::cos(tmp_x*u) * std::cos(tmp_y*v);
        }
    }
    return ret / 4;
}

// 逆 DCT 变换
JpegBlock IDCT8X8(const JpegBlock& dct_matrix) {
    JpegBlock ret;
    for (int x = 0; x < 8; ++x) {
        for (int y = 0; y < 8; ++y) {
            ret[x][y] = (int16_t)round(f1(x, y, dct_matrix));
        }
    }
    return ret;
}

// 每项 + 128
void ShiftSubImage(JpegBlock& sub_image) {
    for (auto& row : sub_image)
        for (auto& pixel : row)
            pixel += 128;
}

//R = Y + 1.402(Cr - 128)
//G = Y - 0.34414(Cb - 128) - 0.71414(Cr - 128)
//B = Y + 1.772(Cb - 128)
std::array<int16_t, 3> YCbCrToRBG(int16_t Y, int16_t Cb, int16_t Cr) {
    std::array<int16_t, 3> rgb;
    rgb[0] = int16_t(Y + 1.402*(Cr - 128));
    rgb[1] = int16_t(Y - 0.34414*(Cb - 128) - 0.71414*(Cr - 128));
    rgb[2] = int16_t(Y + 1.772*(Cb - 128));
    return rgb;
}

std::array<std::array<RgbPixel, 8>, 8> BlockYCbCrToRBG(const UnitBlock& block) {
    const JpegBlock& Y = block[0];
    const JpegBlock& Cb = block[1];
    const JpegBlock& Cr = block[2];
    std::array<std::array<RgbPixel, 8>, 8> ret;
    for (int i = 0; i < 8; ++i) {
        for (int j = 0; j < 8; ++j) {
            ret[i][j] = YCbCrToRBG(Y[i][j], Cb[i][j], Cr[i][j]);
        }
    }
    return ret;
}

void DeHuffmanEntropyCode(const uint8_t* code, int code_len, const std::vector<DefHuffmanTable::HuffmanTable>& huffman_tables) {
    uint8_t sample_precision = SOF0_.P; // 基线顺序模式下，为 8 bit
    uint16_t vertical_size = SOF0_.Y;
    uint16_t horizontal_size = SOF0_.X;
    uint8_t component_count = SOS_.Ns;

    // get quant tables
    DefQuantTable::QuantTable y_quant_table;
    DefQuantTable::QuantTable c_quant_table;
    auto quant_tables = GetQuantTables(DQT_vec_);
    for (auto& table : quant_tables) {
        if (table.Tq == 0)
            y_quant_table = table;
        else
            c_quant_table = table;
    }
    
    // generate huffman trees
    std::shared_ptr<HuffmanNode> y_dc_root(nullptr);
    std::shared_ptr<HuffmanNode> y_ac_root(nullptr);
    std::shared_ptr<HuffmanNode> c_dc_root(nullptr);
    std::shared_ptr<HuffmanNode> c_ac_root(nullptr);
    for (auto& table : huffman_tables) {
        if (table.Tc == 0) { // DC
            if (table.Th == 0)
                y_dc_root = HuffmanTableToHuffmanTree(table); // Y component
            else
                c_dc_root = HuffmanTableToHuffmanTree(table); // CbCr component
        }
        else { // AC
            if (table.Th == 0)
                y_ac_root = HuffmanTableToHuffmanTree(table); // Y component
            else
                c_ac_root = HuffmanTableToHuffmanTree(table); // CbCr component
        }
    }
    //PrintHuffmanTree(y_dc_root.get());

    // generate quantized DCT
    std::vector<UnitBlock> unit_blocks;
    int bit_id = 0;
    while (bit_id < code_len * 8) {
        if ((bit_id + 7) / 8 >= code_len) {
            auto bits = GetBits(code, bit_id, code_len * 8 - bit_id);
            if (!ExistHuffmanCode(y_dc_root.get(), bits))
                break;
        }
        UnitBlock unit_block;
        unit_block[0] = EntropyCodeToQuantizedDCTMatrix(code, bit_id, y_dc_root, y_ac_root);
        unit_block[1] = EntropyCodeToQuantizedDCTMatrix(code, bit_id, c_dc_root, c_ac_root);
        unit_block[2] = EntropyCodeToQuantizedDCTMatrix(code, bit_id, c_dc_root, c_ac_root);
        unit_blocks.push_back(unit_block);
    }
    AmendQuantizedDCTMatrix(unit_blocks);

    // 对 quantized DCT 反量化，获得原始的 DCT
    std::vector<UnitBlock> origin_dct_sets(unit_blocks.size());
    for (int i = 0; i < (int)unit_blocks.size(); ++i) {
        origin_dct_sets[i][0] = DeQuantize(unit_blocks[i][0], y_quant_table);
        origin_dct_sets[i][1] = DeQuantize(unit_blocks[i][1], c_quant_table);
        origin_dct_sets[i][2] = DeQuantize(unit_blocks[i][2], c_quant_table);
    }

    // 对 DCT 进行 IDCT 操作，获得降采样的源像素点集合
    std::vector<UnitBlock> idct_sets(origin_dct_sets.size());
    for (int i = 0; i < (int)origin_dct_sets.size(); ++i) {
        idct_sets[i][0] = IDCT8X8(origin_dct_sets[i][0]);
        idct_sets[i][1] = IDCT8X8(origin_dct_sets[i][1]);
        idct_sets[i][2] = IDCT8X8(origin_dct_sets[i][2]);
        for (int j = 0; j < 3; ++j)
            ShiftSubImage(idct_sets[i][j]);
        PrintBlock(idct_sets[i][0]); std::cout << std::endl;
        PrintBlock(idct_sets[i][1]); std::cout << std::endl;
        PrintBlock(idct_sets[i][2]); std::cout << std::endl;
    }

    // 反降采样
    // 现在只支持 1:1 采样
    std::vector<UnitBlock> YCbCr_blocks = idct_sets;

    // YCbCr -> RGB
    std::vector<std::array<std::array<RgbPixel, 8>, 8>> RGB_blocks;
    for (auto& YCbCr_block : YCbCr_blocks) {
        RGB_blocks.push_back(BlockYCbCrToRBG(YCbCr_block));
    }
    std::cout << "RGB:\n";
    for (auto& rgb_block : RGB_blocks) {
        PrintRgbBlock(rgb_block);
        std::cout << std::endl;
    }
}

void DecodeJpeg() {

}

void ParseJpeg(const std::string& file_path) {
    std::vector<uint8_t> bytes = FileToVector(file_path);
    const uint8_t* data = bytes.data();
    data_begin_ = data;
    data_end_ = data + bytes.size();
    const uint8_t* p = data;
    
    p = FindSOI(data, bytes.size());
    if (!p) {
        std::cout << "The file is corrupt" << std::endl;
        return;
    }
    p += 2;

    APPn_data_vec_ = GetAPPnData(p);

    while (p[0] == 0xFF && p[1] == DQT) {
        DQT_vec_.push_back(GetDQT(p));
    }

    SOF0_ = GetSOF(p);

    while (p[0] == 0xFF && p[1] == DHT) {
        DHT_vec_.push_back(GetDHT(p));
    }
    for (const DefHuffmanTable& dht : DHT_vec_) {
        for (const DefHuffmanTable::HuffmanTable& huffman_table : dht.tables) {
            huffman_tables_.push_back(huffman_table);
        }
    }

    SOS_ = GetSOS(p);

    entropy_code_ = GetEntropyCode(p);
    printf("0x%x\n", (unsigned int)(p - data_begin_));

    //DecodeJpeg();

    DeHuffmanEntropyCode(entropy_code_.data(), (int)entropy_code_.size(), huffman_tables_);
}