#pragma once

#include <array>
#include <cstdint>
#include <vector>

// 8 x 8 矩阵
typedef std::array<std::array<int16_t, 8>, 8> JpegBlock;

void PrintBlock(const std::array<std::array<int16_t, 8>, 8>& block);

// 8 x 8 像素块，含有 3 个通道，一次为 Y，Cb，Cr
typedef std::array<JpegBlock, 3> UnitBlock;

// SOI: start of image

// APPn: Reserved for application segments
struct APPn {
    uint16_t Lp;
    std::vector<uint8_t> Api;
};

// SOFn: Start Of Frame markers
struct FrameHeader {
    FrameHeader() : Lf(0), P(0), Y(0), X(0), Nf(0) {}

    // Frame component-specification parameters
    struct ComponentParam {
        ComponentParam() : Ci(0), Hi(0), Vi(0), Tqi(0) {}
        uint8_t Ci;  // Component identifier
        uint8_t Hi;  // Horizontal sampling factor
        uint8_t Vi;  // Vertical sampling factor
        uint8_t Tqi; // Quantization table destination selector
    };

    uint16_t Lf; // Frame header length
    uint8_t P;   // Sample precision
    uint16_t Y;  // Number of lines
    uint16_t X;  // Number of samples per line
    uint8_t Nf;  // Number of image components in frame
    std::vector<ComponentParam> com_params; // Frame component-specification parameters
};

// SOS: scan header
struct ScanHeader {
    struct ComponentParam {
        uint8_t Csj; // Scan component selector
        uint8_t Tdj; // DC entropy coding table destination selector
        uint8_t Taj; // AC entropy coding table destination selector
    };

    uint16_t Ls; // Scan header length
    uint8_t Ns;  // Number of image components in scan
    std::vector<ComponentParam> com_params; // Scan component-specification parameters
    uint8_t Ss;  // Start of spectral or predictor selection
    uint8_t Se;  // End of spectral selection
    uint8_t Ah;  // Successive approximation bit position high
    uint8_t Al;  // Successive approximation bit position low or point transform
};

// DQT: Define quantization table(s)
struct DefQuantTable {
    struct QuantTable {
        uint8_t Pq;
        uint8_t Tq;
        std::array<int16_t, 64 > Qk;
    };

    uint16_t Lq;
    std::vector<QuantTable> tables;
};

// DHT: Define Huffman table(s)
struct DefHuffmanTable {
    struct HuffmanTable {
        uint8_t Tc;                                // Table class (0: DC; 1: AC)
        uint8_t Th;                                // Huffman table destination identifier (0: Y; 1: C)
        std::array<uint8_t, 16> Li;                // Number of Huffman codes of length i
        std::array<std::vector<uint8_t>, 16> Vij;  // Value associated with each Huffman code
    };

    uint16_t Lh;                      // Huffman table definition length
    std::vector<HuffmanTable> tables; // Huffman tables
};

// binary tree node
struct HuffmanNode {
    HuffmanNode() : val(0), left(nullptr), right(nullptr) {}
    HuffmanNode(uint8_t value) : val(value), left(nullptr), right(nullptr) {}
    ~HuffmanNode() { 
        if (left) 
            delete left;
        if (right)
            delete right;
    }

    uint8_t val;
    HuffmanNode* left;
    HuffmanNode* right;
};

// print Huffman tree
void PrintHuffmanTree(HuffmanNode* root);

// runLength encode
struct RLE {
    uint8_t runlength; // 之前的 0 的个数
    uint8_t size;      // amplitude 的位数
    int16_t amplitude; // 非 0 整数
};

// 返回矩阵在 zigzag 顺序下的坐标序列
std::vector<std::pair<int, int>> ZigzigOrderOfMatrix(int n);