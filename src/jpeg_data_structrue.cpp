#include <iostream>

#include "jpeg_data_structrue.h"

void HuffmanTreeBacktrack(HuffmanNode* node, std::vector<bool>& bits) {
    if (!node)
        return;
    if (!node->left && !node->right) {
        std::cout << (int)node->val << ": ";
        for (bool bit : bits) {
            if (bit)
                std::cout << 1;
            else
                std::cout << 0;
        }
        std::cout << std::endl;
    }

    bits.push_back(false);
    HuffmanTreeBacktrack(node->left, bits);
    bits.pop_back();
    bits.push_back(true);
    HuffmanTreeBacktrack(node->right, bits);
    bits.pop_back();
}

void PrintHuffmanTree(HuffmanNode* root) {
    std::vector<bool> bits;
    HuffmanTreeBacktrack(root, bits);
}

void PrintBlock(const std::array<std::array<int16_t, 8>, 8>& block) {
    for (auto& row : block) {
        for (int16_t val : row)
            std::cout << val << " ";
        std::cout << std::endl;
    }
}

void PrintRgbBlock(const std::array<std::array<RgbPixel, 8>, 8>& rgb_block) {
    for (int i = 0; i < 8; ++i) {
        for (int j = 0; j < 8; ++j) {
            for (int k = 0; k < 3; ++k) {
                std::cout << rgb_block[i][j][k] << " ";
            }
            std::cout << "\t";
        }
        std::cout << std::endl;
    }
}

std::vector<std::pair<int, int>> ZigzigOrderOfMatrix(int n) {
    std::vector<std::pair<int, int>> zigzag;
    zigzag.push_back({ 0, 0 });

    int i = 0;
    int j = 1;

    for (int k = 1; k < n; ++k) {
        if (i == 0) {
            int count = j - i + 1;
            while (count-- > 0)
                zigzag.push_back({ i++, j-- });
            j = 0;
        }
        else {
            int count = i - j + 1;
            while (count-- > 0)
                zigzag.push_back({ i--, j++ });
            i = 0;
        }
    }

    if (i == 0) {
        i = 1;
        j = n - 1;
    }
    else {
        i = n - 1;
        j = 1;
    }

    for (int k = 1; k < n; ++k) {
        if (i == n - 1) {
            int count = i - j + 1;
            while (count-- > 0)
                zigzag.push_back({ i--, j++ });
            i += 2;
            j = n - 1;
        }
        else {
            int count = j - i + 1;
            while (count-- > 0)
                zigzag.push_back({ i++, j-- });
            j += 2;
            i = n - 1;
        }
    }
    return zigzag;
}