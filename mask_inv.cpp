#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <sstream>
#include <iomanip>

int main() {
    const std::string INPUT_FILE = "C:/Users/User/c++1/mask-mura-61.txt";
    const std::string OUTPUT_FILE = "C:/Users/User/c++1/mask-mura-expanded_inv.txt";
    const int BLOCK_SIZE = 32;          // Увеличенный размер блока
    const double RADIUS = 12.0;         // 3.0 * (32/8) = 12.0
    const double CENTER = 15.5;         // Центр блока 32×32: (31/2 = 15.5)

    // Чтение исходной маски
    std::ifstream inFile(INPUT_FILE);
    if (!inFile.is_open()) {
        std::cerr << "Ошибка: не удалось открыть входной файл " << INPUT_FILE << std::endl;
        return 1;
    }

    std::vector<std::vector<int>> mask;
    std::string line;
    while (std::getline(inFile, line)) {
        if (line.empty()) continue;
        std::istringstream iss(line);
        std::vector<int> row;
        int val;
        while (iss >> val) row.push_back(val);
        if (!row.empty()) mask.push_back(row);
    }
    inFile.close();

    int origH = mask.size();
    int origW = (origH > 0) ? mask[0].size() : 0;
    std::cout << "Исходный размер маски: " << origH << " x " << origW << std::endl;

    // Генерация блока с окружностью (для элементов '1')
    std::vector<std::vector<int>> circleBlock(BLOCK_SIZE, std::vector<int>(BLOCK_SIZE, 0));
    for (int i = 0; i < BLOCK_SIZE; ++i) {
        for (int j = 0; j < BLOCK_SIZE; ++j) {
            double dx = i - CENTER;
            double dy = j - CENTER;
            if (dx * dx + dy * dy <= RADIUS * RADIUS) {
                circleBlock[i][j] = 1;
            }
        }
    }

    // Генерация блока нулей (для элементов '0') — один раз для эффективности
    std::vector<std::vector<int>> zeroBlock(BLOCK_SIZE, std::vector<int>(BLOCK_SIZE, 0));

    // Формирование результата
    int newH = origH * BLOCK_SIZE;
    int newW = origW * BLOCK_SIZE;
    std::vector<std::vector<int>> result(newH, std::vector<int>(newW, 0));

    for (int i = 0; i < origH; ++i) {
        for (int j = 0; j < origW; ++j) {
            const auto& block = (mask[i][j] == 1) ? circleBlock : zeroBlock;
            int base_i = i * BLOCK_SIZE;
            int base_j = j * BLOCK_SIZE;
            for (int bi = 0; bi < BLOCK_SIZE; ++bi) {
                for (int bj = 0; bj < BLOCK_SIZE; ++bj) {
                    result[base_i + bi][base_j + bj] = block[bi][bj];
                }
            }
        }
    }

    // Запись в файл
    std::ofstream outFile(OUTPUT_FILE);
    if (!outFile.is_open()) {
        std::cerr << "Ошибка: не удалось создать выходной файл " << OUTPUT_FILE << std::endl;
        return 1;
    }

    for (int i = 0; i < newH; ++i) {
        for (int j = 0; j < newW; ++j) {
            outFile << result[i][j];
            if (j < newW - 1) outFile << " ";
        }
        outFile << "\n";
    }
    outFile.close();

    std::cout << "Готово! Результат сохранён в " << OUTPUT_FILE << std::endl;
    std::cout << "Новый размер: " << newH << " x " << newW << " (ожидаемо: " 
              << origH * BLOCK_SIZE << " x " << origW * BLOCK_SIZE << ")" << std::endl;
    
    // Дополнительная информация о структуре блока
    int circle_pixels = 0;
    for (const auto& row : circleBlock)
        for (int val : row) circle_pixels += val;
    std::cout << "Блок для '1': " << circle_pixels << " пикселей со значением 1 из " 
              << BLOCK_SIZE * BLOCK_SIZE << " (радиус=" << RADIUS << ", центр=" << CENTER << ")" << std::endl;

    return 0;
}