
#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <sstream>
#include <iomanip>

int main() {
    const std::string INPUT_FILE = "C:/Users/User/c++1/mask-mura-61.txt";
    const std::string OUTPUT_FILE = "C:/Users/User/c++1/mask-mura-expanded.txt";
    const int BLOCK_SIZE = 32;          // Размер блока 32x32
    const double RADIUS = 12.0;         // Масштабированный радиус окружности
    const double CENTER = 15.5;         // Центр блока 32x32

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

    // Блок для элементов '1': полностью заполнен 1
    std::vector<std::vector<int>> fullBlock(BLOCK_SIZE, std::vector<int>(BLOCK_SIZE, 1));

    // Блок для элементов '0': внутри окружности 0, снаружи 1
    std::vector<std::vector<int>> holeBlock(BLOCK_SIZE, std::vector<int>(BLOCK_SIZE, 1));
    for (int i = 0; i < BLOCK_SIZE; ++i) {
        for (int j = 0; j < BLOCK_SIZE; ++j) {
            double dx = i - CENTER;
            double dy = j - CENTER;
            if (dx * dx + dy * dy <= RADIUS * RADIUS) {
                holeBlock[i][j] = 0; // Внутри окружности — 0
            }
        }
    }

    // Формирование результата
    int newH = origH * BLOCK_SIZE;
    int newW = origW * BLOCK_SIZE;
    std::vector<std::vector<int>> result(newH, std::vector<int>(newW, 0));

    for (int i = 0; i < origH; ++i) {
        for (int j = 0; j < origW; ++j) {
            const auto& block = (mask[i][j] == 1) ? fullBlock : holeBlock;
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
    
    // Информация о структуре блоков
    int hole_zeros = 0;
    for (const auto& row : holeBlock)
        for (int val : row) if (val == 0) hole_zeros++;
    std::cout << "Блок для '1': полностью заполнен 1 (" << BLOCK_SIZE << "x" << BLOCK_SIZE << " = " << BLOCK_SIZE*BLOCK_SIZE << " ед.)" << std::endl;
    std::cout << "Блок для '0': " << hole_zeros << " нулей внутри окружности (радиус=" << RADIUS << ", центр=" << CENTER << ")" << std::endl;

    return 0;
}