#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <cmath>
#include <filesystem>  // Для создания папок

// Геометрические параметры из DetectorConstruction.cc
const double f = 17.2; // расстояние от источника до коллиматора
const double d = 4.2;  // расстояние от коллиматора до детектора
const double hole_step = 0.036129; // шаг отверстий в маске (в см)
const double detector_width_mm = 14.08; // ширина детектора (в мм)

// === Размеры (без изменений) ===
const int DETECTOR_SIZE_Y = 256;
const int DETECTOR_SIZE_X = 256;
const int SOURCE_SIZE_Y = 128;
const int SOURCE_SIZE_X = 128;
// ==============================================================================

// === ПАРАМЕТРЫ ГЕНЕРАЦИИ ДАТАСЕТА ===
const double source_radius = 10.0;      // Радиус круга (в индексах источника)
const int GRID_STEP = 10;               // Шаг сетки для перемещения источника
const int MARGIN = 15;                  // Отступ от края (чтобы круг не обрезался)
// ==============================================================================

// Функция генерации тенеграммы для источника в заданной позиции
std::vector<std::vector<double>> generate_tengram_at_position(
    int center_y, int center_x,
    const std::vector<std::vector<int>>& mask,
    int MASK_SIZE_Y, int MASK_SIZE_X)
{
    std::vector<std::vector<double>> tengram(DETECTOR_SIZE_Y, std::vector<double>(DETECTOR_SIZE_X, 0.0));
    
    for (int sy = 0; sy < SOURCE_SIZE_Y; ++sy) {
        for (int sx = 0; sx < SOURCE_SIZE_X; ++sx) {
            double dy = sy - center_y;
            double dx = sx - center_x;
            double distance_squared = dy * dy + dx * dx;
            
            if (distance_squared <= source_radius * source_radius) {
                for (int mx = 0; mx < MASK_SIZE_Y; ++mx) {
                    for (int my = 0; my < MASK_SIZE_X; ++my) {
                        if (mask[mx][my] == 0) {
                            double mx_centered = 1.7 * (mx - MASK_SIZE_Y / 2);
                            double my_centered = 1.7 * (my - MASK_SIZE_X / 2);
                            
                            int sy_centered = sy - SOURCE_SIZE_Y / 2;
                            int sx_centered = sx - SOURCE_SIZE_X / 2;
                            
                            double proj_y = (-1) * (sy_centered - (sy_centered - 0.25 * mx_centered) * (f + d) / f);
                            double proj_x = (+1) * (sx_centered - (sx_centered - 0.25 * my_centered) * (f + d) / f);
                            
                            int det_y = static_cast<int>(round(proj_y + DETECTOR_SIZE_Y / 2));
                            int det_x = static_cast<int>(round(proj_x + DETECTOR_SIZE_X / 2));
                            
                            if (det_y >= 0 && det_y < DETECTOR_SIZE_Y && det_x >= 0 && det_x < DETECTOR_SIZE_X) {
                                tengram[det_y][det_x] += 1.0;
                            }
                        }
                    }
                }
            }
        }
    }
    return tengram;
}

// Функция сохранения источника в файл
void save_source(const std::vector<std::vector<int>>& source, const std::string& filename) {
    std::ofstream file(filename);
    if (!file.is_open()) return;
    for (int i = 0; i < SOURCE_SIZE_Y; ++i) {
        for (int j = 0; j < SOURCE_SIZE_X; ++j) {
            file << source[i][j];
            if (j < SOURCE_SIZE_X - 1) file << "\t";
        }
        file << "\n";
    }
    file.close();
}

// Функция сохранения тенеграммы в файл
void save_tengram(const std::vector<std::vector<double>>& tengram, const std::string& filename) {
    std::ofstream file(filename);
    if (!file.is_open()) return;
    for (int i = 0; i < DETECTOR_SIZE_Y; ++i) {
        for (int j = 0; j < DETECTOR_SIZE_X; ++j) {
            file << tengram[i][j];
            if (j < DETECTOR_SIZE_X - 1) file << "\t";
        }
        file << "\n";
    }
    file.close();
}

int main()
{
    // Создаём папки для вывода
    std::filesystem::create_directories("./output/sources");
    std::filesystem::create_directories("./output/tenegrams");
    
    std::cout << "=== Generating dataset ===\n";
    std::cout << "Source plane: " << SOURCE_SIZE_Y << "x" << SOURCE_SIZE_X << "\n";
    std::cout << "Detector size: " << DETECTOR_SIZE_Y << "x" << DETECTOR_SIZE_X << "\n";
    std::cout << "Circle radius: " << source_radius << ", Grid step: " << GRID_STEP << "\n";
    std::cout << "Margin from edge: " << MARGIN << "\n\n";

    // Загрузка маски
    std::ifstream mask_file("./mask-mura-expanded_16.txt");
    if (!mask_file.is_open()) {
        std::cerr << "Cannot open mask-mura-expanded_16.txt\n";
        return 1;
    }
    std::vector<std::vector<int>> mask;
    std::string line;
    while (std::getline(mask_file, line)) {
        std::istringstream iss(line);
        std::vector<int> row;
        int val;
        while (iss >> val) row.push_back(val);
        mask.push_back(row);
    }
    mask_file.close();

    int MASK_SIZE_Y = mask.size();
    int MASK_SIZE_X = (MASK_SIZE_Y > 0) ? mask[0].size() : 0;
    if (MASK_SIZE_X == 0 || MASK_SIZE_Y == 0) {
        std::cerr << "Loaded mask file is empty\n";
        return 1;
    }
    std::cout << "Mask loaded: " << MASK_SIZE_Y << " x " << MASK_SIZE_X << "\n\n";

    // Генерация пар источник-тенеграмма
    int pair_index = 0;
    
    // Проходим по сетке позиций для центра круга
    for (int cy = MARGIN; cy < SOURCE_SIZE_Y - MARGIN; cy += GRID_STEP) {
        for (int cx = MARGIN; cx < SOURCE_SIZE_X - MARGIN; cx += GRID_STEP) {
            
            // 1. Создаём источник-круг в текущей позиции
            std::vector<std::vector<int>> true_source(SOURCE_SIZE_Y, std::vector<int>(SOURCE_SIZE_X, 0));
            for (int sy = 0; sy < SOURCE_SIZE_Y; ++sy) {
                for (int sx = 0; sx < SOURCE_SIZE_X; ++sx) {
                    double dy = sy - cy;
                    double dx = sx - cx;
                    if (dy*dy + dx*dx <= source_radius * source_radius) {
                        true_source[sy][sx] = 1;
                    }
                }
            }
            
            // 2. Генерируем тенеграмму для этого источника
            auto tengram = generate_tengram_at_position(cy, cx, mask, MASK_SIZE_Y, MASK_SIZE_X);
            
            // 3. Сохраняем пару файлов
            std::string source_file = "./output/sources/circus_2D_" + std::to_string(pair_index) + ".txt";
            std::string tengram_file = "./output/tenegrams/" + std::to_string(pair_index) + ".txt";
            
            save_source(true_source, source_file);
            save_tengram(tengram, tengram_file);
            
            pair_index++;
            
            // Прогресс
            if (pair_index % 20 == 0) {
                std::cout << "  Generated " << pair_index << " pairs...\n";
            }
        }
    }
    
    std::cout << "\n=== Dataset generation completed ===\n";
    std::cout << "Total pairs generated: " << pair_index << "\n";
    std::cout << "Sources saved to: ./output/sources/circus_2D_*.txt\n";
    std::cout << "Teneagrams saved to: ./output/tenegrams/*.txt\n";
    
    return 0;
}