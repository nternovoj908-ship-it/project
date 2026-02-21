#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <cmath>
#include <filesystem>
#include <chrono>
#include <omp.h>

// Геометрические параметры из DetectorConstruction.cc
const double f_base = 17.2; // базовое расстояние от источника до коллиматора
const double d = 4.2;       // расстояние от коллиматора до детектора
const double hole_step = 0.036129; // шаг отверстий в маске (в см)
const double detector_width_mm = 14.08; // ширина детектора (в мм)

// === РАЗМЕРЫ (без изменений) ===
const int DETECTOR_SIZE_Y = 256;
const int DETECTOR_SIZE_X = 256;
const int SOURCE_SIZE_Y = 128;  // размер среза источника
const int SOURCE_SIZE_X = 128;
// ==============================================================================

// ============================================================================
// === НОВАЯ ФУНКЦИЯ: Точечный источник =======================================
// ============================================================================

// Создаёт матрицу источника с единицей в позиции (py, px)
std::vector<std::vector<int>> create_point_source(int py, int px) {
    std::vector<std::vector<int>> source(SOURCE_SIZE_Y, std::vector<int>(SOURCE_SIZE_X, 0));
    if (py >= 0 && py < SOURCE_SIZE_Y && px >= 0 && px < SOURCE_SIZE_X) {
        source[py][px] = 1;
    }
    return source;
}

// ============================================================================
// === ФУНКЦИЯ ГЕНЕРАЦИИ ТЕНЕГРАММЫ (без изменений логики) ====================
// ============================================================================

std::vector<std::vector<double>> generate_tengram_for_source(
    const std::vector<std::vector<int>>& source,
    const std::vector<std::vector<int>>& mask,
    double f_current)
{
    std::vector<std::vector<double>> tengram(DETECTOR_SIZE_Y, std::vector<double>(DETECTOR_SIZE_X, 0.0));
    int MASK_SIZE_Y = mask.size();
    int MASK_SIZE_X = (MASK_SIZE_Y > 0) ? mask[0].size() : 0;

    for (int sy = 0; sy < SOURCE_SIZE_Y; ++sy) {
        for (int sx = 0; sx < SOURCE_SIZE_X; ++sx) {
            if (source[sy][sx] == 0) continue;

            for (int mx = 0; mx < MASK_SIZE_Y; ++mx) {
                for (int my = 0; my < MASK_SIZE_X; ++my) {
                    if (mask[mx][my] == 0) { // 0 = отверстие
                        double mx_centered = 1.7 * (mx - MASK_SIZE_Y / 2);
                        double my_centered = 1.7 * (my - MASK_SIZE_X / 2);

                        int sy_centered = sy - SOURCE_SIZE_Y / 2;
                        int sx_centered = sx - SOURCE_SIZE_X / 2;

                        double proj_y = (-1) * (sy_centered - (sy_centered - 0.125 * mx_centered) * (f_current + d) / f_current);
                        double proj_x = (+1) * (sx_centered - (sx_centered - 0.125 * my_centered) * (f_current + d) / f_current);

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
    return tengram;
}

// ============================================================================
// === ФУНКЦИИ СОХРАНЕНИЯ =====================================================
// ============================================================================

void save_tengram(const std::vector<std::vector<double>>& tengram, const std::string& filename) {
    std::ofstream file(filename + ".txt");
    if (!file.is_open()) return;
    file << std::fixed;
    file.precision(6);
    for (int i = 0; i < DETECTOR_SIZE_Y; ++i) {
        for (int j = 0; j < DETECTOR_SIZE_X; ++j) {
            file << tengram[i][j];
            if (j < DETECTOR_SIZE_X - 1) file << "\t";
        }
        file << "\n";
    }
    file.close();
}

// ============================================================================
// === ОСНОВНАЯ ПРОГРАММА =====================================================
// ============================================================================

int main()
{
    // Настройка OpenMP
    int max_threads = omp_get_max_threads();
    omp_set_num_threads(max_threads);
    omp_set_nested(0);  // Отключаем вложенную параллелизацию для стабильности
    
    std::cout << "=== 3D Point Source Dataset Generation (OpenMP enabled) ===\n";
    std::cout << "Available CPU threads: " << max_threads << "\n";
    std::cout << "Source slice size: " << SOURCE_SIZE_Y << "x" << SOURCE_SIZE_X << "\n";
    std::cout << "Detector size: " << DETECTOR_SIZE_Y << "x" << DETECTOR_SIZE_X << "\n";
    std::cout << "Total teneagrams: " << (3 * SOURCE_SIZE_Y * SOURCE_SIZE_X) << "\n";
    std::cout << "  Slice f (0):   " << SOURCE_SIZE_Y << "x" << SOURCE_SIZE_X << " positions\n";
    std::cout << "  Slice f+1 (1): " << SOURCE_SIZE_Y << "x" << SOURCE_SIZE_X << " positions\n";
    std::cout << "  Slice f+2 (2): " << SOURCE_SIZE_Y << "x" << SOURCE_SIZE_X << " positions\n";
    std::cout << "Base f = " << f_base << ", d = " << d << "\n\n";

    // Создаём папку для вывода
    std::filesystem::create_directories("./3D/tengram_point");
    std::cout << "Output directory: ./3D/tengram_point/\n\n";

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

    // Счётчики
    int total_generated = 0;
    auto global_start = std::chrono::high_resolution_clock::now();

    std::cout << "Starting generation...\n\n";

    // ========================================================================
    // ГЛАВНЫЙ ЦИКЛ: 3 среза × 128×128 позиций = 49,152 тенеграммы
    // ========================================================================
    
    // Массив значений f для каждого среза
    const double f_values[3] = {f_base, f_base + 1, f_base + 2};
    
    for (int slice = 0; slice < 3; ++slice) {
        double f_current = f_values[slice];
        std::cout << "Processing slice " << slice << " (f = " << f_current << ")...\n";
        
        auto slice_start = std::chrono::high_resolution_clock::now();
        
        // === ПАРАЛЛЕЛЬНЫЙ ЦИКЛ ПО ПОЗИЦИЯМ ИСТОЧНИКА ===
        #pragma omp parallel for collapse(2) schedule(dynamic, 64)
        for (int i = 0; i < SOURCE_SIZE_Y; ++i) {
            for (int j = 0; j < SOURCE_SIZE_X; ++j) {
                
                // 1. Создаём точечный источник в позиции (i, j)
                auto point_source = create_point_source(i, j);
                
                // 2. Генерируем тенеграмму
                auto tengram = generate_tengram_for_source(point_source, mask, f_current);
                
                // 3. Формируем имя файла: fig_{i}_{j}_{slice}.txt
                std::string filename = "./3D/tengram_point/fig_" + 
                                       std::to_string(i) + "_" + 
                                       std::to_string(j) + "_" + 
                                       std::to_string(slice);
                
                // 4. Сохраняем тенеграмму
                save_tengram(tengram, filename);
                
                // Атомарное увеличение счётчика
                #pragma omp atomic
                total_generated++;
            }
        }
        
        auto slice_end = std::chrono::high_resolution_clock::now();
        auto slice_time = std::chrono::duration_cast<std::chrono::seconds>(slice_end - slice_start).count();
        
        std::cout << "  Slice " << slice << " completed: " 
                  << (SOURCE_SIZE_Y * SOURCE_SIZE_X) << " teneagrams in " 
                  << slice_time << "s (" 
                  << (static_cast<double>(slice_time) / (SOURCE_SIZE_Y * SOURCE_SIZE_X) * 1000) 
                  << " ms/tenegram)\n\n";
        
        // Прогресс
        auto current = std::chrono::high_resolution_clock::now();
        auto elapsed = std::chrono::duration_cast<std::chrono::seconds>(current - global_start).count();
        double avg = static_cast<double>(elapsed) / total_generated;
        double remaining = avg * (3 * SOURCE_SIZE_Y * SOURCE_SIZE_X - total_generated);
        
        std::cout << "Overall progress: " << total_generated << "/" 
                  << (3 * SOURCE_SIZE_Y * SOURCE_SIZE_X) 
                  << " (" << (100.0 * total_generated / (3 * SOURCE_SIZE_Y * SOURCE_SIZE_X)) << "%)\n";
        std::cout << "  Elapsed: " << elapsed << "s, ETA: " 
                  << static_cast<int>(remaining / 60) << "m " 
                  << static_cast<int>(remaining) % 60 << "s\n\n";
    }
    
    auto global_end = std::chrono::high_resolution_clock::now();
    auto total_time = std::chrono::duration_cast<std::chrono::seconds>(global_end - global_start).count();
    
    std::cout << "\n=== Generation Completed ===\n";
    std::cout << "Total teneagrams generated: " << total_generated << "\n";
    std::cout << "Total time: " << total_time / 60 << "m " << total_time % 60 << "s\n";
    std::cout << "Average time per teneagram: " 
              << (static_cast<double>(total_time) / total_generated * 1000) << " ms\n";
    std::cout << "\nDataset structure:\n";
    std::cout << "  3D/\n";
    std::cout << "  └── tengram_point/\n";
    std::cout << "      ├── fig_0_0_0.txt      (point at [0][0], slice f)\n";
    std::cout << "      ├── fig_0_1_0.txt      (point at [0][1], slice f)\n";
    std::cout << "      ├── ...\n";
    std::cout << "      ├── fig_127_127_0.txt  (point at [127][127], slice f)\n";
    std::cout << "      ├── fig_0_0_1.txt      (point at [0][0], slice f+1)\n";
    std::cout << "      ├── ...\n";
    std::cout << "      └── fig_127_127_2.txt  (point at [127][127], slice f+2)\n";
    
    return 0;
}