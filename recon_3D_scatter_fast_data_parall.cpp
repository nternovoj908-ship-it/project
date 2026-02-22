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

// === РАЗМЕРЫ: отдельные для источника и детектора ===
const int DETECTOR_SIZE_Y = 256;  // размер детектора (не менять)
const int DETECTOR_SIZE_X = 256;

const int SOURCE_SIZE_Y = 128;    // размер плоскости реконструкции (источника)
const int SOURCE_SIZE_X = 128;
// ==============================================================================

// ============================================================================
// === ФУНКЦИЯ РЕКОНСТРУКЦИИ ОДНОГО СРЕЗА (с OpenMP и новым принципом) =======
// ============================================================================

void reconstruct_slice(
    const std::vector<std::vector<double>>& real_tengram,
    const std::vector<std::vector<int>>& mask,
    double f_current,
    std::vector<std::vector<double>>& source_reconstruction)
{
    int MASK_SIZE_Y = mask.size();
    int MASK_SIZE_X = (MASK_SIZE_Y > 0) ? mask[0].size() : 0;

    // === ПАРАЛЛЕЛЬНЫЙ ЦИКЛ РЕКОНСТРУКЦИИ ===
    #pragma omp parallel for collapse(2) schedule(dynamic)
    for (int sx = 0; sx < SOURCE_SIZE_Y; ++sx) {
        for (int sy = 0; sy < SOURCE_SIZE_X; ++sy) {
            
            // Локальная вспомогательная тенеграмма (без сохранения)
            std::vector<std::vector<int>> aux_tengram(DETECTOR_SIZE_Y, std::vector<int>(DETECTOR_SIZE_X, 0));

            // Проходим по всем отверстиям маски (0 = дырка)
            for (int mx = 0; mx < MASK_SIZE_Y; ++mx) {
                for (int my = 0; my < MASK_SIZE_X; ++my) {
                    if (mask[mx][my] == 0) {
                        
                        // Координаты относительно центра маски
                        double mx_centered = 1.7 * (mx - MASK_SIZE_Y / 2);
                        double my_centered = 1.7 * (my - MASK_SIZE_X / 2);

                        // Координаты источника центрируются относительно SOURCE_SIZE
                        int sx_centered = sx - SOURCE_SIZE_Y / 2;
                        int sy_centered = sy - SOURCE_SIZE_X / 2;

                        // Формула проекции (без изменений: коэффициенты 0.125, (-1), (+1))
                        double dx_real = (-1) * (sx_centered - (sx_centered - 0.125 * mx_centered) * (f_current + d) / f_current);
                        double dy_real = (+1) * (sy_centered - (sy_centered - 0.125 * my_centered) * (f_current + d) / f_current);

                        // Преобразуем в целые координаты детектора
                        int dx = static_cast<int>(round(dx_real + DETECTOR_SIZE_Y / 2));
                        int dy = static_cast<int>(round(dy_real + DETECTOR_SIZE_X / 2));

                        // Проверяем границы детектора
                        if (dx >= 0 && dx < DETECTOR_SIZE_Y && dy >= 0 && dy < DETECTOR_SIZE_X) {
                            aux_tengram[dx][dy]++;
                        }
                    }
                }
            }

            // === НОВЫЙ ПРИНЦИП СООТВЕТСТВИЯ: накопление значений ===
            for (int i = 0; i < DETECTOR_SIZE_Y; ++i) {
                for (int j = 0; j < DETECTOR_SIZE_X; ++j) {
                    if (aux_tengram[i][j] > 0 && real_tengram[i][j] > 0) {
                        // Добавляем ЗНАЧЕНИЕ из реальной тенеграммы к источнику
                        source_reconstruction[sx][sy] += real_tengram[i][j];
                    }
                }
            }
        }
    }
}

// ============================================================================
// === ФУНКЦИЯ СОХРАНЕНИЯ СРЕЗА ===============================================
// ============================================================================

void save_slice(const std::vector<std::vector<double>>& source, const std::string& filename) {
    std::ofstream file(filename + ".txt");
    if (!file.is_open()) return;
    
    file << std::fixed;
    file.precision(2);  // <<< precision(2) как просили
    
    for (int i = 0; i < SOURCE_SIZE_Y; ++i) {
        for (int j = 0; j < SOURCE_SIZE_X; ++j) {
            file << source[i][j];
            if (j < SOURCE_SIZE_X - 1) file << "\t";
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
    // === OpenMP: Информация о потоках ===
    int max_threads = omp_get_max_threads();
    omp_set_num_threads(max_threads);
    
    std::cout << "=== 3D Batch Reconstruction (OpenMP enabled) ===\n";
    std::cout << "Available CPU threads: " << max_threads << "\n";
    std::cout << "Source slice size: " << SOURCE_SIZE_Y << "x" << SOURCE_SIZE_X << "\n";
    std::cout << "Detector size: " << DETECTOR_SIZE_Y << "x" << DETECTOR_SIZE_X << "\n";
    std::cout << "Input tengrams: ./3D/tengram/fig_*.txt\n";
    std::cout << "Output: ./3D/reconstruction/{i}/fig_i_slice_{0,1,2}.txt\n\n";

    // Создаём папку для результатов
    std::filesystem::create_directories("./3D/reconstruction");
    
    // Загрузка маски (одна для всех реконструкций)
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

    // === Цикл обработки всех тенеграмм (~3000 файлов) ===
    int processed_count = 0;
    int max_index = 3000;  // Ожидаемое количество файлов
    auto total_start = std::chrono::high_resolution_clock::now();
    
    for (int idx = 227; idx < max_index; ++idx) {
        // Формируем путь к файлу тенеграммы
        std::string tengram_path = "./3D/tengram/fig_" + std::to_string(idx) + ".txt";
        
        // Проверяем существование файла
        if (!std::filesystem::exists(tengram_path)) {
            if (processed_count > 0) {
                std::cout << "No more tengram files found after index " << idx << "\n";
            }
            break;
        }
        
        // === Чтение тенеграммы ===
        std::ifstream tengram_file(tengram_path);
        if (!tengram_file.is_open()) {
            std::cerr << "Cannot open " << tengram_path << "\n";
            continue;
        }
        std::vector<std::vector<double>> real_tengram;
        while (std::getline(tengram_file, line)) {
            std::istringstream iss(line);
            std::vector<double> row;
            double val;
            while (iss >> val) row.push_back(val);
            real_tengram.push_back(row);
        }
        tengram_file.close();

        // Проверка размера детектора
        int file_det_y = real_tengram.size();
        int file_det_x = (file_det_y > 0) ? real_tengram[0].size() : 0;
        if (file_det_x != DETECTOR_SIZE_X || file_det_y != DETECTOR_SIZE_Y) {
            std::cerr << "Warning: Tengram " << idx << " size mismatch (" 
                      << file_det_y << "x" << file_det_x << ")\n";
            continue;
        }

        // === Создаём папку для результатов этого объёмного источника ===
        std::string output_dir = "./3D/reconstruction/" + std::to_string(idx);
        std::filesystem::create_directories(output_dir);
        
        // === Реконструкция трёх срезов (f, f+1, f+2) ===
        std::vector<std::vector<double>> recon_slice_0(SOURCE_SIZE_Y, std::vector<double>(SOURCE_SIZE_X, 0.0));
        std::vector<std::vector<double>> recon_slice_1(SOURCE_SIZE_Y, std::vector<double>(SOURCE_SIZE_X, 0.0));
        std::vector<std::vector<double>> recon_slice_2(SOURCE_SIZE_Y, std::vector<double>(SOURCE_SIZE_X, 0.0));
        
        auto slice_start = std::chrono::high_resolution_clock::now();
        
        // Срез 0: f = f_base
        reconstruct_slice(real_tengram, mask, f_base, recon_slice_0);
        
        // Срез 1: f = f_base + 1
        reconstruct_slice(real_tengram, mask, f_base + 1, recon_slice_1);
        
        // Срез 2: f = f_base + 2
        reconstruct_slice(real_tengram, mask, f_base + 2, recon_slice_2);
        
        auto slice_end = std::chrono::high_resolution_clock::now();
        auto slice_time = std::chrono::duration_cast<std::chrono::milliseconds>(slice_end - slice_start).count();
        
        // === Сохранение трёх срезов ===
        std::string base_name = output_dir + "/fig_" + std::to_string(idx) + "_slice_";
        save_slice(recon_slice_0, base_name + "0");  // f
        save_slice(recon_slice_1, base_name + "1");  // f+1
        save_slice(recon_slice_2, base_name + "2");  // f+2
        
        processed_count++;
        
        // Прогресс каждые 100 итераций
        if (processed_count % 100 == 0) {
            auto current_time = std::chrono::high_resolution_clock::now();
            auto elapsed = std::chrono::duration_cast<std::chrono::seconds>(current_time - total_start).count();
            double avg_time = static_cast<double>(elapsed) / processed_count;
            double remaining = avg_time * (max_index - processed_count);
            
            std::cout << "Progress: " << processed_count << "/" << max_index 
                      << " (" << (100.0 * processed_count / max_index) << "%)\n";
            std::cout << "  Slice time: " << slice_time << "ms, Avg: " << avg_time << "s/iter\n";
            std::cout << "  ETA: " << static_cast<int>(remaining / 60) << "m " 
                      << static_cast<int>(remaining) % 60 << "s\n\n";
        }
    }
    
    auto total_end = std::chrono::high_resolution_clock::now();
    auto total_time = std::chrono::duration_cast<std::chrono::seconds>(total_end - total_start).count();
    
    std::cout << "\n=== 3D Reconstruction Completed ===\n";
    std::cout << "Total volumes processed: " << processed_count << "\n";
    std::cout << "Total time: " << total_time / 60 << "m " << total_time % 60 << "s\n";
    std::cout << "Average time per volume: " << (static_cast<double>(total_time) / processed_count) << "s\n";
    std::cout << "\nDataset structure:\n";
    std::cout << "  3D/\n";
    std::cout << "  ├── tengram/           (input: fig_0.txt ... fig_" << (processed_count-1) << ".txt)\n";
    std::cout << "  └── reconstruction/\n";
    if (processed_count > 0) {
        std::cout << "      ├── 0/\n";
        std::cout << "      │   ├── fig_0_slice_0.txt  (f=" << f_base << ")\n";
        std::cout << "      │   ├── fig_0_slice_1.txt  (f=" << (f_base+1) << ")\n";
        std::cout << "      │   └── fig_0_slice_2.txt  (f=" << (f_base+2) << ")\n";
        std::cout << "      ├── 1/\n";
        std::cout << "      └── ...\n";
    }
    
    return 0;
}
