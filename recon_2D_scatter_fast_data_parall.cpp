#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <cmath>
#include <filesystem>  // Для работы с директориями
#include <omp.h>       // <<<--- ДОБАВЛЕНО: OpenMP для параллелизации

// Геометрические параметры из DetectorConstruction.cc
const double f = 17.2; // расстояние от источника до коллиматора (изм. строка 6)
const double d = 4.2;  // расстояние от коллиматора до детектора (изм. строка 7)
const double hole_step = 0.036129; // шаг отверстий в маске (в см)
const double detector_width_mm = 14.08; // ширина детектора (в мм)
// Порог совпадения (оставлен для совместимости, но не используется в новой логике)
const double MATCH_THRESHOLD = 0.15; // 80%

// === Размеры (без изменений) ===
const int DETECTOR_SIZE_Y = 256;
const int DETECTOR_SIZE_X = 256;
const int SOURCE_SIZE_Y = 128;  // Размер плоскости реконструкции
const int SOURCE_SIZE_X = 128;
// ==============================================================================

int main()
{
    // === OpenMP: Информация о потоках ===
    int max_threads = omp_get_max_threads();
    std::cout << "=== Batch reconstruction (OpenMP enabled) ===\n";
    std::cout << "Available CPU threads: " << max_threads << "\n";
    std::cout << "Detector size: " << DETECTOR_SIZE_Y << "x" << DETECTOR_SIZE_X << "\n";
    std::cout << "Source reconstruction grid: " << SOURCE_SIZE_Y << "x" << SOURCE_SIZE_X << "\n";
    std::cout << "Input tengrams: ./output/tenegrams/{i}.txt\n";
    std::cout << "Output reconstructions: ./output/reconstruction/{i}.txt\n\n";

    // Создаём папку для результатов реконструкции
    std::filesystem::create_directories("./output/reconstruction");

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
        while (iss >> val) {
            row.push_back(val);
        }
        mask.push_back(row);
    }
    mask_file.close();

    int MASK_SIZE_Y = mask.size();
    int MASK_SIZE_X = (MASK_SIZE_Y > 0) ? mask[0].size() : 0;
    if (MASK_SIZE_X == 0 || MASK_SIZE_Y == 0) {
        std::cerr << "Empty mask file\n";
        return 1;
    }
    std::cout << "Mask loaded: " << MASK_SIZE_Y << " x " << MASK_SIZE_X << "\n\n";

    // === Цикл обработки всех тенеграмм ===
    int processed_count = 0;
    int max_index = 1000;  // Ожидаемое максимальное количество файлов
    
    for (int idx = 0; idx < max_index; ++idx) {
        // Формируем путь к файлу тенеграммы
        std::string tengram_path = "./output/tenegrams/" + std::to_string(idx) + ".txt";
        
        // Проверяем существование файла
        if (!std::filesystem::exists(tengram_path)) {
            // Если файл не найден — возможно, это конец датасета
            if (processed_count > 0) {
                std::cout << "No more tengram files found after index " << idx << "\n";
            }
            break;
        }
        
        // === Чтение тенеграммы (логика без изменений) ===
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
            while (iss >> val) {
                row.push_back(val);
            }
            real_tengram.push_back(row);
        }
        tengram_file.close();

        // Проверка размера
        int file_det_y = real_tengram.size();
        int file_det_x = (file_det_y > 0) ? real_tengram[0].size() : 0;
        if (file_det_x != DETECTOR_SIZE_X || file_det_y != DETECTOR_SIZE_Y) {
            std::cerr << "Warning: Tengram " << idx << " size mismatch (" 
                      << file_det_y << "x" << file_det_x << ")\n";
            continue;
        }

        // === Массив реконструкции (размер SOURCE_SIZE) ===
        std::vector<std::vector<double>> source_reconstruction(SOURCE_SIZE_Y, std::vector<double>(SOURCE_SIZE_X, 0.0));

        // === ПАРАЛЛЕЛЬНЫЙ ЦИКЛ РЕКОНСТРУКЦИИ ===
        #pragma omp parallel for collapse(2) schedule(dynamic)
        for (int sx = 0; sx < SOURCE_SIZE_Y; ++sx) {
            for (int sy = 0; sy < SOURCE_SIZE_X; ++sy) {
                // Создаём вспомогательную тенеграмму (локальная для каждого потока)
                std::vector<std::vector<int>> aux_tengram(DETECTOR_SIZE_Y, std::vector<int>(DETECTOR_SIZE_X, 0));

                // Проходим по всем отверстиям маски
                for (int mx = 0; mx < MASK_SIZE_Y; ++mx) {
                    for (int my = 0; my < MASK_SIZE_X; ++my) {
                        if (mask[mx][my] == 0) { // 0 = дырка
                            // Координаты относительно центра маски
                            int mx_centered = 1.7*(mx - MASK_SIZE_Y / 2);
                            int my_centered = 1.7*(my - MASK_SIZE_X / 2);

                            // Координаты источника центрируются относительно SOURCE_SIZE
                            int sx_centered = sx - SOURCE_SIZE_Y / 2;
                            int sy_centered = sy - SOURCE_SIZE_X / 2;

                            // Формула проекции (без изменений)
                            double dx_real = (-1)*(sx_centered - (sx_centered - 0.25*mx_centered) * (f + d) / f);
                            double dy_real = (+1)*(sy_centered - (sy_centered - 0.25*my_centered) * (f + d) / f);

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

                // === Логика соответствия (без изменений) ===
                for (int i = 0; i < DETECTOR_SIZE_Y; ++i) {
                    for (int j = 0; j < DETECTOR_SIZE_X; ++j) {
                        if (aux_tengram[i][j] > 0 && real_tengram[i][j] > 0) {
                            source_reconstruction[sx][sy] += real_tengram[i][j];
                        }
                    }
                }
            }
        }
        // === Конец параллельного цикла ===

        // === Сохранение результата в папку reconstruction ===
        std::string output_path = "./output/reconstruction/" + std::to_string(idx) + ".txt";
        std::ofstream out(output_path);
        if (!out.is_open()) {
            std::cerr << "Cannot open " << output_path << " for writing\n";
            continue;
        }
        for (int i = 0; i < SOURCE_SIZE_Y; ++i) {
            for (int j = 0; j < SOURCE_SIZE_X; ++j) {
                out << source_reconstruction[i][j];
                if (j < SOURCE_SIZE_X - 1) out << "\t";
            }
            out << "\n";
        }
        out.close();
        
        processed_count++;
        
        // Прогресс (защищён от гонок вывода)
        if (processed_count % 50 == 0) {
            #pragma omp critical
            {
                std::cout << "  Processed " << processed_count << " reconstructions...\n";
            }
        }
    }
    
    std::cout << "\n=== Batch reconstruction completed ===\n";
    std::cout << "Total reconstructions saved: " << processed_count << "\n";
    std::cout << "Results saved to: ./output/reconstruction/{i}.txt\n";
    
    return 0;
}
