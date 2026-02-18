#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <cmath>

// Геометрические параметры из DetectorConstruction.cc
const double f = 17.2; // расстояние от источника до коллиматора
const double d = 4.2;  // расстояние от коллиматора до детектора
const double hole_step = 0.036129; // шаг отверстий в маске (в см)
const double detector_width_mm = 14.08; // ширина детектора (в мм)

// === ИЗМЕНЕНИЕ: Раздельные размеры для детектора и источника ===
// Размеры детектора (не изменяются)
const int DETECTOR_SIZE_Y = 256;
const int DETECTOR_SIZE_X = 256;

// Размеры плоскости источника (уменьшены в 2 раза по каждой оси = в 4 раза по площади)
const int SOURCE_SIZE_Y = 128;  // 256 / 2
const int SOURCE_SIZE_X = 128;  // 256 / 2
// ==============================================================================

int main()
{
    // Параметры кругового источника
    const double source_radius = 10.0; // Радиус круга в индексах источника
    // === ИЗМЕНЕНИЕ: Центр вычисляется относительно SOURCE_SIZE ===
    const int center_source_y = SOURCE_SIZE_Y / 2; // Центр круга Y (в индексах источника)
    const int center_source_x = SOURCE_SIZE_X / 2; // Центр круга X (в индексах источника)
    // ==============================================================================

    std::cout << "Creating tengram for a circular source.\n";
    std::cout << "Detector Size: " << DETECTOR_SIZE_Y << " x " << DETECTOR_SIZE_X << "\n";
    std::cout << "Source Plane Size: " << SOURCE_SIZE_Y << " x " << SOURCE_SIZE_X << "\n";
    std::cout << "Circular source center: (" << center_source_y << ", " << center_source_x << ")\n";
    std::cout << "Circular source radius: " << source_radius << " (in source indices)\n";

    // --- Шаг 1: Загрузка маски (динамически определяем размеры) ---
    //std::ifstream mask_file("C:/Users/User/c++1/mask-mura-expanded.txt");
    std::ifstream mask_file("./mask-mura-expanded_16.txt");
    if (!mask_file.is_open()) {
        std::cerr << "Cannot open mask-mura-expanded.txt\n";
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
        std::cerr << "Loaded mask file is empty\n";
        return 1;
    }
    std::cout << "Mask loaded successfully. Size: " << MASK_SIZE_Y << " x " << MASK_SIZE_X << "\n";

    // Подсчитаем отверстия (0) и закрытые (1) для информации
    int holes = 0, closed = 0;
    for (int i = 0; i < MASK_SIZE_Y; ++i) {
        for (int j = 0; j < MASK_SIZE_X; ++j) {
            if (mask[i][j] == 0) holes++;
            else closed++;
        }
    }
    std::cout << "Mask stats: " << holes << " holes (0), " << closed << " closed (1)\n";

    // === ДОБАВЛЕНО: Массив для хранения истинного источника ===
    std::vector<std::vector<int>> true_source(SOURCE_SIZE_Y, std::vector<int>(SOURCE_SIZE_X, 0));
    // ===========================================================

    // --- Шаг 2: Создание тенеграммы ---
    // Инициализация тенеграммы нулями (размер детектора — НЕ изменяется)
    std::vector<std::vector<double>> generated_tengram(DETECTOR_SIZE_Y, std::vector<double>(DETECTOR_SIZE_X, 0.0));

    // Подсчитаем количество точек внутри круга для информации
    int points_in_circle = 0;
    // === ИЗМЕНЕНИЕ: Итерация по плоскости источника (SOURCE_SIZE) ===
    for (int sy = 0; sy < SOURCE_SIZE_Y; ++sy) {
        for (int sx = 0; sx < SOURCE_SIZE_X; ++sx) {
            // Проверяем, находится ли точка (sy, sx) внутри круга
            double dy = sy - center_source_y;
            double dx = sx - center_source_x;
            double distance_squared = dy * dy + dx * dx;
            if (distance_squared <= source_radius * source_radius) {
                points_in_circle++;
                // === ДОБАВЛЕНО: Сохраняем точку в истинный источник ===
                true_source[sy][sx] = 1;
                // =====================================================

                // --- Проекция этой точки на детектор ---
                // Используем ТОЧНУЮ копию логики из реконструкции (эталон первых 3 тенеграмм)
                for (int mx = 0; mx < MASK_SIZE_Y; ++mx) {
                    for (int my = 0; my < MASK_SIZE_X; ++my) {
                        if (mask[mx][my] == 0) { // 0 = отверстие (дырка) — ЭТАЛОННОЕ определение
                            // Координаты относительно центра маски (в единицах шага) — ЭТАЛОННЫЙ масштаб 1.7
                            double mx_centered = 1.7 * (mx - MASK_SIZE_Y / 2);
                            double my_centered = 1.7 * (my - MASK_SIZE_X / 2);

                            // === ИЗМЕНЕНИЕ: Координаты источника центрируются относительно SOURCE_SIZE ===
                            int sy_centered = sy - SOURCE_SIZE_Y / 2; // Y источника
                            int sx_centered = sx - SOURCE_SIZE_X / 2; // X источника
                            // ==============================================================================

                            // === Формула проекции (без изменений, использует DETECTOR_SIZE для выхода) ===
                            double proj_y = (-1) * (sy_centered - (sy_centered - 0.25 * mx_centered) * (f + d) / f);
                            double proj_x = (+1) * (sx_centered - (sx_centered - 0.25 * my_centered) * (f + d) / f);

                            // Преобразуем в целые координаты детектора — смещение по DETECTOR_SIZE (без изменений)
                            int det_y = static_cast<int>(round(proj_y + DETECTOR_SIZE_Y / 2));
                            int det_x = static_cast<int>(round(proj_x + DETECTOR_SIZE_X / 2));
                            // ==============================================================================

                            // Проверяем, входят ли координаты в размер детектора
                            if (det_y >= 0 && det_y < DETECTOR_SIZE_Y && det_x >= 0 && det_x < DETECTOR_SIZE_X) {
                                generated_tengram[det_y][det_x] += 1.0;
                            }
                        }
                    }
                }
            }
        }
    }
    // ==============================================================================

    std::cout << "Number of points inside the circular source: " << points_in_circle << "\n";
    std::cout << "Tengram generation completed.\n";

    // === ДОБАВЛЕНО: Сохранение истинного источника в файлы ===
    // 1. Сохранение в .txt
    std::ofstream source_txt("./true_source.txt");
    if (source_txt.is_open()) {
        for (int i = 0; i < SOURCE_SIZE_Y; ++i) {
            for (int j = 0; j < SOURCE_SIZE_X; ++j) {
                source_txt << true_source[i][j];
                if (j < SOURCE_SIZE_X - 1) source_txt << "\t";
            }
            source_txt << "\n";
        }
        source_txt.close();
        std::cout << "True source saved to true_source.txt\n";
    } else {
        std::cerr << "Warning: Cannot open true_source.txt for writing\n";
    }

    // 2. Сохранение в .pgm
    std::ofstream source_pgm("./true_source.pgm");
    if (source_pgm.is_open()) {
        source_pgm << "P2\n" << SOURCE_SIZE_X << " " << SOURCE_SIZE_Y << "\n255\n";
        for (int i = 0; i < SOURCE_SIZE_Y; ++i) {
            for (int j = 0; j < SOURCE_SIZE_X; ++j) {
                unsigned char gray_val = (true_source[i][j] > 0) ? 255 : 0;
                source_pgm << static_cast<int>(gray_val) << " ";
            }
            source_pgm << "\n";
        }
        source_pgm.close();
        std::cout << "True source saved to true_source.pgm\n";
    } else {
        std::cerr << "Warning: Cannot open true_source.pgm for writing\n";
    }
    // =============================================================

    // Подсчитаем и выведем статистику по полученной тенеграмме (размер детектора)
    int nonzero_count_gen = 0;
    double max_val_gen = 0.0;
    for (int i = 0; i < DETECTOR_SIZE_Y; ++i) {
        for (int j = 0; j < DETECTOR_SIZE_X; ++j) {
            if (generated_tengram[i][j] != 0) {
                nonzero_count_gen++;
                if (generated_tengram[i][j] > max_val_gen) {
                    max_val_gen = generated_tengram[i][j];
                }
            }
        }
    }
    std::cout << "Generated tengram stats:\n";
    std::cout << "  Non-zero elements: " << nonzero_count_gen << "\n";
    std::cout << "  Max value: " << max_val_gen << "\n";

    // --- Шаг 3: Сохранение тенеграммы (размер детектора — без изменений) ---
    // 3.1. Сохранение в .txt
    std::ofstream tengram_txt_file("./generated_circular_tengram.txt");
    if (!tengram_txt_file.is_open()) {
        std::cerr << "Cannot open generated_circular_tengram.txt for writing.\n";
        return 1;
    }
    for (int i = 0; i < DETECTOR_SIZE_Y; ++i) {
        for (int j = 0; j < DETECTOR_SIZE_X; ++j) {
            tengram_txt_file << generated_tengram[i][j];
            if (j < DETECTOR_SIZE_X - 1) tengram_txt_file << "\t";
        }
        tengram_txt_file << "\n";
    }
    tengram_txt_file.close();
    std::cout << "Generated tengram saved to generated_circular_tengram.txt\n";

    // 3.2. Сохранение в .pgm
    std::ofstream tengram_pgm_file("./generated_circular_tengram.pgm");
    if (!tengram_pgm_file.is_open()) {
        std::cerr << "Cannot open generated_circular_tengram.pgm for writing.\n";
        return 1;
    }
    tengram_pgm_file << "P2\n" << DETECTOR_SIZE_X << " " << DETECTOR_SIZE_Y << "\n255\n";
    for (int i = 0; i < DETECTOR_SIZE_Y; ++i) {
        for (int j = 0; j < DETECTOR_SIZE_X; ++j) {
            unsigned char gray_val = (generated_tengram[i][j] > 0) ? 255 : 0;
            tengram_pgm_file << static_cast<int>(gray_val) << " ";
        }
        tengram_pgm_file << "\n";
    }
    tengram_pgm_file.close();
    std::cout << "Generated tengram saved to generated_circular_tengram.pgm\n";

    std::cout << "Program finished successfully.\n";
    return 0;
}
