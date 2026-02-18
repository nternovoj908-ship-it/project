#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <cmath>

// Геометрические параметры из DetectorConstruction.cc
const double f_base = 17.2; // базовое расстояние от источника до коллиматора
const double d = 4.2;       // расстояние от коллиматора до детектора
const double hole_step = 0.036129; // шаг отверстий в маске (в см)
const double detector_width_mm = 14.08; // ширина детектора (в мм)

// Размеры детектора (желательно загружать из файла, но для примера зададим жёстко)
const int DETECTOR_SIZE_Y = 256;
const int DETECTOR_SIZE_X = 256;

// ============================================================================
// === ФУНКЦИИ (без изменений, кроме добавления save_source) ==================
// ============================================================================

// --- Функция генерации тенеграммы для заданного источника и значения f ---
std::vector<std::vector<double>> generate_tengram_for_source(
    const std::vector<std::vector<int>>& source,
    const std::vector<std::vector<int>>& mask,
    double f_current)
{
    std::vector<std::vector<double>> tengram(DETECTOR_SIZE_Y, std::vector<double>(DETECTOR_SIZE_X, 0.0));
    int MASK_SIZE_Y = mask.size();
    int MASK_SIZE_X = (MASK_SIZE_Y > 0) ? mask[0].size() : 0;

    for (int sy = 0; sy < DETECTOR_SIZE_Y; ++sy) {
        for (int sx = 0; sx < DETECTOR_SIZE_X; ++sx) {
            if (source[sy][sx] == 0) continue; // Пропускаем пустые точки источника

            for (int mx = 0; mx < MASK_SIZE_Y; ++mx) {
                for (int my = 0; my < MASK_SIZE_X; ++my) {
                    if (mask[mx][my] == 0) { // 0 = отверстие (дырка)
                        double mx_centered = 1.7 * (mx - MASK_SIZE_Y / 2);
                        double my_centered = 1.7 * (my - MASK_SIZE_X / 2);

                        int sy_centered = sy - DETECTOR_SIZE_Y / 2;
                        int sx_centered = sx - DETECTOR_SIZE_X / 2;

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

// --- Функция создания кругового источника в центре ---
std::vector<std::vector<int>> create_circle_source(double radius, int center_y, int center_x) {
    std::vector<std::vector<int>> source(DETECTOR_SIZE_Y, std::vector<int>(DETECTOR_SIZE_X, 0));
    for (int sy = 0; sy < DETECTOR_SIZE_Y; ++sy) {
        for (int sx = 0; sx < DETECTOR_SIZE_X; ++sx) {
            double dy = sy - center_y;
            double dx = sx - center_x;
            if (dy*dy + dx*dx <= radius*radius) {
                source[sy][sx] = 1;
            }
        }
    }
    return source;
}

// --- Функция создания квадратного источника в левом верхнем углу ---
std::vector<std::vector<int>> create_square_source(int side, int start_y, int start_x) {
    std::vector<std::vector<int>> source(DETECTOR_SIZE_Y, std::vector<int>(DETECTOR_SIZE_X, 0));
    for (int sy = start_y; sy < start_y + side && sy < DETECTOR_SIZE_Y; ++sy) {
        for (int sx = start_x; sx < start_x + side && sx < DETECTOR_SIZE_X; ++sx) {
            source[sy][sx] = 1;
        }
    }
    return source;
}

// --- Функция создания треугольного источника в правом нижнем углу ---
// Прямоугольный треугольник, катеты направлены вверх и влево от угла
std::vector<std::vector<int>> create_triangle_source(int leg, int corner_y, int corner_x) {
    std::vector<std::vector<int>> source(DETECTOR_SIZE_Y, std::vector<int>(DETECTOR_SIZE_X, 0));
    for (int dy = 0; dy < leg; ++dy) {
        for (int dx = 0; dx < leg - dy; ++dx) {
            int sy = corner_y - dy;
            int sx = corner_x - dx;
            if (sy >= 0 && sy < DETECTOR_SIZE_Y && sx >= 0 && sx < DETECTOR_SIZE_X) {
                source[sy][sx] = 1;
            }
        }
    }
    return source;
}

// --- Функция сохранения тенеграммы в .txt и .pgm ---
void save_tengram(const std::vector<std::vector<double>>& tengram, const std::string& base_filename) {
    // Сохранение в .txt
    std::ofstream txt_file(base_filename + ".txt");
    if (!txt_file.is_open()) {
        std::cerr << "Cannot open " << base_filename << ".txt for writing.\n";
        return;
    }
    for (int i = 0; i < DETECTOR_SIZE_Y; ++i) {
        for (int j = 0; j < DETECTOR_SIZE_X; ++j) {
            txt_file << tengram[i][j];
            if (j < DETECTOR_SIZE_X - 1) txt_file << "\t";
        }
        txt_file << "\n";
    }
    txt_file.close();

    // Сохранение в .pgm
    std::ofstream pgm_file(base_filename + ".pgm");
    if (!pgm_file.is_open()) {
        std::cerr << "Cannot open " << base_filename << ".pgm for writing.\n";
        return;
    }
    pgm_file << "P2\n" << DETECTOR_SIZE_X << " " << DETECTOR_SIZE_Y << "\n255\n";
    for (int i = 0; i < DETECTOR_SIZE_Y; ++i) {
        for (int j = 0; j < DETECTOR_SIZE_X; ++j) {
            unsigned char gray_val = (tengram[i][j] > 0) ? 255 : 0;
            pgm_file << static_cast<int>(gray_val) << " ";
        }
        pgm_file << "\n";
    }
    pgm_file.close();
    std::cout << "Saved: " << base_filename << ".txt/.pgm\n";
}

// === НОВАЯ ФУНКЦИЯ: Сохранение источника (массив int) в .txt и .pgm ===
void save_source(const std::vector<std::vector<int>>& source, const std::string& base_filename) {
    // Сохранение в .txt
    std::ofstream txt_file(base_filename + ".txt");
    if (!txt_file.is_open()) {
        std::cerr << "Cannot open " << base_filename << ".txt for writing.\n";
        return;
    }
    for (int i = 0; i < DETECTOR_SIZE_Y; ++i) {
        for (int j = 0; j < DETECTOR_SIZE_X; ++j) {
            txt_file << source[i][j];
            if (j < DETECTOR_SIZE_X - 1) txt_file << "\t";
        }
        txt_file << "\n";
    }
    txt_file.close();

    // Сохранение в .pgm
    std::ofstream pgm_file(base_filename + ".pgm");
    if (!pgm_file.is_open()) {
        std::cerr << "Cannot open " << base_filename << ".pgm for writing.\n";
        return;
    }
    pgm_file << "P2\n" << DETECTOR_SIZE_X << " " << DETECTOR_SIZE_Y << "\n255\n";
    for (int i = 0; i < DETECTOR_SIZE_Y; ++i) {
        for (int j = 0; j < DETECTOR_SIZE_X; ++j) {
            unsigned char gray_val = (source[i][j] > 0) ? 255 : 0;
            pgm_file << static_cast<int>(gray_val) << " ";
        }
        pgm_file << "\n";
    }
    pgm_file.close();
    std::cout << "Saved source: " << base_filename << ".txt/.pgm\n";
}
// ======================================================================

// --- Функция поэлементного сложения двух тенеграмм ---
std::vector<std::vector<double>> sum_tengrams(
    const std::vector<std::vector<double>>& t1,
    const std::vector<std::vector<double>>& t2)
{
    std::vector<std::vector<double>> result(DETECTOR_SIZE_Y, std::vector<double>(DETECTOR_SIZE_X, 0.0));
    for (int i = 0; i < DETECTOR_SIZE_Y; ++i) {
        for (int j = 0; j < DETECTOR_SIZE_X; ++j) {
            result[i][j] = t1[i][j] + t2[i][j];
        }
    }
    return result;
}

// ============================================================================
// === ОСНОВНАЯ ПРОГРАММА (без изменений, кроме добавления сохранения источников) ===
// ============================================================================

int main()
{
    std::cout << "=== Multi-layer tengram generation ===\n";
    std::cout << "Detector size: " << DETECTOR_SIZE_Y << " x " << DETECTOR_SIZE_X << "\n";
    std::cout << "Base f = " << f_base << ", d = " << d << "\n\n";

    // --- Шаг 1: Загрузка маски ---
    std::ifstream mask_file("./mask-mura-expanded.txt");
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

    // --- Шаг 2: Создание источников ---
    // Слой 1: Круг в центре, f = f_base
    auto source_circle = create_circle_source(10.0, DETECTOR_SIZE_Y/2, DETECTOR_SIZE_X/2);
    // Слой 2: Квадрат в левом верхнем углу, f = f_base + 1
    auto source_square = create_square_source(20, 10, 10); // сторона 20, начало с (10,10)
    // Слой 3: Треугольник в правом нижнем углу, f = f_base + 2
    auto source_triangle = create_triangle_source(20, DETECTOR_SIZE_Y - 10, DETECTOR_SIZE_X - 10); // катет 20, угол в (246,246)

    std::cout << "Sources created:\n";
    std::cout << "  Layer 1: Circle (r=10) at center, f = " << f_base << "\n";
    std::cout << "  Layer 2: Square (side=20) at top-left, f = " << f_base + 1 << "\n";
    std::cout << "  Layer 3: Triangle (leg=20) at bottom-right, f = " << f_base + 2 << "\n\n";

    // === НОВОЕ: Сохранение каждого источника отдельно ===
    std::cout << "Saving individual sources...\n";
    save_source(source_circle,  "./circular_source");
    save_source(source_square,  "./square_source");
    save_source(source_triangle,"./triangle_source");
    std::cout << "\n";
    // ================================================

    // --- Шаг 3: Генерация тенеграмм для каждого слоя ---
    std::cout << "Generating tengrams...\n";

    auto tengram_f     = generate_tengram_for_source(source_circle,  mask, f_base);
    std::cout << "  Layer 1 (f=" << f_base << ") done.\n";

    auto tengram_f1    = generate_tengram_for_source(source_square,  mask, f_base + 1);
    std::cout << "  Layer 2 (f=" << f_base + 1 << ") done.\n";

    auto tengram_f2    = generate_tengram_for_source(source_triangle, mask, f_base + 2);
    std::cout << "  Layer 3 (f=" << f_base + 2 << ") done.\n\n";

    // --- Шаг 4: Сохранение индивидуальных тенеграмм ---
    std::cout << "Saving individual tengrams...\n";
    save_tengram(tengram_f,  "./tengram_layer1_f");
    save_tengram(tengram_f1, "./tengram_layer2_f+1");
    save_tengram(tengram_f2, "./tengram_layer3_f+2");
    std::cout << "\n";

    // --- Шаг 5: Сохранение комбинированных тенеграмм ---
    std::cout << "Saving combined tengrams...\n";

    // f + (f+1)
    auto tengram_f_f1 = sum_tengrams(tengram_f, tengram_f1);
    save_tengram(tengram_f_f1, "./tengram_combined_f+f+1");

    // f + (f+2)
    auto tengram_f_f2 = sum_tengrams(tengram_f, tengram_f2);
    save_tengram(tengram_f_f2, "./tengram_combined_f+f+2");

    // f + (f+1) + (f+2)
    auto tengram_all = sum_tengrams(tengram_f_f1, tengram_f2);
    save_tengram(tengram_all, "./tengram_combined_all");

    std::cout << "\n=== Program finished successfully ===\n";
    return 0;
}
