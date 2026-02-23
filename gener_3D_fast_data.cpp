#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <cmath>
#include <filesystem>
#include <random>
#include <omp.h>

// Геометрические параметры из DetectorConstruction.cc
const double f_base = 17.2; // базовое расстояние от источника до коллиматора
const double d = 4.2;       // расстояние от коллиматора до детектора
const double hole_step = 0.036129; // шаг отверстий в маске (в см)
const double detector_width_mm = 14.08; // ширина детектора (в мм)

// === РАЗМЕРЫ (ИЗМЕНЕНО) ===
const int DETECTOR_SIZE_Y = 256;
const int DETECTOR_SIZE_X = 256;

const int SOURCE_SIZE_Y = 128;  // НОВОЕ: размер среза источника
const int SOURCE_SIZE_X = 128;  // НОВОЕ: размер среза источника
// ==============================================================================

// Параметры фигур
const int CIRCLE_RADIUS = 10;
const int SQUARE_SIDE = 20;
const int TRIANGLE_LEG = 20;

// Количество итераций
const int TOTAL_ITERATIONS = 6000;

// ============================================================================
// === ФУНКЦИИ СОЗДАНИЯ ФИГУР (ИЗМЕНЕНО: размер SOURCE_SIZE) =================
// ============================================================================

// Круг в центре среза
std::vector<std::vector<int>> create_circle_source(int center_y, int center_x) {
    std::vector<std::vector<int>> source(SOURCE_SIZE_Y, std::vector<int>(SOURCE_SIZE_X, 0));
    for (int sy = 0; sy < SOURCE_SIZE_Y; ++sy) {
        for (int sx = 0; sx < SOURCE_SIZE_X; ++sx) {
            double dy = sy - center_y;
            double dx = sx - center_x;
            if (dy*dy + dx*dx <= CIRCLE_RADIUS * CIRCLE_RADIUS) {
                source[sy][sx] = 1;
            }
        }
    }
    return source;
}

// Квадрат в центре среза
std::vector<std::vector<int>> create_square_source(int center_y, int center_x) {
    std::vector<std::vector<int>> source(SOURCE_SIZE_Y, std::vector<int>(SOURCE_SIZE_X, 0));
    int half = SQUARE_SIDE / 2;
    for (int sy = center_y - half; sy <= center_y + half && sy < SOURCE_SIZE_Y; ++sy) {
        for (int sx = center_x - half; sx <= center_x + half && sx < SOURCE_SIZE_X; ++sx) {
            if (sy >= 0 && sx >= 0) {
                source[sy][sx] = 1;
            }
        }
    }
    return source;
}

// Треугольник в центре среза (прямоугольный, катеты вправо и вниз)
std::vector<std::vector<int>> create_triangle_source(int center_y, int center_x) {
    std::vector<std::vector<int>> source(SOURCE_SIZE_Y, std::vector<int>(SOURCE_SIZE_X, 0));
    // Смещаем центр, чтобы треугольник был примерно центрирован
    int start_y = center_y - TRIANGLE_LEG / 2;
    int start_x = center_x - TRIANGLE_LEG / 2;
    
    for (int dy = 0; dy < TRIANGLE_LEG; ++dy) {
        for (int dx = 0; dx < TRIANGLE_LEG - dy; ++dx) {
            int sy = start_y + dy;
            int sx = start_x + dx;
            if (sy >= 0 && sy < SOURCE_SIZE_Y && sx >= 0 && sx < SOURCE_SIZE_X) {
                source[sy][sx] = 1;
            }
        }
    }
    return source;
}

// ============================================================================
// === ФУНКЦИЯ ГЕНЕРАЦИИ ТЕНЕГРАММЫ (ИЗМЕНЕНО: SOURCE_SIZE) ===================
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

                        double proj_y = (-1) * (sy_centered - (sy_centered - 0.25 * mx_centered) * (f_current + d) / f_current);
                        double proj_x = (+1) * (sx_centered - (sx_centered - 0.25 * my_centered) * (f_current + d) / f_current);

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

void save_source_txt(const std::vector<std::vector<int>>& source, const std::string& filename) {
    std::ofstream file(filename + ".txt");
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

void save_source_pgm(const std::vector<std::vector<int>>& source, const std::string& filename) {
    std::ofstream file(filename + ".pgm");
    if (!file.is_open()) return;
    file << "P2\n" << SOURCE_SIZE_X << " " << SOURCE_SIZE_Y << "\n255\n";
    for (int i = 0; i < SOURCE_SIZE_Y; ++i) {
        for (int j = 0; j < SOURCE_SIZE_X; ++j) {
            unsigned char gray_val = (source[i][j] > 0) ? 255 : 0;
            file << static_cast<int>(gray_val) << " ";
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
    // Включаем OpenMP
    omp_set_num_threads(omp_get_num_procs());
    int max_threads = omp_get_max_threads();
    
    std::cout << "=== 3D Dataset Generation (OpenMP enabled) ===\n";
    std::cout << "Available CPU threads: " << max_threads << "\n";
    std::cout << "Source slice size: " << SOURCE_SIZE_Y << "x" << SOURCE_SIZE_X << "\n";
    std::cout << "Detector size: " << DETECTOR_SIZE_Y << "x" << DETECTOR_SIZE_X << "\n";
    std::cout << "Total iterations: " << TOTAL_ITERATIONS << "\n";
    std::cout << "Base f = " << f_base << ", d = " << d << "\n\n";

    // Создаём папки
    std::filesystem::create_directories("./3D/tengram");
    std::filesystem::create_directories("./3D/true_sources");
    
    std::cout << "Output directories created.\n\n";

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

    // Генератор случайных чисел
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<> dist_y(0, SOURCE_SIZE_Y - 1);
    std::uniform_int_distribution<> dist_x(0, SOURCE_SIZE_X - 1);

    // Счётчики для прогресса
    int completed = 0;
    auto start_time = std::chrono::high_resolution_clock::now();

    std::cout << "Starting generation...\n\n";

    // ========================================================================
    // ГЛАВНЫЙ ЦИКЛ ГЕНЕРАЦИИ (3000 итераций)
    // ========================================================================
    for (int iter = 0; iter < TOTAL_ITERATIONS; ++iter) {
        
        // === Случайные координаты для каждого объекта ===
        int circle_y = dist_y(gen);
        int circle_x = dist_x(gen);
        
        int square_y = dist_y(gen);
        int square_x = dist_x(gen);
        
        int triangle_y = dist_y(gen);
        int triangle_x = dist_x(gen);
        
        // === Создание источников со случайными позициями ===
        auto source_circle = create_circle_source(circle_y, circle_x);
        auto source_square = create_square_source(square_y, square_x);
        auto source_triangle = create_triangle_source(triangle_y, triangle_x);
        
        // === Генерация тенеграмм для каждого слоя (параллельно) ===
        std::vector<std::vector<double>> tengram_f, tengram_f1, tengram_f2;
        
        #pragma omp parallel sections
        {
            #pragma omp section
            {
                tengram_f = generate_tengram_for_source(source_circle, mask, f_base);
            }
            #pragma omp section
            {
                tengram_f1 = generate_tengram_for_source(source_square, mask, f_base + 1);
            }
            #pragma omp section
            {
                tengram_f2 = generate_tengram_for_source(source_triangle, mask, f_base + 2);
            }
        }
        
        // === Суммирование тенеграмм ===
        std::vector<std::vector<double>> combined_tengram(DETECTOR_SIZE_Y, std::vector<double>(DETECTOR_SIZE_X, 0.0));
        for (int i = 0; i < DETECTOR_SIZE_Y; ++i) {
            for (int j = 0; j < DETECTOR_SIZE_X; ++j) {
                combined_tengram[i][j] = tengram_f[i][j] + tengram_f1[i][j] + tengram_f2[i][j];
            }
        }
        
        // === Сохранение тенеграммы ===
        std::string tengram_path = "./3D/tengram/fig_" + std::to_string(iter);
        save_tengram(combined_tengram, tengram_path);
        
        // === Создание папки для истинных источников ===
        std::string sources_dir = "./3D/true_sources/" + std::to_string(iter);
        std::filesystem::create_directories(sources_dir);
        
        // === Сохранение срезов ===
        // Срез 0 (f) - круг
        std::string slice0_path = sources_dir + "/fig_" + std::to_string(iter) + "_slice_0";
        save_source_txt(source_circle, slice0_path);
        
        // Срез 1 (f+1) - квадрат
        std::string slice1_path = sources_dir + "/fig_" + std::to_string(iter) + "_slice_1";
        save_source_txt(source_square, slice1_path);
        
        // Срез 2 (f+2) - треугольник
        std::string slice2_path = sources_dir + "/fig_" + std::to_string(iter) + "_slice_2";
        save_source_txt(source_triangle, slice2_path);
        
        // === Для первых 3 итераций сохраняем .pgm визуализации ===
        if (iter < 3) {
            save_source_pgm(source_circle, slice0_path);
            save_source_pgm(source_square, slice1_path);
            save_source_pgm(source_triangle, slice2_path);
            
            std::cout << "Iteration " << iter << ":\n";
            std::cout << "  Circle:    (" << circle_y << ", " << circle_x << ")\n";
            std::cout << "  Square:    (" << square_y << ", " << square_x << ")\n";
            std::cout << "  Triangle:  (" << triangle_y << ", " << triangle_x << ")\n";
            std::cout << "  Saved PGM visualizations for all 3 slices.\n\n";
        }
        
        completed++;
        
        // Прогресс каждые 100 итераций
        if (completed % 100 == 0) {
            auto current_time = std::chrono::high_resolution_clock::now();
            auto elapsed = std::chrono::duration_cast<std::chrono::seconds>(current_time - start_time).count();
            double avg_time = static_cast<double>(elapsed) / completed;
            double remaining = avg_time * (TOTAL_ITERATIONS - completed);
            
            std::cout << "Progress: " << completed << "/" << TOTAL_ITERATIONS 
                      << " (" << (100.0 * completed / TOTAL_ITERATIONS) << "%)\n";
            std::cout << "  Elapsed: " << elapsed << "s, Avg: " << avg_time << "s/iter\n";
            std::cout << "  ETA: " << static_cast<int>(remaining / 60) << "m " 
                      << static_cast<int>(remaining) % 60 << "s\n\n";
        }
    }
    
    auto end_time = std::chrono::high_resolution_clock::now();
    auto total_time = std::chrono::duration_cast<std::chrono::seconds>(end_time - start_time).count();
    
    std::cout << "\n=== Generation Completed ===\n";
    std::cout << "Total iterations: " << completed << "\n";
    std::cout << "Total time: " << total_time / 60 << "m " << total_time % 60 << "s\n";
    std::cout << "Average time per iteration: " << (static_cast<double>(total_time) / completed) << "s\n";
    std::cout << "\nDataset structure:\n";
    std::cout << "  3D/\n";
    std::cout << "  ├── tengram/           (" << TOTAL_ITERATIONS << " files: fig_0.txt ... fig_" << (TOTAL_ITERATIONS-1) << ".txt)\n";
    std::cout << "  └── true_sources/\n";
    for (int i = 0; i < std::min(3, TOTAL_ITERATIONS); ++i) {
        std::cout << "      ├── " << i << "/\n";
        std::cout << "      │   ├── fig_" << i << "_slice_0.txt  (circle, f=" << f_base << ")\n";
        std::cout << "      │   ├── fig_" << i << "_slice_1.txt  (square, f=" << (f_base+1) << ")\n";
        std::cout << "      │   └── fig_" << i << "_slice_2.txt  (triangle, f=" << (f_base+2) << ")\n";
    }
    std::cout << "      └── ...\n";
    
    return 0;
}
