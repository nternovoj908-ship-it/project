#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <cmath>
#include <filesystem>

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
const int MARGIN = 20;                  // Отступ от края плоскости источника
const int GRID_N = 14;                  // Количество позиций на ось для каждой фигуры
// Размеры фигур (в индексах источника)
const int CIRCLE_RADIUS = 10;
const int SQUARE_SIDE = 18;
const int TRIANGLE_LEG = 20;
const int RECT_WIDTH = 24, RECT_HEIGHT = 12;
const int ELLIPSE_RX = 14, ELLIPSE_RY = 8;
// ==============================================================================

// Типы геометрических фигур
enum class ShapeType { CIRCLE, SQUARE, TRIANGLE, RECTANGLE, ELLIPSE, COUNT };

// ============================================================================
// === ФУНКЦИИ СОЗДАНИЯ ФИГУР =================================================
// ============================================================================

// Круг
void create_circle(std::vector<std::vector<int>>& source, int cy, int cx, int radius) {
    for (int sy = 0; sy < SOURCE_SIZE_Y; ++sy) {
        for (int sx = 0; sx < SOURCE_SIZE_X; ++sx) {
            double dy = sy - cy;
            double dx = sx - cx;
            if (dy*dy + dx*dx <= radius * radius) {
                source[sy][sx] = 1;
            }
        }
    }
}

// Квадрат (сторона side, центр в cy,cx)
void create_square(std::vector<std::vector<int>>& source, int cy, int cx, int side) {
    int half = side / 2;
    for (int sy = cy - half; sy <= cy + half && sy < SOURCE_SIZE_Y; ++sy) {
        for (int sx = cx - half; sx <= cx + half && sx < SOURCE_SIZE_X; ++sx) {
            if (sy >= 0 && sx >= 0) {
                source[sy][sx] = 1;
            }
        }
    }
}

// Прямоугольный треугольник (прямой угол в центре, катеты направлены вправо и вниз)
void create_triangle(std::vector<std::vector<int>>& source, int cy, int cx, int leg) {
    for (int dy = 0; dy < leg; ++dy) {
        for (int dx = 0; dx < leg - dy; ++dx) {
            int sy = cy + dy;
            int sx = cx + dx;
            if (sy >= 0 && sy < SOURCE_SIZE_Y && sx >= 0 && sx < SOURCE_SIZE_X) {
                source[sy][sx] = 1;
            }
        }
    }
}

// Прямоугольник (ширина w, высота h, центр в cy,cx)
void create_rectangle(std::vector<std::vector<int>>& source, int cy, int cx, int w, int h) {
    int hw = w / 2, hh = h / 2;
    for (int sy = cy - hh; sy <= cy + hh && sy < SOURCE_SIZE_Y; ++sy) {
        for (int sx = cx - hw; sx <= cx + hw && sx < SOURCE_SIZE_X; ++sx) {
            if (sy >= 0 && sx >= 0) {
                source[sy][sx] = 1;
            }
        }
    }
}

// Эллипс (радиусы rx, ry, центр в cy,cx)
void create_ellipse(std::vector<std::vector<int>>& source, int cy, int cx, int rx, int ry) {
    for (int sy = 0; sy < SOURCE_SIZE_Y; ++sy) {
        for (int sx = 0; sx < SOURCE_SIZE_X; ++sx) {
            double dy = (sy - cy) / static_cast<double>(ry);
            double dx = (sx - cx) / static_cast<double>(rx);
            if (dy*dy + dx*dx <= 1.0) {
                source[sy][sx] = 1;
            }
        }
    }
}

// Фабрика фигур
void create_shape(std::vector<std::vector<int>>& source, ShapeType type, int cy, int cx) {
    // Очищаем источник
    for (auto& row : source) std::fill(row.begin(), row.end(), 0);
    
    switch (type) {
        case ShapeType::CIRCLE:
            create_circle(source, cy, cx, CIRCLE_RADIUS);
            break;
        case ShapeType::SQUARE:
            create_square(source, cy, cx, SQUARE_SIDE);
            break;
        case ShapeType::TRIANGLE:
            create_triangle(source, cy, cx, TRIANGLE_LEG);
            break;
        case ShapeType::RECTANGLE:
            create_rectangle(source, cy, cx, RECT_WIDTH, RECT_HEIGHT);
            break;
        case ShapeType::ELLIPSE:
            create_ellipse(source, cy, cx, ELLIPSE_RX, ELLIPSE_RY);
            break;
        default:
            create_circle(source, cy, cx, CIRCLE_RADIUS);
    }
}

// ============================================================================
// === ФУНКЦИИ ГЕНЕРАЦИИ И СОХРАНЕНИЯ (без изменений логики) ==================
// ============================================================================

std::vector<std::vector<double>> generate_tengram_at_position(
    const std::vector<std::vector<int>>& source,
    const std::vector<std::vector<int>>& mask,
    int MASK_SIZE_Y, int MASK_SIZE_X)
{
    std::vector<std::vector<double>> tengram(DETECTOR_SIZE_Y, std::vector<double>(DETECTOR_SIZE_X, 0.0));
    
    for (int sy = 0; sy < SOURCE_SIZE_Y; ++sy) {
        for (int sx = 0; sx < SOURCE_SIZE_X; ++sx) {
            if (source[sy][sx] == 0) continue;
            
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
    return tengram;
}

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

// ============================================================================
// === ОСНОВНАЯ ПРОГРАММА =====================================================
// ============================================================================

int main()
{
    // Создаём папки для вывода
    std::filesystem::create_directories("./output/sources");
    std::filesystem::create_directories("./output/tenegrams");
    
    // Расчёт параметров сетки
    const int usable_size = SOURCE_SIZE_Y - 2 * MARGIN;  // 128 - 40 = 88
    const double grid_step = static_cast<double>(usable_size) / GRID_N;  // 88 / 14 ≈ 6.29
    
    std::cout << "=== Generating diverse dataset (GRID MODE) ===\n";
    std::cout << "Source plane: " << SOURCE_SIZE_Y << "x" << SOURCE_SIZE_X << "\n";
    std::cout << "Detector size: " << DETECTOR_SIZE_Y << "x" << DETECTOR_SIZE_X << "\n";
    std::cout << "Margin: " << MARGIN << ", Usable area: " << usable_size << "x" << usable_size << "\n";
    std::cout << "Grid: " << GRID_N << "x" << GRID_N << " = " << (GRID_N*GRID_N) << " positions per shape\n";
    std::cout << "Grid step: " << grid_step << "\n";
    std::cout << "Shapes: circle, square, triangle, rectangle, ellipse (5 total)\n";
    std::cout << "Expected pairs: " << (GRID_N*GRID_N*5) << "\n\n";

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

    // Предварительный расчёт позиций сетки (одинаково для всех фигур)
    std::vector<int> grid_positions;
    for (int g = 0; g < GRID_N; ++g) {
        int pos = MARGIN + static_cast<int>(round(g * grid_step));
        // Ограничиваем в допустимых границах
        pos = std::max(MARGIN, std::min(SOURCE_SIZE_Y - MARGIN - 1, pos));
        grid_positions.push_back(pos);
    }
    
    std::cout << "Grid positions (first 5): ";
    for (int i = 0; i < std::min(5, GRID_N); ++i) {
        std::cout << grid_positions[i] << " ";
    }
    std::cout << "...\n\n";

    // Генерация пар источник-тенеграмма
    int pair_index = 0;
    
    // Циклически перебираем фигуры
    for (int shape_idx = 0; shape_idx < static_cast<int>(ShapeType::COUNT); ++shape_idx) {
        ShapeType current_shape = static_cast<ShapeType>(shape_idx);
        std::cout << "Processing shape " << (shape_idx + 1) << "/" << static_cast<int>(ShapeType::COUNT) << "... ";
        
        // Проходим по детерминированной сетке позиций
        for (int gy = 0; gy < GRID_N; ++gy) {
            for (int gx = 0; gx < GRID_N; ++gx) {
                
                // Позиция центра фигуры на детерминированной сетке
                int cy = grid_positions[gy];
                int cx = grid_positions[gx];
                
                // 1. Создаём источник с текущей фигурой
                std::vector<std::vector<int>> true_source(SOURCE_SIZE_Y, std::vector<int>(SOURCE_SIZE_X, 0));
                create_shape(true_source, current_shape, cy, cx);
                
                // 2. Генерируем тенеграмму
                auto tengram = generate_tengram_at_position(true_source, mask, MASK_SIZE_Y, MASK_SIZE_X);
                
                // 3. Сохраняем пару файлов
                std::string source_file = "./output/sources/circus_2D_" + std::to_string(pair_index) + ".txt";
                std::string tengram_file = "./output/tenegrams/" + std::to_string(pair_index) + ".txt";
                
                save_source(true_source, source_file);
                save_tengram(tengram, tengram_file);
                
                pair_index++;
            }
        }
        std::cout << "generated " << (GRID_N*GRID_N) << " pairs\n";
    }
    
    std::cout << "\n=== Dataset generation completed ===\n";
    std::cout << "Total pairs generated: " << pair_index << "\n";
    std::cout << "Sources saved to: ./output/sources/circus_2D_*.txt\n";
    std::cout << "Teneagrams saved to: ./output/tenegrams/*.txt\n";
    std::cout << "\nShape distribution (exact):\n";
    int per_shape = GRID_N * GRID_N;
    std::cout << "  Circle:     " << per_shape << " (indices 0-" << per_shape-1 << ")\n";
    std::cout << "  Square:     " << per_shape << " (indices " << per_shape << "-" << 2*per_shape-1 << ")\n";
    std::cout << "  Triangle:   " << per_shape << " (indices " << 2*per_shape << "-" << 3*per_shape-1 << ")\n";
    std::cout << "  Rectangle:  " << per_shape << " (indices " << 3*per_shape << "-" << 4*per_shape-1 << ")\n";
    std::cout << "  Ellipse:    " << per_shape << " (indices " << 4*per_shape << "-" << 5*per_shape-1 << ")\n";
    
    return 0;
}
    

