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
// Порог совпадения
const double MATCH_THRESHOLD = 0.15; // 80%

// ============================================================================
// === НОВАЯ ФУНКЦИЯ: Реконструкция для заданного значения f ==================
// ============================================================================
void reconstruct_source(
    const std::vector<std::vector<double>>& real_tengram,
    const std::vector<std::vector<int>>& mask,
    double f_current,
    const std::string& output_basename)
{
    int DETECTOR_SIZE_Y = real_tengram.size();
    int DETECTOR_SIZE_X = (DETECTOR_SIZE_Y > 0) ? real_tengram[0].size() : 0;
    int MASK_SIZE_Y = mask.size();
    int MASK_SIZE_X = (MASK_SIZE_Y > 0) ? mask[0].size() : 0;

    std::cout << "  Reconstructing with f = " << f_current << "...\n";

    // Массив реконструкции
    std::vector<std::vector<int>> source_reconstruction(DETECTOR_SIZE_Y, std::vector<int>(DETECTOR_SIZE_X, 0));

    int aux_counter = 0;
    bool center_aux_saved = false;
    bool found_nonzero = false; // для вывода первого ненулевого элемента

    // Проходим по каждому элементу в массиве реконструкции
    for (int sx = 0; sx < DETECTOR_SIZE_Y; ++sx) {
        for (int sy = 0; sy < DETECTOR_SIZE_X; ++sy) {
            // Создаём вспомогательную тенеграмму (изначально всё 0)
            std::vector<std::vector<int>> aux_tengram(DETECTOR_SIZE_Y, std::vector<int>(DETECTOR_SIZE_X, 0));

            // Для первой вспомогательной тенеграммы будем выводить параметры
            bool printed_params = false;

            // Проходим по всем отверстиям маски
            for (int mx = 0; mx < MASK_SIZE_Y; ++mx) {
                for (int my = 0; my < MASK_SIZE_X; ++my) {
                    if (mask[mx][my] == 0) { // 0 = дырка, 1 = закрыто
                        // Координаты относительно центра маски (в единицах шага)
                        int mx_centered = 1.7*(mx - MASK_SIZE_Y / 2);
                        int my_centered = 1.7*(my - MASK_SIZE_X / 2);

                        // Координаты относительно центра массива реконструкции (в единицах шага)
                        int sx_centered = sx - DETECTOR_SIZE_Y / 2;
                        int sy_centered = sy - DETECTOR_SIZE_X / 2;

                        // Для первой вспомогательной тенеграммы (sx=0, sy=0) и первой дырки в маске
                        if (sx == 0 && sy == 0 && !printed_params) {
                            std::cout << "  First aux tengram params for first hole (mx=" << mx << ", my=" << my << "):\n";
                            std::cout << "    mx_centered = " << mx_centered << "\n";
                            std::cout << "    my_centered = " << my_centered << "\n";
                            std::cout << "    sx_centered = " << sx_centered << "\n";
                            std::cout << "    sy_centered = " << sy_centered << "\n";

                            // Формула проекции (исправленная, с множителями -1 и +1)
                            double dx_real = (-1)*(sx_centered - (sx_centered - 0.125*mx_centered) * (f_current + d) / f_current);
                            double dy_real = (+1)*(sy_centered - (sy_centered - 0.125*my_centered) * (f_current + d) / f_current);
                            std::cout << "    dx_real = " << dx_real << "\n";
                            std::cout << "    dy_real = " << dy_real << "\n";

                            // Преобразуем в целые координаты
                            int dx = static_cast<int>(round(dx_real + DETECTOR_SIZE_Y / 2));
                            int dy = static_cast<int>(round(dy_real + DETECTOR_SIZE_X / 2));
                            std::cout << "    dx = " << dx << "\n";
                            std::cout << "    dy = " << dy << "\n";
                            printed_params = true;
                        }

                        // === ЕДИНАЯ ФОРМУЛА ДЛЯ ВСЕХ ИСТОЧНИКОВ ===
                        double dx_real = (-1)*(sx_centered - (sx_centered - 0.125*mx_centered) * (f_current + d) / f_current);
                        double dy_real = (+1)*(sy_centered - (sy_centered - 0.125*my_centered) * (f_current + d) / f_current);

                        // Преобразуем в целые координаты
                        int dx = static_cast<int>(round(dx_real + DETECTOR_SIZE_Y / 2));
                        int dy = static_cast<int>(round(dy_real + DETECTOR_SIZE_X / 2));

                        // Проверяем, входят ли координаты в размер детектора
                        if (dx >= 0 && dx < DETECTOR_SIZE_Y && dy >= 0 && dy < DETECTOR_SIZE_X) {
                            aux_tengram[dx][dy]++;
                        }
                    }
                }
            }

            // Сохраняем вспомогательную тенеграмму, если это центр (только для первого среза)
            if (sx == DETECTOR_SIZE_X / 2 && sy == DETECTOR_SIZE_Y / 2 && !center_aux_saved && output_basename.find("f") != std::string::npos && output_basename.find("+") == std::string::npos) {
                std::string filename_txt = "C:/Users/User/c++1/aux_tengram_center_f.txt";
                std::ofstream aux_out_txt(filename_txt);
                for (int i = 0; i < DETECTOR_SIZE_Y; ++i) {
                    for (int j = 0; j < DETECTOR_SIZE_X; ++j) {
                        aux_out_txt << aux_tengram[i][j];
                        if (j < DETECTOR_SIZE_X - 1) aux_out_txt << "\t";
                    }
                    aux_out_txt << "\n";
                }
                aux_out_txt.close();

                std::string filename_pgm = "C:/Users/User/c++1/aux_tengram_center_f.pgm";
                std::ofstream aux_out_pgm(filename_pgm);
                aux_out_pgm << "P2\n" << DETECTOR_SIZE_X << " " << DETECTOR_SIZE_Y << "\n255\n";
                for (int i = 0; i < DETECTOR_SIZE_Y; ++i) {
                    for (int j = 0; j < DETECTOR_SIZE_X; ++j) {
                        int gray = aux_tengram[i][j] ? 255 : 0; // белый = проекция
                        aux_out_pgm << gray << " ";
                    }
                    aux_out_pgm << "\n";
                }
                aux_out_pgm.close();
                center_aux_saved = true;
                std::cout << "  Center aux tengram saved as aux_tengram_center_f.txt/.pgm\n";
            }

            // Выводим первые 3 вспомогательные тенеграммы (только для первого среза)
            if (aux_counter < 3 && output_basename.find("f") != std::string::npos && output_basename.find("+") == std::string::npos) {
                std::string filename_txt = "C:/Users/User/c++1/aux_tengram_" + std::to_string(aux_counter + 1) + "_f.txt";
                std::ofstream aux_out_txt(filename_txt);
                for (int i = 0; i < DETECTOR_SIZE_Y; ++i) {
                    for (int j = 0; j < DETECTOR_SIZE_X; ++j) {
                        aux_out_txt << aux_tengram[i][j];
                        if (j < DETECTOR_SIZE_X - 1) aux_out_txt << "\t";
                    }
                    aux_out_txt << "\n";
                }
                aux_out_txt.close();

                std::string filename_pgm = "C:/Users/User/c++1/aux_tengram_" + std::to_string(aux_counter + 1) + "_f.pgm";
                std::ofstream aux_out_pgm(filename_pgm);
                aux_out_pgm << "P2\n" << DETECTOR_SIZE_X << " " << DETECTOR_SIZE_Y << "\n255\n";
                for (int i = 0; i < DETECTOR_SIZE_Y; ++i) {
                    for (int j = 0; j < DETECTOR_SIZE_X; ++j) {
                        int gray = aux_tengram[i][j] ? 255 : 0; // белый = проекция
                        aux_out_pgm << gray << " ";
                    }
                    aux_out_pgm << "\n";
                }
                aux_out_pgm.close();

                // Также выводим информацию о первой вспомогательной тенеграмме
                if (aux_counter == 0) {
                    std::cout << "  Aux tengram 1 size: " << DETECTOR_SIZE_Y << " x " << DETECTOR_SIZE_X << "\n";
                    int aux_nonzero = 0;
                    for (int i = 0; i < DETECTOR_SIZE_Y; ++i) {
                        for (int j = 0; j < DETECTOR_SIZE_X; ++j) {
                            if (aux_tengram[i][j] != 0) aux_nonzero++;
                        }
                    }
                    std::cout << "  Non-zero elements in aux tengram 1: " << aux_nonzero << "\n";

                    found_nonzero = false;
                    for (int i = 0; i < DETECTOR_SIZE_Y; ++i) {
                        for (int j = 0; j < DETECTOR_SIZE_X; ++j) {
                            if (aux_tengram[i][j] != 0) {
                                std::cout << "  First non-zero element in aux tengram 1 at: (" << i << ", " << j << ") = " << aux_tengram[i][j] << "\n";
                                found_nonzero = true;
                                break;
                            }
                        }
                        if (found_nonzero) break;
                    }
                }
                aux_counter++;
            }

            // Считаем количество совпадений
            int matches = 0;
            int total_aux = 0;
            for (int i = 0; i < DETECTOR_SIZE_Y; ++i) {
                for (int j = 0; j < DETECTOR_SIZE_X; ++j) {
                    if (aux_tengram[i][j] > 0) {
                        total_aux++;
                        if (real_tengram[i][j] > 0) {
                            matches++;
                        }
                    }
                }
            }

            // Проверяем, удовлетворяет ли порогу
            if (total_aux > 0) {
                double ratio = static_cast<double>(matches) / total_aux;
                // Новое правило для обозначения источника:
                if (ratio >= 0.99) {
                    source_reconstruction[sx][sy]++;
                }
                else if (ratio >= 0.89 && matches >= 5) {
                    source_reconstruction[sx][sy]++;
                }
                else if (ratio >= 0.89 && matches >= 7) {
                    source_reconstruction[sx][sy]++;
                }
            }
        }
    }

    // Сохраняем результат в .txt
    std::ofstream out(output_basename + ".txt");
    for (int i = 0; i < DETECTOR_SIZE_Y; ++i) {
        for (int j = 0; j < DETECTOR_SIZE_X; ++j) {
            out << source_reconstruction[i][j];
            if (j < DETECTOR_SIZE_X - 1) out << "\t";
        }
        out << "\n";
    }
    out.close();

    // Создаём PGM-изображение
    std::ofstream img(output_basename + ".pgm");
    img << "P2\n" << DETECTOR_SIZE_X << " " << DETECTOR_SIZE_Y << "\n255\n";
    for (int i = 0; i < DETECTOR_SIZE_Y; ++i) {
        for (int j = 0; j < DETECTOR_SIZE_X; ++j) {
            int gray = source_reconstruction[i][j] ? 255 : 0; // белый = источник
            img << gray << " ";
        }
        img << "\n";
    }
    img.close();

    std::cout << "  Reconstructed source saved to " << output_basename << ".txt/.pgm\n";
}
// ============================================================================

int main()
{
    std::cout << "=== Multi-slice reconstruction ===\n";
    std::cout << "Base f = " << f_base << ", d = " << d << "\n\n";

    // --- Читаем тенеграмму (одна для всех срезов) ---
    std::ifstream tengram_file("C:/Users/User/c++1/generated_circular_tengram.txt");
    if (!tengram_file.is_open()) {
        std::cerr << "Cannot open generated_circular_tengram.txt\n";
        return 1;
    }
    std::vector<std::vector<double>> real_tengram;
    std::string line;
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

    int DETECTOR_SIZE_Y = real_tengram.size();
    int DETECTOR_SIZE_X = (DETECTOR_SIZE_Y > 0) ? real_tengram[0].size() : 0;

    if (DETECTOR_SIZE_X == 0 || DETECTOR_SIZE_Y == 0) {
        std::cerr << "Empty tengram file\n";
        return 1;
    }
    std::cout << "Real tengram size: " << DETECTOR_SIZE_Y << " x " << DETECTOR_SIZE_X << "\n";

    // Подсчитаем ненулевые элементы
    int nonzero_count = 0;
    for (int i = 0; i < DETECTOR_SIZE_Y; ++i) {
        for (int j = 0; j < DETECTOR_SIZE_X; ++j) {
            if (real_tengram[i][j] != 0) {
                nonzero_count++;
            }
        }
    }
    std::cout << "Non-zero elements in real tengram: " << nonzero_count << "\n\n";

    // --- Читаем маску (одна для всех срезов) ---
    std::ifstream mask_file("C:/Users/User/c++1/mask-mura-expanded.txt");
    if (!mask_file.is_open()) {
        std::cerr << "Cannot open mask-mura-expanded.txt\n";
        return 1;
    }
    std::vector<std::vector<int>> mask;
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
    std::cout << "Mask size: " << MASK_SIZE_Y << " x " << MASK_SIZE_X << "\n\n";

    // --- Выполняем реконструкцию для трёх значений f ---
    std::cout << "Starting reconstruction passes...\n\n";

    // Срез 1: f = f_base
    reconstruct_source(real_tengram, mask, f_base, "C:/Users/User/c++1/reconstructed_source_f");
    std::cout << "\n";

    // Срез 2: f = f_base + 1
    reconstruct_source(real_tengram, mask, f_base + 1, "C:/Users/User/c++1/reconstructed_source_f+1");
    std::cout << "\n";

    // Срез 3: f = f_base + 2
    reconstruct_source(real_tengram, mask, f_base + 2, "C:/Users/User/c++1/reconstructed_source_f+2");
    std::cout << "\n";

    std::cout << "=== All reconstructions completed ===\n";
    return 0;
}