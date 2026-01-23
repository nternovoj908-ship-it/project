#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <algorithm>

// Параметры матрицы
const int MATRIX_SIZE = 256; // размер квадратной матрицы
const double WORLD_SIZE_MM = 14.08; // размер детектора в мм (например, 14.08 мм)

int main()
{
    // Открываем файл MyTree.txt
    std::ifstream file("MyTree.txt");
    if (!file.is_open())
    {
        std::cerr << "Cannot open MyTree.txt\n";
        return 1;
    }

    // Пропускаем заголовок
    std::string line;
    std::getline(file, line);

    // Матрица для хранения данных
    std::vector<std::vector<double>> matrix(MATRIX_SIZE, std::vector<double>(MATRIX_SIZE, 0.0));

    double x, y, edep, ekin, fullene;

    while (file >> x >> y >> edep >> ekin >> fullene)
    {
        // Преобразуем координаты x, y в индексы массива
        int ix = static_cast<int>((x / WORLD_SIZE_MM + 0.5) * MATRIX_SIZE);
        int iy = static_cast<int>((y / WORLD_SIZE_MM + 0.5) * MATRIX_SIZE);

        // Проверяем, что индекс в пределах матрицы
        if (ix >= 0 && ix < MATRIX_SIZE && iy >= 0 && iy < MATRIX_SIZE)
        {
            matrix[iy][ix] += edep; // суммируем энергию в одной ячейке
        }
    }

    file.close();

    // Найдём максимальное значение для нормализации
    double max_val = 0.0;
    for (int i = 0; i < MATRIX_SIZE; ++i)
        for (int j = 0; j < MATRIX_SIZE; ++j)
            if (matrix[i][j] > max_val)
                max_val = matrix[i][j];

    // Сохраняем матрицу в файл
    std::ofstream out("matrix_output.txt");
    for (int i = 0; i < MATRIX_SIZE; ++i)
    {
        for (int j = 0; j < MATRIX_SIZE; ++j)
        {
            out << matrix[i][j];
            if (j < MATRIX_SIZE - 1) out << "\t";
        }
        out << "\n";
    }
    out.close();

    // Создаём PGM-изображение (серые оттенки)
    std::ofstream img("matrix_image.pgm");
    img << "P2\n" << MATRIX_SIZE << " " << MATRIX_SIZE << "\n255\n";
    for (int i = 0; i < MATRIX_SIZE; ++i)
    {
        for (int j = 0; j < MATRIX_SIZE; ++j)
        {
            int gray = 0;
            if (max_val > 0)
                gray = static_cast<int>((matrix[i][j] / max_val) * 255);
            img << gray << " ";
        }
        img << "\n";
    }
    img.close();

    std::cout << "Matrix saved to matrix_output.txt\n";
    std::cout << "Image saved to matrix_image.pgm\n";

    return 0;
}
