/**
 * Симуляция потока жидкости с использованием клеточного газового автомата (LGA) по модели FHP-II
 *
 * Описание:
 * Реализация двумерного клеточного газового автомата на квадратной решётке N×N.
 * Каждая клетка содержит до 8 частиц, движущихся в заданных направлениях.
 * Поддерживаются правила столкновений FHP-II, отражение от препятствий,
 * периодические граничные условия и источник частиц слева.
 *
 * Функционал:
 * - Инициализация решётки с заданной плотностью
 * - Распространение и столкновения частиц
 * - Обработка препятствий и граничных условий
 * - Визуализация состояния решётки
 * - Сохранение результата в бинарный файл grid.bin
 * - Расчёт плотности и скорости потока
 */

#include <stdio.h>
#include <stdlib.h>
#include <time.h>

/**
 * @brief Перечисление направлений движения частиц в модели FHP-II.
 *
 * Определяет 8 возможных направлений, в которых могут двигаться частицы.
 */
enum direction {
    North=0,
    NE=1,
    E=2,
    SE=3,
    S=4,
    SW=5,
    W=6,
    NW=7
};

/**
 * @brief Структура, представляющая одну ячейку решетки клеточного газового автомата.
 *
 * Каждая ячейка может содержать частицы, движущиеся в 8 направлениях (FHP-II модель),
 * а также флаг, указывающий, является ли ячейка препятствием.
 */
struct Cell {
    unsigned char particles; /**< Битовая маска: 1 бит на каждое направление. */
    unsigned char obstacle; /**< флаг: 0 - свободно, 1 - препятствие. */
};

struct Cell **board = NULL;

/**
 * @brief Инициализирует решетку случайным распределением частиц и устанавливает прямоугольное препятствие.
 *
 * @param N Размер решетки (N x N)
 * @param board Двумерный массив структур Cell, представляющий решетку.
 * @param p Вероятность наличия частицы в каждом направлении.
 * @param x Начальная координата X препятствия.
 * @param y Начальная координата Y препятствия.
 * @param w Ширина препятствия.
 * @param h Высота препятствия.
 */
void init(const int N, struct Cell **board, const double p, const int x, const int y, const int w, const int h) {
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            board[i][j].obstacle = 0;
            board[i][j].particles = 0;

            if ((j >= y) && (i >= x) && (j < y + w) && (i < x + h))
                board[i][j].obstacle = 1;
            else {
                for (int k = 0; k < 8; k++) {
                    const double particle_p = (double)rand() / RAND_MAX;
                    if (particle_p <= p)
                        board[i][j].particles |= (1 << k);
                }
            }
        }
    }
}

/**
 * @brief Обнуляет временную решетку на каждом шаге симуляции.
 *
 * @param N Размер решетки (N x N).
 * @param temp_board Временная решетка, которую нужно очистить
 * @param board Основная решетка для копирования информации о препятствиях.
 */
void reset(const int N, struct Cell **temp_board, struct Cell **board) {
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            temp_board[i][j].obstacle = board[i][j].obstacle;
            temp_board[i][j].particles = 0;
        }
    }
}

/**
 * @brief Отражает направление движения частиц на противоположное при столкновении с препятствием.
 *
 * @param dir Исходное направление движения частиц.
 * @return Противоположное направление после отражения или -1, если направление не определено.
 */
int reflection(const int dir) {
    int ref_dir = -1;

    if (dir == North)
        ref_dir = S;
    else if (dir == NE)
        ref_dir = SW;
    else if (dir == E)
        ref_dir = W;
    else if (dir == SE)
        ref_dir = NW;
    else if (dir == S)
        ref_dir = North;
    else if (dir == SW)
        ref_dir = NE;
    else if (dir == W)
        ref_dir = E;
    else if (dir == NW)
        ref_dir = SE;

    return ref_dir;
}

/**
 * @brief Перемещает частицу в указанном направлении с учетом границ, препятствий и периодических условий.
 *
 * @param N Размер решетки (N x N).
 * @param board Текущее состояние решетки.
 * @param temp_board Временная решетка для хранения промежуточного состояния.
 * @param dir Направление движения частицы.
 * @param i Координата X текущей ячейки.
 * @param j Координата Y текущей ячейки.
 */
void direction(const int N, struct Cell **board, struct Cell **temp_board, const int dir, const int i, const int j) {
    if (dir == North) {
        if (i - 1 >= 0) {
            if (board[i - 1][j].obstacle == 0)
                temp_board[i - 1][j].particles |= (1<<dir);
            else
                temp_board[i][j].particles |= (1<<reflection(dir));
        }
        else {
            const int new_i = (i + N - 1) % N;
            if (board[new_i][j].obstacle == 0)
                temp_board[new_i][j].particles |= (1<<dir);
            else
                temp_board[i][j].particles |= (1<<reflection(dir));
        }
    }

    else if (dir == NE) {
        const int new_i = (i - 1 + N) % N;
        const int new_j = (j + 1) % N;

        if (board[new_i][new_j].obstacle == 0)
            temp_board[new_i][new_j].particles |= (1<<dir);
        else
            temp_board[i][j].particles |= (1<<reflection(dir));
    }

    else if (dir == E) {
        if (j + 1 < N) {
            if (board[i][j + 1].obstacle == 0)
                temp_board[i][j + 1].particles |= (1<<dir);
            else
                temp_board[i][j].particles |= (1<<reflection(dir));
        }
        else
            temp_board[i][j].particles |= (1<<dir);
    }

    else if (dir == SE) {
        const int new_i = (i + 1) % N;
        const int new_j = (j + 1) % N;

        if (board[new_i][new_j].obstacle == 0)
            temp_board[new_i][new_j].particles |= (1<<dir);
        else
            temp_board[i][j].particles |= (1<<reflection(dir));
    }

    else if (dir == S) {
        if (i + 1 < N) {
            if (board[i + 1][j].obstacle == 0)
                temp_board[i + 1][j].particles |= (1<<dir);
            else
                temp_board[i][j].particles |= (1<<reflection(dir));
        }
        else {
            const int new_i = (i + 1 - N) % N;
            if (board[new_i][j].obstacle == 0)
                temp_board[new_i][j].particles |= (1<<dir);
            else
                temp_board[i][j].particles |= (1<<reflection(dir));
        }
    }

    else if (dir == SW) {
        const int new_i = (i + 1) % N;
        const int new_j = (j - 1 + N) % N;

        if (board[new_i][new_j].obstacle == 0)
            temp_board[new_i][new_j].particles |= (1<<dir);
        else
            temp_board[i][j].particles |= (1<<reflection(dir));
    }

    else if (dir == W) {
        if (j - 1 >= 0) {
            if (board[i][j - 1].obstacle == 0)
                temp_board[i][j - 1].particles |= (1<<dir);
            else
                temp_board[i][j].particles |= (1<<reflection(dir));
        }
        else
            temp_board[i][j].particles |= (1<<dir);
    }

    else if (dir == NW) {
        const int new_i = (i - 1 + N) % N;
        const int new_j = (j - 1 + N) % N;

        if (board[new_i][new_j].obstacle == 0)
            temp_board[new_i][new_j].particles |= (1<<dir);
        else
            temp_board[i][j].particles |= (1<<reflection(dir));
    }
}

/**
 * @brief Перемещает все частицы из текущего состояния решетки в новое положение.
 *
 * @param N Размер решетки (N x N).
 * @param board Текущее состояние решетки.
 * @param temp_board Временная решетка для хранения нового состояния после распространения.
 */
void spreading(const int N, struct Cell **board, struct Cell **temp_board) {
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            if (board[i][j].obstacle == 1)
                continue;
            for (int k = 0; k < 8; k++) {
                if ((board[i][j].particles >> k) & 1)
                    direction(N, board, temp_board, k, i, j);
            }
        }
    }
}

/**
 * @brief Обрабатывает столкновения двух частиц согласно правилам FHP-II.
 *
 * @param N Размер решетеки (N x N).
 * @param board Решетка, в которой происходит столкновение.
 * @param pos Массив из 8 элементов, обозначающих наличие частиц в каждом направлении.
 * @param i Координата X текущей ячейки.
 * @param j Координата Y текущей ячейки.
 */
void two_particles(int N, struct Cell **board, int pos[8], const int i, const int j) {
    if (pos[0] == 1 && pos[4] == 1) {
        board[i][j].particles |= (1<<2);
        board[i][j].particles |= (1<<6);
        pos[0] = pos[4] = 0;
    }
    if (pos[1] == 1 && pos[5] == 1) {
        board[i][j].particles |= (1<<0);
        board[i][j].particles |= (1<<4);
        pos[1] = pos[5] = 0;
    }
    if (pos[3] == 1 && pos[7] == 1) {
        board[i][j].particles |= (1<<4);
        board[i][j].particles |= (1<<0);
        pos[3] = pos[7] = 0;
    }
    if (pos[2] == 1 && pos[6] == 1) {
        board[i][j].particles |= (1<<1);
        board[i][j].particles |= (1<<5);
        pos[2] = pos[6] = 0;
    }
}

/**
 * @brief Обрабатывает столкновения трех частиц согласно правилам FHP-II.
 *
 * @param N Размер решетки (N x N).
 * @param board Решетка, в которой происходит столкновение.
 * @param pos Массив из 8 элементов, обозначающих наличие частиц в каждом направлении.
 * @param i Координата X текущей ячейки.
 * @param j Координата Y текущей ячейки.
 */
void three_particles(int N, struct Cell **board, int pos[8], const int i, const int j) {
    if (pos[0] == 1 && pos[3] == 1 && pos[5] == 1) {
        board[i][j].particles |= (1<<4);
        board[i][j].particles |= (1<<1);
        board[i][j].particles |= (1<<7);
        pos[0] = pos[3] = pos[5] = 0;
    }
    if (pos[0] == 1 && pos[2] == 1 && pos[3] == 1) {
        board[i][j].particles |= (1<<3);
        board[i][j].particles |= (1<<4);
        board[i][j].particles |= (1<<5);
        pos[0] = pos[2] = pos[3] = 0;
    }
    if (pos[2] == 1 && pos[4] == 1 && pos[5] == 1) {
        board[i][j].particles |= (1<<0);
        board[i][j].particles |= (1<<1);
        board[i][j].particles |= (1<<6);
        pos[2] = pos[4] = pos[5] = 0;
    }
    if (pos[0] == 1 && pos[4] == 1 && pos[7] == 1) {
        board[i][j].particles |= (1<<0);
        board[i][j].particles |= (1<<1);
        board[i][j].particles |= (1<<2);
        pos[0] = pos[4] = pos[7] = 0;
    }
}

/**
 * @brief Обрабатывает столкновения между частицами в каждой ячейке решетки.
 *
 * @param N Размер решетки (N x N).
 * @param board Решетка после распространения частиц.
 * @param temp_board Временная решетка с новым состояние после столкновения.
 */
void collision(const int N, struct Cell **board, struct Cell **temp_board) {
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            board[i][j].particles = 0;
            if (temp_board[i][j].obstacle == 1) {
                board[i][j] = temp_board[i][j];
                continue;
            }
            int pos[8] = {0, 0, 0, 0, 0, 0, 0, 0};
            for (int k = 0; k < 8; k++) {
                pos[k] = (temp_board[i][j].particles >> k) & 1;
            }

            two_particles(N, board, pos, i, j);
            three_particles(N, board, pos, i , j);

            for (int k = 0; k < 8; k++) {
                if (pos[k] != 0)
                    board[i][j].particles |= (1<<k);
            }
        }
    }
}

/**
 * @brief Рассчитывает среднюю плотность и скорость потока жидкости.
 *
 * @param N Размер решетки (N x N).
 * @param board Текущее состояние решетки.
 * @param density Указатель на переменную для сохранения средней плотности.
 * @param vx Указатель на переменную для горизонтальной компоненты скорости.
 * @param vy Указатель на переменную для вертикальной компоненты скорости.
 */
void density_flow_rate(const int N, struct Cell **board, double *density, double *vx, double *vy) {
    int total_particles = 0;
    const double jx[8] = {0.0, 0.5, 1.0, 0.5, 0.0, -0.5, -1.0, -0.5};
    const double jy[8] = {1.0, 0.5, 0.0, -0.5, -1.0, -0.5, 0.0, 0.5};
    double total_x = 0, total_y = 0;


    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            if (board[i][j].obstacle == 0) {
                for (int k = 0; k < 8; k++) {
                    if ((board[i][j].particles >> k) & 1) {
                        total_particles++;
                        total_x += jx[k];
                        total_y += jy[k];
                    }
                }
            }
        }
    }

    *density = total_particles / (double)(N * N);
    *vx = total_x / total_particles;
    *vy = total_y / total_particles;
}

/**
 * @brief Визуализирует текущее состояние решетки в консоли.
 *
 * @param N Размер решетки (N x N).
 * @param board Текущее состояние решетки.
 */
void print_board(const int N, struct Cell **board) {
    int count_particles = 0;

    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            count_particles = 0;

            if (board[i][j].obstacle == 0) {
                for (int k = 0; k < 8; k++) {
                    if ((board[i][j].particles >> k) & 1)
                        count_particles++;
                }
                if (count_particles == 0)
                    printf("\033[90m.\033[0m ");
                else if (count_particles <= 2)
                    printf("\033[32m*\033[0m ");
                else
                    printf("\033[33m#\033[0m ");
            }
            else
                printf("\033[41mX\033[0m ");
        }
        printf("\n");
    }
}

/**
 * @brief Сохраняет текущее состояние решетки в бинарный файл.
 *
 * @param N Размер решетки (N x N).
 * @param filename Имя файла для сохранения данных.
 * @param board Текущее состояние решетки.
 */
void save_bin(const int N, const char *filename, struct Cell **board){
    FILE *file = fopen(filename, "wb");
    if (!file) {
        for (int i = 0; i < N; i++) {
            free(board[i]);
        }
        free(board);
        perror("Error opening file for writing");
        exit(EXIT_FAILURE);
    }

    fwrite(&N, sizeof(int), 1, file);

    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            fwrite(&board[i][j].particles, sizeof(unsigned char), 1, file);
            fwrite(&board[i][j].obstacle, sizeof(unsigned char), 1, file);
        }
    }

    fclose(file);
}

/**
 * @brief Загружает состояние решетки из бинарного файла.
 *
 * @param filename Имя файла для чтения данных.
 * @param N_out Указатель для возврата размера решетки.
 * @return Указатель на двумерный массив структур Cell, представляющий загруженную решетку.
 */
struct Cell** load_bin(const char *filename, int *N_out) {
    FILE *file = fopen(filename, "rb");
    if (!file) {
        perror("Error opening file for reading");
        exit(EXIT_FAILURE);
    }

    fread(N_out, sizeof(int), 1, file);
    int N = *N_out;

    struct Cell **board = malloc(N * sizeof(struct Cell *));
    if (board == NULL) {
        perror("Error allocating memory");
        exit(EXIT_FAILURE);
    }
    for (int i = 0; i < N; i++) {
        board[i] = malloc(N * sizeof(struct Cell));
        if (board[i] == NULL) {
            perror("Error allocating memory");
            exit(EXIT_FAILURE);
        }
    }

    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            fread(&board[i][j].particles, sizeof(unsigned char), 1, file);
            fread(&board[i][j].obstacle, sizeof(unsigned char), 1, file);
        }
    }

    fclose(file);
    return board;
}

int main() {
    int N, x, y, w, h, iteration, ans;
    double p;
    srand(time(NULL));

    printf("Generate new (0). Load from file (1): ");
    scanf("%d", &ans);

    if (ans == 0) {
        // создание новой решетки
        printf("Enter N p x y w h iteration: ");
        scanf("%d %lf %d %d %d %d %d", &N, &p, &x, &y, &w, &h, &iteration);

        board = malloc(N * sizeof(struct Cell *));
        if (board == NULL) {
            free(board);
            perror("Error allocating memory");
            exit(EXIT_FAILURE);
        }
        for (int i = 0; i < N; i++) {
            board[i] = malloc(N * sizeof(struct Cell));
            if (board[i] == NULL) {
                for (int j = 0; j < i; j++) {
                    free(board[i]);
                }
                free(board);
                perror("Error allocating memory");
                exit(EXIT_FAILURE);
            }
        }

        init(N, board, p, x, y, w, h);
        printf("Initial state\n");
        print_board(N, board);
        printf("\n");
    }
    else {
        // загрузка существующей решетки из файла
        board = load_bin("grid.bin", &N);

        printf("Enter the number of iterations: ");
        scanf("%d", &iteration);

        printf("Loaded state\n");
        print_board(N, board);
        printf("\n");
    }

    // Создание временной решетки
    struct Cell **temp_board = malloc(N * sizeof(struct Cell *));
    if (temp_board == NULL) {
        for (int i = 0; i < N; i++) {
            free(board[i]);
        }
        free(board);
        perror("Error allocating memory");
        exit(EXIT_FAILURE);
    }
    for (int i = 0; i < N; i++) {
        temp_board[i] = malloc(N * sizeof(struct Cell));
        if (temp_board[i] == NULL) {
            for (int j = 0; j < N; j++) {
                free(board[j]);
            }
            free(board);
            perror("Error allocating memory");
            exit(EXIT_FAILURE);
        }
    }

    // Основной цикл
    for (int step = 0; step < iteration; step++) {
        reset(N, temp_board, board);
        spreading(N, board, temp_board); // распространение
        collision(N, board, temp_board); // отражение

        // источник частиц слева (движение на восток)
        for (int i = 0; i < N; i++) {
            if (board[i][0].obstacle == 0) {
                board[i][0].particles |= (1<<2);
            }
        }
    }

    // Вывод конечного состояния решетки
    printf("Final state\n");
    print_board(N, board);
    save_bin(N, "grid.bin", board);
    printf("\n");

    // Вывод плотности и среднего потока
    double density, vx, vy;
    density_flow_rate(N, board, &density, &vx, &vy);
    printf("Average density: %4f\nFlow rate: (%4f, %4f)\n", density, vx, vy);

    // Освобождение памяти
    for (int i = 0; i < N; i++) {
        free(board[i]);
        free(temp_board[i]);
    }
    free(board);
    free(temp_board);
    free(load_bin("grid.bin", &N));

    return 0;
}
