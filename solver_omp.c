#define _GNU_SOURCE
#define _POSIX_C_SOURCE 200809L
#include <stdio.h>
#include <errno.h>
#include <float.h>
#include <limits.h>
#include <stdint.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>
#ifdef __linux__
#include <sched.h>
#endif

typedef struct {
    double x;
    double y;
} Point;

/* Вариант 5: A(0,0), B(3,0), C(3,2), D(1,3) */
static const Point DOMAIN_POLYGON[4] = {
    {0.0, 0.0},
    {3.0, 0.0},
    {3.0, 2.0},
    {1.0, 3.0}
};

#define DOMAIN_VERTICES 4
#define MAX_POLYGON_VERTICES 16

static const double X_MIN = -1.0;
static const double X_MAX = 4.0;
static const double Y_MIN = -1.0;
static const double Y_MAX = 4.0;


static int parse_positive_int(const char *text)
{
    char *end;
    errno = 0;
    long value = strtol(text, &end, 10);
    if (errno || end == text || *end || value < 1 || value > INT_MAX) {
        fprintf(stderr, "Invalid positive integer: %s\n", text);
        exit(EXIT_FAILURE);
    }
    return (int)value;
}

static double parse_tolerance(const char *text)
{
    char *end;
    errno = 0;
    double value = strtod(text, &end);
    if (errno || end == text || *end || !isfinite(value) || value <= 0.0) {
        fprintf(stderr, "Tolerance must be finite and positive\n");
        exit(EXIT_FAILURE);
    }
    return value;
}

static void *xcalloc(size_t n, size_t size)
{
    if (size && n > SIZE_MAX / size) {
        fprintf(stderr, "Allocation size overflow\n");
        exit(EXIT_FAILURE);
    }
    void *ptr = calloc(n, size);
    if (ptr == NULL) {
        fprintf(stderr, "Memory allocation error\n");
        exit(EXIT_FAILURE);
    }
    return ptr;
}

/* Индекс внутреннего узла в одномерном массиве */
static inline size_t idx(int i, int j, int N)
{
    return (size_t)(i - 1) * (size_t)(N - 1) + (size_t)(j - 1);
}

static inline double f_function(double x, double y)
{
    (void)x;
    (void)y;
    return 1.0;
}

/* Продолжение граничной функции на всю Pi */
static inline double v0_function(double x, double y)
{
    return sin(x + y);
}

/* Отсечение многоугольника полуплоскостью a*x+b*y<=c */
static int clip_halfplane(const Point *input, int n, Point *output,
                          double a, double b, double c)
{
    if (n <= 0)
        return 0;
    int out_n = 0;

    for (int k = 0; k < n; ++k) {
        Point S = input[k];
        Point E = input[(k + 1) % n];
        double fs = a * S.x + b * S.y - c;
        double fe = a * E.x + b * E.y - c;
        int Sin = (fs <= 1.0e-13);
        int Ein = (fe <= 1.0e-13);

        if (Sin && Ein) {
            output[out_n++] = E;
        } else if (Sin && !Ein) {
            double dx = E.x - S.x;
            double dy = E.y - S.y;
            double denominator = a * dx + b * dy;
            if (fabs(denominator) > 1.0e-15) {
                double t = (c - a * S.x - b * S.y) / denominator;
                output[out_n].x = S.x + t * dx;
                output[out_n].y = S.y + t * dy;
                ++out_n;
            }
        } else if (!Sin && Ein) {
            double dx = E.x - S.x;
            double dy = E.y - S.y;
            double denominator = a * dx + b * dy;
            if (fabs(denominator) > 1.0e-15) {
                double t = (c - a * S.x - b * S.y) / denominator;
                output[out_n].x = S.x + t * dx;
                output[out_n].y = S.y + t * dy;
                ++out_n;
            }
            output[out_n++] = E;
        }
    }
    return out_n;
}

/* Площадь многоугольника */
static double polygon_area(const Point *polygon, int n)
{
    if (n < 3)
        return 0.0;
    double area = 0.0;
    /* Translate to a polygon vertex to avoid cancellation of large x*y
       products for tiny cells (important because Q is scaled by h^-4). */
    for (int i = 1; i + 1 < n; ++i) {
        double ax = polygon[i].x - polygon[0].x;
        double ay = polygon[i].y - polygon[0].y;
        double bx = polygon[i + 1].x - polygon[0].x;
        double by = polygon[i + 1].y - polygon[0].y;
        area += ax * by - bx * ay;
    }
    return 0.5 * fabs(area);
}

/* Доля площади ячейки S_ij/(h1*h2), лежащая внутри D */
static double cell_fraction_inside_domain(double xmin, double xmax,
                                          double ymin, double ymax)
{
    /* Локальные массивы: у каждой OpenMP-нити собственные p1 и p2 */
    Point p1[MAX_POLYGON_VERTICES];
    Point p2[MAX_POLYGON_VERTICES];
    for (int k = 0; k < DOMAIN_VERTICES; ++k)
        p1[k] = DOMAIN_POLYGON[k];

    int n = DOMAIN_VERTICES;
    n = clip_halfplane(p1, n, p2, -1.0, 0.0, -xmin);
    n = clip_halfplane(p2, n, p1, 1.0, 0.0, xmax);
    n = clip_halfplane(p1, n, p2, 0.0, -1.0, -ymin);
    n = clip_halfplane(p2, n, p1, 0.0, 1.0, ymax);

    double cell_area = (xmax - xmin) * (ymax - ymin);
    if (cell_area <= 0.0)
        return 0.0;

    double fraction = polygon_area(p1, n) / cell_area;
    if (fraction < 0.0)
        fraction = 0.0;
    if (fraction > 1.0)
        fraction = 1.0;
    return fraction;
}

static int point_inside_domain(double x, double y)
{
    Point P = {x, y};
    for (int k = 0; k < DOMAIN_VERTICES; ++k) {
        Point A = DOMAIN_POLYGON[k];
        Point B = DOMAIN_POLYGON[(k + 1) % DOMAIN_VERTICES];
        double cross = (B.x - A.x) * (P.y - A.y) -
                       (B.y - A.y) * (P.x - A.x);
        if (cross < -1.0e-12)
            return 0;
    }
    return 1;
}

/* Формирование Q_ij, F_ij, диагонали A и правой части */
static double build_system(int M, int N, double h1, double h2,
                           double epsilon, double *q,
                           double *diag, double *b)
{
    const double inv_h1_2 = 1.0 / (h1 * h1);
    const double inv_h2_2 = 1.0 / (h2 * h2);
    const double inv_epsilon_2 = 1.0 / (epsilon * epsilon);
    double reconstructed_area = 0.0;

#pragma omp parallel for reduction(+:reconstructed_area) schedule(static)
    for (int i = 1; i <= M - 1; ++i) {
        for (int j = 1; j <= N - 1; ++j) {
            size_t k = idx(i, j, N);
            double x = X_MIN + i * h1;
            double y = Y_MIN + j * h2;
            double xmin = x - 0.5 * h1;
            double xmax = x + 0.5 * h1;
            double ymin = y - 0.5 * h2;
            double ymax = y + 0.5 * h2;

            double fraction = cell_fraction_inside_domain(xmin, xmax, ymin, ymax);
            reconstructed_area += fraction * h1 * h2;

            q[k] = (1.0 - fraction) * inv_epsilon_2;
            double F = fraction * f_function(x, y);
            b[k] = F + q[k] * v0_function(x, y);

            /* Перенос условий Дирихле с Gamma в правую часть */
            if (i == 1)
                b[k] += inv_h1_2 * v0_function(X_MIN, y);
            if (i == M - 1)
                b[k] += inv_h1_2 * v0_function(X_MAX, y);
            if (j == 1)
                b[k] += inv_h2_2 * v0_function(x, Y_MIN);
            if (j == N - 1)
                b[k] += inv_h2_2 * v0_function(x, Y_MAX);

            diag[k] = 2.0 * inv_h1_2 + 2.0 * inv_h2_2 + q[k];
        }
    }
    return reconstructed_area;
}

/* Параллельное матрично-свободное вычисление Aw */
static void apply_A(const double *w, double *Aw, const double *diag,
                    int M, int N, double h1, double h2)
{
    const double inv_h1_2 = 1.0 / (h1 * h1);
    const double inv_h2_2 = 1.0 / (h2 * h2);

#pragma omp parallel for schedule(static)
    for (int i = 1; i <= M - 1; ++i) {
        for (int j = 1; j <= N - 1; ++j) {
            size_t k = idx(i, j, N);
            double value = diag[k] * w[k];
            if (i > 1)
                value -= inv_h1_2 * w[idx(i - 1, j, N)];
            if (i < M - 1)
                value -= inv_h1_2 * w[idx(i + 1, j, N)];
            if (j > 1)
                value -= inv_h2_2 * w[idx(i, j - 1, N)];
            if (j < N - 1)
                value -= inv_h2_2 * w[idx(i, j + 1, N)];
            Aw[k] = value;
        }
    }
}

static void apply_A_team(const double *w, double *Aw, const double *diag,
                    int M, int N, double h1, double h2)
{
    const double inv_h1_2 = 1.0 / (h1 * h1);
    const double inv_h2_2 = 1.0 / (h2 * h2);

#pragma omp for schedule(static)
    for (int i = 1; i <= M - 1; ++i) {
        for (int j = 1; j <= N - 1; ++j) {
            size_t k = idx(i, j, N);
            double value = diag[k] * w[k];
            if (i > 1)
                value -= inv_h1_2 * w[idx(i - 1, j, N)];
            if (i < M - 1)
                value -= inv_h1_2 * w[idx(i + 1, j, N)];
            if (j > 1)
                value -= inv_h2_2 * w[idx(i, j - 1, N)];
            if (j < N - 1)
                value -= inv_h2_2 * w[idx(i, j + 1, N)];
            Aw[k] = value;
        }
    }
}

typedef struct {
    size_t iterations;
    size_t restarts;
    int converged;
    double step_norm;
    double residual_norm;
} SolverResult;

/* Метод сопряжённых градиентов с диагональным предобуславливанием */
static SolverResult conjugate_gradient(int M, int N, double h1, double h2,
                                       const double *diag, const double *b,
                                       double *w, double tolerance,
                                       double *r, double *z, double *p, double *Ap,
                                       double *previous_w)
{
    const size_t n = (size_t)(M - 1) * (size_t)(N - 1);
    const double gamma_n = ((double)n * DBL_EPSILON) /
                           (1.0 - (double)n * DBL_EPSILON);
    double rz = 0.0, pAp = 0.0, pp = 0.0, alpha = 0.0;
    double rz_new = 0.0, beta = 0.0, rr = 0.0;
    double previous_J = 0.0, previous_J_abs = 0.0;
    double candidate_J = 0.0, candidate_J_abs = 0.0;
    double step_norm = INFINITY;
    int stop = 0, restart = 0;
    SolverResult result = {0, 0, 0, INFINITY, INFINITY};

    /* One team for the entire solve. All threads follow the same branches.
       Implicit barriers on for/single publish scalars and completed vectors. */
#pragma omp parallel
    {
#pragma omp for reduction(+:rz) schedule(static)
        for (size_t k = 0; k < n; ++k) {
            w[k] = 0.0;
            r[k] = b[k];
            z[k] = r[k] / diag[k];
            p[k] = z[k];
            rz += r[k] * z[k];
        }
#pragma omp single
        {
            if (rz == 0.0) {
                result.converged = 1;
                result.step_norm = 0.0;
                stop = 1;
            }
        }
        for (size_t iteration = 0; iteration < n; ++iteration) {
            if (stop) break;
#pragma omp single
            {
                pAp = pp = candidate_J = candidate_J_abs = rz_new = 0.0;
                restart = 0;
            }
            apply_A_team(p, Ap, diag, M, N, h1, h2);
#pragma omp for reduction(+:pAp,pp) schedule(static)
            for (size_t k = 0; k < n; ++k) {
                pAp += p[k] * Ap[k];
                pp += p[k] * p[k];
            }
#pragma omp single
            {
                if (!isfinite(pAp) || pAp <= 0.0 || !isfinite(rz) || rz <= 0.0)
                    stop = 1;
                else {
                    alpha = rz / pAp;
                    step_norm = fabs(alpha) * sqrt(h1 * h2 * pp);
                }
            }
            if (stop) break;
#pragma omp for reduction(+:candidate_J,candidate_J_abs) schedule(static)
            for (size_t k = 0; k < n; ++k) {
                previous_w[k] = w[k];
                w[k] += alpha * p[k];
                r[k] -= alpha * Ap[k];
                double term = (b[k] + r[k]) * w[k];
                candidate_J += term;
                candidate_J_abs += fabs(term);
            }
#pragma omp single
            {
                candidate_J *= -0.5 * h1 * h2;
                candidate_J_abs *= 0.5 * h1 * h2;
                double slack = 8.0 * gamma_n * fmax(1.0, previous_J_abs + candidate_J_abs);
                if (!isfinite(candidate_J) || !isfinite(step_norm)) stop = 1;
                else if (candidate_J > previous_J + slack) {
                    restart = 1;
                    rz = previous_J = previous_J_abs = 0.0;
                } else {
                    previous_J = candidate_J;
                    previous_J_abs = candidate_J_abs;
                    result.iterations = iteration + 1;
                    result.step_norm = step_norm;
                    if (step_norm < tolerance) {
                        result.converged = 1;
                        stop = 1;
                    }
                }
            }
            if (stop) break;
            if (restart) {
#pragma omp for schedule(static)
                for (size_t k = 0; k < n; ++k) w[k] = previous_w[k];
                apply_A_team(w, Ap, diag, M, N, h1, h2);
#pragma omp for reduction(+:rz,previous_J,previous_J_abs) schedule(static)
                for (size_t k = 0; k < n; ++k) {
                    r[k] = b[k] - Ap[k];
                    z[k] = r[k] / diag[k];
                    p[k] = z[k];
                    rz += r[k] * z[k];
                    double term = (b[k] + r[k]) * w[k];
                    previous_J += term;
                    previous_J_abs += fabs(term);
                }
#pragma omp single
                {
                    previous_J *= -0.5 * h1 * h2;
                    previous_J_abs *= 0.5 * h1 * h2;
                    ++result.restarts;
                    if (result.restarts >= 32) stop = 1;
                    else {
                        result.iterations = iteration + 1;
                        result.step_norm = INFINITY;
                    }
                }
                continue;
            }
#pragma omp for reduction(+:rz_new) schedule(static)
            for (size_t k = 0; k < n; ++k) {
                z[k] = r[k] / diag[k];
                rz_new += r[k] * z[k];
            }
#pragma omp single
            {
                if (!isfinite(rz_new) || rz_new < 0.0) stop = 1;
                else { beta = rz_new / rz; rz = rz_new; }
            }
            if (stop) break;
#pragma omp for schedule(static)
            for (size_t k = 0; k < n; ++k) p[k] = z[k] + beta * p[k];
        }
#pragma omp for reduction(+:rr) schedule(static)
        for (size_t k = 0; k < n; ++k) rr += r[k] * r[k];
    }
    result.residual_norm = sqrt(h1 * h2 * rr);
    return result;
}

/* Запись решения последовательная и не входит в Algorithm time */
static int write_solution(const char *filename, const double *w,
                           int M, int N, double h1, double h2)
{
    FILE *fp = fopen(filename, "w");
    if (fp == NULL) {
        fprintf(stderr, "Cannot open %s\n", filename);
        return 0;
    }

    fprintf(fp, "i,j,x,y,value,inside_D\n");
    for (int i = 0; i <= M; ++i) {
        double x = X_MIN + i * h1;
        for (int j = 0; j <= N; ++j) {
            double y = Y_MIN + j * h2;
            double value;
            if (i == 0 || i == M || j == 0 || j == N)
                value = v0_function(x, y);
            else
                value = w[idx(i, j, N)];

            fprintf(fp, "%d,%d,%.17g,%.17g,%.17g,%d\n",
                    i, j, x, y, value, point_inside_domain(x, y));
        }
    }
    int ok = !ferror(fp);
    if (fclose(fp) != 0) ok = 0;
    return ok;
}

int main(int argc, char **argv)
{
    int M = 80;
    int N = 80;
    int threads = 1;
    double tolerance = 1.0e-8;
    const char *output_file = "solution_omp.csv";

    if (argc == 2 || argc > 6) {
        fprintf(stderr, "Usage: %s M N threads tolerance output.csv (use - to skip CSV)\n", argv[0]);
        return EXIT_FAILURE;
    }

    /* Запуск: ./solver_omp M N threads tolerance output.csv */
    if (argc >= 3) {
        M = parse_positive_int(argv[1]);
        N = parse_positive_int(argv[2]);
    }
    if (argc >= 4)
        threads = parse_positive_int(argv[3]);
    if (argc >= 5)
        tolerance = parse_tolerance(argv[4]);
    if (argc >= 6)
        output_file = argv[5];

    if (M < 2 || N < 2) {
        fprintf(stderr, "M and N must be >= 2\n");
        return EXIT_FAILURE;
    }
    if (threads < 1) {
        fprintf(stderr, "Number of threads must be >= 1\n");
        return EXIT_FAILURE;
    }

    omp_set_dynamic(0);
    omp_set_num_threads(threads);

    int actual_threads = 1;
#pragma omp parallel
    {
#pragma omp single
        actual_threads = omp_get_num_threads();
#ifdef __linux__
        /* Placement diagnostics outside the measured solver. */
#pragma omp critical
        printf("Placement     : thread=%d cpu=%d\n", omp_get_thread_num(), sched_getcpu());
#endif
    }

    double h1 = (X_MAX - X_MIN) / M;
    double h2 = (Y_MAX - Y_MIN) / N;
    double h = (h1 > h2) ? h1 : h2;
    double epsilon = h * h; /* По условию epsilon=h^2 */
    if ((size_t)(M - 1) > SIZE_MAX / (size_t)(N - 1))
        return EXIT_FAILURE;
    size_t unknowns = (size_t)(M - 1) * (size_t)(N - 1);

    double *q = xcalloc(unknowns, sizeof(double));
    double *diag = xcalloc(unknowns, sizeof(double));
    double *b = xcalloc(unknowns, sizeof(double));
    double *w = xcalloc(unknowns, sizeof(double));
    double *r = xcalloc(unknowns, sizeof(double));
    double *z = xcalloc(unknowns, sizeof(double));
    double *p = xcalloc(unknowns, sizeof(double));
    double *Ap = xcalloc(unknowns, sizeof(double));
    double *previous_w = xcalloc(unknowns, sizeof(double));


    printf("Program       : OpenMP\n");
    printf("Threads       : %d\n", actual_threads);
    printf("Variant       : 5\n");
    printf("Domain        : A(0,0) B(3,0) C(3,2) D(1,3)\n");
    printf("Pi            : [%g,%g] x [%g,%g]\n", X_MIN, X_MAX, Y_MIN, Y_MAX);
    printf("Grid          : %d x %d\n", M, N);
    printf("Unknowns      : %zu\n", unknowns);
    printf("h1            : %.12e\n", h1);
    printf("h2            : %.12e\n", h2);
    printf("epsilon       : %.12e\n", epsilon);
    printf("Tolerance     : %.12e\n", tolerance);

    /* Построение системы вне измеряемого итерационного решателя. */
    double reconstructed_area = build_system(M, N, h1, h2, epsilon, q, diag, b);
    /* Internal solve-only wall time; setup, allocation and output excluded. */
    double start_time = omp_get_wtime();
    SolverResult result = conjugate_gradient(M, N, h1, h2, diag, b, w, tolerance, r, z, p, Ap, previous_w);
    double algorithm_time = omp_get_wtime() - start_time;

    apply_A(w, Ap, diag, M, N, h1, h2);
    double true_residual2 = 0.0, solution2 = 0.0, solution_max = 0.0;
    double functional = 0.0;
    for (size_t k = 0; k < unknowns; ++k) {
        double residual = b[k] - Ap[k];
        true_residual2 += residual * residual;
        solution2 += w[k] * w[k];
        solution_max = fmax(solution_max, fabs(w[k]));
        functional += (0.5 * Ap[k] - b[k]) * w[k];
    }
    printf("True residual : %.12e\n", sqrt(h1 * h2 * true_residual2));
    printf("Solution E norm: %.12e\n", sqrt(h1 * h2 * solution2));
    printf("Solution C norm: %.12e\n", solution_max);
    printf("Functional J  : %.12e\n", h1 * h2 * functional);

    printf("Area(D)       : %.12f (expected 6.500000000000)\n", reconstructed_area);
    printf("Iterations    : %zu\n", result.iterations);
    printf("Restarts      : %zu\n", result.restarts);
    printf("Converged     : %s\n", result.converged ? "yes" : "no");
    printf("Step norm     : %.12e\n", result.step_norm);
    printf("Residual norm : %.12e\n", result.residual_norm);
    printf("Timing scope  : iterative_solver_only_v2\n");
    printf("Algorithm time: %.6f s\n", algorithm_time);

    int output_ok = 1;
    if (strcmp(output_file, "-") != 0)
        output_ok = write_solution(output_file, w, M, N, h1, h2);
    printf("Output        : %s\n", output_file);

    /* Формат: M,N,threads,iterations,time,step_norm,residual_norm */
    printf("RESULT,%d,%d,%d,%zu,%.9f,%.12e,%.12e\n",
           M, N, actual_threads, result.iterations, algorithm_time,
           result.step_norm, result.residual_norm);

    free(q);
    free(diag);
    free(b);
    free(w);
    free(r);
    free(z);
    free(p);
    free(Ap);
    free(previous_w);
    return result.converged && output_ok ? EXIT_SUCCESS : EXIT_FAILURE;
}
