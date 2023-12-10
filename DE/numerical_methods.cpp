#include <iostream>
#include <cmath>
#include <vector>
#include <algorithm>
#include <functional>

template <class T>
using results = std::vector<std::pair<T, double>>;

template <class T>
using system_results = std::vector<std::tuple<T, T, double>>;

using y_x = std::function<double(double)>;
using f_xy = std::function<double(double, double)>;
using f_txy = std::function<double(double, double, double)>;

using numerical_method = std::function<results<double>(int, const f_xy&, double, double, double)>;

[[nodiscard]] constexpr double y_impl(const double x) noexcept {
    const auto c = -1 / M_E;
    const auto sol = std::sqrt(std::pow(M_E, 1 / x) * c + 2);
    return sol;
}

[[nodiscard]] constexpr double f_impl(const double x, const double y) noexcept {
    return (2 - y * y) / (2 * x * x * y);
}

[[nodiscard]] constexpr double f_sys_impl(
        const double t,
        const double x,
        const double y
) noexcept {
    return 2 + y - x * x;
}

[[nodiscard]] constexpr double g_sys_impl(
        const double t,
        const double x,
        const double y
) noexcept {
    return 2 * (x * x - x * y);
}

[[nodiscard]] constexpr inline double ith(
        const double start,
        const double h,
        const double i
) noexcept { return start + h * i; }

[[nodiscard]] auto exact(
        const int n,
        const y_x& y,
        const double x0,
        const double b
) noexcept {
    const auto h = (b - x0) / (n - 1);
    results<double> sols(n);

    for (int i = 0; i < n; ++i) {
        const auto x = ith(x0, h, i);
        sols[i] = { x, y(x) };
    }

    return sols;
}

[[nodiscard]] auto euler(
        const int n,
        const f_xy& f,
        const double x0,
        const double y0,
        const double b
) noexcept {
    const auto h = (b - x0) / (n - 1);
    results<double> sols(n, { x0, y0 });

    for (int i = 1; i < n; ++i) {
        const auto [x_prev, y_prev] = sols[i - 1];
        sols[i] = { ith(x0, h, i), y_prev + h * (f(x_prev, y_prev)) };
    }

    return sols;
}

[[nodiscard]] auto euler_sys(
        const int n,
        const f_txy& f,
        const f_txy& g,
        const double t0,
        const double x0,
        const double y0,
        const double b
) noexcept {
    const auto h = (b - t0) / (n - 1);
    system_results<double> sols(n, { t0, x0, y0 });

    for (int i = 1; i < n; ++i) {
        const auto [t_prev, x_prev, y_prev] = sols[i - 1];

        sols[i] = {
                ith(t0, h, i),
                x_prev + h * f(t_prev, x_prev, y_prev),
                y_prev + h * g(t_prev, x_prev, y_prev)
        };
    }

    return sols;
}

[[nodiscard]] auto improved_euler(
        const int n,
        const f_xy& f,
        const double x0,
        const double y0,
        const double b
) noexcept {
    const auto h = (b - x0) / (n - 1);
    results<double> sols(n, { x0, y0 });

    for (int i = 1; i < n; ++i) {
        const auto x = ith(x0, h, i);
        const auto [x_prev, y_prev] = sols[i - 1];
        const auto k1 = f(x_prev, y_prev);
        const auto k2 = f(x, y_prev + h * k1);
        sols[i] = { x, y_prev + h / 2 * (k1 + k2) };
    }

    return sols;
}

[[nodiscard]] auto improved_euler_sys(
        const int n,
        const f_txy& f,
        const f_txy& g,
        const double t0,
        const double x0,
        const double y0,
        const double b
) noexcept {
    const auto h = (b - t0) / (n - 1);
    system_results<double> sols(n, { t0, x0, y0 });

    for (int i = 1; i < n; ++i) {
        const auto t = ith(t0, h, i);
        const auto [t_prev, x_prev, y_prev] = sols[i - 1];

        const auto k1x = f(t_prev, x_prev, y_prev);
        const auto k1y = g(t_prev, x_prev, y_prev);

        const auto k2x = f(t, x_prev + h * k1x, y_prev + h * k1y);
        const auto k2y = g(t, x_prev + h * k1x, y_prev + h * k1y);

        const auto x = x_prev + h / 2 * (k1x + k2x);
        const auto y = y_prev + h / 2 * (k1y + k2y);

        sols[i] = { t, x, y };
    }

    return sols;
}

[[nodiscard]] auto runge_kutta(
        const int n,
        const f_xy& f,
        const double x0,
        const double y0,
        const double b
) noexcept {
    const auto h = (b - x0) / (n - 1);
    results<double> sols(n, { x0, y0 });

    for (int i = 1; i < n; ++i) {
        const auto x = ith(x0, h, i);
        const auto [x_prev, y_prev] = sols[i - 1];
        const auto k1 = f(x_prev, y_prev);
        const auto k2 = f(x_prev + h / 2, y_prev + h / 2 * k1);
        const auto k3 = f(x_prev + h / 2, y_prev + h / 2 * k2);
        const auto k4 = f(x_prev + h, y_prev + h * k3);
        sols[i] = { x, y_prev + h / 6 * (k1 + 2 * k2 + 2 * k3 + k4) };
    }

    return sols;
}

[[nodiscard]] auto runge_kutta_sys(
        const int n,
        const f_txy& f,
        const f_txy& g,
        const double t0,
        const double x0,
        const double y0,
        const double b
) noexcept {
    const auto h = (b - t0) / (n - 1);
    system_results<double> sols(n, { t0, x0, y0 });

    for (int i = 1; i < n; ++i) {
        const auto t = ith(t0, h, i);
        const auto [t_prev, x_prev, y_prev] = sols[i - 1];

        const auto k1x = f(t_prev, x_prev, y_prev);
        const auto k1y = g(t_prev, x_prev, y_prev);

        const auto k2x = f(t, x_prev + h / 2 * k1x, y_prev + h / 2 * k1y);
        const auto k2y = g(t, x_prev + h / 2 * k1x, y_prev + h / 2 * k1y);

        const auto k3x = f(t, x_prev + h / 2 * k2x, y_prev + h / 2 * k2y);
        const auto k3y = g(t, x_prev + h / 2 * k2x, y_prev + h / 2 * k2y);

        const auto k4x = f(t, x_prev + h * k3x, y_prev + h * k3y);
        const auto k4y = g(t, x_prev + h * k3x, y_prev + h * k3y);

        const auto x = x_prev + h / 6 * (k1x + 2 * k2x + 2 * k3x + k4x);
        const auto y = y_prev + h / 6 * (k1y + 2 * k2y + 2 * k3y + k4y);

        sols[i] = { t, x, y };
    }

    return sols;
}

[[nodiscard]] auto local_errors(const results<double>& result, const y_x& f) noexcept {
    results<double> errs(result.size());

    for (int i = 0; i < errs.size(); ++i) {
        const auto [x, y] = result[i];
        errs[i] = { x, std::abs(y - f(x)) };
    }

    return errs;
}

[[nodiscard]] auto global_errors(
        const int n1,
        const int n2,
        const numerical_method& method,
        const y_x& yx,
        const f_xy& f,
        const double x0,
        const double y0,
        const double b
) noexcept {
    results<int> errs(n2 - n1 + 1);

    for (int n = n1; n <= n2; ++n) {
        const auto local_errs = local_errors(method(n, f, x0, y0, b), yx);

        const auto [_, max] = *std::max_element(
                local_errs.begin(),
                local_errs.end(),
                [](const auto x, const auto y) {
                    return x.second < y.second;
                }
        );

        errs[n - n1] = { n, max };
    }

    return errs;
}

template <class T> void print_results(
        const results<T>& res,
        const char* first_msg,
        const char* const second_msg,
        const char* const first_arg_fmt = "%.5lf "
) noexcept {
    std::printf("%s=\n", first_msg);

    for (auto&& [x, y] : res)
        std::printf(first_arg_fmt, x);

    std::printf("\n%s=\n", second_msg);

    for (auto&& [x, y] : res)
        std::printf("%.5lf ", y);
}

template <class T> void print_system_results(
        const system_results<T>& res,
        const char* x_msg,
        const char* const y_msg
) noexcept {
    std::puts("ti=");

    for (auto&& [t, x, y] : res)
        std::printf("%.5f ", t);

    std::printf("\n%s=\n", x_msg);

    for (auto&& [t, x, y] : res)
        std::printf("%.5lf ", x);

    std::printf("\n%s=\n", y_msg);

    for (auto&& [t, x, y] : res)
        std::printf("%.5lf ", y);
}

void print_solutions(const results<double>& res, const char* const y_msg) noexcept {
    print_results(res, "xi", y_msg);
}

void print_local_errors(const results<double>& errs, const char* const err_msg) noexcept {
    print_results(errs, "xi", err_msg);
}

void print_global_errors(const results<int>& errs, const char* const ge_msg) {
    print_results(errs, "ni", ge_msg, "%d ");
}

void print_exact(const int n, const double x0, const double b) {
    print_solutions(exact(n, y_impl, x0, b), "y(xi)");
}

void print_euler(const int n, const double x0, const double y0, const double b) {
    print_solutions(euler(n, f_impl, x0, y0, b), "Euler_yi");
}

void print_euler_sys(
        const int n,
        const double t0,
        const double x0,
        const double y0,
        const double b
) {
    print_system_results(
            euler_sys(n, f_sys_impl, g_sys_impl, t0, x0, y0, b),
            "Euler_xi",
            "Euler_yi"
    );
}

void print_improved_euler(const int n, const double x0, const double y0, const double b) {
    print_solutions(improved_euler(n, f_impl, x0, y0, b),"iEuler_yi");
}

void print_improved_euler_sys(
        const int n,
        const double t0,
        const double x0,
        const double y0,
        const double b
) {
    print_system_results(
            improved_euler_sys(n, f_sys_impl, g_sys_impl, t0, x0, y0, b),
            "iEuler_xi",
            "iEuler_yi"
    );
}

void print_runge_kutta(
        const int n,
        const double x0,
        const double y0,
        const double b
) {
    print_solutions(runge_kutta(n, f_impl, x0, y0, b), "RK4_yi");
}

void print_runge_kutta_sys(
        const int n,
        const double t0,
        const double x0,
        const double y0,
        const double b
) {
    print_system_results(
            runge_kutta_sys(n, f_sys_impl, g_sys_impl, t0, x0, y0, b),
            "RK4_xi",
            "RK4_yi"
    );
}

void print_euler_le(const int n, const double x0, const double y0, const double b) {
    print_local_errors(
            local_errors(euler(n, f_impl, x0, y0, b), y_impl),
            "Euler_LE(xi)"
    );
}

void print_improved_euler_le(
        const int n,
        const double x0,
        const double y0,
        const double b
) {
    print_local_errors(
            local_errors(improved_euler(n, f_impl, x0, y0, b), y_impl),
            "iEuler_LE(xi)"
    );
}

void print_runge_kutta_le(
        const int n,
        const double x0,
        const double y0,
        const double b
) {
    print_local_errors(
            local_errors(runge_kutta(n, f_impl, x0, y0, b), y_impl),
            "RK4_LE(xi)"
    );
}

void print_euler_ge(
        const int n1,
        const int n2,
        const double x0,
        const double y0,
        const double b
) {
    print_global_errors(
            global_errors(n1, n2, euler, y_impl, f_impl, x0, y0, b),
            "Euler_GE(n)"
    );
}

void print_improved_euler_ge(
        const int n1,
        const int n2,
        const double x0,
        const double y0,
        const double b
) {
    print_global_errors(
            global_errors(n1, n2, improved_euler, y_impl, f_impl, x0, y0, b),
            "iEuler_GE(n)"
    );
}

void print_runge_kutta_ge(
        const int n1,
        const int n2,
        const double x0,
        const double y0,
        const double b
) {
    print_global_errors(
            global_errors(n1, n2, runge_kutta, y_impl, f_impl, x0, y0, b),
            "RK4_GE(n)"
    );
}

auto read_de_params() {
    int n = 0, n1 = 0, n2 = 0;
    std::scanf("%d%d%d", &n, &n1, &n2);
    return std::make_tuple(n, n1, n2);
}

void solve_de() {
    const auto [n, n1, n2] = read_de_params();

    int task = 0;
    std::scanf("%d", &task);

    const double x0 = 1, y0 = 1, b = 2;

    switch (task) {
        case 1:  print_exact(n, x0, b);                      break;
        case 2:  print_euler(n, x0, y0, b);                  break;
        case 3:  print_improved_euler(n, x0, y0, b);         break;
        case 4:  print_runge_kutta(n, x0, y0, b);            break;
        case 5:  print_euler_le(n, x0, y0, b);               break;
        case 6:  print_improved_euler_le(n, x0, y0, b);      break;
        case 7:  print_runge_kutta_le(n, x0, y0, b);         break;
        case 8:  print_euler_ge(n1, n2, x0, y0, b);          break;
        case 9:  print_improved_euler_ge(n1, n2, x0, y0, b); break;
        case 10: print_runge_kutta_ge(n1, n2, x0, y0, b);    break;
        default:                                             break;
    }
}

void solve_de_sys() {
    int n = 0, task = 0;
    std::scanf("%d%d", &n, &task);

    const double t0 = 0, x0 = -2, y0 = -2, b = 0.55;

    switch (task) {
        case 1:  print_euler_sys(n, t0, x0, y0, b);          break;
        case 2:  print_improved_euler_sys(n, t0, x0, y0, b); break;
        case 3:  print_runge_kutta_sys(n, t0, x0, y0, b);    break;
        default:                                             break;
    }
}

int main() {
    solve_de_sys();
    return 0;
}
