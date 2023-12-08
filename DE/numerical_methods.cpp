#include <iostream>
#include <cmath>
#include <vector>
#include <algorithm>
#include <functional>

template <class T>
using results = std::vector<std::pair<T, double>>;

template <class T>
using system_results = std::vector<std::tuple<T, T, double>>;

[[nodiscard]] constexpr inline double ith(
        const double h,
        const double i
) noexcept { return h * i; }

[[nodiscard]] auto exact(const int n) noexcept {
    const auto h = M_PI / (n - 1);
    results<double> sols(n, { 0, 0 });

    for (int i = 1; i < n; ++i) {
        const auto x = ith(h, i);
        sols[i] = { x, std::sin(x) };
    }

    return sols;
}

[[nodiscard]] auto euler(const int n) noexcept {
    const auto h = M_PI / (n - 1);
    results<double> sols(n, { 0, 0 });

    for (int i = 1; i < n; ++i) {
        const auto [x_prev, y_prev] = sols[i - 1];
        sols[i] = { ith(h, i), y_prev + h * std::cos(x_prev) };
    }

    return sols;
}

[[nodiscard]] auto euler_sys(const int n) noexcept {
    const auto h = M_PI / (n - 1);
    system_results<double> sols(n, { 0, 1, 2 });

    for (int i = 1; i < n; ++i) {
        const auto [t_prev, x_prev, y_prev] = sols[i - 1];
        sols[i] = { ith(h, i), x_prev + h * y_prev, y_prev + h * -4 * x_prev };
    }

    return sols;
}

[[nodiscard]] auto improved_euler(const int n) noexcept {
    const auto h = M_PI / (n - 1);
    results<double> sols(n, { 0, 0 });

    for (int i = 1; i < n; ++i) {
        const auto x = ith(h, i);
        const auto [x_prev, y_prev] = sols[i - 1];
        const auto k1 = std::cos(x_prev);
        const auto k2 = std::cos(x);
        sols[i] = { x, y_prev + h / 2 * (k1 + k2) };
    }

    return sols;
}

[[nodiscard]] auto improved_euler_sys(const int n) noexcept {
    const auto h = M_PI / (n - 1);
    system_results<double> sols(n, { 0, 1, 2 });

    for (int i = 1; i < n; ++i) {
        const auto t = ith(h, i);
        const auto [t_prev, x_prev, y_prev] = sols[i - 1];

        const auto k1x = y_prev;
        const auto k1y = -4 * x_prev;

        const auto k2x = y_prev + h * k1y;
        const auto k2y = -4 * (x_prev + h * k1x);

        const auto x = x_prev + h / 2 * (k1x + k2x);
        const auto y = y_prev + h / 2 * (k1y + k2y);

        sols[i] = { t, x, y };
    }

    return sols;
}

[[nodiscard]] auto runge_kutta(const int n) noexcept {
    const auto h = M_PI / (n - 1);
    results<double> sols(n, { 0, 0 });

    for (int i = 1; i < n; ++i) {
        const auto x = ith(h, i);
        const auto [x_prev, y_prev] = sols[i - 1];
        const auto k1 = std::cos(x_prev);
        const auto k2 = std::cos(x_prev + h / 2);
        const auto k3 = k2;
        const auto k4 = std::cos(x_prev + h);
        sols[i] = { x, y_prev + h / 6 * (k1 + 2 * k2 + 2 * k3 + k4) };
    }

    return sols;
}

[[nodiscard]] auto runge_kutta_sys(const int n) noexcept {
    const auto h = M_PI / (n - 1);
    system_results<double> sols(n, { 0, 1, 2 });

    for (int i = 1; i < n; ++i) {
        const auto t = ith(h, i);
        const auto [t_prev, x_prev, y_prev] = sols[i - 1];

        const auto k1x = y_prev;
        const auto k1y = -4 * x_prev;

        const auto k2x = y_prev + h / 2 * k1y;
        const auto k2y = -4 * (x_prev + h / 2 * k1x);

        const auto k3x = y_prev + h / 2 * k2y;
        const auto k3y = -4 * (x_prev + h / 2 * k2x);

        const auto k4x = y_prev + h * k3y;
        const auto k4y = -4 * (x_prev + h * k3x);

        const auto x = x_prev + h / 6 * (k1x + 2 * k2x + 2 * k3x + k4x);
        const auto y = y_prev + h / 6 * (k1y + 2 * k2y + 2 * k3y + k4y);

        sols[i] = { t, x, y };
    }

    return sols;
}

[[nodiscard]] auto local_errors(results<double>&& result) noexcept {
    results<double> errs(result.size());

    for (int i = 0; i < errs.size(); ++i) {
        const auto [x, y] = result[i];
        errs[i] = { x, std::abs(y - std::sin(x)) };
    }

    return errs;
}

[[nodiscard]] auto global_errors(
        const int n1,
        const int n2,
        std::function<results<double>(int)>&& method
) noexcept {
    results<int> errs(n2 - n1 + 1);

    for (int n = n1; n <= n2; ++n) {
        const auto local_errs = local_errors(method(n));

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
        results<T>&& res,
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
        system_results<T>&& res,
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

void print_solutions(results<double>&& res, const char* const y_msg) noexcept {
    print_results(std::move(res), "ith", y_msg);
}

void print_errors(results<double>&& res, const char* const err_msg) noexcept {
    print_results(local_errors(std::move(res)), "ith", err_msg);
}

void print_global_errors(results<int>&& errs, const char* const ge_msg) {
    print_results(std::move(errs), "ni", ge_msg, "%d ");
}

void print_exact(const int n) {
    print_solutions(exact(n), "y(ith)");
}

void print_euler(const int n) {
    print_solutions(euler(n), "Euler_yi");
}

void print_euler_sys(const int n) {
    print_system_results(euler_sys(n), "Euler_xi", "Euler_yi");
}

void print_improved_euler(const int n) {
    print_solutions(improved_euler(n), "iEuler_yi");
}

void print_improved_euler_sys(const int n) {
    print_system_results(improved_euler_sys(n), "iEuler_xi", "iEuler_yi");
}

void print_runge_kutta(const int n) {
    print_solutions(runge_kutta(n), "RK4_yi");
}

void print_runge_kutta_sys(const int n) {
    print_system_results(runge_kutta_sys(n), "RK4_xi", "RK4_yi");
}

void print_euler_errs(const int n) {
    print_errors(euler(n), "Euler_LE(ith)");
}

void print_improved_euler_errs(const int n) {
    print_errors(improved_euler(n), "iEuler_LE(ith)");
}

void print_runge_kutta_errs(const int n) {
    print_errors(runge_kutta(n), "RK4_LE(ith)");
}

void print_euler_ge(const int n1, const int n2) {
    print_global_errors(
            std::move(global_errors(n1, n2, euler)),
            "Euler_GE(n)"
    );
}

void print_improved_euler_ge(const int n1, const int n2) {
    print_global_errors(
            std::move(global_errors(n1, n2, improved_euler)),
            "iEuler_GE(n)"
    );
}

void print_runge_kutta_ge(const int n1, const int n2) {
    print_global_errors(
            std::move(global_errors(n1, n2, runge_kutta)),
            "RK4_GE(n)"
    );
}

int solve_du() {
    int n = 0, n1 = 0, n2 = 0, task = 0;
    std::scanf("%d%d%d%d", &n, &n1, &n2, &task);

    switch (task) {
        case 1:  print_exact(n);                  break;
        case 2:  print_euler(n);                  break;
        case 3:  print_improved_euler(n);         break;
        case 4:  print_runge_kutta(n);            break;
        case 5:  print_euler_errs(n);             break;
        case 6:  print_improved_euler_errs(n);    break;
        case 7:  print_runge_kutta_errs(n);       break;
        case 8:  print_euler_ge(n1, n2);          break;
        case 9:  print_improved_euler_ge(n1, n2); break;
        case 10: print_runge_kutta_ge(n1, n2);    break;
        default:                                  break;
    }

    return 0;
}

int solve_system_du() {
    int n = 0, task = 0;
    std::scanf("%d%d", &n, &task);

    switch (task) {
        case 1:  print_euler_sys(n);          break;
        case 2:  print_improved_euler_sys(n); break;
        case 3:  print_runge_kutta_sys(n);    break;
        default:                              break;
    }

    return 0;
}

int main() {
    return solve_system_du();
}
