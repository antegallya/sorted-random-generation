#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <random>
#include <vector>
#include <map>

template < class Generator >
void naive(std::vector<double> & x, double sample(Generator&), Generator& g) {
    for (size_t i = 0; i < x.size(); ++i) {
        x[i] = sample(g);
    }

    std::sort(x.begin(), x.end());
}

template < class Generator >
using method_type = decltype(naive<Generator>);

template < class Generator >
void exp_dist(std::vector<double> & x, double sample(Generator&),
              Generator& g) {
    double sum = 0.;
    for (size_t i = 0; i < x.size(); ++i) {
        /* sample is sampling [0, 1), but we need (0, 1]. */
        sum -= std::log(1. - sample(g));
        x[i] = sum;
    }
    sum -= std::log(1. - sample(g));
    for (size_t i = 0; i < x.size(); ++i) {
        x[i] = x[i] / sum;
    }
}

template < class Generator >
void bentley_saxe(std::vector<double> & x, double sample(Generator&),
                  Generator& g) {
    double curmax = 1.;
    for (size_t i = x.size(); i > 0; --i) {
        curmax *= std::pow(1. - sample(g), 1. / i);
        x[i - 1] = curmax;
    }
}

template < class Generator >
void bentley_saxe_expln(std::vector<double> & x, double sample(Generator&),
                  Generator& g) {
    double curmax = 1.;
    for (size_t i = x.size(); i > 0; --i) {
        /* exp(log(x)/n) is faster than pow(x, 1/n). */
        curmax *= std::exp(std::log(1. - sample(g)) / i);
        x[i - 1] = curmax;
    }
}

void check_distrib(std::vector<double> x) {
    std::string sneg = std::is_sorted(x.begin(), x.end()) ? "" : "un";
    std::cout << "Output is " << sneg << "ordered" << std::endl;
    /* Add a uniformity test ? */

    auto x_limits = std::minmax_element(x.begin(), x.end());
    std::cout << "X in [" << *(x_limits.first) << ','
              << (double) *(x_limits.second) << ']' << std::endl;

    std::array<unsigned int, 40> hist;
    hist.fill(0);
    for (auto p : x) {
        ++hist[std::floor(hist.size() * p)];
    }

    auto hist_limits = std::minmax_element(hist.begin(), hist.end());
    unsigned int bmin = *hist_limits.first;
    unsigned int bmax = *hist_limits.second;
    unsigned int bspan = bmax - bmin;
    const unsigned short w = 40;
    std::cout << hist.size() << " buckets. Expected bucket fill ratio "
              << 1. / hist.size() << '.' << std::endl;
    std::cout << std::fixed << std::setprecision(5)
              << "fill ratio:" << std::setw(2) << '|'
              << (double) bmin / x.size() << std::setw(w - 5 - 2)
              << (double) bmax / x.size() << '|' << std::endl;
    for (size_t i = 0; i < hist.size(); ++i) {
        std::cout << std::fixed << std::setprecision(2) << std::setw(2)
                  << '[' << ((double) i) / hist.size() << ','
                  << ((double) (i + 1)) / hist.size() << ')' << ' '
                  << std::string(w * (hist[i] - bmin) / bspan, '*')
                  << std::endl;
    }
}

template < class Generator >
method_type<Generator>& method_fun(const std::string& method) {
    if (method == "naive") return naive<Generator>;
    if (method == "exp_dist") return exp_dist<Generator>;
    if (method == "bentley_saxe") return bentley_saxe<Generator>;
    if (method == "bentley_saxe_expln") return bentley_saxe_expln<Generator>;
    throw std::invalid_argument("Unknown generation method");
}

void usage(const std::string& arg0) {
    std::cerr << "usage: " << arg0 << " <METHOD> <N>" << std::endl
              << "  METHOD is one of naive, exp_dist" << std::endl;
    exit(1);
}

int main(int argc, char* argv[]) {
    if (argc != 3) {
        usage(argv[0]);
    }
    /* I'm too lazy to check this int. */
    size_t n = std::stoi(argv[2]);
    std::random_device rnd_dev;
    std::default_random_engine rnd_engine(rnd_dev());
    method_type<typeof rnd_engine>* method;
    try {
        method = method_fun<typeof rnd_engine>(argv[1]);
    }
    catch (const std::invalid_argument& ia) {
        usage(argv[0]);
    }

    std::vector<double> x(n);
    auto sample = [](auto& g){
        return std::generate_canonical<double, 10>(g);
    };

    auto start = std::chrono::steady_clock::now();
    method(x, sample, rnd_engine);
    auto end = std::chrono::steady_clock::now();

    std::chrono::duration<double> elapsed_seconds = end-start;
    std::cout << "elapsed time: " << elapsed_seconds.count() << 's'
              << std::endl;
    check_distrib(x);
}
