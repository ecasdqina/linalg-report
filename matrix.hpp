#include <stdlib.h>
#include <cstring>
#include <algorithm>
#include <stdio.h>
#include <initializer_list>
#include <vector>
#include <numeric>
#include <cmath>
#include <random>

constexpr double EPS = 1e-3;

bool doublecomp(double x, double y = 0) {
	return std::abs(x - y) <= EPS;
}

// ベクトルの正規化
void normalize(std::vector<double>& x) {
	double abs = 0;
	for (size_t i = 0; i < x.size(); ++i) abs += x[i] * x[i];
	abs = sqrt(abs);

	for (size_t i = 0; i < x.size(); ++i) x[i] /= abs;
}

struct matrix_t {
	size_t n;
	double* a;
	double* tmp;

	using value_type = double;

	matrix_t(): n(0), a(nullptr), tmp(nullptr) {}
	matrix_t(size_t n)
	: n(n), a(new double[n * n]), tmp(new double[n]) {}
	matrix_t(size_t n, double v): matrix_t(n) {
		std::memset(a, v, n * n * sizeof(double));
	}
	matrix_t(const matrix_t& b): matrix_t(b.size()) {
		std::memcpy(a, b.a, n * n * sizeof(double));
	}
	matrix_t(const std::vector<std::vector<double>>& init): matrix_t(init.size()) {
		for (size_t i = 0; i < n; ++i) for (size_t j = 0; j < n; ++j) a[i * n + j] = init[i][j];
	}
	~matrix_t() {
		delete[] a;
		delete[] tmp;
	}

	matrix_t operator=(const matrix_t& b) {
		std::memcpy(a, b.a, n * n * sizeof(double));
		return (*this);
	}

	matrix_t operator+(const matrix_t& b) const {
		matrix_t x(*this);
		for (size_t i = 0; i < n * n; ++i) x[i] += b[i];
		return x;
	}
	matrix_t operator-(const matrix_t& b) const {
		matrix_t x(*this);
		for (size_t i = 0; i < n * n; ++i) x[i] -= b[i];
		return x;
	}
	matrix_t operator*(const matrix_t& b) const {
		matrix_t x(size(), 0);
		for (size_t i = 0; i < n; ++i) {
			for (size_t k = 0; k < n; ++k) {
				for (size_t j = 0; j < n; ++j) {
					x[i * n + j] += a[i * n + k] * b[k * n + j];
				}
			}
		}
		return x;
	}
	std::vector<double> operator*(const std::vector<double>& b) const {
		std::vector<double> x(n, 0);
		for (size_t i = 0; i < n; ++i) {
			for (size_t j = 0; j < n; ++j) {
				x[i] += a[i * n + j] * b[j];
			}
		}
		return x;
	}
	matrix_t operator*(double b) const {
		matrix_t x(*this);
		for (size_t i = 0; i < n * n; ++i) x[i] *= b;
		return x;
	}
	matrix_t operator/(double b) const {
		matrix_t x(*this);
		for (size_t i = 0; i < n * n; ++i) x[i] /= b;
		return x;
	}

	// 逆行列
	matrix_t inverse() const {
		matrix_t x(*this), y(n, 0);
		for (size_t i = 0; i < n; ++i) y[i * n + i] = 1;
		for (size_t i = 0; i < n; ++i) {
			if (doublecomp(x[i * n + i])) {
				for (size_t j = i + 1; j < n; ++j) {
					if (doublecomp(x[j * n + i])) continue;

					x.swap_row(i, j);
					y.swap_row(i, j);
					break;
				}
			}
			if (doublecomp(x[i * n + i])) return matrix_t();

			const double div = x[i * n + i];
			for (size_t j = 0; j < n; ++j) {
				x[i * n + j] /= div;
				y[i * n + j] /= div;
			}

			for (size_t j = i + 1; j < n; ++j) {
				const double scal = x[j * n + i];
				for (size_t k = 0; k < n; ++k) {
					x[j * n + k] -= x[i * n + k] * scal;
					y[j * n + k] -= y[i * n + k] * scal;
				}
			}
		}
		for (int i = n - 1; i >= 0; --i) {
			for (size_t j = 0; j < i; ++j) {
				const double scal = x[j * n + i];
				for (size_t k = 0; k < n; ++k) y[j * n + k] -= y[i * n + k] * scal;
				x[j * n + i] = 0;
			}
		}
		return y;
	}

	// 行列式
	double determinant() const {
		double det = 1;
		matrix_t x(*this);
		for (size_t i = 0; i < n; ++i) {
			if (doublecomp(x[i * n + i])) {
				for (size_t j = i + 1; j < n; ++j) {
					if (doublecomp(x[j * n + i])) continue;

					x.swap_row(i, j);
					det *= -1;
					break;
				}
			}
			if (doublecomp(x[i * n + i])) return 0;

			const double div = x[i * n + i];
			for (size_t j = i; j < n; ++j) x[i * n + j] /= div;
			det *= div;

			for (size_t j = i + 1; j < n; ++j) {
				const double scal = x[j * n + i];
				for (size_t k = i + 1; k < n; ++k) {
					x[j * n + k] -= x[i * n + k] * scal;
				}
			}
		}
		return det;
	}

	// LU 分解
	std::pair<matrix_t, matrix_t> lu_decomposition() const {
		matrix_t x(*this), l(n, 0), u(n, 0);
		for (size_t i = 0; i + 1 < n; ++i) {
			const double div = x[i * n + i];
			if (doublecomp(div)) {
				return {matrix_t(), matrix_t()};
			}

			l[i * n + i] = 1;
			u[i * n + i] = div;
			for (size_t j = i + 1; j < n; ++j) {
				l[j * n + i] = x[j * n + i] / div;
				u[i * n + j] = x[i * n + j];
			}

			for (size_t j = i + 1; j < n; ++j) {
				for (size_t k = i + 1; k < n; ++k) {
					x[j * n + k] -= l[j * n + i] * u[i * n + k];
				}
			}
		}
		l[n * n - 1] = 1;
		u[n * n - 1] = x[n * n - 1];

		return {l, u};
	}
	// LU 法
	std::vector<double> lutishauser() const {
		matrix_t x(*this);
		while (!x.is_upper_trimatrix()) {
			auto [l, u] = x.lu_decomposition();

			x = u * l;
		}

		std::vector<double> evs(n);
		for (size_t i = 0; i < n; ++i) evs[i] = x[i * n + i];
		return evs;
	}

	// LU 分解を用いた線形方程式ソルバ
	std::vector<double> linsolve_lu(const std::vector<double>& b) {
		auto [l, u] = lu_decomposition();

		std::vector<double> y(n);
		for (size_t i = 0; i < n; ++i) {
			y[i] = b[i];
			for (size_t j = 0; j < i; ++j) y[i] -= l[i * n + j] * y[j];
//			y[i] /= l[i * n + i]; // bacause l[i * n + i] = 1;
		}

		std::vector<double> x(n);
		for (int i = n - 1; i >= 0; --i) {
			x[i] = y[i];
			for (size_t j = i + 1; j < n; ++j) x[i] -= u[i * n + j] * x[j];
			x[i] /= u[i * n + i];
		}
		return x;
	}

	// 逆反復法
	std::vector<double> inviter(double ev) const {
		matrix_t b(*this);
		for (size_t i = 0; i < n; ++i) b[i * n + i] -= ev;

		std::vector<double> x(n, 1. / sqrt((double)n)), y(n);
		while (true) {
			y = b.linsolve_lu(x);
			x.swap(y);
			normalize(x);

			if (x[0] < 0) {
				for (size_t i = 0; i < n; ++i) x[i] *= -1;
			}

			bool cont = false;
			for (size_t i = 0; i < n; ++i) {
				if (doublecomp(x[i], y[i])) continue;

				cont = true;
			}
			if (!cont) break;
		}
		return x;
	}

	// 上三角行列かどうか
	bool is_upper_trimatrix() const {
		for (size_t i = 0; i < n; ++i) {
			for (size_t j = 0; j < i; ++j) {
				if (doublecomp(a[i * n + j])) continue;

				return false;
			}
		}
		return true;
	}

	void swap_row(size_t i, size_t j) {
		std::memcpy(tmp, a + i * n, n * sizeof(double));
		std::memcpy(a + i * n, a + j * n, n * sizeof(double));
		std::memcpy(a + j * n, tmp, n * sizeof(double));
	}

	size_t size() const { return n; }
	double get(size_t i, size_t j) const { return a[i * n + j]; }
	double& operator[](size_t k) { return a[k]; }
	double operator[](size_t k) const { return a[k]; }

	void print() const {
		if (n == 0) {
			puts("[]");
			fflush(stdout);
		} else {
			for (size_t i = 0; i < n; ++i) {
				printf(i == 0 ? "[[" : " [");
				for (size_t j = 0; j < n; ++j) {
					printf("%lf%s", a[i * n + j], (j + 1 == n ? "]" : " "));
				}
				puts("]");
				fflush(stdout);
			}
		}
	}
};

matrix_t make_random_matrix(size_t n, bool symmetric = false, bool regular = true) {
	std::random_device seed_gen;
	std::default_random_engine engine(seed_gen());
	std::uniform_real_distribution<> dist(0, 1.0);

	matrix_t a(n);
	do {
		for (size_t i = 0; i < n * n; ++i) a[i] = dist(engine);
	} while (regular and doublecomp(a.determinant()));

	if (symmetric) {
		for (size_t i = 0; i < n; ++i) {
			for (size_t j = 0; j < i; ++j) {
				a[i * n + j] = a[j * n + i];
			}
		}
	}
	return a;
}

std::vector<double> make_random_vector(size_t n) {
	std::random_device seed_gen;
	std::default_random_engine engine(seed_gen());
	std::uniform_real_distribution<> dist(-1.0, 1.0);

	std::vector<double> a(n);
	for (double& v: a) v = dist(engine);
	return a;
}
