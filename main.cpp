#include <chrono>

#include "matrix.hpp"

// 計測時の実行回数
#define ATTEMPT 1

int main(int argc, char** argv) {
	if (argv[argc - 1][0] == 'd') {
		for (size_t n = 100; n <= 1000; ++n) {
			matrix_t a = make_random_matrix(n);

			auto start = std::chrono::system_clock::now();
			for (size_t i = 0; i < ATTEMPT; ++i) {
				double d = a.determinant();
			}
			auto end = std::chrono::system_clock::now();

			double mili = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
			printf("%zu, %lf\n", n, mili);
		}
	} else if (argv[argc - 1][0] == 'l') {
		for (size_t n = 100; n <= 1000; ++n) {
			matrix_t a = make_random_matrix(n);
			std::vector<double> b = make_random_vector(n);

			auto start = std::chrono::system_clock::now();
			for (size_t i = 0; i < ATTEMPT; ++i) {
				std::vector<double> x = a.inverse() * b; // 逆行列で解く
//				std::vector<double> x = a.linsolve_lu(b); // LU 分解で解く。
			}
			auto end = std::chrono::system_clock::now();

			// 誤差評価
			std::vector<double> x = a.inverse() * b; // 逆行列で解く
//			std::vector<double> x = a.linsolve_lu(b); // LU 分解で解く。

			std::vector<double> y = a * x;
			double eps = 0; // 二乗和誤差
			for (size_t i = 0; i < n; ++i) eps += (y[i] - b[i]) * (y[i] - b[i]);

			size_t mili = (size_t)std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
			printf("%zu, %zu %.12lf\n", n, mili, eps);
		}
	} else if (argv[argc - 1][0] == 'e') {
		for (size_t n = 10; n <= 100; n += 10) {
			matrix_t a = make_random_matrix(n, true);

			auto start = std::chrono::system_clock::now();
			for (size_t i = 0; i < ATTEMPT; ++i) {
				auto x = a.lutishauser(); // 固有値

				for (const double& lambda: x) {
					auto v = a.inviter(lambda); // 固有ベクトル

					auto w = a * v;
				}
			}
			auto end = std::chrono::system_clock::now();

			// トレースと固有値の総和の誤差
			auto x = a.lutishauser();
			double eps = std::accumulate(x.begin(), x.end(), 0.);
			for (size_t i = 0; i < n; ++i) eps -= a[i * n + i];

			size_t mili = (size_t)std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
			printf("%zu, %zu %.12lf\n", n, mili, eps);
		}
	}
	return 0;
}
