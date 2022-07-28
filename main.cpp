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

			double mili = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
			printf("%zu, %lf\n", n, mili);
		}
	} else if (argv[argc - 1][0] == 'e') {
		// 流石に遅すぎるので範囲を狭めている
		for (size_t n = 10; n <= 30; n += 5) {
			matrix_t a = make_random_matrix(n, true);

			auto start = std::chrono::system_clock::now();
			for (size_t i = 0; i < ATTEMPT; ++i) {
				auto x = a.lutishauser(); // 固有値

				for (const double& lambda: x) {
					auto v = a.inviter(lambda); // 固有ベクトル
				}
			}
			auto end = std::chrono::system_clock::now();

			double mili = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
			printf("%zu, %lf\n", n, mili);
		}
	}
	return 0;
}
