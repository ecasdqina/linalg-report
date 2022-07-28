#include "matrix.hpp"

#include <cassert>

bool test(const matrix_t& a, const matrix_t& b) {
	for (size_t i = 0; i < a.size() * a.size(); ++i) {
		if (doublecomp(a[i], b[i])) continue;

		return false;
	}
	return true;
}

bool test(const std::vector<double>& a, const std::vector<double>& b) {
	for (size_t i = 0; i < a.size(); ++i) {
		if (doublecomp(a[i], b[i])) continue;

		return false;
	}
	return true;
}

int main() {
	{
		matrix_t a({{0, 1, 1}, {1, 0, 1}, {1, 1, 0}});

		assert (doublecomp(a.determinant(), 2));
		assert (test(a.inverse(), matrix_t({{-0.5, 0.5, 0.5}, {0.5, -0.5, 0.5}, {0.5, 0.5, -0.5}})));
	}

	{
		matrix_t a({{8, 16, 24, 32}, {2, 7, 12, 17}, {6, 17, 32, 59}, {7, 22, 46, 105}});
		auto [l, u] = a.lu_decomposition();
		assert (test(a, l * u));
	}

	{
		matrix_t a({{1, 0, -1}, {0, -1, 0}, {-1, 0, 1}});
		assert (test(a.lutishauser(), {2, -1, 0}));
		assert (test(a * std::vector<double>({1, 2, 3}), {-2, -2, 2}));
	}

	{
		matrix_t a({{8, 16, 24, 32}, {2, 7, 12, 17}, {6, 17, 32, 59}, {7, 22, 46, 105}});

		auto x = a.linsolve_lu({1, 2, 3, 4});
		assert (test(a * x, {1, 2, 3, 4}));
	}

	{
		matrix_t a({{8, 16, 24, 32}, {2, 7, 12, 17}, {6, 17, 32, 59}, {7, 22, 46, 105}});

		auto x = a.lutishauser();
		for (double v: x) {
			auto y = a.inviter(v);

			auto z = a * y;
			for (double& e: z) e /= v;

			assert (test(z, y));
		}
	}

	printf("TEST = Ok\n");
	return 0;
}
