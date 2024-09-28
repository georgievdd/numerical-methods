#include "matrix.cpp"

namespace matrix {
    template <typename T>
    class SLAE {
    private:
        Matrix<T> _coefficients;
        Matrix<T> _results;     
    public:
        SLAE(const Matrix<T>& coefficients, const Matrix<T>& results)
            : _coefficients(coefficients), _results(results) {
            if (coefficients.rowCount() != results.rowCount()) {
                throw std::invalid_argument("Number of rows in coefficients and results must match.");
            }
        }

        Matrix<T> solve() {
            size_t n = _coefficients.rowCount();
            Matrix<T> augmentedMatrix(n, n + 1);

            for (size_t i = 0; i < n; ++i) {
                for (size_t j = 0; j < n; ++j) {
                    augmentedMatrix[i][j] = _coefficients[i][j];
                }
                augmentedMatrix[i][n] = _results[i][0];
            }

            utils::triangulateMatrix(augmentedMatrix);

            Matrix<T> solution(1, n);
            for (int i = n - 1; i >= 0; --i) {
                if (utils::equals(augmentedMatrix.row(i)[i], 0)) {
                    throw std::runtime_error("Matrix is singular or nearly singular.");
                }
                solution[0][i] = augmentedMatrix.row(i)[n] / augmentedMatrix.row(i)[i];
                for (size_t k = i + 1; k < n; ++k) {
                    solution[0][i] -= augmentedMatrix.row(i)[k] * solution[0][k] / augmentedMatrix.row(i)[i];
                }
            }

            return utils::transponate(solution);
        }
    };
}


void testSLAE() {
    using namespace matrix::syntactic;
    using namespace matrix;
    for (int i{2}; i < 21; ++i) {
        const auto coefficients = matrix::Matrix<float>::random(i, i);
        const auto results = matrix::Matrix<float>::random(i, 1);
        if (utils::rank(coefficients) != utils::rank(utils::addColumn(coefficients, results))) {
            std::cout << "Matrix is singular\n";
            continue;
        }
        matrix::SLAE slae(coefficients, results);
        const auto solution = slae.solve();
        if ((coefficients * solution) != results) {
            std::cout << "results:\n";
            std::cout << results << std::endl;
            std::cout << "coefficients * solution:\n";
            std::cout << coefficients * solution << std::endl;
        }
    }
}