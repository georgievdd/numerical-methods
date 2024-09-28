#include <vector>
#include <stdexcept>
#include <iterator>
#include <iostream>
#include <random>
#include <cassert>

using std::vector;
using std::out_of_range;

namespace matrix {
    template <typename T>
    class Matrix {
    private:
        vector<vector<T>> _data;
    public:
        Matrix() = default;
        Matrix(size_t n, size_t m) : _data(n, vector<T>(m)) {}
        Matrix(const Matrix& other) : _data(other._data) {}
        Matrix(Matrix&& other) noexcept : _data(std::move(other._data)) {}
        virtual ~Matrix() = default;
        Matrix(std::initializer_list<std::initializer_list<T>> init) {
            for (const auto& iterable1 : init) {
                vector<T> row;
                for (const T& item : iterable1) {
                    row.emplace_back(item);
                }
                _data.emplace_back(std::move(row));
            }
        }

        static Matrix random(size_t n, size_t m) {
            Matrix matrix(n, m);
            std::random_device rd;
            std::mt19937 gen(rd());
            std::uniform_real_distribution<T> dis(0.0, 10000.0);

            for (size_t i = 0; i < n; ++i) {
                for (size_t j = 0; j < m; ++j) {
                    matrix._data[i][j] = dis(gen);
                }
            }
            return std::move(matrix);
        }

        static Matrix single(size_t n, size_t m) {
            Matrix matrix(n, m);
            for (size_t i = 0; i < n; ++i) {
                matrix._data[i][i] = 1;
            }
            return std::move(matrix);
        }

        Matrix& operator=(const Matrix& other) {
            if (this != &other) {
                _data.clear();
                _data.resize(other.rowCount(), vector<T>(other.colCount()));
                for (size_t i = 0; i < other.rowCount(); ++i) {
                    for (size_t j = 0; j < other.colCount(); ++j) {
                        _data[i][j] = other._data[i][j];
                    }
                }
            }
            return *this;
        }
        Matrix& operator=(Matrix&& other) noexcept {
            if (this != &other) {
                _data.clear();
                _data = std::move(other._data);
                other._data.clear();
            }
            return *this;
        }
        template <typename Iterable2>
        Matrix& operator=(std::initializer_list<std::initializer_list<T>> init) {
            return *this = Matrix(init);
        }
        const vector<T>& row(size_t index) const {
            if (index >= _data.size()) throw out_of_range("Row index out of range");
            return _data[index];
        }
        vector<T>& row(size_t index) {
            if (index >= _data.size()) throw out_of_range("Row index out of range");
            return _data[index];
        }
        vector<T>& operator [](size_t index) {
            return row(index);
        }
        const vector<T>& operator [](size_t index) const {
            return row(index);
        }
        size_t rowCount() const { return _data.size(); }
        size_t colCount() const { return _data.empty() ? 0 : _data[0].size(); }

        class ColumnIterator {
        private:
            Matrix<T>& _matrix;
            size_t _col;
            size_t _row;

        public:
            using iterator_category = std::forward_iterator_tag;
            using difference_type = std::ptrdiff_t;
            using value_type = T;
            using pointer = T*;
            using reference = T&;

            ColumnIterator(Matrix<T>& matrix, size_t col, size_t row = 0)
                : _matrix(matrix), _col(col), _row(row) {}

            reference operator*() {
                return _matrix._data[_row][_col];
            }

            pointer operator->() {
                return &_matrix._data[_row][_col];
            }

            // Префиксный инкремент
            ColumnIterator& operator++() {
                ++_row;
                return *this;
            }

            // Постфиксный инкремент
            ColumnIterator operator++(int) {
                ColumnIterator tmp = *this;
                ++(*this);
                return tmp;
            }

            bool operator==(const ColumnIterator& other) const {
                return _row == other._row && _col == other._col;
            }

            bool operator!=(const ColumnIterator& other) const {
                return !(*this == other);
            }
        };

        class ConstColumnIterator {
        private:
            const Matrix<T>& _matrix;
            size_t _col;
            size_t _row;

        public:
            using iterator_category = std::forward_iterator_tag;
            using difference_type = std::ptrdiff_t;
            using value_type = const T;
            using pointer = const T*;
            using reference = const T&;

            ConstColumnIterator(const Matrix<T>& matrix, size_t col, size_t row = 0)
                : _matrix(matrix), _col(col), _row(row) {}

            reference operator*() const {
                return _matrix._data[_row][_col];
            }

            pointer operator->() const {
                return &_matrix._data[_row][_col];
            }

            ConstColumnIterator& operator++() {
                ++_row;
                return *this;
            }

            ConstColumnIterator operator++(int) {
                ConstColumnIterator tmp = *this;
                ++(*this);
                return tmp;
            }

            bool operator==(const ConstColumnIterator& other) const {
                return _row == other._row && _col == other._col;
            }

            bool operator!=(const ConstColumnIterator& other) const {
                return !(*this == other);
            }
        };

        ColumnIterator beginColumn(size_t col) {
            if (col >= colCount()) throw out_of_range("Column index out of range");
            return ColumnIterator(*this, col, 0);
        }

        ColumnIterator endColumn(size_t col) {
            if (col >= colCount()) throw out_of_range("Column index out of range");
            return ColumnIterator(*this, col, rowCount());
        }

        ConstColumnIterator beginColumn(size_t col) const {
            if (col >= colCount()) throw out_of_range("Column index out of range");
            return ConstColumnIterator(*this, col, 0);
        }

        ConstColumnIterator endColumn(size_t col) const {
            if (col >= colCount()) throw out_of_range("Column index out of range");
            return ConstColumnIterator(*this, col, rowCount());
        }
    };

    namespace utils {

        constexpr float epsilon = 0.9;

        template <typename L, typename R>
        bool equals(L l, R r) {
            return std::abs((float)l - (float)r) < epsilon;
        }

        template <typename T>
        Matrix<T> transponate(const Matrix<T>& mat) {
            size_t rows = mat.rowCount();
            size_t cols = mat.colCount();
            Matrix<T> transposed(cols, rows);
            for (size_t i = 0; i < rows; ++i) {
                for (size_t j = 0; j < cols; ++j) {
                    transposed[j][i] = mat[i][j];
                }
            }
            return std::move(transposed);
        }

        template <typename T>
        Matrix<T> addColumn(const Matrix<T>& a, const Matrix<T>& b) {
            if (a.rowCount() != b.rowCount()) {
                throw std::invalid_argument("Number of rows must match for adding a column.");
            }
            Matrix<T> result(a.rowCount(), a.colCount() + b.colCount());
            for (size_t i = 0; i < a.rowCount(); ++i) {
                for (size_t j = 0; j < a.colCount(); ++j) {
                    result.row(i)[j] = a.row(i)[j];
                }
            }
            for (size_t i = 0; i < b.rowCount(); ++i) {
                for (size_t j = 0; j < b.colCount(); ++j) {
                    result.row(i)[a.colCount() + j] = b.row(i)[j];
                }
            }
            return std::move(result);
        }

        template <typename T>
        Matrix<T> multiply(const Matrix<T>& m, const Matrix<T>& other) {
            if (m.colCount() != other.rowCount()) {
                throw std::invalid_argument("Matrix dimensions do not allow multiplication");
            }

            Matrix<T> result(m.rowCount(), other.colCount());
            for (size_t i = 0; i < m.rowCount(); ++i) {
                for (size_t j = 0; j < other.colCount(); ++j) {
                    T sum = 0;
                    for (size_t k = 0; k < m.colCount(); ++k) {
                        sum += m.row(i)[k] * other.row(k)[j];
                    }
                    result.row(i)[j] = sum;
                }
            }
            return std::move(result);
        }

        template <typename T, typename Number>
        void multiply(Matrix<T>& m, const Number alpha) {
            for (size_t i = 0; i < m.rowCount(); ++i) {
                for (size_t j = 0; j < m.colCount(); ++j) {
                    m.row(i)[j] *= alpha;
                }
            }
        }

        template <typename T>
        void add(Matrix<T>& m, const Matrix<T>& other) {
            if (m.rowCount() != other.rowCount() || m.colCount() != other.colCount()) {
                throw std::invalid_argument("Matrix dimensions must match for addition");
            }
            for (size_t i = 0; i < m.rowCount(); ++i) {
                for (size_t j = 0; j < m.colCount(); ++j) {
                    m.row(i)[j] += other.row(i)[j];
                }
            }
        }

        template <typename T>
        void triangulateMatrix(Matrix<T>& mat) {
            size_t n = mat.rowCount();
            size_t m = mat.colCount();
            if (n > m) {
                throw std::invalid_argument("Matrix must have at least as many columns as rows.");
            }
            for (size_t i = 0; i < n; ++i) {
                size_t maxRow = i;
                for (size_t k = i + 1; k < n; ++k) {
                    if (std::abs(mat.row(k)[i]) > std::abs(mat.row(maxRow)[i])) {
                        maxRow = k;
                    }
                }
                if (i != maxRow) {
                    std::swap(mat.row(i), mat.row(maxRow));
                }
                T pivot = mat.row(i)[i];
                if (pivot == 0) {
                    throw std::runtime_error("Matrix is singular and cannot be triangulated.");
                }
                for (size_t j = i; j < m; ++j) {
                    mat.row(i)[j] /= pivot;
                }
                for (size_t k = i + 1; k < n; ++k) {
                    T factor = mat.row(k)[i];
                    for (size_t j = i; j < m; ++j) {
                        mat.row(k)[j] -= factor * mat.row(i)[j];
                    }
                }
            }
        }
        
        template <typename T>
        Matrix<T> trianguleMatrix(const Matrix<T>& mat) {
            Matrix<T> copy = mat;
            triangulateMatrix(copy);
            return std::move(copy);
        }
    
        template <typename T>
        float determinant(const Matrix<T>& mat) {
            size_t n = mat.rowCount();
            size_t m = mat.colCount();
            if (n != m) {
                throw std::invalid_argument("Matrix must be square to calculate determinant.");
            }
            
            Matrix<T> tempMat = mat;
            float det = 1;

            for (size_t i = 0; i < n; ++i) {
                size_t maxRow = i;
                for (size_t k = i + 1; k < n; ++k) {
                    if (std::abs(tempMat.row(k)[i]) > std::abs(tempMat.row(maxRow)[i])) {
                        maxRow = k;
                    }
                }
                if (tempMat.row(maxRow)[i] == 0) {
                    return 0;
                }
                if (i != maxRow) {
                    std::swap(tempMat.row(i), tempMat.row(maxRow));
                    det = -det;
                }
                T pivot = tempMat.row(i)[i];
                for (size_t j = i; j < n; ++j) {
                    tempMat.row(i)[j] = (float)tempMat.row(i)[j] / pivot;
                }
                for (size_t k = i + 1; k < n; ++k) {
                    T factor = tempMat.row(k)[i];
                    for (size_t j = i; j < n; ++j) {
                        tempMat.row(k)[j] -= factor * tempMat.row(i)[j];
                    }
                }
                det *= pivot;
            }
            return det;
        }

        template <typename T>
        bool isSingular(const Matrix<T>& mat) {
            return determinant(mat) == 0;
        }
    
        template <typename T>
        int rank(const Matrix<T>& mat) {
            Matrix<T> copy = utils::trianguleMatrix(mat);
            size_t n = copy.rowCount();
            int rank = 0;

            for (size_t i = 0; i < n; ++i) {
                bool nonZeroRow = false;
                for (size_t j = 0; j < copy.colCount(); ++j) {
                    if (!utils::equals(copy.row(i)[j], 0)) {
                        nonZeroRow = true;
                        break;
                    }
                }
                if (nonZeroRow) {
                    ++rank;
                }
            }

            return rank;
        }

    }

    namespace syntactic {
        template <typename T>
        bool operator==(const Matrix<T>& lhs, const Matrix<T>& rhs) {
            if (lhs.rowCount() != rhs.rowCount() || lhs.colCount() != rhs.colCount()) {
                return false;
            }

            for (size_t i = 0; i < lhs.rowCount(); ++i) {
                for (size_t j = 0; j < lhs.colCount(); ++j) {
                    if (!utils::equals(lhs[i][j], rhs[i][j])) {
                        return false;
                    }
                }
            }
            return true;
        }

        template <typename T>
        bool operator!=(const Matrix<T>& lhs, const Matrix<T>& rhs) {
            return !(lhs == rhs);
        }
        template <typename T>
        Matrix<T> operator+(const Matrix<T>& lhs, const Matrix<T>& rhs) {
            Matrix<T> result = lhs;
            utils::add(result, rhs);
            return result;
        }

        template <typename T>
        Matrix<T> operator*(const Matrix<T>& lhs, const Matrix<T>& rhs) {
            return utils::multiply(lhs, rhs);
        }

        template <typename T, typename Number>
        Matrix<T> operator*(const Matrix<T>& m, const Number alpha) {
            Matrix<T> result = m;
            utils::multiply(result, alpha);
            return result;
        }

        template <typename T, typename Number>
        Matrix<T> operator*(const Number alpha, const Matrix<T>& m) {
            return m * alpha;
        }
        template <typename T>
        std::ostream& operator<<(std::ostream& os, const Matrix<T>& matrix) {
            for (size_t i = 0; i < matrix.rowCount(); ++i) {
                for (size_t j = 0; j < matrix.colCount(); ++j) {
                    os << matrix.row(i)[j];
                    if (matrix.colCount() != j + 1)
                        os << " ";
                }
                if (i != matrix.rowCount() - 1)
                    os << "\n";
            }
            return os;
        }
    }
}

void test() {
    using namespace matrix::syntactic;
    using matrix::utils::equals;
    using matrix::Matrix;
    {
        Matrix<float> matrix4 = {{1, 1, 7, 3},
                                {6, 5, 3, 2},
                                {2, 2, 4, 3},
                                {2, 2, 4, 3}};
        assert(matrix::utils::isSingular(matrix4));
    }
    {
        Matrix<float> matrix2 = {{1, 2},
                                 {3, 4}};
        auto determinant = matrix::utils::determinant(matrix2);
        assert(equals(determinant, -2));
    }
}
