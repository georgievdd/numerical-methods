#include <iostream>
#include <vector>
#include <random>
#include <cmath>
using namespace std;

template <typename T>
class Matrix {
private:
    vector<vector<T>> data_;
public:

    Matrix(int n, int m) {
        data_.resize(n, vector<T>(m));
    }

    Matrix(int n, int m, bool random) {
        data_.resize(n, vector<T>(m));
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < m; ++j) {
                data_[i][j] = rand();
            }
        }
    }

    friend std::istream& operator>>(std::istream& in, Matrix<T>& other) {
        for (int i = 0; i < other.size(); ++i) {
            for (int j = 0; j < other[i].size(); ++j) {
                in >> other[i][j];
            }
        }
        return in;
    }
    const vector<T>& operator[](int i) const {
        return data_[i];
    }

    size_t size() const {
        return data_.size();
    }

    const vector<vector<T>>& payload() const {
        return data_;
    }
};


class Math {
public:
    template <typename T>
    static T evklNorm(const Matrix<T>& matrix) {
        T result = 0;
        for (const auto& v: matrix.payload()) {
            for (const auto& e: v) {
                result += e;
            }
        }
        return sqrt(result);
    }
    template <typename T>
    static T kubNorm(const Matrix<T>& matrix) {
        T maxRowSum = 0;
        for (const auto& row : matrix.payload()) {
            T rowSum = 0.0;
            for (T val : row) {
                rowSum += fabs(val);
            }
            maxRowSum = max(maxRowSum, rowSum);
        }
        return maxRowSum;
    }
};

int main() {

    Matrix<float> matrix(2, 2, true);

    cout << Math::evklNorm(matrix) << "\n"
         << Math::kubNorm(matrix) << endl;


    return 0;
}