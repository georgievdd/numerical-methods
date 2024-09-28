#include "slae.cpp"
using namespace std;
using namespace matrix::syntactic;
using namespace matrix;
#define NM(x, y) "(" << x << "," << y << ")"
int main() {
    freopen("output.md", "w", stdout);
    for (int n : {10, 20}) {
        const auto coefficients = Matrix<float>::random(n, n);
        const auto results = matrix::Matrix<float>::random(n, 1);
        const auto expandedMatrix = utils::addColumn(results, coefficients);
        if (utils::rank(coefficients) != utils::rank(expandedMatrix)) {
            cout << "## matrix is singular\n";
        }
        cout << "## random coefficients of matrix " << NM(n, n) << ": \n```" << coefficients << "\n```\n";
        cout << "## random results of matrix " << NM(n, n) << ": \n```" << results << "\n```\n";
        cout << "## rank: " << utils::rank(expandedMatrix) << "\n";
        cout << "## triangled matrix:\n```" << utils::trianguleMatrix(expandedMatrix) << "\n```\n";
        cout << "## solution:\n```" << SLAE(coefficients, results).solve() << "\n```\n";
    }

    return 0;
}