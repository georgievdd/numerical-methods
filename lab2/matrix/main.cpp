#include <sstream>
#include <iomanip>
#include "slae.cpp"
#include "markup.cpp"
using namespace std;
using namespace matrix::syntactic;
using namespace matrix;
#define NM(x, y) "(" << x << "," << y << ")"
int main() {
    stringstream ss;
    for (int n : {6, 12}) {
        const auto coefficients = Matrix<float>::random(n, n);
        freopen(((string)"output" + to_string(n) + ".md").c_str(), "w", stdout);
        const auto results = matrix::Matrix<float>::random(n, 1);
        const auto expandedMatrix = utils::addColumn(results, coefficients);
        if (utils::rank(coefficients) != utils::rank(expandedMatrix)) {
            header << "matrix is singular" << endl;
        }
        header << "random coefficients of matrix" << NM(n, n) << ":" << endl;
        ss = stringstream();
        ss << coefficients;
        code << ss.str() << endl;
        header << "random results of matrix " << NM(n, n) << ":" << endl;
        ss = stringstream();
        ss << results;
        code << ss.str() << endl;
        header << "determinant: " << endl;
        paragraph << utils::determinant(coefficients) << endl;
        header << "rank: " << utils::rank(expandedMatrix) << endl;
        header << "triangled matrix:" << endl;
        ss = stringstream();
        ss << utils::trianguleMatrix(expandedMatrix);
        code << ss.str() << endl;
        header << "solution: " << endl;
        ss = stringstream();
        ss << SLAE(coefficients, results).solve();
        code << ss.str() << endl;
    }

    return 0;
}