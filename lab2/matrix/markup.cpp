#include <iostream>
#include <string>
#include <sstream>
using namespace std;
namespace markup {
    class MarkupInstance {
        const string prefix, suffix;
    public:
        MarkupInstance() = default;
        MarkupInstance(const string& prefix, const string& suffix = "")
            : prefix(prefix), suffix(suffix) {}
        virtual ~MarkupInstance() = default;
        virtual ostream& operator<<(const string& s) const {
            const auto output = prefix + s + suffix;
            return cout << output;
        }
        virtual ostream& operator<<(const float s) const {
            const auto output = prefix + to_string(s) + suffix;
            return cout << output;
        }
    };
    struct MarkupHeader : MarkupInstance {
        MarkupHeader() : MarkupInstance("# ", "") {}
    };
    struct MarkupCode : MarkupInstance {
        MarkupCode() : MarkupInstance("```\n", "\n```") {}
    };
    struct MarkupParagraph : MarkupInstance {
        MarkupParagraph() : MarkupInstance("", "") {}
    };
}
const markup::MarkupHeader header;
const markup::MarkupCode code;
const markup::MarkupParagraph paragraph;

