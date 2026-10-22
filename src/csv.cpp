#include "dsp/csv.h"

#include <filesystem>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <stdexcept>

namespace dsp {

namespace {

std::string escape(const std::string& s) {
    const bool needQuote = s.find_first_of(",\"\n\r") != std::string::npos;
    if (!needQuote) return s;
    std::string out = "\"";
    for (char c : s) {
        if (c == '"') out += '"';
        out += c;
    }
    out += '"';
    return out;
}

void writeLine(std::ostream& os, const CsvRow& row) {
    for (std::size_t i = 0; i < row.size(); ++i) {
        if (i) os << ',';
        os << escape(row[i]);
    }
    os << '\n';
}

}  // namespace

void writeCsv(const std::string& path, const CsvRow& header, const std::vector<CsvRow>& rows) {
    const std::filesystem::path p(path);
    if (p.has_parent_path()) std::filesystem::create_directories(p.parent_path());

    std::ofstream os(p, std::ios::binary);
    if (!os) throw std::runtime_error("无法写入 CSV: " + path);

    writeLine(os, header);
    for (const auto& r : rows) writeLine(os, r);
}

std::string fmt(double v, int precision) {
    std::ostringstream ss;
    ss << std::fixed << std::setprecision(precision) << v;
    return ss.str();
}

std::string fmtInt(long long v) { return std::to_string(v); }

}  // namespace dsp
