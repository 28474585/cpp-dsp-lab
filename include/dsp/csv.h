#pragma once

#include <string>
#include <vector>

namespace dsp {

using CsvRow = std::vector<std::string>;

// 写 UTF-8 CSV（带表头）。目录不存在会自动创建。
void writeCsv(const std::string& path, const CsvRow& header, const std::vector<CsvRow>& rows);

// 把 double 格式化成固定小数位的字符串，避免 CSV 里出现 1e+06 之类的科学计数法。
std::string fmt(double v, int precision = 6);

std::string fmtInt(long long v);

}  // namespace dsp
