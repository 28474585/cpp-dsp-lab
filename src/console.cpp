#include "dsp/console.h"

#ifdef _WIN32
#include <windows.h>
#endif

namespace dsp {

void enableUtf8Console() {
#ifdef _WIN32
    // 切到 UTF-8 代码页，与 -fexec-charset=UTF-8 的输出保持一致。
    // 同时尝试切输入代码页，避免读取到 GBK 编码的参数时出错。
    SetConsoleOutputCP(CP_UTF8);
    SetConsoleCP(CP_UTF8);
#endif
}

}  // namespace dsp
