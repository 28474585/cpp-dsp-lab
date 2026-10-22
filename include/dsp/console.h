#pragma once

// 控制台编码设置。
//
// 背景：源码按 UTF-8 编译，字符串字面量在可执行文件里就是 UTF-8 字节。
// 但 Windows 控制台默认代码页是 936(GBK)，直接把 UTF-8 字节打出去会显示乱码。
// 这里把控制台输出代码页切成 UTF-8，让中文正常显示。非 Windows 平台为空操作。

namespace dsp {

void enableUtf8Console();

}  // namespace dsp
