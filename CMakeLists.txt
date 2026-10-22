cmake_minimum_required(VERSION 3.16)
project(cpp_dsp_lab LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

# 源码统一按 UTF-8 编译，保证中文输出不乱码
add_compile_options(-finput-charset=UTF-8 -fexec-charset=UTF-8)
add_compile_options(-Wall)

include_directories(include third_party)

# 随阶段推进逐步扩充
set(CORE_SOURCES
  src/fft.cpp
  src/window.cpp
  src/filter.cpp
  src/qpsk.cpp
  src/csv.cpp
  src/console.cpp
)

add_executable(fft_demo apps/fft_demo.cpp ${CORE_SOURCES})
add_executable(filter_demo apps/filter_demo.cpp ${CORE_SOURCES})
add_executable(qpsk_demo apps/qpsk_demo.cpp ${CORE_SOURCES})
add_executable(tests tests/test_all.cpp ${CORE_SOURCES})
