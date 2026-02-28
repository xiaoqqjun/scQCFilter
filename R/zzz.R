# zzz.R - 包初始化文件
# 这个文件在包加载时自动执行

# 确保所有依赖包都被加载
.onLoad <- function(libname, pkgname) {
  # 依赖包会自动被加载（通过NAMESPACE中的import）
  # 这里无需额外操作
}

.onAttach <- function(libname, pkgname) {
  packageStartupMessage("scQCFilter loaded successfully!")
}
