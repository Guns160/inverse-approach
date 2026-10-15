#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""GSVD 滤波解法示例：与逐 alpha 正规方程解法的对比（Python 版）。

与同目录 ../matlab/gsvd_example.m 构造同一个合成反问题，输出可逐项对照。

对应仓库正式实现（同一套算术，本文件为精简自包含版）：
- ``inverse_common/gsvd.py``            LAPACKE dggsvd3 紧凑 GSVD（本文件内嵌精简版）
- ``inverse_common/factors.py``         端点消去 + Cholesky 罚项装配（第 1 节同结构）
- ``inverse_common/reconstruction.py``  全 alpha 网格的广义滤波（第 3 节核心公式）

被替换的原始解法见下方引用（原来的代码，解方程部分原文）。

运行：python3 gsvd_example.py
依赖：NumPy。GSVD 需要 LP64 LAPACKE_dggsvd3（自动探测系统 / NumPy / SciPy 的
OpenBLAS，也可用环境变量 JPSI_LAPACKE_LIBRARY 指定库路径）。SciPy 没有暴露
gsvd 接口，所以正式实现与本示例都直接用 ctypes 调 LAPACKE。
"""

# ---------------------------------------------------------------- 原始解法（引用）
# 原来的代码的解方程部分（原文，逐个 alpha 重新求解正规方程）：
#   for i=1:n+mm
#       for j=1:n+mm
#          A2(i,j)=trapz(y,Aphi_New(:,i).*Aphi_New(:,j));
#       end
#       D2(i,1)=trapz(y,Aphi_New(:,i).*g_delta);
#   end
#   count_alphak=1000;
#   alpha_range=logspace(log10(alpha_min),log10(alpha_max),count_alphak);
#   for k=1:count_alphak
#       alpha2(k)=alpha_range(count_alphak-k+1);
#      T=(A2+alpha2(k)*H2);
#      Final= T\ [D2, e1, en];    % 每个 alpha 都重新做一次完整分解
#   end
#   % 更早版本（显式逆矩阵，精度更差）：
#   % FA(:,k)=inv(A2+alpha*H2)*D2;
#   % FA5(:,k)=inv(A2+alpha2(k)*H2)*(D2-Lag);
#
# 正规方程 A2=A'*W*A 把条件数平方；选 alpha（L-curve / 准最优准则）需要对整条
# 网格逐点求解，1000 个点就是 1000 次完整分解。

import ctypes
import ctypes.util
import glob
import os
import sys
import time
from functools import lru_cache
from pathlib import Path

import numpy as np

_LAPACK_COL_MAJOR = 102
_LAPACK_INT = ctypes.c_int
_DOUBLE_PTR = ctypes.POINTER(ctypes.c_double)
_INT_PTR = ctypes.POINTER(_LAPACK_INT)


# ------------------------------------------- 改进核心 1：紧凑 GSVD（LAPACKE）
def _candidate_libraries():
    """返回可能导出 LAPACKE_dggsvd3 的动态库（精简自 inverse_common/gsvd.py）。"""
    candidates = []
    for env in ("JPSI_LAPACKE_LIBRARY", "ETA_C_LAPACKE_LIBRARY"):
        explicit = os.environ.get(env)
        if explicit:
            candidates.append(explicit)
    for name in ("lapacke", "lapack", "openblas", "blas"):
        found = ctypes.util.find_library(name)
        if found:
            candidates.append(found)
    for prefix in (Path(sys.prefix), Path("/usr"), Path("/usr/local"),
                   Path("/opt/homebrew/opt/lapack")):
        for directory in (prefix / "lib", prefix / "lib64"):
            for pattern in ("liblapacke.so*", "liblapack.so*", "libopenblas.so*",
                            "libblas.so*", "liblapacke*.dylib", "liblapack*.dylib",
                            "libopenblas*.dylib"):
                candidates.extend(sorted(glob.glob(str(directory / pattern))))
    numpy_root = Path(np.__file__).resolve().parent
    for directory in (numpy_root.parent / "numpy.libs", numpy_root / ".dylibs",
                      numpy_root.parent / ".dylibs"):
        if directory.is_dir():
            for pattern in ("*lapack*", "*openblas*", "*blas*"):
                candidates.extend(sorted(glob.glob(str(directory / pattern))))
    for entry in sys.path:
        for directory in (Path(entry) / "scipy" / ".dylibs", Path(entry) / "scipy.libs"):
            if directory.is_dir():
                for pattern in ("*openblas*.dylib", "*openblas*.so*", "*lapack*.so*"):
                    candidates.extend(sorted(glob.glob(str(directory / pattern))))
    seen = set()
    return [item for item in candidates
            if item and not (item in seen or seen.add(item))]


@lru_cache(maxsize=1)
def _load_dggsvd3():
    """定位并配置 LAPACKE dggsvd3；跳过 ILP64（64 位整数）的 OpenBLAS。"""
    for candidate in _candidate_libraries():
        try:
            library = ctypes.CDLL(candidate)
            routine = library.LAPACKE_dggsvd3
        except (OSError, AttributeError):
            continue
        config = getattr(library, "openblas_get_config", None)
        if config is not None:
            config.restype = ctypes.c_char_p
            if b"USE64BITINT" in (config() or b""):
                continue  # ILP64 接口与 LP64 不兼容
        routine.restype = _LAPACK_INT
        routine.argtypes = [
            _LAPACK_INT, ctypes.c_char, ctypes.c_char, ctypes.c_char,
            _LAPACK_INT, _LAPACK_INT, _LAPACK_INT, _INT_PTR, _INT_PTR,
            _DOUBLE_PTR, _LAPACK_INT, _DOUBLE_PTR, _LAPACK_INT,
            _DOUBLE_PTR, _DOUBLE_PTR, _DOUBLE_PTR,
            _LAPACK_INT, _DOUBLE_PTR, _LAPACK_INT, _DOUBLE_PTR,
            _LAPACK_INT, _INT_PTR,
        ]
        return routine, candidate
    return None


def _compact_gsvd_square(r, b_mat):
    """对方阵对 (R, B) 调 LAPACKE dggsvd3，按 Regularization Tools cgsvd 的
    排序返回（正则化模在前、不受罚模在后）。算术同 inverse_common/gsvd.py。"""
    loaded = _load_dggsvd3()
    if loaded is None:
        raise RuntimeError(
            "找不到 LP64 LAPACKE_dggsvd3。可用 JPSI_LAPACKE_LIBRARY 指定库路径；"
            "正式实现同样不设正规方程回退。")
    routine, library_path = loaded
    m, n = r.shape
    p, n_b = b_mat.shape
    if n_b != n or m < n or p >= n:
        raise ValueError("compact M3 GSVD 需要 m>=n>p")
    a_work = np.array(r, dtype=np.float64, order="F", copy=True)
    b_work = np.array(b_mat, dtype=np.float64, order="F", copy=True)
    alpha_out = np.empty(n, dtype=np.float64)
    beta_out = np.empty(n, dtype=np.float64)
    u = np.empty((m, m), dtype=np.float64, order="F")
    v = np.empty((p, p), dtype=np.float64, order="F")
    q = np.empty((n, n), dtype=np.float64, order="F")
    iwork = np.empty(n, dtype=np.intc)
    k = _LAPACK_INT()
    ell = _LAPACK_INT()
    info = routine(
        _LAPACK_COL_MAJOR, b"U", b"V", b"Q", m, n, p,
        ctypes.byref(k), ctypes.byref(ell),
        a_work.ctypes.data_as(_DOUBLE_PTR), m,
        b_work.ctypes.data_as(_DOUBLE_PTR), p,
        alpha_out.ctypes.data_as(_DOUBLE_PTR),
        beta_out.ctypes.data_as(_DOUBLE_PTR),
        u.ctypes.data_as(_DOUBLE_PTR), m,
        v.ctypes.data_as(_DOUBLE_PTR), p,
        q.ctypes.data_as(_DOUBLE_PTR), n,
        iwork.ctypes.data_as(_INT_PTR),
    )
    if info != 0:
        raise RuntimeError(f"LAPACKE_dggsvd3 失败，INFO={info}")
    if k.value + ell.value != n or ell.value != p or k.value != n - p:
        raise RuntimeError(f"意外的秩结构：k={k.value}, l={ell.value}")
    if not np.all(np.isfinite(alpha_out)) or not np.all(np.isfinite(beta_out)):
        raise RuntimeError("dggsvd3 返回了非有限的广义奇异值")
    raw_r = np.triu(a_work[:n, :n])
    if np.any(np.abs(np.diag(raw_r)) <= np.finfo(float).tiny):
        raise RuntimeError("dggsvd3 返回了奇异的上三角因子")
    compact_raw = np.arange(k.value, n)   # LAPACK 把不受罚模放在最前 k 列
    null_raw = np.arange(k.value)
    ordering = np.r_[compact_raw, null_raw]
    x_raw = q @ np.linalg.solve(raw_r, np.eye(n))
    return {
        "U": u[:, :n][:, ordering],       # 对应 (R, B) 的 U，正则化模在前
        "X": x_raw[:, ordering],
        "sigma": alpha_out[compact_raw].copy(),
        "mu": beta_out[compact_raw].copy(),
        "sigma_unreg": alpha_out[null_raw].copy(),
        "library": library_path,
    }


def compact_gsvd(a, b):
    """economy QR + LAPACKE dggsvd3（同 inverse_common/gsvd.py 的 qr_reduce 路径）。

    高矩阵 A 先做 economy QR：min ||Ax-b|| 与 min ||Rx-Q'b|| 只差一个与 alpha
    无关的正交残差，因此 (R, B) 的紧凑 GSVD 正是 Tikhonov 解所需的分解。
    返回的 U 已乘回 Q，因此 beta = U.T @ b 直接可用。
    """
    q_a, r_a = np.linalg.qr(a, mode="reduced")
    dec = _compact_gsvd_square(r_a, b)
    dec["U"] = q_a @ dec["U"]
    return dec


# ---------------------------------------------------------------- 1) 合成反问题
def trapezoid_weights(t):
    """均匀网格梯形权（正式程序 inverse_common/numerics.py trapw 的均匀情形）。"""
    dt = t[1] - t[0]
    w = np.full(t.size, dt)
    w[0] = w[-1] = dt / 2.0
    return w


def mass_plus_stiffness(s):
    """帽基质量 Gram（M）+ 一阶差分刚度（M1），H = M + M1。

    构造与原来的代码中的 M、M1、H1 相同（h 为帽基步长）。
    """
    h = s[1] - s[0]
    n_hats = s.size
    big_m = np.zeros((n_hats, n_hats))
    big_m1 = np.zeros((n_hats, n_hats))
    big_m[0, 0] = big_m[-1, -1] = 4.0
    big_m[0, 1] = big_m[1, 0] = big_m[-1, -2] = big_m[-2, -1] = 2.0
    big_m1[0, 0] = big_m1[-1, -1] = 1.0
    big_m1[0, 1] = big_m1[1, 0] = big_m1[-1, -2] = big_m1[-2, -1] = -1.0
    for j in range(1, n_hats - 1):
        big_m[j, j] = 8.0
        big_m[j, j - 1] = big_m[j, j + 1] = 2.0
        big_m1[j, j] = 2.0
        big_m1[j, j - 1] = big_m1[j, j + 1] = -1.0
    return big_m * (h / 12.0) + big_m1 * (1.0 / h)


def build_problem():
    """构造与正式程序同结构的合成反问题（同 factors.py / factorsOneWindow.m）。

    - A = sqrt(w)·K_free：极点列（不受罚）+ 内部帽基核列，两端帽列固定；
    - B = [0, chol(H_ii)]：谱质量+刚度矩阵内部块的 Cholesky 因子作罚项行。
    本示例边界取 0；真实管道中非零边界通过 q = H_ii^{-1}(H_if·boundary) 移位
    并入 b0（见 factors.py 的 b0 表达式）。
    """
    m2_min, m2_max, dm2 = 5.0, 20.0, 0.05      # Borel 网格 t（数据点）
    s_h, lam, ds = 4.0, 20.0, 0.1              # 谱 s 帽基网格
    masses = np.array([3.0969, 3.6861])        # 示例极点：J/psi、psi(2S)
    t = np.arange(m2_min, m2_max + 1e-9, dm2)
    s = np.arange(s_h, lam + 1e-9, ds)
    w = trapezoid_weights(t)
    poles = np.exp(-np.square(masses) / t[:, None]) / t[:, None]
    khat = np.exp(-np.outer(t, s) / 60.0) * ds  # 光滑示意核（正式程序用解析 hatKernel）
    h_mat = mass_plus_stiffness(s)
    n_hats = s.size
    interior = n_hats - 2
    h_ii = h_mat[1:-1, 1:-1]
    r_b = np.linalg.cholesky(h_ii).T           # 上三角，同 MATLAB chol(...,'upper')
    b_mat = np.hstack([np.zeros((interior, masses.size)), r_b])
    k_full = np.hstack([poles, khat])           # 极点列在前，帽列在后
    n_full = k_full.shape[1]
    free = np.r_[:masses.size, masses.size + 1:n_full - 1]  # 极点 + 内部帽
    a_mat = np.sqrt(w)[:, None] * k_full[:, free]
    x_true = np.concatenate([
        np.array([1.2, 0.45]),
        0.8 * np.exp(-np.square(s[1:-1] - 12.0) / 12.0)])
    b0 = a_mat @ x_true                         # 无噪声数据（重点在解法对比）
    return a_mat, b_mat, b0, h_ii


def main():
    # 个别 BLAS（如 Apple Silicon 上的 OpenBLAS）会在 dgemv/dgemm 结果完全正常时
    # 留下陈旧的浮点异常标志，触发 numpy 的虚假 RuntimeWarning；统一屏蔽。
    # 关键结果另有显式对照（第 4 节与直接解、解析范数对比）。
    np.seterr(divide="ignore", over="ignore", invalid="ignore")
    a_mat, b_mat, b0, h_ii = build_problem()
    m2 = a_mat.shape[0]
    n = a_mat.shape[1]
    p = b_mat.shape[0]
    k_poles = n - p
    alpha = np.logspace(-16.0, -2.0, 1001)

    # ---- 2) 基线：原来的代码的逐 alpha 正规方程解法 ------------------------------
    a2 = a_mat.T @ a_mat
    d2 = a_mat.T @ b0
    h2 = np.zeros((n, n))
    h2[k_poles:, k_poles:] = h_ii              # 与 ||B x||^2 完全相同的罚项
    np.linalg.solve(a2 + alpha[0] * h2, d2)     # 预热，使计时处于稳态
    t0 = time.perf_counter()
    x_old = np.empty((n, alpha.size))
    for j in range(alpha.size):
        x_old[:, j] = np.linalg.solve(a2 + alpha[j] * h2, d2)
    t_old = time.perf_counter() - t0

    # ---- 3) 改进：一次 compact GSVD + 全 alpha 网格滤波 --------------------
    compact_gsvd(a_mat, b_mat)                 # 预热（库加载 / BLAS 线程池）
    t0 = time.perf_counter()
    dec = compact_gsvd(a_mat, b_mat)
    t_gsvd = time.perf_counter() - t0
    sigma = dec["sigma"]
    mu = dec["mu"]
    beta_all = dec["U"].T @ b0                 # U 已含 economy-QR 的 Q
    beta_c = beta_all[:p]
    beta_n = beta_all[p:]
    recon = dec["U"] @ beta_all
    null2 = float((b0 - recon) @ (b0 - recon))  # 与 alpha 无关的正交残差

    t0 = time.perf_counter()
    coef = (sigma[:, None] * beta_c[:, None]) / (
        sigma[:, None] ** 2 + alpha[None, :] * mu[:, None] ** 2)  # 广义滤波
    un = beta_n / dec["sigma_unreg"]
    x_new = dec["X"][:, :p] @ coef + (dec["X"][:, p:] @ un)[:, None]
    res_an = np.sqrt(np.sum((sigma[:, None] * coef - beta_c[:, None]) ** 2,
                            axis=0) + null2)    # 解析 ||A x - b0||_w
    pen_an = np.sqrt(np.sum((mu[:, None] * coef) ** 2, axis=0))  # 解析 ||B x||
    t_filter = time.perf_counter() - t0

    # ---- 4) 验证与对比 ----------------------------------------------------
    check_idx = np.linspace(0, alpha.size - 1, 5).astype(int)
    dev_r = dev_p = 0.0
    print(f"问题规模: m2={m2} 数据点, n={n} 未知量({k_poles} 极点 + {n - k_poles} 内部帽), "
          f"p={p} 罚项行, alpha 网格 {alpha.size} 点")
    print(f"GSVD 库: {dec['library']}")
    print(f"cond(A)={np.linalg.cond(a_mat):.2e}; 正规方程把条件数平方到≈"
          f"{np.linalg.cond(a_mat) ** 2:.1e}，浮点下 A'A 已分辨不出小奇异方向"
          f"（实测 cond(A'A)={np.linalg.cond(a2):.2e}）")
    print(f"耗时: 旧法 {1e3 * t_old:.1f} ms; 新法 = GSVD {1e3 * t_gsvd:.1f} ms + "
          f"全网格滤波 {1e3 * t_filter:.1f} ms; 加速比 "
          f"{t_old / (t_gsvd + t_filter):.1f} 倍")
    print()
    print("  alpha        旧法相对误差   新法相对误差（相对增广最小二乘直接解）")
    for j in check_idx:
        a = alpha[j]
        stacked = np.vstack([a_mat, np.sqrt(a) * b_mat])
        rhs = np.concatenate([b0, np.zeros(p)])
        x_ref = np.linalg.lstsq(stacked, rhs, rcond=None)[0]
        scale = max(np.linalg.norm(x_ref), np.finfo(float).tiny)
        e_old = np.linalg.norm(x_old[:, j] - x_ref) / scale
        e_new = np.linalg.norm(x_new[:, j] - x_ref) / scale
        r_direct = np.linalg.norm(a_mat @ x_new[:, j] - b0)
        p_direct = np.linalg.norm(b_mat @ x_new[:, j])
        dev_r = max(dev_r, abs(r_direct - res_an[j])
                    / max(r_direct, np.finfo(float).tiny))
        dev_p = max(dev_p, abs(p_direct - pen_an[j])
                    / max(p_direct, np.finfo(float).tiny))
        print(f"{a:8.1e}   {e_old:10.3e}   {e_new:10.3e}")
    rel_grid = (np.linalg.norm(x_new - x_old, axis=0)
                / np.maximum(np.linalg.norm(x_new, axis=0), np.finfo(float).tiny))
    j_max = int(np.argmax(rel_grid))
    print(f"解析范数自检: 残差范数最大相对偏差 {dev_r:.2e}, "
          f"罚项范数最大相对偏差 {dev_p:.2e}")
    print(f"全网格新旧解最大相对差 {rel_grid[j_max]:.2e}（出现在 alpha={alpha[j_max]:.1e}；"
          f"小 alpha 处旧法失真）")


if __name__ == "__main__":
    main()
