/*
 * liion_cycle.c — 锂离子电池恒流充放电(GCD)数据分析（C 语言实现，自用实测数据）
 *
 * 输入：一个目录，里面是若干 CSV 曲线文件，格式：
 *       第一行表头 "Capacity_mAh_g,Voltage_V"，其后每行 "比容量,电压"
 *       文件名形如 "Sheet1_curve02_CH.csv" / "Sheet2_curve01_DI.csv"
 *         CH = 充电曲线（电压随容量上升），DI = 放电曲线（电压随容量下降）
 *         curve 后面的两位数字是该曲线在第几列（用于按原始顺序配对）
 *
 * 输出：逐循环的比容量、库仑效率、放电容量保持率、极化电压(充放电中点电压差)
 *
 * 与 408 的对应（面试可讲）：
 *   - 动态数组 realloc            -> 顺序表（线性表）
 *   - qsort + 函数指针比较器       -> 排序算法 / 回调
 *   - 线性插值求 50% 容量处电压     -> 数值方法
 *   - 均值 / 保持率 / 极值统计      -> 基础算法
 *   - CSV 解析（strtod / 字符串）  -> 文件与字符串处理
 *
 * 编译： gcc -O2 -Wall -o liion_cycle.exe liion_cycle.c -lm
 * 运行： liion_cycle.exe liion_data/cell01
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <dirent.h>

#define MAX_PATH 512
#define INIT_CAP 16

/* 单条充放电曲线 */
typedef struct {
    double *cap;      /* 比容量 mAh/g */
    double *vol;      /* 电压 V */
    int      n;
    int      index;   /* 解析自文件名，用于按原始列序配对 */
    int      is_charge; /* 1=充电, 0=放电 */
    char     name[128];
} Curve;

/* 曲线集合（动态数组） */
typedef struct {
    Curve *a;
    int    n;
    int    cap;
} CurveList;

/* 单个循环的统计结果 */
typedef struct {
    double q_charge;
    double q_discharge;
    double ce;   /* 库仑效率 %，无对应充电时为 -1 */
    double dv;   /* 极化电压 = 充电中点电压 - 放电中点电压，无效时为 NAN */
} CycleStat;

/* ---------- 顺序表：动态扩容 ---------- */
static void curve_list_push(CurveList *L, Curve c) {
    if (L->n == L->cap) {
        L->cap = L->cap ? L->cap * 2 : INIT_CAP;
        L->a = (Curve *)realloc(L->a, L->cap * sizeof(Curve));
    }
    L->a[L->n++] = c;
}

/* ---------- 读一条曲线 ---------- */
static int read_curve(const char *path, Curve *out) {
    FILE *fp = fopen(path, "r");
    if (!fp) return 0;
    /* 跳过表头 */
    char line[256];
    if (!fgets(line, sizeof(line), fp)) { fclose(fp); return 0; }

    int cap_alloc = 256, cap_n = 0;
    double *ca = (double *)malloc(cap_alloc * sizeof(double));
    double *vo = (double *)malloc(cap_alloc * sizeof(double));
    while (fgets(line, sizeof(line), fp)) {
        double c, v;
        if (sscanf(line, "%lf,%lf", &c, &v) == 2) {
            if (cap_n == cap_alloc) {
                cap_alloc *= 2;
                ca = (double *)realloc(ca, cap_alloc * sizeof(double));
                vo = (double *)realloc(vo, cap_alloc * sizeof(double));
            }
            ca[cap_n] = c;
            vo[cap_n] = v;
            cap_n++;
        }
    }
    fclose(fp);
    if (cap_n < 5) { free(ca); free(vo); return 0; }

    out->cap = ca; out->vol = vo; out->n = cap_n;
    out->is_charge = (vo[cap_n - 1] > vo[0]) ? 1 : 0;
    return 1;
}

/* ---------- 从文件名解析列序号 "xxx_curve07_CH.csv" ---------- */
static int parse_index(const char *name) {
    const char *p = strstr(name, "curve");
    if (!p) return 0;
    p += 5;
    int idx = 0, got = 0;
    while (*p >= '0' && *p <= '9') { idx = idx * 10 + (*p - '0'); p++; got = 1; }
    return got ? idx : 0;
}

/* ---------- qsort 比较器：按列序号升序 ---------- */
static int cmp_index(const void *a, const void *b) {
    int ia = ((const Curve *)a)->index;
    int ib = ((const Curve *)b)->index;
    return ia - ib;
}

/* ---------- 线性插值：求 50% 最大容量处的电压（中点电压） ---------- */
static double midpoint_voltage(const Curve *c) {
    double qmax = c->cap[0];
    for (int i = 1; i < c->n; i++) if (c->cap[i] > qmax) qmax = c->cap[i];
    double target = qmax / 2.0;
    for (int i = 1; i < c->n; i++) {
        double lo = c->cap[i - 1], hi = c->cap[i];
        if ((lo <= target && target <= hi) || (hi <= target && target <= lo)) {
            if (lo == hi) continue;
            double f = (target - lo) / (hi - lo);
            return c->vol[i - 1] + f * (c->vol[i] - c->vol[i - 1]);
        }
    }
    return c->vol[c->n / 2];
}

/* ---------- 读目录内所有 CSV ---------- */
static void read_dir(const char *dir, CurveList *L) {
    DIR *d = opendir(dir);
    if (!d) { fprintf(stderr, "无法打开目录: %s\n", dir); return; }
    struct dirent *e;
    while ((e = readdir(d)) != NULL) {
        size_t len = strlen(e->d_name);
        if (len < 4 || strcmp(e->d_name + len - 4, ".csv") != 0) continue;
        char full[MAX_PATH];
        snprintf(full, sizeof(full), "%s/%s", dir, e->d_name);
        Curve c;
        if (read_curve(full, &c)) {
            c.index = parse_index(e->d_name);
            snprintf(c.name, sizeof(c.name), "%s", e->d_name);
            curve_list_push(L, c);
        }
    }
    closedir(d);
    qsort(L->a, L->n, sizeof(Curve), cmp_index);
}

/* ---------- 分析一个电池（一个目录） ---------- */
static void analyze(const char *label, CurveList *L) {
    printf("\n========================================\n");
    printf("  电池: %s   (曲线数=%d)\n", label, L->n);
    printf("========================================\n");

    /* 按列序配对：遇到充电先暂存，遇到放电则合成一个循环 */
    CycleStat stats[64];
    int nstat = 0;
    const Curve *pending = NULL;
    for (int i = 0; i < L->n; i++) {
        const Curve *c = &L->a[i];
        if (c->is_charge) {
            pending = c;
        } else {
            CycleStat s;
            s.q_discharge = c->cap[0];
            for (int k = 1; k < c->n; k++) if (c->cap[k] > s.q_discharge) s.q_discharge = c->cap[k];
            s.q_charge = -1; s.ce = -1; s.dv = NAN;
            if (pending) {
                double v_ch_end   = pending->vol[pending->n - 1];
                double v_di_start = c->vol[0];
                /* 同一循环判定：放电起点电压不应高于充电终点电压（容差 0.2 V）。
                   若不满足，说明该放电与暂存充电不属于同一循环，不配对，
                   放电按“无对应充电”处理，pending 保留等待后续匹配。 */
                if (v_di_start <= v_ch_end + 0.2) {
                    s.q_charge = pending->cap[0];
                    for (int k = 1; k < pending->n; k++) if (pending->cap[k] > s.q_charge) s.q_charge = pending->cap[k];
                    s.ce = s.q_discharge / s.q_charge * 100.0;
                    double vmch = midpoint_voltage(pending);
                    double vmdi = midpoint_voltage(c);
                    s.dv = vmch - vmdi;
                    pending = NULL;
                }
            }
            if (nstat < 64) stats[nstat++] = s;
        }
    }

    printf("\n  %-6s %-10s %-10s %-10s %-10s\n", "循环", "充电容量", "放电容量", "库仑效率", "极化ΔV");
    printf("  -----------------------------------------------------------\n");
    double first_dis = 0, last_dis = 0;
    int dis_cnt = 0;
    double ce_sum = 0; int ce_cnt = 0;
    double dv_sum = 0; int dv_cnt = 0;
    for (int i = 0; i < nstat; i++) {
        CycleStat s = stats[i];
        char ce_buf[16], dv_buf[16];
        if (s.ce > 0) { snprintf(ce_buf, sizeof(ce_buf), "%.1f%%", s.ce); ce_sum += s.ce; ce_cnt++; }
        else strcpy(ce_buf, "—(无充电)");
        if (!isnan(s.dv)) { snprintf(dv_buf, sizeof(dv_buf), "%.3f V", s.dv); dv_sum += s.dv; dv_cnt++; }
        else strcpy(dv_buf, "—");
        printf("  %-6d %-10.2f %-10.2f %-10s %-10s\n",
               i + 1, s.q_charge > 0 ? s.q_charge : 0, s.q_discharge, ce_buf, dv_buf);
        if (s.q_discharge > 0) {
            if (dis_cnt == 0) first_dis = s.q_discharge;
            last_dis = s.q_discharge;
            dis_cnt++;
        }
    }

    printf("\n  --- 汇总 ---\n");
    if (ce_cnt > 0)
        printf("  平均库仑效率(有效循环): %.1f%%\n", ce_sum / ce_cnt);
    if (dv_cnt > 0)
        printf("  平均极化电压ΔV      : %.3f V\n", dv_sum / dv_cnt);
    if (dis_cnt > 0 && first_dis > 0)
        printf("  放电容量保持率(末/首): %.1f%%  (首放 %.2f -> 末放 %.2f mAh/g, 共 %d 次放电)\n",
               last_dis / first_dis * 100.0, first_dis, last_dis, dis_cnt);
}

int main(int argc, char **argv) {
    if (argc < 2) {
        printf("用法: %s <数据目录> [<数据目录2> ...]\n", argv[0]);
        printf("例:   %s liion_data/cell01\n", argv[0]);
        return 1;
    }
    for (int a = 1; a < argc; a++) {
        CurveList L = {0};
        read_dir(argv[a], &L);
        if (L.n == 0) {
            printf("\n[警告] 目录 %s 中没有读到任何曲线 CSV。\n", argv[a]);
            continue;
        }
        analyze(argv[a], &L);
        for (int i = 0; i < L.n; i++) { free(L.a[i].cap); free(L.a[i].vol); }
        free(L.a);
    }
    return 0;
}
