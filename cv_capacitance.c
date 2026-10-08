/*
 * cv_capacitance.c — Supercapacitor CV Analysis in Pure C
 * ---------------------------------------------------------------
 * Reads Cyclic Voltammetry (CV) files exported by a CHI600E
 * electrochemical workstation, integrates the current-potential
 * curve with the trapezoidal rule, and computes the double-layer
 * capacitance C = Q / dV, where Q = (1/v) * INTEGRAL(I dE).
 *
 * Two experiments are analysed:
 *   (A) scan-rate series   -> rate capability (capacitance retention)
 *   (B) electrolyte series -> electrolyte screening at fixed scan rate
 *
 * Pure C99, no third-party library. Build:
 *     gcc -O2 -Wall -o cv_capacitance.exe cv_capacitance.c -lm
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define MAX_PATH   512
#define MAX_LINE   1024
#define MAX_FILES  32
#define INIT_CAP   4096      /* initial capacity of the dynamic array */

/* One measured point: potential (V), current (A) */
typedef struct {
    double e;                /* potential, V */
    double i;                /* current,   A */
} Point;

/* Parsed result of one CV file */
typedef struct {
    char   name[64];
    double e_high;           /* upper potential limit, V   */
    double e_low;            /* lower potential limit, V   */
    double scan_rate;        /* V/s                        */
    double sensitivity;      /* A/V                        */
    int    npoints;
    double q_charge;         /* charge of forward scan, C  */
    double capacitance;      /* F                          */
} CVRun;

/* ---------------------------------------------------------------
 * grab_double: pull a numeric value out of a CHI header line.
 * e.g. "Scan Rate (V/s) = 0.1"  -> 0.1
 * Returns 1 on success, 0 if the key was not found.
 * --------------------------------------------------------------- */
static int grab_double(const char *line, const char *key, double *out)
{
    if (strncmp(line, key, strlen(key)) != 0) return 0;
    const char *eq = strchr(line, '=');
    if (!eq) return 0;
    *out = strtod(eq + 1, NULL);
    return 1;
}

/* ---------------------------------------------------------------
 * read_cv_file: parse one CHI600E CV text export.
 * Header carries the experiment parameters; data section starts
 * after the "Potential/V, Current/A" column title.
 * Points are stored in a dynamic array (realloc doubling) — this is
 * exactly the "growable sequential list" from Data Structures.
 * --------------------------------------------------------------- */
static int read_cv_file(const char *path, const char *label, CVRun *run)
{
    FILE *fp = fopen(path, "r");
    if (!fp) { fprintf(stderr, "  [skip] cannot open %s\n", path); return 0; }

    char line[MAX_LINE];
    int  in_data = 0, n = 0, cap = INIT_CAP;
    Point *pts = (Point *)malloc(sizeof(Point) * cap);
    if (!pts) { fclose(fp); return 0; }

    memset(run, 0, sizeof(*run));
    snprintf(run->name, sizeof(run->name), "%s", label);
    run->scan_rate = 0.0;

    while (fgets(line, sizeof(line), fp)) {
        if (!in_data) {
            grab_double(line, "High E (V)",      &run->e_high);
            grab_double(line, "Low E (V)",       &run->e_low);
            grab_double(line, "Scan Rate (V/s)", &run->scan_rate);
            grab_double(line, "Sensitivity (A/V)",&run->sensitivity);
            if (strstr(line, "Potential/V") && strstr(line, "Current/A")) in_data = 1;
            continue;
        }
        /* data line: "-0.001, -8.834e-3" */
        double e, i;
        if (sscanf(line, "%lf, %lf", &e, &i) != 2) continue;
        if (n >= cap) {                       /* array full -> double it */
            cap *= 2;
            Point *tmp = (Point *)realloc(pts, sizeof(Point) * cap);
            if (!tmp) { free(pts); fclose(fp); return 0; }
            pts = tmp;
        }
        pts[n].e = e;
        pts[n].i = i;
        n++;
    }
    fclose(fp);

    if (n < 10 || run->scan_rate <= 0.0) { free(pts); return 0; }
    run->npoints = n;

    /* -----------------------------------------------------------
     * Locate the forward (anodic) scan: from the lowest potential
     * to the highest. Segment order of a CHI 3-segment CV is
     *   Init E -> Low E -> High E -> Init E,
     * so the forward half-cycle is [idx_min, idx_max].
     * ----------------------------------------------------------- */
    int idx_min = 0, idx_max = 0;
    for (int k = 1; k < n; k++) {
        if (pts[k].e < pts[idx_min].e) idx_min = k;
        if (pts[k].e > pts[idx_max].e) idx_max = k;
    }
    int lo = (idx_min < idx_max) ? idx_min : idx_max;
    int hi = (idx_min < idx_max) ? idx_max : idx_min;

    /* Trapezoidal rule: INTEGRAL(I dE) over the forward scan */
    double area = 0.0;
    for (int k = lo; k < hi; k++) {
        double dE = pts[k + 1].e - pts[k].e;
        double i_avg = 0.5 * (pts[k].i + pts[k + 1].i);
        area += i_avg * dE;
    }

    /* Q = (1/v) * INTEGRAL(I dE);   C = Q / dV */
    double dV = run->e_high - run->e_low;
    run->q_charge     = fabs(area) / run->scan_rate;
    run->capacitance  = (dV > 0.0) ? run->q_charge / dV : 0.0;

    free(pts);
    return 1;
}

/* qsort comparator: sort by capacitance, descending */
static int cmp_cap_desc(const void *a, const void *b)
{
    const CVRun *x = (const CVRun *)a, *y = (const CVRun *)b;
    if (x->capacitance > y->capacitance) return -1;
    if (x->capacitance < y->capacitance) return  1;
    return 0;
}

int main(void)
{
    /* (A) scan-rate series — 2M H2SO4, eight scan rates */
    const char *rate_files[MAX_FILES] = {
        "cv_data/rate_0.002Vps.txt","cv_data/rate_0.005Vps.txt",
        "cv_data/rate_0.01Vps.txt", "cv_data/rate_0.02Vps.txt",
        "cv_data/rate_0.05Vps.txt", "cv_data/rate_0.1Vps.txt",
        "cv_data/rate_0.2Vps.txt",  "cv_data/rate_0.5Vps.txt"
    };
    const char *rate_labels[MAX_FILES] = {
        "0.002 V/s","0.005 V/s","0.01 V/s","0.02 V/s",
        "0.05 V/s","0.1 V/s","0.2 V/s","0.5 V/s"
    };
    /* (B) electrolyte series — fixed 0.1 V/s */
    const char *elec_files[MAX_FILES] = {
        "cv_data/elec_0.5M_H2SO4.txt","cv_data/elec_1M_HCl.txt",
        "cv_data/elec_1M_KOH.txt",    "cv_data/elec_2M_H2SO4.txt"
    };
    const char *elec_labels[MAX_FILES] = {
        "0.5M H2SO4","1M HCl","1M KOH","2M H2SO4"
    };

    printf("=====================================================\n");
    printf(" Supercapacitor CV Analysis  (CHI600E data, pure C)\n");
    printf(" C = Q/dV,  Q = (1/v) * INTEGRAL(I dE)  [trapezoid]\n");
    printf("=====================================================\n\n");

    /* ---------- (A) rate capability ---------- */
    CVRun runs[MAX_FILES];
    int nr = 0;
    printf("[A] Scan-rate series  (2M H2SO4, window %.1f V)\n",
           (runs[0].e_high - runs[0].e_low) > 0 ? (runs[0].e_high - runs[0].e_low) : 1.6);
    printf("    %-12s %10s %12s %12s %10s\n",
           "scan rate", "points", "Q (C)", "C (F)", "retention");
    printf("    ------------------------------------------------------\n");

    for (int k = 0; k < 8; k++) {
        CVRun r;
        if (!read_cv_file(rate_files[k], rate_labels[k], &r)) continue;
        runs[nr++] = r;
    }
    /* Baseline = the slowest scan, i.e. the closest to equilibrium.
     * Faster scans leave ions less time to reach the pores, so the
     * measured capacitance drops — that drop is the rate capability. */
    int    i_base = 0;
    double base   = 0.0;
    for (int k = 0; k < nr; k++)
        if (runs[k].scan_rate < runs[i_base].scan_rate) i_base = k;
    base = runs[i_base].capacitance;

    FILE *fout = fopen("cv_rate_capacitance.csv", "w");
    if (fout) fprintf(fout, "scan_rate_Vps,points,Q_C,C_F,retention_pct\n");
    for (int k = 0; k < nr; k++) {
        double ret = (base > 0.0) ? 100.0 * runs[k].capacitance / base : 0.0;
        printf("    %-12s %10d %12.6f %12.6f %9.1f%%\n",
               runs[k].name, runs[k].npoints, runs[k].q_charge,
               runs[k].capacitance, ret);
        if (fout) fprintf(fout, "%s,%d,%.8f,%.8f,%.2f\n",
                          runs[k].name, runs[k].npoints,
                          runs[k].q_charge, runs[k].capacitance, ret);
    }
    if (fout) fclose(fout);
    printf("    -> written: cv_rate_capacitance.csv\n\n");

    /* ---------- (B) electrolyte screening ---------- */
    CVRun eruns[MAX_FILES];
    int ne = 0;
    printf("[B] Electrolyte screening  (fixed 0.1 V/s)\n");
    printf("    %-14s %10s %12s %12s\n", "electrolyte", "points", "Q (C)", "C (F)");
    printf("    ------------------------------------------------------\n");

    for (int k = 0; k < 4; k++) {
        CVRun r;
        if (!read_cv_file(elec_files[k], elec_labels[k], &r)) continue;
        eruns[ne++] = r;
    }
    qsort(eruns, ne, sizeof(CVRun), cmp_cap_desc);   /* rank by capacitance */

    fout = fopen("cv_electrolyte.csv", "w");
    if (fout) fprintf(fout, "electrolyte,points,Q_C,C_F\n");
    for (int k = 0; k < ne; k++) {
        printf("    %-14s %10d %12.6f %12.6f\n",
               eruns[k].name, eruns[k].npoints, eruns[k].q_charge, eruns[k].capacitance);
        if (fout) fprintf(fout, "%s,%d,%.8f,%.8f\n",
                          eruns[k].name, eruns[k].npoints,
                          eruns[k].q_charge, eruns[k].capacitance);
    }
    if (fout) fclose(fout);
    printf("    -> written: cv_electrolyte.csv\n\n");

    if (ne > 0) {
        printf("    best electrolyte: %s  (%.6f F)\n", eruns[0].name, eruns[0].capacitance);
    }
    if (nr >= 2) {
        double fast = runs[nr - 1].capacitance;
        printf("    rate capability : %.1f%% retained at %s vs %s\n",
               100.0 * fast / base, runs[nr - 1].name, runs[0].name);
    }
    printf("\nNote: active-mass loading was not recorded in the raw files,\n");
    printf("      so results are capacitance C (F), not specific capacitance (F/g).\n");
    return 0;
}
