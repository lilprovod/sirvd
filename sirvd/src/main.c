/*Main proc*/
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <ctype.h>

#include "rk.h"
#include "sirvd.h"

typedef enum {
    METHOD_RK2,
    METHOD_RK3,
    METHOD_RK4,
    METHOD_RK56_ADAPTIVE,
} Method;

static int streq(const char* a, const char* b)
{
    /**
     * @brief Сравнивает строки
     */
    while (*a && *b) {
        if (tolower( (unsigned char)*a ) != tolower( (unsigned char)*b )) return 0;
        a++;
        b++;
    }
    return (*a == '\0' && *b == '\0');
}

static int parse_method(const char* s, Method* out)
{
    if (streq(s, "rk2"))   { *out = METHOD_RK2;   return 1; }
    if (streq(s, "rk3"))   { *out = METHOD_RK3;   return 1; }
    if (streq(s, "rk4"))   { *out = METHOD_RK4;   return 1; }
    if (streq(s, "rk56"))  { *out = METHOD_RK56_ADAPTIVE; return 1; }
    
    return 0;
}

static int is_kv(const char* s) 
{
    /**
     * @brief Проверяет, является ли строка парой ключ-значение
     */
    return strchr(s, '=') != NULL;
}

static int parse_kv(const char* s, char* key, size_t ksz, char* val, size_t vsz)
{
    /**
     * @brief Парсит пару ключ-значение
     */
    const char* eq = strchr(s, '=');
    if (!eq) return 0;

    size_t key_len = (size_t)(eq - s);
    size_t val_len = strlen(eq + 1);

    if (key_len == 0 || val_len == 0) return 0;
    if (key_len >= ksz || val_len >= vsz) return 0;

    memcpy(key, s, key_len);
    key[key_len] = '\0';

    memcpy(val, eq + 1, val_len + 1);

    return 1;
}


int main(int argc, char** argv)
{
    /// Default params
    Method method = METHOD_RK4;
    double h = 0.1;
    double t_end = 1000.0;

    double atol = 1.0;
    double rtol = 1e-6;

    const char* out_name = "output.csv";

    /// SIRVD model
    // Parameters
    SIRVDParams p;
    p.N = 1000000.0;
    p.alpha = 0.0015;
    p.beta  = 0.25;
    p.gamma = 0.10;
    p.delta = 0.0005;
    p.sigma = 0.005;

    // At start
    double x[SIRVD_DIM];

    x[0] = p.N - 10.0;
    x[1] = 10.0;            // I0
    x[2] = 0.0;             // R0
    x[3] = 0.0;             // V0
    x[4] = 0.0;             // D0

    /// Парсинг аргументов
    for (int i = 1; i < argc; i++) {
        const char* arg = argv[i];

        if (is_kv(arg)) {
            char key[64], val[256];

            if (!parse_kv(arg, key, sizeof(key), val, sizeof(val))) {
                printf("Bad argument: %s\n", arg);
                return 1;
            }

            if (streq(key, "method") || streq(key, "m")) {
                Method m;

                if (!parse_method(val, &m)) {
                    printf("Unknown method: %s\nUse: rk2 | rk3 | rk4 | rk56\n", val);
                    return 1;
                }

                method = m;
            } else if (streq(key, "h")) {
                h = atof(val);
            } else if (streq(key, "t_end") || streq(key, "time") || streq(key, "t")) {
                t_end = atof(val);
            } else if (streq(key, "o") || streq(key, "out") || streq(key, "output") || streq(key, "output_file")) {
                out_name = argv[i] + (strchr(arg, '=') - arg) + 1;

            } else if (streq(key, "atol")) {
                atol = atof(val);
            } else if (streq(key, "rtol")) {
                rtol = atof(val);

            } else if (streq(key, "N")) {
                p.N = atof(val);
                
            } else if (streq(key, "alpha") || streq(key, "a")) {
                double alpha = atof(val);
                p.alpha = clamp(alpha, 0.0, 1.0);
            } else if (streq(key, "beta") || streq(key, "b")) {
                double beta = atof(val);
                p.beta = clamp(beta, 0.0, 1.0);
            } else if (streq(key, "gamma") || streq(key, "g")) {
                double gamma = atof(val);
                p.gamma = clamp(gamma, 0.0, 1.0);
            } else if (streq(key, "delta") || streq(key, "d")) {
                double delta = atof(val);
                p.delta = clamp(delta, 0.0, 1.0);
            } else if (streq(key, "sigma") || streq(key, "s")) {
                double sigma = atof(val);
                p.sigma = clamp(sigma, 0.0, 1.0);

            } else if (streq(key, "i0") || streq(key, "i_start")) {
                double I0 = atof(val);
                x[1] = clamp(I0, 0.0, p.N);
            } else if (streq(key, "r0") || streq(key, "r_start")) {
                double R0 = atof(val);
                x[2] = clamp(R0, 0.0, p.N);
            } else if (streq(key, "v0") || streq(key, "v_start")) {
                double V0 = atof(val);
                x[3] = clamp(V0, 0.0, p.N);
            } else if (streq(key, "d0") || streq(key, "d_start")) {
                double D0 = atof(val);
                x[4] = clamp(D0, 0.0, p.N);
            } else {
                printf("Unknown option: %s\n", key);
                return 1;
            }
        }
    }

    double from_all_groups = x[1] + x[2] + x[3] + x[4];
    if (from_all_groups > p.N) {
        printf("Bad initialize: I0 + R0 + V0 + D0 = %.10f > N = %.10f\n", from_all_groups, p.N);
        return 1;
    }
    x[0] = p.N - from_all_groups;       // S0

    /// Валидация
    if (h <= 0.0)     { printf("Bad h: %f\nMust be positive\n", h); return 1; }
    if (t_end <= 0.0) { printf("Bad t_end: %f\nMust be positive\n", t_end); return 1;   }
    if (atol <= 0.0)  { printf("Bad atol: %f\nMust be positive\n", atol); return 1;     }
    if (rtol < 0.0)   { printf("Bad rtol: %f\nMust be non-negative\n", rtol); return 1; }


    /// File
    FILE* f = fopen(out_name, "w");
    if (!f) {
        printf("Cannot open output file: %s\n", out_name);
        return 1;
    }
    fprintf(f, "t,h,h_min,h_max,S,I,R,V,D\n");


    /// Основной цикл (прогноз)
    double t = 0.0;

    while (t < t_end)
    {
        fprintf(f, "%.10f,%.10f,%.10f,%.10f,%.10f,%.10f,%.10f,%.10f,%.10f\n",
                t, h, H_MIN, H_MAX, x[0], x[1], x[2], x[3], x[4]);
        
        switch (method) {
            case METHOD_RK2:
                rk2_step(t, h, x, SIRVD_DIM, sirvd_rhs_wrap, &p);   break;
            case METHOD_RK3:
                rk3_step(t, h, x, SIRVD_DIM, sirvd_rhs_wrap, &p);   break;
            case METHOD_RK4:
                rk4_step(t, h, x, SIRVD_DIM, sirvd_rhs_wrap, &p);   break;
            case METHOD_RK56_ADAPTIVE:
                while (!rk56_adaptive_step(
                    &t, &h, H_MIN, H_MAX, x, SIRVD_DIM, sirvd_rhs_wrap, &p, t_end, atol, rtol))
                {
                    //
                }
                continue;
        }

        t += h;
    }

    fclose(f);

    printf("Job done! Output written to: %s\n", out_name);
    return 0;
}