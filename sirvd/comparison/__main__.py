import sys
import subprocess
from pathlib import Path

import pandas as pd
import numpy as np


COMP = ['S', 'I', 'R', 'V', 'D']

DEFAULT_PARAMS  = {
    'n': 1_000_000,
    't_end': 1_000,
    'a':     0.0015,
    'b':     0.25,
    'g':     0.10,
    'd':     0.0005,
    's':     0.005,
    'i0':    10,
    'atol':  0.01,
    'rtol':  1e-10
}

DEFAULT_METHODS = "rk2,rk3,rk4"
DEFAULT_HS      = "1.0,0.5,0.1"


def build_common_params(args: dict) -> list[str]:
    """Собирает общие параметры системы в список"""
    return [
        f"t_end={args.get('t_end', DEFAULT_PARAMS['t_end'])}",
        f"n={args.get('n', DEFAULT_PARAMS['n'])}",
        f"a={args.get('a', DEFAULT_PARAMS['a'])}",
        f"b={args.get('b', DEFAULT_PARAMS['b'])}",
        f"g={args.get('g', DEFAULT_PARAMS['g'])}",
        f"d={args.get('d', DEFAULT_PARAMS['d'])}",
        f"s={args.get('s', DEFAULT_PARAMS['s'])}",
        f"i0={args.get('i0', DEFAULT_PARAMS['i0'])}",
    ]


def parse_kv(argv: list) -> dict:
    """Парсит аргумент вида ключ=значение."""
    d = dict()
    for a in argv:
        if '=' in a:
            k, v = a.split('=', 1)
            d[k.strip().lower()] = v.strip()
    
    return d


def parse_list_string(lst: str) -> list:
    """Парсит строку-список вида '`str`,`str`,`str`...'."""
    return [x.strip() for x in lst.split(',') if x.strip()]


def parse_float_list(lst: str) -> list:
    """Парсит строку-список вида '`float`,`float`,`float`...'"""
    raw = parse_list_string(lst)
    result = [float(x.strip()) for x in raw]
    
    return sorted(result, reverse=True)


def run_cmd(cmd):
    print(" ".join(cmd))
    subprocess.run(cmd, check=True)


def find_sirvd() -> str:
    """Ищет `sirvd.exe` для автоматического запуска моделирования."""
    if Path("./sirvd.exe").exists():
        return "./sirvd.exe"
    if Path("./sirvd").exists():
        return "./sirvd"
    
    raise FileNotFoundError("Unable to find modeling binary file (/sirvd)")


def make_reference(sirvd: str, args: dict, out_dir: Path) -> Path:
    """Запускает эталонный метод RK56 и создает `ref.csv`"""
    ref_csv = out_dir / "rk56_ref.csv"
    cmd = [
        sirvd,
        "method=rk56",
        f"atol={args.get('atol', DEFAULT_PARAMS['atol'])}",
        f"rtol={args.get('rtol', DEFAULT_PARAMS['rtol'])}",
        f"out={ref_csv}",
        *build_common_params(args)
    ]
    
    run_cmd(cmd) # запускаем моделирование эталона
    return ref_csv


def make_methods_runs(sirvd: str, args: dict, out_dir: Path, methods: list[str], hs: list[float]) -> list:
    """ Запускает моделирование для всех методов `methods` и шагов `hs`.\n
        Возвращает список результатов `.csv`
    """
    processed = []
    common = build_common_params(args)
    
    def _format_h(h: float) -> str:
        return str(h).replace('.', '_')
    
    for m in methods:
        if m == 'rk56':
            continue
        for h in hs:
            out_csv = out_dir / f"{m}_{_format_h(h)}.csv"
            cmd = [
                sirvd,
                f"method={m}",
                f"h={h}",
                f"out={out_csv}",
                *common
            ]
            
            run_cmd(cmd)
            processed.append((m, h, out_csv))
    
    return processed


def linear_interpolation(t_ref, y_ref, t_query):
    """
        Линейная интерполяция.
        
        :param t_ref: Массив точек эталона
        :param y_ref: Массив значений эталона
        :param t_query: Точки, в которых нужно вычислить значение
    """
    return np.interp(t_query, t_ref, y_ref)


def compare_methods(ref_csv, method_csv, epsilon=1e-14):
    """Сравнивает метод с эталоном по данным из `.csv` и возвращает накопившуюся ошибку."""
    reference = pd.read_csv(ref_csv)
    method = pd.read_csv(method_csv)
    
    # эталон
    t_ref = reference['t'].to_numpy()
    
    # метод для сравнения
    t_method = method['t'].to_numpy()

    abs_sum, rel_sum = 0.0, 0.0
    for c in COMP:
        y_ref = reference[c].to_numpy()
        y_method = method[c].to_numpy()
        
        y_ref_interp = linear_interpolation(t_ref, y_ref, t_method)
        
        # абсолютная
        abs_err = np.abs(y_method - y_ref_interp)
        abs_sum += float(np.sum(abs_err))
        
        # относительная
        rel_err = abs_err / (np.abs(y_ref_interp) + epsilon)
        rel_sum += float(np.sum(rel_err))
    
    return abs_sum, rel_sum, len(t_method)


def compare_all(files: list, ref_csv: Path):
    print(f"{'Method':^10} {'h':^10} {'Points':^8} {'Absolute sum error (e_abs_sum)':^35} {'Relative sum error (e_rel_sum)':^35}")
    
    for file in files:
        method_name, h_value, filename = file
        e_abs_sum, e_rel_sum, points = compare_methods(ref_csv, filename)

        print(f"{method_name.upper():^10} {h_value:^10} {points:^8} {e_abs_sum:^35.15e} {e_rel_sum:^35.15e}")


def clear_dir(dir: Path):
    for p in dir.glob('*.csv'):
        p.unlink()


def main():
    if any(x in ('-h', '--help') for x in sys.argv[1:]):
        print(
            "Usage:\n"
            "\tpython -m comparison [methods] [hs] (out_dir) (model params)\n\n"
            "Example:\n"
            "\tpython -m comparison methods=rk2,rk4 hs=1.5,1.0,0.5,0.01 out_dir=cmp_out n=500000 t_end=365 g=0.05"
        )
        sys.exit(0)
    
    args = parse_kv(sys.argv[1:])
    
    methods = parse_list_string(args.get("methods", DEFAULT_METHODS))
    hs      = parse_float_list(args.get("hs", DEFAULT_HS))
    
    out_dir = Path(args.get("out_dir", "comparison_out"))
    out_dir.mkdir(exist_ok=True)
    
    sirvd = find_sirvd()
    
    print("\nModeling reference...")
    ref_csv = make_reference(sirvd, args, out_dir)
    
    print("\nModeling your methods...")
    files = make_methods_runs(sirvd, args, out_dir, methods, hs)
    
    print("\nComparing all methods...")
    compare_all(files, ref_csv)
    
    clear_dir(out_dir)


if __name__ == "__main__":
    main()
