import os
import sys

import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt

mpl.rcParams["font.family"] = "DejaVu Sans"


def draw_from_csv(filename: str = "output.csv"):
    if not os.path.isfile(filename):
        raise FileNotFoundError(f"Cannot find file for plot: {filename}\n")
    
    df = pd.read_csv(filename)
    
    backgr = "#181818"
    accent = "#A855F7"
    grid_c = "#3A3A3A"
    text_c = "#E6E6E6"
    grey_c = "#8A8A8A"
    suspct = "#38BDF8"
    infect = "#EF4444"
    recovd = "#22C55E"
    
    fig, axes = plt.subplots(
        nrows=3, ncols=1,
        figsize=(12, 9),
        sharex=True,
        gridspec_kw={'height_ratios': [3, 1.5, 1.5]}
    )
    
    fig.canvas.manager.set_window_title(title="Моделирование заражения (SIRVD)")
    
    fig.set_facecolor(backgr)
    fig.tight_layout(pad=4)
    fig.suptitle("SIRVD", color=text_c, fontsize=16, y=0.99, fontweight='bold')
    
    def stylize(ax):
        ax.set_facecolor(backgr)
        
        ax.grid(True, which='both', linestyle='--', linewidth=0.6, alpha=0.5, color=grid_c)
        
        ax.tick_params(colors=text_c)
        ax.yaxis.label.set_color(text_c)
        ax.xaxis.label.set_color(text_c)
        ax.title.set_color(text_c)
        
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.spines["bottom"].set_color(grid_c)
        ax.spines["left"].set_color(grid_c)
        
        ax.ticklabel_format(style='plain', axis='y')
        
        legend = ax.legend(frameon=False)
        for text in legend.get_texts():
            text.set_color(text_c)
    
    # Основной график: S, I, R
    ax = axes[0]
    ax.plot(df["t"], df["S"], label="S — восприимчивые", color=suspct, linewidth=2.0, alpha=0.95)
    ax.plot(df["t"], df["I"], label="I — зараженные",    color=infect, linewidth=2.0, alpha=0.95)
    ax.plot(df["t"], df["R"], label="R — выздоровевшие", color=recovd, linewidth=2.0, alpha=0.95)
    
    ax.set_ylabel("Люди", color=text_c)
    ax.set_title("Динамика распространения", color=text_c, fontsize=12)
    
    stylize(ax)
    
    
    # График суммы населения
    ax = axes[1]
    df["sum"] = df["S"] + df["I"] + df["R"] + df["V"] + df["D"]
    
    ax.plot(df["t"], df["sum"], label="Население (S + I + R + V + D)", color=accent, linewidth=2.0, alpha=0.95)
    
    ax.set_ylabel("Люди", color=text_c)
    ax.set_title("Проверка сохранения населения", color=text_c, fontsize=12)
    
    stylize(ax)
    
    
    # График изменения шага (h)
    ax = axes[2]
    
    ax.plot(df["t"], df["h"], label="Шаг h", color=accent, linewidth=2.0, alpha=0.95)
    ax.plot(df["t"], df["h_max"], linestyle=':', color=grey_c, linewidth=1.2, alpha=0.9)
    ax.plot(df["t"], df["h_min"], linestyle=':', color=grey_c, linewidth=1.2, alpha=0.9)
    
    ax.set_ylabel("Дни", color=text_c)
    ax.set_title("Изменение адаптивного шага")
    
    stylize(ax)
    
    ax.set_xlabel("Время, дни", color=text_c)
    
    plt.show()


if __name__ == "__main__":
    out_name = "output.csv"
    
    if len(sys.argv) >= 2:
        out_name = sys.argv[1]
    
    draw_from_csv(filename=out_name)