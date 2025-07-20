#!/usr/bin/env python


# Take a list of integers, and make a line plot showing dots where the values are.

import sys
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator, ScalarFormatter
from matplotlib import collections as mc


def stats(l):
    h = {}
    for x in l:
        if x not in h:
            h[x] = 0
        h[x] += 1
    h = list(h.items())
    h.sort()
    print(h)


def plot_samples(samples, d=False):
    # two figs in 1 plot
    fig, axs = plt.subplots(1, 1, figsize=(5, 10), sharex=True)

    ax = axs

    # secondary y-axis for tree depth
    # names = ["pos-", "pos+", "lex-", "lex+", "colex-", "colex+", "rand"]
    names = ["pos-", "colex-"]
    for i, (sample, name) in enumerate(zip(samples, names)):
        print(name)

        xs = [x[0] for x in sample]
        ymins = [x[1][0] for x in sample]
        ymaxs = [x[1][1] for x in sample]

        # ymax stats
        stats(ymins)
        stats(ymaxs)
        stats([x[0] - y[0] for x, y in zip(sample[1:], sample)])

        # For each x, plot a line between ymin and ymax
        lines = []
        lines2 = []
        for j, (x, ymin, ymax) in enumerate(zip(xs, ymins, ymaxs)):
            # x = x + i / len(samples) * 3
            lines.append([(x - ymin + 1, j), (x - ymax, j)])
            if ymax > 7:
                lines2.append([(x, j), (x, j + 4)])
        # Next colour of the default colour cycle
        icolor = plt.rcParams["axes.prop_cycle"].by_key()["color"][i]
        lc = mc.LineCollection(lines, alpha=0.5, color=icolor, lw=2)
        lc2 = mc.LineCollection(lines2, alpha=0.5, color="black", lw=0.5)
        if i == 0 or True:
            ax.add_collection(lc)
            ax.add_collection(lc2)

            # ax2.plot(
            #     xs,
            #     ymins,
            #     marker="o",
            #     linestyle="None",
            #     color=icolor,
            #     alpha=0.5,
            #     label=name,
            #     markersize=5,
            # )
            # ax2.plot(
            #     xs,
            #     ymaxs,
            #     marker="o",
            #     linestyle="None",
            #     color=icolor,
            #     alpha=0.5,
            #     label=name,
            #     markersize=5,
            # )
        # ax2.autoscale()

        ax.plot(
            xs,
            range(len(xs)),
            marker="o",
            linestyle="None",
            alpha=0.5,
            label=name,
            markersize=5,
        )

    ax.set_title("STPD samples")
    # ax2.set_xlabel("Text pos")
    ax.set_ylabel("Sample index")
    ax.legend()
    # gridlines every 100
    ax.xaxis.set_major_locator(MultipleLocator(100))

    # ax.set_ylim(-100, 200)
    ax.yaxis.set_major_locator(MultipleLocator(20))
    ax.yaxis.set_major_formatter(ScalarFormatter())

    # ax2.set_ylabel("Tree depth range")
    # ax2.set_ylim(0, 200)
    # ax2.yaxis.set_major_locator(MultipleLocator(20))
    # ax2.yaxis.set_major_formatter(ScalarFormatter())

    # plt.yscale("log", base=2)
    # ax.yaxis.set_major_formatter(ScalarFormatter())
    # yticks at powers of 2
    # plt.yticks([2**i for i in range(0, 5)])

    ax.grid(True, alpha=0.5)
    # ax2.grid(True)
    plt.show()


# Read lines from stdin and plot them
samples = []
for line in sys.stdin:
    samples.append(eval(line))

plot_samples(samples)
