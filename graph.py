from typing import NoReturn
from pathlib import Path
from matplotlib import pyplot as plt


def load(filename: str):
    xs, ys, zs = list(), list(), list()
    with Path(filename).open() as f:
        for datum in f.readlines():
            buf = datum.split(',')

            xs.append(int(buf[0]))
            ys.append(int(buf[1]))
            zs.append(float(buf[2]) if len(buf) == 3 else -1)
    return xs, ys, zs


def plot(title, xlabel, ylabel, xs, ys):
    fig, ax = plt.subplots()
    ax.plot(xs, ys)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.savefig(title + ".png")


def main() -> NoReturn:
    xs, ys, zs = load("det.dat")
    plot("Det_NT", "Matrix size", "Elapsed time[ms]", xs, ys)

    xs, ys, zs = load("lin.dat")
    plot("Lin_NT", "Matrix size", "Elapsed time[ms]", xs, ys)
    plot("Lin_NE", "Matrix size", "Error", xs, zs)

    xs, ys, zs = load("eig.dat")
    plot("Eig_NT", "Matrix size", "Elapsed time[ms]", xs, ys)
    plot("Eig_NE", "Matrix size", "Error", xs, zs)


if __name__ == '__main__':
    main()
