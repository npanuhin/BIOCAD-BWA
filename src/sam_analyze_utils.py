import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Polygon
# from threading import Thread, Lock
# from functools import wraps


# def threded():
#     def decorator(funtion):
#         @wraps(funtion)
#         def run(self, *args, **kwargs):
#             self.thread = Thread(target=funtion, args=[self] + list(args), kwargs=kwargs)
#             self.thread.start()
#         return run
#     return decorator


# def threadSave():
#     def decorator(funtion):
#         @wraps(funtion)
#         def run(self, *args, **kwargs):
#             if self.thread is not None:
#                 self.thread.join()

#             funtion(self, *args, **kwargs)
#         return run
#     return decorator


class Plot:
    def __init__(self, title, fontsize, grid_size=None, figsize=None, nameX=None, nameY=None):
        self.fig = plt.figure(title, figsize, tight_layout=True)
        self.ax = self.fig.add_subplot()

        self.legend = None

        # self.fig.rc("font", size=fontsize)               # controls default text sizes
        # self.fig.rc("axes", titlesize=fontsize)          # fontsize of the axes title
        # self.ax.rc("axes", labelsize=8)                  # fontsize of the x and y labels
        # self.ax.rc("xtick", labelsize=fontsize)         # fontsize of the tick labels
        # self.ax.rc("ytick", labelsize=fontsize)         # fontsize of the tick labels
        # self.ax.rc("legend", fontsize=fontsize)         # legend fontsize
        # self.ax.rc("figure", titlesize=fontsize)        # fontsize of the figure title
        self.ax.ticklabel_format(style="plain")
        if grid_size is not None:
            self.ax.set_xticks(list(range(0, int(grid_size * 1000), int(grid_size))))
            self.ax.set_yticks(list(range(0, int(grid_size * 1000), int(grid_size))))
        self.ax.tick_params(axis="x", which="major", labelsize=fontsize, rotation=30)
        self.ax.tick_params(axis="y", which="major", labelsize=fontsize)
        self.ax.grid(which="major", linestyle='-', linewidth="1", alpha=0.1, color="black")

        # self.ax.minorticks_on()
        # self.ax.grid(which="minor", linestyle=':', linewidth="1", color="black")

        if nameX is not None:
            self.ax.set_xlabel(nameX)
        if nameY is not None:
            self.ax.set_ylabel(nameY)

        # self.thread = None

    def __del__(self):
        self.close()

    def scatter(self, dots, dotsize=None, *args, **kwargs):
        self.ax.scatter(*zip(*dots), s=dotsize, *args, **kwargs)
        # self.ax.plot(*zip(*dots), color='none')  # Walkaround for relim() to work

    def line(self, x1, y1, x2, y2, *args, **kwargs):
        self.ax.plot([x1, x2], [y1, y2], *args, **kwargs)

    def legendLine(self, legend_dict, fontsize=None, *line_args, **line_kwargs):
        legend_objects, legend_names = [], []
        for name in legend_dict:
            legend_names.append(name)
            legend_objects.append(Line2D([0], [0], color=legend_dict[name], *line_args, **line_kwargs))

        self.legend = self.ax.legend(legend_objects, legend_names, fontsize=fontsize)

    def poligon(self, dots, *args, **kwargs):
        self.ax.add_patch(Polygon(dots, *args, **kwargs))

    def tight(self):
        self.ax.ignore_existing_data_limits = True
        # self.ax.update_datalim(self.scatter.get_datalim(self.ax.transData))
        # self.ax.relim()
        self.ax.autoscale_view()

    def save(self, path, *args, **kwargs):
        self.fig.savefig(path, dpi=400, *args, **kwargs)

    def clear(self):
        for artist in self.ax.lines + self.ax.collections:
            artist.remove()

        if self.legend is not None:
            self.legend.remove()

    def show(self):
        plt.show()

    def close(self):
        plt.close(self.fig)


# class Threaded(Thread):
#     def __init__(self, *args, **kwargs):
#         self.target = kwargs.pop('target')
#         self.finished = False
#         super(Threaded, self).__init__(target=self.saveTarget, *args, **kwargs)

#     def saveTarget(self):
#         self.target()
#         self.finished = True
