from os.path import normpath as os_normpath, join as os_join, exists as os_exists
from json import load as json_load
# from threading import Thread, Lock
# from functools import wraps


def mkpath(*paths):
    return os_normpath(os_join(*paths))


def prtNum(num):
    return "{:,}".format(num).replace(',', "'")


def equalE(value1, value2, E):
    '''Epsilon comparison'''
    return value1 - E < value2 < value1 + E


def distance2(x1, y1, x2, y2):
    return (x1 - x2) ** 2 + (y1 - y2) ** 2


def linearApproxDots(dots):
    n, sumx, sumy, sumx2, sumxy = len(dots), 0, 0, 0, 0
    for x, y in dots:
        sumx += x
        sumy += y
        sumx2 += x ** 2
        sumxy += x * y

    k = (n * sumxy - (sumx * sumy)) / (n * sumx2 - sumx * sumx)
    b = (sumy - k * sumx) / n
    return k, b


def linearApproxLines(lines):
    dots = []
    for line in lines:
        dots += line.dots
    return linearApproxDots(dots)


def YCoordOnLine(x1, y1, x2, y2, target_x):
    return y1 + (y2 - y1) * ((target_x - x1) / (x2 - x1))


def setSettings(settings, alternative_settings_path=None):
    if alternative_settings_path is not None and os_exists(alternative_settings_path):
        with open(alternative_settings_path, 'r', encoding="utf-8") as settings_file:
            settings.update(json_load(settings_file))

    counted = {}

    for key, value in settings.items():
        if not isinstance(value, str) or value[0] != '$':
            if isinstance(value, float) and value % 1 == 0:
                settings[key] = int(value)
            counted[key] = value

    while True:
        for key, value in settings.items():
            if key not in counted:
                try:
                    settings[key] = eval(value[1:], None, counted)
                except NameError:
                    continue
                if isinstance(value, float) and value % 1 == 0:
                    settings[key] = int(value)
                counted[key] = value
                break
        else:
            break


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


# class Threaded(Thread):
#     def __init__(self, *args, **kwargs):
#         self.target = kwargs.pop('target')
#         self.finished = False
#         super(Threaded, self).__init__(target=self.saveTarget, *args, **kwargs)

#     def saveTarget(self):
#         self.target()
#         self.finished = True
