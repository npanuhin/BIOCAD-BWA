class Line:
    """
    Properties:
        start_x  {0}
        start_y  {1}
        end_x    {2}
        end_y    {3}
        dots = [dot1, ..., dotN] {4}
        coords = (start_x, start_y, end_x, end_y)
    """

    def __init__(self, start_x=None, start_y=None, end_x=None, end_y=None, dots=[]):
        self.start_x = start_x
        self.start_y = start_y
        self.end_x = end_x
        self.end_y = end_y
        self.dots = dots

    @property
    def coords(self):
        return self.start_x, self.start_y, self.end_x, self.end_y

    # @property
    # def x1(self):
    #     return self.start_x

    # @property
    # def y1(self):
    #     return self.start_y

    # @property
    # def x2(self):
    #     return self.end_x

    # @property
    # def y2(self):
    #     return self.end_y
