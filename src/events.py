# --------------------------------------------------------------------------------> Rotation
class Rotation:
    def __init__(self, start_line, end_line, rotation_center):
        self.type = "Rotation"
        self.start_line = start_line
        self.end_line = end_line
        self.rotation_center = rotation_center

    def __repr__(self):
        return "Rotation(start_line={}, end_line={}, rotation_center={})".format(
            self.start_line, self.end_line, self.rotation_center
        )


# --------------------------------------------------------------------------------> Deletion
class Deletion:
    def __init__(self, start_x, start_y, length):
        self.type = "Deletion"
        self.start_x = start_x
        self.start_y = start_y
        self.length = length

    def __repr__(self):
        return "Deletion(start_x={}, start_y={}, length={})".format(
            self.start_x, self.start_y, self.length
        )

    @property
    def size(self):
        return self.length


# --------------------------------------------------------------------------------> Insertion
class Insertion:
    def __init__(self, start_x, start_y, height):
        self.type = "Insertion"
        self.start_x = start_x
        self.start_y = start_y
        self.height = height

    def __repr__(self):
        return "Insertion(start_x={}, start_y={}, height={})".format(
            self.start_x, self.start_y, self.height
        )

    @property
    def size(self):
        return self.height


# --------------------------------------------------------------------------------> Translocation
class Translocation:
    def __init__(self, start_x, start_y, height):
        self.type = "Translocation"
        self.start_x = start_x
        self.start_y = start_y
        self.height = height

    def __repr__(self):
        return "Translocation(start_x={}, start_y={}, height={})".format(
            self.start_x, self.start_y, self.height
        )

    @property
    def size(self):
        return self.height


# --------------------------------------------------------------------------------> Duplication
class Duplication:
    def __init__(self, start_x, start_y, length, height, line_index):
        self.type = "Duplication"
        self.start_x = start_x
        self.start_y = start_y
        self.length = length
        self.height = height
        self.line_index = line_index

    def __repr__(self):
        return "Duplication(start_x={}, start_y={}, length={}, height={}, line_index={})".format(
            self.start_x, self.start_y, self.length, self.height, self.line_index
        )

    @property
    def size(self):
        return self.length


# --------------------------------------------------------------------------------> Pass (Nothing)
class Pass:
    def __init__(self):
        pass
