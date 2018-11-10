class ooi:
  def __init__(self, snaps):
    self.objects = dict()
    self.snap_nums = snaps

  def add(self, obj_name, obj_list):
    if obj_name in self.objects.keys():
      print obj_name, "already exists. Choose another key."
      return

    self.objects[obj_name] = obj_list

  def snap_nums(self):
      return self.snap_nums

  def get(self, obj_name):
    return self.objects[obj_name]

evo = ooi([10,
           11,
           12,
           13,
           14,
           15,
           16,
           17])
evo.add("f1",
        [ [86],
          [102],
          [89],
          [52],
          [38],
          [41] ])

evo.add("f2",
        [ [64],
          [72,94],
          [72,80],
          [38,48],
          [34],
          [33],
          [28],
          [38] ])

evo.add("f3",
        [ [138],
          [116],
          [97,126],
          [55],
          [52],
          [48],
          [45],
          [51,61] ])

evo.add("f4",
        [ [67],
          [71],
          [69],
          [39],
          [35,41],
          [35],
          [27,35],
          [39,55] ])

evo.add("f5",
        [ [91],
          [115] ])

z = range(10,28+1)
z.reverse()
devo = ooi(z)
devo.add("f1",
         [ [31],
           [27],
           [31],
           [32],
           [31],
           [28],
           [22],
           [26],
           [18],
           [17] ])

devo.add("f2",
         [ [78],
           [79],
           [82],
           [84],
           [81],
           [78],
           [56],
           [51,81],
           [47,75],
           [51],
           [35],
           [34],
           [24] ])

devo.add("f3",
         [ [47],
           [45],
           [53],
           [49],
           [44],
           [41],
           [39],
           [36],
           [33],
           [39] ])

devo.add("f4",
         [ [142],
           [144],
           [144],
           [162],
           [160],
           [122],
           [134,126],
           [133,104],
           [123,140,100,115,104,108,109,114],
           [116,124,98,104,98,106],
           [108, 81, 91, 92, 94],
           [91, 82],
           [73, 67],
           [80, 75],
           [98, 91, 88],
           [94],
           [144],
           [155],
           [174, 159, 127] ])

devo.add("f5",
         [ [99],
           [103],
           [107, 117],
           [120] ])
