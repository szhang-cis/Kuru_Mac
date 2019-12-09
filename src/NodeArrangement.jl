function NodeArrangementHex(C)

   linear_bases_idx = [1,(C+1),(C+1)^2,(C+1)^2-C]
   quad_aranger = vcat(linear_bases_idx,deleteat!(collect(1:(C+1)^2),sort(linear_bases_idx)))
   element_numbering = copy(quad_aranger)
   for i=1:C
      faces_z = quad_aranger.+maximum(element_numbering)
      element_numbering = vcat(element_numbering,faces_z)
      if i==C
            element_numbering = element_numbering[!(element_numbering.==faces_z[1:4])]
            element_numbering = vcat(element_numbering[1:4],faces_z[1:4],element_numbering[4:end])
      end
   end
   traversed_edge_numbering_hex = nothing


   # GET FACE NUMBERING ORDER FROM TETRAHEDRAL ELEMENT
   face_0,face_1,face_2,face_3 = [],[],[],[]
   if C==0
      face_0 = [0,1,2,3]  # constant Z =-1 plane
      face_1 = [4,5,6,7]  # constant Z = 1 plane
      face_2 = [0,1,5,4]  # constant Y =-1 plane
      face_3 = [3,2,6,7]  # constant Y = 1 plane
      face_4 = [0,3,7,4]  # constant X =-1 plane
      face_5 = [1,2,6,5]  # constant X = 1 plane

      face_numbering = np.array([face_0,face_1,face_2,face_3,face_4,face_5])

   elseif C==1
      face_numbering = np.array([[ 0,  1,  2,  3,  8,  9, 10, 11, 12],
       [ 4,  5,  6,  7, 22, 23, 24, 25, 26],
       [ 0,  1,  5,  4,  8, 13, 17, 14, 22],
       [ 3,  2,  6,  7, 12, 16, 21, 15, 26],
       [ 0,  3,  7,  4,  9, 13, 18, 16, 23],
       [ 1,  2,  6,  5, 11, 14, 20, 15, 25]])

   elseif C==2
      face_numbering = np.array([[ 0,  1,  2,  3,  8,  9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19],
           [ 4,  5,  6,  7, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63],
           [ 0,  1,  5,  4,  8,  9, 20, 24, 25, 21, 36, 40, 41, 37, 52, 53],
           [ 3,  2,  6,  7, 18, 19, 23, 34, 35, 22, 39, 50, 51, 38, 62, 63],
           [ 0,  3,  7,  4, 10, 14, 20, 26, 30, 23, 36, 42, 46, 39, 54, 58],
           [ 1,  2,  6,  5, 13, 17, 21, 29, 33, 22, 37, 45, 49, 38, 57, 61]])

   elseif C==3
      face_numbering = np.array([[  0,   1,   2,   3,   8,   9,  10,  11,  12,  13,  14,  15,  16,
         17,  18,  19,  20,  21,  22,  23,  24,  25,  26,  27,  28],
       [  4,   5,   6,   7, 104, 105, 106, 107, 108, 109, 110, 111, 112,
        113, 114, 115, 116, 117, 118, 119, 120, 121, 122, 123, 124],
       [  0,   1,   5,   4,   8,   9,  10,  29,  33,  34,  35,  30,  54,
         58,  59,  60,  55,  79,  83,  84,  85,  80, 104, 105, 106],
       [  3,   2,   6,   7,  26,  27,  28,  32,  51,  52,  53,  31,  57,
         76,  77,  78,  56,  82, 101, 102, 103,  81, 122, 123, 124],
       [  0,   3,   7,   4,  11,  16,  21,  29,  36,  41,  46,  32,  54,
         61,  66,  71,  57,  79,  86,  91,  96,  82, 107, 112, 117],
       [  1,   2,   6,   5,  15,  20,  25,  30,  40,  45,  50,  31,  55,
         65,  70,  75,  56,  80,  90,  95, 100,  81, 111, 116, 121]])

   elseif C==4
      face_numbering = np.array([[  0,   1,   2,   3,   8,   9,  10,  11,  12,  13,  14,  15,  16,
         17,  18,  19,  20,  21,  22,  23,  24,  25,  26,  27,  28,  29,
         30,  31,  32,  33,  34,  35,  36,  37,  38,  39],
       [  4,   5,   6,   7, 184, 185, 186, 187, 188, 189, 190, 191, 192,
        193, 194, 195, 196, 197, 198, 199, 200, 201, 202, 203, 204, 205,
        206, 207, 208, 209, 210, 211, 212, 213, 214, 215],
       [  0,   1,   5,   4,   8,   9,  10,  11,  40,  44,  45,  46,  47,
         41,  76,  80,  81,  82,  83,  77, 112, 116, 117, 118, 119, 113,
        148, 152, 153, 154, 155, 149, 184, 185, 186, 187],
       [  3,   2,   6,   7,  36,  37,  38,  39,  43,  72,  73,  74,  75,
         42,  79, 108, 109, 110, 111,  78, 115, 144, 145, 146, 147, 114,
        151, 180, 181, 182, 183, 150, 212, 213, 214, 215],
       [  0,   3,   7,   4,  12,  18,  24,  30,  40,  48,  54,  60,  66,
         43,  76,  84,  90,  96, 102,  79, 112, 120, 126, 132, 138, 115,
        148, 156, 162, 168, 174, 151, 188, 194, 200, 206],
       [  1,   2,   6,   5,  17,  23,  29,  35,  41,  53,  59,  65,  71,
         42,  77,  89,  95, 101, 107,  78, 113, 125, 131, 137, 143, 114,
        149, 161, 167, 173, 179, 150, 193, 199, 205, 211]])

   elseif C==5
      face_numbering = np.array([[  0,   1,   2,   3,   8,   9,  10,  11,  12,  13,  14,  15,  16,
         17,  18,  19,  20,  21,  22,  23,  24,  25,  26,  27,  28,  29,
         30,  31,  32,  33,  34,  35,  36,  37,  38,  39,  40,  41,  42,
         43,  44,  45,  46,  47,  48,  49,  50,  51,  52],
       [  4,   5,   6,   7, 298, 299, 300, 301, 302, 303, 304, 305, 306,
        307, 308, 309, 310, 311, 312, 313, 314, 315, 316, 317, 318, 319,
        320, 321, 322, 323, 324, 325, 326, 327, 328, 329, 330, 331, 332,
        333, 334, 335, 336, 337, 338, 339, 340, 341, 342],
       [  0,   1,   5,   4,   8,   9,  10,  11,  12,  53,  57,  58,  59,
         60,  61,  54, 102, 106, 107, 108, 109, 110, 103, 151, 155, 156,
        157, 158, 159, 152, 200, 204, 205, 206, 207, 208, 201, 249, 253,
        254, 255, 256, 257, 250, 298, 299, 300, 301, 302],
       [  3,   2,   6,   7,  48,  49,  50,  51,  52,  56,  97,  98,  99,
        100, 101,  55, 105, 146, 147, 148, 149, 150, 104, 154, 195, 196,
        197, 198, 199, 153, 203, 244, 245, 246, 247, 248, 202, 252, 293,
        294, 295, 296, 297, 251, 338, 339, 340, 341, 342],
       [  0,   3,   7,   4,  13,  20,  27,  34,  41,  53,  62,  69,  76,
         83,  90,  56, 102, 111, 118, 125, 132, 139, 105, 151, 160, 167,
        174, 181, 188, 154, 200, 209, 216, 223, 230, 237, 203, 249, 258,
        265, 272, 279, 286, 252, 303, 310, 317, 324, 331],
       [  1,   2,   6,   5,  19,  26,  33,  40,  47,  54,  68,  75,  82,
         89,  96,  55, 103, 117, 124, 131, 138, 145, 104, 152, 166, 173,
        180, 187, 194, 153, 201, 215, 222, 229, 236, 243, 202, 250, 264,
        271, 278, 285, 292, 251, 309, 316, 323, 330, 337]])

   elseif C==6
      face_numbering =  np.array([[  0,   1,   2,   3,   8,   9,  10,  11,  12,  13,  14,  15,  16,
         17,  18,  19,  20,  21,  22,  23,  24,  25,  26,  27,  28,  29,
         30,  31,  32,  33,  34,  35,  36,  37,  38,  39,  40,  41,  42,
         43,  44,  45,  46,  47,  48,  49,  50,  51,  52,  53,  54,  55,
         56,  57,  58,  59,  60,  61,  62,  63,  64,  65,  66,  67],
       [  4,   5,   6,   7, 452, 453, 454, 455, 456, 457, 458, 459, 460,
        461, 462, 463, 464, 465, 466, 467, 468, 469, 470, 471, 472, 473,
        474, 475, 476, 477, 478, 479, 480, 481, 482, 483, 484, 485, 486,
        487, 488, 489, 490, 491, 492, 493, 494, 495, 496, 497, 498, 499,
        500, 501, 502, 503, 504, 505, 506, 507, 508, 509, 510, 511],
       [  0,   1,   5,   4,   8,   9,  10,  11,  12,  13,  68,  72,  73,
         74,  75,  76,  77,  69, 132, 136, 137, 138, 139, 140, 141, 133,
        196, 200, 201, 202, 203, 204, 205, 197, 260, 264, 265, 266, 267,
        268, 269, 261, 324, 328, 329, 330, 331, 332, 333, 325, 388, 392,
        393, 394, 395, 396, 397, 389, 452, 453, 454, 455, 456, 457],
       [  3,   2,   6,   7,  62,  63,  64,  65,  66,  67,  71, 126, 127,
        128, 129, 130, 131,  70, 135, 190, 191, 192, 193, 194, 195, 134,
        199, 254, 255, 256, 257, 258, 259, 198, 263, 318, 319, 320, 321,
        322, 323, 262, 327, 382, 383, 384, 385, 386, 387, 326, 391, 446,
        447, 448, 449, 450, 451, 390, 506, 507, 508, 509, 510, 511],
       [  0,   3,   7,   4,  14,  22,  30,  38,  46,  54,  68,  78,  86,
         94, 102, 110, 118,  71, 132, 142, 150, 158, 166, 174, 182, 135,
        196, 206, 214, 222, 230, 238, 246, 199, 260, 270, 278, 286, 294,
        302, 310, 263, 324, 334, 342, 350, 358, 366, 374, 327, 388, 398,
        406, 414, 422, 430, 438, 391, 458, 466, 474, 482, 490, 498],
       [  1,   2,   6,   5,  21,  29,  37,  45,  53,  61,  69,  85,  93,
        101, 109, 117, 125,  70, 133, 149, 157, 165, 173, 181, 189, 134,
        197, 213, 221, 229, 237, 245, 253, 198, 261, 277, 285, 293, 301,
        309, 317, 262, 325, 341, 349, 357, 365, 373, 381, 326, 389, 405,
        413, 421, 429, 437, 445, 390, 465, 473, 481, 489, 497, 505]])

   elseif C==8
      face_numbering = np.array([[  0,   1,   2,   3,   8,   9,  10,  11,  12,  13,  14,  15,  16,
         17,  18,  19,  20,  21,  22,  23,  24,  25,  26,  27,  28,  29,
         30,  31,  32,  33,  34,  35,  36,  37,  38,  39,  40,  41,  42,
         43,  44,  45,  46,  47,  48,  49,  50,  51,  52,  53,  54,  55,
         56,  57,  58,  59,  60,  61,  62,  63,  64,  65,  66,  67,  68,
         69,  70,  71,  72,  73,  74,  75,  76,  77,  78,  79,  80,  81,
         82,  83,  84],
       [  4,   5,   6,   7, 652, 653, 654, 655, 656, 657, 658, 659, 660,
        661, 662, 663, 664, 665, 666, 667, 668, 669, 670, 671, 672, 673,
        674, 675, 676, 677, 678, 679, 680, 681, 682, 683, 684, 685, 686,
        687, 688, 689, 690, 691, 692, 693, 694, 695, 696, 697, 698, 699,
        700, 701, 702, 703, 704, 705, 706, 707, 708, 709, 710, 711, 712,
        713, 714, 715, 716, 717, 718, 719, 720, 721, 722, 723, 724, 725,
        726, 727, 728],
       [  0,   1,   5,   4,   8,   9,  10,  11,  12,  13,  14,  85,  89,
         90,  91,  92,  93,  94,  95,  86, 166, 170, 171, 172, 173, 174,
        175, 176, 167, 247, 251, 252, 253, 254, 255, 256, 257, 248, 328,
        332, 333, 334, 335, 336, 337, 338, 329, 409, 413, 414, 415, 416,
        417, 418, 419, 410, 490, 494, 495, 496, 497, 498, 499, 500, 491,
        571, 575, 576, 577, 578, 579, 580, 581, 572, 652, 653, 654, 655,
        656, 657, 658],
       [  3,   2,   6,   7,  78,  79,  80,  81,  82,  83,  84,  88, 159,
        160, 161, 162, 163, 164, 165,  87, 169, 240, 241, 242, 243, 244,
        245, 246, 168, 250, 321, 322, 323, 324, 325, 326, 327, 249, 331,
        402, 403, 404, 405, 406, 407, 408, 330, 412, 483, 484, 485, 486,
        487, 488, 489, 411, 493, 564, 565, 566, 567, 568, 569, 570, 492,
        574, 645, 646, 647, 648, 649, 650, 651, 573, 722, 723, 724, 725,
        726, 727, 728],
       [  0,   3,   7,   4,  15,  24,  33,  42,  51,  60,  69,  85,  96,
        105, 114, 123, 132, 141, 150,  88, 166, 177, 186, 195, 204, 213,
        222, 231, 169, 247, 258, 267, 276, 285, 294, 303, 312, 250, 328,
        339, 348, 357, 366, 375, 384, 393, 331, 409, 420, 429, 438, 447,
        456, 465, 474, 412, 490, 501, 510, 519, 528, 537, 546, 555, 493,
        571, 582, 591, 600, 609, 618, 627, 636, 574, 659, 668, 677, 686,
        695, 704, 713],
       [  1,   2,   6,   5,  23,  32,  41,  50,  59,  68,  77,  86, 104,
        113, 122, 131, 140, 149, 158,  87, 167, 185, 194, 203, 212, 221,
        230, 239, 168, 248, 266, 275, 284, 293, 302, 311, 320, 249, 329,
        347, 356, 365, 374, 383, 392, 401, 330, 410, 428, 437, 446, 455,
        464, 473, 482, 411, 491, 509, 518, 527, 536, 545, 554, 563, 492,
        572, 590, 599, 608, 617, 626, 635, 644, 573, 667, 676, 685, 694,
        703, 712, 721]])

   elseif C==9
      face_numbering = np.array([[  0,   1,   2,   3,   8,   9,  10,  11,  12,  13,  14,  15,  16,
         17,  18,  19,  20,  21,  22,  23,  24,  25,  26,  27,  28,  29,
         30,  31,  32,  33,  34,  35,  36,  37,  38,  39,  40,  41,  42,
         43,  44,  45,  46,  47,  48,  49,  50,  51,  52,  53,  54,  55,
         56,  57,  58,  59,  60,  61,  62,  63,  64,  65,  66,  67,  68,
         69,  70,  71,  72,  73,  74,  75,  76,  77,  78,  79,  80,  81,
         82,  83,  84,  85,  86,  87,  88,  89,  90,  91,  92,  93,  94,
         95,  96,  97,  98,  99, 100, 101, 102, 103],
       [  4,   5,   6,   7, 904, 905, 906, 907, 908, 909, 910, 911, 912,
        913, 914, 915, 916, 917, 918, 919, 920, 921, 922, 923, 924, 925,
        926, 927, 928, 929, 930, 931, 932, 933, 934, 935, 936, 937, 938,
        939, 940, 941, 942, 943, 944, 945, 946, 947, 948, 949, 950, 951,
        952, 953, 954, 955, 956, 957, 958, 959, 960, 961, 962, 963, 964,
        965, 966, 967, 968, 969, 970, 971, 972, 973, 974, 975, 976, 977,
        978, 979, 980, 981, 982, 983, 984, 985, 986, 987, 988, 989, 990,
        991, 992, 993, 994, 995, 996, 997, 998, 999],
       [  0,   1,   5,   4,   8,   9,  10,  11,  12,  13,  14,  15, 104,
        108, 109, 110, 111, 112, 113, 114, 115, 105, 204, 208, 209, 210,
        211, 212, 213, 214, 215, 205, 304, 308, 309, 310, 311, 312, 313,
        314, 315, 305, 404, 408, 409, 410, 411, 412, 413, 414, 415, 405,
        504, 508, 509, 510, 511, 512, 513, 514, 515, 505, 604, 608, 609,
        610, 611, 612, 613, 614, 615, 605, 704, 708, 709, 710, 711, 712,
        713, 714, 715, 705, 804, 808, 809, 810, 811, 812, 813, 814, 815,
        805, 904, 905, 906, 907, 908, 909, 910, 911],
       [  3,   2,   6,   7,  96,  97,  98,  99, 100, 101, 102, 103, 107,
        196, 197, 198, 199, 200, 201, 202, 203, 106, 207, 296, 297, 298,
        299, 300, 301, 302, 303, 206, 307, 396, 397, 398, 399, 400, 401,
        402, 403, 306, 407, 496, 497, 498, 499, 500, 501, 502, 503, 406,
        507, 596, 597, 598, 599, 600, 601, 602, 603, 506, 607, 696, 697,
        698, 699, 700, 701, 702, 703, 606, 707, 796, 797, 798, 799, 800,
        801, 802, 803, 706, 807, 896, 897, 898, 899, 900, 901, 902, 903,
        806, 992, 993, 994, 995, 996, 997, 998, 999],
       [  0,   3,   7,   4,  16,  26,  36,  46,  56,  66,  76,  86, 104,
        116, 126, 136, 146, 156, 166, 176, 186, 107, 204, 216, 226, 236,
        246, 256, 266, 276, 286, 207, 304, 316, 326, 336, 346, 356, 366,
        376, 386, 307, 404, 416, 426, 436, 446, 456, 466, 476, 486, 407,
        504, 516, 526, 536, 546, 556, 566, 576, 586, 507, 604, 616, 626,
        636, 646, 656, 666, 676, 686, 607, 704, 716, 726, 736, 746, 756,
        766, 776, 786, 707, 804, 816, 826, 836, 846, 856, 866, 876, 886,
        807, 912, 922, 932, 942, 952, 962, 972, 982],
       [  1,   2,   6,   5,  25,  35,  45,  55,  65,  75,  85,  95, 105,
        125, 135, 145, 155, 165, 175, 185, 195, 106, 205, 225, 235, 245,
        255, 265, 275, 285, 295, 206, 305, 325, 335, 345, 355, 365, 375,
        385, 395, 306, 405, 425, 435, 445, 455, 465, 475, 485, 495, 406,
        505, 525, 535, 545, 555, 565, 575, 585, 595, 506, 605, 625, 635,
        645, 655, 665, 675, 685, 695, 606, 705, 725, 735, 745, 755, 765,
        775, 785, 795, 706, 805, 825, 835, 845, 855, 865, 875, 885, 895,
        806, 921, 931, 941, 951, 961, 971, 981, 991]])

   elseif C==10
      face_numbering = np.array([[   0,    1,    2,    3,    8,    9,   10,   11,   12,   13,   14,
          15,   16,   17,   18,   19,   20,   21,   22,   23,   24,   25,
          26,   27,   28,   29,   30,   31,   32,   33,   34,   35,   36,
          37,   38,   39,   40,   41,   42,   43,   44,   45,   46,   47,
          48,   49,   50,   51,   52,   53,   54,   55,   56,   57,   58,
          59,   60,   61,   62,   63,   64,   65,   66,   67,   68,   69,
          70,   71,   72,   73,   74,   75,   76,   77,   78,   79,   80,
          81,   82,   83,   84,   85,   86,   87,   88,   89,   90,   91,
          92,   93,   94,   95,   96,   97,   98,   99,  100,  101,  102,
         103,  104,  105,  106,  107,  108,  109,  110,  111,  112,  113,
         114,  115,  116,  117,  118,  119,  120,  121,  122,  123,  124,
         125,  126,  127,  128,  129,  130,  131,  132,  133,  134,  135,
         136,  137,  138,  139,  140,  141,  142,  143,  144,  145,  146,
         147],
       [   4,    5,    6,    7, 1588, 1589, 1590, 1591, 1592, 1593, 1594,
        1595, 1596, 1597, 1598, 1599, 1600, 1601, 1602, 1603, 1604, 1605,
        1606, 1607, 1608, 1609, 1610, 1611, 1612, 1613, 1614, 1615, 1616,
        1617, 1618, 1619, 1620, 1621, 1622, 1623, 1624, 1625, 1626, 1627,
        1628, 1629, 1630, 1631, 1632, 1633, 1634, 1635, 1636, 1637, 1638,
        1639, 1640, 1641, 1642, 1643, 1644, 1645, 1646, 1647, 1648, 1649,
        1650, 1651, 1652, 1653, 1654, 1655, 1656, 1657, 1658, 1659, 1660,
        1661, 1662, 1663, 1664, 1665, 1666, 1667, 1668, 1669, 1670, 1671,
        1672, 1673, 1674, 1675, 1676, 1677, 1678, 1679, 1680, 1681, 1682,
        1683, 1684, 1685, 1686, 1687, 1688, 1689, 1690, 1691, 1692, 1693,
        1694, 1695, 1696, 1697, 1698, 1699, 1700, 1701, 1702, 1703, 1704,
        1705, 1706, 1707, 1708, 1709, 1710, 1711, 1712, 1713, 1714, 1715,
        1716, 1717, 1718, 1719, 1720, 1721, 1722, 1723, 1724, 1725, 1726,
        1727],
       [   0,    1,    5,    4,    8,    9,   10,   11,   12,   13,   14,
          15,   16,   17,  148,  152,  153,  154,  155,  156,  157,  158,
         159,  160,  161,  149,  292,  296,  297,  298,  299,  300,  301,
         302,  303,  304,  305,  293,  436,  440,  441,  442,  443,  444,
         445,  446,  447,  448,  449,  437,  580,  584,  585,  586,  587,
         588,  589,  590,  591,  592,  593,  581,  724,  728,  729,  730,
         731,  732,  733,  734,  735,  736,  737,  725,  868,  872,  873,
         874,  875,  876,  877,  878,  879,  880,  881,  869, 1012, 1016,
        1017, 1018, 1019, 1020, 1021, 1022, 1023, 1024, 1025, 1013, 1156,
        1160, 1161, 1162, 1163, 1164, 1165, 1166, 1167, 1168, 1169, 1157,
        1300, 1304, 1305, 1306, 1307, 1308, 1309, 1310, 1311, 1312, 1313,
        1301, 1444, 1448, 1449, 1450, 1451, 1452, 1453, 1454, 1455, 1456,
        1457, 1445, 1588, 1589, 1590, 1591, 1592, 1593, 1594, 1595, 1596,
        1597],
       [   3,    2,    6,    7,  138,  139,  140,  141,  142,  143,  144,
         145,  146,  147,  151,  282,  283,  284,  285,  286,  287,  288,
         289,  290,  291,  150,  295,  426,  427,  428,  429,  430,  431,
         432,  433,  434,  435,  294,  439,  570,  571,  572,  573,  574,
         575,  576,  577,  578,  579,  438,  583,  714,  715,  716,  717,
         718,  719,  720,  721,  722,  723,  582,  727,  858,  859,  860,
         861,  862,  863,  864,  865,  866,  867,  726,  871, 1002, 1003,
        1004, 1005, 1006, 1007, 1008, 1009, 1010, 1011,  870, 1015, 1146,
        1147, 1148, 1149, 1150, 1151, 1152, 1153, 1154, 1155, 1014, 1159,
        1290, 1291, 1292, 1293, 1294, 1295, 1296, 1297, 1298, 1299, 1158,
        1303, 1434, 1435, 1436, 1437, 1438, 1439, 1440, 1441, 1442, 1443,
        1302, 1447, 1578, 1579, 1580, 1581, 1582, 1583, 1584, 1585, 1586,
        1587, 1446, 1718, 1719, 1720, 1721, 1722, 1723, 1724, 1725, 1726,
        1727],
       [   0,    3,    7,    4,   18,   30,   42,   54,   66,   78,   90,
         102,  114,  126,  148,  162,  174,  186,  198,  210,  222,  234,
         246,  258,  270,  151,  292,  306,  318,  330,  342,  354,  366,
         378,  390,  402,  414,  295,  436,  450,  462,  474,  486,  498,
         510,  522,  534,  546,  558,  439,  580,  594,  606,  618,  630,
         642,  654,  666,  678,  690,  702,  583,  724,  738,  750,  762,
         774,  786,  798,  810,  822,  834,  846,  727,  868,  882,  894,
         906,  918,  930,  942,  954,  966,  978,  990,  871, 1012, 1026,
        1038, 1050, 1062, 1074, 1086, 1098, 1110, 1122, 1134, 1015, 1156,
        1170, 1182, 1194, 1206, 1218, 1230, 1242, 1254, 1266, 1278, 1159,
        1300, 1314, 1326, 1338, 1350, 1362, 1374, 1386, 1398, 1410, 1422,
        1303, 1444, 1458, 1470, 1482, 1494, 1506, 1518, 1530, 1542, 1554,
        1566, 1447, 1598, 1610, 1622, 1634, 1646, 1658, 1670, 1682, 1694,
        1706],
       [   1,    2,    6,    5,   29,   41,   53,   65,   77,   89,  101,
         113,  125,  137,  149,  173,  185,  197,  209,  221,  233,  245,
         257,  269,  281,  150,  293,  317,  329,  341,  353,  365,  377,
         389,  401,  413,  425,  294,  437,  461,  473,  485,  497,  509,
         521,  533,  545,  557,  569,  438,  581,  605,  617,  629,  641,
         653,  665,  677,  689,  701,  713,  582,  725,  749,  761,  773,
         785,  797,  809,  821,  833,  845,  857,  726,  869,  893,  905,
         917,  929,  941,  953,  965,  977,  989, 1001,  870, 1013, 1037,
        1049, 1061, 1073, 1085, 1097, 1109, 1121, 1133, 1145, 1014, 1157,
        1181, 1193, 1205, 1217, 1229, 1241, 1253, 1265, 1277, 1289, 1158,
        1301, 1325, 1337, 1349, 1361, 1373, 1385, 1397, 1409, 1421, 1433,
        1302, 1445, 1469, 1481, 1493, 1505, 1517, 1529, 1541, 1553, 1565,
        1577, 1446, 1609, 1621, 1633, 1645, 1657, 1669, 1681, 1693, 1705,
        1717]])

   else
        # THIS IS A FLOATING POINT BASED ALGORITHM
        from Florence.QuadratureRules.NumericIntegrator import GaussLobattoQuadrature
        xs = GaussLobattoQuadrature(C+2)[0]
        x,y,z = np.meshgrid(xs,xs,xs)
        fekete = np.concatenate((y.T.flatten()[:,None],x.T.flatten()[:,None],z.T.flatten()[:,None]),axis=1)
        fekete = fekete[element_numbering,:]
        tol=1e-12
        nsize = int((C+2)**3)

        x,y = np.meshgrid(xs,xs)
        points = np.concatenate((x.flatten()[:,None],y.flatten()[:,None]),axis=1)
        facer = points[quad_aranger]

        all_nodes = np.arange(nsize)

        # CONSTANT Z PLANES
        facer3d = np.zeros((facer.shape[0],3))
        facer3d[:,:2] = facer

        from Florence.Tensor import in2d_unsorted

        # facer3d[:,2] = xs[0]
        # face_0 = __FaceArrangementHex__(fekete,facer3d)
        # facer3d[:,2] = xs[-1]
        # face_1 = __FaceArrangementHex__(fekete,facer3d)

        # # CONSTANT Y PLANES
        # facer3d = np.zeros((facer.shape[0],3))
        # facer3d[:,[0,2]] = facer

        # facer3d[:,1] = xs[0]
        # face_2 = __FaceArrangementHex__(fekete,facer3d)
        # facer3d[:,1] = xs[-1]
        # face_3 = __FaceArrangementHex__(fekete,facer3d)

        # # CONSTANT X PLANES
        # facer3d = np.zeros((facer.shape[0],3))
        # facer3d[:,[1,2]] = facer

        # facer3d[:,0] = xs[0]
        # face_4 = __FaceArrangementHex__(fekete,facer3d)
        # facer3d[:,0] = xs[-1]
        # face_5 = __FaceArrangementHex__(fekete,facer3d)


        facer3d[:,2] = xs[0]
        face_0 = in2d_unsorted(fekete,facer3d)
        facer3d[:,2] = xs[-1]
        face_1 = in2d_unsorted(fekete,facer3d)

        # CONSTANT Y PLANES
        facer3d = np.zeros((facer.shape[0],3))
        facer3d[:,[0,2]] = facer

        facer3d[:,1] = xs[0]
        face_2 = in2d_unsorted(fekete,facer3d)
        facer3d[:,1] = xs[-1]
        face_3 = in2d_unsorted(fekete,facer3d)

        # CONSTANT X PLANES
        facer3d = np.zeros((facer.shape[0],3))
        facer3d[:,[1,2]] = facer

        facer3d[:,0] = xs[0]
        face_4 = in2d_unsorted(fekete,facer3d)
        facer3d[:,0] = xs[-1]
        face_5 = in2d_unsorted(fekete,facer3d)

        # SANITY CHECK
        assert len(face_0) == len(face_1)
        assert len(face_1) == len(face_2)
        assert len(face_2) == len(face_3)
        assert len(face_3) == len(face_4)
        assert len(face_4) == len(face_5)

        face_numbering = np.array([face_0,face_1,face_2,face_3,face_4,face_5])
        # print("face_numbering = ", repr(face_numbering))

   return face_numbering, traversed_edge_numbering_hex, element_numbering
end