 face_numbering = [ 0,  1,  2,  3,  8,  9, 10, 11, 12,
        4,  5,  6,  7, 22, 23, 24, 25, 26,
        0,  1,  5,  4,  8, 13, 17, 14, 22,
        3,  2,  6,  7, 12, 16, 21, 15, 26,
        0,  3,  7,  4,  9, 13, 18, 16, 23,
        1,  2,  6,  5, 11, 14, 20, 15, 25]
      face_numbering = transpose(reshape(face_numbering.+1,9,6))
      face_numbering = copy(face_numbering)
println(typeof(face_numbering))
println(face_numbering)
