import numpy as np

face_0 = [0,1,2,3]  # constant Z =-1 plane
face_1 = [4,5,6,7]  # constant Z = 1 plane
face_2 = [0,1,5,4]  # constant Y =-1 plane
face_3 = [3,2,6,7]  # constant Y = 1 plane
face_4 = [0,3,7,4]  # constant X =-1 plane
face_5 = [1,2,6,5]  # constant X = 1 plane

face_numbering = np.array([face_0,face_1,face_2,face_3,face_4,face_5])
print(face_numbering)
