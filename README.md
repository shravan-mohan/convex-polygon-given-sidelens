# Generate a convex polygon with given side lengths. 
This code generate a convex polygon with given side lengths, using a convex optimization routine. There are added 
functionalities to check if the determined coordinates adhere to the given side lengths, the polygon is simple and it 
is convex. 

# Use
    A, r = getConvexPolygon(r)
    print(checkDists(A, r))
    print(checkSimple(A, r))
    print(checkConvexity(A, r))
    
# Packages
1. Numpy
2. CVXPY
3. SymPy Geometry
4. MatplotLib
