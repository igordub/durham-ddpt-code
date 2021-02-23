# Cylinder #1 on left.
x1,y1,z1 = 10, 0, 0 # start point
r1,g1,b1 = 1,0,0 # color (red)
x2,y2,z2 = 0.1, 0, 0 # end point
r2,g2,b2 = 1,1,0 # color (yellow)
radius = 10
cmd.load_cgo( [ 9.0, x1, y1, z1, x2, y2, z2, radius, r1, g1, b1, r2, g2, b2 ], "cylinder1" )

# Cylinder #2 (the circle)
x1,y1,z1 = -0.1, 0, 0 # start point
r1,g1,b1 = 1,0,0 # color (red)
x2,y2,z2 = 0.1, 0, 0 # end point
r2,g2,b2 = 1,1,0 # color (yellow)
radius = 10
cmd.load_cgo( [ 9.0, x1, y1, z1, x2, y2, z2, radius, r1, g1, b1, r2, g2, b2 ], "cylinder2" )

# to make the cylinder transparent, use
cmd.load_cgo( [ 25.0, 0.25, 9.0, x1, y1, z1, x2, y2, z2, radius, r1, g1, b1, r2, g2, b2 ], "cylinder1" )