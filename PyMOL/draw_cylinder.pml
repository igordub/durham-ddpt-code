cmd.delete("ENM")

# Cylinder
x1,y1,z1 = 10, 0, 0 # start point
r1,g1,b1 = 0,0,0 # color (black)
x2,y2,z2 = 0.1, 0, 0 # end point

spring_strength = 1
# Spring streght scales as square of srping's radius r**2 
radius = spring_strength**(0.5) * 0.25
cmd.load_cgo( [ 9.0, x1, y1, z1, x2, y2, z2, radius, r1, g1, b1, r1, g1, b1 ], "ENM" )

cmd.reset()
cmd.zoom("ENM")