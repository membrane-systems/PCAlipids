python

for i in range(1,2):
	cmd.reinitialize()
	cmd.set("orthoscopic","on")
	cmd.set("bg_rgb", [1,1,1])
	cmd.set("ray_shadow","false")
	cmd.viewport(800,400)

	cmd.load("average.pdb","temp")
	cmd.load("extreme"+str(i)+".pdb", "temp")

	cmd.select("NC3", "name NC3")
	cmd.select("PO4", "name PO4")
	cmd.select("GL1", "name GL1")
	cmd.select("GL2", "name GL2")
	cmd.select("C1A", "name C1A")
	cmd.select("C2A", "name C2A")
	cmd.select("C3A", "name C3A")
	cmd.select("C4A", "name C4A")
	cmd.select("C1B", "name C1B")
	cmd.select("C2B", "name C2B")
	cmd.select("C3B", "name C3B")
	cmd.select("C4B", "name C4B")
	
        cmd.unbond("all","all")
	cmd.bond("C4A","C3A")
	cmd.bond("C4B","C3B")
	cmd.bond("C3A","C2A")
	cmd.bond("C3B","C2B")
	cmd.bond("C2A","C1A")
	cmd.bond("C2B","C1B")
	cmd.bond("C1A","GL1")
	cmd.bond("C1B","GL2")
	cmd.bond("GL1","GL2")
	cmd.bond("GL1","PO4")
	cmd.bond("PO4","NC3")

	cmd.split_states("temp")
	cmd.set_name("temp_0001","average")
	cmd.delete("temp")

	cmd.set_view ("	-0.263482451,    0.010398334,    0.964607179,\
     			0.964657962,   -0.000748686,    0.263503969,\
     			0.003461469,    0.999943972,   -0.009831118,\
    			-0.000000638,   -0.000003096,  -65.826843262,\
    			-0.028863965,    0.106127433,    0.871515691,\
			42.953697205,   88.699935913,   20.000000000")

	cmd.hide("everything")
	cmd.show("sticks","all")
	cmd.set("stick_radius",.1)
	cmd.set("stick_transparency", 0.2)
	cmd.set("stick_transparency", 0, "average")

	cmd.set_color("avCol", [0.7, 0.7, 0.7])

	ffCol  = [0.05, 0.0, 1.0]
	lfCol  = [1.0, 0.0, 0.05]
	subCol = [y-x for x,y in zip(ffCol,lfCol)]
	stepCol= [x/19 for x in subCol]

	for x in range(2,10,1):
	    addCol = [y*(x-2) for y in stepCol]
	    mCol = [y+z for y,z in zip(ffCol,addCol)]
	    cmd.set_color("midCol"+str(x),mCol)
	    cmd.color("midCol"+str(x),"temp_000"+str(x))

	for x in range(10,22,1):
	    addCol = [y*(x-2) for y in stepCol]
 	    mCol = [y+z for y,z in zip(ffCol,addCol)]
 	    cmd.set_color("midCol"+str(x),mCol)
  	    cmd.color("midCol"+str(x),"temp_00"+str(x))


	cmd.color("avCol","average")

	cmd.ray(1400,700)
	cmd.png("fig"+str(i)+"_front.png")

python end


