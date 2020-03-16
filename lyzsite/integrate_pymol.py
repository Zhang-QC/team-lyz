def val_to_color(col):
	'''
	Taking in a float value between zero and one and generate
	RGB color list [#, #, #] that can be read by pyMol in the 
	blue-red spectrum.

	Input:
		col: a float between 0 and 1

	Output:
		a list of third floats
	'''
	return [col, 0, 1 - col]


def normalize(lst):
	'''
	Taking in a list of floats, normalized them by dividing them
	based on the largest value.
	'''
	l = []
	m = max(lst)
	for i in lst:
		l.append(i / m)
	return l


def create(filename, save_to_one, save_to_two, pdb_id,col_lst):
	file = open(filename,"w+")
	lst = [
	"import __main__", 
	"__main__.pymol_argv = [ 'pymol', '-qei']",
	"import pymol", 
	"pymol.finish_launching", 
	"from pymol import cmd",
	#"cmd.stereo('walleye')",
	#"cmd.set('stereo_shift',0.23)",
	#"cmd.set('stereo_angle',1.0)",
	"cmd.fetch(" + "'"+pdb_id +"'" +")",
	"cmd.show('surface', 'all')",
	"cmd.remove('solvent')",
	"cmd.set_color('[0.00 , 0.00 , 1.00]', [0.00, 0.00 , 1.00])",
	"cmd.color('[0.00 , 0.00 , 1.00]','all')",
	"cmd.set('transparency', 0.5)",
	"cmd.set_color('[1.00 , 0.00 , 0.00]', [1.00 , 0.00 , 0.00])"]
	l = normalize(col_lst)
	print(l)
	for index, val in enumerate(l):
		col = str(val_to_color(val))
		st1 = "cmd.show('sphere', 'resi " + str(index) +"')"
		st2 = "cmd.set_color('"+col+"', "+col+")"
		st3 = "cmd.color('"+col+"','resi "+ str(index) +"')"
		st4 = "cmd.set('transparency', "+ str(0.5*val)+")"
		lst += [st1, st2, st3, st4]
	lst += [
	"cmd.png(" + save_to_one + " , ray = 1, quiet = 1)",
	"cmd.rotate('z', 180, 'all')",
	"cmd.png(" + save_to_two + " , ray = 1, quiet = 1)",
	"cmd.quit()"
	]
	for x in lst:
		file.write(x + "\r\n")
