def create(filename,save_to, pdb_id,lst_name = [],color_list = []):
	file = open(filename,"w+")
	lst = [
	"import __main__", 
	"__main__.pymol_argv = [ 'pymol', '-qei']",
	"import pymol", 
	"pymol.finish_launching", 
	"from pymol import cmd",
	"cmd.stereo('walleye')",
	"cmd.set('stereo_shift',0.23)",
	"cmd.set('stereo_angle',1.0)",
	"cmd.fetch(" + "'"+pdb_id +"'" +")",
	#"cmd.set_color(" + "'" +lst_name+
	#"'"+','+str(color_list)+')'
	"cmd.png(" + save_to + " , ray = 1, quiet = 1)",
	"cmd.quit()"
	]
	for x in lst:
		file.write(x + "\r\n")
