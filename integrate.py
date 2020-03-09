
def create(filename,protein,lst_name,color_list):
	file = open(filename + ".py","w+")
	lst = ["import __main__", "__main__.pymol_argv = [ 'pymol', '-qei']",
	"import pymol", "pymol.finish_launching", "from pymol import cmd","cmd.stereo('walleye')",
	"cmd.set('stereo_shift',0.23)","cmd.set('stereo_angle',1.0)","cmd.fetch(" + "'"+protein +"'" +")","cmd.set_color(" + "'" +lst_name+
	"'"+','+str(color_list)+')']
	for x in lst:
		file.write(x + "\r\n")








'''
import __main__
__main__.pymol_argv = [ 'pymol', '-qei' ]
 

import pymol
 
# Call the function below before using any PyMOL modules.
pymol.finish_launching()
 
from pymol import cmd
cmd.stereo('walleye')
cmd.set('stereo_shift', 0.23)
cmd.set('stereo_angle', 1.0)
'''

