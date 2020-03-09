import __main__
__main__.pymol_argv = [ 'pymol', '-qei']
import pymol
pymol.finish_launching
from pymol import cmd
cmd.stereo('walleye')
cmd.set('stereo_shift',0.23)
cmd.set('stereo_angle',1.0)
cmd.fetch('1J6Z')
cmd.png(5 , ray = 1, quiet = 1)
cmd.quit()
