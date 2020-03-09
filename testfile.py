import __main__
__main__.pymol_argv = [ 'pymol', '-qei']
import pymol
pymol.finish_launching
from pymol import cmd
cmd.stereo('walleye')
cmd.set('stereo_shift',0.23)
cmd.set('stereo_angle',1.0)
cmd.fetch('1J6Z')
cmd.set_color('list',[0.1, 0.2, 0.3])
