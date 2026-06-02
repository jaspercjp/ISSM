from checkfield import checkfield
from fielddisplay import fielddisplay
#from project3d import project3d # Uncomment if/when extrude is implemented
from WriteData import WriteData
import numpy as np


class calvingcrevassedepth(object):
    """CALVINCREVASSEDEPTH class definition

    Usage:
        calvingcrevassedepth = calvingcrevassedepth()
    """

    def __init__(self):  # {{{
        self.crevasse_opening_stress = 1
        self.critical_stress         = 125e3
        self.water_height            = 0.
        self.meltwater_mask          = 1.
        self.timescale               = 1.
        self.min_iceberg_size        = 8000.
        self.advect_icefront         = 1
        self.propagate_from_front    = 1
        self.remove_only_marked      = 0
        self.hydrofracture_stabilization = 0
        self.hydrofracture_weakening_factor = 1e-3
        self.hydrofracture_min_ocean_levelset = 0.

        #self.setdefaultparameters() # Uncomment if/when setdefaultparameters is used
    # }}}
    def __repr__(self):  # {{{
        s = '   Calving Pi parameters:'
        s += '{}\n'.format(fielddisplay(self, 'crevasse_opening_stress', '0: stress only in the ice-flow direction, 1: max principal'))
        s += '{}\n'.format(fielddisplay(self, 'critical_stress', 'absolute hydrofracture critical stress threshold r_xx [Pa]'))
        s += '{}\n'.format(fielddisplay(self, 'water_height', 'water height in the crevasse [m]'))
        s += '{}\n'.format(fielddisplay(self, 'meltwater_mask', 'binary nodal mask gating crevasse-depth calving (0: disabled, 1: enabled)'))
        s += '{}\n'.format(fielddisplay(self, 'timescale', 'how often to apply CD calving criterion [yr]. set to 0.0 or negative if want calving at every timestep.'))
        s += '{}\n'.format(fielddisplay(self, 'min_iceberg_size', 'minimum iceberg size [m] for calving to occur'))
        s += '{}\n'.format(fielddisplay(self, 'advect_icefront', 'whether to advect ice front (0: false, 1: true)'))
        s += '{}\n'.format(fielddisplay(self, 'propagate_from_front', 'whether to seed calving propagation from the current ice front (0: check all floating ice directly, 1: propagate from front)'))
        s += '{}\n'.format(fielddisplay(self, 'remove_only_marked', 'whether to remove only marked critical/propagated nodes (1) instead of all downstream floating ice (0)'))
        s += '{}\n'.format(fielddisplay(self, 'hydrofracture_stabilization', 'whether to regularize hydrofracture by weakening marked ice instead of immediately removing it (0: disabled, 1: enabled)'))
        s += '{}\n'.format(fielddisplay(self, 'hydrofracture_weakening_factor', 'local multiplier applied to rheology B in hydrofracture-regularized weak ice'))
        s += '{}\n'.format(fielddisplay(self, 'hydrofracture_min_ocean_levelset', 'minimum MaskOceanLevelset value eligible for hydrofracture [m]; set negative to protect a floating buffer downstream of the grounding line'))
        return s
    # }}}
    def setdefaultparameters(self):  # {{{
        return self
    # }}}
    def extrude(self, md):  # {{{
        return self
    # }}}
    def checkconsistency(self, md, solution, analyses):  # {{{
        #Early return
        if solution != 'TransientSolution' or not md.transient.ismovingfront:
            return md

        md = checkfield(md, 'fieldname', 'calving.crevasse_opening_stress', 'numel', [1], 'values', [0,1,2,3,4])
        md = checkfield(md, 'fieldname', 'calving.critical_stress', 'numel', [1], '>', 0.)
        md = checkfield(md, 'fieldname', 'calving.water_height', 'NaN', 1, 'Inf', 1, 'timeseries', 1, '>=', 0) 
        meltwater_mask = np.asarray(self.meltwater_mask)
        if meltwater_mask.size == 1:
            md = checkfield(md, 'fieldname', 'calving.meltwater_mask', 'NaN', 1, 'Inf', 1, 'values', [0, 1])
        else:
            md = checkfield(md, 'fieldname', 'calving.meltwater_mask', 'NaN', 1, 'Inf', 1, 'timeseries', 1, '>=', 0, '<=', 1)
            maskvalues = meltwater_mask
            if maskvalues.ndim > 1 and maskvalues.shape[0] == md.mesh.numberofvertices + 1:
                maskvalues = maskvalues[:-1, :]
            md = checkfield(md, 'fieldname', 'calving.meltwater_mask', 'field', maskvalues, 'values', [0, 1])
        md = checkfield(md, 'fieldname', 'calving.min_iceberg_size', 'numel', [1], '>=', 0)
        md = checkfield(md, 'fieldname', 'calving.advect_icefront', 'numel', [1], 'values', [0, 1])
        md = checkfield(md, 'fieldname', 'calving.propagate_from_front', 'numel', [1], 'values', [0, 1])
        md = checkfield(md, 'fieldname', 'calving.remove_only_marked', 'numel', [1], 'values', [0, 1])
        md = checkfield(md, 'fieldname', 'calving.hydrofracture_stabilization', 'numel', [1], 'values', [0, 1])
        md = checkfield(md, 'fieldname', 'calving.hydrofracture_weakening_factor', 'numel', [1], '>', 0., '<=', 1.)
        md = checkfield(md, 'fieldname', 'calving.hydrofracture_min_ocean_levelset', 'numel', [1], '<=', 0.)

        return md
    # }}}
    def marshall(self, prefix, md, fid):  # {{{
        yts = md.constants.yts
        WriteData(fid, prefix, 'name', 'md.calving.law', 'data', 6, 'format', 'Integer')
        WriteData(fid, prefix, 'object', self, 'fieldname', 'crevasse_opening_stress', 'format', 'Integer')
        WriteData(fid, prefix, 'object', self, 'fieldname', 'critical_stress', 'format', 'Double')
        WriteData(fid, prefix, 'object', self, 'fieldname', 'timescale', 'format', 'Double')
        WriteData(fid, prefix, 'object', self, 'fieldname', 'min_iceberg_size', 'format', 'Double')
        WriteData(fid, prefix, 'object', self, 'fieldname', 'advect_icefront', 'format', 'Boolean')
        WriteData(fid, prefix, 'object', self, 'fieldname', 'propagate_from_front', 'format', 'Boolean')
        WriteData(fid, prefix, 'object', self, 'fieldname', 'remove_only_marked', 'format', 'Boolean')
        WriteData(fid, prefix, 'object', self, 'fieldname', 'hydrofracture_stabilization', 'format', 'Boolean')
        WriteData(fid, prefix, 'object', self, 'fieldname', 'hydrofracture_weakening_factor', 'format', 'Double')
        WriteData(fid, prefix, 'object', self, 'fieldname', 'hydrofracture_min_ocean_levelset', 'format', 'Double')
        WriteData(fid, prefix, 'object', self, 'fieldname', 'water_height', 'format', 'DoubleMat', 'mattype', 1, 'timeserieslength', md.mesh.numberofvertices + 1, 'yts', md.constants.yts)
        WriteData(fid, prefix, 'object', self, 'fieldname', 'meltwater_mask', 'format', 'DoubleMat', 'mattype', 1, 'timeserieslength', md.mesh.numberofvertices + 1, 'yts', md.constants.yts)
    # }}}
