%CALVINCREVASSEDEPTH class definition
%
%   Usage:
%      calvingcrevassedepth=calvingcrevassedepth();

classdef calvingcrevassedepth
	properties (SetAccess=public)
		crevasse_opening_stress=1;
		critical_stress = 125e3;
		water_height = 0.;
		meltwater_mask = 1.;
		timescale=1.;
		min_iceberg_size=8000.;
		advect_icefront=1;
		propagate_from_front=1;
		remove_only_marked=0;
		hydrofracture_stabilization=0;
		hydrofracture_weakening_factor=1e-3;
		hydrofracture_min_ocean_levelset=0.;
	end
	methods
		function self = calvingcrevassedepth(varargin) % {{{
			switch nargin
				case 0
					self=setdefaultparameters(self);
				case 1
					inputstruct=varargin{1};
					list1 = properties('calvingcrevassedepth');
					list2 = fieldnames(inputstruct);
					for i=1:length(list1)
						fieldname = list1{i};
						if ismember(fieldname,list2),
							self.(fieldname) = inputstruct.(fieldname);
						end
					end
				otherwise
					error('constructor not supported');
			end
		end % }}}
		function self = extrude(self,md) % {{{
		end % }}}
		function self = setdefaultparameters(self) % {{{

			self.critical_stress         = 125e3;
			self.crevasse_opening_stress = 1;
         	self.water_height       = 0.;
			self.meltwater_mask     = 1.;
			self.timescale			= 1.;
			self.min_iceberg_size   = 8000.;
			self.advect_icefront    = 1;
			self.propagate_from_front = 1;
			self.remove_only_marked = 0;
			self.hydrofracture_stabilization = 0;
			self.hydrofracture_weakening_factor = 1e-3;
			self.hydrofracture_min_ocean_levelset = 0.;
		end % }}}
		function md = checkconsistency(self,md,solution,analyses) % {{{
			%Early return
			if (~strcmp(solution,'TransientSolution') | md.transient.ismovingfront==0), return; end
			md = checkfield(md,'fieldname','calving.crevasse_opening_stress','numel',[1],'values',[0,1,2,3,4,5]);
			md = checkfield(md,'fieldname','calving.critical_stress','numel',[1],'>',0);
			md = checkfield(md,'fieldname','calving.water_height','NaN',1,'Inf',1,'timeseries',1,'>=',0);
			if numel(self.meltwater_mask)==1
				md = checkfield(md,'fieldname','calving.meltwater_mask','NaN',1,'Inf',1,'values',[0 1]);
			else
				md = checkfield(md,'fieldname','calving.meltwater_mask','NaN',1,'Inf',1,'timeseries',1,'>=',0,'<=',1);
				maskvalues=self.meltwater_mask;
				if size(maskvalues,1)==md.mesh.numberofvertices+1
					maskvalues=maskvalues(1:end-1,:);
				end
				if any(~ismember(maskvalues(:),[0 1]))
					md = checkmessage(md,'field ''calving.meltwater_mask'' should contain only binary values 0 or 1');
				end
			end
			md = checkfield(md,'fieldname','calving.min_iceberg_size','numel',[1],'>=',0);
			md = checkfield(md,'fieldname','calving.advect_icefront','numel',[1],'values',[0,1]);
			md = checkfield(md,'fieldname','calving.propagate_from_front','numel',[1],'values',[0,1]);
			md = checkfield(md,'fieldname','calving.remove_only_marked','numel',[1],'values',[0,1]);
			md = checkfield(md,'fieldname','calving.hydrofracture_stabilization','numel',[1],'values',[0,1]);
			md = checkfield(md,'fieldname','calving.hydrofracture_weakening_factor','numel',[1],'>',0,'<=',1);
			md = checkfield(md,'fieldname','calving.hydrofracture_min_ocean_levelset','numel',[1],'<=',0);
		end % }}}
		function disp(self) % {{{
			disp(sprintf('   Calving Pi parameters:'));
			fielddisplay(self,'crevasse_opening_stress','0: stress only in the ice-flow direction, 1: max principal, 2: buttressing based');
			fielddisplay(self,'critical_stress','absolute hydrofracture critical stress threshold r_xx [Pa]');
			fielddisplay(self,'water_height','water height in the crevasse [m]');
			fielddisplay(self,'meltwater_mask','binary nodal mask gating crevasse-depth calving (0: disabled, 1: enabled)');
			fielddisplay(self,'timescale','how often to apply CD calving criterion [yr]. set to 0.0 or negative if want calving at every timestep.');
			fielddisplay(self,'min_iceberg_size','minimum iceberg size [m] for calving to occur');
			fielddisplay(self,'advect_icefront','whether to advect ice front (0: false, 1: true)');
			fielddisplay(self,'propagate_from_front','whether to seed calving propagation from the current ice front (0: check all floating ice directly, 1: propagate from front)');
			fielddisplay(self,'remove_only_marked','whether to remove only marked critical/propagated nodes (1) instead of all downstream floating ice (0)');
			fielddisplay(self,'hydrofracture_stabilization','whether to regularize hydrofracture by weakening marked ice instead of immediately removing it (0: disabled, 1: enabled)');
			fielddisplay(self,'hydrofracture_weakening_factor','local multiplier applied to rheology B in hydrofracture-regularized weak ice');
			fielddisplay(self,'hydrofracture_min_ocean_levelset','minimum MaskOceanLevelset value eligible for hydrofracture [m]; set negative to protect a floating buffer downstream of the grounding line');

		end % }}}
		function marshall(self,prefix,md,fid) % {{{
			yts=md.constants.yts;
			WriteData(fid,prefix,'name','md.calving.law','data',6,'format','Integer');
			WriteData(fid,prefix,'object',self,'fieldname','crevasse_opening_stress','format','Integer');
			WriteData(fid,prefix,'object',self,'fieldname','critical_stress','format','Double');
         	WriteData(fid,prefix,'object',self,'fieldname','timescale','format','Double');
			WriteData(fid,prefix,'object',self,'fieldname','min_iceberg_size','format','Double');
			WriteData(fid,prefix,'object',self,'fieldname','advect_icefront','format','Boolean');
			WriteData(fid,prefix,'object',self,'fieldname','propagate_from_front','format','Boolean');
			WriteData(fid,prefix,'object',self,'fieldname','remove_only_marked','format','Boolean');
			WriteData(fid,prefix,'object',self,'fieldname','hydrofracture_stabilization','format','Boolean');
			WriteData(fid,prefix,'object',self,'fieldname','hydrofracture_weakening_factor','format','Double');
			WriteData(fid,prefix,'object',self,'fieldname','hydrofracture_min_ocean_levelset','format','Double');
			WriteData(fid,prefix,'object',self,'fieldname','water_height','format','DoubleMat','mattype',1,'timeserieslength',md.mesh.numberofvertices+1,'yts',md.constants.yts);
			WriteData(fid,prefix,'object',self,'fieldname','meltwater_mask','format','DoubleMat','mattype',1,'timeserieslength',md.mesh.numberofvertices+1,'yts',md.constants.yts);
		end % }}}
	end
end
