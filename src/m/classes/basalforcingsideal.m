%BASAL FORCINGS class definition
%
%   Usage:
%      basalforcings=basalforcings();

classdef basalforcingsideal
	properties (SetAccess=public) 
		alpha = NaN;
		beta  = NaN;
		gamma = NaN;
		m0    = NaN;
		groundedice_melting_rate = NaN;
		geothermalflux = NaN;
	end
	methods
		function self = basalforcingsideal(varargin) % {{{
			switch nargin
				case 0
					self=setdefaultparameters(self);
				case 1
					self =structtoobj(basalforcingsideal(),varargin{1});
				otherwise
					error('constructor not supported');
			end
		end % }}}
		function disp(self) % {{{
			disp(sprintf('   ideal basal forcings parameters:'));

			fielddisplay(self,'alpha','exponent controlling the concentration of melt in the GZ. higher alpha is higher concentration');
			fielddisplay(self,'beta','exponent controlling the spread of melt over the ice shelf. higher beta is more evenly spread');
			fielddisplay(self,'gamma','exponent controlling the spread of melt across the width. higher gamma is more concentrated towards the north.');
			fielddisplay(self,'m0','constant that scales the melt');
			fielddisplay(self,'groundedice_melting_rate','melt rate applied on grounded ice');
			fielddisplay(self,'geothermalflux','melt rate applied on grounded ice');

		end % }}}
		function self = extrude(self,md) % {{{
			self.groundedice_melting_rate=project3d(md,'vector',self.groundedice_melting_rate,'type','node','layer',1);
			self.geothermalflux=project3d(md,'vector',self.geothermalflux,'type','element','layer',1); %bedrock only gets geothermal flux
		end % }}}
		function self = initialize(self,md) % {{{

		end % }}}
		function self = setdefaultparameters(self) % {{{
			self.alpha = 1/3.0;
			self.beta = 1.1;
			self.gamma = 0;
			self.m0 = 3.0;
		end % }}}
		function md = checkconsistency(self,md,solution,analyses) % {{{
			md = checkfield(md,'fieldname','basalforcings.alpha','NaN',1,'Inf',1,'>=', 0.0);
			md = checkfield(md,'fieldname','basalforcings.beta','NaN',1,'Inf',1,'>=', 1.0);
			md = checkfield(md,'fieldname','basalforcings.gamma','NaN',1,'Inf',1,'>=', 0.0);
			md = checkfield(md,'fieldname','basalforcings.m0','NaN',1,'Inf',1,'>=', 0.0);
			md = checkfield(md,'fieldname','basalforcings.groundedice_melting_rate','NaN',1,'Inf',1,'timeseries',1);
			md = checkfield(md,'fieldname','basalforcings.geothermalflux','NaN',1,'Inf',1,'timeseries',1);
		end % }}}
		function marshall(self,prefix,md,fid) % {{{
			WriteData(fid,prefix,'name','md.basalforcings.model','data',11,'format','Integer');
			WriteData(fid,prefix,'object',self,'fieldname','alpha','format','Double');
			WriteData(fid,prefix,'object',self,'fieldname','beta','format','Double');
			WriteData(fid,prefix,'object',self,'fieldname','gamma','format','Double');
			WriteData(fid,prefix,'object',self,'fieldname','m0','format','Double','scale',1./md.constants.yts);
			WriteData(fid,prefix,'object',self,'fieldname','groundedice_melting_rate','format','DoubleMat','mattype',1,'scale',1./md.constants.yts,'timeserieslength',md.mesh.numberofvertices+1);
			WriteData(fid,prefix,'object',self,'fieldname','geothermalflux','format','DoubleMat','name','md.basalforcings.geothermalflux','mattype',1,'timeserieslength',md.mesh.numberofelements+1,'yts',md.constants.yts);
		end % }}}
		function savemodeljs(self,fid,modelname) % {{{
		
			writejs1Darray(fid,[modelname '.basalforcings.alpha'],self.alpha);
			writejs1Darray(fid,[modelname '.basalforcings.beta'],self.beta);
			writejs1Darray(fid,[modelname '.basalforcings.gamma'],self.gamma);
			writejs1Darray(fid,[modelname '.basalforcings.m0'],self.m0);
			writejs1Darray(fid,[modelname '.basalforcings.groundedice_melting_rate'],self.groundedice_melting_rate);
			writejs1Darray(fid,[modelname '.basalforcings.geothermalflux'],self.geothermalflux);

		end % }}}
	end
end
