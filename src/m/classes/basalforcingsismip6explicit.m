%ISMIP6 EXPLICIT BASAL FORCINGS class definition
%
%   Usage:
%      basalforcingsismip6explicit=basalforcingsismip6explicit();

classdef basalforcingsismip6explicit
	properties (SetAccess=public) 
		num_basins                = 1;
		basin_id                  = NaN;
		gamma_0                   = 0.;
		ocean_temperature         = NaN;
		ocean_salinity            = NaN;
		lambda_1                  = 0.;
		lambda_2                  = 0.;
		lambda_3                  = 0.;
		tf_depths                 = NaN;
		delta_t                   = NaN;
		islocal                   = 0;
		isslope                   = 0;
		geothermalflux            = NaN;
		groundedice_melting_rate  = NaN;
		melt_anomaly              = NaN;
	end
	methods
		function self = extrude(self,md) % {{{
			self.basin_id=project3d(md,'vector',self.basin_id,'type','element','layer',1);
			
			self.geothermalflux=project3d(md,'vector',self.geothermalflux,'type','element','layer',1); %bedrock only gets geothermal flux
			self.groundedice_melting_rate=project3d(md,'vector',self.groundedice_melting_rate,'type','node','layer',1);
			self.melt_anomaly=project3d(md,'vector',self.melt_anomaly,'type','element','layer',1); %bedrock only gets geothermal flux
		end % }}}
		function self = basalforcingsismip6explicit(varargin) % {{{
			switch nargin
				case 0
					self=setdefaultparameters(self);
				case 1
					self=setdefaultparameters(self);
					self=structtoobj(self,varargin{1});
				otherwise
					error('constructor not supported');
			end
		end % }}}
		function self = initialize(self,md) % {{{
			if self.gamma_0 == 0,
				self.gamma_0 = 14477;
				disp('      no basalforcings.gamma_0 specified: value set to 14477 m/yr');
			end
			if isnan(self.groundedice_melting_rate),
				self.groundedice_melting_rate=zeros(md.mesh.numberofvertices,1);
				disp('      no basalforcings.groundedice_melting_rate specified: values set as zero');
			end

		end % }}}
		function self = setdefaultparameters(self) % {{{
			self.gamma_0 = 14477; %m/yr
			self.lambda_1 = -0.0573; % degC / PSU
			self.lambda_2 = 0.0832; % degC
			self.lambda_3 = 7.61e-4; % degC / m
			self.islocal = false;
			self.isslope = 0;
		end % }}}
		function md = checkconsistency(self,md,solution,analyses) % {{{

			md = checkfield(md,'fieldname','basalforcings.num_basins','numel',1,'NaN',1,'Inf',1,'>',0);
			md = checkfield(md,'fieldname','basalforcings.basin_id','Inf',1,'>=',0,'<=',md.basalforcings.num_basins,'size',[md.mesh.numberofelements 1]);
			md = checkfield(md,'fieldname','basalforcings.gamma_0','numel',1,'NaN',1,'Inf',1,'>=',0);
			md = checkfield(md,'fieldname','basalforcings.lambda_1','numel',1,'NaN',1,'Inf',1);
			md = checkfield(md,'fieldname','basalforcings.lambda_2','numel',1,'NaN',1,'Inf',1);
			md = checkfield(md,'fieldname','basalforcings.lambda_3','numel',1,'NaN',1,'Inf',1);
			md = checkfield(md,'fieldname','basalforcings.tf_depths','NaN',1,'Inf',1,'size',[1,NaN],'<=',0);
			md = checkfield(md,'fieldname','basalforcings.delta_t','NaN',1,'Inf',1,'numel',md.basalforcings.num_basins,'size',[1,md.basalforcings.num_basins]);
			md = checkfield(md,'fieldname','basalforcings.islocal','values',[0 1]);
			md = checkfield(md,'fieldname','basalforcings.isslope','values',[0 1]);
			md = checkfield(md,'fieldname','basalforcings.geothermalflux','NaN',1,'Inf',1,'>=',0,'timeseries',1);
			md = checkfield(md,'fieldname','basalforcings.groundedice_melting_rate','NaN',1,'Inf',1,'timeseries',1);
			if length(md.basalforcings.melt_anomaly)>1,
				md = checkfield(md,'fieldname','basalforcings.melt_anomaly','NaN',1,'Inf',1,'timeseries',1);
			end

			% md = checkfield(md,'fieldname','basalforcings.ocean_temperature','size',[1,1,numel(md.basalforcings.tf_depths)]);
			% md = checkfield(md,'fieldname','basalforcings.ocean_salinity','size',[1,1,numel(md.basalforcings.tf_depths)]);
			% for i=1:numel(md.basalforcings.tf_depths)
			% 	md = checkfield(md,'fieldname',['basalforcings.ocean_temperature{' num2str(i) '}'],'field',md.basalforcings.ocean_temperature{i},'size',[md.mesh.numberofvertices+1 NaN],'NaN',1,'Inf',1,'timeseries',1);
			% 	md = checkfield(md,'fieldname',['basalforcings.ocean_salinity{' num2str(i) '}'],'field',md.basalforcings.ocean_salinity{i},'size',[md.mesh.numberofvertices+1 NaN],'NaN',1,'Inf',1,'>=',0,'timeseries',1);
			% end

		end % }}}
		function disp(self) % {{{
			disp(sprintf('   ISMIP6 EXPLICIT basal melt rate parameterization:'));
			fielddisplay(self,'num_basins','number of basins the model domain is partitioned into [unitless]');
			fielddisplay(self,'basin_id','basin number assigned to each node (unitless)');
			fielddisplay(self,'gamma_0','melt rate coefficient (m/yr)');
			fielddisplay(self,'lambda_1','liquidus constant 1 (degree C / PSU)');
			fielddisplay(self,'lambda_2','liquidus constant 2 (degree C)');
			fielddisplay(self,'lambda_3','liquidus constant 3 (degree C / m)');
			fielddisplay(self,'tf_depths','elevation of vertical layers in ocean thermal forcing dataset');
			fielddisplay(self,'ocean_temperature','ocean temperature (degrees C)');
			fielddisplay(self,'ocean_salinity','ocean salinity (PSU)');
			fielddisplay(self,'delta_t','Ocean temperature correction per basin (degrees C)');
			fielddisplay(self,'islocal','boolean to use the local version of the ISMIP6 melt rate parameterization (default false)');
			fielddisplay(self,'isslope','boolean to use slope dependent melting (default 0)');
			fielddisplay(self,'geothermalflux','geothermal heat flux (W/m^2)');
			fielddisplay(self,'groundedice_melting_rate','basal melting rate (positive if melting) (m/yr)');
			fielddisplay(self,'melt_anomaly','floating ice basal melt anomaly (m/yr)');

		end % }}}
		function marshall(self,prefix,md,fid) % {{{

			yts=md.constants.yts;

			WriteData(fid,prefix,'name','md.basalforcings.model','data',10,'format','Integer');
			WriteData(fid,prefix,'object',self,'fieldname','num_basins','format','Integer');
			WriteData(fid,prefix,'object',self,'fieldname','basin_id','data',self.basin_id-1,'name','md.basalforcings.basin_id','format','IntMat','mattype',2);   %0-indexed
			WriteData(fid,prefix,'object',self,'fieldname','gamma_0','format','Double','scale',1./yts);
			WriteData(fid,prefix,'object',self,'fieldname','lambda_1','format','Double');
			WriteData(fid,prefix,'object',self,'fieldname','lambda_2','format','Double');
			WriteData(fid,prefix,'object',self,'fieldname','lambda_3','format','Double');
			WriteData(fid,prefix,'object',self,'fieldname','tf_depths','format','DoubleMat','name','md.basalforcings.tf_depths');
			WriteData(fid,prefix,'object',self,'fieldname','ocean_temperature','format','MatArray','name','md.basalforcings.ocean_temperature','timeserieslength',md.mesh.numberofvertices+1,'yts',md.constants.yts);
			WriteData(fid,prefix,'object',self,'fieldname','ocean_salinity','format','MatArray','name','md.basalforcings.ocean_salinity','timeserieslength',md.mesh.numberofvertices+1,'yts',md.constants.yts);
			WriteData(fid,prefix,'object',self,'fieldname','delta_t','format','DoubleMat','name','md.basalforcings.delta_t','timeserieslength',md.mesh.numberofvertices+1,'yts',md.constants.yts);
			WriteData(fid,prefix,'object',self,'fieldname','islocal','format','Boolean');
			WriteData(fid,prefix,'object',self,'fieldname','isslope','format','Integer');
			WriteData(fid,prefix,'object',self,'fieldname','geothermalflux','format','DoubleMat','name','md.basalforcings.geothermalflux','mattype',1,'timeserieslength',md.mesh.numberofelements+1,'yts',md.constants.yts);
			WriteData(fid,prefix,'object',self,'fieldname','groundedice_melting_rate','format','DoubleMat','mattype',1,'scale',1./yts,'timeserieslength',md.mesh.numberofvertices+1,'yts',md.constants.yts);
			WriteData(fid,prefix,'object',self,'fieldname','melt_anomaly','format','DoubleMat','mattype',1,'scale',1./yts,'timeserieslength',md.mesh.numberofvertices+1,'yts',md.constants.yts);

		end % }}}
	end
end
