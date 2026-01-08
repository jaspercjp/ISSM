/*!\file FloatingiceMeltingRatex
 * \brief: calculates Floating ice melting rate
 */

#include "./FloatingiceMeltingRatex.h"
#include "../../shared/shared.h"
#include "../../toolkits/toolkits.h"
#include "./../../classes/Inputs/DatasetInput.h"
#include "../InputDuplicatex/InputDuplicatex.h"

void FloatingiceMeltingRatex(FemModel* femmodel){/*{{{*/

	/*Intermediaties*/
	int basalforcing_model;
	int melt_style;

	/*First, get melt_interpolation model from parameters*/
	femmodel->parameters->FindParam(&melt_style,GroundinglineMeltInterpolationEnum);
	if(melt_style==IntrusionMeltEnum){
		InputDuplicatex(femmodel,MaskOceanLevelsetEnum,DistanceToGroundinglineEnum);//FIXME Duplicate first so that it can preserve the sign
		femmodel->DistanceToFieldValue(MaskOceanLevelsetEnum,0.,DistanceToGroundinglineEnum);
	}

	/*First, get BMB model from parameters*/
	femmodel->parameters->FindParam(&basalforcing_model,BasalforcingsEnum);
	if(basalforcing_model==8) basalforcing_model=BasalforcingsIsmip6ExplicitEnum;

	/*branch to correct module*/
	switch(basalforcing_model){
		case FloatingMeltRateEnum:
		case MantlePlumeGeothermalFluxEnum:
			/*Nothing to be done*/
			break;
		case LinearFloatingMeltRateEnum:
			if(VerboseSolution())_printf0_("	  call Linear Floating melting rate module\n");
			LinearFloatingiceMeltingRatex(femmodel);
			break;
		case MismipFloatingMeltRateEnum:
			if(VerboseSolution())_printf0_("	  call Mismip Floating melting rate module\n");
			MismipFloatingiceMeltingRatex(femmodel);
			break;
		case SpatialLinearFloatingMeltRateEnum:
			if(VerboseSolution())_printf0_("	  call Spatial Linear Floating melting rate module\n");
			SpatialLinearFloatingiceMeltingRatex(femmodel);
			break;
		case BasalforcingsPicoEnum:
			if(VerboseSolution())_printf0_("   call Pico Floating melting rate module\n");
			FloatingiceMeltingRatePicox(femmodel);
			break;
		case BasalforcingsIsmip6Enum:
			if(VerboseSolution())_printf0_("   call ISMIP 6 Floating melting rate module\n");
			FloatingiceMeltingRateIsmip6x(femmodel);
			break;
		case BasalforcingsIsmip6ExplicitEnum:
			if(VerboseSolution())_printf0_("   call ISMIP 6 Explicit Floating melting rate module\n");
			FloatingiceMeltingRateIsmip6Explicitx(femmodel);
			break;
		case BeckmannGoosseFloatingMeltRateEnum:
			if(VerboseSolution())_printf0_("        call BeckmannGoosse Floating melting rate module\n");
			BeckmannGoosseFloatingiceMeltingRatex(femmodel);
			break;
		case LinearFloatingMeltRatearmaEnum:
			if(VerboseSolution())_printf0_("        call Linear Floating melting rate ARMA module\n");
			LinearFloatingiceMeltingRatearmax(femmodel);
			break;
		default:
			_error_("Basal forcing model "<<EnumToStringx(basalforcing_model)<<" not supported yet");
	}

}/*}}}*/

void LinearFloatingiceMeltingRatex(FemModel* femmodel){/*{{{*/

	for(Object* & object : femmodel->elements->objects){
      Element* element = xDynamicCast<Element*>(object);
		element->LinearFloatingiceMeltingRate();
	}

}/*}}}*/
void SpatialLinearFloatingiceMeltingRatex(FemModel* femmodel){/*{{{*/

	for(Object* & object : femmodel->elements->objects){
      Element* element = xDynamicCast<Element*>(object);
		element->SpatialLinearFloatingiceMeltingRate();
	}

}/*}}}*/
void MismipFloatingiceMeltingRatex(FemModel* femmodel){/*{{{*/

	for(Object* & object : femmodel->elements->objects){
      Element* element = xDynamicCast<Element*>(object);
		element->MismipFloatingiceMeltingRate();
	}
}
/*}}}*/
void FloatingiceMeltingRateIsmip6x(FemModel* femmodel){/*{{{*/

	int         num_basins, basinid, num_depths, domaintype;
	IssmDouble  area, tf, base, time;
	bool        islocal;
	IssmDouble* tf_depths = NULL;

	femmodel->parameters->FindParam(&num_basins,BasalforcingsIsmip6NumBasinsEnum);
	femmodel->parameters->FindParam(&tf_depths,&num_depths,BasalforcingsIsmip6TfDepthsEnum); _assert_(tf_depths);
	femmodel->parameters->FindParam(&islocal,BasalforcingsIsmip6IsLocalEnum);

	/*Binary search works for vectors that are sorted in increasing order only, make depths positive*/
	for(int i=0;i<num_depths;i++) tf_depths[i] = -tf_depths[i];

	IssmDouble* tf_weighted_avg     = xNewZeroInit<IssmDouble>(num_basins);
	IssmDouble* tf_weighted_avg_cpu = xNewZeroInit<IssmDouble>(num_basins);
	IssmDouble* areas_summed        = xNewZeroInit<IssmDouble>(num_basins);
	IssmDouble* areas_summed_cpu    = xNewZeroInit<IssmDouble>(num_basins);

	/*Get TF at each ice shelf point - linearly intepolate in depth and time*/
	for(Object* & object : femmodel->elements->objects){
      Element* element = xDynamicCast<Element*>(object);
		int      numvertices = element->GetNumberOfVertices();

		/*Set melt to 0 if non floating*/
		if(!element->IsIceInElement() || !element->IsAllFloating() || !element->IsOnBase()){
			IssmDouble* values = xNewZeroInit<IssmDouble>(numvertices);
			element->AddInput(BasalforcingsFloatingiceMeltingRateEnum,values,P1DGEnum);
			element->AddInput(BasalforcingsIsmip6TfShelfEnum,values,P1DGEnum);
			xDelete<IssmDouble>(values);
			continue;
		}

		/*Get TF on all vertices using linear interpolation*/
		IssmDouble*    tf_test        = xNew<IssmDouble>(numvertices);
		IssmDouble*    depth_vertices = xNew<IssmDouble>(numvertices);
		DatasetInput* tf_input = element->GetDatasetInput(BasalforcingsIsmip6TfEnum); _assert_(tf_input);

		element->GetInputListOnVertices(&depth_vertices[0],BaseEnum);

		Gauss* gauss=element->NewGauss();
		for(int iv=0;iv<numvertices;iv++){
			gauss->GaussVertex(iv);

			/*Find out where the ice shelf base is within tf_depths*/
			IssmDouble depth = -depth_vertices[iv]; /*NOTE: make sure we are dealing with depth>0*/
			int offset;
			int found=binary_search(&offset,depth,tf_depths,num_depths);
			if(!found) _error_("depth not found");

			if (offset==-1){
				/*get values for the first depth: */
				_assert_(depth<=tf_depths[0]);
				tf_input->GetInputValue(&tf_test[iv],gauss,0);
			}
			else if(offset==num_depths-1){
				/*get values for the last time: */
				_assert_(depth>=tf_depths[num_depths-1]);
				tf_input->GetInputValue(&tf_test[iv],gauss,num_depths-1);
			}
			else {
				/*get values between two times [offset:offset+1], Interpolate linearly*/
				_assert_(depth>=tf_depths[offset] && depth<tf_depths[offset+1]);
				IssmDouble deltaz=tf_depths[offset+1]-tf_depths[offset];
				IssmDouble alpha2=(depth-tf_depths[offset])/deltaz;
				IssmDouble alpha1=(1.-alpha2);
				IssmDouble tf1,tf2;
				tf_input->GetInputValue(&tf1,gauss,offset);
				tf_input->GetInputValue(&tf2,gauss,offset+1);
				tf_test[iv] = alpha1*tf1 + alpha2*tf2;
			}
		}

		element->AddInput(BasalforcingsIsmip6TfShelfEnum,tf_test,P1DGEnum);
		xDelete<IssmDouble>(tf_test);
		xDelete<IssmDouble>(depth_vertices);
		delete gauss;
	}

	if(!islocal) {
		/*Compute sums of tf*area and shelf-area per cpu*/
		for(Object* & object : femmodel->elements->objects){
			Element* element = xDynamicCast<Element*>(object);
			if(!element->IsOnBase() || !element->IsIceInElement() || !element->IsAllFloating()) continue;
			/*Spawn basal element if on base to compute element area*/
			Element* basalelement = element->SpawnBasalElement();
			Input* tf_input=basalelement->GetInput(BasalforcingsIsmip6TfShelfEnum); _assert_(tf_input);
			basalelement->GetInputValue(&basinid,BasalforcingsIsmip6BasinIdEnum);
			Gauss* gauss=basalelement->NewGauss(1); gauss->GaussPoint(0);
			tf_input->GetInputValue(&tf,gauss);
			delete gauss;
			area=basalelement->GetHorizontalSurfaceArea();
			tf_weighted_avg[basinid]+=tf*area;
			areas_summed[basinid]   +=area;
			/*Delete spawned element if we are in 3D*/
			basalelement->FindParam(&domaintype,DomainTypeEnum);
			if(basalelement->IsSpawnedElement()){basalelement->DeleteMaterials(); delete basalelement;};
		}

		/*Syncronize across cpus*/
		ISSM_MPI_Allreduce(tf_weighted_avg,tf_weighted_avg_cpu,num_basins,ISSM_MPI_DOUBLE,ISSM_MPI_SUM,IssmComm::GetComm());
		ISSM_MPI_Allreduce(areas_summed,areas_summed_cpu,num_basins,ISSM_MPI_DOUBLE,ISSM_MPI_SUM,IssmComm::GetComm());

		/*Make sure Area is not zero to avoid dividing by 0 if a basin is not present in the model*/
		for(int k=0;k<num_basins;k++) if(areas_summed_cpu[k]==0.) areas_summed_cpu[k] = 1.;

		/*Compute weighted means and save*/
		for(int k=0;k<num_basins;k++){tf_weighted_avg_cpu[k] = tf_weighted_avg_cpu[k]/areas_summed_cpu[k];}
		femmodel->parameters->AddObject(new DoubleVecParam(BasalforcingsIsmip6AverageTfEnum,tf_weighted_avg_cpu,num_basins));
	}

   /*Compute meltrates*/
	for(Object* & object : femmodel->elements->objects){
      Element* element = xDynamicCast<Element*>(object);
		element->Ismip6FloatingiceMeltingRate();
	}

	/*Cleanup and return */
	xDelete<IssmDouble>(tf_weighted_avg);
	xDelete<IssmDouble>(tf_weighted_avg_cpu);
	xDelete<IssmDouble>(areas_summed);
	xDelete<IssmDouble>(areas_summed_cpu);
	xDelete<IssmDouble>(tf_depths);
}
/*{{{*/
void FloatingiceMeltingRateIsmip6Explicitx(FemModel* femmodel){

	int         num_basins, basinid, num_depths, domaintype;
	IssmDouble  area, tf, base, time;
	bool        islocal;
	IssmDouble* tf_depths = NULL;
	IssmDouble lambda1, lambda2, lambda3;

	femmodel->parameters->FindParam(&num_basins,BasalforcingsIsmip6NumBasinsEnum);
	femmodel->parameters->FindParam(&tf_depths,&num_depths,BasalforcingsIsmip6TfDepthsEnum); _assert_(tf_depths);
	femmodel->parameters->FindParam(&islocal,BasalforcingsIsmip6IsLocalEnum);
	femmodel->parameters->FindParam(&lambda1,BasalforcingsIsmip6ExplicitLambda1Enum);
	femmodel->parameters->FindParam(&lambda2,BasalforcingsIsmip6ExplicitLambda2Enum);
	femmodel->parameters->FindParam(&lambda3,BasalforcingsIsmip6ExplicitLambda3Enum);

	/*Copy to local array to avoid modifying global parameters*/
	IssmDouble* tf_depths_local = xNew<IssmDouble>(num_depths);
	xMemCpy<IssmDouble>(tf_depths_local, tf_depths, num_depths);

	/*Binary search works for vectors that are sorted in increasing order only, make depths positive*/
	for(int i=0;i<num_depths;i++) tf_depths_local[i] = -tf_depths_local[i];

	bool decreasing = false;
	if (num_depths > 1 && tf_depths_local[0] > tf_depths_local[num_depths-1]) decreasing = true;

	if(decreasing){
		/*Reverse array so that it is sorted increasing for binary_search*/
		for(int i=0; i<num_depths/2; i++){
			IssmDouble temp = tf_depths_local[i];
			tf_depths_local[i] = tf_depths_local[num_depths-1-i];
			tf_depths_local[num_depths-1-i] = temp;
		}
	}

	IssmDouble* tf_weighted_avg     = xNewZeroInit<IssmDouble>(num_basins);
	IssmDouble* tf_weighted_avg_cpu = xNewZeroInit<IssmDouble>(num_basins);
	IssmDouble* areas_summed        = xNewZeroInit<IssmDouble>(num_basins);
	IssmDouble* areas_summed_cpu    = xNewZeroInit<IssmDouble>(num_basins);

	/*Get TF at each ice shelf point - linearly intepolate in depth and time*/

	_printf0_("\tcomputing thermal forcing at each element\n");
	for(Object* & object : femmodel->elements->objects){
      Element* element = xDynamicCast<Element*>(object);
		int      numvertices = element->GetNumberOfVertices();

		/*Set melt to 0 if non floating*/
		if(!element->IsIceInElement() || !element->IsAllFloating() || !element->IsOnBase()){
			IssmDouble* values = xNewZeroInit<IssmDouble>(numvertices);
			element->AddInput(BasalforcingsFloatingiceMeltingRateEnum,values,P1DGEnum);
			element->AddInput(BasalforcingsIsmip6TfShelfEnum,values,P1DGEnum);
			xDelete<IssmDouble>(values);
			continue;
		}

		/*Get TF on all vertices using linear interpolation*/
		IssmDouble*    tf_test        = xNew<IssmDouble>(numvertices);
		IssmDouble*    depth_vertices = xNew<IssmDouble>(numvertices);
		DatasetInput* temp_input = element->GetDatasetInput(BasalforcingsIsmip6ExplicitOceanTemperatureInputEnum); _assert_(temp_input);
		DatasetInput* sal_input = element->GetDatasetInput(BasalforcingsIsmip6ExplicitOceanSalinityInputEnum); _assert_(sal_input);
		
		/*We must serve the DatasetInput manually as it is not part of the standard element inputs served in Element::Serve*/
		int* vertexlids = xNew<int>(numvertices);
		element->GetVerticesLidList(vertexlids);
		temp_input->Serve(numvertices,vertexlids);
		sal_input->Serve(numvertices,vertexlids);
		xDelete<int>(vertexlids);
		
		element->GetInputListOnVertices(&depth_vertices[0],BaseEnum);

		Gauss* gauss=element->NewGauss();
		for(int iv=0;iv<numvertices;iv++){
			gauss->GaussVertex(iv);

			/*Find out where the ice shelf base is within tf_depths*/
			IssmDouble depth = -depth_vertices[iv]; /*NOTE: make sure we are dealing with depth>0*/
			int offset;
			if(xIsNan<IssmDouble>(depth)) _error_("depth is NaN");
			int found=binary_search(&offset,depth,tf_depths_local,num_depths);
			if(!found) _error_("depth not found");

			IssmDouble temp, sal;
			int idx1, idx2;

			if(decreasing){
				if (offset==-1){
					/*get values for the first depth: (originally last)*/
					_assert_(depth<=tf_depths_local[0]);
					idx1 = num_depths - 1 - 0;
					temp_input->GetInputValue(&temp,gauss,idx1);
					sal_input->GetInputValue(&sal,gauss,idx1);
				}
				else if(offset==num_depths-1){
					/*get values for the last time: */
					_assert_(depth>=tf_depths_local[num_depths-1]);
					idx1 = num_depths - 1 - (num_depths-1);
					temp_input->GetInputValue(&temp,gauss,idx1);
					sal_input->GetInputValue(&sal,gauss,idx1);
				}
				else {
					/*get values between two times [offset:offset+1], Interpolate linearly*/
					_assert_(depth>=tf_depths_local[offset] && depth<tf_depths_local[offset+1]);
					IssmDouble deltaz=tf_depths_local[offset+1]-tf_depths_local[offset];
					IssmDouble alpha2=(depth-tf_depths_local[offset])/deltaz;
					IssmDouble alpha1=(1.-alpha2);
					IssmDouble temp1, temp2, sal1, sal2;

					idx1 = num_depths - 1 - offset;
					idx2 = num_depths - 1 - (offset+1);
					
					temp_input->GetInputValue(&temp1,gauss,idx1);
					temp_input->GetInputValue(&temp2,gauss,idx2);
					temp = alpha1*temp1 + alpha2*temp2;

					sal_input->GetInputValue(&sal1,gauss,idx1);
					sal_input->GetInputValue(&sal2,gauss,idx2);
					sal = alpha1*sal1 + alpha2*sal2;
				}
			}
			else{
				if (offset==-1){
					/*get values for the first depth: */
					_assert_(depth<=tf_depths_local[0]);
					temp_input->GetInputValue(&temp,gauss,0);
					sal_input->GetInputValue(&sal,gauss,0);
				}
				else if(offset==num_depths-1){
					/*get values for the last time: */
					_assert_(depth>=tf_depths_local[num_depths-1]);
					temp_input->GetInputValue(&temp,gauss,num_depths-1);
					sal_input->GetInputValue(&sal,gauss,num_depths-1);
				}
				else {
					/*get values between two times [offset:offset+1], Interpolate linearly*/
					_assert_(depth>=tf_depths_local[offset] && depth<tf_depths_local[offset+1]);
					IssmDouble deltaz=tf_depths_local[offset+1]-tf_depths_local[offset];
					IssmDouble alpha2=(depth-tf_depths_local[offset])/deltaz;
					IssmDouble alpha1=(1.-alpha2);
					IssmDouble temp1, temp2, sal1, sal2;
					
					temp_input->GetInputValue(&temp1,gauss,offset);
					temp_input->GetInputValue(&temp2,gauss,offset+1);
					temp = alpha1*temp1 + alpha2*temp2;

					sal_input->GetInputValue(&sal1,gauss,offset);
					sal_input->GetInputValue(&sal2,gauss,offset+1);
					sal = alpha1*sal1 + alpha2*sal2;
				}
			}
			
			if(xIsNan<IssmDouble>(temp) || xIsNan<IssmDouble>(sal)){
				int my_rank = IssmComm::GetRank();
				_printf_("Rank " << my_rank << " Element " << element->Id() << " Node " << iv << ": temp or sal is NaN.\n");
				_printf_("depth=" << depth << " offset=" << offset << " idx1=" << idx1 << " idx2=" << idx2 << "\n");
				_printf_("temp=" << temp << " sal=" << sal << "\n");
				_error_("NaN in input data");
			}

			/* Calculate TF = T - (lambda1*S + lambda2 + lambda3*z) */
			/* depth_vertices is usually negative (elevation). */
			tf_test[iv] = temp - (lambda1 * sal + lambda2 + lambda3 * depth_vertices[iv] + 273.0);
			
			if(xIsNan<IssmDouble>(tf_test[iv])){
				int my_rank = IssmComm::GetRank();
				_printf_("Rank " << my_rank << " Element " << element->Id() << " Node " << iv << ": tf_test is NaN.\n");
				_printf_("temp=" << temp << " sal=" << sal << " depth_vertices=" << depth_vertices[iv] << "\n");
				_error_("NaN in tf_test");
			}
		}

		element->AddInput(BasalforcingsIsmip6TfShelfEnum,tf_test,P1DGEnum);
		xDelete<IssmDouble>(tf_test);
		xDelete<IssmDouble>(depth_vertices);
		delete gauss;
	}

	/*Compute sums of tf*area and shelf-area per cpu*/
	if(!islocal) {
		for(Object* & object : femmodel->elements->objects){
			Element* element = xDynamicCast<Element*>(object);
			if(!element->IsOnBase() || !element->IsIceInElement() || !element->IsAllFloating()) continue;
			/*Spawn basal element if on base to compute element area*/
			Element* basalelement = element->SpawnBasalElement();
			Input* tf_input=basalelement->GetInput(BasalforcingsIsmip6TfShelfEnum); _assert_(tf_input);
			basalelement->GetInputValue(&basinid,BasalforcingsIsmip6BasinIdEnum);
			Gauss* gauss=basalelement->NewGauss(1); gauss->GaussPoint(0);
			tf_input->GetInputValue(&tf,gauss);
			delete gauss;
			area=basalelement->GetHorizontalSurfaceArea();
			tf_weighted_avg[basinid]+=tf*area;
			areas_summed[basinid]   +=area;
			/*Delete spawned element if we are in 3D*/
			basalelement->FindParam(&domaintype,DomainTypeEnum);
			if(basalelement->IsSpawnedElement()){basalelement->DeleteMaterials(); delete basalelement;};
		}

		/*Syncronize across cpus*/
		ISSM_MPI_Allreduce(tf_weighted_avg,tf_weighted_avg_cpu,num_basins,ISSM_MPI_DOUBLE,ISSM_MPI_SUM,IssmComm::GetComm());
		ISSM_MPI_Allreduce(areas_summed,areas_summed_cpu,num_basins,ISSM_MPI_DOUBLE,ISSM_MPI_SUM,IssmComm::GetComm());

		/*Make sure Area is not zero to avoid dividing by 0 if a basin is not present in the model*/
		for(int k=0;k<num_basins;k++) if(areas_summed_cpu[k]==0.) areas_summed_cpu[k] = 1.;

		/*Compute weighted means and save*/
		for(int k=0;k<num_basins;k++){tf_weighted_avg_cpu[k] = tf_weighted_avg_cpu[k]/areas_summed_cpu[k];}
		femmodel->parameters->AddObject(new DoubleVecParam(BasalforcingsIsmip6AverageTfEnum,tf_weighted_avg_cpu,num_basins));
	}

   /*Compute meltrates*/
	for(Object* & object : femmodel->elements->objects){
      Element* element = xDynamicCast<Element*>(object);
		element->Ismip6ExplicitFloatingiceMeltingRate();
		
		/*Check for NaN in meltrates? Ismip6ExplicitFloatingiceMeltingRate is void and adds input using AddInput*/
		/*We can't easily check inside the element function call from here without modifying Element::Ismip6ExplicitFloatingiceMeltingRate*/
	}
	/*Wait, I can modify Element::Ismip6ExplicitFloatingiceMeltingRate in Element.cpp code if needed. */
	/*But for now let's focus on inputs in this file.*/

	/*Cleanup and return */
	xDelete<IssmDouble>(tf_weighted_avg);
	xDelete<IssmDouble>(tf_weighted_avg_cpu);
	xDelete<IssmDouble>(areas_summed);
	xDelete<IssmDouble>(areas_summed_cpu);
	xDelete<IssmDouble>(tf_depths_local);
}
/*}}}*/
void BeckmannGoosseFloatingiceMeltingRatex(FemModel* femmodel){/*{{{*/

	for(Object* & object : femmodel->elements->objects){
      Element* element = xDynamicCast<Element*>(object);
		element->BeckmannGoosseFloatingiceMeltingRate();
	}
}
/*}}}*/
void LinearFloatingiceMeltingRatearmax(FemModel* femmodel){/*{{{*/

	/*Get time parameters*/
   IssmDouble time,dt,starttime,tstep_arma;
   femmodel->parameters->FindParam(&time,TimeEnum);
   femmodel->parameters->FindParam(&dt,TimesteppingTimeStepEnum);
   femmodel->parameters->FindParam(&starttime,TimesteppingStartTimeEnum);
   femmodel->parameters->FindParam(&tstep_arma,BasalforcingsARMATimestepEnum);

   /*Determine if this is a time step for the ARMA model*/
   bool isstepforarma = false;

   #ifndef _HAVE_AD_
   if((fmod(time,tstep_arma)<fmod((time-dt),tstep_arma)) || (time<=starttime+dt) || tstep_arma==dt) isstepforarma = true;
   #else
   _error_("not implemented yet");
   #endif

   /*Load parameters*/
   bool isstochastic;
   bool isdeepmeltingstochastic = false;
   int M,N,arorder,maorder,numbasins,numparams,numbreaks,my_rank;
   femmodel->parameters->FindParam(&numbasins,BasalforcingsLinearNumBasinsEnum);
	femmodel->parameters->FindParam(&arorder,BasalforcingsARMAarOrderEnum);
	femmodel->parameters->FindParam(&maorder,BasalforcingsARMAmaOrderEnum);
	femmodel->parameters->FindParam(&numparams,BasalforcingsLinearNumParamsEnum);
   femmodel->parameters->FindParam(&numbreaks,BasalforcingsLinearNumBreaksEnum);
   IssmDouble* datebreaks     = NULL;
	IssmDouble* arlagcoefs     = NULL;
   IssmDouble* malagcoefs     = NULL;
   IssmDouble* polyparams     = NULL;
	IssmDouble* deepwaterel    = NULL;
   IssmDouble* upperwaterel   = NULL;
   IssmDouble* upperwatermelt = NULL;
   IssmDouble* perturbation   = NULL;

	/*Get autoregressive parameters*/
   femmodel->parameters->FindParam(&datebreaks,&M,&N,BasalforcingsARMAdatebreaksEnum);  _assert_(M==numbasins); _assert_(N==max(numbreaks,1));
   femmodel->parameters->FindParam(&polyparams,&M,&N,BasalforcingsARMApolyparamsEnum);  _assert_(M==numbasins); _assert_(N==(numbreaks+1)*numparams);
	femmodel->parameters->FindParam(&arlagcoefs,&M,&N,BasalforcingsARMAarlagcoefsEnum);  _assert_(M==numbasins); _assert_(N==arorder);
   femmodel->parameters->FindParam(&malagcoefs,&M,&N,BasalforcingsARMAmalagcoefsEnum);  _assert_(M==numbasins); _assert_(N==maorder);

	/*Get basin-specific parameters*/
   femmodel->parameters->FindParam(&deepwaterel,&M,BasalforcingsDeepwaterElevationEnum);            _assert_(M==numbasins);
   femmodel->parameters->FindParam(&upperwaterel,&M,BasalforcingsUpperwaterElevationEnum);          _assert_(M==numbasins);
   femmodel->parameters->FindParam(&upperwatermelt,&M,BasalforcingsUpperwaterMeltingRateEnum);      _assert_(M==numbasins);

	/*Evaluate whether stochasticity on DeepwaterMeltingRate is requested*/
	femmodel->parameters->FindParam(&isstochastic,StochasticForcingIsStochasticForcingEnum);
   if(isstochastic){
      int  numstochasticfields;
      int* stochasticfields;
      femmodel->parameters->FindParam(&numstochasticfields,StochasticForcingNumFieldsEnum);
      femmodel->parameters->FindParam(&stochasticfields,&N,StochasticForcingFieldsEnum); _assert_(N==numstochasticfields);
      for(int i=0;i<numstochasticfields;i++){
         if(stochasticfields[i]==BasalforcingsDeepwaterMeltingRatearmaEnum) isdeepmeltingstochastic = true;
      }
      xDelete<int>(stochasticfields);
   }

	/*Loop over each element to compute FloatingiceMeltingRate at vertices*/
   for(Object* &object:femmodel->elements->objects){
      Element* element = xDynamicCast<Element*>(object);
      /*Compute ARMA*/
      element->ArmaProcess(isstepforarma,arorder,maorder,numparams,numbreaks,tstep_arma,polyparams,arlagcoefs,malagcoefs,datebreaks,isdeepmeltingstochastic,BasalforcingsDeepwaterMeltingRatearmaEnum);
		element->BasinLinearFloatingiceMeltingRate(deepwaterel,upperwatermelt,upperwaterel,perturbation);
	}

	/*Cleanup*/
	xDelete<IssmDouble>(arlagcoefs);
	xDelete<IssmDouble>(malagcoefs);
	xDelete<IssmDouble>(polyparams);
	xDelete<IssmDouble>(datebreaks);
	xDelete<IssmDouble>(deepwaterel);
	xDelete<IssmDouble>(upperwaterel);
	xDelete<IssmDouble>(upperwatermelt);
	xDelete<IssmDouble>(perturbation);
}/*}}}*/
