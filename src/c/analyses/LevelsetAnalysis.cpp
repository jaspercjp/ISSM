#ifdef HAVE_CONFIG_H
#include <config.h>
#else
#error "Cannot compile with HAVE_CONFIG_H symbol! run configure first!"
#endif
#include "../classes/Inputs/TriaInput.h"
#include "../classes/classes.h"
#include "../modules/modules.h"
#include "../shared/shared.h"
#include "../solutionsequences/solutionsequences.h"
#include "../toolkits/toolkits.h"
#include "./LevelsetAnalysis.h"
#include <algorithm>
#include <limits>
#include <math.h>
#include <vector>

namespace {

static const IssmDouble kCalvedLevelsetValue = 400.0;

struct CrevasseCalvingConfig {
  IssmDouble time_yr;
  IssmDouble dt_yr;
  IssmDouble timescale_yr;
  IssmDouble min_iceberg_size;
  IssmDouble critical_stress;
  IssmDouble hydrofracture_min_ocean_levelset;
  bool advect_icefront;
  bool propagate_from_front;
  bool remove_only_marked;
  bool hydrofracture_stabilization;
  IssmDouble hydrofracture_weakening_factor;
};

struct CrevasseNodalFields {
  IssmDouble *averaged_principal_stress_1;
  IssmDouble *averaged_principal_stress_2;
  IssmDouble *averaged_thickness;
  IssmDouble *averaged_groundingline_distance;
  IssmDouble *floating_nodes;
  IssmDouble *seed_nodes;
  IssmDouble *meltwater_mask;
  IssmDouble *node_x;
  IssmDouble ice_front_x;

  CrevasseNodalFields()
      : averaged_principal_stress_1(NULL), averaged_principal_stress_2(NULL),
        averaged_thickness(NULL), averaged_groundingline_distance(NULL),
        floating_nodes(NULL), seed_nodes(NULL), meltwater_mask(NULL),
        node_x(NULL),
        ice_front_x(std::numeric_limits<IssmDouble>::max()) {}
};

struct CrevassePropagationResult {
  IssmDouble *propagated_nodes;
  IssmDouble *applied_nodes;
  IssmDouble propagated_min_x;
  IssmDouble propagated_area;
  bool has_critical_region;

  CrevassePropagationResult()
      : propagated_nodes(NULL), applied_nodes(NULL),
        propagated_min_x(std::numeric_limits<IssmDouble>::max()),
        propagated_area(0.0), has_critical_region(false) {}
};

struct CrevasseStressSummary {
  int count;
  IssmDouble ratio_sum;
  IssmDouble ratio_min;
  IssmDouble ratio_max;
  IssmDouble threshold_sum;
  IssmDouble threshold_min;
  IssmDouble threshold_max;

  CrevasseStressSummary()
      : count(0), ratio_sum(0.0),
        ratio_min(std::numeric_limits<IssmDouble>::max()),
        ratio_max(-std::numeric_limits<IssmDouble>::max()),
        threshold_sum(0.0),
        threshold_min(std::numeric_limits<IssmDouble>::max()),
        threshold_max(-std::numeric_limits<IssmDouble>::max()) {}
};

struct CrevasseStressDiagnostics {
  CrevasseStressSummary floating;
  CrevasseStressSummary seed;
  CrevasseStressSummary propagated;
};

IssmDouble ComputeLEFMRxx(IssmDouble averaged_stress_1,
                          IssmDouble averaged_stress_2);

/* Return true when the node already has a user-prescribed level-set SPC. */
bool HasStaticLevelsetSpc(Element *element, int node_index, Gauss *gauss) {
  Input *spc_input = element->GetInput(SpcLevelsetEnum);
  if (!spc_input)
    return false;

  gauss->GaussVertex(node_index);
  IssmDouble spc_value;
  spc_input->GetInputValue(&spc_value, gauss);
  return !xIsNan<IssmDouble>(spc_value);
}

/* Count unique master nodes flagged for calving across all MPI ranks. */
int CountMarkedMasterNodes(FemModel *femmodel, const IssmDouble *marked_nodes) {
  int local_marked_nodes = 0;

  for (Object *&object : femmodel->nodes->objects) {
    Node *node = xDynamicCast<Node *>(object);
    if (node->IsClone())
      continue;
    if (marked_nodes[node->Lid()] > 0.5)
      local_marked_nodes++;
  }

  int global_marked_nodes = 0;
  ISSM_MPI_Allreduce(&local_marked_nodes, &global_marked_nodes, 1, ISSM_MPI_INT,
                     ISSM_MPI_SUM, IssmComm::GetComm());
  return global_marked_nodes;
}

void UpdateCrevasseStressSummary(CrevasseStressSummary *summary,
                                 IssmDouble r_xx,
                                 IssmDouble stress_threshold) {
  summary->count++;
  summary->ratio_sum += r_xx;
  summary->ratio_min = fmin(summary->ratio_min, r_xx);
  summary->ratio_max = fmax(summary->ratio_max, r_xx);
  summary->threshold_sum += stress_threshold;
  summary->threshold_min = fmin(summary->threshold_min, stress_threshold);
  summary->threshold_max = fmax(summary->threshold_max, stress_threshold);
}

void ReduceCrevasseStressSummary(CrevasseStressSummary *summary) {
  int global_count = 0;
  IssmDouble global_ratio_sum = 0.0;
  IssmDouble global_ratio_min = summary->ratio_min;
  IssmDouble global_ratio_max = summary->ratio_max;
  IssmDouble global_threshold_sum = 0.0;
  IssmDouble global_threshold_min = summary->threshold_min;
  IssmDouble global_threshold_max = summary->threshold_max;

  ISSM_MPI_Allreduce(&summary->count, &global_count, 1, ISSM_MPI_INT,
                     ISSM_MPI_SUM, IssmComm::GetComm());
  ISSM_MPI_Allreduce(&summary->ratio_sum, &global_ratio_sum, 1,
                     ISSM_MPI_DOUBLE, ISSM_MPI_SUM, IssmComm::GetComm());
  ISSM_MPI_Allreduce(&summary->ratio_min, &global_ratio_min, 1,
                     ISSM_MPI_DOUBLE, ISSM_MPI_MIN, IssmComm::GetComm());
  ISSM_MPI_Allreduce(&summary->ratio_max, &global_ratio_max, 1,
                     ISSM_MPI_DOUBLE, ISSM_MPI_MAX, IssmComm::GetComm());
  ISSM_MPI_Allreduce(&summary->threshold_sum, &global_threshold_sum, 1,
                     ISSM_MPI_DOUBLE, ISSM_MPI_SUM, IssmComm::GetComm());
  ISSM_MPI_Allreduce(&summary->threshold_min, &global_threshold_min, 1,
                     ISSM_MPI_DOUBLE, ISSM_MPI_MIN, IssmComm::GetComm());
  ISSM_MPI_Allreduce(&summary->threshold_max, &global_threshold_max, 1,
                     ISSM_MPI_DOUBLE, ISSM_MPI_MAX, IssmComm::GetComm());

  summary->count = global_count;
  summary->ratio_sum = global_ratio_sum;
  summary->ratio_min = global_ratio_min;
  summary->ratio_max = global_ratio_max;
  summary->threshold_sum = global_threshold_sum;
  summary->threshold_min = global_threshold_min;
  summary->threshold_max = global_threshold_max;
}

CrevasseStressDiagnostics ComputeCrevasseStressDiagnostics(
    FemModel *femmodel, const CrevasseNodalFields &fields,
    const CrevassePropagationResult &propagation) {
  CrevasseStressDiagnostics diagnostics;
  IssmDouble critical_stress;
  femmodel->parameters->FindParam(&critical_stress, CalvingCriticalStressEnum);

  for (Object *&object : femmodel->nodes->objects) {
    Node *node = xDynamicCast<Node *>(object);
    if (node->IsClone())
      continue;

    int lid = node->Lid();
    if (fields.floating_nodes[lid] <= 0.5)
      continue;

    IssmDouble r_xx = ComputeLEFMRxx(
        fields.averaged_principal_stress_1[lid],
        fields.averaged_principal_stress_2[lid]);

    UpdateCrevasseStressSummary(&diagnostics.floating, r_xx,
                                critical_stress);
    if (fields.seed_nodes[lid] > 0.5)
      UpdateCrevasseStressSummary(&diagnostics.seed, r_xx,
                                  critical_stress);
    if (propagation.propagated_nodes && propagation.propagated_nodes[lid] > 0.5)
      UpdateCrevasseStressSummary(&diagnostics.propagated, r_xx,
                                  critical_stress);
  }

  ReduceCrevasseStressSummary(&diagnostics.floating);
  ReduceCrevasseStressSummary(&diagnostics.seed);
  ReduceCrevasseStressSummary(&diagnostics.propagated);
  return diagnostics;
}

void PrintCrevasseStressSummary(const char *label,
                                const CrevasseStressSummary &summary) {
  if (summary.count == 0) {
    _printf0_("\t" << label << " r_xx/critical stress = n/a\n");
    return;
  }

  _printf0_("\t" << label << " r_xx[min/mean/max] = " << summary.ratio_min
            << ", " << summary.ratio_sum / reCast<IssmDouble>(summary.count)
            << ", " << summary.ratio_max << "\n");
  _printf0_("\t" << label << " crit.[min/mean/max] = "
            << summary.threshold_min << ", "
            << summary.threshold_sum / reCast<IssmDouble>(summary.count)
            << ", " << summary.threshold_max << "\n");
}

/* Build a nodal calving mask for the floating region seaward of the inland
 * calving cutoff. This removes the entire disconnected seaward patch once a
 * valid calving position has been identified, rather than only the subset of
 * nodes that happened to be visited during the stress-propagation flood fill. */
void BuildAppliedCalvingMask(FemModel *femmodel,
                             const CrevasseNodalFields &fields,
                             CrevassePropagationResult *result) {
  int numnodes = femmodel->nodes->NumberOfNodes();
  int localmasters = femmodel->nodes->NumberOfNodesLocal();

  Vector<IssmDouble> *vec_applied =
      new Vector<IssmDouble>(localmasters, numnodes);

  if (result->propagated_min_x == std::numeric_limits<IssmDouble>::max()) {
    vec_applied->Assemble();
    xDelete<IssmDouble>(result->applied_nodes);
    femmodel->GetLocalVectorWithClonesNodes(&result->applied_nodes, vec_applied);
    delete vec_applied;
    return;
  }

  for (Object *&object : femmodel->nodes->objects) {
    Node *node = xDynamicCast<Node *>(object);
    if (node->IsClone())
      continue;
    int lid = node->Lid();
    if (fields.floating_nodes[lid] > 0.5 &&
        fields.node_x[lid] >= result->propagated_min_x) {
      vec_applied->SetValue(node->Pid(), 1.0, INS_VAL);
    }
  }

  vec_applied->Assemble();
  xDelete<IssmDouble>(result->applied_nodes);
  femmodel->GetLocalVectorWithClonesNodes(&result->applied_nodes, vec_applied);
  delete vec_applied;
}

/* Remove transient calving SPCs while preserving any static user constraints. */
void ClearDynamicLevelsetConstraintsPreservingStaticSpc(FemModel *femmodel) {
  for (Object *&object : femmodel->elements->objects) {
    Element *element = xDynamicCast<Element *>(object);
    int numnodes = element->GetNumberOfNodes();
    Gauss *gauss = element->NewGauss();

    for (int in = 0; in < numnodes; in++) {
      Node *node = element->GetNode(in);
      if (!node->IsActive())
        continue;
      if (!HasStaticLevelsetSpc(element, in, gauss))
        node->DofInFSet(0);
    }
    delete gauss;
  }
}

/* Freeze the present ice front by constraining front nodes to their current level-set values. */
void FreezeIceFrontAdvectionIfNeeded(FemModel *femmodel) {
  for (Object *&object : femmodel->elements->objects) {
    Element *element = xDynamicCast<Element *>(object);
    if (!element->IsIcefront())
      continue;

    int numnodes = element->GetNumberOfNodes();
    Input *levelset_input = element->GetInput(MaskIceLevelsetEnum);
    _assert_(levelset_input);
    Gauss *gauss = element->NewGauss();

    for (int in = 0; in < numnodes; in++) {
      gauss->GaussVertex(in);
      Node *node = element->GetNode(in);
      if (!node->IsActive())
        continue;

      IssmDouble current_levelset;
      levelset_input->GetInputValue(&current_levelset, gauss);
      node->ApplyConstraint(0, current_levelset);
    }
    delete gauss;
  }
}

/* Compute the nondimensional crevasse stress ratio used for propagation checks in the HFB framework. */
IssmDouble ComputeHFBCrevasseStressRatio(IssmDouble averaged_s_xx,
                                      IssmDouble averaged_s_yy,
                                      IssmDouble averaged_thickness,
                                      IssmDouble rho_ice,
                                      IssmDouble rho_seawater,
                                      IssmDouble constant_g) {
  IssmDouble r_it = rho_ice * constant_g * averaged_thickness *
                    (1.0 - rho_ice / rho_seawater) / 2.0;
  if (r_it <= 0.0)
    return -std::numeric_limits<IssmDouble>::max();

  IssmDouble r_xx = 2.0 * averaged_s_xx + averaged_s_yy;
  return fmin(r_xx / r_it, 1.5);
}

/* Absolute LEFM opening stress used for hydrofracture checks [Pa]. */
IssmDouble ComputeLEFMRxx(IssmDouble averaged_stress_1,
                          IssmDouble averaged_stress_2) {
  return 2.0 * averaged_stress_1 + averaged_stress_2;
}

void ComputePrincipalStresses2D(IssmDouble *pstress_1, IssmDouble *pstress_2,
                                IssmDouble s_xx, IssmDouble s_yy,
                                IssmDouble s_xy) {
  IssmDouble mean_stress = 0.5 * (s_xx + s_yy);
  IssmDouble radius =
      sqrt(0.25 * (s_xx - s_yy) * (s_xx - s_yy) + s_xy * s_xy);

  *pstress_1 = mean_stress + radius;
  *pstress_2 = mean_stress - radius;
}

bool IsCriticalMeltwaterNode(const CrevasseCalvingConfig &config,
                             const CrevasseNodalFields &fields, int lid) {
  if (fields.floating_nodes[lid] <= 0.5)
    return false;
  if (fields.meltwater_mask[lid] < 0.5)
    return false;
  if (fields.averaged_groundingline_distance[lid] >
      config.hydrofracture_min_ocean_levelset)
    return false;

  IssmDouble r_xx = ComputeLEFMRxx(
      fields.averaged_principal_stress_1[lid],
      fields.averaged_principal_stress_2[lid]);

  return r_xx >= config.critical_stress;
}

/* Store the instantaneous hydrofracture prediction as a P1 input so users can
 * request it as an output. The prediction is deliberately the local criterion
 * only: floating ice, meltwater_mask == 1, and r_xx >= critical_stress.
 * It does not include front-connectivity propagation or minimum-iceberg-size
 * filtering, so it exposes the raw hydrofracture susceptibility mask. */
void StoreHydrofracturePredictionInput(FemModel *femmodel,
                                       const CrevasseCalvingConfig &config,
                                       const CrevasseNodalFields &fields) {
  int numvertices = femmodel->vertices->NumberOfVertices();
  int numvertices_local = femmodel->vertices->NumberOfVerticesLocal();
  Vector<IssmDouble> *prediction =
      new Vector<IssmDouble>(numvertices_local, numvertices);

  _printf0_("\texecuting hydrofracture prediction check\n");

  int local_predicted_nodes = 0;
  for (Object *&object : femmodel->nodes->objects) {
    Node *node = xDynamicCast<Node *>(object);
    if (node->IsClone())
      continue;

    if (IsCriticalMeltwaterNode(config, fields, node->Lid()))
      local_predicted_nodes++;
  }

  int global_predicted_nodes = 0;
  ISSM_MPI_Allreduce(&local_predicted_nodes, &global_predicted_nodes, 1,
                     ISSM_MPI_INT, ISSM_MPI_SUM, IssmComm::GetComm());
  _printf0_("\tchecking hydrofracture criterion: "
            << global_predicted_nodes
            << " floating meltwater nodes satisfy r_xx >= "
            << config.critical_stress
            << " Pa and signed grounding-line distance <= "
            << config.hydrofracture_min_ocean_levelset << " m\n");

  for (Object *&object : femmodel->elements->objects) {
    Element *element = xDynamicCast<Element *>(object);
    int numnodes = element->GetNumberOfNodes();

    for (int in = 0; in < numnodes; in++) {
      Node *node = element->GetNode(in);
      if (!node->IsActive())
        continue;

      IssmDouble predicted =
          IsCriticalMeltwaterNode(config, fields, node->Lid())
              ? 1.0
              : 0.0;
      prediction->SetValue(element->vertices[in]->Pid(), predicted, INS_VAL);
    }
  }

  prediction->Assemble();
  InputUpdateFromVectorx(femmodel, prediction, HydrofracturePredictedEnum,
                         VertexPIdEnum);
  delete prediction;
}

/* Store the regularized weak-ice mask. This mask is consumed by material
 * routines as a local multiplier on rheology B, which keeps hydrofractured
 * regions mechanically weak without creating interior no-ice holes in the
 * active mesh. */
void StoreHydrofractureWeakIceInput(FemModel *femmodel,
                                    const IssmDouble *weak_nodes) {
  int numvertices = femmodel->vertices->NumberOfVertices();
  int numvertices_local = femmodel->vertices->NumberOfVerticesLocal();
  Vector<IssmDouble> *weak =
      new Vector<IssmDouble>(numvertices_local, numvertices);

  int local_weak_nodes = 0;
  for (Object *&object : femmodel->nodes->objects) {
    Node *node = xDynamicCast<Node *>(object);
    if (node->IsClone())
      continue;
    if (weak_nodes && weak_nodes[node->Lid()] > 0.5)
      local_weak_nodes++;
  }

  int global_weak_nodes = 0;
  ISSM_MPI_Allreduce(&local_weak_nodes, &global_weak_nodes, 1, ISSM_MPI_INT,
                     ISSM_MPI_SUM, IssmComm::GetComm());
  _printf0_("\thydrofracture weak-ice nodes = " << global_weak_nodes << "\n");

  for (Object *&object : femmodel->elements->objects) {
    Element *element = xDynamicCast<Element *>(object);
    int numnodes = element->GetNumberOfNodes();

    for (int in = 0; in < numnodes; in++) {
      Node *node = element->GetNode(in);
      if (!node->IsActive())
        continue;

      IssmDouble weak_value =
          (weak_nodes && weak_nodes[node->Lid()] > 0.5) ? 1.0 : 0.0;
      weak->SetValue(element->vertices[in]->Pid(), weak_value, INS_VAL);
    }
  }

  weak->Assemble();
  InputUpdateFromVectorx(femmodel, weak, HydrofractureWeakIceEnum,
                         VertexPIdEnum);
  delete weak;
}

/* Release temporary nodal arrays assembled for crevasse propagation. */
void CleanupCrevasseNodalFields(CrevasseNodalFields *fields) {
  xDelete<IssmDouble>(fields->averaged_principal_stress_1);
  xDelete<IssmDouble>(fields->averaged_principal_stress_2);
  xDelete<IssmDouble>(fields->averaged_thickness);
  xDelete<IssmDouble>(fields->averaged_groundingline_distance);
  xDelete<IssmDouble>(fields->floating_nodes);
  xDelete<IssmDouble>(fields->seed_nodes);
  xDelete<IssmDouble>(fields->meltwater_mask);
  xDelete<IssmDouble>(fields->node_x);
}

/* Recompute a signed distance-to-grounding-line field from the current ocean
 * level-set sign. MaskOceanLevelset is updated by grounding-line migration as
 * a flotation criterion, so its magnitude is not a persistent distance in
 * meters. The hydrofracture buffer must therefore use this derived distance,
 * not the raw MaskOceanLevelset value. */
void PrepareHydrofractureGroundinglineDistance(FemModel *femmodel) {
  InputDuplicatex(femmodel, MaskOceanLevelsetEnum, DistanceToGroundinglineEnum);
  femmodel->DistanceToFieldValue(MaskOceanLevelsetEnum, 0.,
                                 DistanceToGroundinglineEnum);
}

/* Assemble nodal stress, thickness, flotation, grounding-line distance, and
 * seed-front fields for calving propagation. */
void BuildCrevasseNodalFields(FemModel *femmodel, CrevasseNodalFields *fields) {
  int numnodes = femmodel->nodes->NumberOfNodes();
  int localmasters = femmodel->nodes->NumberOfNodesLocal();
  int localsize = femmodel->nodes->NumberOfNodesLocalAll();

  Vector<IssmDouble> *vec_principal_stress_1 =
      new Vector<IssmDouble>(localmasters, numnodes);
  Vector<IssmDouble> *vec_principal_stress_2 =
      new Vector<IssmDouble>(localmasters, numnodes);
  Vector<IssmDouble> *vec_thickness =
      new Vector<IssmDouble>(localmasters, numnodes);
  Vector<IssmDouble> *vec_groundingline_distance =
      new Vector<IssmDouble>(localmasters, numnodes);
  Vector<IssmDouble> *vec_count =
      new Vector<IssmDouble>(localmasters, numnodes);
  Vector<IssmDouble> *vec_floating =
      new Vector<IssmDouble>(localmasters, numnodes);
  Vector<IssmDouble> *vec_seed = new Vector<IssmDouble>(localmasters, numnodes);
  Vector<IssmDouble> *vec_meltwater_mask =
      new Vector<IssmDouble>(localmasters, numnodes);
  Vector<IssmDouble> *vec_x = new Vector<IssmDouble>(localmasters, numnodes);

  IssmDouble local_ice_front_x = std::numeric_limits<IssmDouble>::max();

  for (Object *&object : femmodel->elements->objects) {
    Element *element = xDynamicCast<Element *>(object);
    if (!element->IsIceInElement())
      continue;

    int numnodes_element = element->GetNumberOfNodes();
    Input *s_xx_input = element->GetInput(DeviatoricStressxxEnum);
    Input *s_yy_input = element->GetInput(DeviatoricStressyyEnum);
    Input *s_xy_input = element->GetInput(DeviatoricStressxyEnum);
    Input *thickness_input = element->GetInput(ThicknessEnum);
    Input *ice_levelset_input = element->GetInput(MaskIceLevelsetEnum);
    Input *ocean_levelset_input = element->GetInput(MaskOceanLevelsetEnum);
    Input *groundingline_distance_input =
        element->GetInput(DistanceToGroundinglineEnum);
    Input *meltwater_mask_input = element->GetInput(CalvingMeltwaterMaskEnum);
    _assert_(s_xx_input);
    _assert_(s_yy_input);
    _assert_(s_xy_input);
    _assert_(thickness_input);
    _assert_(ice_levelset_input);
    _assert_(ocean_levelset_input);
    _assert_(groundingline_distance_input);
    _assert_(meltwater_mask_input);

    Gauss *gauss = element->NewGauss();
    bool is_icefront = element->IsIcefront();

    for (int in = 0; in < numnodes_element; in++) {
      gauss->GaussVertex(in);
      Node *node = element->GetNode(in);
      if (!node->IsActive())
        continue;

      IssmDouble s_xx, s_yy, s_xy, thickness, ice_levelset, ocean_levelset,
                 groundingline_distance, meltwater_mask, principal_stress_1,
                 principal_stress_2;
      s_xx_input->GetInputValue(&s_xx, gauss);
      s_yy_input->GetInputValue(&s_yy, gauss);
      s_xy_input->GetInputValue(&s_xy, gauss);
      thickness_input->GetInputValue(&thickness, gauss);
      ice_levelset_input->GetInputValue(&ice_levelset, gauss);
      ocean_levelset_input->GetInputValue(&ocean_levelset, gauss);
      groundingline_distance_input->GetInputValue(&groundingline_distance,
                                                  gauss);
      meltwater_mask_input->GetInputValue(&meltwater_mask, gauss);

      ComputePrincipalStresses2D(&principal_stress_1, &principal_stress_2,
                                 s_xx, s_yy, s_xy);

      vec_principal_stress_1->SetValue(node->Pid(), principal_stress_1,
                                       ADD_VAL);
      vec_principal_stress_2->SetValue(node->Pid(), principal_stress_2,
                                       ADD_VAL);
      vec_thickness->SetValue(node->Pid(), thickness, ADD_VAL);
      vec_groundingline_distance->SetValue(node->Pid(),
                                           groundingline_distance, ADD_VAL);
      vec_meltwater_mask->SetValue(node->Pid(), meltwater_mask, ADD_VAL);
      vec_count->SetValue(node->Pid(), 1.0, ADD_VAL);
      vec_x->SetValue(node->Pid(), element->vertices[in]->x, INS_VAL);

      if (ice_levelset <= 0.0 && ocean_levelset < 0.0) {
        vec_floating->SetValue(node->Pid(), 1.0, INS_VAL);
        if (is_icefront) {
          vec_seed->SetValue(node->Pid(), 1.0, INS_VAL);
          local_ice_front_x = fmin(local_ice_front_x, element->vertices[in]->x);
        }
      }
    }
    delete gauss;
  }

  vec_principal_stress_1->Assemble();
  vec_principal_stress_2->Assemble();
  vec_thickness->Assemble();
  vec_groundingline_distance->Assemble();
  vec_count->Assemble();
  vec_floating->Assemble();
  vec_seed->Assemble();
  vec_meltwater_mask->Assemble();
  vec_x->Assemble();

  IssmDouble *local_count = NULL;
  femmodel->GetLocalVectorWithClonesNodes(&fields->averaged_principal_stress_1,
                                          vec_principal_stress_1);
  femmodel->GetLocalVectorWithClonesNodes(&fields->averaged_principal_stress_2,
                                          vec_principal_stress_2);
  femmodel->GetLocalVectorWithClonesNodes(&fields->averaged_thickness,
                                          vec_thickness);
  femmodel->GetLocalVectorWithClonesNodes(
      &fields->averaged_groundingline_distance, vec_groundingline_distance);
  femmodel->GetLocalVectorWithClonesNodes(&fields->floating_nodes, vec_floating);
  femmodel->GetLocalVectorWithClonesNodes(&fields->seed_nodes, vec_seed);
  femmodel->GetLocalVectorWithClonesNodes(&fields->meltwater_mask,
                                          vec_meltwater_mask);
  femmodel->GetLocalVectorWithClonesNodes(&fields->node_x, vec_x);
  femmodel->GetLocalVectorWithClonesNodes(&local_count, vec_count);

  for (int i = 0; i < localsize; i++) {
    if (local_count[i] > 0.0) {
      fields->averaged_principal_stress_1[i] /= local_count[i];
      fields->averaged_principal_stress_2[i] /= local_count[i];
      fields->averaged_thickness[i] /= local_count[i];
      fields->averaged_groundingline_distance[i] /= local_count[i];
      fields->meltwater_mask[i] /= local_count[i];
    }
  }

  ISSM_MPI_Allreduce(&local_ice_front_x, &fields->ice_front_x, 1,
                     ISSM_MPI_DOUBLE, ISSM_MPI_MIN, IssmComm::GetComm());

  xDelete<IssmDouble>(local_count);
  delete vec_principal_stress_1;
  delete vec_principal_stress_2;
  delete vec_thickness;
  delete vec_groundingline_distance;
  delete vec_count;
  delete vec_floating;
  delete vec_seed;
  delete vec_meltwater_mask;
  delete vec_x;
}

/* Flood-fill the floating critical region starting from the current ice-front seed nodes. */
void PropagateCriticalCrevasseRegion(FemModel *femmodel,
                                     const CrevasseCalvingConfig &config,
                                     const CrevasseNodalFields &fields,
                                     CrevassePropagationResult *result) {
  int numnodes = femmodel->nodes->NumberOfNodes();
  int localmasters = femmodel->nodes->NumberOfNodesLocal();
  int localsize = femmodel->nodes->NumberOfNodesLocalAll();
  _printf0_("\texecuting hydrofracture propagation checks from ice front\n");

  Vector<IssmDouble> *vec_propagated =
      new Vector<IssmDouble>(localmasters, numnodes);

  /* Seed the propagation with the present ice-front nodes. The region grown
   * below is therefore the connected floating patch that can be reached from
   * the current terminus while satisfying the crevasse criterion. */
  for (Object *&object : femmodel->elements->objects) {
    Element *element = xDynamicCast<Element *>(object);
    int numnodes_element = element->GetNumberOfNodes();
    for (int in = 0; in < numnodes_element; in++) {
      Node *node = element->GetNode(in);
      if (!node->IsActive())
        continue;
      if (fields.seed_nodes[node->Lid()] > 0.5)
        vec_propagated->SetValue(node->Pid(), 1.0, INS_VAL);
    }
  }

  vec_propagated->Assemble();
  femmodel->GetLocalVectorWithClonesNodes(&result->propagated_nodes,
                                          vec_propagated);

  int global_new = 1;
  while (global_new > 0) {
    /* Snapshot of the propagated mask at the beginning of this sweep.
     * previous[i] > 0.5 means node i was already part of the connected
     * critical region before any new nodes are added in the current pass.
     * previous[i] <= 0.5 means node i is still outside that region.
     *
     * Keeping this frozen copy prevents one sweep from chaining through
     * several neighbors at once; instead, the region grows one edge-neighbor
     * layer per iteration and then synchronizes across MPI ranks. */
    std::vector<IssmDouble> previous(result->propagated_nodes,
                                     result->propagated_nodes + localsize);

    /* Examine fully icy triangles and test whether the propagated region
     * should spread across any of their edges. */
    for (Object *&object : femmodel->elements->objects) {
      Element *element = xDynamicCast<Element *>(object);
      if (!element->IsIceInElement())
        continue;

      int numnodes_element = element->GetNumberOfNodes();
      _assert_(numnodes_element == 3);

      IssmDouble ls[3];
      element->GetInputListOnVertices(&ls[0], MaskIceLevelsetEnum);
      int nrice = 0;
      for (int i = 0; i < 3; i++) {
        if (ls[i] < 0.)
          nrice++;
      }
      if (nrice < 3)
        continue;

      static const int edge_pairs[3][2] = {{0, 1}, {1, 2}, {2, 0}};
      for (int edge = 0; edge < 3; edge++) {
        int from = edge_pairs[edge][0];
        int to = edge_pairs[edge][1];

        Node *from_node = element->GetNode(from);
        Node *to_node = element->GetNode(to);
        if (from_node->IsActive() && to_node->IsActive()) {
          int from_lid = from_node->Lid();
          int to_lid = to_node->Lid();

          /* Edge direction 1: grow from "from" to "to". This branch is taken
           * only when the from-node was already in the propagated region at the
           * start of the sweep and the to-node was not. If the to-node is
           * floating and its local stress state exceeds the critical stress, the
           * calving region spreads across this edge and marks the to-node. */
          if (previous[from_lid] > 0.5 && previous[to_lid] <= 0.5 &&
              IsCriticalMeltwaterNode(config, fields, to_lid)) {
              vec_propagated->SetValue(to_node->Pid(), 1.0, INS_VAL);
          }

          /* Edge direction 2: same test in the opposite orientation. The edge
           * list has an arbitrary local ordering, but the propagation itself is
           * undirected: if either endpoint was previously marked and the other
           * satisfies the flotation and stress checks, the region should spread
           * across that edge. The two branches simply cover both possible
           * orderings of the same edge. */
          if (previous[to_lid] > 0.5 && previous[from_lid] <= 0.5 &&
              IsCriticalMeltwaterNode(config, fields, from_lid)) {
              vec_propagated->SetValue(from_node->Pid(), 1.0, INS_VAL);
          }
        }
      }
    }

    vec_propagated->Assemble();
    xDelete<IssmDouble>(result->propagated_nodes);
    femmodel->GetLocalVectorWithClonesNodes(&result->propagated_nodes,
                                            vec_propagated);

    /* Stop once no rank contributes any newly marked nodes, meaning the
     * connected critical region has reached its final extent. */
    int local_new = 0;
    for (int i = 0; i < localsize; i++) {
      if (result->propagated_nodes[i] > 0.5 && previous[i] <= 0.5)
        local_new++;
    }

    ISSM_MPI_Allreduce(&local_new, &global_new, 1, ISSM_MPI_INT,
                       ISSM_MPI_SUM, IssmComm::GetComm());
    if (global_new > 0)
      result->has_critical_region = true;
  }

  delete vec_propagated;
}

/* Directly mark all floating nodes that satisfy the gated LEFM criterion. */
void MarkCriticalMeltwaterNodes(FemModel *femmodel,
                                const CrevasseCalvingConfig &config,
                                const CrevasseNodalFields &fields,
                                CrevassePropagationResult *result) {
  int numnodes = femmodel->nodes->NumberOfNodes();
  int localmasters = femmodel->nodes->NumberOfNodesLocal();
  int localsize = femmodel->nodes->NumberOfNodesLocalAll();

  Vector<IssmDouble> *vec_marked =
      new Vector<IssmDouble>(localmasters, numnodes);

  _printf0_("\texecuting hydrofracture checks on all floating ice\n");

  for (Object *&object : femmodel->nodes->objects) {
    Node *node = xDynamicCast<Node *>(object);
    if (node->IsClone())
      continue;

    int lid = node->Lid();
    if (IsCriticalMeltwaterNode(config, fields, lid)) {
      vec_marked->SetValue(node->Pid(), 1.0, INS_VAL);
    }
  }

  vec_marked->Assemble();
  xDelete<IssmDouble>(result->propagated_nodes);
  femmodel->GetLocalVectorWithClonesNodes(&result->propagated_nodes,
                                          vec_marked);

  xDelete<IssmDouble>(result->applied_nodes);
  result->applied_nodes = xNew<IssmDouble>(localsize);
  for (int i = 0; i < localsize; i++)
    result->applied_nodes[i] = result->propagated_nodes[i];

  result->has_critical_region =
      CountMarkedMasterNodes(femmodel, result->propagated_nodes) > 0;

  delete vec_marked;
}

/* Diagnose the inland extent and area of the fully propagated floating calving region. */
void ComputePropagatedRegionDiagnostics(FemModel *femmodel,
                                        const CrevasseNodalFields &fields,
                                        CrevassePropagationResult *result) {
  IssmDouble local_min_x = std::numeric_limits<IssmDouble>::max();
  IssmDouble local_area = 0.0;

  for (Object *&object : femmodel->elements->objects) {
    Element *element = xDynamicCast<Element *>(object);
    int numnodes_element = element->GetNumberOfNodes();
    _assert_(numnodes_element == 3);

    bool all_propagated = true;
    bool all_floating = true;
    for (int in = 0; in < numnodes_element; in++) {
      Node *node = element->GetNode(in);
      if (!node->IsActive()) {
        all_propagated = false;
        all_floating = false;
        continue;
      }

      int lid = node->Lid();
      if (result->propagated_nodes[lid] <= 0.5)
        all_propagated = false;
      if (fields.floating_nodes[lid] <= 0.5)
        all_floating = false;
      if (result->propagated_nodes[lid] > 0.5)
        local_min_x = fmin(local_min_x, fields.node_x[lid]);
    }

    if (all_propagated && all_floating) {
      IssmDouble x1 = element->vertices[0]->x;
      IssmDouble y1 = element->vertices[0]->y;
      IssmDouble x2 = element->vertices[1]->x;
      IssmDouble y2 = element->vertices[1]->y;
      IssmDouble x3 = element->vertices[2]->x;
      IssmDouble y3 = element->vertices[2]->y;
      local_area +=
          0.5 * fabs(x1 * (y2 - y3) + x2 * (y3 - y1) + x3 * (y1 - y2));
    }
  }

  ISSM_MPI_Allreduce(&local_min_x, &result->propagated_min_x, 1,
                     ISSM_MPI_DOUBLE, ISSM_MPI_MIN, IssmComm::GetComm());
  ISSM_MPI_Allreduce(&local_area, &result->propagated_area, 1,
                     ISSM_MPI_DOUBLE, ISSM_MPI_SUM, IssmComm::GetComm());
}

/* Apply level-set SPCs that remove the propagated calving region during the current solve. */
void ApplyPropagatedRegionConstraints(FemModel *femmodel,
                                      const CrevasseCalvingConfig &config,
                                      const CrevasseNodalFields &fields,
                                      const CrevassePropagationResult &result) {
  for (Object *&object : femmodel->elements->objects) {
    Element *element = xDynamicCast<Element *>(object);
    int numnodes = element->GetNumberOfNodes();
    Gauss *gauss = element->NewGauss();
    Input *levelset_input = element->GetInput(MaskIceLevelsetEnum);
    _assert_(levelset_input);

    for (int in = 0; in < numnodes; in++) {
      gauss->GaussVertex(in);
      Node *node = element->GetNode(in);
      if (!node->IsActive())
        continue;
      if (HasStaticLevelsetSpc(element, in, gauss))
        continue;

      int lid = node->Lid();
      if (result.applied_nodes && result.applied_nodes[lid] > 0.5) {
        node->ApplyConstraint(0, kCalvedLevelsetValue);
      } else if (!config.advect_icefront && fields.seed_nodes[lid] > 0.5) {
        IssmDouble current_levelset;
        levelset_input->GetInputValue(&current_levelset, gauss);
        node->ApplyConstraint(0, current_levelset);
      } else {
        node->DofInFSet(0);
      }
    }
    delete gauss;
  }
}

} // namespace

void LevelsetAnalysis::CreateConstraints(Constraints *constraints,
                                         IoModel *iomodel) { /*{{{*/

  /*intermediary: */
  int finiteelement;
  int code, vector_layout;
  IssmDouble *spcdata = NULL;
  int M, N;

  /*Get finite element type for this analysis*/
  // cout << "[LevelsetAnalysis::CreateConstraints]: Finding md.levelset.fe" <<
  // endl;
  iomodel->FindConstant(&finiteelement, "md.levelset.fe");

  /*First of, find the record for the enum, and get code  of data type: */
  iomodel->SetFilePointerToData(&code, &vector_layout,
                                "md.levelset.spclevelset");
  if (code != 7)
    _error_("expecting a IssmDouble vector for constraints "
            "md.levelset.spclevelset");
  if (vector_layout != 1)
    _error_("expecting a nodal vector for constraints md.levelset.spclevelset");

  /*Fetch vector:*/
  iomodel->FetchData(&spcdata, &M, &N, "md.levelset.spclevelset");

  /*Call IoModelToConstraintsx*/
  if (N > 1) {
    /*If it is a time series, most likely we are forcing the ice front position
     * and do not want to have a Dynamic Constraint*/
    _assert_(M == iomodel->numberofvertices + 1);
    IoModelToConstraintsx(constraints, iomodel, spcdata, M, N,
                          LevelsetAnalysisEnum, finiteelement);
  } else {
    /*This is not a time series, we probably have calving on, we need the
     * levelset constraints to update as the levelset moves*/
    _assert_(N == 1);
    _assert_(M == iomodel->numberofvertices);
    IoModelToDynamicConstraintsx(constraints, iomodel, spcdata, M, N,
                                 LevelsetAnalysisEnum, finiteelement);
  }

  /*Clean up*/
  xDelete<IssmDouble>(spcdata);
}
/*}}}*/
void LevelsetAnalysis::CreateLoads(Loads *loads, IoModel *iomodel) { /*{{{*/
  return;
} /*}}}*/
void LevelsetAnalysis::CreateNodes(Nodes *nodes, IoModel *iomodel,
                                   bool isamr) { /*{{{*/
  int finiteelement;
  // cout << "[LevelsetAnalysis::CreateNodes]: Finding md.levelset.fe" << endl;
  if (!isamr)
    iomodel->FindConstant(&finiteelement, "md.levelset.fe");
  else
    // [TODO]: Somehow, using the buttress based cd law requires me to hard code
    // finiteelement here... Otherwise, it returns an error IoModel.cpp could
    // not find "md.levelset.fe"...
    finiteelement = P1Enum;
  // cout << "[LevelsetAnalysis::CreateNodes]: md.levelset.fe=" << finiteelement
  // << endl;
  if (iomodel->domaintype != Domain2DhorizontalEnum)
    iomodel->FetchData(2, "md.mesh.vertexonbase", "md.mesh.vertexonsurface");
  ::CreateNodes(nodes, iomodel, LevelsetAnalysisEnum, finiteelement);
  iomodel->DeleteData(2, "md.mesh.vertexonbase", "md.mesh.vertexonsurface");
}
/*}}}*/
int LevelsetAnalysis::DofsPerNode(int **doflist, int domaintype,
                                  int approximation) { /*{{{*/
  return 1;
}
/*}}}*/
void LevelsetAnalysis::UpdateElements(Elements *elements, Inputs *inputs,
                                      IoModel *iomodel, int analysis_counter,
                                      int analysis_type) { /*{{{*/

  /*Finite element type*/
  int finiteelement;
  // cout << "[UpdateElements]: FindingConstant" << endl;
  iomodel->FindConstant(&finiteelement, "md.levelset.fe");

  /*Update elements: */
  int counter = 0;
  for (int i = 0; i < iomodel->numberofelements; i++) {
    if (iomodel->my_elements[i]) {
      Element *element = (Element *)elements->GetObjectByOffset(counter);
      element->Update(inputs, i, iomodel, analysis_counter, analysis_type,
                      finiteelement);
      counter++;
    }
  }

  iomodel->FetchDataToInput(inputs, elements, "md.mask.ice_levelset",
                            MaskIceLevelsetEnum);
  iomodel->FetchDataToInput(inputs, elements, "md.mask.ocean_levelset",
                            MaskOceanLevelsetEnum);
  iomodel->FetchDataToInput(inputs, elements, "md.levelset.spclevelset",
                            SpcLevelsetEnum);
  iomodel->FetchDataToInput(inputs, elements, "md.initialization.vx", VxEnum);
  iomodel->FetchDataToInput(inputs, elements, "md.initialization.vy", VyEnum);

  /*Get moving front parameters*/
  bool isstochastic;
  int calvinglaw;
  bool advect_icefront = true; // default to true
  iomodel->FindConstant(&calvinglaw, "md.calving.law");
  iomodel->FindConstant(&isstochastic,
                        "md.stochasticforcing.isstochasticforcing");

  // Check if we should disable ice front advection for crevasse depth calving
  //    if(calvinglaw == CalvingCrevasseDepthEnum){
  //        if(iomodel->FindConstant(&advect_icefront,"md.calving.advect_icefront")){
  //            advect_icefront =
  //            iomodel->FindConstant(&advect_icefront,"md.calving.advect_icefront");
  //        }
  //    }
  switch (calvinglaw) {

    /*"Continuous" calving laws*/
  case DefaultCalvingEnum:
    iomodel->FetchDataToInput(inputs, elements, "md.calving.calvingrate",
                              CalvingCalvingrateEnum);
    if (isstochastic) {
      iomodel->FetchDataToInput(inputs, elements,
                                "md.stochasticforcing.default_id",
                                StochasticForcingDefaultIdEnum);
      iomodel->FetchDataToInput(inputs, elements, "md.calving.calvingrate",
                                BaselineCalvingCalvingrateEnum);
    }
    break;
  case CalvingLevermannEnum:
    iomodel->FetchDataToInput(inputs, elements, "md.calving.coeff",
                              CalvinglevermannCoeffEnum);
    break;
  case CalvingVonmisesEnum:
    iomodel->FetchDataToInput(inputs, elements,
                              "md.calving.stress_threshold_groundedice",
                              CalvingStressThresholdGroundediceEnum);
    iomodel->FetchDataToInput(inputs, elements,
                              "md.calving.stress_threshold_floatingice",
                              CalvingStressThresholdFloatingiceEnum);
    iomodel->FetchDataToInput(inputs, elements, "md.geometry.bed", BedEnum);
    break;
  case CalvingVonmisesADEnum:
    iomodel->FetchDataToInput(inputs, elements, "md.calving.basin_id",
                              CalvingBasinIdEnum);
    iomodel->FetchDataToInput(inputs, elements, "md.geometry.bed", BedEnum);
    break;
  case CalvingDev2Enum:
    iomodel->FetchDataToInput(inputs, elements,
                              "md.calving.stress_threshold_groundedice",
                              CalvingStressThresholdGroundediceEnum);
    iomodel->FetchDataToInput(inputs, elements,
                              "md.calving.stress_threshold_floatingice",
                              CalvingStressThresholdFloatingiceEnum);
    break;
  case CalvingTestEnum:
    break;
  case CalvingParameterizationEnum:
    iomodel->FetchDataToInput(inputs, elements, "md.geometry.bed", BedEnum);
    break;
  case CalvingCalvingMIPEnum:
    break;

  /*"Discrete" calving laws (need to specify rate as 0 so that we can still
   * solve the level set equation)*/
  case CalvingMinthicknessEnum:
    iomodel->FetchDataToInput(inputs, elements, "md.geometry.bed", BedEnum);
    iomodel->ConstantToInput(inputs, elements, 0., CalvingratexEnum, P1Enum);
    iomodel->ConstantToInput(inputs, elements, 0., CalvingrateyEnum, P1Enum);
    break;
  case CalvingHabEnum:
    iomodel->FetchDataToInput(inputs, elements, "md.calving.flotation_fraction",
                              CalvingHabFractionEnum);
    iomodel->ConstantToInput(inputs, elements, 0., CalvingratexEnum, P1Enum);
    iomodel->ConstantToInput(inputs, elements, 0., CalvingrateyEnum, P1Enum);
    break;
  case CalvingCrevasseDepthEnum:
    iomodel->FetchDataToInput(inputs, elements, "md.calving.water_height",
                              WaterheightEnum);
    iomodel->FetchDataToInput(inputs, elements, "md.calving.meltwater_mask",
                              CalvingMeltwaterMaskEnum);
    iomodel->ConstantToInput(inputs, elements, 0.,
                             HydrofracturePredictedEnum, P1Enum);
    iomodel->ConstantToInput(inputs, elements, 0.,
                             HydrofractureWeakIceEnum, P1Enum);
    iomodel->ConstantToInput(inputs, elements, 0., CalvingratexEnum, P1Enum);
    iomodel->ConstantToInput(inputs, elements, 0., CalvingrateyEnum, P1Enum);
    break;
  case CalvingPollardEnum:
    break;

  default:
    _error_("Calving law " << EnumToStringx(calvinglaw)
                           << " not supported yet");
  }

  /*Get frontal melt parameters*/
  int melt_parameterization;
  iomodel->FindConstant(&melt_parameterization,
                        "md.frontalforcings.parameterization");
  switch (melt_parameterization) {
  case FrontalForcingsDefaultEnum:
    iomodel->FetchDataToInput(inputs, elements,
                              "md.frontalforcings.meltingrate",
                              CalvingMeltingrateEnum);
    if ((calvinglaw == CalvingParameterizationEnum) ||
        (calvinglaw == CalvingCalvingMIPEnum)) {
      iomodel->FetchDataToInput(inputs, elements,
                                "md.frontalforcings.ablationrate",
                                CalvingAblationrateEnum);
    }
    break;
  case FrontalForcingsRignotEnum:
    /*Retrieve thermal forcing only in the case of non-arma
     * FrontalForcingsRignot*/
    iomodel->FetchDataToInput(inputs, elements,
                              "md.frontalforcings.thermalforcing",
                              ThermalForcingEnum);
    iomodel->FetchDataToInput(inputs, elements, "md.frontalforcings.basin_id",
                              FrontalForcingsBasinIdEnum);
    iomodel->FetchDataToInput(inputs, elements,
                              "md.frontalforcings.subglacial_discharge",
                              FrontalForcingsSubglacialDischargeEnum);
    break;
  case FrontalForcingsRignotarmaEnum:
    bool isdischargearma;
    iomodel->FindConstant(&isdischargearma,
                          "md.frontalforcings.isdischargearma");
    iomodel->FetchDataToInput(inputs, elements, "md.frontalforcings.basin_id",
                              FrontalForcingsBasinIdEnum);
    if (isdischargearma == false)
      iomodel->FetchDataToInput(inputs, elements,
                                "md.frontalforcings.subglacial_discharge",
                                FrontalForcingsSubglacialDischargeEnum);
    break;
  default:
    _error_("Frontal forcings" << EnumToStringx(melt_parameterization)
                               << " not supported yet");
  }
}
/*}}}*/
void LevelsetAnalysis::UpdateParameters(Parameters *parameters,
                                        IoModel *iomodel, int solution_enum,
                                        int analysis_enum) { /*{{{*/

  parameters->AddObject(iomodel->CopyConstantObject("md.levelset.stabilization",
                                                    LevelsetStabilizationEnum));
  parameters->AddObject(iomodel->CopyConstantObject(
      "md.levelset.reinit_frequency", LevelsetReinitFrequencyEnum));
  parameters->AddObject(iomodel->CopyConstantObject("md.levelset.kill_icebergs",
                                                    LevelsetKillIcebergsEnum));
  parameters->AddObject(iomodel->CopyConstantObject("md.levelset.migration_max",
                                                    MigrationMaxEnum));

  int calvinglaw;
  IssmDouble *transparam = NULL;
  IssmDouble yts;
  int N, M;
  bool interp, cycle;

  iomodel->FindConstant(&calvinglaw, "md.calving.law");
  switch (calvinglaw) {
  case DefaultCalvingEnum:
  case CalvingLevermannEnum:
    break;
  case CalvingVonmisesEnum:
    parameters->AddObject(iomodel->CopyConstantObject(
        "md.calving.min_thickness", CalvingMinthicknessEnum));
    break;
  case CalvingVonmisesADEnum:
    parameters->AddObject(iomodel->CopyConstantObject(
        "md.calving.min_thickness", CalvingMinthicknessEnum));
    parameters->AddObject(iomodel->CopyConstantObject(
        "md.calving.num_basins", CalvingNumberofBasinsEnum));

    iomodel->FetchData(&transparam, &M, &N,
                       "md.calving.stress_threshold_groundedice");
    _assert_(M >= 1 && N >= 1);
    parameters->AddObject(new DoubleVecParam(
        CalvingADStressThresholdGroundediceEnum, transparam, M));
    xDelete<IssmDouble>(transparam);

    iomodel->FetchData(&transparam, &M, &N,
                       "md.calving.stress_threshold_floatingice");
    _assert_(M >= 1 && N >= 1);
    parameters->AddObject(new DoubleVecParam(
        CalvingADStressThresholdFloatingiceEnum, transparam, M));
    xDelete<IssmDouble>(transparam);

    break;
  case CalvingMinthicknessEnum:
    parameters->AddObject(iomodel->CopyConstantObject(
        "md.calving.min_thickness", CalvingMinthicknessEnum));
    break;
  case CalvingHabEnum:
    break;
  case CalvingCrevasseDepthEnum:
    parameters->AddObject(iomodel->CopyConstantObject(
        "md.calving.crevasse_opening_stress", CalvingCrevasseDepthEnum));
    parameters->AddObject(iomodel->CopyConstantObject(
        "md.calving.critical_stress", CalvingCriticalStressEnum));
    parameters->AddObject(iomodel->CopyConstantObject(
        "md.calving.timescale", CrevasseDepthCalvingTimescaleEnum));
    parameters->AddObject(iomodel->CopyConstantObject(
        "md.calving.min_iceberg_size", CalvingMinIcebergSizeEnum));
    parameters->AddObject(iomodel->CopyConstantObject(
        "md.calving.advect_icefront", CalvingAdvectIcefrontEnum));
    parameters->AddObject(iomodel->CopyConstantObject(
        "md.calving.propagate_from_front", CalvingPropagateFromFrontEnum));
    parameters->AddObject(iomodel->CopyConstantObject(
        "md.calving.remove_only_marked", CalvingRemoveOnlyMarkedEnum));
    parameters->AddObject(iomodel->CopyConstantObject(
        "md.calving.hydrofracture_stabilization",
        CalvingHydrofractureStabilizationEnum));
    parameters->AddObject(iomodel->CopyConstantObject(
        "md.calving.hydrofracture_weakening_factor",
        CalvingHydrofractureWeakeningFactorEnum));
    parameters->AddObject(iomodel->CopyConstantObject(
        "md.calving.hydrofracture_min_ocean_levelset",
        CalvingHydrofractureMinOceanLevelsetEnum));
    /* Initialize the scheduled calving clock. The runtime calving routine
     * resets this to md.timestepping.start_time before the first check so
     * restarted/nonzero-start transients calve exactly one timescale after
     * their own start time instead of being measured from model year zero. */
    parameters->AddObject(new DoubleParam(LastCalvingTimeEnum, 0.0));
    // Initialize calving occurred flag to 0
    parameters->AddObject(new DoubleParam(CalvingOccurredEnum, 0.0));
    // Initialize propagated calving area to 0
    parameters->AddObject(new DoubleParam(CalvingPropagatedAreaEnum, 0.0));
    // Initialize propagated minimum x-coordinate to large value
    parameters->AddObject(new DoubleParam(CalvingPropagatedMinXEnum, 1e20));
    break;
  case CalvingDev2Enum:
    parameters->AddObject(
        iomodel->CopyConstantObject("md.calving.height_above_floatation",
                                    CalvingHeightAboveFloatationEnum));
    break;
  case CalvingTestEnum:
    iomodel->FindConstant(&interp, "md.timestepping.interp_forcing");
    iomodel->FindConstant(&cycle, "md.timestepping.cycle_forcing");
    iomodel->FetchData(&transparam, &N, &M, "md.calving.speedfactor");
    if (N == 1) {
      _assert_(M == 1);
      parameters->AddObject(
          new DoubleParam(CalvingTestSpeedfactorEnum, transparam[0]));
    } else {
      _assert_(N == 2);
      parameters->AddObject(new TransientParam(CalvingTestSpeedfactorEnum,
                                               &transparam[0], &transparam[M],
                                               interp, cycle, M));
    }
    xDelete<IssmDouble>(transparam);
    iomodel->FetchData(&transparam, &N, &M, "md.calving.independentrate");
    if (N == 1) {
      _assert_(M == 1);
      parameters->AddObject(
          new DoubleParam(CalvingTestIndependentRateEnum, transparam[0]));
    } else {
      _assert_(N == 2);
      parameters->AddObject(new TransientParam(CalvingTestIndependentRateEnum,
                                               &transparam[0], &transparam[M],
                                               interp, cycle, M));
    }
    xDelete<IssmDouble>(transparam);
    break;
  case CalvingParameterizationEnum:
    parameters->AddObject(iomodel->CopyConstantObject("md.calving.use_param",
                                                      CalvingUseParamEnum));
    parameters->AddObject(iomodel->CopyConstantObject(
        "md.calving.min_thickness", CalvingMinthicknessEnum));
    parameters->AddObject(
        iomodel->CopyConstantObject("md.calving.theta", CalvingThetaEnum));
    parameters->AddObject(
        iomodel->CopyConstantObject("md.calving.alpha", CalvingAlphaEnum));
    parameters->AddObject(
        iomodel->CopyConstantObject("md.calving.xoffset", CalvingXoffsetEnum));
    parameters->AddObject(
        iomodel->CopyConstantObject("md.calving.yoffset", CalvingYoffsetEnum));
    parameters->AddObject(iomodel->CopyConstantObject(
        "md.calving.vel_lowerbound", CalvingVelLowerboundEnum));
    parameters->AddObject(iomodel->CopyConstantObject(
        "md.calving.vel_threshold", CalvingVelThresholdEnum));
    parameters->AddObject(iomodel->CopyConstantObject(
        "md.calving.vel_upperbound", CalvingVelUpperboundEnum));
    break;
  case CalvingPollardEnum:
    parameters->AddObject(
        iomodel->CopyConstantObject("md.calving.rc", CalvingRcEnum));
    break;
  case CalvingCalvingMIPEnum:
    parameters->AddObject(iomodel->CopyConstantObject("md.calving.experiment",
                                                      CalvingUseParamEnum));
    parameters->AddObject(iomodel->CopyConstantObject(
        "md.calving.min_thickness", CalvingMinthicknessEnum));
    break;
  default:
    _error_("Calving law " << EnumToStringx(calvinglaw)
                           << " not supported yet");
  }

  /*Get frontal melt parameters*/
  int melt_parameterization;
  iomodel->FindConstant(&melt_parameterization,
                        "md.frontalforcings.parameterization");
  switch (melt_parameterization) {
  case FrontalForcingsDefaultEnum:
    break;
  case FrontalForcingsRignotarmaEnum:
    parameters->AddObject(iomodel->CopyConstantObject(
        "md.frontalforcings.num_basins", FrontalForcingsNumberofBasinsEnum));
    /*Retrieve thermal forcing parameters*/
    parameters->AddObject(iomodel->CopyConstantObject(
        "md.frontalforcings.num_params", FrontalForcingsNumberofParamsEnum));
    parameters->AddObject(iomodel->CopyConstantObject(
        "md.frontalforcings.num_breaks", FrontalForcingsNumberofBreaksEnum));
    parameters->AddObject(
        iomodel->CopyConstantObject("md.frontalforcings.monthlyvals_numbreaks",
                                    FrontalForcingsNumberofMonthBreaksEnum));
    parameters->AddObject(iomodel->CopyConstantObject(
        "md.frontalforcings.ar_order", FrontalForcingsARMAarOrderEnum));
    parameters->AddObject(iomodel->CopyConstantObject(
        "md.frontalforcings.ma_order", FrontalForcingsARMAmaOrderEnum));
    parameters->AddObject(iomodel->CopyConstantObject(
        "md.frontalforcings.arma_timestep", FrontalForcingsARMATimestepEnum));
    iomodel->FetchData(&transparam, &M, &N,
                       "md.frontalforcings.polynomialparams");
    parameters->AddObject(new DoubleMatParam(FrontalForcingsARMApolyparamsEnum,
                                             transparam, M, N));
    xDelete<IssmDouble>(transparam);
    iomodel->FetchData(&transparam, &M, &N, "md.frontalforcings.datebreaks");
    parameters->AddObject(new DoubleMatParam(FrontalForcingsARMAdatebreaksEnum,
                                             transparam, M, N));
    xDelete<IssmDouble>(transparam);
    iomodel->FetchData(&transparam, &M, &N, "md.frontalforcings.arlag_coefs");
    parameters->AddObject(new DoubleMatParam(FrontalForcingsARMAarlagcoefsEnum,
                                             transparam, M, N));
    xDelete<IssmDouble>(transparam);
    iomodel->FetchData(&transparam, &M, &N, "md.frontalforcings.malag_coefs");
    parameters->AddObject(new DoubleMatParam(FrontalForcingsARMAmalagcoefsEnum,
                                             transparam, M, N));
    xDelete<IssmDouble>(transparam);
    iomodel->FetchData(&transparam, &M, &N,
                       "md.frontalforcings.monthlyvals_datebreaks");
    parameters->AddObject(new DoubleMatParam(
        FrontalForcingsARMAmonthdatebreaksEnum, transparam, M, N));
    xDelete<IssmDouble>(transparam);
    iomodel->FetchData(&transparam, &M, &N,
                       "md.frontalforcings.monthlyvals_intercepts");
    parameters->AddObject(new DoubleMatParam(
        FrontalForcingsARMAmonthinterceptsEnum, transparam, M, N));
    xDelete<IssmDouble>(transparam);
    iomodel->FetchData(&transparam, &M, &N,
                       "md.frontalforcings.monthlyvals_trends");
    parameters->AddObject(new DoubleMatParam(FrontalForcingsARMAmonthtrendsEnum,
                                             transparam, M, N));
    xDelete<IssmDouble>(transparam);
    parameters->AddObject(
        iomodel->CopyConstantObject("md.frontalforcings.isdischargearma",
                                    FrontalForcingsIsDischargeARMAEnum));
    /*Retrieve subglacial discharge parameters */
    bool isdischargearma;
    parameters->FindParam(&isdischargearma, FrontalForcingsIsDischargeARMAEnum);
    if (isdischargearma == true) {
      parameters->AddObject(
          iomodel->CopyConstantObject("md.frontalforcings.sd_num_params",
                                      FrontalForcingsSdNumberofParamsEnum));
      parameters->AddObject(
          iomodel->CopyConstantObject("md.frontalforcings.sd_num_breaks",
                                      FrontalForcingsSdNumberofBreaksEnum));
      parameters->AddObject(iomodel->CopyConstantObject(
          "md.frontalforcings.sd_ar_order", FrontalForcingsSdarOrderEnum));
      parameters->AddObject(iomodel->CopyConstantObject(
          "md.frontalforcings.sd_ma_order", FrontalForcingsSdmaOrderEnum));
      parameters->AddObject(
          iomodel->CopyConstantObject("md.frontalforcings.sd_arma_timestep",
                                      FrontalForcingsSdARMATimestepEnum));
      iomodel->FetchData(&transparam, &M, &N,
                         "md.frontalforcings.sd_polynomialparams");
      parameters->AddObject(new DoubleMatParam(FrontalForcingsSdpolyparamsEnum,
                                               transparam, M, N));
      xDelete<IssmDouble>(transparam);
      iomodel->FetchData(&transparam, &M, &N,
                         "md.frontalforcings.sd_datebreaks");
      parameters->AddObject(new DoubleMatParam(FrontalForcingsSddatebreaksEnum,
                                               transparam, M, N));
      xDelete<IssmDouble>(transparam);
      iomodel->FetchData(&transparam, &M, &N,
                         "md.frontalforcings.sd_arlag_coefs");
      parameters->AddObject(new DoubleMatParam(FrontalForcingsSdarlagcoefsEnum,
                                               transparam, M, N));
      xDelete<IssmDouble>(transparam);
      iomodel->FetchData(&transparam, &M, &N,
                         "md.frontalforcings.sd_malag_coefs");
      parameters->AddObject(new DoubleMatParam(FrontalForcingsSdmalagcoefsEnum,
                                               transparam, M, N));
      xDelete<IssmDouble>(transparam);
      iomodel->FetchData(&transparam, &M, &N,
                         "md.frontalforcings.sd_monthlyfrac");
      parameters->AddObject(new DoubleMatParam(FrontalForcingsSdMonthlyFracEnum,
                                               transparam, M, N));
      xDelete<IssmDouble>(transparam);
    }
    break;
  case FrontalForcingsRignotEnum:
    parameters->AddObject(iomodel->CopyConstantObject(
        "md.frontalforcings.num_basins", FrontalForcingsNumberofBasinsEnum));
    break;
  default:
    _error_("Frontal forcings " << EnumToStringx(melt_parameterization)
                                << " not supported yet");
  }
}
/*}}}*/

/*Finite element Analysis*/
void LevelsetAnalysis::Core(FemModel *femmodel) { /*{{{*/

  /*parameters: */
  int stabilization;
  femmodel->parameters->FindParam(&stabilization, LevelsetStabilizationEnum);

  /*activate formulation: */
  femmodel->SetCurrentConfiguration(LevelsetAnalysisEnum);

  if (VerboseSolution())
    _printf0_("   call computational core:\n");
  if (stabilization == 4) {
    solutionsequence_fct(femmodel);
  } else {
    solutionsequence_linear(femmodel);
  }
} /*}}}*/
void LevelsetAnalysis::PreCore(FemModel *femmodel) { /*{{{*/
  _error_("not implemented");
} /*}}}*/
ElementVector *LevelsetAnalysis::CreateDVector(Element *element) { /*{{{*/
  /*Default, return NULL*/
  return NULL;
} /*}}}*/
ElementMatrix *
LevelsetAnalysis::CreateJacobianMatrix(Element *element) { /*{{{*/
  /* Jacobian required for the Newton solver */
  _error_("not implemented yet");
} /*}}}*/
ElementMatrix *LevelsetAnalysis::CreateKMatrix(Element *element) { /*{{{*/

  if (!element->IsOnBase())
    return NULL;
  Element *basalelement = element->SpawnBasalElement();

  /*Intermediaries */
  int stabilization, dim, domaintype;
  int i, j, k, row, col;
  IssmDouble kappa, factor;
  IssmDouble Jdet, dt, D_scalar;
  IssmDouble h, hx, hy, hz;
  IssmDouble vel, w[3];
  IssmDouble migrationmax;
  IssmDouble *xyz_list = NULL;

  /*Get problem dimension and whether there is moving front or not*/
  basalelement->FindParam(&domaintype, DomainTypeEnum);
  basalelement->FindParam(&stabilization, LevelsetStabilizationEnum);
  switch (domaintype) {
  case Domain2DverticalEnum:
    dim = 1;
    break;
  case Domain2DhorizontalEnum:
    dim = 2;
    break;
  case Domain3DEnum:
    dim = 2;
    break;
  default:
    _error_("mesh " << EnumToStringx(domaintype) << " not supported yet");
  }
  /*Fetch number of nodes and dof for this finite element*/
  int numnodes = basalelement->GetNumberOfNodes();

  /*Initialize Element vector and other vectors*/
  ElementMatrix *Ke = basalelement->NewElementMatrix();
  IssmDouble *basis = xNew<IssmDouble>(numnodes);
  IssmDouble *dbasis = xNew<IssmDouble>(2 * numnodes);

  /*Retrieve all inputs and parameters*/
  basalelement->GetVerticesCoordinates(&xyz_list);
  basalelement->FindParam(&dt, TimesteppingTimeStepEnum);
  basalelement->FindParam(&migrationmax, MigrationMaxEnum);

  h = basalelement->CharacteristicLength();

  Input *mf_vx_input = NULL;
  Input *mf_vy_input = NULL;

  /*Load velocities*/
  switch (domaintype) {
  case Domain2DverticalEnum:
    mf_vx_input = basalelement->GetInput(MovingFrontalVxEnum);
    _assert_(mf_vx_input);
    break;
  case Domain2DhorizontalEnum:
    mf_vx_input = basalelement->GetInput(MovingFrontalVxEnum);
    _assert_(mf_vx_input);
    mf_vy_input = basalelement->GetInput(MovingFrontalVyEnum);
    _assert_(mf_vy_input);
    break;
  case Domain3DEnum:
    mf_vx_input = basalelement->GetInput(MovingFrontalVxEnum);
    _assert_(mf_vx_input);
    mf_vy_input = basalelement->GetInput(MovingFrontalVyEnum);
    _assert_(mf_vy_input);
    break;
  default:
    _error_("mesh " << EnumToStringx(domaintype) << " not supported yet");
  }

  /* Start  looping on the number of gaussian points: */
  Gauss *gauss = basalelement->NewGauss(2);
  while (gauss->next()) {

    basalelement->JacobianDeterminant(&Jdet, xyz_list, gauss);
    basalelement->NodalFunctions(basis, gauss);
    basalelement->NodalFunctionsDerivatives(dbasis, xyz_list, gauss);
    D_scalar = gauss->weight * Jdet;

    /* Transient */
    if (dt != 0.) {
      for (i = 0; i < numnodes; i++) {
        for (j = 0; j < numnodes; j++) {
          Ke->values[i * numnodes + j] += D_scalar * basis[j] * basis[i];
        }
      }
      D_scalar = D_scalar * dt;
    }

    /* Levelset speed */
    mf_vx_input->GetInputValue(&w[0], gauss);
    mf_vy_input->GetInputValue(&w[1], gauss);

    /* Apply limiter to the migration rate */
    vel = 0.;
    for (i = 0; i < dim; i++)
      vel += w[i] * w[i];
    vel = sqrt(vel) + 1e-14;
    /* !!NOTE: This is different from the previous version 25838 (and before).
     * The current threshold restrict both advance and retreat velocity. */
    if (vel > migrationmax) {
      for (i = 0; i < dim; i++)
        w[i] = w[i] / vel * migrationmax;
    }

    /*Compute D*/
    for (i = 0; i < numnodes; i++) {
      for (j = 0; j < numnodes; j++) {
        for (k = 0; k < dim; k++) {
          Ke->values[i * numnodes + j] +=
              D_scalar * w[k] * dbasis[k * numnodes + j] * basis[i];
        }
      }
    }

    /* Stabilization */
    vel = 0.;
    for (i = 0; i < dim; i++)
      vel += w[i] * w[i];
    vel = sqrt(vel) + 1.e-14;
    switch (stabilization) {
    case 0:
      /*Nothing to be done*/
      break;
    case 1:
      /* Artificial Diffusion */
      basalelement->ElementSizes(&hx, &hy, &hz);
      h = sqrt(pow(hx * w[0] / vel, 2) + pow(hy * w[1] / vel, 2));
      kappa = h * vel / 2.;
      for (i = 0; i < numnodes; i++) {
        for (j = 0; j < numnodes; j++) {
          for (k = 0; k < dim; k++) {
            Ke->values[i * numnodes + j] += D_scalar * kappa *
                                            dbasis[k * numnodes + j] *
                                            dbasis[k * numnodes + i];
          }
        }
      }
      break;
    case 2: {
      /* Streamline Upwinding */
      mf_vx_input->GetInputAverage(&w[0]);
      mf_vy_input->GetInputAverage(&w[1]);
      vel = sqrt(w[0] * w[0] + w[1] * w[1]) + 1.e-8;
      IssmDouble tau = h / (2 * vel);
      factor = dt * gauss->weight * Jdet * tau;
      for (int i = 0; i < numnodes; i++) {
        for (int j = 0; j < numnodes; j++) {
          Ke->values[i * numnodes + j] += factor *
                                          (w[0] * dbasis[0 * numnodes + i] +
                                           w[1] * dbasis[1 * numnodes + i]) *
                                          (w[0] * dbasis[0 * numnodes + j] +
                                           w[1] * dbasis[1 * numnodes + j]);
        }
      }
    } break;
    case 5: {
      /*SUPG*/
      IssmDouble vx, vy;
      mf_vx_input->GetInputAverage(&vx);
      mf_vy_input->GetInputAverage(&vy);
      vel = sqrt(vx * vx + vy * vy) + 1.e-8;
      IssmPDouble xi = 1.;
      IssmDouble tau = xi * h / (2 * vel);

      /*Mass matrix - part 2*/
      factor = gauss->weight * Jdet * tau;
      for (int i = 0; i < numnodes; i++) {
        for (int j = 0; j < numnodes; j++) {
          Ke->values[i * numnodes + j] +=
              factor * basis[j] *
              (vx * dbasis[0 * numnodes + i] + vy * dbasis[1 * numnodes + i]);
        }
      }

      /*Advection matrix - part 2, A*/
      factor = dt * gauss->weight * Jdet * tau;
      for (int i = 0; i < numnodes; i++) {
        for (int j = 0; j < numnodes; j++) {
          Ke->values[i * numnodes + j] +=
              factor *
              (vx * dbasis[0 * numnodes + j] + vy * dbasis[1 * numnodes + j]) *
              (vx * dbasis[0 * numnodes + i] + vy * dbasis[1 * numnodes + i]);
        }
      }

      break;
    }
    case 6: {
      /*SUPG*/
      IssmDouble vx, vy;
      mf_vx_input->GetInputAverage(&vx);
      mf_vy_input->GetInputAverage(&vy);
      vel = sqrt(vx * vx + vy * vy) + 1.e-8;
      IssmDouble ECN, K;
      ECN = vel * dt / h;
      K = 1. / tanh(ECN) - 1. / ECN;
      //				if (ECN<1e-6) K = ECN /3.0;

      /*According to Hilmar, xi=K is too large*/
      IssmDouble xi = 0.1 * K;

      IssmDouble tau = xi * h / (2 * vel);
      Input *levelset_input = NULL;

      IssmDouble kappa;
      IssmDouble p = 4, q = 4;
      IssmDouble phi[3];

      levelset_input = basalelement->GetInput(MaskIceLevelsetEnum);
      _assert_(levelset_input);
      levelset_input->GetInputValue(&phi[0], gauss);

      IssmDouble dphidx = 0., dphidy = 0.;
      IssmDouble nphi;

      for (int i = 0; i < numnodes; i++) {
        dphidx += phi[i] * dbasis[0 * numnodes + i];
        dphidy += phi[i] * dbasis[1 * numnodes + i];
      }
      nphi = sqrt(dphidx * dphidx + dphidy * dphidy);

      if (nphi >= 1) {
        kappa = 1 - 1.0 / nphi;
      } else {
        kappa = 0.5 / M_PI * sin(2 * M_PI * nphi) / nphi;
      }

      kappa = kappa * vel / h;

      /*Mass matrix - part 2*/
      factor = gauss->weight * Jdet * tau;
      for (int i = 0; i < numnodes; i++) {
        for (int j = 0; j < numnodes; j++) {
          Ke->values[i * numnodes + j] +=
              factor * basis[j] *
              (vx * dbasis[0 * numnodes + i] + vy * dbasis[1 * numnodes + i]);
        }
      }

      /*Advection matrix - part 2, A*/
      factor = dt * gauss->weight * Jdet * tau;
      for (int i = 0; i < numnodes; i++) {
        for (int j = 0; j < numnodes; j++) {
          Ke->values[i * numnodes + j] +=
              factor *
              (vx * dbasis[0 * numnodes + j] + vy * dbasis[1 * numnodes + j]) *
              (vx * dbasis[0 * numnodes + i] + vy * dbasis[1 * numnodes + i]);
        }
      }
      /*Add the pertubation term \nabla\cdot(\kappa*\nabla\phi)*/
      factor = dt * gauss->weight * Jdet * kappa;
      for (int i = 0; i < numnodes; i++) {
        for (int j = 0; j < numnodes; j++) {
          for (int k = 0; k < dim; k++) {
            Ke->values[i * numnodes + j] +=
                factor * dbasis[k * numnodes + j] * dbasis[k * numnodes + i];
          }
        }
      }

      break;
    }
    default:
      _error_("unknown type of stabilization in LevelsetAnalysis.cpp");
    }
  }

  /*Clean up and return*/
  xDelete<IssmDouble>(xyz_list);
  xDelete<IssmDouble>(basis);
  xDelete<IssmDouble>(dbasis);
  delete gauss;
  if (basalelement->IsSpawnedElement()) {
    basalelement->DeleteMaterials();
    delete basalelement;
  };
  return Ke;
} /*}}}*/
ElementVector *LevelsetAnalysis::CreatePVector(Element *element) { /*{{{*/

  if (!element->IsOnBase())
    return NULL;
  Element *basalelement = element->SpawnBasalElement();

  /*Intermediaries */
  int domaintype, stabilization;
  IssmDouble Jdet, dt;
  IssmDouble lsf;
  IssmDouble *xyz_list = NULL;

  /*Fetch number of nodes and dof for this finite element*/
  int numnodes = basalelement->GetNumberOfNodes();
  basalelement->FindParam(&stabilization, LevelsetStabilizationEnum);

  /*Initialize Element vector*/
  ElementVector *pe = basalelement->NewElementVector();
  basalelement->FindParam(&dt, TimesteppingTimeStepEnum);
  _assert_(dt > 0.);

  /*Initialize basis vector*/
  IssmDouble *basis = xNew<IssmDouble>(numnodes);
  IssmDouble *dbasis = NULL;
  if ((stabilization == 5) | (stabilization == 6))
    dbasis = xNew<IssmDouble>(2 * numnodes);

  /*Retrieve all inputs and parameters*/
  basalelement->GetVerticesCoordinates(&xyz_list);
  Input *levelset_input = basalelement->GetInput(MaskIceLevelsetEnum);
  _assert_(levelset_input);
  Input *mf_vx_input = basalelement->GetInput(MovingFrontalVxEnum);
  _assert_(mf_vx_input);
  Input *mf_vy_input = basalelement->GetInput(MovingFrontalVyEnum);
  _assert_(mf_vy_input);

  IssmDouble h = basalelement->CharacteristicLength();

  /* Start  looping on the number of gaussian points: */
  Gauss *gauss = basalelement->NewGauss(2);
  while (gauss->next()) {
    basalelement->JacobianDeterminant(&Jdet, xyz_list, gauss);
    basalelement->NodalFunctions(basis, gauss);

    /* old function value */
    levelset_input->GetInputValue(&lsf, gauss);
    IssmDouble factor = Jdet * gauss->weight * lsf;
    for (int i = 0; i < numnodes; i++)
      pe->values[i] += factor * basis[i];

    if (stabilization == 5) { /*SUPG*/
      IssmDouble vx, vy, vel;
      basalelement->NodalFunctionsDerivatives(dbasis, xyz_list, gauss);
      mf_vx_input->GetInputAverage(&vx);
      mf_vy_input->GetInputAverage(&vy);
      vel = sqrt(vx * vx + vy * vy) + 1.e-8;
      IssmPDouble xi = 1.;
      IssmDouble tau = xi * h / (2 * vel);

      /*Force vector - part 2*/
      factor = Jdet * gauss->weight * lsf;
      for (int i = 0; i < numnodes; i++) {
        pe->values[i] += factor * (tau * vx * dbasis[0 * numnodes + i] +
                                   tau * vy * dbasis[1 * numnodes + i]);
      }
    } else if (stabilization == 6) {
      IssmDouble vx, vy, vel;
      basalelement->NodalFunctionsDerivatives(dbasis, xyz_list, gauss);
      mf_vx_input->GetInputAverage(&vx);
      mf_vy_input->GetInputAverage(&vy);
      vel = sqrt(vx * vx + vy * vy) + 1.e-8;

      IssmDouble ECN, K;
      ECN = vel * dt / h;
      K = 1. / tanh(ECN) - 1. / ECN;
      //		if (ECN<1e-6) K = ECN /3.0;

      /*According to Hilmar, xi=K is too large*/
      IssmDouble xi = 0.1 * K;

      IssmDouble tau = xi * h / (2 * vel);

      /*Force vector - part 2*/
      factor = Jdet * gauss->weight * lsf * tau;
      for (int i = 0; i < numnodes; i++) {
        pe->values[i] += factor * (vx * dbasis[0 * numnodes + i] +
                                   vy * dbasis[1 * numnodes + i]);
      }
    }
  }

  /*Clean up and return*/
  xDelete<IssmDouble>(xyz_list);
  xDelete<IssmDouble>(basis);
  xDelete<IssmDouble>(dbasis);
  basalelement->FindParam(&domaintype, DomainTypeEnum);
  if (basalelement->IsSpawnedElement()) {
    basalelement->DeleteMaterials();
    delete basalelement;
  };
  delete gauss;

  return pe;
} /*}}}*/
IssmDouble LevelsetAnalysis::GetDistanceToStraight(IssmDouble *q,
                                                   IssmDouble *s0,
                                                   IssmDouble *s1) { /*{{{*/
  // returns distance d of point q to straight going through points s0, s1
  // d=|a x b|/|b|
  // with a=q-s0, b=s1-s0

  /* Intermediaries */
  const int dim = 2;
  int i;
  IssmDouble a[dim], b[dim];
  IssmDouble norm_b;

  for (i = 0; i < dim; i++) {
    a[i] = q[i] - s0[i];
    b[i] = s1[i] - s0[i];
  }

  norm_b = 0.;
  for (i = 0; i < dim; i++)
    norm_b += b[i] * b[i];
  norm_b = sqrt(norm_b);
  _assert_(norm_b > 0.);

  return fabs(a[0] * b[1] - a[1] * b[0]) / norm_b;
} /*}}}*/
void LevelsetAnalysis::GetSolutionFromInputs(Vector<IssmDouble> *solution,
                                             Element *element) { /*{{{*/
  element->GetSolutionFromInputsOneDof(solution, MaskIceLevelsetEnum);
} /*}}}*/
void LevelsetAnalysis::GradientJ(Vector<IssmDouble> *gradient, Element *element,
                                 int control_type, int control_interp,
                                 int control_index) { /*{{{*/
  _error_("Not implemented yet");
} /*}}}*/
void LevelsetAnalysis::InputUpdateFromSolution(IssmDouble *solution,
                                               Element *element) { /*{{{*/

  int domaintype;
  element->FindParam(&domaintype, DomainTypeEnum);
  switch (domaintype) {
  case Domain2DhorizontalEnum:
    element->InputUpdateFromSolutionOneDof(solution, MaskIceLevelsetEnum);
    break;
  case Domain3DEnum:
    element->InputUpdateFromSolutionOneDofCollapsed(solution,
                                                    MaskIceLevelsetEnum);
    break;
  default:
    _error_("mesh " << EnumToStringx(domaintype) << " not supported yet");
  }
} /*}}}*/
void LevelsetAnalysis::PostProcess(FemModel *femmodel) { /*{{{*/

  /*This function is only used by "discrete calving laws" for which we change
   * the value of the levelset after the advection step (level set equation
   * solve) based on the law*/

  /*Intermediaries*/
  int calvinglaw;
  IssmDouble newlevelset[6];
  femmodel->parameters->FindParam(&calvinglaw, CalvingLawEnum);

  /*Apply minimum thickness criterion*/
  if (calvinglaw == CalvingMinthicknessEnum ||
      calvinglaw == CalvingVonmisesEnum ||
      calvinglaw == CalvingParameterizationEnum ||
      calvinglaw == CalvingVonmisesADEnum ||
      calvinglaw == CalvingCalvingMIPEnum) {

    IssmDouble mig_max = femmodel->parameters->FindParam(MigrationMaxEnum);
    IssmDouble dt = femmodel->parameters->FindParam(TimesteppingTimeStepEnum);

    /*Get current distance to terminus*/
    InputDuplicatex(femmodel, MaskIceLevelsetEnum, DistanceToCalvingfrontEnum);
    femmodel->DistanceToFieldValue(MaskIceLevelsetEnum, 0,
                                   DistanceToCalvingfrontEnum);

    /*Intermediaries*/
    IssmDouble thickness, bed, sealevel, distance, levelset;
    IssmDouble min_thickness =
        femmodel->parameters->FindParam(CalvingMinthicknessEnum);

    /*Loop over all elements of this partition*/
    for (Object *&object : femmodel->elements->objects) {
      Element *element = xDynamicCast<Element *>(object);

      /*no need to postprocess an ice free element*/
      if (!element->IsIceInElement())
        continue;

      int numnodes = element->GetNumberOfNodes();
      _assert_(numnodes < 7);
      Gauss *gauss = element->NewGauss();
      Input *H_input = element->GetInput(ThicknessEnum);
      _assert_(H_input);
      Input *b_input = element->GetInput(BedEnum);
      _assert_(b_input);
      Input *sl_input = element->GetInput(SealevelEnum);
      _assert_(sl_input);
      Input *dis_input = element->GetInput(DistanceToCalvingfrontEnum);
      _assert_(dis_input);
      Input *levelset_input = element->GetInput(MaskIceLevelsetEnum);
      _assert_(levelset_input);

      /*Potentially constrain nodes of this element*/
      for (int in = 0; in < numnodes; in++) {
        gauss->GaussNode(element->GetElementType(), in);

        levelset_input->GetInputValue(&levelset, gauss);
        H_input->GetInputValue(&thickness, gauss);
        b_input->GetInputValue(&bed, gauss);
        sl_input->GetInputValue(&sealevel, gauss);
        dis_input->GetInputValue(&distance, gauss);

        if (thickness < min_thickness && bed < sealevel &&
            fabs(distance) < mig_max * dt && levelset < 0) {
          newlevelset[in] =
              +400.; // Arbitrary > 0 number (i.e. deactivate this node)
        } else {
          newlevelset[in] = levelset;
        }
      }
      element->AddInput(MaskIceLevelsetEnum, &newlevelset[0],
                        element->GetElementType());
      delete gauss;
    }
  }
} /*}}}*/
void LevelsetAnalysis::UpdateConstraints(FemModel *femmodel) { /*{{{*/

  /*Intermediaries*/
  int calvinglaw;
  IssmDouble yts;
  femmodel->parameters->FindParam(&yts, ConstantsYtsEnum);
  femmodel->parameters->FindParam(&calvinglaw, CalvingLawEnum);
  IssmDouble mig_max = femmodel->parameters->FindParam(MigrationMaxEnum);
  IssmDouble dt = femmodel->parameters->FindParam(TimesteppingTimeStepEnum);

  /* Get current distance to terminus
   * Only do this if necessary, PostProcess is already doing it for a few
   * calving law Do not repeat the process is this function is particularly
   * slow*/
  bool computedistance = true;
  if (calvinglaw == CalvingMinthicknessEnum ||
      calvinglaw == CalvingVonmisesEnum ||
      calvinglaw == CalvingParameterizationEnum ||
      calvinglaw == CalvingVonmisesADEnum ||
      calvinglaw == CalvingCalvingMIPEnum) {
    int step;
    femmodel->parameters->FindParam(&step, StepEnum);
    if (step > 1) {
      computedistance = false;
    }
  }
  if (computedistance) {
    InputDuplicatex(femmodel, MaskIceLevelsetEnum, DistanceToCalvingfrontEnum);
    femmodel->DistanceToFieldValue(MaskIceLevelsetEnum, 0,
                                   DistanceToCalvingfrontEnum);
  }

  if (calvinglaw == CalvingHabEnum) {

    /*Intermediaries*/
    IssmDouble thickness, water_depth, distance, hab_fraction;

    /*Loop over all elements of this partition*/
    for (Object *&object : femmodel->elements->objects) {
      Element *element = xDynamicCast<Element *>(object);

      IssmDouble rho_ice = element->FindParam(MaterialsRhoIceEnum);
      IssmDouble rho_water = element->FindParam(MaterialsRhoSeawaterEnum);

      int numnodes = element->GetNumberOfNodes();
      Gauss *gauss = element->NewGauss();
      Input *H_input = element->GetInput(ThicknessEnum);
      _assert_(H_input);
      Input *bed_input = element->GetInput(BedEnum);
      _assert_(bed_input);
      Input *hab_fraction_input = element->GetInput(CalvingHabFractionEnum);
      _assert_(hab_fraction_input);
      Input *dis_input = element->GetInput(DistanceToCalvingfrontEnum);
      _assert_(dis_input);

      /*Potentially constrain nodes of this element*/
      for (int in = 0; in < numnodes; in++) {
        gauss->GaussNode(element->GetElementType(), in);
        Node *node = element->GetNode(in);
        if (!node->IsActive())
          continue;

        H_input->GetInputValue(&thickness, gauss);
        bed_input->GetInputValue(&water_depth, gauss);
        dis_input->GetInputValue(&distance, gauss);
        hab_fraction_input->GetInputValue(&hab_fraction, gauss);

        if (thickness <
                ((rho_water / rho_ice) * (1 + hab_fraction) * -water_depth) &&
            fabs(distance) < mig_max * dt) {
          node->ApplyConstraint(0, +1.);
        } else {
          /* no ice, set no spc */
          node->DofInFSet(0);
        }
      }
      delete gauss;
    }
  } else if (calvinglaw == CalvingCrevasseDepthEnum) {
    CrevasseCalvingConfig config;
    femmodel->parameters->FindParam(&config.time_yr, TimeEnum);
	    femmodel->parameters->FindParam(&config.dt_yr,
	                                    TimesteppingTimeStepEnum);
    IssmDouble start_time_yr;
    femmodel->parameters->FindParam(&start_time_yr, TimesteppingStartTimeEnum);
	    config.time_yr /= yts;
	    config.dt_yr /= yts;
    start_time_yr /= yts;
	    config.min_iceberg_size = femmodel->parameters->FindParam(CalvingMinIcebergSizeEnum);
	    config.critical_stress = femmodel->parameters->FindParam(CalvingCriticalStressEnum);

      config.timescale_yr = 0.0;
    if (femmodel->parameters->Exist(CrevasseDepthCalvingTimescaleEnum))
      femmodel->parameters->FindParam(&config.timescale_yr,
                                      CrevasseDepthCalvingTimescaleEnum);
    if (femmodel->parameters->Exist(CalvingAdvectIcefrontEnum))
      femmodel->parameters->FindParam(&config.advect_icefront,
                                      CalvingAdvectIcefrontEnum);
    else
      config.advect_icefront = true;
    if (femmodel->parameters->Exist(CalvingPropagateFromFrontEnum))
      femmodel->parameters->FindParam(&config.propagate_from_front,
                                      CalvingPropagateFromFrontEnum);
    else
      config.propagate_from_front = true;
    if (femmodel->parameters->Exist(CalvingRemoveOnlyMarkedEnum))
      femmodel->parameters->FindParam(&config.remove_only_marked,
                                      CalvingRemoveOnlyMarkedEnum);
    else
      config.remove_only_marked = false;
    if (femmodel->parameters->Exist(CalvingHydrofractureStabilizationEnum))
      femmodel->parameters->FindParam(&config.hydrofracture_stabilization,
                                      CalvingHydrofractureStabilizationEnum);
    else
      config.hydrofracture_stabilization = false;
    if (femmodel->parameters->Exist(CalvingHydrofractureWeakeningFactorEnum))
      femmodel->parameters->FindParam(&config.hydrofracture_weakening_factor,
                                      CalvingHydrofractureWeakeningFactorEnum);
    else
      config.hydrofracture_weakening_factor = 1.e-3;
    if (femmodel->parameters->Exist(CalvingHydrofractureMinOceanLevelsetEnum))
      femmodel->parameters->FindParam(&config.hydrofracture_min_ocean_levelset,
                                      CalvingHydrofractureMinOceanLevelsetEnum);
    else
      config.hydrofracture_min_ocean_levelset = 0.0;

    IssmDouble last_calving_time = 0.0;
    if (femmodel->parameters->Exist(LastCalvingTimeEnum)) {
      femmodel->parameters->FindParam(&last_calving_time, LastCalvingTimeEnum);
      if (last_calving_time < start_time_yr) {
        last_calving_time = start_time_yr;
        femmodel->parameters->SetParam(last_calving_time, LastCalvingTimeEnum);
      }
    } else {
      last_calving_time = start_time_yr;
      femmodel->parameters->AddObject(new DoubleParam(LastCalvingTimeEnum,
                                                      last_calving_time));
    }

    bool perform_calving = false;
    IssmDouble scheduled_calving_time = config.time_yr;
    if (!config.advect_icefront)
      FreezeIceFrontAdvectionIfNeeded(femmodel);

    _printf0_("   Crevasse Depth Calving"
	              << " (timescale=" << config.timescale_yr << "yrs, "
	              << "min. iceberg size=" << config.min_iceberg_size / 1000
	              << "km, critical stress=" << config.critical_stress / 1000
	              << "kPa, propagate from front=" << config.propagate_from_front
	              << ", remove only marked=" << config.remove_only_marked
	              << ", stabilization=" << config.hydrofracture_stabilization
	              << ", weakening factor="
	              << config.hydrofracture_weakening_factor
	              << ", min ocean levelset="
	              << config.hydrofracture_min_ocean_levelset
	              << ")"
	              << ":\n");
    _printf0_("[LevelsetAnalysis] :: Using LEFM-style crevasse propagation with locally computed r_xx and prescribed critical stress.\n");
    _printf0_("\tstart time=" << start_time_yr << " yr\n");
    _printf0_("\tlast check time=" << last_calving_time << " yr\n");
    _printf0_("\tcurrent time=" << config.time_yr << " yr\n");
    _printf0_("\ttime step=" << config.dt_yr << " yr\n");

    if (config.timescale_yr <= 0.0) {
      perform_calving = true;
    } else {
      IssmDouble elapsed_since_check = config.time_yr - last_calving_time;
      IssmDouble schedule_tolerance =
          1.e-12 * fmax(1.0, fabs(config.time_yr));
      if (elapsed_since_check + schedule_tolerance >= config.timescale_yr) {
        IssmDouble intervals_elapsed =
            floor((elapsed_since_check + schedule_tolerance) /
                  config.timescale_yr);
        intervals_elapsed = fmax(1.0, intervals_elapsed);
        scheduled_calving_time =
            last_calving_time + intervals_elapsed * config.timescale_yr;
        perform_calving = true;
      }
    }

    if (!perform_calving) {
      CrevasseNodalFields nodal_fields;
      PrepareHydrofractureGroundinglineDistance(femmodel);
      BuildCrevasseNodalFields(femmodel, &nodal_fields);
      StoreHydrofracturePredictionInput(femmodel, config, nodal_fields);
      CleanupCrevasseNodalFields(&nodal_fields);

      if (config.advect_icefront)
        ClearDynamicLevelsetConstraintsPreservingStaticSpc(femmodel);
      return;
    }

    femmodel->parameters->SetParam(scheduled_calving_time,
                                   LastCalvingTimeEnum);

    PrepareHydrofractureGroundinglineDistance(femmodel);
    InputDuplicatex(femmodel, MaskIceLevelsetEnum, DistanceToCalvingfrontEnum);
    femmodel->DistanceToFieldValue(MaskIceLevelsetEnum, 0,
                                   DistanceToCalvingfrontEnum);

    CrevasseNodalFields nodal_fields;
    BuildCrevasseNodalFields(femmodel, &nodal_fields);
    StoreHydrofracturePredictionInput(femmodel, config, nodal_fields);

    CrevassePropagationResult propagation;
    if (config.propagate_from_front) {
      PropagateCriticalCrevasseRegion(femmodel, config, nodal_fields,
                                      &propagation);
    } else {
      MarkCriticalMeltwaterNodes(femmodel, config, nodal_fields,
                                 &propagation);
    }
    ComputePropagatedRegionDiagnostics(femmodel, nodal_fields, &propagation);
    if (config.propagate_from_front && !config.remove_only_marked)
      BuildAppliedCalvingMask(femmodel, nodal_fields, &propagation);
    else if (config.propagate_from_front && config.remove_only_marked) {
      int localsize = femmodel->nodes->NumberOfNodesLocalAll();
      xDelete<IssmDouble>(propagation.applied_nodes);
      propagation.applied_nodes = xNew<IssmDouble>(localsize);
      for (int i = 0; i < localsize; i++)
        propagation.applied_nodes[i] = propagation.propagated_nodes[i];
    }
    CrevasseStressDiagnostics stress_diagnostics =
        ComputeCrevasseStressDiagnostics(femmodel, nodal_fields, propagation);

    bool has_front_extent =
        nodal_fields.ice_front_x < std::numeric_limits<IssmDouble>::max() &&
        propagation.propagated_min_x < std::numeric_limits<IssmDouble>::max();
    bool meets_size_threshold =
        !config.propagate_from_front ||
        (has_front_extent &&
         nodal_fields.ice_front_x > propagation.propagated_min_x &&
         (nodal_fields.ice_front_x - propagation.propagated_min_x >=
          config.min_iceberg_size));

    femmodel->parameters->SetParam(propagation.propagated_min_x,
                                   CalvingPropagatedMinXEnum);
    femmodel->parameters->SetParam(propagation.propagated_area,
                                   CalvingPropagatedAreaEnum);

    int total_nodes = femmodel->nodes->NumberOfNodes();
    int marked_nodes =
        meets_size_threshold
            ? CountMarkedMasterNodes(femmodel, propagation.applied_nodes)
            : 0;
    IssmDouble removed_area =
        meets_size_threshold ? propagation.propagated_area : 0.0;

    if (nodal_fields.ice_front_x < std::numeric_limits<IssmDouble>::max())
      _printf0_("\tice front          = "
                << nodal_fields.ice_front_x / 1000 << "km.\n");
    else
      _printf0_("\tice front          = n/a\n");

    if (propagation.propagated_min_x < std::numeric_limits<IssmDouble>::max())
      _printf0_("\tpropagated min x   = "
                << propagation.propagated_min_x / 1000
                << "km\n");
    else
      _printf0_("\tpropagated min x   = n/a\n");
    _printf0_("\tstress-connected region = "
              << propagation.has_critical_region << "\n");
    _printf0_("\tsize criterion met = " << meets_size_threshold << "\n");
    _printf0_("\tmarked nodes       = " << marked_nodes << "/"
              << total_nodes << "\n");
    _printf0_("\tremoved area       = " << removed_area / 1.e6
              << "km^2.\n");
    PrintCrevasseStressSummary("seed nodes", stress_diagnostics.seed);
    PrintCrevasseStressSummary("floating nodes", stress_diagnostics.floating);
    PrintCrevasseStressSummary("propagated", stress_diagnostics.propagated);
    _printf0_("\tcalving occurred   = "
              << (propagation.has_critical_region && meets_size_threshold)
              << "\n");

    if (propagation.has_critical_region && meets_size_threshold) {
      if (config.hydrofracture_stabilization) {
        StoreHydrofractureWeakIceInput(femmodel, propagation.applied_nodes);
        if (config.advect_icefront)
          ClearDynamicLevelsetConstraintsPreservingStaticSpc(femmodel);
        femmodel->parameters->SetParam(0.0, CalvingOccurredEnum);
        _printf0_("\tRegularized hydrofracture as weak ice; no level-set removal executed.\n");
      } else {
        StoreHydrofractureWeakIceInput(femmodel, NULL);
        ApplyPropagatedRegionConstraints(femmodel, config, nodal_fields,
                                         propagation);
        femmodel->parameters->SetParam(1.0, CalvingOccurredEnum);
        _printf0_("\tExecuted calving.\n");
      }
    } else {
      StoreHydrofractureWeakIceInput(femmodel, NULL);
      if (config.advect_icefront)
        ClearDynamicLevelsetConstraintsPreservingStaticSpc(femmodel);
      femmodel->parameters->SetParam(0.0, CalvingOccurredEnum);
      _printf0_("\tDid not execute calving.\n");
    }

    CleanupCrevasseNodalFields(&nodal_fields);
    xDelete<IssmDouble>(propagation.applied_nodes);
    xDelete<IssmDouble>(propagation.propagated_nodes);
  }

  /*Default, do nothing*/
  return;
}
