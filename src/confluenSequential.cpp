#include "confluen.h"
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::interfaces(r, cpp)]]


arma::uvec find_in(const arma::uvec& x, const arma::uvec& y) {
  std::vector<arma::uword> matches;
  for (arma::uword i = 0; i < x.n_elem; ++i) {
    if (arma::any(y == x(i))) {
      matches.push_back(i);
    }
  }
  return arma::uvec(matches);
}

arma::vec inflow_add(const arma::vec& num_Outflow_LastStep, const arma::umat& int_InflowCell) {
  arma::vec num_Inflow_m3(int_InflowCell.n_rows, arma::fill::zeros);  // One inflow per cell (per row)
  
  for (arma::uword i = 0; i < int_InflowCell.n_rows; ++i) {
    arma::urowvec inflow_cells = int_InflowCell.row(i);  // Get inflow cell indices for current cell
    arma::uvec valid_idx = arma::find(inflow_cells > 0);  // Only consider >0 (valid) inflow cell indices
    if (!valid_idx.empty()) {
      // Convert to zero-based indices
      arma::uvec inflow_indices = arma::conv_to<arma::uvec>::from(inflow_cells.elem(valid_idx)) - 1;
      // Sum the outflows from these inflow cells
      num_Inflow_m3(i) = arma::sum(num_Outflow_LastStep.elem(inflow_indices));
    }
  }
  
  return num_Inflow_m3;
}

//' confluenSequential
//' @name confluenSequential
//' @export
// [[Rcpp::export]]
arma::vec confluen_WaterGAP3U(
    arma::vec& RIVER_water_m3,
    const arma::vec& RIVER_length_km,
    const arma::vec& RIVER_velocity_km,
    const arma::vec& RIVER_inflow_m3,
    const arma::field<arma::uvec>& CELL_cellNumberStep_int,  
    const arma::field<arma::umat>& CELL_inflowCellNumberStep_int,  
    const arma::uvec& Riverlak_cellNumber_int,
    const arma::vec& Riverlak_capacity_m3,
    const arma::uvec& Reservoi_cellNumber_int,
    const arma::vec& Reservoi_demand_m3,
    const arma::vec& Reservoi_capacity_m3,
    const arma::vec& Reservoi_meanInflow_m3,
    const arma::vec& Reservoi_meanDemand_m3,
    const arma::uvec& Reservoi_isIrrigate_01,
    const arma::vec& param_Riverlak_lin_storeFactor)
{
  arma::vec RIVER_water_m3_TEMP = RIVER_water_m3;
  arma::vec RIVER_outflow_m3(RIVER_inflow_m3.n_elem, arma::fill::zeros);

  arma::vec Riverlak_water_m3 = RIVER_water_m3_TEMP.elem(Riverlak_cellNumber_int - 1);
  arma::vec Reservoi_water_m3 = RIVER_water_m3_TEMP.elem(Reservoi_cellNumber_int - 1);
  arma::vec Riverlak_inflow_m3 = RIVER_inflow_m3.elem(Riverlak_cellNumber_int - 1);
  arma::vec Reservoi_inflow_m3 = RIVER_inflow_m3.elem(Reservoi_cellNumber_int - 1);

  const int n_Step = CELL_cellNumberStep_int.n_elem;

  arma::uvec idx_Cell_Step = CELL_cellNumberStep_int(0) - 1;

  arma::vec step_RiverOutflow_m3 = riverout_LinearResorvoir(
    RIVER_water_m3_TEMP.elem(idx_Cell_Step),
    RIVER_inflow_m3.elem(idx_Cell_Step),
    RIVER_velocity_km.elem(idx_Cell_Step),
    RIVER_length_km.elem(idx_Cell_Step)
  );

  

  arma::vec step_RIVER_Water_New = RIVER_water_m3_TEMP.elem(idx_Cell_Step) +
    RIVER_inflow_m3.elem(idx_Cell_Step) -
    step_RiverOutflow_m3;

  arma::uvec idx_Riverlak_Step = find_in(Riverlak_cellNumber_int, idx_Cell_Step + 1);
  if (!idx_Riverlak_Step.is_empty()) {
    arma::uvec idx_Step_Riverlak = find_in(idx_Cell_Step + 1, Riverlak_cellNumber_int);

    arma::vec step_RiverlakOutflow_m3 = riverlakout_LinearResorvoir(
      Riverlak_water_m3.elem(idx_Riverlak_Step),
      Riverlak_inflow_m3.elem(idx_Riverlak_Step),
      Riverlak_capacity_m3.elem(idx_Riverlak_Step),
      param_Riverlak_lin_storeFactor.elem(idx_Riverlak_Step)
    );

    step_RiverOutflow_m3.elem(idx_Step_Riverlak) = step_RiverlakOutflow_m3;
    step_RIVER_Water_New.elem(idx_Step_Riverlak) =
      Riverlak_water_m3.elem(idx_Riverlak_Step) +
      RIVER_inflow_m3.elem(idx_Cell_Step.elem(idx_Step_Riverlak)) -
      step_RiverlakOutflow_m3;
  }

  arma::uvec idx_Reservoi_Step = find_in(Reservoi_cellNumber_int, idx_Cell_Step + 1);
  if (!idx_Reservoi_Step.is_empty()) {
    arma::uvec idx_Step_Reservoi = find_in(idx_Cell_Step + 1, Reservoi_cellNumber_int);

    arma::vec step_ReservoiOutflow_m3 = reservoiReleas_Hanasaki(
      Reservoi_water_m3.elem(idx_Reservoi_Step),
      Reservoi_inflow_m3.elem(idx_Reservoi_Step),
      Reservoi_demand_m3.elem(idx_Reservoi_Step),
      Reservoi_capacity_m3.elem(idx_Reservoi_Step),
      Reservoi_meanInflow_m3.elem(idx_Reservoi_Step),
      Reservoi_meanDemand_m3.elem(idx_Reservoi_Step),
      Reservoi_isIrrigate_01.elem(idx_Reservoi_Step)
    );

    step_RiverOutflow_m3.elem(idx_Step_Reservoi) = step_ReservoiOutflow_m3;
    step_RIVER_Water_New.elem(idx_Step_Reservoi) =
      Reservoi_water_m3.elem(idx_Reservoi_Step) +
      Reservoi_inflow_m3.elem(idx_Reservoi_Step) -
      step_ReservoiOutflow_m3;
  }

  RIVER_outflow_m3.elem(idx_Cell_Step) = step_RiverOutflow_m3;
  RIVER_water_m3_TEMP.elem(idx_Cell_Step) = step_RIVER_Water_New;

  for (int i_Step = 1; i_Step < n_Step; ++i_Step) {
    idx_Cell_Step = CELL_cellNumberStep_int(i_Step) - 1;

    arma::vec step_UpstreamInflow_m3 = inflow_add(
      RIVER_outflow_m3,
      CELL_inflowCellNumberStep_int(i_Step)
    );
    
    step_UpstreamInflow_m3 += RIVER_inflow_m3.elem(idx_Cell_Step);

    step_RiverOutflow_m3 = riverout_LinearResorvoir(
      RIVER_water_m3_TEMP.elem(idx_Cell_Step),
      step_UpstreamInflow_m3,
      RIVER_velocity_km.elem(idx_Cell_Step),
      RIVER_length_km.elem(idx_Cell_Step)
    );

    step_RIVER_Water_New = RIVER_water_m3_TEMP.elem(idx_Cell_Step) +
      step_UpstreamInflow_m3 -
      step_RiverOutflow_m3;

    idx_Riverlak_Step = find_in(Riverlak_cellNumber_int, idx_Cell_Step + 1);
    if (!idx_Riverlak_Step.is_empty()) {
      arma::uvec idx_Step_Riverlak = find_in(idx_Cell_Step + 1, Riverlak_cellNumber_int);


      
      
      arma::vec step_RiverlakOutflow_m3 = riverlakout_LinearResorvoir(
        Riverlak_water_m3.elem(idx_Riverlak_Step),
        step_UpstreamInflow_m3.elem(idx_Step_Riverlak),
        Riverlak_capacity_m3.elem(idx_Riverlak_Step),
        param_Riverlak_lin_storeFactor.elem(idx_Riverlak_Step)
      );

      step_RiverOutflow_m3.elem(idx_Step_Riverlak) = step_RiverlakOutflow_m3;
      step_RIVER_Water_New.elem(idx_Step_Riverlak) =
        Riverlak_water_m3.elem(idx_Riverlak_Step) +
        step_UpstreamInflow_m3.elem(idx_Step_Riverlak) -
        step_RiverlakOutflow_m3;
    }

    idx_Reservoi_Step = find_in(Reservoi_cellNumber_int, idx_Cell_Step + 1);
    if (!idx_Reservoi_Step.is_empty()) {
      arma::uvec idx_Step_Reservoi = find_in(idx_Cell_Step + 1, Reservoi_cellNumber_int);

      arma::vec step_ReservoiOutflow_m3 = reservoiReleas_Hanasaki(
        Reservoi_water_m3.elem(idx_Reservoi_Step),
        step_UpstreamInflow_m3.elem(idx_Step_Reservoi),
        Reservoi_demand_m3.elem(idx_Reservoi_Step),
        Reservoi_capacity_m3.elem(idx_Reservoi_Step),
        Reservoi_meanInflow_m3.elem(idx_Reservoi_Step),
        Reservoi_meanDemand_m3.elem(idx_Reservoi_Step),
        Reservoi_isIrrigate_01.elem(idx_Reservoi_Step)
      );

      step_RiverOutflow_m3.elem(idx_Step_Reservoi) = step_ReservoiOutflow_m3;
      step_RIVER_Water_New.elem(idx_Step_Reservoi) =
        Reservoi_water_m3.elem(idx_Reservoi_Step) +
        step_UpstreamInflow_m3.elem(idx_Step_Reservoi) -
        step_ReservoiOutflow_m3;
    }

    RIVER_outflow_m3.elem(idx_Cell_Step) = step_RiverOutflow_m3;
    RIVER_water_m3_TEMP.elem(idx_Cell_Step) = step_RIVER_Water_New;
  }

  RIVER_water_m3 = RIVER_water_m3_TEMP;
  return RIVER_outflow_m3;
}

//' @name confluenSequential
//' @export
// [[Rcpp::export]]
arma::vec confluen_WaterGAP3N(
    arma::vec& RIVER_water_m3,
    const arma::vec& RIVER_length_km,
    const arma::vec& RIVER_velocity_km,
    const arma::vec& RIVER_inflow_m3,
    const arma::field<arma::uvec>& CELL_cellNumberStep_int,  
    const arma::field<arma::umat>& CELL_inflowCellNumberStep_int,  
    const arma::uvec& Riverlak_cellNumber_int,
    const arma::vec& Riverlak_capacity_m3,
    const arma::vec& param_Riverlak_lin_storeFactor)
{
  arma::vec RIVER_water_m3_TEMP = RIVER_water_m3;
  arma::vec RIVER_outflow_m3(RIVER_inflow_m3.n_elem, arma::fill::zeros);

  arma::vec Riverlak_water_m3 = RIVER_water_m3_TEMP.elem(Riverlak_cellNumber_int - 1);
  arma::vec Riverlak_inflow_m3 = RIVER_inflow_m3.elem(Riverlak_cellNumber_int - 1);

  const int n_Step = CELL_cellNumberStep_int.n_elem;

  arma::uvec idx_Cell_Step = CELL_cellNumberStep_int(0) - 1;

  arma::vec step_RiverOutflow_m3 = riverout_LinearResorvoir(
    RIVER_water_m3_TEMP.elem(idx_Cell_Step),
    RIVER_inflow_m3.elem(idx_Cell_Step),
    RIVER_velocity_km.elem(idx_Cell_Step),
    RIVER_length_km.elem(idx_Cell_Step)
  );

  

  arma::vec step_RIVER_Water_New = RIVER_water_m3_TEMP.elem(idx_Cell_Step) +
    RIVER_inflow_m3.elem(idx_Cell_Step) -
    step_RiverOutflow_m3;

  arma::uvec idx_Riverlak_Step = find_in(Riverlak_cellNumber_int, idx_Cell_Step + 1);
  if (!idx_Riverlak_Step.is_empty()) {
    arma::uvec idx_Step_Riverlak = find_in(idx_Cell_Step + 1, Riverlak_cellNumber_int);

    arma::vec step_RiverlakOutflow_m3 = riverlakout_LinearResorvoir(
      Riverlak_water_m3.elem(idx_Riverlak_Step),
      Riverlak_inflow_m3.elem(idx_Riverlak_Step),
      Riverlak_capacity_m3.elem(idx_Riverlak_Step),
      param_Riverlak_lin_storeFactor.elem(idx_Riverlak_Step)
    );

    step_RiverOutflow_m3.elem(idx_Step_Riverlak) = step_RiverlakOutflow_m3;
    step_RIVER_Water_New.elem(idx_Step_Riverlak) =
      Riverlak_water_m3.elem(idx_Riverlak_Step) +
      RIVER_inflow_m3.elem(idx_Cell_Step.elem(idx_Step_Riverlak)) -
      step_RiverlakOutflow_m3;
  }


  RIVER_outflow_m3.elem(idx_Cell_Step) = step_RiverOutflow_m3;
  RIVER_water_m3_TEMP.elem(idx_Cell_Step) = step_RIVER_Water_New;

  for (int i_Step = 1; i_Step < n_Step; ++i_Step) {
    idx_Cell_Step = CELL_cellNumberStep_int(i_Step) - 1;

    arma::vec step_UpstreamInflow_m3 = inflow_add(
      RIVER_outflow_m3,
      CELL_inflowCellNumberStep_int(i_Step)
    );
    
    step_UpstreamInflow_m3 += RIVER_inflow_m3.elem(idx_Cell_Step);

    step_RiverOutflow_m3 = riverout_LinearResorvoir(
      RIVER_water_m3_TEMP.elem(idx_Cell_Step),
      step_UpstreamInflow_m3,
      RIVER_velocity_km.elem(idx_Cell_Step),
      RIVER_length_km.elem(idx_Cell_Step)
    );

    step_RIVER_Water_New = RIVER_water_m3_TEMP.elem(idx_Cell_Step) +
      step_UpstreamInflow_m3 -
      step_RiverOutflow_m3;

    idx_Riverlak_Step = find_in(Riverlak_cellNumber_int, idx_Cell_Step + 1);
    if (!idx_Riverlak_Step.is_empty()) {
      arma::uvec idx_Step_Riverlak = find_in(idx_Cell_Step + 1, Riverlak_cellNumber_int);


      
      
      arma::vec step_RiverlakOutflow_m3 = riverlakout_LinearResorvoir(
        Riverlak_water_m3.elem(idx_Riverlak_Step),
        step_UpstreamInflow_m3.elem(idx_Step_Riverlak),
        Riverlak_capacity_m3.elem(idx_Riverlak_Step),
        param_Riverlak_lin_storeFactor.elem(idx_Riverlak_Step)
      );

      step_RiverOutflow_m3.elem(idx_Step_Riverlak) = step_RiverlakOutflow_m3;
      step_RIVER_Water_New.elem(idx_Step_Riverlak) =
        Riverlak_water_m3.elem(idx_Riverlak_Step) +
        step_UpstreamInflow_m3.elem(idx_Step_Riverlak) -
        step_RiverlakOutflow_m3;
    }


    RIVER_outflow_m3.elem(idx_Cell_Step) = step_RiverOutflow_m3;
    RIVER_water_m3_TEMP.elem(idx_Cell_Step) = step_RIVER_Water_New;
  }

  RIVER_water_m3 = RIVER_water_m3_TEMP;
  return RIVER_outflow_m3;
}
