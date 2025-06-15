#include "confluen.h"
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::interfaces(r, cpp)]]



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
    const arma::field<arma::uvec>& NET_cellNumberStep_int,  
    const arma::field<arma::umat>& NET_upstreamCellNumberStep_int,  
    const arma::field<arma::uvec>& NET_riverlakNumberStep_int,  
    const arma::field<arma::uvec>& NET_reservoiNumberStep_int,  
    const arma::uvec& Riverlak_cellNumber_int,
    arma::vec& Riverlak_water_m3,
    const arma::vec& Riverlak_capacity_m3,
    const arma::uvec& Reservoi_cellNumber_int,
    arma::vec& Reservoi_water_m3,
    const arma::vec& Reservoi_demand_m3,
    const arma::vec& Reservoi_capacity_m3,
    const arma::vec& Reservoi_meanInflow_m3,
    const arma::vec& Reservoi_meanDemand_m3,
    const arma::uvec& Reservoi_isIrrigate_01,
    const arma::vec& param_Riverlak_lin_storeFactor)
{
  arma::vec RIVER_outflow_m3(RIVER_inflow_m3.n_elem, arma::fill::zeros);
  
  
  // step loop ----------
  const int n_Step = NET_cellNumberStep_int.n_elem;
  for (int i_Step = 0; i_Step < n_Step; ++i_Step) {
    
    arma::uvec idx_CELL_Step = NET_cellNumberStep_int(i_Step) - 1;
    // inflow
    arma::vec step_RiverInflow_m3 = RIVER_inflow_m3.elem(idx_CELL_Step);
    if (i_Step > 0) {
      step_RiverInflow_m3 += inflow_add(
        RIVER_outflow_m3,
        NET_upstreamCellNumberStep_int(i_Step)
      );
    }
    
    RIVER_outflow_m3.elem(idx_CELL_Step) = riverout_LinearResorvoir(
      RIVER_water_m3.elem(idx_CELL_Step),
      step_RiverInflow_m3,
      RIVER_velocity_km.elem(idx_CELL_Step),
      RIVER_length_km.elem(idx_CELL_Step)
    );
    RIVER_water_m3.elem(idx_CELL_Step) += step_RiverInflow_m3 - RIVER_outflow_m3.elem(idx_CELL_Step);
    
    // Riverlak
    if (!NET_riverlakNumberStep_int(i_Step).is_empty()) {
      arma::uvec idx_Riverlak_Step = NET_riverlakNumberStep_int(i_Step) - 1;
      arma::uvec idx_CELL_Riverlak_Step = Riverlak_cellNumber_int.elem(idx_Riverlak_Step) - 1;
      
      arma::vec step_RiverlakOutflow_m3 = riverlakout_LinearResorvoir(
        Riverlak_water_m3.elem(idx_Riverlak_Step),
        RIVER_outflow_m3.elem(idx_CELL_Riverlak_Step),
        Riverlak_capacity_m3.elem(idx_Riverlak_Step),
        param_Riverlak_lin_storeFactor.elem(idx_Riverlak_Step)
      );
      Riverlak_water_m3.elem(idx_Riverlak_Step) += RIVER_outflow_m3.elem(idx_CELL_Riverlak_Step) - step_RiverlakOutflow_m3;
      RIVER_outflow_m3.elem(idx_CELL_Riverlak_Step) = step_RiverlakOutflow_m3;
    }
    
    // Riverlak
    if (!NET_reservoiNumberStep_int(i_Step).is_empty()) {
      arma::uvec idx_Reservoi_Step = NET_reservoiNumberStep_int(i_Step) - 1;
      arma::uvec idx_CELL_Reservoi_Step = Reservoi_cellNumber_int.elem(idx_Reservoi_Step) - 1;
      
      arma::vec step_ReservoiOutflow_m3 = reservoiReleas_Hanasaki(
        Reservoi_water_m3.elem(idx_Reservoi_Step),
        RIVER_outflow_m3.elem(idx_CELL_Reservoi_Step),
        Reservoi_demand_m3.elem(idx_Reservoi_Step),
        Reservoi_capacity_m3.elem(idx_Reservoi_Step),
        Reservoi_meanInflow_m3.elem(idx_Reservoi_Step),
        Reservoi_meanDemand_m3.elem(idx_Reservoi_Step),
        Reservoi_isIrrigate_01.elem(idx_Reservoi_Step)
      );
      Reservoi_water_m3.elem(idx_Reservoi_Step) += RIVER_outflow_m3.elem(idx_CELL_Reservoi_Step) - step_ReservoiOutflow_m3;
      RIVER_outflow_m3.elem(idx_CELL_Reservoi_Step) = step_ReservoiOutflow_m3;
    }
    
  }
  
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
    const arma::field<arma::uvec>& NET_cellNumberStep_int,  
    const arma::field<arma::umat>& NET_upstreamCellNumberStep_int,  
    const arma::field<arma::uvec>& NET_riverlakNumberStep_int,  
    const arma::uvec& Riverlak_cellNumber_int,
    arma::vec& Riverlak_water_m3,
    const arma::vec& Riverlak_capacity_m3,
    const arma::vec& param_Riverlak_lin_storeFactor)
{
  arma::vec RIVER_outflow_m3(RIVER_inflow_m3.n_elem, arma::fill::zeros);


  // step loop ----------
  const int n_Step = NET_cellNumberStep_int.n_elem;
  for (int i_Step = 0; i_Step < n_Step; ++i_Step) {
    
    arma::uvec idx_CELL_Step = NET_cellNumberStep_int(i_Step) - 1;
    // inflow
    arma::vec step_RiverInflow_m3 = RIVER_inflow_m3.elem(idx_CELL_Step);
    if (i_Step > 0) {
      step_RiverInflow_m3 += inflow_add(
        RIVER_outflow_m3,
        NET_upstreamCellNumberStep_int(i_Step)
      );
    }
    
    RIVER_outflow_m3.elem(idx_CELL_Step) = riverout_LinearResorvoir(
      RIVER_water_m3.elem(idx_CELL_Step),
      step_RiverInflow_m3,
      RIVER_velocity_km.elem(idx_CELL_Step),
      RIVER_length_km.elem(idx_CELL_Step)
    );
    RIVER_water_m3.elem(idx_CELL_Step) += step_RiverInflow_m3 - RIVER_outflow_m3.elem(idx_CELL_Step);
    
    // Riverlak
    if (!NET_riverlakNumberStep_int(i_Step).is_empty()) {
      arma::uvec idx_Riverlak_Step = NET_riverlakNumberStep_int(i_Step) - 1;
      arma::uvec idx_CELL_Riverlak_Step = Riverlak_cellNumber_int.elem(idx_Riverlak_Step) - 1;
      
      arma::vec step_RiverlakOutflow_m3 = riverlakout_LinearResorvoir(
        Riverlak_water_m3.elem(idx_Riverlak_Step),
        RIVER_outflow_m3.elem(idx_CELL_Riverlak_Step),
        Riverlak_capacity_m3.elem(idx_Riverlak_Step),
        param_Riverlak_lin_storeFactor.elem(idx_Riverlak_Step)
      );
      Riverlak_water_m3.elem(idx_Riverlak_Step) += RIVER_outflow_m3.elem(idx_CELL_Riverlak_Step) - step_RiverlakOutflow_m3;
      RIVER_outflow_m3.elem(idx_CELL_Riverlak_Step) = step_RiverlakOutflow_m3;
    }
    
  }
  
  return RIVER_outflow_m3;
}