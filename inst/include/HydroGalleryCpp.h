#ifndef HYDROGALLERYCPP_API_H
#define HYDROGALLERYCPP_API_H

#include <armadillo>
#include <cstdint>
#include <functional>

// Global function declarations from HydroGalleryCpp core
arma::vec atmosSnow_ThresholdT(
    const arma::vec& ATMOS_precipitation_mm,
    const arma::vec& ATMOS_temperature_Cel,
    const arma::vec& param_ATMOS_thr_Ts);

arma::vec intercep_Full(
    const arma::vec& ATMOS_precipitation_mm,
    const arma::vec& LAND_interceptWater_mm,
    const arma::vec& LAND_interceptCapacity_mm);

arma::vec meteo_nettoRadiat_WaterGAP3(
    const arma::vec& ATMOS_solarRadiat_MJ,
    const arma::vec& ATMOS_nettoLongRadiat_MJ,
    const arma::vec& LAND_albedo_1);

arma::vec evatransPotential_PriestleyTaylor(
    const arma::vec& ATMOS_temperature_Cel,
    const arma::vec& ATMOS_netRadiat_MJ,
    const arma::vec& LAND_elevation_m,
    const arma::vec& param_EVATRANS_prt_alpha);

arma::vec evatransActual_SupplyPow(
    const arma::vec& ATMOS_potentialEvatrans_mm,
    const arma::vec& water_mm,
    const arma::vec& capacity_mm,
    const arma::vec& param_EVATRANS_sup_k,
    const arma::vec& param_EVATRANS_sup_gamma);

arma::vec evatransActual_VIC(
    const arma::vec& ATMOS_potentialEvatrans_mm,
    const arma::vec& water_mm,
    const arma::vec& capacity_mm,
    const arma::vec& param_EVATRANS_vic_gamma);

arma::vec evatransActual_WaterGAP3(
    const arma::vec& ATMOS_potentialEvatrans_mm,
    const arma::vec& water_mm,
    const arma::vec& capacity_mm,
    const arma::vec& param_EVATRANS_wat_petmax);

arma::vec snowMelt_Factor(
    const arma::vec& SNOW_ice_mm,
    const arma::vec& ATMOS_temperature_Cel,
    const arma::vec& param_SNOWMELT_fac_f,
    const arma::vec& param_SNOWMELT_fac_Tmelt);

arma::vec infilt_HBV(
    const arma::vec& LAND_water_mm,
    const arma::vec& SOIL_water_mm,
    const arma::vec& SOIL_capacity_mm,
    const arma::vec& param_INFILT_hbv_beta);

arma::vec percola_WaterGAP3(
    const arma::vec& LAND_water_mm,
    const arma::vec& SOIL_potentialPercola_mm,
    const arma::uvec& param_PERCOLA_wat_01,
    const arma::vec& param_PERCOLA_wat_thresh,
    const arma::vec& param_PERCOLA_wat_k);

arma::vec baseflow_SupplyRatio(
    const arma::vec& GROUND_water_mm,
    const arma::vec& param_BASEFLOW_sur_k);

arma::mat landLeafAreaIndex_WaterGAP3(
    const arma::mat& ATMOS_temperature_Cel,
    const arma::mat& ATMOS_precipitation_mm,
    const arma::vec& CELL_latitude_deg,
    const arma::uvec& LAND_growUpDay_d,
    const arma::vec& LAND_leafAreaIndexMin_,
    const arma::vec& LAND_leafAreaIndexMax_,
    const arma::uvec& Time_dayOfYear_d);

arma::vec lakeout_SupplyPow(
    const arma::vec& Lake_water_m3,
    const arma::vec& Lake_capacity_m3,
    const arma::vec& param_Lakeout_sup_storeFactor,
    const arma::vec& param_Lakeout_sup_gamma);

arma::vec riverout_LinearResorvoir(
    const arma::vec& RIVER_water_m3,
    const arma::vec& RIVER_inflow_m3,
    const arma::vec& RIVER_velocity_km,
    const arma::vec& RIVER_length_km);

arma::vec riverlakout_LinearResorvoir(
    const arma::vec& Riverlak_water_m3,
    const arma::vec& Riverlak_inflow_m3,
    const arma::vec& Riverlak_capacity_m3,
    const arma::vec& param_Riverlak_lin_storeFactor);

arma::vec reservoiReleas_Hanasaki(
    const arma::vec& Reservoi_water_m3,
    const arma::vec& Reservoi_inflow_m3,
    const arma::vec& Reservoi_demand_m3,
    const arma::vec& Reservoi_capacity_m3,
    const arma::vec& Reservoi_meanInflow_m3,
    const arma::vec& Reservoi_meanDemand_m3,
    const arma::uvec& Reservoi_isIrrigate_01);

void withdraw_SingleCell(
    arma::vec& CELL_withdrawal_m3,
    arma::vec& CELL_water_m3);

void withdrawSurface_Around(
    arma::vec& CELL_withdrawal_m3,
    arma::vec& RIVER_water_m3,
    const arma::uvec& Lake_cellNumber_int,
    arma::vec& Lake_water_m3,
    const arma::umat& CELL_cellNumberAround_int);

struct DDSOptions {
  int max_iter = 100;
  double r = 0.2;
  std::uint64_t seed = 12345ULL;
  bool verbose = false;
};

struct DDSResult {
  arma::vec x_best;
  double y_best;
  int evaluations;
};

DDSResult cali_DDS(
    const std::function<double(const arma::vec&)>& objective,
    const arma::vec& x_Min,
    const arma::vec& x_Max,
    const arma::vec& x_Init,
    const DDSOptions& opt);

// Backward-compatible namespace expected by existing WaterGAP3Cpp code
namespace HydroGallery {
using ::atmosSnow_ThresholdT;
using ::intercep_Full;
using ::meteo_nettoRadiat_WaterGAP3;
using ::evatransPotential_PriestleyTaylor;
using ::evatransActual_SupplyPow;
using ::evatransActual_VIC;
using ::evatransActual_WaterGAP3;
using ::snowMelt_Factor;
using ::infilt_HBV;
using ::percola_WaterGAP3;
using ::baseflow_SupplyRatio;
using ::landLeafAreaIndex_WaterGAP3;
using ::lakeout_SupplyPow;
using ::riverout_LinearResorvoir;
using ::riverlakout_LinearResorvoir;
using ::reservoiReleas_Hanasaki;
using ::withdraw_SingleCell;
using ::withdrawSurface_Around;
}  // namespace HydroGallery

#endif  // HYDROGALLERYCPP_API_H
