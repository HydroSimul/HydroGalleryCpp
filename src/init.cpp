#include <R_ext/Rdynload.h>
#include "../inst/include/HydroGalleryCpp.h"

extern "C" void R_init_HydroGalleryCpp(DllInfo* dll) {
  R_useDynamicSymbols(dll, FALSE);

  R_RegisterCCallable("HydroGalleryCpp", "atmosSnow_ThresholdT",
                      (DL_FUNC)&atmosSnow_ThresholdT);
  R_RegisterCCallable("HydroGalleryCpp", "intercep_Full",
                      (DL_FUNC)&intercep_Full);
  R_RegisterCCallable("HydroGalleryCpp", "meteo_nettoRadiat_WaterGAP3",
                      (DL_FUNC)&meteo_nettoRadiat_WaterGAP3);
  R_RegisterCCallable("HydroGalleryCpp", "evatransPotential_PriestleyTaylor",
                      (DL_FUNC)&evatransPotential_PriestleyTaylor);
  R_RegisterCCallable("HydroGalleryCpp", "evatransActual_SupplyPow",
                      (DL_FUNC)&evatransActual_SupplyPow);
  R_RegisterCCallable("HydroGalleryCpp", "evatransActual_VIC",
                      (DL_FUNC)&evatransActual_VIC);
  R_RegisterCCallable("HydroGalleryCpp", "evatransActual_WaterGAP3",
                      (DL_FUNC)&evatransActual_WaterGAP3);
  R_RegisterCCallable("HydroGalleryCpp", "snowMelt_Factor",
                      (DL_FUNC)&snowMelt_Factor);
  R_RegisterCCallable("HydroGalleryCpp", "infilt_HBV",
                      (DL_FUNC)&infilt_HBV);
  R_RegisterCCallable("HydroGalleryCpp", "percola_WaterGAP3",
                      (DL_FUNC)&percola_WaterGAP3);
  R_RegisterCCallable("HydroGalleryCpp", "baseflow_SupplyRatio",
                      (DL_FUNC)&baseflow_SupplyRatio);
  R_RegisterCCallable("HydroGalleryCpp", "landLeafAreaIndex_WaterGAP3",
                      (DL_FUNC)&landLeafAreaIndex_WaterGAP3);
  R_RegisterCCallable("HydroGalleryCpp", "lakeout_SupplyPow",
                      (DL_FUNC)&lakeout_SupplyPow);
  R_RegisterCCallable("HydroGalleryCpp", "riverout_LinearResorvoir",
                      (DL_FUNC)&riverout_LinearResorvoir);
  R_RegisterCCallable("HydroGalleryCpp", "riverlakout_LinearResorvoir",
                      (DL_FUNC)&riverlakout_LinearResorvoir);
  R_RegisterCCallable("HydroGalleryCpp", "reservoiReleas_Hanasaki",
                      (DL_FUNC)&reservoiReleas_Hanasaki);
  R_RegisterCCallable("HydroGalleryCpp", "withdraw_SingleCell",
                      (DL_FUNC)&withdraw_SingleCell);
  R_RegisterCCallable("HydroGalleryCpp", "withdrawSurface_Around",
                      (DL_FUNC)&withdrawSurface_Around);
}
