#include "utils.h"

// Declaration from snow.cpp
arma::vec snowMelt_Factor(
    const arma::vec& SNOW_ice_mm,
    const arma::vec& ATMOS_temperature_Cel,
    const arma::vec& param_SNOWMELT_fac_f,
    const arma::vec& param_SNOWMELT_fac_Tmelt
);

//' **wetland snow process**
//' @name wetlandsnow
//' @inheritParams all_vari
//' @description
//' WaterGAP3 wetland freezing and melt partitioning.
//'
//' @param param_Wetland_kup_Tfrozen (Cel) threshold temperature for frozen-day increase
//' @param param_Wetland_kup_f (mm/day/Cel) melt factor parameter for [snowMelt_Factor()]
//' @param param_Wetland_kup_Tmelt (Cel) threshold temperature for melt start
//' @param param_Wetland_kup_MaxFrozenD (day) maximum frozen-day storage threshold
//' @return Wetland_precipittaion_mm (mm) precipitation reaching wetland water
//' @export
arma::vec wetlandsnow_Kupzig(
    const arma::vec& ATMOS_precipitation_mm,
    const arma::vec& ATMOS_temperature_Cel,
    arma::vec& Wetland_frozenDay_d,
    arma::vec& Wetland_iceStorage_mm,
    const arma::vec& param_Wetland_kup_Tfrozen,
    const arma::vec& param_Wetland_kup_f,
    const arma::vec& param_Wetland_kup_Tmelt,
    const arma::vec& param_Wetland_kup_MaxFrozenD
)
{
    // 1) Update frozen-day counter: +1 if below freezing, -1 otherwise, clamped [0, MaxFrozenD]
    arma::vec delta = 2.0 * arma::conv_to<arma::vec>::from(
        ATMOS_temperature_Cel < param_Wetland_kup_Tfrozen
    ) - 1.0;
    Wetland_frozenDay_d = arma::clamp(
        Wetland_frozenDay_d + delta,
        0.0, arma::datum::inf
    );
    Wetland_frozenDay_d = arma::min(Wetland_frozenDay_d, param_Wetland_kup_MaxFrozenD);

    // 2) Partition precipitation: fully frozen -> ice storage, else -> wetland precip
    arma::vec isFull = arma::conv_to<arma::vec>::from(
        Wetland_frozenDay_d >= param_Wetland_kup_MaxFrozenD
    );
    Wetland_iceStorage_mm      +=        isFull  % ATMOS_precipitation_mm;
    arma::vec Wetland_precipitation_mm = (1.0 - isFull) % ATMOS_precipitation_mm;

    // 3) Melt ice storage and add to wetland precipitation
    arma::vec melt = snowMelt_Factor(
        Wetland_iceStorage_mm,
        ATMOS_temperature_Cel,
        param_Wetland_kup_f,
        param_Wetland_kup_Tmelt
    );
    Wetland_iceStorage_mm      -= melt;
    Wetland_precipitation_mm   += melt;

    return Wetland_precipitation_mm;
}