/*
 * SUMMARY:      settings.h - Definition of string, array sizes, etc.
 * USAGE:        Part of DHSVM
 *
 * AUTHOR:       Bart Nijssen
 * ORG:          University of Washington, Department of Civil Engineering
 * E-MAIL:       nijssen@u.washington.edu
 * ORIG-DATE:    Apr-1996
 * DESCRIPTION:  Definition of string, array sizes, etc.
 * DESCRIP-END.
 * FUNCTIONS:    
 * COMMENTS:
 * $Id: settings.h,v 1.18 2004/08/18 01:01:34 colleen Exp $     
 */

#ifndef SETTINGS_H
#define SETTINGS_H

#ifndef _AIX			/* AIX 3.2.5 already defines this */
typedef unsigned char uchar;
#endif
typedef unsigned short unshort;
typedef unsigned int unint;

#define MAX(x,y) ((x) > (y) ? (x) : (y))
#define MIN(x,y) ((x) < (y) ? (x) : (y))
#define INBASIN(x) ((x) != OUTSIDEBASIN)
#ifndef ABSVAL
#define ABSVAL(x)  ( (x) < 0 ? -(x) : (x) )
#endif

#ifndef TRUE
#define TRUE           1
#endif
#ifndef FALSE
#define FALSE          0
#endif

#define DHSVM_HUGE     (1e20)

/* When compiling DHSVM using Microsoft C++ define MSC++ at compile time.
   Microsoft C++ treats all numeric constants as doubles, and issues a
   warning if this constant is assigned to a float without an explicit type
   cast.   This warning can becme somewhat annoying, and the following pragma
   directive disables the warning */
#ifdef MS_C_PLUSPLUS
#pragma warning(disable : 4244)
#endif

/* Default value for not applicable */
#define NOT_APPLICABLE -9999

/* Options for precipitation and wind source */
#define RADAR          1
#define STATION        2
#define MODEL          3

/* Options for flow gradient calculation */
#define TOPOGRAPHY     1
#define WATERTABLE     2

/* Options for meterological interpolation */
#define INVDIST        1
#define NEAREST        2
#define VARCRESS       3

/* Options for model extent */
#define POINT 1
#define BASIN 2

/* Options for temperature and precipitation lapse rates */
#define CONSTANT 1
#define VARIABLE 2
#define MAP 3

/* Options for infiltration */
#define STATIC 1
#define DYNAMIC 2

/* Options for canopy radiation attenuation */
#define FIXED    1
#define VARIABLE 2

/* indicate ICE or GLACIER class */
#define GLACIER -1234

#define TINY       1e-20
#define DEBUG      FALSE

#define HEADERLINES    5
#define BUFSIZE      255
#define MAXUCHAR     255	/* Maximum value of a 1-byte uchar */
#define MAXSTRING    255
#define NAMESIZE     127

#define NDIRS          4	/* Number of directions in which water can flow, must equal 4 */
#define NNEIGHBORS      8 /* Number of directions in which water and  sediment can flow based on fine grid, must equal 8 */

#define NA          -9999	/* Not applicable */

#define N_MM5_MAPS 8

#define MAP_OUTPUT 1
#define IMAGE_OUTPUT 2

enum KEYS {
/* Options *//* list order must match order in InitConstants.c */
  format = 0, extent, gradient, flow_routing, sensible_heat_flux, sediment,
  sed_input_file, routing, infiltration, interpolation, 
  mm5, qpf, prism, canopy_radatt, 
  shading, snotel, outside, rhoverride, precipitation_source, wind_source, 
  temp_lapse, precip_lapse, cressman_radius, cressman_stations, prism_data_path, 
  prism_data_ext, shading_data_path, shading_data_ext, skyview_data_path, 
  /* Area */
  coordinate_system, extreme_north, extreme_west, center_latitude,
  center_longitude, time_zone_meridian, number_of_rows,
  number_of_columns, grid_spacing, point_north, point_east,
  /* Time */
  time_step, model_start, model_end, 
  /* Constants */
  ground_roughness, snow_roughness, rain_threshold, snow_threshold,
  snow_water_capacity, reference_height, rain_lai_multiplier,
  snow_lai_multiplier, min_intercepted_snow, outside_basin,
  temp_lapse_rate, precip_lapse_rate, 
  /* Station information */
  station_name = 0, station_north, station_east, station_elev, station_file,
  /* RADAR information */
  radar_start = 0, radar_file, radar_north, radar_west, radar_rows, radar_cols,
  radar_grid,
  /* Wind model information */
  number_of_maps = 0, wind_map_path, wind_station,
  /* precipitation lapse rate information */
  precip_lapse_rate_file = 0,
  /* MM5 information */
  MM5_start = 0,
  MM5_temperature, MM5_humidity, MM5_wind, MM5_shortwave,
  MM5_longwave, MM5_precip, MM5_terrain, MM5_lapse,
  MM5_rows, MM5_cols, MM5_ext_north, MM5_ext_west, MM5_dy,
  /* Soil information */
  soil_description = 0, lateral_ks, exponent, depth_thresh, max_infiltration, capillary_drive,
  soil_albedo, manning,
  number_of_layers, porosity, pore_size, bubbling_pressure, field_capacity,
  wilting_point, bulk_density, vertical_ks, solids_thermal, thermal_capacity,
  /* Vegetation information */
  veg_description = 0, overstory, understory, fraction, hemifraction, trunk_space,
  aerodynamic_att, radiation_att, clumping_factor, leaf_angle_a, leaf_angle_b,
  scat, snow_int_cap, mass_drip_ratio, snow_int_eff, imperv_frac, detention_frac, detention_decay, height, 
  max_resistance, min_resistance, moisture_threshold, vpd, rpc,  
  number_of_root_zones, root_zone_depth, overstory_fraction,
  understory_fraction, overstory_monlai, understory_monlai, overstory_monalb,
  understory_monalb, 
  /* Sediment vegetation information */
  root_cohesion = 0, rc_min, rc_max, rc_mean, rc_dev, rc_mode, veg_surcharge,
  vs_min, vs_max, vs_mean, vs_dev, vs_mode,
  /* terrain information */
  demfile = 0, maskfile,
  soiltype_file = 0, soildepth_file,
  /* DHSVM channel keys */
  stream_network = 0, stream_map, stream_class,
  road_network, road_map, road_class,
  /* number of each type of output */
  output_path =
    0, initial_state_path, npixels, nstates, nmapvars, nimagevars, ngraphics,
  /* pixel information */
  north = 0, east, name,
  /* state information */
  state_date = 0,
  /* map information */
  map_variable = 0, map_layer, nmaps, map_date,
  /* image information */
  image_variable = 0, image_layer, image_start, image_end, image_interval,
  image_upper, image_lower,
  /* graphics information */
  graphics_variable = 0,
  /* Sediment oconfiguration file */
  /* Sedoptions */ 
  mass_wasting = 0, surface_erosion, road_erosion, channel_routing,
  /* Parameters */
  mass_spacing, max_iterations,
// Channel Parent parameters not currently used
//  channeld50, channeld90,
  debrisd50, debrisd90,
  /* Sedtime*/
  mass_wasting_date = 0, erosion_start = 0, erosion_end,
  /* Sediment information */
  sed_description = 0, kindex, dfifty, cohesion, 
  coh_min, coh_max, coh_mean, coh_dev, coh_mode, friction_angle, fa_min, fa_max, fa_mean, 
  fa_dev, fa_mode
};

#endif
