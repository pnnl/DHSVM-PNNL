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
 * $Id: settings.h,v 3.1.2 2013/10/18 ning Exp $     
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

#define MAXDIRS        8
#define NNEIGHBORS     8    /* Number of directions in which water can flow based on fine grid, must equal 8 */


#define NA          -9999	/* Not applicable */

#define N_MM5_MAPS 8

#define MAP_OUTPUT 1
#define IMAGE_OUTPUT 2

#define MIN_SWE 0.005 

// Canopy type used in canopy gapping option
enum CanopyType {
  Opening,
  Forest
};

enum KEYS {
/* Options *//* list order must match order in InitConstants.c */
  format = 0, extent, gradient, flow_routing, sensible_heat_flux,
  infiltration, interpolation, mm5, qpf, prism, grid, canopy_radatt, 
  shading, snotel, outside, rhoverride, precipitation_source, wind_source, 
  temp_lapse, precip_lapse, cressman_radius, cressman_stations, prism_data_path, 
  prism_data_ext, shading_data_path, shading_data_ext, skyview_data_path, 
  stream_temp, canopy_shading, improv_radiation, gapping, snowslide, sepr, 
  snowstats, routing_neighbors, 
  /* Area */
  coordinate_system, extreme_north, extreme_west, center_latitude,
  center_longitude, time_zone_meridian, number_of_rows,
  number_of_columns, grid_spacing, point_north, point_east,
  /* Time */
  time_step, model_start, model_end, 
  /* Constants */
  ground_roughness, snow_roughness, 
  snow_water_capacity, reference_height, rain_lai_multiplier,
  snow_lai_multiplier, min_intercepted_snow, outside_basin,
  temp_lapse_rate, precip_lapse_rate, 
  max_swe, snowslide_parameter1,
  snowslide_parameter2, gapwind_adj,
  /* Constants that can vary spatially */
  rain_threshold = 0,
  snow_threshold,
  fresh_alb,
  alb_acc_lambda,
  alb_melt_lambda,
  alb_acc_min,
  alb_melt_min,
  multiplier,
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
  MM5_start = 0, MM5_temperature, MM5_humidity, MM5_wind, MM5_shortwave,
  MM5_longwave, MM5_precip, MM5_terrain, MM5_lapse, MM5_lapse_freq,
  MM5_rows, MM5_cols, MM5_ext_north, MM5_ext_west, MM5_dy, MM5_precip_dist, MM5_precip_freq,
  /* grid information */
  grid_ext_north=0, grid_ext_south, grid_ext_east, grid_ext_west, tot_grid, decim,
  grid_met_file, file_prefix, utm_zone,
  /* Soil information */
  soil_description = 0, lateral_ks, exponent, depth_thresh, max_infiltration, capillary_drive,
  soil_albedo, number_of_layers, porosity, pore_size, bubbling_pressure, field_capacity,
  wilting_point, bulk_density, vertical_ks, solids_thermal, thermal_capacity,
  /* Vegetation information */
  veg_description = 0, overstory, understory, fraction, hemifraction, trunk_space,
  aerodynamic_att, beam_attn, diff_attn, clumping_factor, leaf_angle_a, leaf_angle_b,
  scat, snow_int_cap, mass_drip_ratio, snow_int_eff, imperv_frac, detention_frac, 
  detention_decay, height, max_resistance, min_resistance, moisture_threshold, vpd, rpc,
  number_of_root_zones, root_zone_depth, overstory_fraction, understory_fraction, 
  monextn, vf_adj, overstory_monlai, understory_monlai, overstory_monalb, understory_monalb, 
  /* terrain information */
  demfile = 0, maskfile,
  soiltype_file = 0, soildepth_file, kslat_file, porosity_file,fc_file,
  vegtype_file = 0, vegfc_file, veglai_file,
  /* DHSVM channel keys */
  stream_network = 0, stream_map, stream_class, riparian_veg,
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
  graphics_variable = 0
};

#endif
